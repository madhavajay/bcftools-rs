//! Hidden Markov model kernel used by `roh`, `cnv`, and HMM-backed plugins.
//!
//! This is a Rust port of upstream `HMM.c`'s reusable engine: transition
//! matrix precomputation, initial-state handling, snapshots for sliding
//! windows, Viterbi, forward-backward posterior probabilities, and one
//! Baum-Welch transition update.

use std::fmt;

type TransitionHook = Box<dyn FnMut(u32, u32, &mut [f64])>;

#[derive(Debug, Clone, PartialEq)]
pub struct Snapshot {
    snap_at_pos: u32,
    vit_prob: Vec<f64>,
    fwd_prob: Vec<f64>,
}

#[derive(Debug)]
pub enum HmmError {
    InvalidStateCount,
    InvalidTransitionMatrix {
        expected: usize,
        actual: usize,
    },
    InvalidInitialProbabilities {
        expected: usize,
        actual: usize,
    },
    InvalidObservationShape {
        expected: usize,
        sites: usize,
        emissions: usize,
    },
    NonIncreasingSites,
    ZeroNorm,
    TooManyStates(usize),
}

impl fmt::Display for HmmError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidStateCount => f.write_str("HMM requires at least one state"),
            Self::InvalidTransitionMatrix { expected, actual } => write!(
                f,
                "transition matrix has {actual} values, expected {expected}"
            ),
            Self::InvalidInitialProbabilities { expected, actual } => write!(
                f,
                "initial probabilities have {actual} values, expected {expected}"
            ),
            Self::InvalidObservationShape {
                expected,
                sites,
                emissions,
            } => write!(
                f,
                "{sites} sites require {expected} emission values, got {emissions}"
            ),
            Self::NonIncreasingSites => f.write_str("HMM sites must be non-decreasing"),
            Self::ZeroNorm => f.write_str("HMM probability normalization is zero"),
            Self::TooManyStates(n) => write!(f, "HMM supports at most 256 states, got {n}"),
        }
    }
}

impl std::error::Error for HmmError {}

pub struct Hmm {
    nstates: usize,
    ntprob_arr: usize,
    curr_tprob: Vec<f64>,
    tprob_arr: Vec<f64>,
    init_vit_prob: Vec<f64>,
    init_fwd_prob: Vec<f64>,
    init_bwd_prob: Vec<f64>,
    state_snap_at_pos: u32,
    state_vit_prob: Vec<f64>,
    state_fwd_prob: Vec<f64>,
    state_bwd_prob: Vec<f64>,
    vpath: Vec<u8>,
    vpath_sites: usize,
    fwd: Vec<f64>,
    fwd_sites: usize,
    snapshot_request: Option<u32>,
    set_tprob: Option<TransitionHook>,
}

impl Hmm {
    pub fn new(nstates: usize, tprob: &[f64], ntprob: usize) -> Result<Self, HmmError> {
        if nstates == 0 {
            return Err(HmmError::InvalidStateCount);
        }
        if nstates > usize::from(u8::MAX) + 1 {
            return Err(HmmError::TooManyStates(nstates));
        }

        let expected = nstates * nstates;
        if tprob.len() != expected {
            return Err(HmmError::InvalidTransitionMatrix {
                expected,
                actual: tprob.len(),
            });
        }

        let mut hmm = Self {
            nstates,
            ntprob_arr: ntprob,
            curr_tprob: vec![0.0; expected],
            tprob_arr: Vec::new(),
            init_vit_prob: vec![0.0; nstates],
            init_fwd_prob: vec![0.0; nstates],
            init_bwd_prob: vec![1.0; nstates],
            state_snap_at_pos: 0,
            state_vit_prob: vec![0.0; nstates],
            state_fwd_prob: vec![0.0; nstates],
            state_bwd_prob: vec![1.0; nstates],
            vpath: Vec::new(),
            vpath_sites: 0,
            fwd: Vec::new(),
            fwd_sites: 0,
            snapshot_request: None,
            set_tprob: None,
        };
        hmm.set_transition_probabilities(tprob, ntprob)?;
        hmm.init_states(None)?;
        Ok(hmm)
    }

    pub fn state_count(&self) -> usize {
        self.nstates
    }

    pub fn transition_probabilities(&self) -> &[f64] {
        &self.tprob_arr
    }

    pub fn set_transition_probabilities(
        &mut self,
        tprob: &[f64],
        ntprob: usize,
    ) -> Result<(), HmmError> {
        let expected = self.nstates * self.nstates;
        if tprob.len() != expected {
            return Err(HmmError::InvalidTransitionMatrix {
                expected,
                actual: tprob.len(),
            });
        }

        self.ntprob_arr = ntprob;
        let stored = ntprob.max(1);
        self.tprob_arr.resize(expected * stored, 0.0);
        self.tprob_arr[..expected].copy_from_slice(tprob);
        for i in 1..stored {
            let prev = self.tprob_arr[(i - 1) * expected..i * expected].to_vec();
            let base = self.tprob_arr[..expected].to_vec();
            let out = &mut self.tprob_arr[i * expected..(i + 1) * expected];
            multiply_matrix(self.nstates, &base, &prev, out);
        }
        Ok(())
    }

    pub fn init_states(&mut self, probs: Option<&[f64]>) -> Result<(), HmmError> {
        match probs {
            Some(probs) => {
                if probs.len() != self.nstates {
                    return Err(HmmError::InvalidInitialProbabilities {
                        expected: self.nstates,
                        actual: probs.len(),
                    });
                }
                let sum: f64 = probs.iter().sum();
                if sum == 0.0 {
                    return Err(HmmError::ZeroNorm);
                }
                for (dst, src) in self.init_vit_prob.iter_mut().zip(probs) {
                    *dst = *src / sum;
                }
            }
            None => self.init_vit_prob.fill(1.0 / self.nstates as f64),
        }
        self.init_fwd_prob.copy_from_slice(&self.init_vit_prob);
        self.init_bwd_prob.fill(1.0);
        self.reset(None);
        Ok(())
    }

    pub fn set_transition_hook<F>(&mut self, hook: F)
    where
        F: FnMut(u32, u32, &mut [f64]) + 'static,
    {
        self.set_tprob = Some(Box::new(hook));
    }

    pub fn clear_transition_hook(&mut self) {
        self.set_tprob = None;
    }

    pub fn request_snapshot(&mut self, pos: u32) {
        self.snapshot_request = Some(pos);
    }

    pub fn snapshot(&self, pos: u32) -> Snapshot {
        Snapshot {
            snap_at_pos: pos,
            vit_prob: self.state_vit_prob.clone(),
            fwd_prob: self.state_fwd_prob.clone(),
        }
    }

    pub fn restore(&mut self, snapshot: Option<&Snapshot>) {
        match snapshot {
            Some(snapshot)
                if snapshot.snap_at_pos != 0 && snapshot.vit_prob.len() == self.nstates =>
            {
                self.state_snap_at_pos = snapshot.snap_at_pos;
                self.state_vit_prob.copy_from_slice(&snapshot.vit_prob);
                self.state_fwd_prob.copy_from_slice(&snapshot.fwd_prob);
            }
            _ => self.reset(None),
        }
    }

    pub fn reset(&mut self, snapshot: Option<&mut Snapshot>) {
        if let Some(snapshot) = snapshot {
            snapshot.snap_at_pos = 0;
        }
        self.state_snap_at_pos = 0;
        self.state_vit_prob.copy_from_slice(&self.init_vit_prob);
        self.state_fwd_prob.copy_from_slice(&self.init_fwd_prob);
        self.state_bwd_prob.copy_from_slice(&self.init_bwd_prob);
    }

    pub fn run_viterbi(&mut self, emissions: &[f64], sites: &[u32]) -> Result<&[u8], HmmError> {
        self.validate_observations(emissions, sites)?;
        let n = sites.len();
        self.vpath.resize(n * self.nstates, 0);
        self.vpath_sites = n;

        let mut vprob = self.state_vit_prob.clone();
        let mut vprob_tmp = vec![0.0; self.nstates];
        let mut prev_pos = if self.state_snap_at_pos != 0 {
            self.state_snap_at_pos
        } else {
            sites[0]
        };

        for (i, (&pos, eprob)) in sites.iter().zip(emissions.chunks(self.nstates)).enumerate() {
            let pos_diff = pos_diff_forward(prev_pos, pos)?;
            self.set_current_transition_probabilities(pos_diff, prev_pos, pos);
            prev_pos = pos;

            let mut vnorm = 0.0;
            for j in 0..self.nstates {
                let mut vmax = 0.0;
                let mut k_vmax = 0usize;
                for (k, &prob) in vprob.iter().enumerate() {
                    let pval = prob * mat(&self.curr_tprob, self.nstates, j, k);
                    if vmax < pval {
                        vmax = pval;
                        k_vmax = k;
                    }
                }
                self.vpath[i * self.nstates + j] = k_vmax as u8;
                vprob_tmp[j] = vmax * eprob[j];
                vnorm += vprob_tmp[j];
            }
            normalize(&mut vprob_tmp, vnorm)?;
            std::mem::swap(&mut vprob, &mut vprob_tmp);

            if self.snapshot_request == Some(pos) {
                self.state_vit_prob.copy_from_slice(&vprob);
            }
        }

        let mut iptr = (0..self.nstates)
            .max_by(|&a, &b| vprob[a].total_cmp(&vprob[b]))
            .unwrap_or(0);
        for i in (0..n).rev() {
            let iptr_prev = usize::from(self.vpath[i * self.nstates + iptr]);
            self.vpath[i * self.nstates] = iptr as u8;
            iptr = iptr_prev;
        }

        Ok(self.viterbi_path())
    }

    pub fn viterbi_path(&self) -> &[u8] {
        &self.vpath[..self.vpath_sites * self.nstates]
    }

    pub fn viterbi_states(&self) -> Vec<u8> {
        (0..self.vpath_sites)
            .map(|i| self.vpath[i * self.nstates])
            .collect()
    }

    pub fn run_forward_backward(
        &mut self,
        emissions: &[f64],
        sites: &[u32],
    ) -> Result<&[f64], HmmError> {
        self.validate_observations(emissions, sites)?;
        let n = sites.len();
        self.fwd.resize((n + 1) * self.nstates, 0.0);
        self.fwd_sites = n;
        self.fwd[..self.nstates].copy_from_slice(&self.state_fwd_prob);

        let mut bwd = self.state_bwd_prob.clone();
        let mut bwd_tmp = vec![0.0; self.nstates];
        let mut prev_pos = if self.state_snap_at_pos != 0 {
            self.state_snap_at_pos
        } else {
            sites[0]
        };

        for (i, (&pos, eprob)) in sites.iter().zip(emissions.chunks(self.nstates)).enumerate() {
            let pos_diff = pos_diff_forward(prev_pos, pos)?;
            self.set_current_transition_probabilities(pos_diff, prev_pos, pos);
            prev_pos = pos;

            let (left, right) = self.fwd.split_at_mut((i + 1) * self.nstates);
            let fwd_prev = &left[i * self.nstates..(i + 1) * self.nstates];
            let fwd = &mut right[..self.nstates];
            forward_step(self.nstates, &self.curr_tprob, fwd_prev, eprob, fwd)?;

            if self.snapshot_request == Some(pos) {
                self.state_fwd_prob.copy_from_slice(fwd);
            }
        }

        prev_pos = sites[n - 1];
        for i in 0..n {
            let site_index = n - i - 1;
            let pos = sites[site_index];
            let eprob = &emissions[site_index * self.nstates..(site_index + 1) * self.nstates];
            let pos_diff = pos_diff_backward(pos, prev_pos)?;
            self.set_current_transition_probabilities(pos_diff, pos, prev_pos);
            prev_pos = pos;

            backward_step(self.nstates, &self.curr_tprob, &bwd, eprob, &mut bwd_tmp)?;
            let fwd =
                &mut self.fwd[(site_index + 1) * self.nstates..(site_index + 2) * self.nstates];
            let mut norm = 0.0;
            for j in 0..self.nstates {
                fwd[j] *= bwd[j];
                norm += fwd[j];
            }
            normalize(fwd, norm)?;
            std::mem::swap(&mut bwd, &mut bwd_tmp);
        }

        Ok(self.posterior_probabilities())
    }

    pub fn posterior_probabilities(&self) -> &[f64] {
        let start = self.nstates;
        let end = (self.fwd_sites + 1) * self.nstates;
        &self.fwd[start..end]
    }

    pub fn run_baum_welch(
        &mut self,
        emissions: &[f64],
        sites: &[u32],
    ) -> Result<Vec<f64>, HmmError> {
        self.validate_observations(emissions, sites)?;
        let n = sites.len();
        self.fwd.resize((n + 1) * self.nstates, 0.0);
        self.fwd[..self.nstates].copy_from_slice(&self.state_fwd_prob);

        let mut bwd = self.state_bwd_prob.clone();
        let mut bwd_tmp = vec![0.0; self.nstates];
        let mut tmp_xi = vec![0.0; self.nstates * self.nstates];
        let mut tmp_gamma = vec![0.0; self.nstates];
        let mut fwd_bwd = vec![0.0; self.nstates];
        let mut prev_pos = if self.state_snap_at_pos != 0 {
            self.state_snap_at_pos
        } else {
            sites[0]
        };

        for (i, (&pos, eprob)) in sites.iter().zip(emissions.chunks(self.nstates)).enumerate() {
            let pos_diff = pos_diff_forward(prev_pos, pos)?;
            self.set_current_transition_probabilities(pos_diff, prev_pos, pos);
            prev_pos = pos;
            let (left, right) = self.fwd.split_at_mut((i + 1) * self.nstates);
            let fwd_prev = &left[i * self.nstates..(i + 1) * self.nstates];
            let fwd = &mut right[..self.nstates];
            forward_step(self.nstates, &self.curr_tprob, fwd_prev, eprob, fwd)?;
        }

        prev_pos = sites[n - 1];
        for i in 0..n {
            let site_index = n - i - 1;
            let pos = sites[site_index];
            let eprob = &emissions[site_index * self.nstates..(site_index + 1) * self.nstates];
            let pos_diff = pos_diff_backward(pos, prev_pos)?;
            self.set_current_transition_probabilities(pos_diff, pos, prev_pos);
            prev_pos = pos;

            backward_step(self.nstates, &self.curr_tprob, &bwd, eprob, &mut bwd_tmp)?;

            let fwd =
                &mut self.fwd[(site_index + 1) * self.nstates..(site_index + 2) * self.nstates];
            let mut norm = 0.0;
            for j in 0..self.nstates {
                fwd_bwd[j] = fwd[j] * bwd_tmp[j];
                norm += fwd_bwd[j];
            }
            normalize(&mut fwd_bwd, norm)?;
            for j in 0..self.nstates {
                tmp_gamma[j] += fwd_bwd[j];
            }

            for (j, &fwd_j) in fwd.iter().enumerate().take(self.nstates) {
                for k in 0..self.nstates {
                    let value =
                        fwd_j * bwd[k] * mat(&self.tprob_arr, self.nstates, k, j) * eprob[k] / norm;
                    let idx = mat_idx(self.nstates, k, j);
                    tmp_xi[idx] += value;
                }
            }

            fwd.copy_from_slice(&fwd_bwd);
            std::mem::swap(&mut bwd, &mut bwd_tmp);
        }

        let mut out = vec![0.0; self.nstates * self.nstates];
        for j in 0..self.nstates {
            let mut norm = 0.0;
            for k in 0..self.nstates {
                let idx = mat_idx(self.nstates, k, j);
                out[idx] = tmp_xi[idx] / tmp_gamma[j];
                norm += out[idx];
            }
            for k in 0..self.nstates {
                out[mat_idx(self.nstates, k, j)] /= norm;
            }
        }
        self.curr_tprob.copy_from_slice(&out);
        Ok(out)
    }

    fn validate_observations(&self, emissions: &[f64], sites: &[u32]) -> Result<(), HmmError> {
        let expected = sites.len() * self.nstates;
        if sites.is_empty() || emissions.len() != expected {
            return Err(HmmError::InvalidObservationShape {
                expected,
                sites: sites.len(),
                emissions: emissions.len(),
            });
        }
        if sites.windows(2).any(|w| w[1] < w[0]) {
            return Err(HmmError::NonIncreasingSites);
        }
        Ok(())
    }

    fn set_current_transition_probabilities(&mut self, pos_diff: u32, prev_pos: u32, pos: u32) {
        let matrix_len = self.nstates * self.nstates;
        let n = (pos_diff as usize)
            .checked_rem(self.ntprob_arr)
            .unwrap_or(0);
        self.curr_tprob
            .copy_from_slice(&self.tprob_arr[n * matrix_len..(n + 1) * matrix_len]);

        if let Some(full_blocks) = (pos_diff as usize).checked_div(self.ntprob_arr) {
            let block = self.tprob_arr
                [(self.ntprob_arr - 1) * matrix_len..self.ntprob_arr * matrix_len]
                .to_vec();
            for _ in 0..full_blocks {
                let current = self.curr_tprob.clone();
                multiply_matrix(self.nstates, &block, &current, &mut self.curr_tprob);
            }
        }

        if let Some(hook) = self.set_tprob.as_mut() {
            hook(prev_pos, pos, &mut self.curr_tprob);
        }
    }
}

fn pos_diff_forward(prev_pos: u32, pos: u32) -> Result<u32, HmmError> {
    if pos < prev_pos {
        return Err(HmmError::NonIncreasingSites);
    }
    Ok(if pos == prev_pos {
        0
    } else {
        pos - prev_pos - 1
    })
}

fn pos_diff_backward(pos: u32, prev_pos: u32) -> Result<u32, HmmError> {
    if pos > prev_pos {
        return Err(HmmError::NonIncreasingSites);
    }
    Ok(if pos == prev_pos {
        0
    } else {
        prev_pos - pos - 1
    })
}

fn mat_idx(n: usize, i: usize, j: usize) -> usize {
    n * i + j
}

fn mat(matrix: &[f64], n: usize, i: usize, j: usize) -> f64 {
    matrix[mat_idx(n, i, j)]
}

fn multiply_matrix(n: usize, a: &[f64], b: &[f64], dst: &mut [f64]) {
    for i in 0..n {
        for j in 0..n {
            let mut val = 0.0;
            for k in 0..n {
                val += mat(a, n, i, k) * mat(b, n, k, j);
            }
            dst[mat_idx(n, i, j)] = val;
        }
    }
}

fn normalize(values: &mut [f64], norm: f64) -> Result<(), HmmError> {
    if norm == 0.0 || !norm.is_finite() {
        return Err(HmmError::ZeroNorm);
    }
    for value in values {
        *value /= norm;
    }
    Ok(())
}

fn forward_step(
    n: usize,
    tprob: &[f64],
    fwd_prev: &[f64],
    eprob: &[f64],
    fwd: &mut [f64],
) -> Result<(), HmmError> {
    let mut norm = 0.0;
    for j in 0..n {
        let mut pval = 0.0;
        for (k, &prev) in fwd_prev.iter().enumerate().take(n) {
            pval += prev * mat(tprob, n, j, k);
        }
        fwd[j] = pval * eprob[j];
        norm += fwd[j];
    }
    normalize(fwd, norm)
}

fn backward_step(
    n: usize,
    tprob: &[f64],
    bwd: &[f64],
    eprob: &[f64],
    bwd_tmp: &mut [f64],
) -> Result<(), HmmError> {
    let mut norm = 0.0;
    for (j, out) in bwd_tmp.iter_mut().enumerate().take(n) {
        let mut pval = 0.0;
        for k in 0..n {
            pval += bwd[k] * eprob[k] * mat(tprob, n, k, j);
        }
        *out = pval;
        norm += pval;
    }
    normalize(bwd_tmp, norm)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-9, "{a} != {b}");
    }

    #[test]
    fn transition_probabilities_are_precomputed_as_matrix_powers() {
        let hmm = Hmm::new(2, &[0.9, 0.2, 0.1, 0.8], 3).unwrap();
        let tprob = hmm.transition_probabilities();
        assert_eq!(tprob.len(), 12);
        approx_eq(tprob[4], 0.83);
        approx_eq(tprob[5], 0.34);
        approx_eq(tprob[6], 0.17);
        approx_eq(tprob[7], 0.66);
    }

    #[test]
    fn viterbi_tracks_the_most_likely_state_path() {
        let mut hmm = Hmm::new(2, &[0.95, 0.05, 0.05, 0.95], 0).unwrap();
        let sites = [10, 11, 12, 13];
        let emissions = [0.99, 0.01, 0.95, 0.05, 0.05, 0.95, 0.01, 0.99];
        hmm.run_viterbi(&emissions, &sites).unwrap();
        assert_eq!(hmm.viterbi_states(), [0, 0, 1, 1]);
    }

    #[test]
    fn forward_backward_returns_normalized_posteriors() {
        let mut hmm = Hmm::new(2, &[0.9, 0.1, 0.1, 0.9], 0).unwrap();
        let sites = [1, 2, 3];
        let emissions = [0.8, 0.2, 0.5, 0.5, 0.2, 0.8];
        let posterior = hmm.run_forward_backward(&emissions, &sites).unwrap();
        assert_eq!(posterior.len(), 6);
        for probs in posterior.chunks(2) {
            approx_eq(probs[0] + probs[1], 1.0);
        }
        assert!(posterior[0] > posterior[1]);
        assert!(posterior[4] < posterior[5]);
    }

    #[test]
    fn baum_welch_returns_column_normalized_transition_matrix() {
        let mut hmm = Hmm::new(2, &[0.8, 0.3, 0.2, 0.7], 0).unwrap();
        let sites = [1, 2, 3, 4];
        let emissions = [0.9, 0.1, 0.85, 0.15, 0.2, 0.8, 0.1, 0.9];
        let tprob = hmm.run_baum_welch(&emissions, &sites).unwrap();
        for j in 0..2 {
            approx_eq(tprob[mat_idx(2, 0, j)] + tprob[mat_idx(2, 1, j)], 1.0);
        }
    }

    #[test]
    fn snapshot_restore_restarts_from_saved_probabilities() {
        let mut hmm = Hmm::new(2, &[0.9, 0.1, 0.1, 0.9], 0).unwrap();
        hmm.request_snapshot(2);
        let sites = [1, 2, 3];
        let emissions = [0.9, 0.1, 0.8, 0.2, 0.1, 0.9];
        hmm.run_forward_backward(&emissions, &sites).unwrap();
        let snapshot = hmm.snapshot(2);
        hmm.restore(Some(&snapshot));
        assert_eq!(hmm.state_snap_at_pos, 2);
    }
}
