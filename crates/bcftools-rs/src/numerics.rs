//! Numeric helpers shared by analysis subcommands.
//!
//! This module ports the self-contained numeric kernels shared by bcftools:
//! logarithmic histograms (`dist.c`), minimizers (`kmin.c`), complete-linkage
//! clustering (`hclust.c`), peak fitting (`peakfit.c`), single-locus EM
//! (`em.c`), and the consensus-caller allele-frequency posterior (`prob1.c`).

#[derive(Debug, Clone)]
pub struct LogDistribution {
    bins: Vec<u64>,
    nvalues: u64,
    npow: u32,
    nexact: u32,
    nlevel: u32,
}

impl LogDistribution {
    pub fn new(npow: u32) -> Self {
        let nexact = 10_u32.pow(npow);
        let nlevel = nexact - 10_u32.pow(npow.saturating_sub(1));
        Self {
            bins: Vec::new(),
            nvalues: 0,
            npow,
            nexact,
            nlevel,
        }
    }

    pub fn nbins(&self) -> usize {
        self.bins.len()
    }

    pub fn nvalues(&self) -> u64 {
        self.nvalues
    }

    pub fn insert(&mut self, value: u32) -> usize {
        let ibin = self.bin_for_value(value);
        if ibin >= self.bins.len() {
            self.bins.resize(ibin + 1, 0);
        }
        self.bins[ibin] += 1;
        self.nvalues += 1;
        ibin
    }

    pub fn insert_n(&mut self, value: u32, count: u32) -> usize {
        if count == 0 {
            return 0;
        }
        let ibin = self.insert(value);
        self.bins[ibin] += u64::from(count - 1);
        self.nvalues += u64::from(count);
        ibin
    }

    pub fn get(&self, idx: usize) -> Option<DistributionBin> {
        let count = *self.bins.get(idx)?;
        let idx = u32::try_from(idx).ok()?;
        let (beg, end) = if idx < self.nexact {
            (idx, idx + 1)
        } else {
            let level = (idx - self.nexact) / self.nlevel + 1;
            let bin = idx - self.nexact - self.nlevel * (level - 1);
            let step = 10_u32.pow(level);
            let value = 10_u32.pow(level + self.npow - 1) + step * bin;
            (value, value + step)
        };
        Some(DistributionBin { beg, end, count })
    }

    fn bin_for_value(&self, value: u32) -> usize {
        if value <= self.nexact {
            return value as usize;
        }
        let value_f = f64::from(value);
        let npow = value_f.log10() as u32;
        let level = npow - self.npow + 1;
        let step = 10_u32.pow(level);
        let bin = self.nexact + self.nlevel * (level - 1) + (value - 10_u32.pow(npow)) / step;
        bin as usize
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DistributionBin {
    pub beg: u32,
    pub end: u32,
    pub count: u64,
}

pub fn hooke_jeeves<F>(mut func: F, x: &mut [f64], r: f64, eps: f64, max_calls: usize) -> f64
where
    F: FnMut(&[f64]) -> f64,
{
    let n = x.len();
    let mut x1 = vec![0.0; n];
    let mut dx = vec![0.0; n];
    for (step, &value) in dx.iter_mut().zip(x.iter()) {
        *step = value.abs() * r;
        if *step == 0.0 {
            *step = r;
        }
    }

    let mut n_calls = 1usize;
    let mut radius = r;
    let mut fx = func(x);
    let mut fx1;
    loop {
        x1.copy_from_slice(x);
        fx1 = hooke_jeeves_aux(&mut func, &mut x1, fx, &mut dx, &mut n_calls);
        while fx1 < fx {
            for k in 0..n {
                let t = x[k];
                dx[k] = if x1[k] > x[k] {
                    dx[k].abs()
                } else {
                    -dx[k].abs()
                };
                x[k] = x1[k];
                x1[k] = x1[k] + x1[k] - t;
            }
            fx = fx1;
            if n_calls >= max_calls {
                break;
            }
            fx1 = func(&x1);
            n_calls += 1;
            fx1 = hooke_jeeves_aux(&mut func, &mut x1, fx1, &mut dx, &mut n_calls);
            if fx1 >= fx {
                break;
            }
            if x1
                .iter()
                .zip(x.iter())
                .zip(dx.iter())
                .all(|((&x1, &x), &dx)| (x1 - x).abs() <= 0.5 * dx.abs())
            {
                break;
            }
        }
        if radius >= eps {
            if n_calls >= max_calls {
                break;
            }
            radius *= r;
            for step in &mut dx {
                *step *= r;
            }
        } else {
            break;
        }
    }
    fx1
}

fn hooke_jeeves_aux<F>(
    func: &mut F,
    x1: &mut [f64],
    mut fx1: f64,
    dx: &mut [f64],
    n_calls: &mut usize,
) -> f64
where
    F: FnMut(&[f64]) -> f64,
{
    for k in 0..x1.len() {
        x1[k] += dx[k];
        let mut ftmp = func(x1);
        *n_calls += 1;
        if ftmp < fx1 {
            fx1 = ftmp;
        } else {
            dx[k] = -dx[k];
            x1[k] += dx[k] + dx[k];
            ftmp = func(x1);
            *n_calls += 1;
            if ftmp < fx1 {
                fx1 = ftmp;
            } else {
                x1[k] -= dx[k];
            }
        }
    }
    fx1
}

pub fn brent<F>(mut func: F, mut a: f64, mut b: f64, tol: f64) -> (f64, f64)
where
    F: FnMut(f64) -> f64,
{
    let gold1 = 1.618_033_988_7;
    let gold2 = 0.381_966_011_3;
    let tiny = 1e-20;
    let max_iter = 100;

    let mut fa = func(a);
    let mut fb = func(b);
    if fb > fa {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut fa, &mut fb);
    }
    let mut c = b + gold1 * (b - a);
    let mut fc = func(c);
    while fb > fc {
        let bound = b + 100.0 * (c - b);
        let r = (b - a) * (fb - fc);
        let q = (b - c) * (fb - fa);
        let tmp = if (q - r).abs() < tiny {
            if q > r { tiny } else { -tiny }
        } else {
            q - r
        };
        let mut u = b - ((b - c) * q - (b - a) * r) / (2.0 * tmp);
        let mut fu;
        if between(u, b, c) {
            fu = func(u);
            if fu < fc {
                a = b;
                b = u;
                fb = fu;
                break;
            } else if fu > fb {
                c = u;
                break;
            }
            u = c + gold1 * (c - b);
            fu = func(u);
        } else if between(u, c, bound) {
            fu = func(u);
            if fu < fc {
                b = c;
                c = u;
                u = c + gold1 * (c - b);
                fb = fc;
                fc = fu;
                let _ = func(u);
            } else {
                a = b;
                b = c;
                c = u;
                fb = fc;
                break;
            }
        } else if between(bound, c, u) {
            u = bound;
            fu = func(u);
        } else {
            u = c + gold1 * (c - b);
            fu = func(u);
        }
        a = b;
        b = c;
        c = u;
        fa = fb;
        fb = fc;
        fc = fu;
    }
    if a > c {
        std::mem::swap(&mut a, &mut c);
    }

    let mut e: f64 = 0.0;
    let mut d: f64 = 0.0;
    let mut w = b;
    let mut v = b;
    let mut fv = fb;
    let mut fw = fb;
    for _ in 0..max_iter {
        let mid = 0.5 * (a + c);
        let tol1 = tol * b.abs() + tiny;
        let tol2 = 2.0 * tol1;
        if (b - mid).abs() <= tol2 - 0.5 * (c - a) {
            return (b, fb);
        }
        if e.abs() > tol1 {
            let r = (b - w) * (fb - fv);
            let mut q = (b - v) * (fb - fw);
            let mut p = (b - v) * q - (b - w) * r;
            q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            } else {
                q = -q;
            }
            let eold = e;
            e = d;
            if p.abs() >= (0.5 * q * eold).abs() || p <= q * (a - b) || p >= q * (c - b) {
                e = if b >= mid { a - b } else { c - b };
                d = gold2 * e;
            } else {
                d = p / q;
                let u = b + d;
                if u - a < tol2 || c - u < tol2 {
                    d = if mid > b { tol1 } else { -tol1 };
                }
            }
        } else {
            e = if b >= mid { a - b } else { c - b };
            d = gold2 * e;
        }
        let u = if d.abs() >= tol1 {
            b + d
        } else {
            b + if d > 0.0 { tol1 } else { -tol1 }
        };
        let fu = func(u);
        if fu <= fb {
            if u >= b {
                a = b;
            } else {
                c = b;
            }
            v = w;
            w = b;
            b = u;
            fv = fw;
            fw = fb;
            fb = fu;
        } else {
            if u < b {
                a = u;
            } else {
                c = u;
            }
            if fu <= fw || w == b {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if fu <= fv || v == b || v == w {
                v = u;
                fv = fu;
            }
        }
    }
    (b, fb)
}

fn between(u: f64, a: f64, b: f64) -> bool {
    (a > u && u > b) || (a < u && u < b)
}

const EM_ITER_MAX: usize = 50;
const EM_ITER_TRY: usize = 10;
const EM_EPS: f64 = 1e-5;
const MC_DEF_INDEL: f64 = 0.15;
const MC_TINY: f64 = 1e-20;
const CONTRAST_TINY: f64 = 1e-30;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PriorType {
    Full,
    Cond2,
    Flat,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SingleLocusEm {
    pub ref_frequency: f64,
    pub genotype_frequencies: [f64; 3],
    pub hwe_p_value: f64,
    pub group_frequencies: Option<[f64; 2]>,
    pub one_degree_p_value: Option<f64>,
    pub two_degree_p_value: Option<f64>,
}

pub fn single_locus_em(
    pdg: &[[f64; 3]],
    n1: Option<usize>,
    compute_hwe: bool,
    compute_group_lrt: bool,
) -> Option<SingleLocusEm> {
    if pdg.is_empty() {
        return None;
    }
    let mut ref_frequency = est_freq(pdg)?;
    ref_frequency = freqml(ref_frequency, 0, pdg.len(), pdg);

    let mut genotype_frequencies = [
        (1.0 - ref_frequency) * (1.0 - ref_frequency),
        2.0 * ref_frequency * (1.0 - ref_frequency),
        ref_frequency * ref_frequency,
    ];
    let mut hwe_p_value = -1.0;
    if compute_hwe || compute_group_lrt {
        let hwe = genotype_frequencies;
        for _ in 0..EM_ITER_MAX {
            if g3_iter(&mut genotype_frequencies, pdg, 0, pdg.len()) < EM_EPS {
                break;
            }
        }
        if compute_hwe {
            let mut ratio = 1.0;
            for p in pdg {
                ratio *= (p[0] * genotype_frequencies[0]
                    + p[1] * genotype_frequencies[1]
                    + p[2] * genotype_frequencies[2])
                    / (p[0] * hwe[0] + p[1] * hwe[1] + p[2] * hwe[2]);
            }
            hwe_p_value = gammaq(0.5, ratio.ln());
        }
    }

    let n1 = n1.unwrap_or(0);
    let mut group_frequencies = None;
    let mut one_degree_p_value = None;
    let mut two_degree_p_value = None;
    if n1 > 0 && n1 < pdg.len() {
        let g1 = freqml(ref_frequency, 0, n1, pdg);
        let g2 = freqml(ref_frequency, n1, pdg.len(), pdg);
        group_frequencies = Some([g1, g2]);

        if compute_group_lrt {
            let freqs = [ref_frequency, g1, g2];
            let mut f3 = [[0.0; 3]; 3];
            for (out, &f) in f3.iter_mut().zip(freqs.iter()) {
                out[0] = (1.0 - f) * (1.0 - f);
                out[1] = 2.0 * f * (1.0 - f);
                out[2] = f * f;
            }
            let tmp = lk_ratio_test(pdg, n1, &f3).ln().max(0.0);
            one_degree_p_value = Some(gammaq(0.5, tmp));

            let mut g = [genotype_frequencies; 3];
            for _ in 0..EM_ITER_MAX {
                if g3_iter(&mut g[1], pdg, 0, n1) < EM_EPS {
                    break;
                }
            }
            for _ in 0..EM_ITER_MAX {
                if g3_iter(&mut g[2], pdg, n1, pdg.len()) < EM_EPS {
                    break;
                }
            }
            let tmp = lk_ratio_test(pdg, n1, &g).ln().max(0.0);
            two_degree_p_value = Some(gammaq(1.0, tmp));
        }
    }

    Some(SingleLocusEm {
        ref_frequency,
        genotype_frequencies,
        hwe_p_value,
        group_frequencies,
        one_degree_p_value,
        two_degree_p_value,
    })
}

fn est_freq(pdg: &[[f64; 3]]) -> Option<f64> {
    let mut counts = [0usize; 3];
    for p in pdg {
        if p[0] != 1.0 || p[1] != 1.0 || p[2] != 1.0 {
            let which = if p[0] > p[1] { 0 } else { 1 };
            let which = if p[which] > p[2] { which } else { 2 };
            counts[which] += 1;
        }
    }
    let n = counts.iter().sum::<usize>();
    (n != 0).then(|| (0.5 * counts[1] as f64 + counts[2] as f64) / n as f64)
}

fn em_neg_log_likelihood(f: f64, beg: usize, end: usize, pdg: &[[f64; 3]]) -> f64 {
    if !(0.0..=1.0).contains(&f) {
        return 1e300;
    }
    let f3 = [(1.0 - f) * (1.0 - f), 2.0 * f * (1.0 - f), f * f];
    let mut product = 1.0;
    let mut loss = 0.0;
    for p in &pdg[beg..end] {
        product *= p[0] * f3[0] + p[1] * f3[1] + p[2] * f3[2];
        if product < 1e-200 {
            loss -= product.ln();
            product = 1.0;
        }
    }
    loss - product.ln()
}

fn freq_iter(f: &mut f64, pdg: &[[f64; 3]], beg: usize, end: usize) -> f64 {
    let f0 = *f;
    let f3 = [(1.0 - f0) * (1.0 - f0), 2.0 * f0 * (1.0 - f0), f0 * f0];
    let mut next = 0.0;
    for p in &pdg[beg..end] {
        next += (p[1] * f3[1] + 2.0 * p[2] * f3[2]) / (p[0] * f3[0] + p[1] * f3[1] + p[2] * f3[2]);
    }
    next /= (end - beg) as f64 * 2.0;
    let err = (next - *f).abs();
    *f = next;
    err
}

fn freqml(f0: f64, beg: usize, end: usize, pdg: &[[f64; 3]]) -> f64 {
    let mut f = f0;
    let mut converged = false;
    for _ in 0..EM_ITER_TRY {
        if freq_iter(&mut f, pdg, beg, end) < EM_EPS {
            converged = true;
            break;
        }
    }
    if !converged {
        let left = if f0 == f { 0.5 * f0 } else { f0 };
        f = brent(|x| em_neg_log_likelihood(x, beg, end, pdg), left, f, EM_EPS).0;
    }
    f
}

fn g3_iter(g: &mut [f64; 3], pdg: &[[f64; 3]], beg: usize, end: usize) -> f64 {
    let mut next = [0.0; 3];
    for p in &pdg[beg..end] {
        let tmp = [p[0] * g[0], p[1] * g[1], p[2] * g[2]];
        let sum = (tmp[0] + tmp[1] + tmp[2]) * (end - beg) as f64;
        next[0] += tmp[0] / sum;
        next[1] += tmp[1] / sum;
        next[2] += tmp[2] / sum;
    }
    let err = next
        .iter()
        .zip(g.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    *g = next;
    err
}

fn lk_ratio_test(pdg: &[[f64; 3]], n1: usize, f3: &[[f64; 3]; 3]) -> f64 {
    let mut ratio = 1.0;
    for (i, p) in pdg.iter().enumerate() {
        let group = if i < n1 { 1 } else { 2 };
        ratio *= (p[0] * f3[group][0] + p[1] * f3[group][1] + p[2] * f3[group][2])
            / (p[0] * f3[0][0] + p[1] * f3[0][1] + p[2] * f3[0][2]);
    }
    ratio
}

#[derive(Debug, Clone)]
pub struct Prob1 {
    n: usize,
    m: usize,
    n1: Option<usize>,
    ploidy: Option<Vec<u8>>,
    phi: Vec<f64>,
    phi_indel: Vec<f64>,
    phi1: Vec<f64>,
    phi2: Vec<f64>,
    z: Vec<f64>,
    zswap: Vec<f64>,
    z1: Vec<f64>,
    z2: Vec<f64>,
    afs: Vec<f64>,
    afs1: Vec<f64>,
    p_ref_folded: f64,
    p_var_folded: f64,
    t: f64,
    t1: f64,
    t2: f64,
    hypergeom: Option<Vec<Vec<f64>>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Prob1Result {
    pub ac: usize,
    pub f_exp: f64,
    pub f_flat: f64,
    pub p_ref_folded: f64,
    pub p_ref: f64,
    pub p_var_folded: f64,
    pub p_var: f64,
    pub cil: f64,
    pub cih: f64,
    pub cmp: [f64; 3],
    pub p_chi2: f64,
    pub lrt: f64,
}

impl Prob1 {
    pub fn new(n_samples: usize, ploidy: Option<Vec<u8>>) -> Self {
        let mut m = 2 * n_samples;
        let ploidy = ploidy.filter(|p| {
            assert_eq!(p.len(), n_samples);
            m = p.iter().map(|&x| usize::from(x)).sum();
            m != 2 * n_samples
        });
        let mut this = Self {
            n: n_samples,
            m,
            n1: None,
            ploidy,
            phi: vec![0.0; m + 1],
            phi_indel: vec![0.0; m + 1],
            phi1: vec![0.0; m + 1],
            phi2: vec![0.0; m + 1],
            z: vec![0.0; m + 1],
            zswap: vec![0.0; m + 1],
            z1: vec![0.0; m + 1],
            z2: vec![0.0; m + 1],
            afs: vec![0.0; m + 1],
            afs1: vec![0.0; m + 1],
            p_ref_folded: 0.0,
            p_var_folded: 0.0,
            t: 0.0,
            t1: 0.0,
            t2: 0.0,
            hypergeom: None,
        };
        this.init_prior(PriorType::Full, 1e-3);
        this
    }

    pub fn init_prior(&mut self, prior_type: PriorType, theta: f64) {
        init_prior(prior_type, theta, self.m, &mut self.phi);
        self.indel_prior(MC_DEF_INDEL);
    }

    pub fn init_subprior(&mut self, prior_type: PriorType, theta: f64) {
        let Some(n1) = self.n1 else {
            return;
        };
        if n1 == 0 || n1 >= self.n {
            return;
        }
        init_prior(prior_type, theta, 2 * n1, &mut self.phi1);
        init_prior(prior_type, theta, 2 * (self.n - n1), &mut self.phi2);
    }

    pub fn set_n1(&mut self, n1: usize) -> Result<(), &'static str> {
        if n1 == 0 || n1 >= self.n {
            return Err("n1 must split the sample set");
        }
        if self.m != self.n * 2 {
            return Err("n1 cannot be set when haploid samples are present");
        }
        self.n1 = Some(n1);
        self.hypergeom = None;
        Ok(())
    }

    pub fn calculate(
        &mut self,
        pdg: &[[f64; 3]],
        is_indel: bool,
        do_contrast: bool,
    ) -> Option<Prob1Result> {
        if pdg.len() != self.n || self.n == 0 {
            return None;
        }
        let f_exp = self.cal_afs(pdg, is_indel)?;
        let p_ref = self.afs1[self.m];
        let p_var = self.afs1[..self.m].iter().sum();
        let ac = self
            .z
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .map(|(k, _)| self.m - k)
            .unwrap_or(0);

        let sum_z: f64 = self.z.iter().sum();
        let f_flat = self
            .z
            .iter()
            .enumerate()
            .map(|(k, z)| k as f64 * z / sum_z)
            .sum::<f64>()
            / self.m as f64;

        let mut lower = 0usize;
        let mut accum = 0.0;
        while lower <= self.m && accum + self.afs1[lower] <= 0.025 {
            accum += self.afs1[lower];
            lower += 1;
        }
        let mut upper = self.m;
        accum = 0.0;
        while upper > 0 && accum + self.afs1[upper] <= 0.025 {
            accum += self.afs1[upper];
            upper -= 1;
        }
        let cil = (self.m - upper) as f64 / self.m as f64;
        let cih = (self.m - lower) as f64 / self.m as f64;

        let lrt = if let Some(n1) = self.n1 {
            let max0 = self.z.iter().copied().fold(-1.0, f64::max);
            let max1 = self.z1[..=2 * n1].iter().copied().fold(-1.0, f64::max);
            let max2 = self.z2[..=self.m - 2 * n1]
                .iter()
                .copied()
                .fold(-1.0, f64::max);
            let value = (max1 * max2 / max0).ln();
            if value < 0.0 { 1.0 } else { gammaq(0.5, value) }
        } else {
            -1.0
        };

        let mut cmp = [-1.0; 3];
        let p_chi2 = if do_contrast && p_var > 0.5 {
            self.contrast2(&mut cmp)
        } else {
            -1.0
        };

        Some(Prob1Result {
            ac,
            f_exp,
            f_flat,
            p_ref_folded: self.p_ref_folded,
            p_ref,
            p_var_folded: self.p_var_folded,
            p_var,
            cil,
            cih,
            cmp,
            p_chi2,
            lrt,
        })
    }

    pub fn call_gt(&self, pdg: [f64; 3], ref_frequency: f64, sample: usize, is_var: bool) -> u8 {
        let ploidy = self.ploidy.as_ref().map(|p| p[sample]).unwrap_or(2);
        let priors = if ploidy == 2 {
            [
                (1.0 - ref_frequency) * (1.0 - ref_frequency),
                2.0 * ref_frequency * (1.0 - ref_frequency),
                ref_frequency * ref_frequency,
            ]
        } else {
            [1.0 - ref_frequency, 0.0, ref_frequency]
        };
        let mut g = [pdg[0] * priors[0], pdg[1] * priors[1], pdg[2] * priors[2]];
        let sum = g.iter().sum::<f64>();
        let mut max = -1.0;
        let mut max_i = 0usize;
        for (i, value) in g.iter_mut().enumerate() {
            *value /= sum;
            if *value > max {
                max = *value;
                max_i = i;
            }
        }
        if !is_var {
            max_i = 2;
            max = g[2];
        }
        let err = (1.0 - max).max(1e-308);
        let q = (-4.343 * err.ln() + 0.499).floor().min(99.0) as u8;
        (q << 2) | max_i as u8
    }

    pub fn drain_afs(&mut self) -> Vec<f64> {
        let mut out = vec![0.0; self.m + 1];
        for (k, item) in out.iter_mut().enumerate().take(self.m + 1) {
            *item = self.afs[self.m - k];
        }
        self.afs.fill(0.0);
        out
    }

    fn indel_prior(&mut self, x: f64) {
        for i in 0..self.m {
            self.phi_indel[i] = self.phi[i] * x;
        }
        self.phi_indel[self.m] = 1.0 - self.phi[self.m] * x;
    }

    fn cal_afs(&mut self, pdg: &[[f64; 3]], is_indel: bool) -> Option<f64> {
        self.afs1.fill(0.0);
        self.cal_y(pdg);
        let phi = if is_indel { &self.phi_indel } else { &self.phi };
        let sum = (0..=self.m).map(|k| phi[k] * self.z[k]).sum::<f64>();
        for (k, afs1) in self.afs1.iter_mut().enumerate().take(self.m + 1) {
            *afs1 = phi[k] * self.z[k] / sum;
            if !afs1.is_finite() {
                return None;
            }
        }

        let folded_sum = (0..=self.m)
            .map(|k| (phi[k] + phi[self.m - k]) * 0.5 * self.z[k])
            .sum::<f64>();
        self.p_var_folded = (1..self.m)
            .map(|k| (phi[k] + phi[self.m - k]) * 0.5 * self.z[k])
            .sum::<f64>()
            / folded_sum;
        self.p_ref_folded =
            (phi[self.m] + phi[0]) * 0.5 * (self.z[self.m] + self.z[0]) / folded_sum;

        let mut expected = 0.0;
        for k in 0..=self.m {
            self.afs[k] += self.afs1[k];
            expected += k as f64 * self.afs1[k];
        }
        Some(expected / self.m as f64)
    }

    fn cal_y(&mut self, pdg: &[[f64; 3]]) {
        if let Some(n1) = self
            .n1
            .filter(|&n1| n1 > 0 && n1 < self.n && self.m == self.n * 2)
        {
            self.z1[..=2 * n1].fill(0.0);
            self.z2[..=2 * (self.n - n1)].fill(0.0);
            self.t1 = 0.0;
            self.t2 = 0.0;
            self.cal_y_core(pdg, n1);
            self.t2 = self.t;
            self.z2[..=2 * (self.n - n1)].copy_from_slice(&self.z[..=2 * (self.n - n1)]);
            self.cal_y_core(pdg, 0);
            let scale = (self.t - (self.t1 + self.t2)).exp();
            for z in &mut self.z {
                *z *= scale;
            }
        } else {
            self.cal_y_core(pdg, 0);
        }
    }

    fn cal_y_core(&mut self, pdg: &[[f64; 3]], beg: usize) {
        self.z.fill(0.0);
        self.zswap.fill(0.0);
        self.z[0] = 1.0;
        let mut last_min = 0usize;
        let mut last_max = 0usize;
        self.t = 0.0;

        if self.m == self.n * 2 {
            let mut m_seen = 0usize;
            for (sample, sample_pdg) in pdg.iter().enumerate().take(self.n).skip(beg) {
                let j = sample - beg;
                let mut min = last_min;
                let mut max = last_max;
                let m0 = m_seen;
                m_seen += 2;
                let p = [sample_pdg[0], 2.0 * sample_pdg[1], sample_pdg[2]];
                while min < max && self.z[min] < MC_TINY {
                    self.z[min] = 0.0;
                    self.zswap[min] = 0.0;
                    min += 1;
                }
                while max > min && self.z[max] < MC_TINY {
                    self.z[max] = 0.0;
                    self.zswap[max] = 0.0;
                    max -= 1;
                }
                max += 2;
                self.zswap[min..=max].fill(0.0);
                if min == 0 {
                    self.zswap[0] = (m0 + 1) as f64 * (m0 + 2) as f64 * p[0] * self.z[0];
                }
                if min <= 1 {
                    self.zswap[1] = m0 as f64 * (m0 + 1) as f64 * p[0] * self.z[1]
                        + p[1] * self.z[0] * (m0 + 1) as f64;
                }
                for k in min.max(2)..=max {
                    self.zswap[k] =
                        signed_offset(m0, k, 1) * signed_offset(m0, k, 2) * p[0] * self.z[k]
                            + k as f64 * signed_offset(m0, k, 2) * p[1] * self.z[k - 1]
                            + k as f64 * (k - 1) as f64 * p[2] * self.z[k - 2];
                }
                let sum = self.zswap[min..=max].iter().sum::<f64>();
                self.t += (sum / (m_seen as f64 * (m_seen - 1) as f64)).ln();
                for k in min..=max {
                    self.zswap[k] /= sum;
                }
                if min >= 1 {
                    self.zswap[min - 1] = 0.0;
                }
                if min >= 2 {
                    self.zswap[min - 2] = 0.0;
                }
                if j < self.n - 1 {
                    if max < self.m {
                        self.zswap[max + 1] = 0.0;
                    }
                    if max + 2 <= self.m {
                        self.zswap[max + 2] = 0.0;
                    }
                }
                if Some(sample + 1) == self.n1 {
                    self.t1 = self.t;
                    self.z1[..=2 * (sample + 1)].copy_from_slice(&self.zswap[..=2 * (sample + 1)]);
                }
                std::mem::swap(&mut self.z, &mut self.zswap);
                last_min = min;
                last_max = max;
            }
        } else {
            let ploidy = self.ploidy.as_ref().expect("mixed ploidy requires ploidy");
            let mut m_seen = 0usize;
            for sample in 0..self.n {
                let mut min = last_min;
                let mut max = last_max;
                while min < max && self.z[min] < MC_TINY {
                    self.z[min] = 0.0;
                    self.zswap[min] = 0.0;
                    min += 1;
                }
                while max > min && self.z[max] < MC_TINY {
                    self.z[max] = 0.0;
                    self.zswap[max] = 0.0;
                    max -= 1;
                }
                let m0 = m_seen;
                m_seen += usize::from(ploidy[sample]);
                if ploidy[sample] == 1 {
                    let p = [pdg[sample][0], pdg[sample][2]];
                    max += 1;
                    self.zswap[min..=max].fill(0.0);
                    if min == 0 {
                        self.zswap[0] = (m0 + 1) as f64 * p[0] * self.z[0];
                    }
                    for k in min.max(1)..=max {
                        self.zswap[k] = signed_offset(m0, k, 1) * p[0] * self.z[k]
                            + k as f64 * p[1] * self.z[k - 1];
                    }
                    let sum = self.zswap[min..=max].iter().sum::<f64>();
                    self.t += (sum / m_seen as f64).ln();
                    for k in min..=max {
                        self.zswap[k] /= sum;
                    }
                    if min >= 1 {
                        self.zswap[min - 1] = 0.0;
                    }
                    if sample < self.n - 1 && max < self.m {
                        self.zswap[max + 1] = 0.0;
                    }
                } else {
                    let p = [pdg[sample][0], 2.0 * pdg[sample][1], pdg[sample][2]];
                    max += 2;
                    self.zswap[min..=max].fill(0.0);
                    if min == 0 {
                        self.zswap[0] = (m0 + 1) as f64 * (m0 + 2) as f64 * p[0] * self.z[0];
                    }
                    if min <= 1 {
                        self.zswap[1] = m0 as f64 * (m0 + 1) as f64 * p[0] * self.z[1]
                            + p[1] * self.z[0] * (m0 + 1) as f64;
                    }
                    for k in min.max(2)..=max {
                        self.zswap[k] =
                            signed_offset(m0, k, 1) * signed_offset(m0, k, 2) * p[0] * self.z[k]
                                + k as f64 * signed_offset(m0, k, 2) * p[1] * self.z[k - 1]
                                + k as f64 * (k - 1) as f64 * p[2] * self.z[k - 2];
                    }
                    let sum = self.zswap[min..=max].iter().sum::<f64>();
                    self.t += (sum / (m_seen as f64 * (m_seen - 1) as f64)).ln();
                    for k in min..=max {
                        self.zswap[k] /= sum;
                    }
                    if min >= 1 {
                        self.zswap[min - 1] = 0.0;
                    }
                    if min >= 2 {
                        self.zswap[min - 2] = 0.0;
                    }
                    if sample < self.n - 1 {
                        if max < self.m {
                            self.zswap[max + 1] = 0.0;
                        }
                        if max + 2 <= self.m {
                            self.zswap[max + 2] = 0.0;
                        }
                    }
                }
                std::mem::swap(&mut self.z, &mut self.zswap);
                last_min = min;
                last_max = max;
            }
        }
    }

    fn contrast2(&mut self, ret: &mut [f64; 3]) -> f64 {
        let Some(n1) = self.n1 else {
            return 0.0;
        };
        let n2 = self.n - n1;
        if n1 == 0 || n2 == 0 {
            return 0.0;
        }
        if self.hypergeom.is_none() {
            let tmp = ln_gamma((2 * (n1 + n2) + 1) as f64)
                - (ln_gamma((2 * n1 + 1) as f64) + ln_gamma((2 * n2 + 1) as f64));
            let mut hg = vec![vec![0.0; 2 * n2 + 1]; 2 * n1 + 1];
            for (k1, row) in hg.iter_mut().enumerate().take(2 * n1 + 1) {
                for (k2, item) in row.iter_mut().enumerate().take(2 * n2 + 1) {
                    *item = (ln_gamma((k1 + k2 + 1) as f64)
                        + ln_gamma((self.m - k1 - k2 + 1) as f64)
                        - (ln_gamma((k1 + 1) as f64)
                            + ln_gamma((k2 + 1) as f64)
                            + ln_gamma((2 * n1 - k1 + 1) as f64)
                            + ln_gamma((2 * n2 - k2 + 1) as f64)
                            + tmp))
                        .exp();
                }
            }
            self.hypergeom = Some(hg);
        }
        let sum = (0..=self.m).map(|k| self.phi[k] * self.z[k]).sum::<f64>();
        let k10 = (0..=2 * n1)
            .max_by(|&a, &b| (self.phi1[a] * self.z1[a]).total_cmp(&(self.phi1[b] * self.z1[b])))
            .unwrap_or(0);
        let k20 = (0..=2 * n2)
            .max_by(|&a, &b| (self.phi2[a] * self.z2[a]).total_cmp(&(self.phi2[b] * self.z2[b])))
            .unwrap_or(0);

        let mut z = 0.0;
        *ret = [0.0; 3];
        for k1 in (0..=k10).rev() {
            for k2 in (0..=k20).rev() {
                let y = self.contrast2_aux(sum, k1, k2, ret);
                if y < 0.0 {
                    break;
                }
                z += y;
            }
            for k2 in k20 + 1..=2 * n2 {
                let y = self.contrast2_aux(sum, k1, k2, ret);
                if y < 0.0 {
                    break;
                }
                z += y;
            }
        }
        let mut right = [0.0; 3];
        for k1 in k10 + 1..=2 * n1 {
            for k2 in (0..=k20).rev() {
                let y = self.contrast2_aux(sum, k1, k2, &mut right);
                if y < 0.0 {
                    break;
                }
                z += y;
            }
            for k2 in k20 + 1..=2 * n2 {
                let y = self.contrast2_aux(sum, k1, k2, &mut right);
                if y < 0.0 {
                    break;
                }
                z += y;
            }
        }
        for i in 0..3 {
            ret[i] += right[i];
        }
        if ret.iter().sum::<f64>() < 0.95 {
            *ret = [0.0; 3];
            z = 0.0;
            for k1 in 0..=2 * n1 {
                for k2 in 0..=2 * n2 {
                    let y = self.contrast2_aux(sum, k1, k2, ret);
                    if y >= 0.0 {
                        z += y;
                    }
                }
            }
            if ret.iter().sum::<f64>() < 0.95 {
                *ret = [1.0 / 3.0; 3];
                z = 1.0;
            }
        }
        z
    }

    fn contrast2_aux(&self, sum: f64, k1: usize, k2: usize, x: &mut [f64; 3]) -> f64 {
        let hg = self.hypergeom.as_ref().expect("hypergeom initialized");
        let p = self.phi[k1 + k2] * self.z1[k1] * self.z2[k2] / sum * hg[k1][k2];
        let n1 = self.n1.expect("n1 set");
        let n2 = self.n - n1;
        if p < CONTRAST_TINY {
            return -1.0;
        }
        let f1 = 0.5 * k1 as f64 / n1 as f64;
        let f2 = 0.5 * k2 as f64 / n2 as f64;
        if f1 < f2 {
            x[1] += p;
        } else if f1 > f2 {
            x[2] += p;
        } else {
            x[0] += p;
        }
        p * chi2_test(k1, k2, 2 * n1 - k1, 2 * n2 - k2)
    }
}

fn signed_offset(m0: usize, k: usize, add: isize) -> f64 {
    (m0 as isize + add - k as isize) as f64
}

fn init_prior(prior_type: PriorType, theta: f64, m: usize, phi: &mut [f64]) {
    match prior_type {
        PriorType::Cond2 => {
            for (i, item) in phi.iter_mut().enumerate().take(m + 1) {
                *item = 2.0 * (i + 1) as f64 / (m + 1) as f64 / (m + 2) as f64;
            }
        }
        PriorType::Flat => {
            for item in phi.iter_mut().take(m + 1) {
                *item = 1.0 / (m + 1) as f64;
            }
        }
        PriorType::Full => {
            let mut sum = 0.0;
            for (i, item) in phi.iter_mut().enumerate().take(m) {
                *item = theta / (m - i) as f64;
                sum += *item;
            }
            phi[m] = 1.0 - sum;
        }
    }
}

fn chi2_test(a: usize, b: usize, c: usize, d: usize) -> f64 {
    let x = (a + b) as f64 * (c + d) as f64 * (b + d) as f64 * (a + c) as f64;
    if x == 0.0 {
        return 1.0;
    }
    let z = (a * d) as f64 - (b * c) as f64;
    gammaq(0.5, 0.5 * z * z * (a + b + c + d) as f64 / x)
}

pub fn regularized_gamma_p(s: f64, z: f64) -> f64 {
    if z <= 1.0 || z < s {
        gammap_series(s, z)
    } else {
        1.0 - gammaq_continued_fraction(s, z)
    }
}

pub fn regularized_gamma_q(s: f64, z: f64) -> f64 {
    if z <= 0.0 {
        return 1.0;
    }
    if z <= 1.0 || z < s {
        1.0 - gammap_series(s, z)
    } else {
        gammaq_continued_fraction(s, z)
    }
}

fn gammaq(s: f64, z: f64) -> f64 {
    regularized_gamma_q(s, z)
}

fn gammap_series(s: f64, z: f64) -> f64 {
    const EPS: f64 = 1e-14;
    let mut sum = 1.0;
    let mut x = 1.0;
    for k in 1..100 {
        x *= z / (s + k as f64);
        sum += x;
        if x / sum < EPS {
            break;
        }
    }
    (s * z.ln() - z - ln_gamma(s + 1.0) + sum.ln()).exp()
}

fn gammaq_continued_fraction(s: f64, z: f64) -> f64 {
    const EPS: f64 = 1e-14;
    const TINY: f64 = 1e-290;
    let mut f = 1.0 + z - s;
    let mut c = f;
    let mut d = 0.0;
    for j in 1..100 {
        let a = j as f64 * (s - j as f64);
        let b = (2 * j + 1) as f64 + z - s;
        d = (b + a * d).max(TINY);
        c = (b + a / c).max(TINY);
        d = 1.0 / d;
        let delta = c * d;
        f *= delta;
        if (delta - 1.0).abs() < EPS {
            break;
        }
    }
    (s * z.ln() - z - ln_gamma(s) - f.ln()).exp()
}

fn ln_gamma(z: f64) -> f64 {
    const COF: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if z < 0.5 {
        return std::f64::consts::PI.ln()
            - (std::f64::consts::PI * z).sin().ln()
            - ln_gamma(1.0 - z);
    }
    let z = z - 1.0;
    let mut x = COF[0];
    for (i, cof) in COF.iter().enumerate().skip(1) {
        x += cof / (z + i as f64);
    }
    let t = z + 7.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (z + 0.5) * t.ln() - t + x.ln()
}

#[derive(Debug, Clone)]
struct HNode {
    akid: Option<usize>,
    bkid: Option<usize>,
    parent: Option<usize>,
    id: usize,
    idx: usize,
    value: f32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct HCluster {
    pub dist: f32,
    pub members: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct HierarchicalClustering {
    ndat: usize,
    pdist: Vec<f32>,
    nodes: Vec<HNode>,
    root: usize,
    explain: String,
}

impl HierarchicalClustering {
    pub fn new(n: usize, pdist: Vec<f32>) -> Self {
        assert_eq!(pdist.len(), n.saturating_mul(n.saturating_sub(1)) / 2);
        let mut this = Self {
            ndat: n,
            pdist,
            nodes: Vec::with_capacity(n.saturating_mul(2)),
            root: 0,
            explain: String::new(),
        };
        let mut active = Vec::with_capacity(n);
        for idx in 0..n {
            active.push(this.append_leaf(idx));
        }
        while active.len() > 1 {
            let (pos_i, pos_j, min_value) = this.closest_pair(&active);
            let max_pos = pos_i.max(pos_j);
            let min_pos = pos_i.min(pos_j);
            let iclust_id = active.remove(max_pos);
            let jclust_id = active.remove(min_pos);
            let iidx = this.nodes[iclust_id].idx;
            let jidx = this.nodes[jclust_id].idx;

            for &active_id in &active {
                let active_idx = this.nodes[active_id].idx;
                let keep = pdist_get(&this.pdist, active_idx, iidx);
                let candidate = pdist_get(&this.pdist, active_idx, jidx);
                if keep < candidate {
                    pdist_set(&mut this.pdist, active_idx, iidx, candidate);
                }
            }

            let node_id = this.append_parent(iidx, min_value, iclust_id, jclust_id);
            this.nodes[iclust_id].parent = Some(node_id);
            this.nodes[jclust_id].parent = Some(node_id);
            active.push(node_id);
        }
        if let Some(&root) = active.first() {
            this.root = root;
        }
        this
    }

    pub fn create_list(&mut self, min_inter_dist: f32, max_intra_dist: &mut f32) -> Vec<HCluster> {
        let cutoff = self.set_threshold(min_inter_dist, *max_intra_dist);
        *max_intra_dist = cutoff;
        let mut clusters = Vec::new();
        let mut stack = vec![self.root];

        if self.nodes[self.root].value < cutoff {
            clusters.push(self.collect_cluster(self.root));
            return clusters;
        }

        while let Some(node_id) = stack.pop() {
            let node = &self.nodes[node_id];
            let Some(akid) = node.akid else {
                clusters.push(self.collect_cluster(node_id));
                continue;
            };
            let bkid = node.bkid.expect("internal node has both children");

            if node.value >= cutoff && self.nodes[akid].value < cutoff {
                clusters.push(self.collect_cluster(akid));
            } else {
                stack.push(akid);
            }

            if node.value >= cutoff && self.nodes[bkid].value < cutoff {
                clusters.push(self.collect_cluster(bkid));
            } else {
                stack.push(bkid);
            }
        }

        clusters
    }

    pub fn explain(&self) -> Vec<&str> {
        self.explain.lines().collect()
    }

    pub fn create_dot(&self, labels: &[&str], threshold: f32) -> String {
        let mut out = String::from("digraph myGraph {");
        for node in &self.nodes {
            if node.value != 0.0 {
                out.push_str(&format!("\"{}\" [label=\"{}\"];", node.id, node.value));
            } else {
                out.push_str(&format!(
                    "\"{}\" [label=\"{}\"];",
                    node.id, labels[node.idx]
                ));
            }
        }
        for node in &self.nodes {
            if let Some(akid) = node.akid {
                self.push_dot_edge(&mut out, node.id, akid, threshold);
            }
            if let Some(bkid) = node.bkid {
                self.push_dot_edge(&mut out, node.id, bkid, threshold);
            }
        }
        out.push_str("};");
        out
    }

    fn append_leaf(&mut self, idx: usize) -> usize {
        let id = self.nodes.len();
        self.nodes.push(HNode {
            akid: None,
            bkid: None,
            parent: None,
            id,
            idx,
            value: 0.0,
        });
        id
    }

    fn append_parent(&mut self, idx: usize, value: f32, akid: usize, bkid: usize) -> usize {
        let id = self.nodes.len();
        self.nodes.push(HNode {
            akid: Some(akid),
            bkid: Some(bkid),
            parent: None,
            id,
            idx,
            value,
        });
        id
    }

    fn closest_pair(&self, active: &[usize]) -> (usize, usize, f32) {
        let mut min_value = f32::INFINITY;
        let mut min_i = 1usize;
        let mut min_j = 0usize;
        for i in 1..active.len() {
            for j in 0..i {
                let value = pdist_get(
                    &self.pdist,
                    self.nodes[active[i]].idx,
                    self.nodes[active[j]].idx,
                );
                if value < min_value {
                    min_value = value;
                    min_i = i;
                    min_j = j;
                }
            }
        }
        (min_i, min_j, min_value)
    }

    fn set_threshold(&mut self, min_inter_dist: f32, max_intra_dist: f32) -> f32 {
        let mut internal: Vec<_> = (self.ndat..self.nodes.len()).collect();
        internal.sort_by(|&a, &b| self.nodes[a].value.total_cmp(&self.nodes[b].value));
        self.explain.clear();

        let mut threshold = max_intra_dist.abs();
        let mut min_dev = f32::INFINITY;
        let mut imin = None;
        for i in 0..internal.len() {
            let mut dev = 0.0;
            if i > 0 {
                dev += calc_node_dev(&self.nodes, &internal[..i]);
            }
            if i + 1 < internal.len() {
                dev += calc_node_dev(&self.nodes, &internal[i..]);
            }
            let th = self.nodes[internal[i]].value;
            self.explain.push_str(&format!("DEV\t{th}\t{dev}\n"));
            if min_dev > dev && th >= min_inter_dist {
                min_dev = dev;
                imin = Some(i);
            }
        }

        if max_intra_dist > 0.0 {
            threshold = max_intra_dist;
        } else if let Some(i) = imin {
            threshold = self.nodes[internal[i]].value.min(max_intra_dist.abs());
        }
        let max_dist = internal
            .last()
            .map(|&id| self.nodes[id].value)
            .unwrap_or(0.0);
        self.explain.push_str(&format!("TH\t{threshold}\n"));
        self.explain.push_str(&format!("MAX_DIST\t{max_dist}\n"));
        self.explain
            .push_str(&format!("MIN_INTER\t{min_inter_dist}\n"));
        self.explain
            .push_str(&format!("MAX_INTRA\t{}\n", max_intra_dist.abs()));
        threshold
    }

    fn collect_cluster(&self, node_id: usize) -> HCluster {
        let mut members = Vec::new();
        let mut stack = vec![node_id];
        let dist = self.nodes[node_id].value;
        while let Some(id) = stack.pop() {
            let node = &self.nodes[id];
            match (node.akid, node.bkid) {
                (Some(a), Some(b)) => {
                    stack.push(a);
                    stack.push(b);
                }
                _ => members.push(node.id),
            }
        }
        members.sort_unstable();
        HCluster { dist, members }
    }

    fn push_dot_edge(&self, out: &mut String, parent: usize, child: usize, threshold: f32) {
        if self.nodes[parent].value >= threshold && self.nodes[child].value < threshold {
            out.push_str(&format!(
                "\"{}\" -> \"{}\" [color=\"#D43F3A\" penwidth=3];",
                self.nodes[parent].id, self.nodes[child].id
            ));
        } else {
            out.push_str(&format!(
                "\"{}\" -> \"{}\";",
                self.nodes[parent].id, self.nodes[child].id
            ));
        }
    }
}

fn pdist_idx(a: usize, b: usize) -> usize {
    if a > b {
        a * (a - 1) / 2 + b
    } else {
        b * (b - 1) / 2 + a
    }
}

fn pdist_get(pdist: &[f32], a: usize, b: usize) -> f32 {
    pdist[pdist_idx(a, b)]
}

fn pdist_set(pdist: &mut [f32], a: usize, b: usize, value: f32) {
    let idx = pdist_idx(a, b);
    pdist[idx] = value;
}

fn calc_node_dev(nodes: &[HNode], ids: &[usize]) -> f32 {
    let avg = ids.iter().map(|&id| nodes[id].value).sum::<f32>() / ids.len() as f32;
    let dev = ids
        .iter()
        .map(|&id| {
            let delta = nodes[id].value - avg;
            delta * delta
        })
        .sum::<f32>();
    (dev / ids.len() as f32).sqrt()
}

const PEAK_NPARAMS: usize = 5;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PeakKind {
    Gaussian,
    BoundedGaussian,
    Exp,
}

#[derive(Debug, Clone)]
struct MonteCarloParam {
    scan: bool,
    min: f64,
    max: f64,
    best: f64,
}

impl Default for MonteCarloParam {
    fn default() -> Self {
        Self {
            scan: false,
            min: 0.0,
            max: 0.0,
            best: 0.0,
        }
    }
}

#[derive(Debug, Clone)]
struct Peak {
    kind: PeakKind,
    fit_mask: u32,
    params: [f64; PEAK_NPARAMS],
    ori_params: [f64; PEAK_NPARAMS],
    mc: [MonteCarloParam; PEAK_NPARAMS],
}

#[derive(Debug, Default, Clone)]
pub struct PeakFit {
    peaks: Vec<Peak>,
    nmc_iter: usize,
}

impl PeakFit {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn reset(&mut self) {
        self.peaks.clear();
        self.nmc_iter = 0;
    }

    pub fn add_gaussian(&mut self, scale: f64, center: f64, sigma: f64, fit_mask: u32) {
        let mut params = [0.0; PEAK_NPARAMS];
        params[0] = scale;
        params[1] = center;
        params[2] = sigma;
        self.peaks
            .push(Peak::new(PeakKind::Gaussian, params, fit_mask));
    }

    pub fn add_bounded_gaussian(
        &mut self,
        scale: f64,
        center: f64,
        sigma: f64,
        min: f64,
        max: f64,
        fit_mask: u32,
    ) {
        assert!(min < max);
        let mut params = [0.0; PEAK_NPARAMS];
        params[0] = scale;
        params[1] = center;
        params[2] = sigma;
        params[3] = min;
        params[4] = max;
        let mut peak = Peak::new(PeakKind::BoundedGaussian, params, fit_mask);
        peak.ori_params[1] = peak.convert_set(1, center);
        peak.params = peak.ori_params;
        self.peaks.push(peak);
    }

    pub fn add_exp(&mut self, scale: f64, center: f64, sigma: f64, fit_mask: u32) {
        assert_eq!(fit_mask & (1 << 1), 0);
        let mut params = [0.0; PEAK_NPARAMS];
        params[0] = scale;
        params[1] = center;
        params[2] = sigma;
        self.peaks.push(Peak::new(PeakKind::Exp, params, fit_mask));
    }

    pub fn set_mc(&mut self, xmin: f64, xmax: f64, iparam: usize, niter: usize) {
        let peak = self.peaks.last_mut().expect("set_mc requires a peak");
        peak.mc[iparam].scan = true;
        peak.mc[iparam].min = xmin;
        peak.mc[iparam].max = xmax;
        self.nmc_iter = niter;
    }

    pub fn set_params(&mut self, ipk: usize, params: &[f64]) {
        let peak = &mut self.peaks[ipk];
        for (i, value) in params.iter().copied().enumerate().take(PEAK_NPARAMS) {
            peak.params[i] = peak.convert_set(i, value);
        }
    }

    pub fn get_params(&self, ipk: usize) -> Vec<f64> {
        self.peaks[ipk].converted_params().to_vec()
    }

    pub fn evaluate(&self, xvals: &[f64], yvals: &[f64]) -> f64 {
        assert_eq!(xvals.len(), yvals.len());
        let vals = evaluate_peaks(&self.peaks, xvals);
        vals.iter()
            .zip(yvals)
            .map(|(fit, obs)| (fit - obs).abs())
            .sum()
    }

    pub fn run(&mut self, xvals: &[f64], yvals: &[f64]) -> f64 {
        assert_eq!(xvals.len(), yvals.len());
        let nparams = self.fit_param_count();
        if nparams == 0 {
            return self.evaluate(xvals, yvals);
        }

        let mut rng = Lcg::new(0);
        let mut best_fit = f64::INFINITY;
        let mut best_peaks = self.peaks.clone();

        for _ in 0..=self.nmc_iter {
            let mut start_peaks = self.peaks.clone();
            for peak in &mut start_peaks {
                peak.params = peak.ori_params;
                for i in 0..PEAK_NPARAMS {
                    if peak.mc[i].scan {
                        let raw = rng.next_between(peak.mc[i].min, peak.mc[i].max);
                        peak.params[i] = peak.convert_set(i, raw);
                    }
                }
            }

            let mut params = collect_fit_params(&start_peaks);
            let template = start_peaks.clone();
            let fit = hooke_jeeves(
                |params| {
                    let mut peaks = template.clone();
                    scatter_fit_params(&mut peaks, params);
                    evaluate_peak_residuals(&peaks, xvals, yvals)
                },
                &mut params,
                0.5,
                1e-8,
                5000,
            );
            let mut fitted = template;
            scatter_fit_params(&mut fitted, &params);
            let fit = fit.min(evaluate_peak_residuals(&fitted, xvals, yvals));
            if fit < best_fit {
                best_fit = fit;
                best_peaks = fitted;
            }
        }

        for peak in &mut best_peaks {
            for i in 0..PEAK_NPARAMS {
                peak.mc[i].best = peak.params[i];
            }
        }
        self.peaks = best_peaks;
        best_fit
    }

    pub fn sprint_func(&self) -> String {
        self.peaks
            .iter()
            .map(Peak::sprint_func)
            .collect::<Vec<_>>()
            .join(" + ")
    }

    fn fit_param_count(&self) -> usize {
        self.peaks
            .iter()
            .map(|peak| peak.fit_mask.count_ones() as usize)
            .sum()
    }
}

impl Peak {
    fn new(kind: PeakKind, params: [f64; PEAK_NPARAMS], fit_mask: u32) -> Self {
        Self {
            kind,
            fit_mask,
            params,
            ori_params: params,
            mc: std::array::from_fn(|_| MonteCarloParam::default()),
        }
    }

    fn convert_set(&self, iparam: usize, value: f64) -> f64 {
        if self.kind != PeakKind::BoundedGaussian || iparam != 1 {
            return value;
        }
        let min = self.ori_params[3];
        let max = self.ori_params[4];
        let value = value.clamp(min, max);
        (2.0 * (value - min) / (max - min) - 1.0).acos()
    }

    fn converted_params(&self) -> [f64; PEAK_NPARAMS] {
        let mut out = self.params;
        match self.kind {
            PeakKind::Gaussian | PeakKind::Exp => {
                out[0] = self.params[0].abs();
                out[1] = self.params[1].abs();
                out[2] = self.params[2].abs();
            }
            PeakKind::BoundedGaussian => {
                out[0] = self.params[0].abs();
                out[2] = self.params[2].abs();
                out[1] = self.bounded_center();
            }
        }
        out
    }

    fn bounded_center(&self) -> f64 {
        0.5 * (self.params[1].cos() + 1.0) * (self.params[4] - self.params[3]) + self.params[3]
    }

    fn add_values(&self, xvals: &[f64], yvals: &mut [f64]) {
        match self.kind {
            PeakKind::Gaussian => {
                let scale2 = self.params[0] * self.params[0];
                let center = self.params[1];
                let sigma = self.params[2];
                for (&x, y) in xvals.iter().zip(yvals.iter_mut()) {
                    let tmp = (x - center) / sigma;
                    *y += scale2 * (-tmp * tmp).exp();
                }
            }
            PeakKind::BoundedGaussian => {
                let scale2 = self.params[0] * self.params[0];
                let center = self.bounded_center();
                let sigma = self.params[2];
                for (&x, y) in xvals.iter().zip(yvals.iter_mut()) {
                    let tmp = (x - center) / sigma;
                    *y += scale2 * (-tmp * tmp).exp();
                }
            }
            PeakKind::Exp => {
                let scale2 = self.params[0] * self.params[0];
                let center = self.params[1];
                let sigma = self.params[2];
                for (&x, y) in xvals.iter().zip(yvals.iter_mut()) {
                    *y += scale2 * ((x - center) / sigma / sigma).exp();
                }
            }
        }
    }

    fn sprint_func(&self) -> String {
        match self.kind {
            PeakKind::Gaussian => format!(
                "{}**2 * exp(-(x-{})**2/{}**2)",
                self.params[0].abs(),
                self.params[1],
                self.params[2].abs()
            ),
            PeakKind::BoundedGaussian => format!(
                "{}**2 * exp(-(x-{})**2/{}**2)",
                self.params[0].abs(),
                self.bounded_center(),
                self.params[2].abs()
            ),
            PeakKind::Exp => format!(
                "{}**2 * exp((x-{})/{}**2)",
                self.params[0].abs(),
                self.params[1],
                self.params[2].abs()
            ),
        }
    }
}

fn collect_fit_params(peaks: &[Peak]) -> Vec<f64> {
    let mut out = Vec::new();
    for peak in peaks {
        for i in 0..PEAK_NPARAMS {
            if peak.fit_mask & (1 << i) != 0 {
                out.push(peak.params[i]);
            }
        }
    }
    out
}

fn scatter_fit_params(peaks: &mut [Peak], params: &[f64]) {
    let mut iparam = 0;
    for peak in peaks {
        for i in 0..PEAK_NPARAMS {
            if peak.fit_mask & (1 << i) != 0 {
                peak.params[i] = params[iparam];
                iparam += 1;
            }
        }
    }
}

fn evaluate_peaks(peaks: &[Peak], xvals: &[f64]) -> Vec<f64> {
    let mut vals = vec![0.0; xvals.len()];
    for peak in peaks {
        peak.add_values(xvals, &mut vals);
    }
    vals
}

fn evaluate_peak_residuals(peaks: &[Peak], xvals: &[f64], yvals: &[f64]) -> f64 {
    evaluate_peaks(peaks, xvals)
        .iter()
        .zip(yvals)
        .map(|(fit, obs)| (fit - obs).abs())
        .sum()
}

#[derive(Debug, Clone)]
struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_f64(&mut self) -> f64 {
        self.state = self.state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
        ((self.state / 65_536) % 32_768) as f64 / 32_767.0
    }

    fn next_between(&mut self, min: f64, max: f64) -> f64 {
        self.next_f64() * (max - min) + min
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_distribution_keeps_exact_bins_then_log_bins() {
        let mut dist = LogDistribution::new(2);
        assert_eq!(dist.insert(0), 0);
        assert_eq!(dist.insert(42), 42);
        assert_eq!(dist.insert(100), 100);
        assert_eq!(dist.insert(101), 100);
        assert_eq!(dist.nbins(), 101);
        assert_eq!(dist.nvalues(), 4);
        assert_eq!(
            dist.get(100),
            Some(DistributionBin {
                beg: 100,
                end: 110,
                count: 2
            })
        );
    }

    #[test]
    fn log_distribution_insert_n_matches_upstream_counter_behavior() {
        let mut dist = LogDistribution::new(1);
        let bin = dist.insert_n(12, 3);
        assert_eq!(bin, 10);
        assert_eq!(dist.get(bin).unwrap().count, 3);
        assert_eq!(dist.nvalues(), 4);
    }

    #[test]
    fn hooke_jeeves_minimizes_quadratic_bowl() {
        let mut x = [0.0, 0.0];
        let f = hooke_jeeves(
            |x| (x[0] - 2.0).powi(2) + (x[1] + 3.0).powi(2),
            &mut x,
            0.5,
            1e-8,
            10_000,
        );
        assert!(f < 1e-10, "f={f}, x={x:?}");
        assert!((x[0] - 2.0).abs() < 1e-4, "x={x:?}");
        assert!((x[1] + 3.0).abs() < 1e-4, "x={x:?}");
    }

    #[test]
    fn brent_minimizes_univariate_quadratic() {
        let (xmin, fmin) = brent(|x| (x - 4.0).powi(2) + 2.0, 0.0, 1.0, 1e-10);
        assert!((xmin - 4.0).abs() < 1e-6, "xmin={xmin}");
        assert!((fmin - 2.0).abs() < 1e-9, "fmin={fmin}");
    }

    #[test]
    fn incomplete_gamma_helpers_match_known_values() {
        let p = regularized_gamma_p(1.0, 2.0);
        let q = regularized_gamma_q(1.0, 2.0);
        assert!((p - (1.0 - (-2.0_f64).exp())).abs() < 1e-12, "p={p}");
        assert!((q - (-2.0_f64).exp()).abs() < 1e-12, "q={q}");
        assert!((p + q - 1.0).abs() < 1e-12);
    }

    #[test]
    fn single_locus_em_estimates_ref_frequency_and_group_tests() {
        let pdg = [
            [0.001, 0.01, 1.0],
            [0.001, 1.0, 0.001],
            [1.0, 0.01, 0.001],
            [1.0, 0.01, 0.001],
        ];
        let result = single_locus_em(&pdg, Some(2), true, true).unwrap();
        assert!((result.ref_frequency - 0.375).abs() < 0.05, "{result:?}");
        assert!(result.hwe_p_value.is_finite());
        assert!(result.group_frequencies.unwrap()[0] > result.group_frequencies.unwrap()[1]);
        assert!(result.one_degree_p_value.unwrap().is_finite());
        assert!(result.two_degree_p_value.unwrap().is_finite());
    }

    #[test]
    fn prob1_posterior_favors_reference_for_reference_likelihoods() {
        let pdg = [[0.001, 0.01, 1.0], [0.001, 0.01, 1.0]];
        let mut p1 = Prob1::new(2, None);
        let result = p1.calculate(&pdg, false, false).unwrap();
        assert_eq!(result.ac, 0, "{result:?}");
        assert!(result.p_ref > 0.9, "{result:?}");
        assert!(result.p_var < 0.1, "{result:?}");
        assert!(result.f_exp > 0.9, "{result:?}");
        let gt = p1.call_gt(pdg[0], result.f_exp, 0, false);
        assert_eq!(gt & 3, 2);
    }

    #[test]
    fn prob1_supports_group_split_and_mixed_ploidy() {
        let pdg = [
            [0.001, 0.01, 1.0],
            [0.001, 1.0, 0.001],
            [1.0, 0.01, 0.001],
            [1.0, 0.01, 0.001],
        ];
        let mut split = Prob1::new(4, None);
        split.set_n1(2).unwrap();
        split.init_subprior(PriorType::Full, 1e-3);
        let result = split.calculate(&pdg, false, true).unwrap();
        assert!(result.lrt.is_finite(), "{result:?}");
        assert!(result.p_chi2.is_finite(), "{result:?}");

        let mixed_pdg = [[0.001, 0.01, 1.0], [1.0, 0.01, 0.001]];
        let mut mixed = Prob1::new(2, Some(vec![1, 2]));
        assert!(mixed.set_n1(1).is_err());
        let result = mixed.calculate(&mixed_pdg, false, false).unwrap();
        assert!(result.f_exp.is_finite(), "{result:?}");
        assert_eq!(mixed.call_gt(mixed_pdg[0], result.f_exp, 0, true) & 3, 2);
    }

    #[test]
    fn hierarchical_clustering_uses_complete_linkage() {
        let pdist = vec![
            0.1, // d(1,0)
            0.8, 0.9, // d(2,0), d(2,1)
            0.85, 0.95, 0.2, // d(3,0), d(3,1), d(3,2)
        ];
        let mut clust = HierarchicalClustering::new(4, pdist);
        let mut cutoff = 0.5;
        let mut clusters = clust.create_list(0.0, &mut cutoff);
        clusters.sort_by_key(|cluster| cluster.members[0]);
        assert_eq!(cutoff, 0.5);
        assert_eq!(
            clusters,
            [
                HCluster {
                    dist: 0.1,
                    members: vec![0, 1]
                },
                HCluster {
                    dist: 0.2,
                    members: vec![2, 3]
                }
            ]
        );
    }

    #[test]
    fn hierarchical_clustering_reports_threshold_explanation_and_dot() {
        let pdist = vec![0.1, 0.8, 0.9, 0.85, 0.95, 0.2];
        let mut clust = HierarchicalClustering::new(4, pdist);
        let mut cutoff = -0.6;
        let clusters = clust.create_list(0.0, &mut cutoff);
        assert!(!clusters.is_empty());
        assert!(clust.explain().iter().any(|line| line.starts_with("TH\t")));

        let dot = clust.create_dot(&["a", "b", "c", "d"], cutoff);
        assert!(dot.starts_with("digraph myGraph {"));
        assert!(dot.contains("\"0\" [label=\"a\"]"));
        assert!(dot.ends_with("};"));
    }

    #[test]
    fn peakfit_evaluates_and_formats_peak_sum() {
        let mut fit = PeakFit::new();
        fit.add_gaussian(2.0, 1.0, 0.5, 0);
        fit.add_exp(1.0, 0.0, 2.0, 0);
        let x = [1.0];
        let y = [4.0 + 0.25_f64.exp()];
        assert!(fit.evaluate(&x, &y) < 1e-12);
        let func = fit.sprint_func();
        assert!(func.contains("exp(-(x-1"));
        assert!(func.contains(" + "));
    }

    #[test]
    fn peakfit_bounded_gaussian_converts_center_back_to_requested_bounds() {
        let mut fit = PeakFit::new();
        fit.add_bounded_gaussian(1.0, 0.75, 0.1, 0.0, 1.0, 0);
        let params = fit.get_params(0);
        assert!((params[1] - 0.75).abs() < 1e-12, "{params:?}");
        fit.set_params(0, &[1.0, 1.5, 0.1, 0.0, 1.0]);
        let params = fit.get_params(0);
        assert!((params[1] - 1.0).abs() < 1e-12, "{params:?}");
    }

    #[test]
    fn peakfit_native_solver_recovers_simple_gaussian_parameters() {
        let x: Vec<_> = (-10..=10).map(|i| f64::from(i) / 5.0).collect();
        let y: Vec<_> = x
            .iter()
            .map(|&x| 3.0_f64.powi(2) * (-(x - 0.4).powi(2) / 0.7_f64.powi(2)).exp())
            .collect();
        let mut fit = PeakFit::new();
        fit.add_gaussian(2.5, 0.0, 1.0, 0b111);
        let before = fit.evaluate(&x, &y);
        let after = fit.run(&x, &y);
        assert!(after < before, "before={before}, after={after}");
        let params = fit.get_params(0);
        assert!((params[0] - 3.0).abs() < 0.2, "{params:?}");
        assert!((params[1] - 0.4).abs() < 0.2, "{params:?}");
        assert!((params[2] - 0.7).abs() < 0.2, "{params:?}");
    }
}
