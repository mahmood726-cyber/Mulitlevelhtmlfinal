// tests.js — pure Node tests for engine.js (786-MIII Multilevel & Multivariate).
// Every expected value is hand-derived INDEPENDENTLY (not by running the engine).
// Run: node tests.js   ->  exits 0 if all pass, 1 otherwise.

const { Engine, Matrix, Metric } = require('./engine.js');

let N = 0, M = 0;
function approx(label, got, exp, tol) {
    tol = (tol === undefined) ? 1e-6 : tol;
    if (typeof got === 'number' && Math.abs(got - exp) <= tol) { N++; console.log('PASS ' + label + '  (' + got + ')'); }
    else { M++; console.log('FAIL ' + label + '  got=' + got + ' expected=' + exp + ' tol=' + tol); }
}
function ok(label, cond, info) {
    if (cond) { N++; console.log('PASS ' + label); }
    else { M++; console.log('FAIL ' + label + (info ? '  ' + info : '')); }
}

// ---------------------------------------------------------------------------
// 1. Normal-tail p (two-sided). getP is a p-value, NOT a CDF.
//    At z=0 the two-sided tail probability is 1.0 (NOT 0.5).  [This engine has
//    no normalCDF; the spec's "CDF(0)=0.5" check does not apply — getP is a p.]
//    At z=1.96 the two-sided p is ~0.05.  getP uses |z| so it is symmetric.
// ---------------------------------------------------------------------------
approx('getP(0) two-sided p = 1.0', Engine.getP(0), 1.0, 1e-3);
approx('getP(1.96) ~= 0.05', Engine.getP(1.96), 0.05, 5e-4);
approx('getP(-1.96) == getP(1.96)', Engine.getP(-1.96), Engine.getP(1.96), 1e-12);

// ---------------------------------------------------------------------------
// 2. t-approximation getT(df) = 1.96 + 2.37/df  (df>0); 1.96 for df<=0.
//    Hand: getT(1)=1.96+2.37=4.33 ; getT(10)=1.96+0.237=2.197 ; getT(0)=1.96.
// ---------------------------------------------------------------------------
approx('getT(1) = 4.33', Engine.getT(1), 4.33, 1e-9);
approx('getT(10) = 2.197', Engine.getT(10), 2.197, 1e-9);
approx('getT(0) = 1.96 (floor)', Engine.getT(0), 1.96, 1e-12);

// ---------------------------------------------------------------------------
// 3. Matrix algebra: inverse and product (the engine's IGLS depends on these).
//    A = [[4,3],[6,3]] ; det = 4*3 - 3*6 = -6 ; A^-1 = 1/det*[[3,-3],[-6,4]]
//       = [[-0.5, 0.5],[1, -0.6666667]] ; A^-1 * A must equal I.
// ---------------------------------------------------------------------------
{
    const Ai = Matrix.inv([[4, 3], [6, 3]]);
    approx('Matrix.inv [0][0] = -0.5', Ai[0][0], -0.5, 1e-9);
    approx('Matrix.inv [0][1] =  0.5', Ai[0][1], 0.5, 1e-9);
    approx('Matrix.inv [1][0] =  1.0', Ai[1][0], 1.0, 1e-9);
    approx('Matrix.inv [1][1] = -0.6666667', Ai[1][1], -0.6666666667, 1e-9);
    const II = Matrix.dot(Ai, [[4, 3], [6, 3]]);
    ok('Matrix.inv*A = I (A.Ainv.A=A identity)',
        Math.abs(II[0][0] - 1) < 1e-9 && Math.abs(II[1][1] - 1) < 1e-9 &&
        Math.abs(II[0][1]) < 1e-9 && Math.abs(II[1][0]) < 1e-9);
    // transpose
    const T = Matrix.transpose([[1, 2, 3], [4, 5, 6]]);
    ok('Matrix.transpose shape 3x2', T.length === 3 && T[0].length === 2 && T[2][0] === 3 && T[2][1] === 6);
}

// ---------------------------------------------------------------------------
// 4. FULLY HAND-WORKED FIXED-EFFECT POOLING (Continuous MD, 2 effects, 1 cluster).
//    Fixed => tau_b^2 = tau_w^2 = 0, so Sigma is diagonal( vi ) and the IGLS
//    reduces to the inverse-variance weighted mean.
//    Effect 1: m1=10,s1=2,n1=10 ; m2=5,s2=2,n2=10
//        MD_1 = 10-5 = 5 ; vi_1 = 4/10 + 4/10 = 0.8 ; w_1 = 1/0.8 = 1.25
//    Effect 2: m1=12,s1=3,n1=20 ; m2=6,s2=3,n2=20
//        MD_2 = 12-6 = 6 ; vi_2 = 9/20 + 9/20 = 0.9 ; w_2 = 1/0.9 = 1.111111...
//    XtWX = sum(w) = 2.3611111 ; CovBeta = 1/2.3611111 = 0.4235294
//    beta = (1.25*5 + 1.111111*6)/2.3611111 = 12.9166667/2.3611111 = 5.4705882
//    pse(model) = sqrt(0.4235294) = 0.6507914
//    Only 1 cluster => tVal = 1.96 (the m>1 t-correction does not engage).
//    displayVal = beta (MD linear) = 5.4705882
//    lo = 5.4705882 - 1.96*0.6507914 = 4.1950371
//    hi = 5.4705882 + 1.96*0.6507914 = 6.7461393
//    Q (fixed branch never updates Q) = 0 ; tau2 = 0 ; I2 = 0.
//    k_studies (clusters) = 1 ; k_effects = 2.
// ---------------------------------------------------------------------------
{
    const eff = Engine.calcEffects(
        [{ id: 'A', m1: 10, s1: 2, n1: 10, m2: 5, s2: 2, n2: 10 },
         { id: 'A', m1: 12, s1: 3, n1: 20, m2: 6, s2: 3, n2: 20 }], 'cont', 'MD');
    approx('MD effect1 es', eff[0].es, 5, 1e-12);
    approx('MD effect1 vi', eff[0].vi, 0.8, 1e-12);
    approx('MD effect2 es', eff[1].es, 6, 1e-12);
    approx('MD effect2 vi', eff[1].vi, 0.9, 1e-12);
    const r = Engine.pool(eff, 'fixed', null, null, null);
    approx('IGLS fixed pooled es', r.es, 5.4705882353, 1e-9);
    approx('IGLS fixed model se', r.se, 0.6507913735, 1e-9);
    approx('IGLS fixed displayVal == es (MD linear)', r.displayVal, r.es, 1e-12);
    approx('IGLS fixed lo', r.lo, 4.1950371433, 1e-7);
    approx('IGLS fixed hi', r.hi, 6.7461393273, 1e-7);
    approx('IGLS fixed tau2 = 0', r.tau2, 0, 1e-12);
    approx('IGLS fixed I2 = 0', r.I2, 0, 1e-12);
    ok('k_studies (clusters) = 1', r.k_studies === 1);
    ok('k_effects = 2', r.k_effects === 2);
}

// ---------------------------------------------------------------------------
// 5. POOL ON LOG SCALE then back-transform (OR). Two identical effects in two
//    separate clusters => no heterogeneity => tau_b^2 = tau_w^2 = 0, I2 = 0.
//    OR: e1=10,n1=100,e2=20,n2=100 (no zero cells => cc=0).
//    OR  = (10*80)/(20*90) = 800/1800 = 0.4444444 ; logOR = ln(0.4444444) = -0.8109302
//    vi  = 1/10 + 1/90 + 1/20 + 1/80 = 0.1 + 0.0111111 + 0.05 + 0.0125 = 0.1736111
//    Two identical, independent clusters => pooled logOR = -0.8109302, disp = 0.4444444.
// ---------------------------------------------------------------------------
{
    const eff = Engine.calcEffects(
        [{ id: '1', e1: 10, n1: 100, e2: 20, n2: 100 },
         { id: '2', e1: 10, n1: 100, e2: 20, n2: 100 }], 'binary', 'OR');
    approx('OR each es (logOR)', eff[0].es, -0.8109302162, 1e-9);
    approx('OR each vi', eff[0].vi, 0.1736111111, 1e-9);
    const r = Engine.pool(eff, 'random', null, null, null);
    approx('OR identical tau2_b = 0', r.tau2_b, 0, 1e-9);
    approx('OR identical tau2_w = 0', r.tau2_w, 0, 1e-9);
    approx('OR identical I2 = 0', r.I2, 0, 1e-9);
    approx('OR pooled logOR', r.es, -0.8109302162, 1e-9);
    approx('OR pooled displayVal (back-transform)', r.displayVal, 0.4444444444, 1e-9);
    ok('OR pooled on LOG scale then exp()-back-transformed',
        Math.abs(Math.exp(r.es) - r.displayVal) < 1e-12);
    ok('OR pooled Q ~ 0 (identical)', Math.abs(r.Q) < 1e-9, 'Q=' + r.Q);
}

// ---------------------------------------------------------------------------
// 6. TRUE MULTILEVEL STRUCTURE: variance decomposition is non-negative and
//    sums coherently. 3 studies / 4 effects (study S1 contributes 2 effects).
//    Cheung's decomposition: I2_total = (tau_b^2 + tau_w^2)/(tau_b^2+tau_w^2+typVi)*100,
//    with I2_3 (between) and I2_2 (within) summing to I2_total. All in [0,100].
// ---------------------------------------------------------------------------
{
    const eff = Engine.calcEffects(
        [{ id: 'S1', m1: 10, s1: 2, n1: 30, m2: 5, s2: 2, n2: 30 },
         { id: 'S1', m1: 11, s1: 2, n1: 30, m2: 5, s2: 2, n2: 30 },
         { id: 'S2', m1: 8,  s1: 2, n1: 30, m2: 7, s2: 2, n2: 30 },
         { id: 'S3', m1: 9,  s1: 2, n1: 30, m2: 6, s2: 2, n2: 30 }], 'cont', 'MD');
    const r = Engine.pool(eff, 'random', null, null, null);
    ok('multilevel k_studies = 3 clusters', r.k_studies === 3);
    ok('multilevel k_effects = 4', r.k_effects === 4);
    ok('tau2_b >= 0', r.tau2_b >= 0, 'tau2_b=' + r.tau2_b);
    ok('tau2_w >= 0', r.tau2_w >= 0, 'tau2_w=' + r.tau2_w);
    ok('tau2 == tau2_b + tau2_w', Math.abs(r.tau2 - (r.tau2_b + r.tau2_w)) < 1e-12);
    ok('I2_total in [0,100]', r.I2 >= 0 && r.I2 <= 100, 'I2=' + r.I2);
    ok('I2 components non-negative', r.I2_2 >= 0 && r.I2_3 >= 0);
    ok('I2_3 + I2_2 == I2_total (Cheung decomposition)',
        Math.abs((r.I2_3 + r.I2_2) - r.I2) < 1e-9, 'sum=' + (r.I2_3 + r.I2_2) + ' I2=' + r.I2);
    ok('pVal in [0,1]', r.pVal >= 0 && r.pVal <= 1, 'pVal=' + r.pVal);
    ok('model se finite & > 0', isFinite(r.se) && r.se > 0);
    // RVE sandwich: with >= 2 clusters, robust SE is finite and > 0.
    ok('RVE se_r finite & > 0 (m>=2 clusters)', isFinite(r.se_r) && r.se_r > 0, 'se_r=' + r.se_r);
}

// ---------------------------------------------------------------------------
// 7. EDGE: k=1 single effect passthrough (1 cluster, 1 effect, fixed).
//    es = 5 ; pooled es passes through = 5 ; displayVal = 5.
//    RVE small-sample correction sqrt(m/(m-1)) with m=1 => sqrt(1/0) = Inf,
//    so se_r is NaN (RVE is undefined for a single cluster) — the MODEL-based
//    estimate is still finite. Documented caveat; test asserts that contract.
// ---------------------------------------------------------------------------
{
    const eff = Engine.calcEffects(
        [{ id: 'X', m1: 10, s1: 2, n1: 10, m2: 5, s2: 2, n2: 10 }], 'cont', 'MD');
    ok('k=1 produces one effect', eff.length === 1);
    const r = Engine.pool(eff, 'fixed', null, null, null);
    approx('k=1 pooled es passthrough', r.es, 5, 1e-12);
    approx('k=1 displayVal passthrough', r.displayVal, 5, 1e-12);
    ok('k=1 model se finite', isFinite(r.se) && r.se > 0);
    ok('k=1 RVE se_r is NaN (undefined for single cluster)', Number.isNaN(r.se_r));
}

// ---------------------------------------------------------------------------
// 8. EDGE: empty guards.
//    calcEffects([]) -> [] ; pool([]) -> null ; pool(all-excluded) -> null.
// ---------------------------------------------------------------------------
ok('calcEffects([]) returns empty array',
    Array.isArray(Engine.calcEffects([], 'binary', 'OR')) && Engine.calcEffects([], 'binary', 'OR').length === 0);
ok('pool([]) returns null', Engine.pool([], 'fixed', null, null, null) === null);
ok('pool(all-excluded) returns null',
    Engine.pool([{ es: 1, vi: 1, excluded: true, id: 'a', metric: 'MD' }], 'fixed', null, null, null) === null);

// ---------------------------------------------------------------------------
// 9. SMD (continuous) uses Hedges' g small-sample correction.
//    A: m1=10,s1=2,n1=20 ; m2=8,s2=2,n2=20. df=38.
//    sp = sqrt(((19*4)+(19*4))/38) = sqrt(152/38) = sqrt(4) = 2.
//    d  = (10-8)/2 = 1.0 ; J = 1 - 3/(4*38-1) = 1 - 3/151 = 0.9801325.
//    g  = 1.0 * 0.9801325 = 0.9801325.
//    vi = (n1+n2)/(n1*n2) + g^2/(2*(n1+n2)) = 40/400 + 0.9606597/80 = 0.1 + 0.0120082 = 0.1120082.
// ---------------------------------------------------------------------------
{
    const eff = Engine.calcEffects(
        [{ id: 'A', m1: 10, s1: 2, n1: 20, m2: 8, s2: 2, n2: 20 }], 'cont', 'SMD');
    approx('SMD Hedges g', eff[0].es, 0.9801324503, 1e-9);
    approx('SMD vi', eff[0].vi, 0.1120082, 1e-6);
}

// ---------------------------------------------------------------------------
// 10. Survival HR: log scale, se from CI width / 3.92.
//     hr=0.8, lci=0.6, uci=0.95 => es=ln(0.8)=-0.2231436.
//     se = (ln(0.95)-ln(0.6))/3.92 = (-0.0512933+0.5108256)/3.92 = 0.4595324/3.92 = 0.1172276.
//     vi = se^2 = 0.0137423.
// ---------------------------------------------------------------------------
{
    const eff = Engine.calcEffects([{ id: '1', hr: 0.8, lci: 0.6, uci: 0.95, n1: 100 }], 'survival', 'HR');
    approx('HR es = ln(0.8)', eff[0].es, -0.2231435513, 1e-9);
    approx('HR vi', eff[0].vi, 0.0137423, 1e-6);
}

// ---------------------------------------------------------------------------
// 11. Proportion (logit). e=30,n=100 (no CC). p=0.3.
//     es = ln(0.3/0.7) = ln(0.4285714) = -0.8472979.
//     vi = 1/(n*p*(1-p)) = 1/(100*0.3*0.7) = 1/21 = 0.0476190.
//     displayVal = e/n = 0.3.
// ---------------------------------------------------------------------------
{
    const eff = Engine.calcEffects([{ id: '1', e: 30, n: 100 }], 'prop', 'Prop');
    approx('Prop logit es', eff[0].es, -0.8472978604, 1e-9);
    approx('Prop vi', eff[0].vi, 0.0476190476, 1e-9);
    approx('Prop displayVal', eff[0].displayVal, 0.3, 1e-12);
}

// ---------------------------------------------------------------------------
// 12. Multilevel Egger's test: guard (<3 effects) and output sanity.
//     k<3 => {p:null, intercept:0}. For k>=3 (asymmetric), p must be in [0,1]
//     and the intercept finite.
// ---------------------------------------------------------------------------
{
    const two = Engine.calcEffects(
        [{ id: '1', e1: 10, n1: 100, e2: 20, n2: 100 },
         { id: '2', e1: 12, n1: 100, e2: 22, n2: 100 }], 'binary', 'OR');
    const g2 = Engine.multilevelEgger(two);
    ok('Egger k<3 returns null p', g2.p === null && g2.intercept === 0);

    const many = Engine.calcEffects(
        [{ id: '1', e1: 5,  n1: 100, e2: 20, n2: 100 },
         { id: '2', e1: 10, n1: 120, e2: 22, n2: 120 },
         { id: '3', e1: 30, n1: 200, e2: 35, n2: 200 },
         { id: '4', e1: 40, n1: 300, e2: 44, n2: 300 }], 'binary', 'OR');
    const g4 = Engine.multilevelEgger(many);
    ok('Egger k>=3 p in [0,1]', g4.p !== null && g4.p >= 0 && g4.p <= 1, 'p=' + g4.p);
    ok('Egger k>=3 intercept finite', isFinite(g4.intercept));
}

// ---------------------------------------------------------------------------
// 13. Binary continuity correction fires ONLY when a cell is zero/full (cc=0.5),
//     not unconditionally. RR with e1=0 must add 0.5; with all cells interior, cc=0.
//     Interior: e1=10,n1=100,e2=20,n2=100 => RR=(10/100)/(20/100)=0.5 => logRR=ln(0.5).
// ---------------------------------------------------------------------------
{
    const interior = Engine.calcEffects([{ id: '1', e1: 10, n1: 100, e2: 20, n2: 100 }], 'binary', 'RR');
    approx('RR interior logRR = ln(0.5) (no cc)', interior[0].es, Math.log(0.5), 1e-12);
    const zero = Engine.calcEffects([{ id: '1', e1: 0, n1: 100, e2: 20, n2: 100 }], 'binary', 'RR');
    // cc=0.5: RR = ((0+0.5)/(100+0.5)) / ((20+0.5)/(100+0.5)) = 0.5/20.5 = 0.0243902
    approx('RR zero-cell logRR with cc=0.5', zero[0].es, Math.log((0.5 / 100.5) / (20.5 / 100.5)), 1e-12);
    ok('RR zero-cell vi finite (cc applied)', isFinite(zero[0].vi) && zero[0].vi > 0);
}

console.log('\n' + N + ' passed, ' + M + ' failed');
process.exit(M === 0 ? 0 : 1);
