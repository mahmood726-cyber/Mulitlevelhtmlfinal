// engine.js — pure statistical core for 786-MIII Multilevel & Multivariate Synthesis.
// Extracted VERBATIM from the inline app script (single source of truth).
// Pure functions/objects only — no DOM / Plotly dependencies:
//   Matrix              : dense matrix algebra (zeros/identity/transpose/dot/inv/scalarMult/add/sub)
//   Metric.isLogRatio   : which metrics are pooled on the log scale
//   Engine.toNum/validNum/getT/getP : numeric helpers + t-approx + normal-tail two-sided p
//   Engine.calcEffects  : per-row effect-size + sampling variance (log scale for ratios/logit prop)
//   Engine.pool         : multilevel IGLS solver (block-diagonal Sigma), RVE sandwich, Cheung's I2
//   Engine.multilevelEgger : multilevel Egger's small-study test (precision-weighted OLS)

        // --- MATRIX ALGEBRA ENGINE ---
        const Matrix = {
            zeros: (r, c) => Array(r).fill(0).map(() => Array(c).fill(0)),
            identity: (n) => Array(n).fill(0).map((_, i) => Array(n).fill(0).map((__, j) => i === j ? 1 : 0)),
            transpose: (m) => m[0].map((_, i) => m.map(r => r[i])),
            dot: (a, b) => {
                const r1 = a.length, c1 = a[0].length, c2 = b[0].length;
                const res = Matrix.zeros(r1, c2);
                for (let i = 0; i < r1; i++) {
                    for (let j = 0; j < c2; j++) {
                        let sum = 0;
                        for (let k = 0; k < c1; k++) sum += a[i][k] * b[k][j];
                        res[i][j] = sum;
                    }
                }
                return res;
            },
            inv: (m) => {
                const n = m.length;
                const A = m.map(row => [...row]);
                const I = Matrix.identity(n);
                for (let i = 0; i < n; i++) {
                    let pivot = A[i][i];
                    for (let j = 0; j < n; j++) { A[i][j] /= pivot; I[i][j] /= pivot; }
                    for (let k = 0; k < n; k++) {
                        if (k !== i) {
                            let f = A[k][i];
                            for (let j = 0; j < n; j++) { A[k][j] -= f * A[i][j]; I[k][j] -= f * I[i][j]; }
                        }
                    }
                }
                return I;
            },
            scalarMult: (m, s) => m.map(r => r.map(x => x * s)),
            add: (a, b) => a.map((r, i) => r.map((x, j) => x + b[i][j])),
            sub: (a, b) => a.map((r, i) => r.map((x, j) => x - b[i][j]))
        };

        const Metric={isLogRatio:m=>['RR','OR','HR'].includes(m)};

        const Engine={
            toNum:v=>{const x=parseFloat(v);return Number.isFinite(x)?x:0;},
            validNum:v=>Number.isFinite(v)&&!isNaN(v),
            getT:df=>(df<=0)?1.96:1.96+2.37/df, // Use for small K correction
            getP:z=>{z=Math.abs(z);const t=1/(1+0.2316419*z);return 2*(0.3989423*Math.exp(-0.5*z*z)*((((1.330274*t-1.821256)*t+1.781478)*t-0.356564)*t+0.319382)*t);},

            calcEffects:(rows,type,met)=>{
                return rows.map((r,i)=>{
                    let es=0,vi=0,dv=0,lo=0,hi=0;
                    const n1=Engine.toNum(r.n1), n2=Engine.toNum(r.n2);

                    if(type==='binary'){
                        const e1=Engine.toNum(r.e1), e2=Engine.toNum(r.e2), cc=(e1===0||e2===0||e1===n1||e2===n2)?.5:0;
                        if(!n1||!n2)return null;
                        if(met==='RR'){ es=Math.log(((e1+cc)/(n1+cc))/((e2+cc)/(n2+cc))); vi=1/(e1+cc)-1/(n1+cc)+1/(e2+cc)-1/(n2+cc); dv=Math.exp(es); }
                        else if(met==='OR'){ es=Math.log(((e1+cc)*(n2-e2+cc))/((e2+cc)*(n1-e1+cc))); vi=1/(e1+cc)+1/(n1-e1+cc)+1/(e2+cc)+1/(n2-e2+cc); dv=Math.exp(es); }
                        else { const p1=(e1+cc)/(n1+2*cc), p2=(e2+cc)/(n2+2*cc); es=p1-p2; vi=(p1*(1-p1))/n1+(p2*(1-p2))/n2; dv=es; }
                    } else if(type==='cont'){
                        const m1=Engine.toNum(r.m1), s1=Engine.toNum(r.s1), m2=Engine.toNum(r.m2), s2=Engine.toNum(r.s2);
                        if(!n1||!n2)return null;
                        if(met==='MD'){ es=m1-m2; vi=(s1*s1)/n1+(s2*s2)/n2; dv=es; }
                        else { const df=n1+n2-2, sp=Math.sqrt(((n1-1)*s1*s1+(n2-1)*s2*s2)/df); es=(m1-m2)/sp*(1-(3/(4*df-1))); vi=(n1+n2)/(n1*n2)+(es*es)/(2*(n1+n2)); dv=es; }
                    } else if(type==='prop'){
                        const n=Engine.toNum(r.n), e=Engine.toNum(r.e); if(!n)return null;
                        const useCC = (e===0 || e===n);
                        const p_raw = e/n;
                        if(useCC){ const p_cc = (e+0.5)/(n+1); es = Math.log(p_cc/(1-p_cc)); vi = 1/((e+0.5)*(1-(e+0.5)/(n+1))); }
                        else { es = Math.log(p_raw/(1-p_raw)); vi = 1/(n*p_raw*(1-p_raw)); }
                        dv=e/n;
                    } else if(type==='survival'){
                        const hr=Engine.toNum(r.hr), l=Engine.toNum(r.lci), u=Engine.toNum(r.uci);
                        if(hr<=0)return null; es=Math.log(hr); const se=(Math.log(u)-Math.log(l))/3.92; vi=se*se; dv=hr;
                    }
                    if(!Engine.validNum(es)||!Engine.validNum(vi))return null;
                    const se=Math.sqrt(vi);
                    if(Metric.isLogRatio(met)||met==='Prop'){
                        if(met==='Prop'){ lo=Math.exp(es-1.96*se)/(1+Math.exp(es-1.96*se)); hi=Math.exp(es+1.96*se)/(1+Math.exp(es+1.96*se)); }
                        else { lo=Math.exp(es-1.96*se); hi=Math.exp(es+1.96*se); }
                    } else { lo=es-1.96*se; hi=es+1.96*se; }

                    return {id:r.id||'Study', subId:i, year:Engine.toNum(r.year), group:r.group, es, vi, se, displayVal:dv, lo, hi, raw:r, metric:met, rowIdx:i, excluded:!!r.excluded};
                }).filter(x=>x);
            },

            // --- MULTILEVEL SOLVER (IGLS + RVE Correction + Cheung's I2) ---
            pool:(eff,mod,meth,est,hk)=>{
                const act=eff.filter(e=>!e.excluded); if(!act.length)return null;
                const k = act.length;

                // 1. Identify Clusters (Studies)
                const clusters = {};
                act.forEach(e => {
                    if (!clusters[e.id]) clusters[e.id] = [];
                    clusters[e.id].push(e);
                });
                const m = Object.keys(clusters).length; // Number of Studies

                // 2. Iterative GLS (IGLS)
                let tb2 = 0.01, tw2 = 0.00;
                if(mod === 'fixed') { tb2 = 0; tw2 = 0; }

                // Vectors
                const Y = act.map(e => [e.es]); // k x 1
                const X = Matrix.zeros(k, 1); for(let i=0; i<k; i++) X[i][0] = 1;

                let beta, CovBeta, Q = 0, W_mat;

                for(let iter=0; iter<10; iter++){
                    // Construct Sigma Matrix (Block Diagonal)
                    const Sigma = Matrix.zeros(k, k);
                    for(let i=0; i<k; i++){
                        for(let j=0; j<k; j++){
                            let val = 0;
                            if(i === j) val += act[i].vi + tw2;
                            if(act[i].id === act[j].id) val += tb2; // Correlated if same study
                            Sigma[i][j] = val;
                        }
                    }

                    W_mat = Matrix.inv(Sigma);
                    const Xt = Matrix.transpose(X);
                    const XtW = Matrix.dot(Xt, W_mat);
                    const XtWX = Matrix.dot(XtW, X);

                    try {
                        CovBeta = Matrix.inv(XtWX);
                        beta = Matrix.dot(Matrix.dot(CovBeta, XtW), Y);
                    } catch(e) {
                        beta = [[act.reduce((a,b)=>a+b.es,0)/k]]; CovBeta = [[0.1]];
                    }

                    if(mod === 'fixed') break;

                    // Cluster-Based Variance Estimation (Moment Method)
                    const res = []; for(let i=0; i<k; i++) res.push(Y[i][0] - beta[0][0]);

                    let Q_new = 0; act.forEach((e,i) => Q_new += (1/e.vi) * (res[i]**2)); Q = Q_new;

                    // Estimate variances
                    const clusterMeans = Object.values(clusters).map(c => c.reduce((sum, e) => sum + e.es, 0) / c.length);
                    const grandMean = clusterMeans.reduce((a,b)=>a+b,0)/clusterMeans.length;
                    const varBetween = clusterMeans.reduce((a,b)=>a+(b-grandMean)**2,0) / (clusterMeans.length-1 || 1);

                    let ssWithin = 0; let dfWithin = 0;
                    Object.values(clusters).forEach(c => {
                        const cMean = c.reduce((sum, e) => sum + e.es, 0) / c.length;
                        c.forEach(e => ssWithin += (e.es - cMean)**2);
                        dfWithin += (c.length - 1);
                    });
                    const varWithin = ssWithin / (dfWithin || 1);

                    const totalVar = Math.max(0, (Q - (k-1))/act.reduce((a,b)=>a+1/b.vi,0));
                    const ratio = totalVar / (varBetween + varWithin || 1);
                    tb2 = varBetween * ratio; tw2 = varWithin * ratio;
                }

                const pooledES = beta[0][0];
                const pooledSE_Model = Math.sqrt(CovBeta[0][0]);

                // --- ADVANCED METHOD 1: ROBUST VARIANCE ESTIMATION (SANDWICH) ---
                // Meat = Sum[ Xj' Wj rj rj' Wj Xj ] with small sample correction
                let Meat = 0;
                const smallSampleCorr = Math.sqrt(m / (m - 1)); // Tipton's simple correction

                // Re-calculate residuals based on final beta
                const residuals = act.map(e => e.es - pooledES);

                // Iterate per cluster for Sandwich construction
                let currentIdx = 0;
                Object.keys(clusters).forEach(cId => {
                    const cSize = clusters[cId].length;
                    // Extract sub-matrices for this cluster
                    const X_j = Matrix.zeros(cSize, 1); for(let i=0; i<cSize; i++) X_j[i][0] = 1;
                    const r_j = Matrix.zeros(cSize, 1); for(let i=0; i<cSize; i++) r_j[i][0] = residuals[currentIdx+i] * smallSampleCorr;

                    // W_j is block from W_mat
                    const W_j = Matrix.zeros(cSize, cSize);
                    for(let i=0; i<cSize; i++){
                        for(let j=0; j<cSize; j++){
                            W_j[i][j] = W_mat[currentIdx+i][currentIdx+j];
                        }
                    }

                    // Calculate term: Xj' Wj rj rj' Wj Xj
                    // Note: Xj' Wj is 1xSize. rj rj' is Size x Size.
                    const XtW_j = Matrix.dot(Matrix.transpose(X_j), W_j); // 1xSize
                    const rrt = Matrix.dot(r_j, Matrix.transpose(r_j));   // Size x Size
                    const term = Matrix.dot(Matrix.dot(Matrix.dot(XtW_j, rrt), W_j), X_j); // 1x1

                    Meat += term[0][0];
                    currentIdx += cSize;
                });

                // Sandwich = V_model * Meat * V_model
                const Var_Robust = CovBeta[0][0] * Meat * CovBeta[0][0];
                const SE_Robust = Math.sqrt(Var_Robust);

                // Calculate CIs (Model vs Robust)
                const tVal = m > 1 ? Engine.getT(m - 1) : 1.96;
                const meta = act[0].metric;
                const isLog = Metric.isLogRatio(meta) || meta==='Prop';

                // Function to transform back
                const trans = (val, se) => {
                    if(isLog){
                        if(meta==='Prop'){
                            return {
                                v: Math.exp(val)/(1+Math.exp(val)),
                                l: Math.exp(val-tVal*se)/(1+Math.exp(val-tVal*se)),
                                h: Math.exp(val+tVal*se)/(1+Math.exp(val+tVal*se))
                            };
                        }
                        return { v: Math.exp(val), l: Math.exp(val-tVal*se), h: Math.exp(val+tVal*se) };
                    }
                    return { v: val, l: val-tVal*se, h: val+tVal*se };
                };

                const resModel = trans(pooledES, pooledSE_Model);
                const resRobust = trans(pooledES, SE_Robust);

                // --- ADVANCED METHOD 2: CHEUNG'S I2 DECOMPOSITION ---
                const typVi = act.reduce((a,b)=>a+b.vi,0) / k; // Simple typical sampling variance
                const denom = tb2 + tw2 + typVi;
                const I2_3 = (tb2 / denom) * 100; // Between-Cluster
                const I2_2 = (tw2 / denom) * 100; // Within-Cluster
                const I2_Total = ((tb2 + tw2) / denom) * 100;

                // Visual Weights
                act.forEach((e) => {
                    e.w = 1 / (e.vi + tw2 + tb2);
                    e.q = (e.es - pooledES)**2 / e.vi;
                    e.imp = e.w * e.q;
                });

                return {
                    es:pooledES, se:pooledSE_Model, se_r: SE_Robust,
                    I2: I2_Total, I2_2, I2_3, Q, tau2:tb2+tw2, tau2_b:tb2, tau2_w:tw2,
                    studies:eff,
                    displayVal:resModel.v, lo:resModel.l, hi:resModel.h, // Model based
                    robVal:resRobust.v, robLo:resRobust.l, robHi:resRobust.h, // RVE based
                    pVal:Engine.getP(pooledES/pooledSE_Model), k_studies:m, k_effects:k
                };
            },

            // --- ADVANCED METHOD 3: MULTILEVEL EGGER'S TEST (Meta-Regression) ---
            multilevelEgger: (eff) => {
                const act = eff.filter(e => !e.excluded);
                if(act.length < 3) return {p:null, intercept:0};

                const k = act.length;
                const Y = act.map(e => [e.es / e.se]); // Standardized Effect
                const X = Matrix.zeros(k, 2);
                for(let i=0; i<k; i++) {
                    X[i][0] = 1 / act[i].se; // Precision (Estimate for slope = original intercept)
                    X[i][1] = 1;        // Intercept (Estimate for Egger's bias term)
                }

                // Simple OLS for Egger's (sufficient for bias detection)
                try {
                    const Xt = Matrix.transpose(X);
                    const XtX = Matrix.dot(Xt, X);
                    const XtY = Matrix.dot(Xt, Y);
                    const beta = Matrix.dot(Matrix.inv(XtX), XtY);

                    // Calculate P for the intercept term (beta[1])
                    const residuals = [];
                    for(let i=0; i<k; i++) residuals.push(Y[i][0] - (beta[0][0]*X[i][0] + beta[1][0]*X[i][1]));
                    const sse = residuals.reduce((a,b)=>a+b*b,0);
                    const mse = sse / (k-2);
                    const varCov = Matrix.scalarMult(Matrix.inv(XtX), mse);
                    const se_intercept = Math.sqrt(varCov[1][1]);
                    const t = beta[1][0] / se_intercept;

                    return { p: Engine.getP(t), intercept: beta[1][0] };
                } catch(e) { return {p:null, intercept:0}; }
            }
        };

if (typeof module!=='undefined'&&module.exports){ module.exports = { Matrix, Metric, Engine }; }
