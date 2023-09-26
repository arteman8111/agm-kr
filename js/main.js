import * as param from "./const.js"
import { atm } from "./atm.js"
import { v_default, ksu, grad, rad } from "./utils.js"
import { p_0, po_0, T_0 } from "./param_0.js"
import { secantMethod } from "./methods.js"
import { p_ksu, po_ksu, T_ksu, v_ksu, a_ksu, M_ksu, um_ksu } from "./param_ksu.js"
import { w, w_next } from "./prandtl.js"
import { tau, pi, eps } from "./gasodynamic.js"
import { cp_koeff, u_koeff, lymbda_koeff } from "./koeff.js"
import { render } from "./html.js"

const atmos = atm(param.h);
let betta = [param.betta_k - rad(param.alfa), param.betta_k + rad(param.alfa), param.betta_k, param.betta_k];
let thet = [secantMethod(ksu(betta[0], param.M_inf), rad(12)), secantMethod(ksu(betta[1], param.M_inf), rad(12))];

let M = [param.M_inf];
let p = [atmos.p];
let po = [atmos.po];
let T = [atmos.T];
let a = [atmos.a]
let v = [v_default(M[0], a[0])];
let um = [];

let p0 = [];
let po0 = [];
let T0 = [];

let cp = [];
let u = [];
let lymbda = [];

for (let i = 0; i <= 2; i++) {
    if (i < 2){
        p.push(p_ksu(M[0], thet[i], p[0]))
        po.push(po_ksu(M[0], thet[i], po[0]))
        v.push(v_ksu(v[0], thet[i], betta[i]))
        T.push(T_ksu(T[0], p[i+1], po[0], p[0], po[i+1]))
        a.push(a_ksu(p[i+1], po[i+1]))
        M.push(M_ksu(v[i+1], a[i+1]))
        um.push(um_ksu(M[i+1]))
    }
    p0.push(p_0(M[i], p[i]));
    po0.push(po_0(M[i], po[i]));
    T0.push(T_0(M[i], T[i]));
}
for (let i = 1; i <= 2; i++) {
    const w0 = w(M[i]);
    const prandtl = (M) => w_next(w0) - w(M);
    const Mt = secantMethod(prandtl, 4);
    const pt = p0[i] / pi(Mt);
    const po_t = po0[i] / eps(Mt);
    const T_t = T0[i] / tau(Mt);
    const a_t = a_ksu(pt, po_t);
    M.push(Mt);
    p.push(pt)
    po.push(po_t)
    T.push(T_t)
    a.push(a_t)
    v.push(v_default(Mt, a_t))
    um.push(um_ksu(Mt))
}
(function posled() {
    let delta = rad(0);
    let p_arr = [1,0];
    let thet_arr = [0,0];
    while(math.abs(p_arr[0]-p_arr[1])/p_arr[0] >= 0.03){
        let init = rad(10)
        const first = ksu(betta[2] + delta, M[3])
        const second = ksu(betta[3] - delta, M[4])
        thet_arr[0] = secantMethod(first, init)
        thet_arr[1] = secantMethod(second, init)
        p_arr[0] = p_ksu(M[3],thet_arr[0],p[3])
        p_arr[1] = p_ksu(M[4],thet_arr[1],p[4])
        delta += rad(0.1)
    }
    for(let i = 0; i <=1; i++){
        p.push(p_arr[i])
        thet.push(thet_arr[i])
    }
    betta.push(betta[2] + delta)
    betta.push(betta[2] - delta)

})()

for (let i = 3; i <= 4; i++){
    const po_t = po_ksu(M[0],thet[0],po[i])
    const T_t = T_ksu(T[i],p[i+2],po[i],p[i],po_t)
    const v_t = v_ksu(v[i],thet[i-1],betta[i+1])
    const a_t = a_ksu(p[i+2],po_t)
    const M_t = M_ksu(v_t,a_t)
    po.push(po_t)
    T.push(T_t)
    v.push(v_t)
    a.push(a_t)
    M.push(M_t)
    um.push(um_ksu(M_t))
}

T.forEach(Tk => {
    cp.push(cp_koeff(Tk,T[0]))
    u.push((u_koeff(Tk,T[0])))
    lymbda.push(lymbda_koeff(Tk,T[0]))
});

const arr = [M,p,T,po,cp,u,lymbda]
const thet__arr = [thet.map(el=>grad(el))]
const angle__arr = [um.map(el=>grad(el)),betta.map(el=>grad(el))]
const data__header = ['M, [-]','p, [Па]','T, [K]','po, [кг/м3]','Cp, [Дж/кг*К]','μ, [Па*с]','λ, [Вт/м*К]']
const thet__header = ['θ, [град]']
const angle__header = ['μ, [град]','β, [град]']
render(arr, data__header,0)
render(thet__arr,thet__header,1)
render(angle__arr,angle__header,1)