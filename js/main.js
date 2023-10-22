import * as param from "./variables/const.js"
import { v_default, ksu, grad, rad } from "./functions/utils.js"
import { p_0, po_0, T_0 } from "./functions/param_0.js"
import { secantMethod } from "./methods/secant-method.js"
import { p_ksu, po_ksu, T_ksu, v_ksu, a_ksu, M_ksu, um_ksu } from "./functions/param_ksu.js"
import { w, w_next } from "./functions/prandtl.js"
import { tau, pi, eps } from "./functions/gasodynamic.js"
import { cp_koeff, u_koeff, lymbda_koeff } from "./functions/koeff.js"
import { betta, thet, M, p, po, T, a, v, um, p0, po0, T0, cp, u, lymbda } from "./variables/arr.js"
import { render } from "./view/html.js";
import { X, Y, Mk, Cxa, Cya, Cx, Cy, mz, Mz } from "./functions/ADH.js"

for (let i = 0; i <= 2; i++) {
    if (i < 2) {
        p.push(p_ksu(M[0], thet[i], p[0]))
        po.push(po_ksu(M[0], thet[i], po[0]))
        v.push(v_ksu(v[0], thet[i], betta[i]))
        T.push(T_ksu(T[0], p[i + 1], po[0], p[0], po[i + 1]))
        a.push(a_ksu(p[i + 1], po[i + 1]))
        M.push(M_ksu(v[i + 1], a[i + 1]))
        um.push(um_ksu(M[i + 1]))
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
    let p_arr = [1, 0];
    let thet_arr = [0, 0];
    while (math.abs(p_arr[0] - p_arr[1]) / p_arr[0] >= math.pow(10, -4)) {
        let init = rad(10)
        const first = ksu(betta[2] + delta, M[3])
        const second = ksu(betta[3] - delta, M[4])
        thet_arr[0] = secantMethod(first, init)
        thet_arr[1] = secantMethod(second, init)
        p_arr[0] = p_ksu(M[3], thet_arr[0], p[3])
        p_arr[1] = p_ksu(M[4], thet_arr[1], p[4])
        delta += rad(math.pow(10, -4))
    }
    for (let i = 0; i <= 1; i++) {
        p.push(p_arr[i])
        thet.push(thet_arr[i])
    }
    betta.push(betta[2] + delta)
    betta.push(betta[2] - delta)

})()

for (let i = 3; i <= 4; i++) {
    const po_t = po_ksu(M[i], thet[i - 1], po[i])
    const T_t = T_ksu(T[i], p[i + 2], po[i], p[i], po_t)
    const v_t = v_ksu(v[i], thet[i - 1], betta[i + 1])
    const a_t = a_ksu(p[i + 2], po_t)
    const M_t = M_ksu(v_t, a_t)
    po.push(po_t)
    T.push(T_t)
    v.push(v_t)
    a.push(a_t)
    M.push(M_t)
    um.push(um_ksu(M_t))
}

T.forEach(Tk => {
    cp.push(cp_koeff(Tk, T[0]))
    u.push((u_koeff(Tk, T[0])))
    lymbda.push(lymbda_koeff(Tk, T[0]))
});

const arr = [M, p, T, po, cp, u, lymbda]
const thet__arr = [thet.map(el => grad(el))]
const angle__arr = [um.map(el => grad(el)), betta.map(el => grad(el))]
const data__header = ['M, [-]', 'p, [Па]', 'T, [K]', 'po, [кг/м3]', 'Cp, [Дж/кг*К]', 'μ, [Па*с]', 'λ, [Вт/м*К]']
const thet__header = ['θ, [град]']
const angle__header = ['μ, [град]', 'β, [град]']
render(arr, data__header, 0)
render(thet__arr, thet__header, 1)
render(angle__arr, angle__header, 1)

const getOprTemp = (flow, T, M) => {
    const Cp_opr = T => 1004.7 * math.pow(T / 288.15, 0.1);
    const u_opr = T => 1.79 * math.pow(10, -5) * math.pow(T / 288.15, 0.76);
    const lymbda_opr = T => 0.0232 * math.pow(T / 288.15, 0.86);

    const r_opr = (Cp, u, lymbda, extent) => math.pow(Cp * u / lymbda, extent);
    const Tr_opr = (T, r, M) => T * (1 + r * (param.k - 1) * M * M / 2);
    const T_opr = (Tr, T) => (param.T_st + T) / 2 + 0.22 * (Tr - T)
    let T_z = (param.T_st + T) / 2;
    let T_z_previous = 0;
    let Tr = 0;
    let rl = 0.83;
    let rt = 0.88;
    let r = rl;
    let extent = 1 / 2;
    if (flow === 'turb') {
        r = rt;
        extent = 1 / 3;
    }
    while (math.abs(T_z - T_z_previous) > (T_z * 0.01)) {
        T_z_previous = T_z;
        Tr = Tr_opr(T, r, M);
        T_z = T_opr(Tr, T);
        r = r_opr(Cp_opr(T_z), u_opr(T_z), lymbda_opr(T_z), extent);
    }
    return [T_z, Tr]
}

let T_opr_laminar = [];
let Tr_laminar = [];
let T_opr_turb = [];
let Tr_turb = [];
for (let i = 1; i <= 4; i++) {
    const [T_iter_laminar_opr, Tr_iter_laminar] = getOprTemp('laminar', T[i], M[i]);
    const [T_iter_turb_opr, Tr_iter_turb] = getOprTemp('turb', T[i], M[i]);
    T_opr_laminar.push(T_iter_laminar_opr);
    T_opr_turb.push(T_iter_turb_opr);
    Tr_turb.push(Tr_iter_turb);
    Tr_laminar.push(Tr_iter_laminar);
}
// const Tr_opr_laminar = T_opr_laminar.map(el => el.at(0));
const T_opr_header = ['Tr (л), [K]', 'T* (л), [K]', 'Tr (т), [K]', 'T* (т), [K]'];
render([Tr_laminar, T_opr_laminar, Tr_turb, T_opr_turb], T_opr_header, 1);

const C_x = Cx(X(p[1], p[0]), X(p[2], p[0]), X(p[3], p[0]), X(p[4], p[0]));
const C_y = Cy(Y(p[1], p[0]), Y(p[2], p[0]), Y(p[3], p[0]), Y(p[4], p[0]));
const m_z = mz(Mk(X(p[1], p[0]), Y(p[1], p[0])), Mk(-X(p[2], p[0]), -Y(p[2], p[0])), Mk(-X(p[3], p[0]), 3 * Y(p[3], p[0])), Mk(X(p[4], p[0]), -3 * Y(p[4], p[0])));
const M_z = Mz(m_z);
const C_xa = Cxa(C_x, C_y);
const C_ya = Cya(C_x, C_y);
const adh_header = ['Cx, [-]', 'Cy, [-]', 'mz, [-]', 'Cxa, [-]', 'Cya, [-]', 'Mz, [-]']
render([[C_x], [C_y], [m_z], [C_xa], [C_ya], [M_z]], adh_header, 0)

const criticalReynoldsNumber = (Tr, M) => {
    const Tst_otn = param.T_st / Tr; // Относительная температура стенки;
    const x = (Tst_otn - 1) / (M * M);
    const y = 1 - 16 * x - 412.5 * math.pow(x, 2) - 35000 * math.pow(x, 3) - 375000 * math.pow(x, 4);
    return 5 * math.pow(10, 6) * y
}
const reynoldsNumber = (po, V, u) => po * V * param.L / u;
const Re_l1 = reynoldsNumber(po[1], v[1], u[1]);
const Re_l2 = reynoldsNumber(po[2], v[2], u[2]);
console.log(Re_l1, Re_l2);

const Re_kr1 = criticalReynoldsNumber(Tr_laminar[0], M[1]);
const Re_kr2 = criticalReynoldsNumber(Tr_laminar[1], M[2]);

const x_kr1 = Re_kr1 * u[1] / (v[1] * po[1]);
const x_kr2 = Re_kr2 * u[2] / (v[2] * po[2]);
render([[Re_kr1, Re_l1], [Re_kr2, Re_l2]], ['Re_кр1 / Re_l1, [-]', 'Re_кр2 / Re_l2, [-]'], 0);
render([[x_kr1], [x_kr2]], ['Xкр1, [м]', 'Xкр2, [м]'], 0)