import * as param from "../variables/const.js"
import {M, T} from "../functions/6Oblastey.js"


const getOprTemp = (flow, T, M) => {
    const Cp_opr = T => 1004.7 * math.pow(T / 288.15, 0.1);
    const u_opr = T => 1.79 * math.pow(10, -5) * math.pow(T / 288.15, 0.76);
    const lymbda_opr = T => 0.0232 * math.pow(T / 261, 0.86);

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

export {
    T_opr_laminar,
    Tr_laminar,
    T_opr_turb,
    Tr_turb
}