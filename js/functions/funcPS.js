import * as param from "../variables/const.js"
import * as ps from "../functions/paramPS.js";
import { M, T, po, v, u} from "../functions/6Oblastey.js"
import { Tr_laminar, T_opr_laminar, T_opr_turb } from "../functions/oprTemp.js";
import { u_koeff } from "../functions/koeff.js";
import { criticalReynoldsNumber } from "../functions/utils.js";

const getSloy = (x_kr, inc) => {
    let xarr = [x_kr / 2, x_kr];
    let delta_ns = [];
    let tau_ns = [];
    let Cfx_ns = [];
    let Cf_ns = [];
    let delta_ns_zv = [];
    let delta_ns_zv_zv = [];
    let delta_ = [];
    let tau = [];
    let Cfx = [];
    let Cf = [];
    let delta_zv = [];
    let delta_zv_zv = [];

    xarr.forEach(
        (x, i) => {
            const Re = ps.Re(v[inc], po[inc], x, u[inc]);
            delta_ns.push(ps.delta_l_ns(x, Re));
            tau_ns.push(ps.tau_l_ns(po[inc], v[inc], Re));
            Cfx_ns.push(ps.Cfx_l_ns(Re));
            Cf_ns.push(ps.Cf_l_ns(Re));
            delta_ns_zv.push(ps.delta_l_ns_zv(Re, x));
            delta_ns_zv_zv.push(ps.delta_l_ns_zv_zv(Re, x));

            delta_.push(ps.delta_l(delta_ns[i], T_opr_laminar[inc - 1], T[inc]));
            tau.push(ps.tau_l(tau_ns[i], T_opr_laminar[inc - 1], T[inc]));
            Cfx.push(ps.Cfx_l(Cfx_ns[i], T_opr_laminar[inc - 1], T[inc]));
            Cf.push(ps.Cf_l(Cf_ns[i], T_opr_laminar[inc - 1], T[inc]));
            delta_zv.push(ps.delta_l_zv(delta_ns_zv[i], delta_[i], delta_ns[i]));
            delta_zv_zv.push(ps.delta_l_zv_zv(delta_ns_zv_zv[i], Cfx[i], Cfx_ns[i]))
        }
    )
    let po_zv = po[inc] * T[inc] / T_opr_laminar[inc - 1];
    let u_zv = u_koeff(T_opr_laminar[inc - 1]);
    let x_kr1 = criticalReynoldsNumber(Tr_laminar[inc - 1], M[inc]) * u[inc] / (po[inc] * v[inc]);
    let Re_krit = po_zv * v[inc] * x_kr1 / u_zv;
    let delta_xk = math.pow((4.64 * x_kr1 * math.pow(po_zv * v[inc] / u_zv, 0.2)) / (math.sqrt(Re_krit) * 0.37), 1.25)
    let b = param.L - x_kr;
    xarr.push(delta_xk, delta_xk + b / 3, delta_xk + 2 * b / 3, delta_xk + b);
    let xarr3_6 = [delta_xk, delta_xk + b / 3, delta_xk + 2 * b / 3, delta_xk + b]
    xarr3_6.forEach(
        (x, i) => {
            const Re = ps.Re(v[inc], po[inc], x, u[inc]);
            delta_ns.push(ps.delta_t_ns(x, Re));
            tau_ns.push(ps.tau_t_ns(po[inc], v[inc], Re));
            Cfx_ns.push(ps.Cfx_t_ns(Re));
            Cf_ns.push(ps.Cf_t_ns(Re));
            delta_ns_zv.push(ps.delta_t_ns_zv(delta_ns[i + 2]));
            delta_ns_zv_zv.push(ps.delta_t_ns_zv_zv(delta_ns[i + 2]));

            delta_.push(ps.delta_t(delta_ns[i + 2], T_opr_turb[inc - 1], T[inc]));
            tau.push(ps.tau_t(tau_ns[i + 2], T_opr_turb[inc - 1], T[inc]));
            Cfx.push(ps.Cfx_t(Cfx_ns[i + 2], T_opr_turb[inc - 1], T[inc]));
            Cf.push(ps.Cf_t(Cf_ns[i + 2], T_opr_turb[inc - 1], T[inc]));
            delta_zv.push(ps.delta_t_zv(delta_ns_zv[i + 2], delta_[i + 2], delta_ns[i + 2]));
            delta_zv_zv.push(ps.delta_t_zv_zv(delta_ns_zv_zv[i + 2], Cfx[i + 2], Cfx_ns[i + 2]))
        }
    )
    let po_zv_T = po[inc + 2] * T[inc + 2] / T_opr_turb[inc + 1];
    let u_zv_T = u_koeff(T_opr_turb[inc + 1]);
    let delta_xk2 = math.pow(math.pow(delta_[5], 5) * po_zv_T * v[inc + 2] / (u_zv_T * math.pow(0.37, 5)), 1 / 4);
    xarr.push(delta_xk2, param.L / 2 + delta_xk2, param.L + delta_xk2);
    let xarr7_9 = [delta_xk2, param.L / 2 + delta_xk2, param.L + delta_xk2];
    xarr7_9.forEach(
        (x, i) => {
            const Re = ps.Re(v[inc + 2], po[inc + 2], x, u[inc + 2]);
            delta_ns.push(ps.delta_t_ns(x, Re));
            tau_ns.push(ps.tau_t_ns(po[inc + 2], v[inc + 2], Re));
            Cfx_ns.push(ps.Cfx_t_ns(Re));
            Cf_ns.push(ps.Cf_t_ns(Re));
            delta_ns_zv.push(ps.delta_t_ns_zv(delta_ns[i + 6]));
            delta_ns_zv_zv.push(ps.delta_t_ns_zv_zv(delta_ns[i + 6]));

            delta_.push(ps.delta_t(delta_ns[i + 6], T_opr_turb[inc + 1], T[inc + 2]));
            tau.push(ps.tau_t(tau_ns[i + 6], T_opr_turb[inc + 1], T[inc + 2]));
            Cfx.push(ps.Cfx_t(Cfx_ns[i + 6], T_opr_turb[inc + 1], T[inc + 2]));
            Cf.push(ps.Cf_t(Cf_ns[i + 6], T_opr_turb[inc + 1], T[inc + 2]));
            delta_zv.push(ps.delta_t_zv(delta_ns_zv[i + 6], delta_[i + 6], delta_ns[i + 6]));
            delta_zv_zv.push(ps.delta_t_zv_zv(delta_ns_zv_zv[i + 6], Cfx[i + 6], Cfx_ns[i + 6]))
        }
    )
    let xarr2 = [
        x_kr / 2,
        x_kr,
        x_kr,
        x_kr + (param.L - x_kr) / 3,
        x_kr + 2 * (param.L - x_kr) / 3,
        param.L,
        0,
        param.L / 2,
        param.L
    ];
    return [
        xarr2,
        xarr,
        delta_ns,
        tau_ns,
        Cfx_ns,
        Cf_ns,
        delta_ns_zv,
        delta_ns_zv_zv,
        delta_,
        tau,
        Cfx,
        Cf,
        delta_zv,
        delta_zv_zv
    ]
}

export {getSloy}