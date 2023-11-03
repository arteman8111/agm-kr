import * as param from "./variables/const.js"
import * as ps from "./functions/paramPS.js";
import { render } from "./view/html.js";
import { grad } from "./functions/utils.js";
import { M, p, T, po, v, a, cp, u, lymbda, T0, p0, po0, thet, um, betta, delta } from "./functions/6Oblastey.js"
import { Tr_laminar, T_opr_laminar, Tr_turb, T_opr_turb } from "./functions/oprTemp.js";
import { Cx, Cy, mz, Cxa, Cya, K, Cd } from "./functions/ADH.js";

// Последний этап АГМ
const criticalReynoldsNumber = (Tr, M) => {
    const Tst_otn = param.T_st / Tr; // Относительная температура стенки;
    const x = (Tst_otn - 1) / (M * M);
    const y = 1 - 16 * x - 412.5 * math.pow(x, 2) - 35000 * math.pow(x, 3) - 375000 * math.pow(x, 4);
    return 5 * math.pow(10, 6) * y
}
const reynoldsNumber = (po, V, u, x = param.L) => po * V * x / u;
const Re_l1 = reynoldsNumber(po[1], v[1], u[1]);
const Re_l2 = reynoldsNumber(po[2], v[2], u[2]);
const Re_kr1 = criticalReynoldsNumber(Tr_laminar[0], M[1]);
const Re_kr2 = criticalReynoldsNumber(Tr_laminar[1], M[2]);
const x_kr1 = Re_kr1 * u[1] / (v[1] * po[1]);
const x_kr2 = Re_kr2 * u[2] / (v[2] * po[2]);

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
    let delta_xk = ps.delta_x(delta_[1], v[inc], po[inc], u[inc]);
    let b = param.L - x_kr1;
    xarr.push(delta_xk, delta_xk + b / 3, delta_xk + 2 * b / 3, delta_xk + b);
    let xarr3_6 = [delta_xk, delta_xk + b / 3, delta_xk + 2 * b / 3, delta_xk + b]
    xarr3_6.forEach(
        (x, i) => {
            const Re = ps.Re(v[inc], po[inc], x, u[inc]);
            delta_ns.push(ps.delta_t_ns(x, Re));
            tau_ns.push(ps.tau_t_ns(po[inc], v[inc], Re));
            Cfx_ns.push(ps.Cfx_t_ns(Re));
            Cf_ns.push(ps.Cf_t_ns(Re));
            delta_ns_zv.push(ps.delta_t_ns_zv(Re, x));
            delta_ns_zv_zv.push(ps.delta_t_ns_zv_zv(Re, x));

            delta_.push(ps.delta_t(delta_ns[i], T_opr_turb[inc - 1], T[inc]));
            tau.push(ps.tau_t(tau_ns[i], T_opr_turb[inc - 1], T[inc]));
            Cfx.push(ps.Cfx_t(Cfx_ns[i], T_opr_turb[inc - 1], T[inc]));
            Cf.push(ps.Cf_t(Cf_ns[i], T_opr_turb[inc - 1], T[inc]));
            delta_zv.push(ps.delta_t_zv(delta_ns_zv[i], delta_[i], delta_ns[i]));
            delta_zv_zv.push(ps.delta_t_zv_zv(delta_ns_zv_zv[i], Cfx[i], Cfx_ns[i]))
        }
    )
    let delta_xk2 = ps.delta_x(delta_[5], v[inc + 2], po[inc + 2], u[inc + 2]);
    xarr.push(delta_xk2, param.L / 2 + delta_xk2, param.L + delta_xk2);
    let xarr7_9 = [delta_xk2, param.L / 2 + delta_xk2, param.L + delta_xk2];
    xarr7_9.forEach(
        (x, i) => {
            const Re = ps.Re(v[inc + 2], po[inc + 2], x, u[inc + 2]);
            delta_ns.push(ps.delta_t_ns(x, Re));
            tau_ns.push(ps.tau_t_ns(po[inc + 2], v[inc + 2], Re));
            Cfx_ns.push(ps.Cfx_t_ns(Re));
            Cf_ns.push(ps.Cf_t_ns(Re));
            delta_ns_zv.push(ps.delta_t_ns_zv(Re, x));
            delta_ns_zv_zv.push(ps.delta_t_ns_zv_zv(Re, x));

            delta_.push(ps.delta_t(delta_ns[i], T_opr_turb[inc + 1], T[inc + 2]));
            tau.push(ps.tau_t(tau_ns[i], T_opr_turb[inc + 1], T[inc + 2]));
            Cfx.push(ps.Cfx_t(Cfx_ns[i], T_opr_turb[inc + 1], T[inc + 2]));
            Cf.push(ps.Cf_t(Cf_ns[i], T_opr_turb[inc + 1], T[inc + 2]));
            delta_zv.push(ps.delta_t_zv(delta_ns_zv[i], delta_[i], delta_ns[i]));
            delta_zv_zv.push(ps.delta_t_zv_zv(delta_ns_zv_zv[i], Cfx[i], Cfx_ns[i]))
        }
    )
    let xarr2 = [x_kr / 2, x_kr, x_kr, param.L - (param.L - x_kr) / 3, x_kr + (param.L - x_kr) / 3, param.L, 0, param.L / 2, param.L];
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

const ps_for_1 = getSloy(x_kr1, 1);
const ps_for_2 = getSloy(x_kr2, 2);
console.log(ps_for_1);
const Cf_sm = (Re_kr, x_kr, delta_x, i) => {
    const ReL = ps.Re(v[i], po[i], param.L, u[i]);
    const Re_dx = ps.Re(v[i], po[i], delta_x, u[i]);
    const Cfsr_l = 1.3 / math.sqrt(Re_kr);
    const Cfsr_t = 0.074 / math.pow(ReL, 1 / 5);
    const Cfsr_t2 = 0.074 / math.pow(Re_dx, 1 / 5)
    const l = param.L - x_kr + delta_x

    return Cfsr_l * x_kr / param.L + Cfsr_t * l / param.L - Cfsr_t2 * delta_x / param.L
}

let delta_x_1_1 = ps_for_1[1][2];
let delta_x_1_2 = ps_for_1[1][6];
let delta_x_2_1 = ps_for_2[1][2];
let delta_x_2_2 = ps_for_2[1][6];

const Cf_sm1_1 = Cf_sm(Re_kr1, x_kr1, delta_x_1_1, 1);
const Cf_sm1_2 = Cf_sm(Re_kr2, x_kr2, delta_x_1_2, 3);
const Cf_sm2_1 = Cf_sm(Re_kr1, 0, delta_x_2_1, 2);
const Cf_sm2_2 = Cf_sm(Re_kr2, 0, delta_x_2_2, 4);


render([M, p, T, po, cp, u, lymbda], ['M, [-]', 'p, [Па]', 'T, [K]', 'po, [кг/м3]', 'Cp, [Дж/кг*К]', 'μ, [Па*с]', 'λ, [Вт/м*К]'], 0, 'Параметры в 6 областях')
render([T0, p0, po0], ['T0, [K]', 'p0, [Па]', 'po, [кг / м3]'], 0, 'Параметры торможения')
render([thet.map(el => grad(el))], ['θ, [град]'], 1, 'Угол наколна СУ')
render([um.map(el => grad(el)), betta.map(el => grad(el))], ['μ, [град]', 'β, [град]'], 1, 'Угол наклона маха и поворота потока')
render([[delta * 180 / math.pi]], ['delta, [град]'], 6, 'Отклонение от оси в 5-6 областях')
render([Tr_laminar, T_opr_laminar, Tr_turb, T_opr_turb], ['Tr (л), [K]', 'T* (л), [K]', 'Tr (т), [K]', 'T* (т), [K]'], 1, 'Определяющая температура');
render([[Cx], [Cy], [mz], [Cxa], [Cya], [K], [Cd]], ['Cx, [-]', 'Cy, [-]', 'mz, [-]', 'Cxa, [-]', 'Cya, [-]', 'K, [-]', 'Cd, [-]'], 0, 'АДХ')
render([[Re_kr1, Re_l1], [Re_kr2, Re_l2]], ['Re_кр1 / Re_l1, [-]', 'Re_кр2 / Re_l2, [-]'], 0, 'Критическое / Обычное число Рейнольдса');
render([[x_kr1], [x_kr2]], ['Xкр1, [м]', 'Xкр2, [м]'], 0, '1-ая и 2-ая критическая точка в смешанном ПС');

render(ps_for_1, ['X', 'Xф', 'δн', 'τн', 'Cfxн', 'Cfн', 'δ*н', 'δ**н', 'δ', 'τ', 'Cfx', 'Cf', 'δ*', 'δ**'], 1, 'Параметры ЛПС и ТПС 1 и 3 грани')
render(ps_for_2, ['X', 'Xф', 'δн', 'τн', 'Cfxн', 'Cfн', 'δ*н', 'δ**н', 'δ', 'τ', 'Cfx', 'Cf', 'δ*', 'δ**'], 1, 'Параметры ЛПС и ТПС 2 и 4 грани')

render([[Cf_sm1_1], [Cf_sm1_2], [Cf_sm2_1], [Cf_sm2_2]], ['1 грань', '2 грань', '3 грань', '4 грань'], 1, "Средние коэфф трения ПС")