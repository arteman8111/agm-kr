import * as param from "./variables/const.js"
import { myId, render } from "./view/html.js";
import { grad } from "./functions/utils.js";
import { M, p, T, po, v, a, cp, u, lymbda, T0, p0, po0, thet, um, betta, delta } from "./functions/6Oblastey.js"
import { Tr_laminar, T_opr_laminar, Tr_turb, T_opr_turb } from "./functions/oprTemp.js";
import { Cx, Cy, mz, Cxa, Cya, K, Cd } from "./functions/ADH.js";
import { criticalReynoldsNumber, reynoldsNumber } from "./functions/utils.js";
import { getSloy } from "./functions/funcPS.js";

// Последний этап АГМ
const Cfsm = (ps) => {
    return (ps[11][1] * ps[1][1] + ps[11][5] * ps[1][5] - ps[11][2] * ps[1][2]) / param.L;
}
const Cfturb = (ps) => {
    return (ps[11][8] * ps[1][8] - ps[11][6] * ps[1][6]) / param.L
}

const Re_l1 = reynoldsNumber(po[1], v[1], u[1]);
const Re_l2 = reynoldsNumber(po[2], v[2], u[2]);
const Re_kr1 = criticalReynoldsNumber(Tr_laminar[0], M[1]);
const Re_kr2 = criticalReynoldsNumber(Tr_laminar[1], M[2]);
const x_kr1 = Re_kr1 * u[1] / (v[1] * po[1]);
const x_kr2 = Re_kr2 * u[2] / (v[2] * po[2]);
const ps13 = getSloy(x_kr1, 1);
const ps24 = getSloy(x_kr2, 2);
const Cf1 = Cfsm(ps13);
const Cf2 = Cfsm(ps24);
const Cf3 = Cfturb(ps13);
const Cf4 = Cfturb(ps24);

const q = (po, v) => {
    return po * v * v / 2;
}

const newADH = () => {
    const Cxf = (Cf1 * q(po[1], v[1]) + Cf2 * q(po[2], v[2]) + Cf3 * q(po[3], v[3]) + Cf4 * q(po[4], v[4])) / (2 * q(po[0], v[0]));
    const Cx_new = Cx + Cxf;
    const Cya_new = Cy * math.cos(param.alfa) - Cx_new * math.sin(param.alfa);
    const Cxa_new = Cy * math.sin(param.alfa) + Cx_new * math.cos(param.alfa);
    const K_new = Cya_new / Cxa_new;

    return [
        [Cx_new], [Cy], [mz], [Cxa_new], [Cya_new], [K_new], [Cd]
    ]
}

function init() {
    const app = document.querySelector('#app');
    myId(app);
    render([M, p, T, po, cp, u, lymbda], ['M, [-]', 'p, [Па]', 'T, [K]', 'po, [кг/м3]', 'Cp, [Дж/кг*К]', 'μ, [Па*с]', 'λ, [Вт/м*К]'], 0, 'Параметры в 6 областях')
    render([T0, p0, po0], ['T0, [K]', 'p0, [Па]', 'po, [кг / м3]'], 0, 'Параметры торможения')
    render([thet.map(el => grad(el))], ['θ, [град]'], 1, 'Угол наколна СУ')
    render([um.map(el => grad(el)), betta.map(el => grad(el))], ['μ, [град]', 'β, [град]'], 1, 'Угол наклона маха и поворота потока')
    render([[delta * 180 / math.pi]], ['delta, [град]'], 6, 'Отклонение от оси в 5-6 областях')
    render([Tr_laminar, T_opr_laminar, Tr_turb, T_opr_turb], ['Tr (л), [K]', 'T* (л), [K]', 'Tr (т), [K]', 'T* (т), [K]'], 1, 'Определяющая температура');
    render([[Re_kr1, Re_l1], [Re_kr2, Re_l2]], ['Re_кр1 / Re_l1, [-]', 'Re_кр2 / Re_l2, [-]'], 0, 'Критическое / Обычное число Рейнольдса');
    render([[x_kr1], [x_kr2]], ['Xкр1, [м]', 'Xкр2, [м]'], 0, '1-ая и 2-ая критическая точка в смешанном ПС');
    render(ps13, ['X', 'Xф', 'δн', 'τн', 'Cfxн', 'Cfн', 'δ*н', 'δ**н', 'δ', 'τ', 'Cfx', 'Cf', 'δ*', 'δ**'], 1, 'Параметры ЛПС и ТПС 1 и 3 грани')
    render(ps24, ['X', 'Xф', 'δн', 'τн', 'Cfxн', 'Cfн', 'δ*н', 'δ**н', 'δ', 'τ', 'Cfx', 'Cf', 'δ*', 'δ**'], 1, 'Параметры ЛПС и ТПС 2 и 4 грани')
    render([[Cf1], [Cf2], [Cf3], [Cf4]], ['1 грань', '2 грань', '3 грань', '4 грань'], 1, "Средние коэфф трения ПС")
    render([[Cx], [Cy], [mz], [Cxa], [Cya], [K], [Cd]], ['Cx, [-]', 'Cy, [-]', 'mz, [-]', 'Cxa, [-]', 'Cya, [-]', 'K, [-]', 'Cd, [-]'], 0, 'АДХ')
    render(newADH(), ['Cx, [-]', 'Cy, [-]', 'mz, [-]', 'Cxa, [-]', 'Cya, [-]', 'K, [-]', 'Cd, [-]'], 0, 'АДХ с учетом трения')
}
init()
