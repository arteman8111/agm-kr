import { InfoDisplay, Renderer } from "./render.js"
import { rad, grad, secantMethod } from "./utils.js";
import { Atmosphere4401 } from "./atmos.js";
import { criticalReynoldsNumber, ksu, v_default, getOprTemp } from "./func.js";
import {
    Rek, delta_l_ns, tau_l_ns, Cfx_l_ns, Cf_l_ns, delta_l_ns_zv, delta_l_ns_zv_zv,
    delta_l, tau_l, Cfx_l, Cf_l, delta_l_zv, delta_l_zv_zv, delta_x,
    delta_t_ns, tau_t_ns, Cfx_t_ns, Cf_t_ns, delta_t_ns_zv, delta_t_ns_zv_zv,
    delta_t, tau_t, Cfx_t, Cf_t, delta_t_zv, delta_t_zv_zv
} from './psloy.js';
import { cp_koeff, u_koeff, lymbda_koeff, pi, eps, tau, p_0, po_0, T_0 } from "./func.js";
import {w, w_next, p_ksu, po_ksu, T_ksu, v_ksu, a_ksu, M_ksu, um_ksu} from './ksu.js';

function init(temperatureWall, alfaDeg, paramC, paramB, height, countMah) {

    // Определение параметров
    const T_st = temperatureWall;
    const alfa = rad(alfaDeg);
    const c = paramC * math.pow(10, -3);
    const b = paramB * math.pow(10, -3);
    const h = height;
    const M_inf = countMah;
    const Re_k = 5 * math.pow(10, 6) // критическое число рейнольдса [-]
    const betta_k = math.atan(c / b) // угол полураствора [рад]
    const L = b / (2 * math.cos(betta_k)); // длина грани [м]
    const betta = [betta_k - alfa, betta_k + alfa, betta_k, betta_k]; // массив углов поворота
    const thet = [secantMethod(ksu(betta[0], M_inf), rad(12)), secantMethod(ksu(betta[1], M_inf), rad(12))]; // массив СУ
    const M = [M_inf]; // мах
    const atmos = new Atmosphere4401(h); // инициализация атмосферы
    const p = [atmos.p]; // давление 
    const po = [atmos.po]; // плотность
    const T = [atmos.T]; // температура
    const a = [atmos.a] // скорость звука
    const v = [v_default(M[0], a[0])]; // скорость потока
    const um = []; // угол маха
    const p0 = []; // давление торможения
    const po0 = []; // плотность торможения
    const T0 = []; // температура торможения
    const cp = []; // удельная теплоемкость 
    const u = []; // динамическая вязкость
    const lymbda = []; // теплопроводность воздуха
    const T_opr_laminar = []; // определяющая температура для ЛПС
    const Tr_laminar = []; // температура восстановления для ЛПС
    const T_opr_turb = []; // определяющая температура для ТПС
    const Tr_turb = []; // температура восстановления для ТПС

    // Число рейнольдса
    const reynoldsNumber = (po, V, u, x = L) => po * V * x / u;

    // Средний коэффициент трения для смешанного ПС
    const Cfsm = ps => (ps[11][1] * ps[1][1] + ps[11][5] * ps[1][5] - ps[11][2] * ps[1][2]) / L;

    // Средний коэффициент трения для ТПС
    const Cfturb = ps => (ps[11][8] * ps[1][8] - ps[11][6] * ps[1][6]) / L

    // Скоростной напор
    const q = (po, v) => po * v * v / 2;

    // Расчет параметров ПС 1-4 грани
    const getSloy = (x_kr, inc) => {

        const xarr = [x_kr / 2, x_kr];
        const delta_ns = [];
        const tau_ns = [];
        const Cfx_ns = [];
        const Cf_ns = [];
        const delta_ns_zv = [];
        const delta_ns_zv_zv = [];
        const delta_ = [];
        const tau = [];
        const Cfx = [];
        const Cf = [];
        const delta_zv = [];
        const delta_zv_zv = [];

        xarr.forEach(
            (x, i) => {
                const Re = Rek(v[inc], po[inc], x, u[inc]);
                delta_ns.push(delta_l_ns(x, Re));
                tau_ns.push(tau_l_ns(po[inc], v[inc], Re));
                Cfx_ns.push(Cfx_l_ns(Re));
                Cf_ns.push(Cf_l_ns(Re));
                delta_ns_zv.push(delta_l_ns_zv(Re, x));
                delta_ns_zv_zv.push(delta_l_ns_zv_zv(Re, x));

                delta_.push(delta_l(delta_ns[i], T_opr_laminar[inc - 1], T[inc]));
                tau.push(tau_l(tau_ns[i], T_opr_laminar[inc - 1], T[inc]));
                Cfx.push(Cfx_l(Cfx_ns[i], T_opr_laminar[inc - 1], T[inc]));
                Cf.push(Cf_l(Cf_ns[i], T_opr_laminar[inc - 1], T[inc]));
                delta_zv.push(delta_l_zv(delta_ns_zv[i], delta_[i], delta_ns[i]));
                delta_zv_zv.push(delta_l_zv_zv(delta_ns_zv_zv[i], Cfx[i], Cfx_ns[i]))
            }
        )
        let po_zv = po[inc] * T[inc] / T_opr_laminar[inc - 1];
        let u_zv = u_koeff(T_opr_laminar[inc - 1]);
        let x_kr1 = criticalReynoldsNumber(Tr_laminar[inc - 1], M[inc], T_st) * u[inc] / (po[inc] * v[inc]);
        let Re_krit = po_zv * v[inc] * x_kr1 / u_zv;
        let delta_xk = math.pow((4.64 * x_kr1 * math.pow(po_zv * v[inc] / u_zv, 0.2)) / (math.sqrt(Re_krit) * 0.37), 1.25)
        let b = L - x_kr;
        xarr.push(delta_xk, delta_xk + b / 3, delta_xk + 2 * b / 3, delta_xk + b);
        let xarr3_6 = [delta_xk, delta_xk + b / 3, delta_xk + 2 * b / 3, delta_xk + b]
        xarr3_6.forEach(
            (x, i) => {
                const Re = Rek(v[inc], po[inc], x, u[inc]);
                delta_ns.push(delta_t_ns(x, Re));
                tau_ns.push(tau_t_ns(po[inc], v[inc], Re));
                Cfx_ns.push(Cfx_t_ns(Re));
                Cf_ns.push(Cf_t_ns(Re));
                delta_ns_zv.push(delta_t_ns_zv(delta_ns[i + 2]));
                delta_ns_zv_zv.push(delta_t_ns_zv_zv(delta_ns[i + 2]));

                delta_.push(delta_t(delta_ns[i + 2], T_opr_turb[inc - 1], T[inc]));
                tau.push(tau_t(tau_ns[i + 2], T_opr_turb[inc - 1], T[inc]));
                Cfx.push(Cfx_t(Cfx_ns[i + 2], T_opr_turb[inc - 1], T[inc]));
                Cf.push(Cf_t(Cf_ns[i + 2], T_opr_turb[inc - 1], T[inc]));
                delta_zv.push(delta_t_zv(delta_ns_zv[i + 2], delta_[i + 2], delta_ns[i + 2]));
                delta_zv_zv.push(delta_t_zv_zv(delta_ns_zv_zv[i + 2], Cfx[i + 2], Cfx_ns[i + 2]))
            }
        )
        let po_zv_T = po[inc + 2] * T[inc + 2] / T_opr_turb[inc + 1];
        let u_zv_T = u_koeff(T_opr_turb[inc + 1]);
        let delta_xk2 = math.pow(math.pow(delta_[5], 5) * po_zv_T * v[inc + 2] / (u_zv_T * math.pow(0.37, 5)), 1 / 4);
        xarr.push(delta_xk2, L / 2 + delta_xk2, L + delta_xk2);
        let xarr7_9 = [delta_xk2, L / 2 + delta_xk2, L + delta_xk2];
        xarr7_9.forEach(
            (x, i) => {
                const Re = Rek(v[inc + 2], po[inc + 2], x, u[inc + 2]);
                delta_ns.push(delta_t_ns(x, Re));
                tau_ns.push(tau_t_ns(po[inc + 2], v[inc + 2], Re));
                Cfx_ns.push(Cfx_t_ns(Re));
                Cf_ns.push(Cf_t_ns(Re));
                delta_ns_zv.push(delta_t_ns_zv(delta_ns[i + 6]));
                delta_ns_zv_zv.push(delta_t_ns_zv_zv(delta_ns[i + 6]));

                delta_.push(delta_t(delta_ns[i + 6], T_opr_turb[inc + 1], T[inc + 2]));
                tau.push(tau_t(tau_ns[i + 6], T_opr_turb[inc + 1], T[inc + 2]));
                Cfx.push(Cfx_t(Cfx_ns[i + 6], T_opr_turb[inc + 1], T[inc + 2]));
                Cf.push(Cf_t(Cf_ns[i + 6], T_opr_turb[inc + 1], T[inc + 2]));
                delta_zv.push(delta_t_zv(delta_ns_zv[i + 6], delta_[i + 6], delta_ns[i + 6]));
                delta_zv_zv.push(delta_t_zv_zv(delta_ns_zv_zv[i + 6], Cfx[i + 6], Cfx_ns[i + 6]))
            }
        )
        let xarr2 = [
            x_kr / 2,
            x_kr,
            x_kr,
            x_kr + (L - x_kr) / 3,
            x_kr + 2 * (L - x_kr) / 3,
            L,
            0,
            L / 2,
            L
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

    // Сила продольная и нормальная
    const X = pk => (pk - p[0]) * c / 2;
    const Y = pk => (pk - p[0]) * b / 2;

    function paramPush() {
        // Параметры в областях 1-2
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
        }
        for (let i = 0; i <= 2; i++) {
            p0.push(p_0(M[i], p[i]));
            po0.push(po_0(M[i], po[i]));
            T0.push(T_0(M[i], T[i]));
        }
        for (let i = 1; i <= 2; i++) {
            const w0 = w(M[i]);
            const prandtl = (M) => w_next(w0, betta_k) - w(M);
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

        // Расчет параметров в конце ромба
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

        // Параметры 5-6 области
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

        // cp, u, lymbda
        T.forEach(Tk => {
            cp.push(cp_koeff(Tk))
            u.push((u_koeff(Tk)))
            lymbda.push(lymbda_koeff(Tk))
        });

        // Параметры торможения
        for (let i = 3; i <= 6; i++) {
            p0.push(p_0(M[i], p[i]));
            po0.push(po_0(M[i], po[i]));
            T0.push(T_0(M[i], T[i]));
        }
    }
    paramPush();

    // АДХ без трения
    const q_inf = q(po[0],v[0]);
    const X1 = X(p[1]);
    const X2 = X(p[2]);
    const X3 = X(p[3]);
    const X4 = X(p[4]);
    const Y1 = Y(p[1]);
    const Y2 = Y(p[2]);
    const Y3 = Y(p[3]);
    const Y4 = Y(p[4]);
    const Cx = (X1 + X2 - X3 - X4) / (b * q_inf);
    const Cy = (-Y1 + Y2 - Y3 + Y4) / (b * q_inf);
    const mz1 = Y1 * (b / 4) + X1 * (c / 4);
    const mz2 = -Y2 * (b / 4) - X2 * (c / 4);
    const mz3 = Y3 * (3 * b / 4) - X3 * (c / 4);
    const mz4 = -Y4 * (3 * b / 4) + X4 * (c / 4);
    const mz = (mz1 + mz2 + mz3 + mz4) / (q_inf * b * b);
    const Cxa = Cy * math.sin(alfa) + Cx * math.cos(alfa);
    const Cya = Cy * math.cos(alfa) - Cx * math.sin(alfa);
    const K = Cya / Cxa;
    const Cd = - (mz / Cy)

    // Определяющая температура и температура восстановления
    for (let i = 1; i <= 4; i++) {
        const [T_iter_laminar_opr, Tr_iter_laminar] = getOprTemp('laminar', T[i], M[i], T_st);
        const [T_iter_turb_opr, Tr_iter_turb] = getOprTemp('turb', T[i], M[i], T_st);
        T_opr_laminar.push(T_iter_laminar_opr);
        T_opr_turb.push(T_iter_turb_opr);
        Tr_turb.push(Tr_iter_turb);
        Tr_laminar.push(Tr_iter_laminar);
    }

    // Число Рейнольдса и критические точки ПС
    const Re_l1 = reynoldsNumber(po[1], v[1], u[1]);
    const Re_l2 = reynoldsNumber(po[2], v[2], u[2]);
    const Re_kr1 = criticalReynoldsNumber(Tr_laminar[0], M[1], T_st);
    const Re_kr2 = criticalReynoldsNumber(Tr_laminar[1], M[2], T_st);
    const x_kr1 = Re_kr1 * u[1] / (v[1] * po[1]);
    const x_kr2 = Re_kr2 * u[2] / (v[2] * po[2]);

    // Расчет ПС
    const ps13 = getSloy(x_kr1, 1);
    const ps24 = getSloy(x_kr2, 2);

    // Расчет среднего коэфф.трения
    const Cf1 = Cfsm(ps13);
    const Cf2 = Cfsm(ps24);
    const Cf3 = Cfturb(ps13);
    const Cf4 = Cfturb(ps24);

    // Перерасчет АДХ с учетом трения
    const newADH = () => {
        const Cxf = (Cf1 * q(po[1], v[1]) + Cf2 * q(po[2], v[2]) + Cf3 * q(po[3], v[3]) + Cf4 * q(po[4], v[4])) / (2 * q(po[0], v[0]));
        const Cx_new = Cx + Cxf;
        const Cya_new = Cy * math.cos(alfa) - Cx_new * math.sin(alfa);
        const Cxa_new = Cy * math.sin(alfa) + Cx_new * math.cos(alfa);
        const K_new = Cya_new / Cxa_new;

        return [
            [Cx_new], [Cxa_new], [Cya_new], [K_new]
        ]
    }

    // Блок отрисовки таблиц
    const getIdTable = new InfoDisplay('#app');
    const getResultTable = new Renderer('#app');
    getIdTable.render(h, M_inf, alfa, T_st, Re_k, c, b);
    getResultTable.render([[M[0]], [p[0]], [T[0]], [po[0]], [a[0]]], ['M, [-]', 'p, [Па]', 'T, [K]', 'po, [кг/м3]', 'a, [м/с]'], 2, 'Параметры атмосферы')
    getResultTable.render([M, p, T, po, a, cp, u, lymbda], ['M, [-]', 'p, [Па]', 'T, [K]', 'po, [кг/м3]', 'a, [м/с]', 'Cp, [Дж/кг*К]', 'μ, [Па*с]', 'λ, [Вт/м*К]'], 0, 'Параметры в 6 областях, где 0 - набегающий поток')
    getResultTable.render([T0, p0, po0], ['T0, [K]', 'p0, [Па]', 'po, [кг / м3]'], 0, 'Параметры торможения')
    getResultTable.render([thet.map(el => grad(el))], ['θ, [град]'], 1, 'Угол наколна СУ')
    getResultTable.render([um.map(el => grad(el)), betta.map(el => grad(el))], ['μ, [град]', 'β, [град]'], 1, 'Угол наклона маха и поворота потока')
    getResultTable.render([Tr_laminar, T_opr_laminar, Tr_turb, T_opr_turb], ['Tr (л), [K]', 'T* (л), [K]', 'Tr (т), [K]', 'T* (т), [K]'], 1, 'Определяющая температура');
    getResultTable.render([[Re_kr1, Re_l1], [Re_kr2, Re_l2]], ['Re_кр1 / Re_l1, [-]', 'Re_кр2 / Re_l2, [-]'], 2, 'Критическое / Обычное число Рейнольдса');
    getResultTable.render([[x_kr1], [x_kr2]], ['Xкр1, [м]', 'Xкр2, [м]'], 2, '1-ая и 2-ая критическая точка в смешанном ПС');
    getResultTable.render(ps13, ['X', 'Xф', 'δн', 'τн', 'Cfxн', 'Cfн', 'δ*н', 'δ**н', 'δ', 'τ', 'Cfx', 'Cf', 'δ*', 'δ**'], 1, 'Параметры ЛПС и ТПС 1 и 3 грани')
    getResultTable.render(ps24, ['X', 'Xф', 'δн', 'τн', 'Cfxн', 'Cfн', 'δ*н', 'δ**н', 'δ', 'τ', 'Cfx', 'Cf', 'δ*', 'δ**'], 1, 'Параметры ЛПС и ТПС 2 и 4 грани')
    getResultTable.render([[Cf1], [Cf2], [Cf3], [Cf4]], ['1 грань', '2 грань', '3 грань', '4 грань'], 2, "Средние коэфф трения ПС")
    getResultTable.render([[Cx], [Cy], [mz], [Cxa], [Cya], [K], [Cd]], ['Cx, [-]', 'Cy, [-]', 'mz, [-]', 'Cxa, [-]', 'Cya, [-]', 'K, [-]', 'Cd, [-]'], 2, 'АДХ')
    getResultTable.render(newADH(), ['Cx, [-]', 'Cxa, [-]', 'Cya, [-]', 'K, [-]'], 2, 'АДХ с учетом трения')
}

const exampleModal = document.getElementById('exampleModal');
exampleModal.addEventListener('show.bs.modal', function (event) {
    event.stopPropagation();
    const divArea = document.querySelector('.area');

    const modal_T_st = this.querySelector('#T_st');
    const modal_alfa = this.querySelector('#alfa');
    const modal_c = this.querySelector('#c');
    const modal_b = this.querySelector('#b');
    const modal_M_inf = this.querySelector('#M_inf');
    const modal_h = this.querySelector('#h');
    const modalFooterButton = this.querySelector('.modal-footer button')

    // Кнопка, запускающая модальное окно
    const button = event.relatedTarget

    modalFooterButton.addEventListener('click', () => {
        console.log("Tst", +modal_T_st.value);
        console.log("alfa", +modal_alfa.value);
        console.log("c", +modal_c.value);
        console.log("b", +modal_b.value);
        console.log("Minf", +modal_M_inf.value);
        console.log("h", +modal_h.value);

        const T_st = +modal_T_st.value;
        const alfa = +modal_alfa.value;
        const c = +modal_c.value;
        const b = +modal_b.value;
        const h = +modal_h.value;
        const M_inf = +modal_M_inf.value;

        button.remove();
        divArea.remove();
        init(T_st, alfa, c, b, h, M_inf);
    })
})