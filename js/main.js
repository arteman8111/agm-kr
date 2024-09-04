import { InfoDisplay, Renderer } from "./render.js"
import { rad, grad, secantMethod } from "./utils.js";
import { Atmosphere4401 } from "./atmos.js";

function init(N, group) {

    // Скорость v
    const v_default = (M, a) => M * a

    // Соотношение для теории косого скачка уплотнения
    const ksu = (betta, M) => {
        return function (thet) {
            const a = 0.5 * (k + 1) * math.pow(M * math.sin(thet), 2);
            const b = 1 + 0.5 * (k - 1) * math.pow(M * math.sin(thet), 2);
            return math.tan(thet) / math.tan(thet - betta) - a / b;
        }
    }

    // Критическое число рейнольдса через аппроксимацию
    const criticalReynoldsNumber = (Tr, M) => {
        const Tst_otn = T_st / Tr; // Относительная температура стенки;
        const x = (Tst_otn - 1) / (M * M);
        const y = 1 - 16 * x - 412.5 * math.pow(x, 2) - 35000 * math.pow(x, 3) - 375000 * math.pow(x, 4);
        return 5 * math.pow(10, 6) * y
    }

    // Число рейнольдса
    const reynoldsNumber = (po, V, u, x = L) => po * V * x / u;

    // w
    const w = (M) => math.sqrt((k + 1) / (k - 1)) * math.atan(math.sqrt((M * M - 1) * (k - 1) / (k + 1))) - math.atan(math.sqrt(M * M - 1))
    const w_next = (w) => w + 2 * betta_k

    // Несжимаемый ЛПС
    const Rek = (v, po, x, u) => v * po * x / u;
    const delta_l_ns = (x, Re) => 4.64 * x / math.sqrt(Re);
    const tau_l_ns = (po, v, Re) => 0.323 * po * math.pow(v, 2) / math.sqrt(Re);
    const Cfx_l_ns = Re => 0.646 / math.sqrt(Re);
    const Cf_l_ns = Re => 1.292 / math.sqrt(Re);
    const delta_l_ns_zv = (Re, x) => 1.74 * x / math.sqrt(Re);
    const delta_l_ns_zv_zv = (Re, x) => 0.646 * x / math.sqrt(Re);

    // Сжимаемый ЛПС
    const delta_l = (delta, T_opr, T) => delta * math.pow(T_opr / T, (n + 1) / 2);
    const tau_l = (tau, T_opr, T) => tau * math.pow(T_opr / T, (n - 1) / 2);
    const Cfx_l = (Cfx, T_opr, T) => Cfx * math.pow(T_opr / T, (n - 1) / 2);
    const Cf_l = (Cf, T_opr, T) => Cf * math.pow(T_opr / T, (n - 1) / 2);
    const delta_l_zv = (delta_ns_zv, delta_l, delta_l_ns) => delta_ns_zv * delta_l / delta_l_ns;
    const delta_l_zv_zv = (delta_ns_zv_zv, Cfx_l, Cfx_l_ns) => delta_ns_zv_zv * Cfx_l / Cfx_l_ns;

    // Расстояние перехода для ПС
    const delta_x = (delta, v, po, u) => math.pow(math.pow(delta, 5) * v * po / (u * math.pow(0.37, 5)), 1 / 4)

    // Несжимаемый ТПС
    const delta_t_ns = (x, Re) => 0.37 * x / math.pow(Re, 1 / 5);
    const tau_t_ns = (po, v, Re) => 0.0255 * po * math.pow(v, 2) / math.pow(Re, 1 / 5);
    const Cfx_t_ns = Re => 0.0578 / math.pow(Re, 1 / 5);
    const Cf_t_ns = Re => 0.074 / math.pow(Re, 1 / 5);
    const delta_t_ns_zv = delta_t_ns => delta_t_ns / 8;
    const delta_t_ns_zv_zv = delta_t_ns => 7 * delta_t_ns / 72;

    // Сжимаемый ТПС
    const delta_t = (delta, T_opr, T) => delta * math.pow(T_opr / T, (n + 1) / 5);
    const tau_t = (tau, T_opr, T) => tau * math.pow(T_opr / T, (n - 4) / 5);
    const Cfx_t = (Cfx, T_opr, T) => Cfx * math.pow(T_opr / T, (n - 4) / 5);
    const Cf_t = (Cf, T_opr, T) => Cf * math.pow(T_opr / T, (n - 4) / 5);
    const delta_t_zv = (delta_t_ns_zv, delta_t, delta_t_ns) => delta_t_ns_zv * delta_t / delta_t_ns;
    const delta_t_zv_zv = (delta_t_ns_zv_zv, Cfx_t, Cfx_t_ns) => delta_t_ns_zv_zv * Cfx_t / Cfx_t_ns;

    const p_ksu = (M, thet, p_prev) => p_prev * (math.pow(M, 2) * math.pow(math.sin(thet), 2) * 2 * k / (k + 1) - (k - 1) / (k + 1));
    const po_ksu = (M, thet, po_prev) => po_prev * ((0.5 * (k + 1) * math.pow(M, 2) * math.pow(math.sin(thet), 2)) / (1 + 0.5 * (k - 1) * math.pow(M, 2) * math.pow(math.sin(thet), 2)));
    const T_ksu = (T_prev, p, po_inf, p_inf, po) => T_prev * p * po_inf / (p_inf * po);
    const v_ksu = (v_prev, thet, betta) => v_prev * math.cos(thet) / math.cos(thet - betta);
    const a_ksu = (p, po) => math.sqrt(k * p / po);
    const M_ksu = (v, a) => v / a;
    const um_ksu = M => math.asin(1 / M);

    // Параметры торможения
    const p_0 = (M, p_prev) => p_prev * pi(M)
    const po_0 = (M, po_prev) => po_prev * eps(M)
    const T_0 = (M, T_prev) => T_prev * tau(M)

    // Удельная теплоемкость / динамическая вязкость / теплопроводность воздуха
    const cp_koeff = T => cp0 * math.pow(T / 288.15, fi)
    const u_koeff = T => u0 * math.pow(T / 288.15, n)
    const lymbda_koeff = T => lymbda0 * math.pow(T / 261, psi)

    // Газодинамические функции
    const pi = M => math.pow(1 + math.pow(M, 2) * (k - 1) / 2, k / (k - 1)); // для отношения pi = p/p0
    const eps = M => math.pow(1 + math.pow(M, 2) * (k - 1) / 2, 1 / (k - 1)); // для отношения eps = po/po0
    const tau = M => (1 + math.pow(M, 2) * (k - 1) / 2); // для отношения tau = T/T0

    // Средний коэффициент трения для смешанного ПС
    const Cfsm = (ps) => {
        return (ps[11][1] * ps[1][1] + ps[11][5] * ps[1][5] - ps[11][2] * ps[1][2]) / L;
    }

    // Средний коэффициент трения для ТПС
    const Cfturb = (ps) => {
        return (ps[11][8] * ps[1][8] - ps[11][6] * ps[1][6]) / L
    }

    // Скоростной напор
    const q = (po, v) => {
        return po * v * v / 2;
    }

    // Расчет определяющей температуры и температуры восстановления 1-4 грани
    const getOprTemp = (flow, T, M) => {
        const Cp_opr = T => 1004.7 * Math.pow(T / 288.15, 0.1);
        const u_opr = T => 1.79e-5 * Math.pow(T / 288.15, 0.76);
        const lymbda_opr = T => 0.0232 * Math.pow(T / 261, 0.86);
    
        const r_opr = (Cp, u, lymbda, extent) => Math.pow(Cp * u / lymbda, extent);
        const Tr_opr = (T, r, M) => T * (1 + r * (k - 1) * M * M / 2);
        const T_opr = (Tr, T) => (T_st + T) / 2 + 0.22 * (Tr - T);
    
        let T_z = (T_st + T) / 2;
        let T_z_previous;
        const r_init = flow === 'turb' ? 0.88 : 0.83;
        const extent = flow === 'turb' ? 1 / 3 : 1 / 2;
    
        while (Math.abs(T_z - (T_z_previous ?? 0)) > (T_z * 0.01)) {
            T_z_previous = T_z;
            const Tr = Tr_opr(T, r_init, M);
            T_z = T_opr(Tr, T);
            r_init = r_opr(Cp_opr(T_z), u_opr(T_z), lymbda_opr(T_z), extent);
        }
    
        return [T_z, Tr_opr(T, r_init, M)];
    }
    

    // Расчет параметров ПС 1-4 грани
    const getSloy = (x_kr, inc) => {
        const xarr = [x_kr / 2, x_kr];
        const results = {
            delta_ns: [], tau_ns: [], Cfx_ns: [], Cf_ns: [], delta_ns_zv: [], delta_ns_zv_zv: [],
            delta_: [], tau: [], Cfx: [], Cf: [], delta_zv: [], delta_zv_zv: []
        };
    
        const computeValues = (x, inc, idxOffset = 0, isTurb = false) => {
            const Re = Rek(v[inc], po[inc], x, u[inc]);
            const deltaFunc = isTurb ? delta_t_ns : delta_l_ns;
            const tauFunc = isTurb ? tau_t_ns : tau_l_ns;
            const CfxFunc = isTurb ? Cfx_t_ns : Cfx_l_ns;
            const CfFunc = isTurb ? Cf_t_ns : Cf_l_ns;
            const delta_zv_Func = isTurb ? delta_t_ns_zv : delta_l_ns_zv;
            const delta_zv_zv_Func = isTurb ? delta_t_ns_zv_zv : delta_l_ns_zv_zv;
            const T_opr = isTurb ? T_opr_turb[inc - 1] : T_opr_laminar[inc - 1];
            results.delta_ns.push(deltaFunc(x, Re));
            results.tau_ns.push(tauFunc(po[inc], v[inc], Re));
            results.Cfx_ns.push(CfxFunc(Re));
            results.Cf_ns.push(CfFunc(Re));
            results.delta_ns_zv.push(delta_zv_Func(results.delta_ns[idxOffset]));
            results.delta_ns_zv_zv.push(delta_zv_zv_Func(results.delta_ns[idxOffset]));
            results.delta_.push(delta_t(results.delta_ns[idxOffset], T_opr, T[inc]));
            results.tau.push(tau_t(results.tau_ns[idxOffset], T_opr, T[inc]));
            results.Cfx.push(Cfx_t(results.Cfx_ns[idxOffset], T_opr, T[inc]));
            results.Cf.push(Cf_t(results.Cf_ns[idxOffset], T_opr, T[inc]));
            results.delta_zv.push(delta_t_zv(results.delta_ns_zv[idxOffset], results.delta_[idxOffset], results.delta_ns[idxOffset]));
            results.delta_zv_zv.push(delta_t_zv_zv(results.delta_ns_zv_zv[idxOffset], results.Cfx[idxOffset], results.Cfx_ns[idxOffset]));
        };
    
        xarr.forEach((x, i) => computeValues(x, inc, i));
        const po_zv = po[inc] * T[inc] / T_opr_laminar[inc - 1];
        const u_zv = u_koeff(T_opr_laminar[inc - 1]);
        const x_kr1 = criticalReynoldsNumber(Tr_laminar[inc - 1], M[inc]) * u[inc] / (po[inc] * v[inc]);
        const Re_krit = po_zv * v[inc] * x_kr1 / u_zv;
        const delta_xk = Math.pow((4.64 * x_kr1 * Math.pow(po_zv * v[inc] / u_zv, 0.2)) / (Math.sqrt(Re_krit) * 0.37), 1.25);
        const xarr3_6 = [delta_xk, delta_xk + (L - x_kr) / 3, delta_xk + 2 * (L - x_kr) / 3, delta_xk + (L - x_kr)];
        xarr3_6.forEach((x, i) => computeValues(x, inc, i + 2, true));
        const po_zv_T = po[inc + 2] * T[inc + 2] / T_opr_turb[inc + 1];
        const u_zv_T = u_koeff(T_opr_turb[inc + 1]);
        const delta_xk2 = Math.pow(Math.pow(results.delta_[5], 5) * po_zv_T * v[inc + 2] / (u_zv_T * Math.pow(0.37, 5)), 1 / 4);
        const xarr7_9 = [delta_xk2, L / 2 + delta_xk2, L + delta_xk2];
        xarr7_9.forEach((x, i) => computeValues(x, inc + 2, i + 6, true));
    
        const xarr2 = [
            x_kr / 2, x_kr, x_kr, x_kr + (L - x_kr) / 3, x_kr + 2 * (L - x_kr) / 3, L, 0, L / 2, L
        ];
    
        return [xarr2, [...xarr, ...xarr3_6, ...xarr7_9], ...Object.values(results)];
    };
    

    // Сила продольная и нормальная
    const X = pk => (pk - p[0]) * c / 2;
    const Y = pk => (pk - p[0]) * b / 2;

    let c, b;
    const cb_get = (N) => {
        if (N <= 6) {
            c = 80 * math.pow(10, -3);
            b = 1000 * math.pow(10, -3);
        } else if (N <= 12) {
            c = 90 * math.pow(10, -3);
            b = 1100 * math.pow(10, -3);
        } else if (N <= 18) {
            c = 100 * math.pow(10, -3);
            b = 1200 * math.pow(10, -3);
        } else if (N <= 24) {
            c = 110 * math.pow(10, -3);
            b = 1300 * math.pow(10, -3);
        }
    }
    cb_get(N);
    let h, M_inf, T_st

    const group_get = (groupN) => {
        if (groupN === 2){
            h = (10 + 0.4 * N) * math.pow(10, 3) // геометрическая высота [м]
            M_inf = 4.5 + 0.1 * N // число маха для набегающего потока [-]
            T_st = 373 // температура стенки [К]
        } else if (groupN === 1) {
            h = (10 + 0.5 * N) * math.pow(10, 3) // геометрическая высота [м]
            M_inf = 4 + 0.15 * N // число маха для набегающего потока [-]
            T_st = 405 // температура стенки [К]
        }
    }
    group_get(group)

    // Определение параметров
    let alfa = rad(2 + 0.1 * N) // угол атаки [град]
    let Re_k = 5 * math.pow(10, 6) // критическое число рейнольдса [-]
    let cp0 = 1004.7; // удельная теплоемкость воздуха [Дж/кг/К]
    let u0 = 1.79 * math.pow(10, -5); // динамическая вязкость [Па*с]
    let lymbda0 = 0.0232; // теплопроводность воздуха [Вт/м/К]
    let k = 1.4 // адиабата коэффициент [-]
    let betta_k = math.atan(c / b) // угол полураствора [рад]
    let L = b / (2 * math.cos(betta_k)); // длина грани [м]
    let fi = 0.1; // коэффициент [-]
    let n = 0.76; // коэффициент [-]
    let psi = 0.85; // коэффициент [-]
    let atmos = new Atmosphere4401(h); // инициализация атмосферы
    let betta = [betta_k - alfa, betta_k + alfa, betta_k, betta_k]; // массив углов поворота
    let thet = [secantMethod(ksu(betta[0], M_inf), rad(12)), secantMethod(ksu(betta[1], M_inf), rad(12))]; // массив СУ
    let M = [M_inf]; // мах
    let p = [atmos.p]; // давление 
    let po = [atmos.po]; // плотность
    let T = [atmos.T]; // температура
    let a = [atmos.a] // скорость звука
    let v = [v_default(M[0], a[0])]; // скорость потока
    let um = []; // угол маха
    let p0 = []; // давление торможения
    let po0 = []; // плотность торможения
    let T0 = []; // температура торможения
    let cp = []; // удельная теплоемкость 
    let u = []; // динамическая вязкость
    let lymbda = []; // теплопроводность воздуха
    let T_opr_laminar = []; // определяющая температура для ЛПС
    let Tr_laminar = []; // температура восстановления для ЛПС
    let T_opr_turb = []; // определяющая температура для ТПС
    let Tr_turb = []; // температура восстановления для ТПС
    let delta; // отклонение от оси в 5-6 областях

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
    // Расчет параметров в конце ромба
    (function posled() {
        delta = rad(0);
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

    // АДХ без трения
    const q_inf = po[0] * v[0] * v[0] / 2; // Скоростной напор набегающего потока
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
        const [T_iter_laminar_opr, Tr_iter_laminar] = getOprTemp('laminar', T[i], M[i]);
        const [T_iter_turb_opr, Tr_iter_turb] = getOprTemp('turb', T[i], M[i]);
        T_opr_laminar.push(T_iter_laminar_opr);
        T_opr_turb.push(T_iter_turb_opr);
        Tr_turb.push(Tr_iter_turb);
        Tr_laminar.push(Tr_iter_laminar);
    }

    // Число Рейнольдса и критические точки ПС
    const Re_l1 = reynoldsNumber(po[1], v[1], u[1]);
    const Re_l2 = reynoldsNumber(po[2], v[2], u[2]);
    const Re_kr1 = criticalReynoldsNumber(Tr_laminar[0], M[1]);
    const Re_kr2 = criticalReynoldsNumber(Tr_laminar[1], M[2]);
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
    getIdTable.render(group, N, h, M_inf, alfa, T_st, Re_k, c, b);
    getResultTable.render([[M[0]], [p[0]], [T[0]], [po[0]], [a[0]]], ['M, [-]', 'p, [Па]', 'T, [K]', 'po, [кг/м3]', 'a, [м/с]'], 2, 'Параметры атмосферы')
    getResultTable.render([M, p, T, po, a, cp, u, lymbda], ['M, [-]', 'p, [Па]', 'T, [K]', 'po, [кг/м3]', 'a, [м/с]', 'Cp, [Дж/кг*К]', 'μ, [Па*с]', 'λ, [Вт/м*К]'], 0, 'Параметры в 6 областях, где 0 - набегающий поток')
    getResultTable.render([T0, p0, po0], ['T0, [K]', 'p0, [Па]', 'po, [кг / м3]'], 0, 'Параметры торможения')
    getResultTable.render([thet.map(el => grad(el))], ['θ, [град]'], 1, 'Угол наколна СУ')
    getResultTable.render([um.map(el => grad(el)), betta.map(el => grad(el))], ['μ, [град]', 'β, [град]'], 1, 'Угол наклона маха и поворота потока')
    getResultTable.render([[delta * 180 / math.pi]], ['delta, [град]'], 2, 'Отклонение от оси в 5-6 областях')
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
const app = document.querySelector('#app')
exampleModal.addEventListener('show.bs.modal', function (event) {
    event.stopPropagation();
    const divArea = document.querySelector('.area');
    const modalBodyInput = this.querySelector('#group-number')
    const modalBodyInputGroup = this.querySelector('#group')
    const modalFooterButton = this.querySelector('.modal-footer button')
    // Кнопка, запускающая модальное окно
    const button = event.relatedTarget
    modalFooterButton.addEventListener('click', () => {
        const value = +modalBodyInput.value;
        const group = +modalBodyInputGroup.value;
        button.remove();
        divArea.remove();
        init(value, group);
    })
})