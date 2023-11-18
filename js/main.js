function init(N) {
    // Блок отрисовки
    const render = (arr, label, index, title1 = `Заголовок`) => {
        const title = `<h2 style="margin: 25px 0 0 0; font-size: 24px; padding:10px 5px; ">${title1}</h3>`
        const app = document.querySelector('#app')
        const box = document.createElement('div')
        const table = document.createElement('table')
        table.className = "data"
        const header = document.createElement('tr')
        const td = (element) => `<td class="data__element">${element}</td>`
        header.insertAdjacentHTML('beforeend', td(''))
        if (index !== 2) {
            for (let i = index; i < arr[0].length + index; i++) {
                header.insertAdjacentHTML('beforeend', td(i))
            }
        }
        for (let i = 0; i < arr.length; i++) {
            const tr = document.createElement('tr');
            tr.insertAdjacentHTML('beforeend', td(label[i]))
            arr[i].forEach(element => tr.insertAdjacentHTML('beforeend', td(math.format(element, { notation: 'exponential', precision: 5 }))))
            table.insertAdjacentElement('beforeend', tr)
        }
        table.insertAdjacentElement('afterbegin', header)
        box.insertAdjacentHTML('beforeend', title)
        box.insertAdjacentElement('beforeend', table)
        app.insertAdjacentElement('beforeend', box)
    }
    const myId = (app) => {
        app.insertAdjacentHTML('afterbegin', `<section class="id">
            <h2 class="id__title">Мои ИД:</h2>
            <dl class="id__list">
                <div class="id__box">
                    <dt class="id__item">Группа:</dt>
                    <dd class="id__item">СМ3-72</dd>
                </div>
                <div class="id__box">
                    <dt class="id__item">Вариант N:</dt>
                    <dd class="id__item">${N}</dd>
                </div>
                <div class="id__box">
                    <dt class="id__item">h, [м]:</dt>
                    <dd class="id__item">${h}</dd>
                </div>
                <div class="id__box">
                    <dt class="id__item">M, [-]:</dt>
                    <dd class="id__item">${M_inf}</dd>
                </div>
                <div class="id__box">
                    <dt class="id__item">α, [град]:</dt>
                    <dd class="id__item">${alfa}</dd>
                </div>
                <div class="id__box">
                    <dt class="id__item">T, [K]:</dt>
                    <dd class="id__item">${T_st}</dd>
                </div>
                <div class="id__box">
                    <dt class="id__item">Reкр, [-]:</dt>
                    <dd class="id__item">${Re_k}</dd>
                </div>
                <div class="id__box">
                    <dt class="id__item">с, [м]:</dt>
                    <dd class="id__item">${c}</dd>
                </div>
                <div class="id__box">
                    <dt class="id__item">b, [м]:</dt>
                    <dd class="id__item">${b}</dd>
                </div>
            </dl>
            </section>`)
    }

    //  Метод секущих для решения нелинейных уравнений
    const secantMethod = (f, init) => {
        let x0 = init;

        const delta = 0.0001;
        const f_toch = (f(x0) - f(x0 - delta)) / delta;
        const x1 = x0 - f(x0) / f_toch;
        const error = 0.0001;


        let x_prev = x0;
        let x_iter = x1;
        let x_next = -1;
        let b = 0;
        function calc(x_prev, x_iter) {
            return (x_iter - x_prev) * f(x_iter) / (f(x_iter) - f(x_prev))
        }
        while (math.abs(x_next - b) > error) {
            x_next = x_iter - calc(x_prev, x_iter);
            b = x_iter;
            x_prev = x_iter;
            x_iter = x_next;
        }
        return x_next
    }

    // Скорость v
    const v_default = (M, a) => {
        return M * a
    }

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

    // Утилита
    const rad = (value) => {
        return value * math.pi / 180;
    }
    const grad = (value) => {
        return value * 180 / math.pi
    }

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
        const Cp_opr = T => 1004.7 * math.pow(T / 288.15, 0.1);
        const u_opr = T => 1.79 * math.pow(10, -5) * math.pow(T / 288.15, 0.76);
        const lymbda_opr = T => 0.0232 * math.pow(T / 261, 0.86);

        const r_opr = (Cp, u, lymbda, extent) => math.pow(Cp * u / lymbda, extent);
        const Tr_opr = (T, r, M) => T * (1 + r * (k - 1) * M * M / 2);
        const T_opr = (Tr, T) => (T_st + T) / 2 + 0.22 * (Tr - T)
        let T_z = (T_st + T) / 2;
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

    // Расчет параметров ПС 1-4 грани
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
        let x_kr1 = criticalReynoldsNumber(Tr_laminar[inc - 1], M[inc]) * u[inc] / (po[inc] * v[inc]);
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
        console.log(xarr2);
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

    // Определение параметров набегающего потока через ГОСТ 4401-81
    const atm = (h) => {
        const r = 6356767;
        const gc = 9.80665;
        const x = 1.4;
        const R = 287.05287;
        const H = r * h / (r + h);
        let bm, Tm, Hm, pm;

        if (h < -2000) {
            return
        }
        if (h < 0) {
            bm = -0.0065; Tm = 301.15; Hm = -2000; pm = 127774;
        } else if (h < 11000) {
            bm = -0.0065; Tm = 288.15; Hm = 0; pm = 101325;
        } else if (h < 20000) {
            bm = 0; Tm = 216.65; Hm = 11000; pm = 22632;
        } else if (h < 32000) {
            bm = 0.001; Tm = 216.65; Hm = 20000; pm = 5474.87;
        } else if (h < 47000) {
            bm = 0.0028; Tm = 228.65; Hm = 32000; pm = 868.014;
        } else if (h < 51000) {
            bm = 0; Tm = 270.65; Hm = 47000; pm = 110.906;
        } else if (h < 71000) {
            bm = -0.0028; Tm = 270.65; Hm = 51000; pm = 66.9384;
        } else if (h < 85000) {
            bm = -0.002; Tm = 214.65; Hm = 71000; pm = 3.95639;
        }

        const T = Tm + bm * (H - Hm);
        const p = bm ? pm * math.exp(-gc * math.log(T / Tm) / (bm * R)) : pm * math.exp(-gc * (H - Hm) / (R * T));
        const po = p / (R * T);
        const a = math.sqrt(x * R * T);
        const g = gc * math.pow(r / (r + h), 2);
        return {
            T,
            p,
            po,
            a,
            g,
            H
        }
    }

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

    // Определение параметров
    let h = (10 + 0.4 * N) * math.pow(10, 3) // геометрическая высота [м]
    let M_inf = 4.5 + 0.1 * N // число маха для набегающего потока [-]
    let alfa = 2 + 0.1 * N // угол атаки [град]
    let T_st = 373 // температура стенки [К]
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
    let atmos = atm(h); // инициализация атмосферы
    let betta = [betta_k - rad(alfa), betta_k + rad(alfa), betta_k, betta_k]; // массив углов поворота
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
    const Cya = Cy * math.cos(alfa) - Cx * math.sin(alfa);
    const Cxa = Cy * math.sin(alfa) + Cx * math.cos(alfa);
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
            [Cx_new], [Cy], [mz], [Cxa_new], [Cya_new], [K_new], [Cd]
        ]
    }

    // Блок отрисовки таблиц
    const app = document.querySelector('#app');
    myId(app);
    render([[M[0]], [p[0]], [T[0]], [po[0]], [a[0]]], ['M, [-]', 'p, [Па]', 'T, [K]', 'po, [кг/м3]', 'a, [м/с]'], 2, 'Параметры атмосферы')
    render([M, p, T, po, a, cp, u, lymbda], ['M, [-]', 'p, [Па]', 'T, [K]', 'po, [кг/м3]', 'a, [м/с]', 'Cp, [Дж/кг*К]', 'μ, [Па*с]', 'λ, [Вт/м*К]'], 0, 'Параметры в 6 областях, где 0 - набегающий поток')
    render([T0, p0, po0], ['T0, [K]', 'p0, [Па]', 'po, [кг / м3]'], 0, 'Параметры торможения')
    render([thet.map(el => grad(el))], ['θ, [град]'], 1, 'Угол наколна СУ')
    render([um.map(el => grad(el)), betta.map(el => grad(el))], ['μ, [град]', 'β, [град]'], 1, 'Угол наклона маха и поворота потока')
    render([[delta * 180 / math.pi]], ['delta, [град]'], 2, 'Отклонение от оси в 5-6 областях')
    render([Tr_laminar, T_opr_laminar, Tr_turb, T_opr_turb], ['Tr (л), [K]', 'T* (л), [K]', 'Tr (т), [K]', 'T* (т), [K]'], 1, 'Определяющая температура');
    render([[Re_kr1, Re_l1], [Re_kr2, Re_l2]], ['Re_кр1 / Re_l1, [-]', 'Re_кр2 / Re_l2, [-]'], 2, 'Критическое / Обычное число Рейнольдса');
    render([[x_kr1], [x_kr2]], ['Xкр1, [м]', 'Xкр2, [м]'], 2, '1-ая и 2-ая критическая точка в смешанном ПС');
    render(ps13, ['X', 'Xф', 'δн', 'τн', 'Cfxн', 'Cfн', 'δ*н', 'δ**н', 'δ', 'τ', 'Cfx', 'Cf', 'δ*', 'δ**'], 1, 'Параметры ЛПС и ТПС 1 и 3 грани')
    render(ps24, ['X', 'Xф', 'δн', 'τн', 'Cfxн', 'Cfн', 'δ*н', 'δ**н', 'δ', 'τ', 'Cfx', 'Cf', 'δ*', 'δ**'], 1, 'Параметры ЛПС и ТПС 2 и 4 грани')
    render([[Cf1], [Cf2], [Cf3], [Cf4]], ['1 грань', '2 грань', '3 грань', '4 грань'], 2, "Средние коэфф трения ПС")
    render([[Cx], [Cy], [mz], [Cxa], [Cya], [K], [Cd]], ['Cx, [-]', 'Cy, [-]', 'mz, [-]', 'Cxa, [-]', 'Cya, [-]', 'K, [-]', 'Cd, [-]'], 2, 'АДХ')
    render(newADH(), ['Cx, [-]', 'Cy, [-]', 'mz, [-]', 'Cxa, [-]', 'Cya, [-]', 'K, [-]', 'Cd, [-]'], 2, 'АДХ с учетом трения')
}

const exampleModal = document.getElementById('exampleModal');
const app = document.querySelector('#app')
exampleModal.addEventListener('show.bs.modal', function (event) {
    event.stopPropagation();
    const divArea = document.querySelector('.area');
    const modalBodyInput = this.querySelector('.modal-body select')
    const modalFooterButton = this.querySelector('.modal-footer button')
    // Кнопка, запускающая модальное окно
    const button = event.relatedTarget
    modalFooterButton.addEventListener('click', () => {
        const value = modalBodyInput.value;
        if (value == 15) {
            alert("Верни долг владик!")
        } else {
            button.remove();
            divArea.remove();
            init(value);
        }
    })
})