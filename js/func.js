const k = 1.4; // k - кэфф адиабаты
const fi = 0.1; // коэффициент [-]
const n = 0.76; // коэффициент [-]
const psi = 0.85; // коэффициент [-]
const cp0 = 1004.7; // удельная теплоемкость воздуха [Дж/кг/К]
const u0 = 1.79 * math.pow(10, -5); // динамическая вязкость [Па*с]
const lymbda0 = 0.0232; // теплопроводность воздуха [Вт/м/К]

// Скорость v
export function v_default(Mya, a1) {
    return Mya * a1;
}

// Соотношение для теории косого скачка уплотнения
export function ksu(betta, M) {
    return function (thet) {
        const a = 0.5 * (k + 1) * math.pow(M * math.sin(thet), 2);
        const b = 1 + 0.5 * (k - 1) * math.pow(M * math.sin(thet), 2);
        return math.tan(thet) / math.tan(thet - betta) - a / b;
    }
}

// Критическое число рейнольдса через аппроксимацию
export const criticalReynoldsNumber = (Tr, M, T_st) => {
    const Tst_otn = T_st / Tr; // Относительная температура стенки;
    const x = (Tst_otn - 1) / (M * M);
    const y = 1 - 16 * x - 412.5 * math.pow(x, 2) - 35000 * math.pow(x, 3) - 375000 * math.pow(x, 4);
    return 5 * math.pow(10, 6) * y
}

// Расчет определяющей температуры и температуры восстановления 1-4 грани
export const getOprTemp = (flow, T, M, T_st) => {
    const Cp_opr = T => 1004.7 * Math.pow(T / 288.15, 0.1);
    const u_opr = T => 1.79e-5 * Math.pow(T / 288.15, 0.76);
    const lymbda_opr = T => 0.0232 * Math.pow(T / 261, 0.86);

    const r_opr = (Cp, u, lymbda, extent) => Math.pow(Cp * u / lymbda, extent);
    const Tr_opr = (T, r, M) => T * (1 + r * (k - 1) * M * M / 2);
    const T_opr = (Tr, T) => (T_st + T) / 2 + 0.22 * (Tr - T);

    let T_z = (T_st + T) / 2;
    let T_z_previous;
    let r_init = flow === 'turb' ? 0.88 : 0.83;
    let extent = flow === 'turb' ? 1 / 3 : 1 / 2;

    while (Math.abs(T_z - (T_z_previous ?? 0)) > (T_z * 0.01)) {
        T_z_previous = T_z;
        const Tr = Tr_opr(T, r_init, M);
        T_z = T_opr(Tr, T);
        r_init = r_opr(Cp_opr(T_z), u_opr(T_z), lymbda_opr(T_z), extent);
    }

    return [T_z, Tr_opr(T, r_init, M)];
}

// Удельная теплоемкость / динамическая вязкость / теплопроводность воздуха
export const cp_koeff = T => cp0 * math.pow(T / 288.15, fi)
export const u_koeff = T => u0 * math.pow(T / 288.15, n)
export const lymbda_koeff = T => lymbda0 * math.pow(T / 261, psi)

// Газодинамические функции
export const pi = M => math.pow(1 + math.pow(M, 2) * (k - 1) / 2, k / (k - 1)); // для отношения pi = p/p0
export const eps = M => math.pow(1 + math.pow(M, 2) * (k - 1) / 2, 1 / (k - 1)); // для отношения eps = po/po0
export const tau = M => (1 + math.pow(M, 2) * (k - 1) / 2); // для отношения tau = T/T0

// Параметры торможения
export const p_0 = (M, p_prev) => p_prev * pi(M)
export const po_0 = (M, po_prev) => po_prev * eps(M)
export const T_0 = (M, T_prev) => T_prev * tau(M)