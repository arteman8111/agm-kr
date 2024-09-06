// Скорость v
export function v_default(Mya, a1) {
    return Mya * a1;
}

// Соотношение для теории косого скачка уплотнения
export function ksu(betta, M, k = 1.4) {
    // k - кэфф адиабаты
    return function (thet) {
        const a = 0.5 * (k + 1) * math.pow(M * math.sin(thet), 2);
        const b = 1 + 0.5 * (k - 1) * math.pow(M * math.sin(thet), 2);
        return math.tan(thet) / math.tan(thet - betta) - a / b;
    }
}
