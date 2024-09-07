const k = 1.4 // коэфф адиабаты

export const w = (M) => math.sqrt((k + 1) / (k - 1)) * math.atan(math.sqrt((M * M - 1) * (k - 1) / (k + 1))) - math.atan(math.sqrt(M * M - 1))
export const w_next = (w, betta_k) => w + 2 * betta_k
export const p_ksu = (M, thet, p_prev) => p_prev * (math.pow(M, 2) * math.pow(math.sin(thet), 2) * 2 * k / (k + 1) - (k - 1) / (k + 1));
export const po_ksu = (M, thet, po_prev) => po_prev * ((0.5 * (k + 1) * math.pow(M, 2) * math.pow(math.sin(thet), 2)) / (1 + 0.5 * (k - 1) * math.pow(M, 2) * math.pow(math.sin(thet), 2)));
export const T_ksu = (T_prev, p, po_inf, p_inf, po) => T_prev * p * po_inf / (p_inf * po);
export const v_ksu = (v_prev, thet, betta) => v_prev * math.cos(thet) / math.cos(thet - betta);
export const a_ksu = (p, po) => math.sqrt(k * p / po);
export const M_ksu = (v, a) => v / a;
export const um_ksu = M => math.asin(1 / M);