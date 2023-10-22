import { k } from "../variables/const.js"

const p_ksu = (M, thet, p_prev) => p_prev * (math.pow(M, 2) * math.pow(math.sin(thet), 2) * 2 * k / (k + 1) - (k - 1) / (k + 1));
const po_ksu = (M, thet, po_prev) => po_prev * ((0.5 * (k + 1) * math.pow(M, 2) * math.pow(math.sin(thet), 2)) / (1 + 0.5 * (k - 1) * math.pow(M, 2) * math.pow(math.sin(thet), 2)));
const T_ksu = (T_prev, p, po_inf, p_inf, po) => T_prev * p * po_inf / (p_inf * po);
const v_ksu = (v_prev, thet, betta) => v_prev * math.cos(thet) / math.cos(thet - betta);
const a_ksu = (p, po) => math.sqrt(k * p / po);
const M_ksu = (v, a) => v / a;
const um_ksu = M => math.asin(1 / M);

export { 
    p_ksu, 
    po_ksu, 
    T_ksu, 
    v_ksu, 
    a_ksu, 
    M_ksu, 
    um_ksu 
}