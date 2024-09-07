import {u_koeff } from "./func.js";
import { criticalReynoldsNumber } from "./func.js";

const n = 0.76; // коэффициент [-]
// Несжимаемый ЛПС
export const Rek = (v, po, x, u) => v * po * x / u;
export const delta_l_ns = (x, Re) => 4.64 * x / math.sqrt(Re);
export const tau_l_ns = (po, v, Re) => 0.323 * po * math.pow(v, 2) / math.sqrt(Re);
export const Cfx_l_ns = Re => 0.646 / math.sqrt(Re);
export const Cf_l_ns = Re => 1.292 / math.sqrt(Re);
export const delta_l_ns_zv = (Re, x) => 1.74 * x / math.sqrt(Re);
export const delta_l_ns_zv_zv = (Re, x) => 0.646 * x / math.sqrt(Re);

// Сжимаемый ЛПС
export const delta_l = (delta, T_opr, T) => delta * math.pow(T_opr / T, (n + 1) / 2);
export const tau_l = (tau, T_opr, T) => tau * math.pow(T_opr / T, (n - 1) / 2);
export const Cfx_l = (Cfx, T_opr, T) => Cfx * math.pow(T_opr / T, (n - 1) / 2);
export const Cf_l = (Cf, T_opr, T) => Cf * math.pow(T_opr / T, (n - 1) / 2);
export const delta_l_zv = (delta_ns_zv, delta_l, delta_l_ns) => delta_ns_zv * delta_l / delta_l_ns;
export const delta_l_zv_zv = (delta_ns_zv_zv, Cfx_l, Cfx_l_ns) => delta_ns_zv_zv * Cfx_l / Cfx_l_ns;

// Расстояние перехода для ПС
export const delta_x = (delta, v, po, u) => math.pow(math.pow(delta, 5) * v * po / (u * math.pow(0.37, 5)), 1 / 4)

// Несжимаемый ТПС
export const delta_t_ns = (x, Re) => 0.37 * x / math.pow(Re, 1 / 5);
export const tau_t_ns = (po, v, Re) => 0.0255 * po * math.pow(v, 2) / math.pow(Re, 1 / 5);
export const Cfx_t_ns = Re => 0.0578 / math.pow(Re, 1 / 5);
export const Cf_t_ns = Re => 0.074 / math.pow(Re, 1 / 5);
export const delta_t_ns_zv = delta_t_ns => delta_t_ns / 8;
export const delta_t_ns_zv_zv = delta_t_ns => 7 * delta_t_ns / 72;

// Сжимаемый ТПС
export const delta_t = (delta, T_opr, T) => delta * math.pow(T_opr / T, (n + 1) / 5);
export const tau_t = (tau, T_opr, T) => tau * math.pow(T_opr / T, (n - 4) / 5);
export const Cfx_t = (Cfx, T_opr, T) => Cfx * math.pow(T_opr / T, (n - 4) / 5);
export const Cf_t = (Cf, T_opr, T) => Cf * math.pow(T_opr / T, (n - 4) / 5);
export const delta_t_zv = (delta_t_ns_zv, delta_t, delta_t_ns) => delta_t_ns_zv * delta_t / delta_t_ns;
export const delta_t_zv_zv = (delta_t_ns_zv_zv, Cfx_t, Cfx_t_ns) => delta_t_ns_zv_zv * Cfx_t / Cfx_t_ns;