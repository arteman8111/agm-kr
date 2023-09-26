import { cp0, u0, lymbda0 } from "./const.js";
const fi = 0.1;
const n = 0.76;
const psi = 0.85;


const cp_koeff = (T, Tk) => cp0 * math.pow(T / Tk, fi)
const u_koeff = (T, Tk) => u0 * math.pow(T / Tk, n)
const lymbda_koeff = (T, Tk) => lymbda0 * math.pow(T / Tk, psi)

export {
    cp_koeff,
    u_koeff,
    lymbda_koeff
}