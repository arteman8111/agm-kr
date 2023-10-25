import { cp0, u0, lymbda0 } from "../variables/const.js";
const fi = 0.1;
const n = 0.76;
const psi = 0.85;


const cp_koeff = T => cp0 * math.pow(T / 288.15, fi)
const u_koeff = T => u0 * math.pow(T / 288.15, n)
const lymbda_koeff = T => lymbda0 * math.pow(T / 261, psi)

export {
    cp_koeff,
    u_koeff,
    lymbda_koeff
}