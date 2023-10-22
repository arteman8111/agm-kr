import { pi, eps, tau } from "./gasodynamic.js";

const p_0 = (M, p_prev) => p_prev * pi(M)
const po_0 = (M, po_prev) => po_prev * eps(M)
const T_0 = (M, T_prev) => T_prev * tau(M)

export {
    p_0,
    po_0,
    T_0
}