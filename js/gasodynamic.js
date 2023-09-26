import { k } from "./const.js";

const pi = M => math.pow(1 + math.pow(M, 2) * (k - 1) / 2, k / (k - 1)); // для отношения pi = p/p0
const eps = M => math.pow(1 + math.pow(M, 2) * (k - 1) / 2, 1 / (k - 1)); // для отношения eps = po/po0
const tau = M => (1 + math.pow(M, 2) * (k - 1) / 2); // для отношения tau = T/T0

export {
    pi,
    eps,
    tau
}