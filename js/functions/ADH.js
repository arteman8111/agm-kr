import * as param from "../variables/const.js"

const X = (pk,p0) => (pk - p0) * param.c / 2;
const Y = (pk,p0) => (pk - p0) * param.b / 2;
const Mk = (Xk, Yk) => Yk * param.b / 4 + Xk * param.c / 4;
const Cxa = (Cx,Cy) => Cx * math.cos(param.alfa) + Cy * math.sin(param.alfa);
const Cya = (Cx,Cy) => -Cx * math.sin(param.alfa) + Cy * math.cos(param.alfa);
const Cx = (X1, X2, X3, X4) => (X1 + X2 - X3 - X4) / (param.L * math.cos(param.betta_k) * param.q_inf);
const Cy = (Y1, Y2, Y3, Y4) => (-Y1 + Y2 + Y3 - Y4) / (param.L * math.cos(param.betta_k) * param.q_inf);
const mz = (M1,M2,M3,M4) => (M1 + M2 + M3 + M4) / (param.q_inf * math.pow(param.L * math.cos(param.betta_k), 2))
const Mz = (mz) => mz * param.q_inf * math.pow(param.b,2)

export {
    X,
    Y,
    Mk,
    Cxa,
    Cya,
    Cx,
    Cy,
    mz,
    Mz
}