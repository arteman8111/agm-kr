import {p, po, v} from "./6Oblastey.js"
import * as param from "../variables/const.js"

const q_inf = po[0] * v[0] * v[0] / 2; // q_inf
const X = pk => (pk - p[0]) * param.c / 2;
const Y = pk => (pk - p[0]) * param.b / 2;

const X1 = X(p[1]);
const X2 = X(p[2]);
const X3 = X(p[3]);
const X4 = X(p[4]);

const Y1 = Y(p[1]);
const Y2 = Y(p[2]);
const Y3 = Y(p[3]);
const Y4 = Y(p[4]);

const Cx = (X1 + X2 - X3 - X4) / (param.b * q_inf);
const Cy = (-Y1 + Y2 - Y3 + Y4) / (param.b * q_inf);
const mz1 = Y1 * (param.b / 4) + X1 * (param.c / 4);
const mz2 = -Y2 * (param.b / 4) - X2 * (param.c / 4);
const mz3 = Y3 * (3 * param.b / 4) - X3 * (param.c / 4);
const mz4 = -Y4 * (3 * param.b / 4) + X4 * (param.c / 4);
const mz = (mz1 + mz2 + mz3 + mz4) / (q_inf * param.b * param.b);
const Cya = Cy * math.cos(param.alfa) - Cx * math.sin(param.alfa);
const Cxa = Cy * math.sin(param.alfa) + Cx * math.cos(param.alfa);
const K = Cya / Cxa;
const Cd = - (mz / Cy)

export {
    Cx, Cy, mz, Cxa, Cya, K, Cd
}