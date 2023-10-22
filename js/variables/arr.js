import * as param from "./const.js"
import { atm } from "../atm.js"
import { rad, ksu, v_default } from "../functions/utils.js";
import { secantMethod } from "../methods/secant-method.js";

const atmos = atm(param.h);

let betta = [param.betta_k - rad(param.alfa), param.betta_k + rad(param.alfa), param.betta_k, param.betta_k];
let thet = [secantMethod(ksu(betta[0], param.M_inf), rad(12)), secantMethod(ksu(betta[1], param.M_inf), rad(12))];
let M = [param.M_inf];
let p = [atmos.p];
let po = [atmos.po];
let T = [atmos.T];
let a = [atmos.a]
let v = [v_default(M[0], a[0])];
let um = [];
let p0 = [];
let po0 = [];
let T0 = [];
let cp = [];
let u = [];
let lymbda = [];


export {
    betta,
    thet,
    M,
    p,
    po,
    T,
    a,
    v,
    um,
    p0,
    po0,
    T0,
    cp,
    u,
    lymbda
}