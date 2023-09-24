import * as param from "./const.js"
import {P,thet,gx,gy} from "./utils.js"
function vx_rk4(t, m, x, y) {
    return P(t) * math.cos(thet(t)) / m - gx(x,y);
}
function vy_rk4(t, m, x, y) {
    return P(t) * math.sin(thet(t)) / m - gy(x,y);
}
function x_rk4(v) {
    return v;
}
function y_rk4(v) {
    return v;
}
function m_rk4(t) {
    return - P(t) / param.W
}
function rk4(){
    let k1_vx = dt * vx_rk4(t, m, x, y);
    let k1_vy = dt * vy_rk4(t, m, x, y);
    let k1_x  = dt * x_rk4(vx);
    let k1_y  = dt * y_rk4(vy); 
    let k1_m  = dt * m_rk4(t)

    let k2_vx = dt * vx_rk4(t + dt/2, m + k1_m/2, x + k1_x/2, y + k1_y/2);
    let k2_vy = dt * vy_rk4(t + dt/2, m + k1_m/2, x + k1_x/2, y + k1_y/2);
    let k2_x  = dt * x_rk4(vx + k1_vx/2);
    let k2_y  = dt * y_rk4(vy + k1_vy/2);
    let k2_m  = dt * m_rk4(t + dt/2);

    let k3_vx = dt * vx_rk4(t + dt/2, m + k2_m/2, x + k2_x/2, y + k2_y/2);
    let k3_vy = dt * vy_rk4(t + dt/2, m + k2_m/2, x + k2_x/2, y + k2_y/2);
    let k3_x  = dt * x_rk4(vx + k2_vx/2);
    let k3_y  = dt * y_rk4(vy + k2_vy/2);
    let k3_m  = dt * m_rk4(t + dt/2);

    let k4_vx = dt * vx_rk4(t + dt, m + k3_m, x + k3_x, y + k3_y);
    let k4_vy = dt * vy_rk4(t + dt, m + k3_m, x + k3_x, y + k3_y);
    let k4_x  = dt * x_rk4(vx + k3_vx);
    let k4_y  = dt * y_rk4(vy + k3_vy);
    let k4_m  = dt * m_rk4(t + dt)


    vx += (k1_vx + 2 * k2_vx + 2 * k3_vx + k4_vx) / 6;
    vy += (k1_vy + 2 * k2_vy + 2 * k3_vy + k4_vy) / 6;
    x  += (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6;
    y  += (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6;
    m  += (k1_m + 2 * k2_m + 2 * k3_m + k4_m) / 6;
    r  = math.sqrt(math.pow(x, 2) + math.pow(y + param.rM, 2));
    v  = math.sqrt(math.pow(vx,2) + math.pow(vy,2))
}
// НУ
let t0 = 0;
let r0 = 0;
let vx = 0;
let vy = 0;
let v = 0;
let x = 0;
let y = 0;
let dt = param.step;
let m = param.m0;

// Радиус-вектор
let r = r0;
let t = t0;
function init() {
    while ( r < param.rk) {
        if (t === 0){
            console.log(`t = ${t.toFixed(6)} || P = ${P(t)} || vx = ${vx}`);
        }
        if (t + dt > param.t1 && t < param.t1){
            dt = param.t1 - t;
            t += dt;
            rk4();
            console.log(`t = ${t.toFixed(6)} || P = ${P(t)} || vx = ${vx}`);
            dt = param.step - dt;
            t += dt;
            rk4();
            console.log(`t = ${t.toFixed(6)} || P = ${P(t)} || vx = ${vx}`);
            dt = param.step
        }
        if (t + dt > param.t2 && t < param.t2){
            dt = param.t2 - t;
            t += dt;
            rk4();
            console.log(`t = ${t.toFixed(6)} || P = ${P(t)} || vx = ${vx}`);
            dt = param.step - dt;
            t += dt;
            rk4();
            console.log(`t = ${t.toFixed(6)} || P = ${P(t)} || vx = ${vx}`);
            dt = param.step
        }

        t  += dt;
        rk4()
        console.log(`t = ${t.toFixed(6)} || P = ${P(t)}  || vx = ${vx}`)
    }
}
init()