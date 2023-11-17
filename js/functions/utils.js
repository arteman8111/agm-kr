import * as param from "../variables/const.js"

const v_default = (M,a) => {
    return M * a
}

const ksu = (betta, M) => {
    return function(thet){
        const a = 0.5 * (param.k + 1) * math.pow(M * math.sin(thet), 2);
        const b = 1 + 0.5 * (param.k - 1) * math.pow(M * math.sin(thet), 2);
        return math.tan(thet) / math.tan(thet - betta) - a / b;
    }
}

const criticalReynoldsNumber = (Tr, M) => {
    const Tst_otn = param.T_st / Tr; // Относительная температура стенки;
    const x = (Tst_otn - 1) / (M * M);
    const y = 1 - 16 * x - 412.5 * math.pow(x, 2) - 35000 * math.pow(x, 3) - 375000 * math.pow(x, 4);
    return 5 * math.pow(10, 6) * y
}
const reynoldsNumber = (po, V, u, x = param.L) => po * V * x / u;

const rad = (value) => {
    return value * math.pi / 180;
}

const grad = (value) => {
    return value * 180 / math.pi
}

export {v_default,ksu,rad,grad,criticalReynoldsNumber,reynoldsNumber}