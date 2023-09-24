import * as param from "./const.js"

const gx = (x, y) => {
    return param.uM * x / math.pow((math.pow(x, 2) + math.pow((param.rM + y), 2)), 1.5)
}
const gy = (x, y) => {
    return param.uM * (param.rM + y) / math.pow((math.pow(x, 2) + math.pow((param.rM + y), 2)), 1.5)
}
const thet = (t) => {
    if (t <= param.tv) {
        return math.pi / 2
    } else if (t <= param.t1) {
        return math.pi / 2 + param.thet_torch * (t - param.tv)
    } else {
        return math.pi / 2 + param.thet_torch * (param.t1 - param.tv)
    }
}
const P = (t) => {
    if (t < param.t2 && t > param.t1){
        return 0;
    } else {
        return param.P2;
    }
}
export {P,thet,gx,gy}
