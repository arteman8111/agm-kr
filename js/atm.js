const atm = (h) => {
    const r = 6356767;
    const gc = 9.80665;
    const x = 1.4;
    const R = 287.05287;
    const H = r*h/(r+h);
    let bm, Tm, Hm, pm;

    if (h < -2000){
        return
    }
    if (h < 0){
        bm = -0.0065;Tm = 301.15;Hm = -2000;pm = 127774;
    } else if (h < 11000){
        bm = -0.0065;Tm = 288.15;Hm = 0;    pm = 101325;
    } else if (h < 20000){
        bm =  0;     Tm = 216.65;Hm = 11000;pm = 22632;
    } else if (h < 32000){
        bm =  0.001; Tm = 216.65;Hm = 20000;pm = 5474.87;
    } else if (h < 47000){
        bm =  0.0028;Tm = 228.65;Hm = 32000;pm = 868.014;
    } else if (h < 51000){
        bm =  0     ;Tm = 270.65;Hm = 47000;pm = 110.906;
    } else if (h < 71000){
        bm = -0.0028;Tm = 270.65;Hm = 51000;pm = 66.9384;
    } else if (h < 85000){
        bm = -0.002 ;Tm = 214.65;Hm = 71000;pm = 3.95639;
    }
    
    const T = Tm + bm * (H - Hm);
    const p = bm ? pm * math.exp(-gc * math.log(T / Tm) / (bm * R)) : pm * math.exp(-gc * (H - Hm) / (R * T)); 
    const po = p/(R * T);
    const a = math.sqrt(x * R * T);
    const g = gc * math.pow(r / (r + h), 2);
    return {
        T,
        p,
        po,
        a,
        g,
        H
    }
}
export {atm}