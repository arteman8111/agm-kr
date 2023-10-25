// ИД для всех
// Мой вариант 4
const N = 4; // Номер варианта из ЭУ
const h = (10 + 0.4*N) * math.pow(10,3) // Геометрическая высота [м]
const M_inf = 4.5 + 0.1*N // Число Маха для набегающего потока [-]
const alfa = 2 + 0.1*N // угол атаки [град]
const T_st = 373 // [К]
const Re_k = 5 * math.pow(10,6) // Критическое число рейнольдса
const cp0 = 1004.7; // удельная теплоемкость воздуха
const u0 = 1.79 * math.pow(10,-5); // динамическая вязкость 
const lymbda0 = 0.0232; // теплопроводность воздуха 
const c = 80 * math.pow(10,-3) // мм 
const b = 1000 * math.pow(10,-3) // мм 
const k = 1.4 // адиабата
const betta_k = math.atan(c/b)
const L = math.sqrt(math.pow(b / 2, 2) + math.pow(c / 2, 2)); // Длина грани


export {
    N,
    h,
    M_inf,
    alfa,
    T_st,
    Re_k,
    cp0,
    u0,
    lymbda0,
    c,
    b,
    k,
    betta_k,
    L
}