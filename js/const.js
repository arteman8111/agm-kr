// ИД для всех
const m0 = 3015 // кг ~ начальная масса СВ
const m_pg = 1550 // кг ~ Предельная масса рабочего запаса топлива
const m_const = 615 // кг ~ Масса конструкции
const tv = 14 // c ~ продолжительность вертикального участка выведения
const step = 0.1 // с ~ шаг инт
const uM = 4903 * math.pow(10,9) // км^3/с^2 ~ гравитационная постоянная Луны
const rM = 1738 * math.pow(10,3) // км ~ радиус Луны

// ИД мои
const W = 3.360 * math.pow(10,3) // м/с ~ Эффективная скорость истечения газов

const P1 = 8.30 * math.pow(10,3) // кН ~ Тяга СВ 1 АУТ
const h_isl_1_1 = 140 * math.pow(10,3)// км ~ высота орбиты для P1
const h_isl_1_2 = 180 * math.pow(10,3)// км ~ высота орбиты для P1

const P2 = 9.90 * math.pow(10,3) // кН ~ Тяга СВ 2 АУТ
const h_isl_2_1 = 180 * math.pow(10,3)// км ~ высота орбиты для P2
const h_isl_2_2 = 310 * math.pow(10,3) // км ~ высота орбиты для P2

// Подгоняемые параметры
const t1 = 5.01; // с
const t2 = 7.12; // с
const thet_torch = -math.pi; // рад

// Конечные параметры
const vk = math.sqrt(uM / (rM + h_isl_2_2));
const rk = h_isl_2_2 + rM;
export {
    m0,
    m_pg,
    m_const,
    tv,
    step,
    uM,
    rM,
    W,
    P1,
    P2,
    h_isl_1_1,
    h_isl_1_2,
    h_isl_2_1,
    h_isl_2_2,
    t1,
    t2,
    thet_torch,
    vk,
    rk
}