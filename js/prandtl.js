import { k, betta_k } from "./const.js";

const w = (M) => math.sqrt((k + 1) / (k - 1)) * math.atan(math.sqrt((M * M - 1) * (k - 1) / (k + 1))) - math.atan(math.sqrt(M * M - 1))
const w_next = (w) => w + 2 * betta_k

export {
    w,
    w_next
}