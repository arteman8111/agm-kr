const secantMethod = (f, init) => {
    let x0 = init;

    const delta = 0.0001;
    const f_toch = (f(x0) - f(x0 - delta)) / delta;
    const x1 = x0 - f(x0) / f_toch;
    const error = 0.0001;


    let x_prev = x0;
    let x_iter = x1;
    let x_next = -1;
    let b = 0;
    function calc(x_prev, x_iter) {
        return (x_iter - x_prev) * f(x_iter) / (f(x_iter) - f(x_prev))
    }
    while (math.abs(x_next - b) > error) {
        x_next = x_iter - calc(x_prev, x_iter);
        b = x_iter;
        x_prev = x_iter;
        x_iter = x_next;
    }
    return x_next
}

export { secantMethod }