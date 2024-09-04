export const secantMethod = (f, x0, tolerance = 1e-4, maxIterations = 100) => {
    let x1 = x0 - f(x0) / ((f(x0) - f(x0 - tolerance)) / tolerance);

    for (let i = 0; i < maxIterations; i++) {
        const f_x1 = f(x1);
        const x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f(x0));

        if (Math.abs(x2 - x1) < tolerance) return x2;

        x0 = x1;
        x1 = x2;
    }

    return x1; // Return last iteration if convergence not reached
};

// Utility functions
export const rad = (value) => value * Math.PI / 180;
export const grad = (value) => value * 180 / Math.PI;
