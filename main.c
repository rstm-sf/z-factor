#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

/*
    P  - pressures,   atm;
    T  - temperature, celsius;
    sg - specific gravity (0.57 < sg < 1.68).
*/
int8_t calcZfactor_DAK(double P, double T, double sg, double *z);

/*
    A9, A10 - coefficients Dranchuk equation;
    Rr      - reduced density;
    Tpr     - pseudo reduced temperature, psia.
*/
double coeff_C3(double A9, double A10, double Rr, double Tpr);

/*
    C  - coefficients Dranchuk equation;
    Rr - reduced density;
    z  - compressibility factor.
*/
double fun_DAK(double *C, double Rr, double z);

int8_t main() {

    const double P  = 3250.0 * 6894.757293168 / 101325;
    const double T  = (213.0 - 32.0) * 5.0 / 9.0;
    const double sg = 0.666;
    double z        = 0.0;
    int8_t err      = calcZfactor_DAK(P, T, sg, &z);
    if(err < 0)
        printf("err = %" PRId16 "\n", err);

    printf("Z = %f\n", z);

    system("pause");

    return 0;
}

int8_t calcZfactor_DAK(const double P, const double T, const double sg, double *z) {

    // Ppc - pseudocritical pressure, psia.
    // Tpc - pseudocritical temperature, K (degrees Rankine, 1(K) = 1*5/9 (°R)).
    // Suttons correlations, B.C. Craft and M.F. Hawkins.
    const double Ppc = 756.8 - 131.0 * sg - 3.60 * sg * sg;
    const double Tpc = (169.2 + 349.5 * sg - 74.0 * sg * sg) * 5.0 / 9.0;

    // Ppr - pseudo reduced pressure (1 (atm) = 1*101325/6894.757293168 (psia)).
    // Tpr - pseudo reduced temperature (1 (°C) = 1+273.15 (K)).
    // Dranchuk-Abbou Kassem: 0.2 < Ppc < 30, 1.0 < Tpc < 3.0.
    const double Ppr = P * 101325 / 6894.757293168 / Ppc;
    const double Tpr = (T + 273.15) / Tpc;
    printf("Ppr = %f, Tpr = %f\n", Ppr, Tpr);

    double A[11];
    A[0] =  0.32650; A[1] = -1.07000; A[2]  = -0.5339; A[3] = 0.01569;
    A[4] = -0.05165; A[5] =  0.54750; A[6]  = -0.7361; A[7] = 0.18440;
    A[8] =  0.10560; A[9] =  0.61340; A[10] =  0.7210;

    const double invTpr = 1.0 / Tpr;
    double tmp          = invTpr*invTpr;

    double C[4];
    C[0] = A[0] + A[1] * invTpr + A[2] * tmp * invTpr + A[3] * tmp*tmp +
           A[4] * tmp*tmp * invTpr;
    tmp  = A[6] * invTpr + A[7] * tmp;
    C[1] = A[5] + tmp;
    C[2] = A[8] * tmp;

    double a    = 1e-2;
    double b    = 4.0;
    double Rra  = 0.27*Ppr * invTpr / a;
    double Rrb  = 0.27*Ppr * invTpr / b;
    C[3]        = coeff_C3(A[9], A[10], Rra, Tpr);
    double fa   = fun_DAK(C, Rra, a);
    C[3]        = coeff_C3(A[9], A[10], Rrb, Tpr);
    double fb   = fun_DAK(C, Rrb, b);

    if (fa * fb > 0.0) {
		printf("f(a)f(b) > 0\n");
		return -1;
	} else if (fa == 0.0) {
        *z = a;
        return 0;
    } else if (fb == 0.0) {
        *z = b;
        return 0;
    }

    uint16_t i             = 0;
    uint16_t const maxIter = 100;
    const double inv2      = 1.0 / 2.0;
    const double epsilon   = 2 * 10e-6;
    double convergence;

    // The method bisection
    for(i; i < maxIter; ++i) {

        *z = (a + b) * inv2;
        convergence = fabs(b - a);
        if(convergence <= epsilon)
            break;

        const double Rr  = 0.27*Ppr * invTpr / *z;
        C[3] = coeff_C3(A[9], A[10], Rr, Tpr);
        const double fz = fun_DAK(C, Rr, *z);

        if (fa * fz < 0) {
            b  = *z;
            fb = fz;
        } else if (fb * fz < 0) {
            a  = *z;
            fa = fz;
        } else if (fz == 0.0) {
            break;
        } else
			return -2;
    }

    printf("Iter = %" PRIu16 "\nConvergence = %e\n", i, convergence);
    if (i == maxIter)
        printf("Warning: max iter = %" PRIu16 "\n", maxIter);

    return 0;

} const;

inline double coeff_C3(const double A9, const double A10, const double Rr,
                 const double Tpr) {
    const double Rr2 = Rr * Rr;
    const double tmp = A10 * Rr2;
    return (A9 * (1.0 + tmp) * Rr2 / (Tpr*Tpr*Tpr) * exp(-tmp));
} const;

inline double fun_DAK(const double *C, const double Rr, const double z) {
    const double Rr2 = Rr * Rr;
    return (z - 1.0 - C[0] * Rr - C[1] * Rr2 + C[2] * Rr2*Rr2*Rr - C[3]);
} const;
