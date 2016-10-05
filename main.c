#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

/*
	Ppr - pseudo reduced pressure, psia;
	Tpr - pseudo reduced temperature, K;
	z   - compressibility factor.
*/
int8_t calcZfactor_DAK(double Ppr, double Tpr, double *z);

/*
	C   - coefficients Dranchuk equation;
	Rr  - reduced density;
	Tpr - pseudo reduced temperature, K;
	z   - compressibility factor.
*/
double fun_DAK(double C, double Rr, double Tpr, double z);

/*
    TEST solve (Applied Petroleum Reservoir Engineering. B.C. Craft, M.F. Hawkins)
*/
void test();

int8_t main() {

    test();

	system("pause");

	return 0;

}

int8_t calcZfactor_DAK(const double Ppr, const double Tpr, double *z) {

	double const invTpr = 1.0 / Tpr;
	double tmp          = invTpr*invTpr;
	double const Rr_z   = 0.27*Ppr * invTpr;

	double C;
	C   = (0.3265 - 1.07 * invTpr - 0.5339 * tmp * invTpr + 0.01569 * tmp*tmp -
		   0.05165 * tmp*tmp * invTpr) * Rr_z;
	tmp = -0.7361 * invTpr + 0.1844 * tmp;
	C  += 1.0 + ((0.5475 + tmp) - (0.1056 * tmp) * Rr_z*Rr_z*Rr_z) * Rr_z*Rr_z;

	uint16_t i             = 0;
	uint16_t const maxIter = 100;
	double const inv2      = 0.5;
	double const epsilon   = 2.0e-6;
	double a               = 0.6;
	double b               = 1.3;
	double convergence;
	double zn;

	// The method bisection
	for(i; i < maxIter; ++i) {

		zn = (a + b) * inv2;
		convergence = fabs(b - a);
		if(convergence <= epsilon)
			break;

		const double Rr  = Rr_z / zn;
		const double fz = fun_DAK(C, Rr, Tpr, zn);

		if (fz > 0) {
			b = zn;
		} else if (fz < 0) {
			a = zn;
		} else if (fz == 0.0) {
			break;
		} else
			return -2;
	}

	*z = zn;

	printf("Iter = %" PRIu16 "\nConvergence = %e\n", i, convergence);
	if (i == maxIter)
		printf("Warning: max iter = %" PRIu16 "\n", maxIter);

	return 0;

}

inline double fun_DAK(const double C, const double Rr, const double Tpr,
	const double z) {

	const double Rr2 = Rr * Rr;
	const double tmp = 0.7210 * Rr2;
	const double C2  = 0.6134 * (1.0 + tmp) * Rr2 / (Tpr*Tpr*Tpr) * exp(-tmp);
	return (z - C - C2);

}

void test() {

    const double P  = 3250.0 * 6894.757293168 / 101325;
    const double T  = (213.0 - 32.0) * 5.0 / 9.0;
    // sg  - specific gravity (0.57 < sg < 1.68).
    const double sg = 0.666;
    double z        = 0.0;

    // Ppc - pseudocritical pressure, psia.
    // Tpc - pseudocritical temperature, K (degrees Rankine, 1(K) = 1*5/9 (°R)).
    // Sutton's correlations, B.C. Craft and M.F. Hawkins.
    const double Ppc = 756.8 - 131.0 * sg - 3.60 * sg * sg;
    const double Tpc = (169.2 + 349.5 * sg - 74.0 * sg * sg) * 5.0 / 9.0;

    // Ppr - pseudo reduced pressure (1 (atm) = 1*101325/6894.757293168 (psia)).
    // Tpr - pseudo reduced temperature (1 (°C) = 1+273.15 (K)).
    // Dranchuk-Abbou Kassem: 0.2 < Ppc < 30, 1.0 < Tpc < 3.0.
    const double Ppr = P * 101325 / 6894.757293168 / Ppc;
    const double Tpr = (T + 273.15) / Tpc;
    printf("Ppr = %f, Tpr = %f\n", Ppr, Tpr);

    int8_t err = calcZfactor_DAK(Ppr, Tpr, &z);
    if(err < 0)
        printf("err = %" PRId16 "\n", err);

    printf("Z = %f\n", z);

}
