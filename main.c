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
	Rr  - reduced density;
	Tpr - pseudo reduced temperature, psia.
*/
double coeff_C3(double Rr, double Tpr);

/*
	C  - coefficients Dranchuk equation;
	Rr - reduced density;
	z  - compressibility factor.
*/
double fun_DAK(double *C, double Rr, double z);

int8_t main() {

	const double P  = 3250.0 * 6894.757293168 / 101325;
	const double T  = (213.0 - 32.0) * 5.0 / 9.0;
	// sg  - specific gravity (0.57 < sg < 1.68).
	const double sg = 0.666;
	double z        = 0.0;

	// Ppc - pseudocritical pressure, psia.
	// Tpc - pseudocritical temperature, K (degrees Rankine, 1(K) = 1*5/9 (°R)).
	// Sutton's correlations, B.C. Craft and M.F. Hawkins.
	const double Ppc = 756.8 - 131.0 * sg - 3.60 * sg * sg;
	const double Tpc = 94.00000752 + 194.1666822 * sg - 41.1111144 * sg * sg;

	// Ppr - pseudo reduced pressure (1 (atm) = 1*101325/6894.757293168 (psia)).
	// Tpr - pseudo reduced temperature (1 (°C) = 1+273.15 (K)).
	// Dranchuk-Abbou Kassem: 0.2 < Ppc < 30, 1.0 < Tpc < 3.0.
	const double Ppr = P * 14.695948775514218902863070110439 / Ppc;
	const double Tpr = (T + 273.15) / Tpc;
	printf("Ppr = %f, Tpr = %f\n", Ppr, Tpr);

	int8_t err      = calcZfactor_DAK(Ppr, Tpr, &z);
	if(err < 0)
		printf("err = %" PRId16 "\n", err);

	printf("Z = %f\n", z);

	system("pause");

	return 0;

}

int8_t calcZfactor_DAK(const double Ppr, const double Tpr, double *z) {

	double const invTpr = 1.0 / Tpr;
	double tmp          = invTpr*invTpr;

	double C[4];
	C[0] = 0.3265 - 1.07 * invTpr - 0.5339 * tmp * invTpr + 0.01569 * tmp*tmp -
		   0.05165 * tmp*tmp * invTpr;
	tmp  = -0.7361 * invTpr + 0.1844 * tmp;
	C[1] = 0.5475 + tmp;
	C[2] = 0.1056 * tmp;

	uint16_t i             = 0;
	uint16_t const maxIter = 100;
	double const inv2      = 0.5;
	double const epsilon   = 2 * 10e-6;
	double const Rr_z      = 0.27*Ppr * invTpr;
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
		C[3] = coeff_C3(Rr, Tpr);
		const double fz = fun_DAK(C, Rr, zn);

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

inline double coeff_C3(const double Rr, const double Tpr) {

	const double Rr2 = Rr * Rr;
	const double tmp = 0.7210 * Rr2;
	return (0.6134 * (1.0 + tmp) * Rr2 / (Tpr*Tpr*Tpr) * exp(-tmp));

}

inline double fun_DAK(const double *C, const double Rr, const double z) {

	const double Rr2 = Rr * Rr;
	return (z - 1.0 - C[0] * Rr - C[1] * Rr2 + C[2] * Rr2*Rr2*Rr - C[3]);

}
