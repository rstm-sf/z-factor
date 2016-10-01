#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

double DAK(double Ppr, double Tpr);

int main() {

    const double Ppr = 2.958;
    const double Tpr = 1.867;
    const double z = DAK(Ppr, Tpr);
    printf("Z = %f\n", z);

    return 0;
}

double DAK(const double Ppr, const double Tpr) {

    double invTpr[5];
    invTpr[0] = 1.0 / Tpr;
    for(uint8_t i = 1; i < 5; ++i)
        invTpr[i] = invTpr[i-1] * invTpr[0];

    double A[11];
    A[0] =  0.32650; A[1] = -1.07000; A[2] = -0.5339; A[3] = 0.01569;
    A[4] = -0.05165; A[5] =  0.54750; A[6] = -0.7361; A[7] = 0.18440;
    A[8] =  0.10560; A[9] =  0.61340; A[10] = 0.7210;

    double tmp;

    double C[4];
    C[0] = A[0] + A[1] * invTpr[0] + A[2] * invTpr[2] +
            A[3] * invTpr[3] + A[4] * invTpr[4];
    tmp = A[6] * invTpr[0] + A[7] * invTpr[1];
    C[1] = A[5] + tmp;
    C[2] = A[8] * tmp;

    double z = 1.0;
    uint16_t i = 0;
    uint16_t const maxLoop = 100;
    double convergence;

    for(i; i < maxLoop; ++i) {

        const double Rr  = 0.27*Ppr * invTpr[0] / z;
        const double Rr2 = Rr * Rr;
        const double Rr5 = Rr2 * Rr2 * Rr;

        tmp = A[10] * Rr2;

        C[3] = A[9] * (1.0 + tmp) * Rr2 * invTpr[2] * exp(-tmp);

        const double zn = 1.0 + C[0] * Rr + C[1] * Rr2 - C[2] * Rr5 + C[3];
        const double invZn = 1.0 / zn;
        const double dfdz = 1.0 + C[0] * Rr * invZn + 2.0*C[1] * Rr2 * invZn +
                             5.0*C[2] * Rr5 * invZn +
                             2.0*A[9] * Rr2 * invTpr[2] * invZn *
                              (1.0 + tmp - tmp * tmp * exp(-tmp));
        const double dz = (zn - z) / dfdz;

        z += dz;
        convergence = fabs(z - zn);
        if(convergence <= 1.0e-6)
            break;
    }

    printf("Loops = %" PRIu16 "\nConvergence = %e\n", i, convergence);
    if (i == maxLoop)
        printf("Warning: max loop = %" PRIu16 "\n", maxLoop);

    return z;

} const;