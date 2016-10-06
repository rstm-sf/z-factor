import numpy as np

'''
	sg  - specific gravity (0.57 < sg < 1.68).
'''
def calcPpc(sg):
	# Ppc - pseudocritical pressure, psia.
	# Sutton's correlations, B.C. Craft and M.F. Hawkins.
	return(756.8 - 131.0 * sg - 3.60 * sg * sg)


'''
	sg  - specific gravity (0.57 < sg < 1.68).
'''
def calcTpc(sg):
	# Tpc - pseudocritical temperature, K (degrees Rankine, 1(K) = 1*5/9 (°R)).
	# Sutton's correlations, B.C. Craft and M.F. Hawkins.
	return((169.2 + 349.5 * sg - 74.0 * sg * sg) * 5.0 / 9.0)


'''
	P   - pressure, atm;
	sg  - specific gravity (0.57 < sg < 1.68).
'''
def calcPpr(P, sg):
	# Ppr - pseudo reduced pressure (1 (atm) = 1*101325/6894.757293168 (psia)).
	# Dranchuk-Abbou Kassem: 0.2 < Ppc < 30.
	return(P * 101325 / 6894.757293168 / calcPpc(sg))


'''
	T  - temperature, °C;
	sg - specific gravity (0.57 < sg < 1.68).
'''
def calcTpr(T, sg):
	# Tpr - pseudo reduced temperature (1 (°C) = 1+273.15 (K)).
	# Dranchuk-Abbou Kassem: 1.0 < Tpc < 3.0.
	return((T + 273.15) / calcTpc(sg))


'''
	C   - coefficients Dranchuk equation;
	Rr  - reduced density;
	Tpr - pseudo reduced temperature, K;
	z   - compressibility factor.
'''
def fun_DAK(C, Rr, Tpr, z):
	Rr2 = Rr * Rr
	tmp = 0.7210 * Rr2
	C2  = 0.6134 * (1.0 + tmp) * Rr2 / (Tpr*Tpr*Tpr) * np.exp(-tmp)
	return (z - C - C2)


'''
	Ppr    - pseudo reduced pressure, psia;
	Tpr    - pseudo reduced temperature, K;
	ra, rb - z locate [ra, rb] (bisection method).
'''
def calcZfactor_DAK(Ppr, Tpr, za = 0.7, zb = 1.1):
	invTpr = 1.0 / Tpr
	tmp    = invTpr*invTpr
	Rr_z   = 0.27*Ppr * invTpr

	C   = (0.3265 - 1.07 * invTpr - 0.5339 * tmp * invTpr + 0.01569 * tmp*tmp -
		   0.05165 * tmp*tmp * invTpr) * Rr_z
	tmp = -0.7361 * invTpr + 0.1844 * tmp
	C  += 1.0 + ((0.5475 + tmp) - (0.1056 * tmp) * Rr_z*Rr_z*Rr_z) * Rr_z*Rr_z

	i       = 0
	maxIter = 100
	inv2    = 0.5
	epsilon = 2.0e-6
	a       = za
	b       = zb
	zn      = 0.0

	# The method bisection
	for i in range(maxIter):

		zn = (a + b) * inv2
		convergence = abs(b - a)
		if(convergence <= epsilon):
			break

		Rr  = Rr_z / zn
		fz = fun_DAK(C, Rr, Tpr, zn)

		if (fz > 0):
			b = zn
		elif (fz < 0):
			a = zn
		elif (fz == 0.0):
			break

	print('Iter =', end=' ')
	print(i)
	print('Convergence =', end=' ')
	print(convergence)

	if (i == maxIter - 1):
		print('Warning: max iter!\n')

	return zn


'''
	TEST solve (Applied Petroleum Reservoir Engineering. B.C. Craft, M.F. Hawkins)
'''
def test():
	P  = 3250.0 * 6894.757293168 / 101325
	T  = (213.0 - 32.0) * 5.0 / 9.0
	# sg  - specific gravity (0.57 < sg < 1.68).
	sg = 0.666
	z  = 0.0

	# Ppr - pseudo reduced pressure.
	# Tpr - pseudo reduced temperature.
	# Dranchuk-Abbou Kassem: 0.2 < Ppc < 30, 1.0 < Tpc < 3.0.
	Ppr = calcPpr(P, sg)
	Tpr = calcTpr(T, sg)
	print('Ppr =', end=' ')
	print(Ppr, end=', ')
	print('Tpr =', end=' ')
	print(Tpr)

	z = calcZfactor_DAK(Ppr, Tpr)

	print('Z =', end=' ')
	print(z)


def main():
	N  = 50
	P  = np.linspace(0, 500, N)
	T  = 0.0
	sg = 0.9
	z  = np.zeros(N)

	for i in range(N):

		Ppr  = calcPpr(P[i], sg)
		Tpr  = calcTpr(T, sg)
		z[i] = calcZfactor_DAK(Ppr, Tpr, 0.2, 2)

	print('Z =',)
	print(z)


main()