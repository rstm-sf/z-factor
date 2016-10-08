import numpy as np
import matplotlib.pyplot as plt
import time

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

	#print('Iter =', end=' ')
	#print(i)
	#print('Convergence =', end=' ')
	#print(convergence)

	if (i == maxIter - 1):
		print('Warning: max iter!\n')

	return zn


'''
	TEST 1: solve (Applied Petroleum Reservoir Engineering. B.C. Craft, M.F. Hawkins)
'''
def test1():
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


'''
	TEST 2: the construction of a family of curves.
'''
def test2():
	N   = 50
	M   = 20
	sg  = 0.661
	P   = np.linspace(0, 500, N)
	T   = np.linspace(-30, 200, M)
	Ppr = calcPpr(P, sg)
	Tpr = calcTpr(T, sg)

	z = np.zeros((M, N))
	for i in range(M):
		tmp = Tpr[i]
		for j in range(N):
			z[i, j] = calcZfactor_DAK(Ppr[j], tmp, 2.5e-4, 6)

	fig  = plt.figure()
	axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])

	for i in range(M):
		axes.plot(Ppr, z[i], label = 'T = ' + str(Tpr[i]))

	handles, labels = axes.get_legend_handles_labels()
	axes.legend(handles, labels, loc = 'upper left', ncol = 2, fontsize = 12)
	axes.set_ylim(z.min(), z.max())
	axes.set_xlim(Ppr.min(), Ppr.max())
	axes.set_ylabel('Compressibility factor Z')
	axes.set_xlabel('Pseudo reduced pressure')

	plt.grid()
	plt.show()


'''
	TEST 3: performance
'''
def test3():
	print('Выберете от чего будет зависеть сжимаемость. Варианты:')
	print('\t' + '1 - псевдо привиденое давление;')
	print('\t' + '2 - псевдо привиденая температура.')
	print('Поле ввода: ', end = '')
	dependence = int(input())
	if (dependence != 1 and dependence != 2):
		print('Incorrectly \'dependence\'!\Выход.')
		return

	M   = 20
	N   = 50
	sg  = 0.661
	z   = np.zeros((M, N))
	za  = 2.5e-4
	zb  = 16

	startTime = 0

	if (dependence == 1):
		startTime = time.time()

		P   = np.linspace(0, 500, N)
		T   = np.linspace(-30, 200, M)
		x = calcPpr(P, sg)
		const = calcTpr(T, sg)

		for i in range(M):
			tmp = const[i]
			for j in range(N):
				z[i, j] = calcZfactor_DAK(x[j], tmp, za, zb)

		str_xyc = ['Pseudo reduced pressure', 'Compressibility factor Z', 'Tpr']

	elif (dependence == 2):
		startTime = time.time()

		P   = np.linspace(0, 500, M)
		T   = np.linspace(-30, 200, N)
		const = calcPpr(P, sg)
		x = calcTpr(T, sg)

		for i in range(M):
			tmp = const[i]
			for j in range(N):
				z[i, j] = calcZfactor_DAK(tmp, x[j], za, zb)

		str_xyc = ['Pseudo reduced pressure', 'Compressibility factor Z', 'Ppr']

	clrs20 = ('#689f38','#009688','#b2dfdb','#e64a19','#00bcd4','#212121',
	          '#757575','#BDBDBD','#fbc02d','#ffeb3b','#0288d1','#03a9f4',
	          '#b3e5fc','#536dfe','#757575','#9C8BBF','#969CE7','#448aff',
	          '#EFBEAA','#c2185b')

	fig  = plt.figure()
	axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])

	str_label = str_xyc[2] + ' = '
	for i in range(M):
		axes.plot(x, z[i], c = clrs20[i], label = str_label + str(const[i]))

	handles, labels = axes.get_legend_handles_labels()
	axes.legend(handles, labels, loc = 'upper left', ncol = 2, fontsize = 10)
	axes.set_ylim(z.min(), z.max())
	axes.set_xlim(x.min(), x.max())
	axes.set_ylabel(str_xyc[1])
	axes.set_xlabel(str_xyc[0])
	plt.grid()

	endTime = time.time()
	print('Elapsed time: {:.3f} sec'.format(endTime - startTime))

	plt.show()


test3()