import numpy as np
import matplotlib.pyplot as plt
import time
import math

'''
	sg  - specific gravity (0.57 < sg < 1.68).
	return: Ppc - pseudocritical pressure, psia.
'''
def calcPpc(sg):
	# Sutton's correlations, B.C. Craft and M.F. Hawkins.
	return(756.8 - 131.0 * sg - 3.60 * sg * sg)


'''
	sg  - specific gravity (0.57 < sg < 1.68).
	return: Tpc - pseudocritical temperature, K.
'''
def calcTpc(sg):
	# 1(K) = 1*5/9 (°R).
	# Sutton's correlations, B.C. Craft and M.F. Hawkins.
	return((169.2 + 349.5 * sg - 74.0 * sg * sg) * 5.0 / 9.0)


'''
	P   - pressure, atm;
	sg  - specific gravity (0.57 < sg < 1.68).
	return: Ppr - pseudo reduced pressure.
'''
def calcPpr(P, sg):
	# 1 (atm) = 1*101325/6894.757293168 (psia).
	# Dranchuk-Abbou Kassem: 0.2 < Ppc < 30.
	return(P * 101325 / 6894.757293168 / calcPpc(sg))


'''
	T  - temperature, °C;
	sg - specific gravity (0.57 < sg < 1.68).
	return: Tpr - pseudo reduced temperature.
'''
def calcTpr(T, sg):
	# 1 (°C) = 1+273.15 (K).
	# Dranchuk-Abbou Kassem: 1.0 < Tpc < 3.0.
	return((T + 273.15) / calcTpc(sg))


'''
	Ppr    - pseudo reduced pressure, psia;
	Tpr    - pseudo reduced temperature, K;
	za, zb - z locate [za, zb] (bisection method).
	return: z - gas compressibility factor based on Dranchuk-Abbou Kassem EoS.
'''
def calcZfactor_DAK(Ppr, Tpr, za = 0.7, zb = 1.1):
	invTpr  = 1.0 / Tpr
	invTpr2 = invTpr*invTpr
	invTpr3 = invTpr2*invTpr
	Rr_z    = 0.27*Ppr * invTpr
	Rr_z2   = Rr_z*Rr_z

	C1  = (0.3265 - 1.07 * invTpr - 0.5339 * invTpr3 +
		  0.01569 * invTpr2*invTpr2 - 0.05165 * invTpr2*invTpr3) * Rr_z
	tmp = -0.7361 * invTpr + 0.1844 * invTpr2
	C2  = (0.5475 + tmp) * Rr_z2
	C3  = 0.1056 * tmp * Rr_z2*Rr_z2*Rr_z
	C4  = 0.6134 * Rr_z2 * invTpr3
	C5  = 0.7210 * Rr_z2

	i       = 0
	maxIter = 100
	inv2    = 0.5
	epsilon = 2.0e-6
	a       = za
	b       = zb
	zn      = 0.0
	one     = 1.0

	# The method bisection
	for i in range(maxIter):

		zn = (a + b) * inv2
		convergence = abs(b - a)
		if(convergence <= epsilon):
			break

		invZn  = one / zn
		invZn2 = invZn*invZn
		tmp = C5 * invZn2
		fz = (zn - one - C1 * invZn - C2 * invZn2 + C3 * invZn2*invZn2*invZn -
			 C4 * invZn2 * (one + tmp) * math.exp(-tmp))

		if (fz > 0):
			b = zn
		elif (fz < 0):
			a = zn
		elif (fz == 0.0):
			break

	if (i == maxIter - 1):
		print('calcZfactor_DAK(). Warning: max iter!\n')

	return zn


'''
	Ppr    - pseudo reduced pressure, psia;
	Tpr    - pseudo reduced temperature, K;
	da, db - dZdT locate [da, db] (bisection method).
	za, zb - z locate [za, zb] (bisection method).
	return: dZ/dTpr,
	Z - gas compressibility factor based on Dranchuk-Abbou Kassem EoS.
'''
def calc_dZdTpr(Ppr, Tpr, da, db, za = 0.7, zb = 1.1):
	z       = calcZfactor_DAK(Ppr, Tpr, za, zb)
	dRrdT   = 0.27*Ppr / (Tpr*Tpr * z)
	i       = 0
	maxIter = 100
	inv2    = 0.5
	epsilon = 2.0e-6
	a       = da
	b       = db
	dZdTn   = 0.0

	for i in range(maxIter):

		dZdTn = (a + b) * inv2
		convergence = abs(b - a)
		if(convergence <= epsilon):
			break

		fz = dZdTn + dRrdT

		if (fz > 0):
			b = dZdTn
		elif (fz < 0):
			a = dZdTn
		elif (fz == 0.0):
			break

	return dZdTn

	if (i == maxIter - 1):
		print('calc_dZdTpr(). Warning: max iter!\n')


'''
	Ppr    - pseudo reduced pressure, psia;
	Tpr    - pseudo reduced temperature, K;
	da, db - dZdPr locate [da, db] (bisection method).
	za, zb - z locate [za, zb] (bisection method).
	return: dZ/dPpr,
	Z   - gas compressibility factor based on Dranchuk-Abbou Kassem EoS.
'''
def calc_dZdPpr(Ppr, Tpr, da, db, za = 0.7, zb = 1.1):
	z       = calcZfactor_DAK(Ppr, Tpr, za, zb)
	dRrdPr  = 0.27 / (Tpr * z)
	i       = 0
	maxIter = 100
	inv2    = 0.5
	epsilon = 2.0e-6
	a       = da
	b       = db
	dZdPrn  = 0.0

	for i in range(maxIter):

		dZdPrn = (a + b) * inv2
		convergence = abs(b - a)
		if(convergence <= epsilon):
			break

		fz = dZdPrn - dRrdPr

		if (fz > 0):
			b = dZdPrn
		elif (fz < 0):
			a = dZdPrn
		elif (fz == 0.0):
			break

	if (i == maxIter - 1):
		print('calc_dZdPpr(). Warning: max iter!\n')

	return dZdPrn


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
		axes.plot(Ppr, z[i], label = 'Tpr = ' + str(Tpr[i]))

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
	print('Выберете один из следующих вариантов:')
	print('\t' + '1. Построить набор из 20 кривых сжимаемости от давления ' +
	             'при заданных значениях температуры по 50 точек на кривой ' +
	             'равномерно распределенных.')
	print('\t' + '2. Построить набор из 20 кривых сжимаемости от температуры '+
	             'при заданных значениях давления по 50 точек на кривой '+
	             ' равномерно распределенных.')
	print('\t' + '3. Построить набор из 20 кривых производной сжимаемости ' +
	             'по давлению как зависимость от давления при заданных значениях '+
	             'температуры по 50 точек на кривой равномерно распределенных.')
	print('\t' + '4. Построить набор из 20 кривых производной сжимаемости '+
	             'по температуре как зависимость от давления при заданных '+
	             'значениях температуры по 50 точек на кривой равномерно '+
	             'распределенных.')
	print('\nПоле ввода: ', end = '')
	dependence = int(input())
	if (dependence < 1 and dependence > 4):
		print('Неверный ввод!\Выход.')
		return

	M   = 20
	N   = 50
	sg  = 0.661
	y   = np.zeros((M, N))
	za  = 2.5e-2
	zb  = 16

	startTime = 0

	if (dependence == 1):
		startTime = time.time()

		P     = np.linspace(0, 500, N)
		T     = np.linspace(-30, 200, M)
		x     = calcPpr(P, sg)
		const = calcTpr(T, sg)

		for i in range(M):
			tmp = const[i]
			for j in range(N):
				y[i, j] = calcZfactor_DAK(x[j], tmp, za, zb)

		str_xyc = ['Pseudo reduced pressure', 'Compressibility factor Z', 'Tpr',
		            'lower right']

	elif (dependence == 2):
		startTime = time.time()

		P     = np.linspace(1, 500, M)
		T     = np.linspace(-30, 200, N)
		const = calcPpr(P, sg)
		x     = calcTpr(T, sg)

		for i in range(M):
			tmp = const[i]
			for j in range(N):
				y[i, j] = calcZfactor_DAK(tmp, x[j], za, zb)

		str_xyc = ['Pseudo reduced temperature', 'Compressibility factor Z', 'Ppr',
		            'lower right']

	elif (dependence == 3):
		startTime = time.time()

		P     = np.linspace(1, 500, M)
		T     = np.linspace(-30, 200, N)
		const = calcPpr(P, sg)
		x     = calcTpr(T, sg)

		for i in range(M):
			tmp = const[i]
			for j in range(N):
				y[i, j] = calc_dZdTpr(tmp, x[j], -zb, -za, za, zb)

		str_xyc = ['Pseudo reduced temperature', 'dZ/dTpr', 'Ppr',
		            'lower right']

	elif (dependence == 4):
		startTime = time.time()

		P     = np.linspace(1, 500, M)
		T     = np.linspace(-30, 200, N)
		const = calcPpr(P, sg)
		x     = calcTpr(T, sg)

		for i in range(M):
			tmp = const[i]
			for j in range(N):
				y[i, j] = calc_dZdPpr(tmp, x[j], za, zb, za, zb)

		str_xyc = ['Pseudo reduced temperature', 'dZ/dPpr', 'Ppr',
		            'upper right']

	clrs20 = ('#689f38','#009688','#b2dfdb','#e64a19','#00bcd4','#212121',
	          '#757575','#BDBDBD','#fbc02d','#ffeb3b','#0288d1','#03a9f4',
	          '#b3e5fc','#536dfe','#757575','#9C8BBF','#969CE7','#448aff',
	          '#EFBEAA','#c2185b')

	fig  = plt.figure()
	axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])

	str_label = str_xyc[2] + ' = '
	for i in range(M):
		axes.plot(x, y[i], c = clrs20[i], label = str_label + str(const[i]))

	handles, labels = axes.get_legend_handles_labels()
	axes.legend(handles, labels, loc = str_xyc[3], ncol = 2, fontsize = 10)
	axes.set_ylim(y.min(), y.max())
	axes.set_xlim(x.min(), x.max())
	axes.set_ylabel(str_xyc[1])
	axes.set_xlabel(str_xyc[0])
	plt.grid()

	endTime = time.time()
	print('Прошло времени: {:.3f} с'.format(endTime - startTime))

	plt.show()


test3()