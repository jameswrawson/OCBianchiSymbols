def get_cpts(K, x):
	a = K.gen()
	tr = a.trace()
	tr2 = (a**2).trace()
	trx = x.trace()
	trax = (a*x).trace()
	k = 0
	l = 0
	if tr == 0:
		k = trx/2
		l = trax/tr2
	else:
		trax = (a*x).trace()
		l = (trx * tr - 2*trax)/(tr**2 - 2*tr2)
		k = (trx - l*tr)/2
	
	return (k, l)

def ctd_frac(K, x):
	quots = []
	while True:
		y = nearest_int(K, x)
		x = x - y
		quots.append(y)
		if x != 0:
			x = 1 / x
		else:
			break
	
	hn = [0, 1]
	kn = [1, 0]
	for quot in quots:
		new_h = quot * hn[-1] + hn[-2]
		new_k = quot * kn[-1] + kn[-2]
		hn.append(new_h)
		kn.append(new_k)
	
	mats = []
	for i in range(1, len(hn) - 1):
		A = hn[i]
		B = hn[i + 1]
		C = kn[i]
		D = kn[i + 1]
		det = A*D - B*C
		if det == 1:
			mat = matrix(K, 2, 2, [A, B, C, D])
			mats.append(mat)
		else:
			mat = matrix(K, 2, 2, [A, B / det, C, D / det])
			mats.append(mat)
	
	return mats

def nearest_int(K, alpha):
	a = K.gen()
	k, l = get_cpts(K, alpha)
	k0 = floor(k)
	l0 = floor(l)
	min_dist = 100
	min_int = 0
	for i in range(0, 2):
		for j in range(0, 2):
			beta = (k0 + i) + (l0 + j)*a
			dist = (alpha - beta).norm()
			if min_dist > dist:
				min_dist = dist
				min_int = beta
	
	return min_int

def bezout(K, alpha, beta):
	u_old = 1
	v_old = 0
	u_new = 0
	v_new = 1
	while beta != 0:
		gamma = nearest_int(K, alpha / beta)
		new_beta = alpha - (gamma * beta)
		u_newnew = u_old - (gamma * u_new)
		v_newnew = v_old - (gamma * v_new)
		alpha = beta
		beta = new_beta
		u_old = u_new
		u_new = u_newnew
		v_old = v_new
		v_new = v_newnew
	
	return (u_old, v_old, alpha)

