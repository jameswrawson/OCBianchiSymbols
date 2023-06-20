class OverconvSymb:
	def __init__(self, parent, values): #values are a dictionary of /tuples/ with coset index, taking values as lists of lists, first index x, second y
		self.parent = parent
		self.values = values
	
	def hecke(self, ell, verbose=False):
		out_vals = {}
		num_csts = len(self.parent.dom.cosets)
		decile = floor(num_csts / 10)
		loop = 1
		for cst in self.parent.dom.cosets:
			if verbose and (loop % decile) == 0:
				print("Hecke {0}% done".format((loop / decile) * 10))
			mat = self.parent.dom.relations[0].get_rep(cst)
			out_vals[tuple(cst)] = [[0 for i in range(self.parent.prec)] for j in range(self.parent.prec)]
			key = tuple(cst)
			ell_id = self.parent.K.ideal(ell)
			for a in ell_id.residues():
				Up = matrix(self.parent.K, 2, 2, [1, a, 0, ell])
				new_mat = Up * mat
				val = self.parent.distAct(Up, self(new_mat))
				for i in range(self.parent.prec):
					for j in range(self.parent.prec):
						out_vals[key][i][j] += val[i][j]
			
			if self.parent.N not in ell_id:
				Up = matrix(self.parent.K, 2, 2, [ell, 0, 0, 1])
				new_mat = Up * mat
				val = self.parent.distAct(Up, self(new_mat))
				for i in range(self.parent.prec):
					for j in range(self.parent.prec):
						out_vals[key][i][j] += val[i][j]
			
			loop += 1
		
		if verbose:
			print("Hecke 100% done")
		
		return OverconvSymb(self.parent, out_vals)
		
		
	
	def __call__(self, matrix):
		out = [[0 for i in range(self.parent.prec)] for j in range(self.parent.prec)]
		path = self.parent.dom.path(matrix)
		for term in path.data:
			summand = self.parent.distAct(term[0].A.inverse(), self.values[tuple(term[1])])
			for i in range(self.parent.prec):
				for j in range(self.parent.prec):
					out[i][j] += self.parent.padics(term[0].n) * summand[i][j]
		
		return out
	
	def __add__(self, other):
		new_vals = {}
		if self.parent != other.parent:
			raise ValueError("Adding modular symbols from different spaces")
		
		for cst in self.parent.dom.cosets:
			k = tuple(cst)
			new_vals[k] = [[0 for i in range(self.parent.prec)] for j in range(self.parent.prec)]
			for i in range(self.parent.prec):
				for j in range(self.parent.prec):
					new_vals[k][i][j] = self.values[k][i][j] + other.values[k][i][j]
		
		return OverconvSymb(self.parent, new_vals)
		
	def __rmul__(self, const):
		new_vals = {}
		
		for cst in self.parent.dom.cosets:
			key = tuple(cst)
			new_vals[key] = [[0 for i in range(self.parent.prec)] for j in range(self.parent.prec)]
			for i in range(self.parent.prec):
				for j in range(self.parent.prec):
					new_vals[key][i][j] = self.parent.convert(const) * self.values[key][i][j]
		
		return OverconvSymb(self.parent, new_vals)
		
	def save(self, filename, append=True):
		f = 0
		if append:
			f = open(filename, "a")
			f.write("#"*10 + "\n")
		else:
			f = open(filename, "w")
		
		for cst in self.values:
			f.write("{0} : {1}\n".format(cst, self.values[cst]))
		
		f.close()
	
	def deformationDirection(self):  #returns lambda s.t. (1, lambda) is the infinitesimal weight it deforms in, False if it is not a smooth point
		errs1 = []
		errs2 = []
		for reln in self.parent.classical.dom.relations:
			cont1 = [0 for i in range((self.parent.wt+1)**2)]
			cont2 = [0 for i in range((self.parent.wt+1)**2)]
			for term in reln.data:
				a, b, c, d = term[0].A[0][0], term[0].A[0][1], term[0].A[1][0], term[0].A[1][1]
				abar, bbar, cbar, dbar = a.conjugate(), b.conjugate(), c.conjugate(), d.conjugate()
				a, b, c, d = self.parent.convert(a), self.parent.convert(b), self.parent.convert(c), self.parent.convert(d)
				abar, bbar, cbar, dbar = self.parent.convert(abar), self.parent.convert(bbar), self.parent.convert(cbar), self.parent.convert(dbar)
				u = (-c / d)
				ubar = (-cbar / dbar)
				
				R.<x, y> = PowerSeriesRing(self.parent.padics, default_prec = self.parent.prec)
				log1 = d.log()
				log2 = dbar.log()
				for i in range(1, self.parent.prec):
					log1 += ((-1)**(i + 1) * (u*x)**i)/i
					log2 += ((-1)**(i + 1) * (ubar*y)**i)/i
				
				f1 = (c + d*x)**self.parent.wt * (cbar + dbar*y)**self.parent.wt * log1
				f2 = (c + d*x)**self.parent.wt * (cbar + dbar*y)**self.parent.wt * log2
				
				for i in range(self.parent.wt + 1):
					for j in range(self.parent.wt + 1):
						err1 = 0
						err2 = 0
						g1 = f1 * ((a*x - b)/(d - c*x))**i * ((abar*y - bbar)/(dbar - cbar*y))**j
						g2 = f2 * ((a*x - b)/(d - c*x))**i * ((abar*y - bbar)/(dbar - cbar*y))**j
						g1coeffs = g1.coefficients()
						g2coeffs = g2.coefficients()
						for k in range(self.parent.prec):
							for l in range(self.parent.prec):
								if x**k * y**l in g1coeffs:
									err1 += g1coeffs[x**k * y**l] * self.values[tuple(term[1])][k][l] 
								
								if x**k * y**l in g2coeffs:
									err2 += g2coeffs[x**k * y**l] * self.values[tuple(term[1])][k][l] 
				
						err1 *= term[0].n
						err2 *= term[0].n
						cont1[(self.parent.wt + 1)*i + j] += err1
						cont2[(self.parent.wt + 1)*i + j] += err1
			
			errs1 = errs1 + cont1
			errs2 = errs2 + cont2
		
		reln_mat = [[self.parent.convert(elt) for elt in row] for row in self.parent.classical.relns]
		p_adic_relns = matrix(self.parent.padics, reln_mat, sparse=True)
		try:
			soln = p_adic_relns.solve_left(vector(self.parent.padics, errs1))
			soln = p_adic_relns.solve_left(vector(self.parent.padics, errs2))
			raise False
		except Exception as inst:
			reln_mat.append(errs2)
			p_adic_relns = matrix(self.parent.padics, reln_mat, sparse=True)
			try:
				soln = p_adic_relns.solve_left(vector(self.parent.padics, errs1))
				return -soln[-1]
			except Exception as inst2:
				return "inf"
		
