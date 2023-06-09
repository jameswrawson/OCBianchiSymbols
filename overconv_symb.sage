class OverconvSymb:
	def __init__(self, parent, values): #values are a dictionary of /tuples/ with coset index, taking values as lists of lists, first index x, second y
		self.parent = parent
		self.values = values
	
	def hecke(self, ell):
		out_vals = {}
		for cst in self.parent.dom.cosets:
			mat = self.parent.dom.relations[0].get_rep(cst)
			out_vals[tuple(cst)] = [[0 for i in range(self.parent.prec)] for j in range(self.parent.prec)]
			key = tuple(cst)
			ell_id = self.parent.K.ideal(ell)
			for a in ell_id.residues():
				Up = matrix(self.parent.K, 2, 2, [1, a, 0, ell])
				new_mat = Up * mat
				val = self.parent.dist_act(Up, self(new_mat))
				for i in range(self.parent.prec):
					for j in range(self.parent.prec):
						out_vals[key][i][j] += val[i][j]
			
			if self.parent.N not in ell_id:
				Up = matrix(self.parent.K, 2, 2, [ell, 0, 0, 1])
				new_mat = Up * mat
				val = self.parent.dist_act(Up, self(new_mat))
				for i in range(self.parent.prec):
					for j in range(self.parent.prec):
						out_vals[key][i][j] += val[i][j]
		return OverconvSymb(self.parent, out_vals)
		
		
	
	def __call__(self, matrix):
		out = [[0 for i in range(self.parent.prec)] for j in range(self.parent.prec)]
		path = self.parent.dom.path(matrix)
		for term in path.data:
			summand = self.parent.dist_act(term[0].A.inverse(), self.values[tuple(term[1])])
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
			err1 = 0
			err2 = 0
			for term in reln.data:
				d, c = term[0].A[1][1], term[0].A[1][0]
				dbar, cbar = term[0].A[1][1].conjugate(), term[0].A[1][0].conjugate()
				d, c = self.parent.convert(d), self.parent.convert(c)
				dbar, cbar = self.parent.convert(dbar), self.parent.convert(cbar)
				u = (-c / d)
				ubar = (-cbar / dbar)
				
				cont1 = 0
				cont2 = 0
				cont1 += d.log() * self.values[tuple(term[1])][0][0]
				cont2 += dbar.log() * self.values[tuple(term[1])][0][0]
				for i in range(1, self.parent.prec):
					cont1 += (((-1)**(i + 1) * u**i)/i) * self.values[tuple(term[1])][i][0]
					cont2 += (((-1)**(i + 1) * ubar**i)/i) * self.values[tuple(term[1])][0][i]
				
				err1 += term[0].n * cont1
				err2 += term[0].n * cont2
			
			errs1.append(err1)
			errs2.append(err2)
		
		reln_mat = [[self.convert(elt) for elt in row] for row in self.parent.classical.relns]
		p_adic_relns = matrix(self.parent.padics, reln_mat, sparse=True)
		try:
			soln = p_adic_relns.solve_left(errs1)
			soln = p_adic_relns.solve_left(errs2)
			raise False
		except:
			reln_mat.append(errs2)
			p_adic_relns = matrix(self.parent.padics, reln_mat, sparse=True)
			try:
				soln = p_adic_relns.solve_left(errs1)
				return -soln[-1]
			except:
				return "inf"
		
