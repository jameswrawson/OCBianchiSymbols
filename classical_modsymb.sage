import fund_dom

class ClassicalModSymb:
	def __init__(self, parent, values): #wt is normalised like a modular symbol, wt = 0 is constant coefficients. values are /tuples/ of coset index dictionary, taking values as lists of lists, first index x, second y
		self.parent = parent
		self.values = values
	
	def hecke(self, ell):
		out_vals = {}
		for cst in self.parent.dom.cosets:
			mat = self.parent.dom.relations[0].get_rep(cst)
			out_vals[tuple(cst)] = [[0 for i in range(self.parent.wt + 1)] for j in range(self.parent.wt + 1)]
			key = tuple(cst)
			ell_id = self.parent.K.ideal(ell)
			for a in ell_id.residues():
				Up = matrix(self.parent.K, 2, 2, [1, a, 0, ell])
				new_mat = Up * mat
				val = self.parent.dist_act(Up, self(new_mat))
				for i in range(self.parent.wt + 1):
					for j in range(self.parent.wt + 1):
						out_vals[key][i][j] += val[i][j]
			
			if self.parent.N not in ell_id:
				Up = matrix(self.parent.K, 2, 2, [ell, 0, 0, 1])
				new_mat = Up * mat
				val = self.parent.dist_act(Up, self(new_mat))
				for i in range(self.parent.wt + 1):
					for j in range(self.parent.wt + 1):
						out_vals[key][i][j] += val[i][j]
		return ClassicalModSymb(self.parent, out_vals)
		
		
	
	def __call__(self, matrix):
		out = [[0 for i in range(self.parent.wt + 1)] for j in range(self.parent.wt + 1)]
		path = self.parent.dom.path(matrix)
		for term in path.data:
			summand = self.parent.dist_act(term[0].A.inverse(), self.values[tuple(term[1])])
			for i in range(self.parent.wt + 1):
				for j in range(self.parent.wt + 1):
					out[i][j] += term[0].n * summand[i][j]
		
		return out
	
	def __add__(self, other):
		new_vals = {}
		if self.parent != other.parent:
			raise ValueError("Adding modular symbols from different spaces")
		
		for cst in self.parent.dom.cosets:
			k = tuple(cst)
			new_vals[k] = [[0 for i in range(self.parent.wt + 1)] for j in range(self.parent.wt + 1)]
			for i in range(self.parent.wt + 1):
				for j in range(self.parent.wt + 1):
					new_vals[k][i][j] = self.values[k][i][j] + other.values[k][i][j]
		
		return ClassicalModSymb(self.parent, new_vals)
		
	def __rmul__(self, const):
		new_vals = {}
		
		for cst in self.parent.dom.cosets:
			key = tuple(cst)
			new_vals[key] = [[0 for i in range(self.parent.wt + 1)] for j in range(self.parent.wt + 1)]
			for i in range(self.parent.wt + 1):
				for j in range(self.parent.wt + 1):
					new_vals[key][i][j] = const * self.values[key][i][j]
		
		return ClassicalModSymb(self.parent, new_vals)
		
