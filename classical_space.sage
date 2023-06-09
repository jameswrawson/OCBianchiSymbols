import fund_dom, classical_modsymb

class ModSymbSpace:
	def __init__(self, field, level, wt):
		self.dom = fund_dom.FundDom(field, level)
		self.K = field
		self.N = level
		self.wt = wt
		
		relns = []
		for cst in self.dom.cosets:
			for i in range(wt + 1):
				for j in range(wt + 1):
					values = []
					for reln in self.dom.relations:
						totals = [[0 for k in range(wt + 1)] for l in range(wt + 1)]
						for term in reln.data:
							if term[1] == cst:
								dist = [[0 for k in range(wt + 1)] for l in range(wt + 1)]
								dist[i][j] = 1
								new_dist = self.distAct(term[0].A.inverse(), dist)
								for k in range(wt + 1):
									for l in range(wt + 1):
										totals[k][l] += term[0].n * new_dist[k][l]
						
						for k in range(wt + 1):
							for l in range(wt + 1):
								values.append(totals[k][l])
					relns.append(values)
		self.relns = matrix(self.K, relns, sparse=True) #transpose of the expected matrix
		self.solved = False
	
	def distAct(self, mat, dist):
		R = PolynomialRing(self.K, 2, 'xy')
		(x, y) = R.gens()
		f = 0
		for i in range(self.wt + 1):
			for j in range(self.wt + 1):
				f += dist[i][j] * binomial(self.wt, i) * binomial(self.wt, j) * x^(self.wt - i) * y^(self.wt - j)
		
		a, b, c, d = mat[0][0], mat[0][1], mat[1][0], mat[1][1]
		abar, bbar, cbar, dbar = a.conjugate(), b.conjugate(), c.conjugate(), d.conjugate()
		g = (c*x + d)**self.wt * (cbar * y + dbar)**self.wt * f((a*x + b)/(c*x + d), (abar*y + bbar)/(cbar*y + dbar))
		denom = g.denominator()
		g = g.numerator()
		
		out = [
			[
				g.coefficient({x : self.wt - i, y : self.wt - j})/ (binomial(self.wt, i) * binomial(self.wt, j) * denom)
				for j in range(self.wt + 1)
			]
			for i in range(self.wt + 1)
		]
		
		return out
	
	def createModularSymbol(self, values):
		return classical_modsymb.ClassicalModSymb(self, values)
	
	def loadModularSymbol(self, filename):
		f = open(filename, "r")
		line = ""
		data = {}
		line = f.readline()
		while line != "#"*10 + "\n":
			if line == "":
				break
			
			key, value = line[:-1].split(" : ")
			cst = []
			reps = key[1:-1].split(", ")
			for rep in reps:
				if rep == "'inf'":
					cst.append('inf')
				else:
					cst.append(self.K(rep))
			
			rows = value[1:-1].split("], ")
			dat = []
			for row in rows:
				vals = row[1:].split(", ")
				out_row = []
				for val in vals:
					out_row.append(self.K(val))
				
				dat.append(out_row)
			
			data[tuple(cst)] = dat
			line = f.readline()
		
		return self.createModularSymbol(data)
	
	def dimension(self):
		if self.solved:
			return len(self.basis)
		else:
			self.basis = self.relns.kernel().get_basis()
			self.solved = True
			return len(self.basis)
	
	def getBasis(self):
		if not self.solved:
			self.space = self.relns.kernel()
			self.basis = self.space.basis()
			self.solved = True
		
		mod_basis = []
		for elt in self.basis:
			data = []
			for i in range(len(self.dom.cosets)):
				data.append(
					(tuple(self.dom.cosets[i]),
					[
					[elt[i*(self.wt + 1)**2  + (self.wt + 1)*j + k] for k in range(self.wt + 1)] for j in range(self.wt + 1)
					]))
			symbol = self.createModularSymbol(dict(data))
			mod_basis.append(symbol)
		return mod_basis
