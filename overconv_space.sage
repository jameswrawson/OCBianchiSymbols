import fund_dom, overconv_symb, classical_space, utilities

class OCModSymbSpace:
	def __init__(self, field, level, wt, prime, prec): #currently assuming p split, and p must divide the level
		self.dom = fund_dom.FundDom(field, level)
		self.K = field
		self.N = level
		self.wt = wt
		self.prec = prec
		self.prime = prime
		if prime.ramification_index() != 1 or prime.residue_class_degree() != 1:
			raise ValueError("Inert or ramified primes are not yet supported")
		p = prime.norm()
		self.padics = Qp(p, prec, type = 'capped-rel')
		gen = prime.gens_reduced()[0]
		a, b = utilities.get_cpts(field, gen)
		poly = field.polynomial()
		c = poly.coefficients()[1]
		d = poly.coefficients()[0]
		disc = c**2 - 4*d
		root = self.padics(disc).square_root(extend=False)
		root = (-c + root)/2
		if p == 2: #re-introduce precision after square roots and dividing by 2
			root = (root**2 + d)/(-c)
			root = (root**2 + d)/(-c)
		target = a + b*root
		if target.valuation() > 0:
			self.gen = root
		else:
			self.gen = -root - c
		
		self.classical = classical_space.ModSymbSpace(field, level, wt)
		self.solved = False
	
	def convert(self, num_elt): #convert number field element into p-adic numbers
		(k, l) = utilities.get_cpts(self.K, self.K(num_elt))
		return (self.padics(k) + self.padics(l)*self.gen).add_bigoh(self.prec)
		
	def distAct(self, mat, dist):
		R.<x, y> = PowerSeriesRing(self.padics, default_prec = self.prec)
		a, b, c, d = mat[0][0], mat[0][1], mat[1][0], mat[1][1]
		abar, bbar, cbar, dbar = a.conjugate(), b.conjugate(), c.conjugate(), d.conjugate()
		a, b, c, d = self.convert(a), self.convert(b), self.convert(c), self.convert(d)
		abar, bbar, cbar, dbar = self.convert(abar), self.convert(bbar), self.convert(cbar), self.convert(dbar)
		out_dist = [[self.padics(0) for i in range(self.prec)] for j in range(self.prec)]
		for i in range(self.prec):
			for j in range(self.prec):
				if i + j >= self.prec:
					break
				f = (a + c*x)**self.wt * (abar + cbar*y)**self.wt * ((d*x + b) / (a + c*x))**i * ((dbar*y + bbar)/(abar + cbar*y))**j
				coeffs = f.coefficients()
				for k in range(self.prec):
					for l in range(self.prec):
						if x**k * y**l in coeffs:
							out_dist[i][j] += dist[k][l] * coeffs[x**k * y**l]
				
				if i > self.wt or j > self.wt:
					res = max(0, self.prec - i - j)
					out_dist[i][j] = out_dist[i][j].add_bigoh(res)
		return out_dist
	
	def createModularSymbol(self, values):
		return overconv_symb.OverconvSymb(self, values)
	
	def getClassicalBasis(self):
		return self.classical.getBasis()
	
	def liftClassical(self, form, e_val):
		data = {}
		for cst in self.classical.dom.cosets:
			dat = [[self.padics(0) for i in range(self.prec)] for j in range(self.prec)]
			for i in range(self.prec):
				for j in range(self.prec):
					if i <= self.wt and j <= self.wt:
						dat[i][j] = self.convert(form.values[tuple(cst)][i][j])
				
			data[tuple(cst)] = dat
		
		symb = self.createModularSymbol(data)
		for i in range(self.prec):
			prim = self.prime.gens_reduced()[0]
			symb = symb.hecke(prim)
			symb = (1 / self.convert(e_val)) * symb.hecke(prim.conjugate())
		
		return symb
	
	def loadModularSymbol(self, filename, appended=True):
		f = open(filename, "r")
		line = ""
		data = {}
		if appended:
			while line != "#"*10 + "\n":
				f.readline()
		
		line = f.readline()
		while line:
			key, value = line[:-1].split(" : ")
			
			dat = []
			rows = value[1:-1].split("[, ")
			for row in rows:
				vals = row[1:].split(", ")
				out_row = []
				for val in vals:
					out_row.append(self.padics(val))
				
				dat.append(out_row)
			
			data[tuple(cst)] = dat
			line = f.readline()
		
		return createModularSymbol(data)
	
	def getBasis(self):
		if self.solved:
			return self.basis
		
		relns = []
		for cst in self.classical.dom.cosets:
			for i in range(self.prec):
				for j in range(self.prec):
					if i + j >= self.prec:
						continue
					
					values = []
					for reln in self.classical.dom.relations:
						totals = [[0 for k in range(self.prec)] for l in range(self.prec)]
						for term in reln.data:
							if term[1] == cst:
								dist = [[0 for k in range(self.prec)] for l in range(self.prec)]
								dist[i][j] = 1.add_bigoh(self.prec - i - j)
								new_dist = self.distAct(term[0].A.inverse(), dist)
								for k in range(self.prec):
									for l in range(self.prec):
										totals[k][l] += term[0].n * new_dist[k][l]
						
						for k in range(self.prec):
							for l in range(self.prec):
								res = self.prec
								if i > self.wt or j > self.wt:
										res = max(0, self.prec - i - j)
								values.append(totals[k][l].add_bigoh(res))
					relns.append(values)
		
		self.relns = matrix(self.K, relns, sparse=True)
		
		basis = self.relns.kernel().basis()
		self.solved = True
		self.basis = []
		
		for elt in basis:
			data = []
			for i in range(len(self.classical.dom.cosets)):
				data.append((tuple(self.dom.cosets[i]), 
				[[elt[i*(self.prec)**2 + self.prec*j + k] for k in range(self.prec)] for j in range(self.prec)]))
			symbol = self.createModularSymbol(dict(data))
			self.basis.append(symbol)
		
		self.solved = True
		
		return self.basis
