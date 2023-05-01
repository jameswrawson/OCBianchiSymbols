import utilities, group_ring
import itertools

class FundDom: #Fundamental domain for Gamma_0(N)
	def __init__(self, field, level): 
		self.K = field
		self.N = level
		disc = self.K.disc()
		prime_facs = field.ideal(level).factor()
		self.primes = []
		for (p, n) in prime_facs:
			self.primes.append(p**n)
		
		identity = matrix(self.K, 2, 2, [1, 0, 0, 1])
		S = matrix(self.K, 2, 2, [0, -1, 1, 0])
		T = matrix(self.K, 2, 2, [1, -1, 1, 0])
		gamma_inf = group_ring.ZGElt(field, level, [[group_ring.ZGCoeff(1, identity), [0] * len(self.primes)]])
		basic_relns = [
			gamma_inf + group_ring.ZGCoeff(1, S) * gamma_inf,
			gamma_inf + group_ring.ZGCoeff(1, T) * gamma_inf + group_ring.ZGCoeff(1, T*T) * gamma_inf
		]
				
		if disc == -4:
			i = self.K(-1).sqrt(extend=False)
			X = matrix(self.K, 2, 2, [i, 1, 1, 0])
			basic_relns.append(gamma_inf + group_ring.ZGCoeff(1, X) * gamma_inf + group_ring.ZGCoeff(1, X*X) * gamma_inf)
			J = matrix(self.K, 2, 2, [1, 0, 0, i])
			basic_relns.append(gamma_inf + group_ring.ZGCoeff(-1, J) * gamma_inf)
		
		elif disc == -8:
			a = self.K(-2).sqrt(extend=False)
			X = matrix(self.K, 2, 2, [a, 1, 1, 0])
			basic_relns.append(gamma_inf + group_ring.ZGCoeff(1, X) * gamma_inf + group_ring.ZGCoeff(1, X*X) * gamma_inf + group_ring.ZGCoeff(1, X*X*X) * gamma_inf)
			J = matrix(self.K, 2, 2, [1, 0, 0, -1])
			basic_relns.append(gamma_inf + group_ring.ZGCoeff(-1, J) * gamma_inf)
		
		
		elif disc == -3:
			a = (1+self.K(-3).sqrt(extend=False))/2
			X = matrix(self.K, 2, 2, [1, a*a, a, 0])
			basic_relns.append(gamma_inf + group_ring.ZGCoeff(1, X) * gamma_inf + group_ring.ZGCoeff(1, X*X) * gamma_inf)
			J = matrix(self.K, 2, 2, [1, 0, 0, a])
			basic_relns.append(gamma_inf + group_ring.ZGCoeff(-1, J) * gamma_inf)
		
		elif disc == -7:
			a = (1+self.K(-7).sqrt(extend=False))/2
			X = matrix(self.K, 2, 2, [-1, a, a-1, 1])
			U = matrix(self.K, 2, 2, [1, a, 0, 1])
			basic_relns.append(gamma_inf + group_ring.ZGCoeff(-1, U) * gamma_inf + group_ring.ZGCoeff(1, X) * gamma_inf + group_ring.ZGCoeff(-1, X*U) * gamma_inf)
			J = matrix(self.K, 2, 2, [1, 0, 0, -1])
			basic_relns.append(gamma_inf + group_ring.ZGCoeff(-1, J) * gamma_inf)
		
		elif disc == -11:
			a = (1+self.K(-11).sqrt(extend=False))/2
			X = matrix(self.K, 2, 2, [-1, a, a-1, 2])
			U = matrix(self.K, 2, 2, [1, a, 0, 1])
			basic_relns.append(gamma_inf + group_ring.ZGCoeff(-1, U) * gamma_inf + group_ring.ZGCoeff(1, X) * gamma_inf + group_ring.ZGCoeff(-1, X*U) * gamma_inf + group_ring.ZGCoeff(1, X*X) * gamma_inf + group_ring.ZGCoeff(-1, X*X*U) * gamma_inf)
			J = matrix(self.K, 2, 2, [1, 0, 0, -1])
			basic_relns.append(gamma_inf + group_ring.ZGCoeff(-1, J) * gamma_inf)
		
		else:
			raise ValueError("Non-Euclidean fields not yet supported")
		
		reses = [list(prime.residues()) + ["inf"] for prime in self.primes]
		all_reps = itertools.product(*reses)
		self.cosets = [list(rep) for rep in all_reps]
		self.relations = []
		for coset in self.cosets:
			mat = basic_relns[0].get_rep(coset)
			for reln in basic_relns:
				self.relations.append(group_ring.ZGCoeff(1, mat) * reln)
	
	def path(self, mat):
		start = 0
		if mat[1][0] == 0:
			start = "inf"
		else:
			start = mat[0][0] / mat[1][0]
		
		end = 0
		if mat[1][1] == 0:
			end = "inf"
		else:
			end = mat[0][1] / mat[1][1]
		
		terms = group_ring.ZGElt(self.K, self.N, [])
		identity = matrix(self.K, 2, 2, [1, 0, 0, 1])
		gamma_inf = group_ring.ZGElt(self.K, self.N, [[group_ring.ZGCoeff(1, identity), [0] * len(self.primes)]])
		if start != "inf":
			fracs = utilities.ctd_frac(self.K, start)
			for frac in fracs:
				terms = terms + group_ring.ZGCoeff(-1, frac) * gamma_inf
		
		if end != "inf":
			fracs = utilities.ctd_frac(self.K, end)
			for frac in fracs:
				terms = terms + group_ring.ZGCoeff(1, frac) * gamma_inf
		
		return terms
