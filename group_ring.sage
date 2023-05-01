import copy
import itertools
import utilities

class ZGCoeff:
	def __init__(self, n, A):
		self.n = n
		self.A = A
	
class ZGElt:
	def __init__(self, K, level, data):
		self.K = K
		self.N = level
		self.prime_fac = K.ideal(level).factor()
		self.primes = []
		for (p, n) in self.prime_fac:
			self.primes.append(p**n)
		
		self.data = data
	
	def __add__(self, other):
		new_data = copy.deepcopy(self.data) + copy.deepcopy(other.data)
		out = ZGElt(self.K, self.N, new_data)
		out.simplify()
		return out
	
	def __rmul__(self, coeff):
		new_data = []
		for term in self.data:
			mat = self.get_rep(term[1])
			
			new_mat = coeff.A * term[0].A * mat
			new_coset = self.get_coset(new_mat)
			new_rep =  self.get_rep(new_coset)
			new_coeff = new_mat * new_rep.inverse()
			
			new_data.append([ZGCoeff(term[0].n * coeff.n, new_coeff), new_coset])
			
		out = ZGElt(self.K, self.N, new_data)
		out.simplify()
		return out
	
	def __repr__(self):
		output = ""
		for term in self.data:
			output += "{0} * {1} * D{2} + ".format(term[0].n, term[0].A, term[1])
		
		if self.data == []:
			output = "0"
		else:
			output = output[:-3]
		
		return output
	
	def simplify(self):
		i = 0
		while i < len(self.data):
			j = i + 1
			while j < len(self.data):
				if self.data[i][1] == self.data[j][1]:
					if self.data[i][0].A == self.data[j][0].A or self.data[i][0].A == -self.data[j][0].A:
						self.data[i][0].n += self.data[j][0].n
						del self.data[j]
					else:
						j += 1
				else:
					j += 1
				
			i += 1
		
		i = 0
		while i < len(self.data):
			if self.data[i][0].n == 0:
				del self.data[i]
			else:
				i += 1
	
	def get_rep(self, coset):
		cs = []
		ds = []
		for c in coset:
			if c != "inf":
				cs.append(c)
				ds.append(1)
			else:
				cs.append(1)
				ds.append(0)
		
		c = self.K.idealchinese(self.primes, cs)
		d = self.K.idealchinese(self.primes, ds)
		
		R = self.K.maximal_order()
		level_ideal = self.K.ideal(self.N)
		while R(c).gcd(R(d)) != 1:
			if c in level_ideal:
				d += self.N
			else:
				c += self.N
		
		u, v, gcd = utilities.bezout(self.K, c, d)
		if gcd != 1:
			u = u / gcd
			v = v / gcd
		
		return matrix(self.K, 2, 2, [v, -u, c, d])
		
	def get_coset(self, elt):
		coset = []
		level_ideal = self.K.ideal(self.N)
		for prime in self.primes:
			if elt[1][1] in prime:
				coset.append("inf")
			else:
				coset.append(prime.reduce(elt[1][0] * elt[1][1].inverse_mod(prime)))
		
		return coset
