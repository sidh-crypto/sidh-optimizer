# Copyright (c) 2016 Luca De Feo.

from collections import OrderedDict, namedtuple

class Cost():
    '''
    A class representing operation counts
    '''
    def __init__(self, from_cost=None, **kwds):
        self.costs = OrderedDict()
        if from_cost is not None:
            self.costs.update(from_cost.costs)
        for k, v in kwds.items():
            if k in self.costs:
                self.costs[k] += v
            else:
                self.costs[k] = v

    def __str__(self):
        return " + ".join("%d%s" % (v,k) for k, v in self.costs.items())
    
    __repr__ = __str__

    def __add__(self, other):
        assert type(other) == type(self)
        return Cost(other, **self.costs)

    def __mul__(self, scalar):
        assert type(scalar) == int
        return Cost(**{ k:scalar*v for (k,v) in self.costs.items() })

    __rmul__ = __mul__
    
    def weigh(self, *args, **kwds):
        '''
        Get a single numerical weight for this Cost by weighting with
        the given weights.
        '''
        weights = { k:1 for k in args }
        weights.update(kwds)
        return sum(self.costs.get(k, 0)*v for k, v in weights.items())

Party = namedtuple('Party', ['mul', 'isogeny', 'next_curve'])

# Arithmetic in GF(p^2), in terms of number of ops in GF(p).
#
# We use the following elementary operations:
#
# - i: inversion mod p of an integer of size log(p)
# - m: multiplication of integers of size log(p)
# - a: addition of integers of size log(p)
# - mod: reduction mod p of integers of size 2log(p)
#
# Assuming -1 is not a square in GF(p), we have the following costs
# for operations in GF(p^2)
A = Cost(a=2)                   # Addition
M = Cost(m=3, a=4, mod=2)       # Multiplication
m = Cost(m=2, mod=2)            # Scalar multipliccation
S = Cost(m=2, a=3, mod=2)       # Squaring
I = Cost(i=1, m=4, a=1, mod=2)  # Inversion

# DJP paper (original)

double = 3*M + 2*S + 4*A
triple = double + 4*M + 2*S + 4*A

DJP = {
    2: Party(double,
                 2*M + S + A,
                 I + 4*M + S + A),
    3: Party(triple,
                 4*M + 2*S + 2*A,
                 I + 4*M +   S + 4*A + m),
    4: Party(2*double,
                 6*M + S + 7*A,
                 double + 2*I + 6*M + 6*A + m),
}

# CLN paper (Microsoft)

CLN = {
    3: Party(triple + M,
                 6*M + 2*S + 2*A,
                 3*M + 3*S + 8*A),
    4: Party(2*double + 2*M,
                 9*M + S + 6*A,
                 5*S + 7*A),
}
