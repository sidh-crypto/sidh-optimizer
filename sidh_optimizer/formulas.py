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

# Arithmetic in GF(p^2), in terms of number of ops in GF(p)
A = Cost(a=2, mod=2)
M = Cost(m=3, a=4, mod=2)
S = Cost(m=2, a=3, mod=2)
I = Cost(i=1, m=4, a=1, mod=2)

# DJP paper (original)

DJP = {
    2: Party(3*M + 2*S, 2*M + S, Cost()),
    3: Party(7*M + 4*S, 4*M + 2*S, Cost()),
    4: Party(6*M + 4*S, 6*M + S, Cost()),
}

# CLN paper (Microsoft)
CLN = {}
