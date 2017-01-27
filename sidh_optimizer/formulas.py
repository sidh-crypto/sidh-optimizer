# coding=utf-8

# Copyright (c) 2016-2017 Luca De Feo.

from collections import OrderedDict, namedtuple
import uuid
import sympy

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
        if not isinstance(other, Cost):
            return NotImplemented
        return Cost(other, **self.costs)

    def __neg__(self):
        return -1*self

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, scalar):
        if not isinstance(scalar, int):
            return NotImplemented
        return Cost(**{ k:int(scalar)*v for (k,v) in self.costs.items() })

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
m = Cost(m=2, mod=2)            # Scalar multiplication
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
    2: Party(double + M,
                 3*M + S + A,
                 3*M + S + 3*A),
    3: Party(triple + M,
                 6*M + 2*S + 2*A,
                 3*M + 3*S + 8*A),
    4: Party(2*double + 2*M,
                 9*M + S + 6*A,
                 5*S + 7*A),
}


########################
# Generic formulas

class Formula():
    def __init__(self, op, children, formula):
       self.op = op
       self.children = children
       self.formula = formula
       self._tag = None

    def __add__(self, other):
        if not isinstance(other, Formula):
            return NotImplemented
        return Formula('+', [self, other], self.formula + other.formula)

    def __sub__(self, other):
        if not isinstance(other, Formula):
            return NotImplemented
        return Formula('-', [self, other], self.formula - other.formula)
    
    def __mul__(self, other):
        if other is self:
            return Formula('^2', [self, other], self.formula**2)
        elif isinstance(other, Formula):
            return Formula('*', [self, other], self.formula * other.formula)
        elif type(other) == int:
            if other < 0:
                return (-other)*self
            elif other == 0:
                raise RuntimeError('No no-ops, please')
            elif other == 1:
                return self
            else:
                val = (self + self)*(other // 2)
                if other % 2:
                    val = val + self
                return val
        else:
            return NotImplemented

    __rmul__ = __mul__
        
    def __truediv__(self, other):
        if not isinstance(other, Formula):
            return NotImplemented
        return Formula('/', [self, other], self.formula / other.formula)

    def __pow__(self, exp):
        if not isinstance(exp, int):
            return NotImplemented
        if exp <= 0:
            raise NotImplementedError
        if exp == 1:
            return self
        else:
            val = (self * self)^(exp // 2)
            if exp % 2:
                val = val*self
            return val

    __xor__ = __pow__

    def __neg__(self):
        return Formula('u-', [self], -self.formula)

    def __repr__(self):
        return repr(self.formula)

    def cost(self, costs=None, tag=None):
        if costs is None:
            costs = { '+': Cost(A=1), '-': Cost(A=1), 'u-': Cost(),
                      '*': Cost(M=1), '^2': Cost(S=1), '/': Cost(M=1, I=1) }
        # Topological sort: count this node only if it has not been counted yet.
        # Soooo hacky!
        if tag is None:
            tag = uuid.uuid4()
        if self._tag == tag:
            return Cost()
        self._tag = tag
        return sum((c.cost(costs, tag) for c in self.children), costs[self.op])


class Var(Formula):
    def __init__(self, name, comment=''):
        self.formula = sympy.symbols(name)
        self.comment = comment

    def cost(self, *args, **kwds):
        return Cost()
    
class Const(Formula):
    def __init__(self, value):
        self.formula = sympy.numbers.sympify(value)

    def cost(self, *args, **kwds):
        return Cost()

def cost(*args, costs=None, tag=None):
    '''
    Compute cost of formulas, sharing common subexpressions.
    '''
    if tag is None:
        tag = uuid.uuid4()
    return sum((f.cost(costs, tag) for f in args), Cost())
    

##########################################
# Isogenies of Edwards curves
#
# To avoid inversions, we use the following generalized Edwards curves
#
#    (x² + y²)E = E + Dx²y²
#
# Setting d = D/E we recover the usual form for E. curves.
##########################################
def ed_double(X, Y, Z):
    '''
    <https://hyperelliptic.org/EFD/g1p/auto-edwards-projective.html#doubling-dbl-2007-bl>

    Opcount: 3M + 4S + 6A
    '''
    b = (X + Y)^2
    c = X^2
    d = Y^2
    e = c + d
    h = Z^2
    j = e - 2*h
    return (b-e)*j, e*(c-d), e*j

def ed_add(X1, Y1, Z1, X2, Y2, Z2, D, E):
    '''
    <https://hyperelliptic.org/EFD/g1p/auto-edwards-projective.html#addition-add-2007-bl>
    plus minor modifications

    Opcount: 13M + 1S + 7A
    '''
    a = Z1*Z2
    b = a^2
    ea = E*a
    eb = E*b
    c = X1*X2
    d = Y1*Y2
    e = D*c*d
    f = eb - e
    g = eb + e
    return ea*f*((X1 + Y1)*(X2 + Y2) - c - d), ea*g*(d-c), f*g

def edwards(ell, D, E, KX, KY, KZ, X, Y, Z):
    '''
    Compute generic formula for isogenies of Edwards curves of
    degree ell.

    Uses Sec. 4.3 of
    <https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/edwardsisogenies.pdf>

    D, E: curve parameters
    KX, KY, KZ: kernel generator
    X, Y, Z: function coordinates
    '''
    assert ell % 2 == 1
    
    ### Next curve
    
    # compute kernel points
    # Depending on ell, some additions may be replaced with doublings
    kernel = [(KX, KY, KZ)]
    for i in range(1, ell//2):
        tmpX, tmpY, tmpZ = kernel[-1]
        kernel.append(ed_add(tmpX, tmpY, tmpZ, KX, KY, KZ, D, E))

    # compute constants appearing in isogeny formula (3)
    D4, E4 = 1, 1
    constants = []
    for PX, PY, PZ in kernel:
        x2, y2, z2 = PX^2, PY^2, PZ^2
        y4, z4 = y2^2, z2^2
        D4, E4 = D4*y4, E4*z4
        tmp = z4 * E
        constants.append((y2*tmp, x2*tmp, x2*y4*D))

    ### Isogeny evaluation
    # This uses the constants defined in the previous phase,
    # named α, β, γ
    x, y, z = X, Y, Z
    X2, Y2, Z2 = X^2, Y^2, Z^2
    Z4 = Z2^2
    X2Y2 = X2 + Y2 - Z4
    sub = lambda n: ''.join(chr(ord(c)+ord('₀')-ord('0')) for c in str(n))
    for i, _ in enumerate(constants):
        a = c = Var('β' + sub(i))
        b = Var('α' + sub(i))
        d = Var('γ' + sub(i))
        e, f, g = (a - b)*(X2 + Y2), a*X2, b*Y2
        x *= f - g
        y *= e - f + g
        z *= c*Z4 - d*X2Y2
    Z2S = Z2^(ell // 2)
    x *= Z2S
    y *= Z2S

    return (D4**2 * D**ell, E4**2 * E**ell), constants, (x, y, z)

##########################################
# Isogenies of Montgomery curves
# <http://eprint.iacr.org/2011/506.pdf>
##########################################
def montgomery(ell):
    pass
