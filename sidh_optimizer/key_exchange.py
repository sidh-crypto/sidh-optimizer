# coding=utf-8

# Copyright (c) 2016 Luca De Feo.

from .paths import optimal_paths

def opcount(n, costs, weights=None):
    """
    Determine the total operation count for the optimal strategy
    determined by the given costs with the given weights.
    """
    if weights is None:
        w = (1,1)
    else:
        w = (costs.mul.weigh(**weights), costs.isogeny.weigh(**weights))
    cost, path = optimal_paths(n, *w)[-1]
    l, r =  path.count()
    total = l*costs.mul + r*costs.isogeny + n*costs.next_curve
        
    return path, total, cost
