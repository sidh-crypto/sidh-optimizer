# coding=utf-8

# Copyright (c) 2016 Luca De Feo.

import math

class Path():
    r""" A path represents what is called a "strategy" in the SIDH
    paper. In practice, it is a tree-like structure looking like this:

       ╱╲   
      ╱╲ ╲  
     ╱ ╱  ╲ 
    ╱╲ ╲ `╱╲
    """
 
    def __init__(self, path):
        r"""
        To create a path,
        
        >>> Path(list_of_integers)
        
        where `list_of_integers` is a list coding each level of the
        tree on the bits of an integer. This constructor is meant for
        internal use.

        Caveat: the bits of the integers, read from right to left,
        determine the edges of the path from left to right.
        """
        n = 4
        for floor in path:
            if floor >= n:
                raise ValueError("Malformed path")
            n <<= 2
        self.path = path

    def ascii_art(self, ascii=False, chars=None):
        "ASCII-art representation of the path."
        n = len(self.path)
        s = ''
        if chars is not None:
            l, r, lr, nil = chars
        elif ascii:
            l, r, lr, nil = ('/ ', ' \\', '/\\', ' `')
        else:
            l, r, lr, nil = ('╱ ', ' ╲', '╱╲', ' `')
        for i, floor in enumerate(self.path):
            s += (n - i - 1) * ' '
            for j in range(i+1):
                edges = floor & 3
                if edges == 0:
                    s += nil
                elif edges == 1:
                    s += l
                elif edges == 2:
                    s += r
                else:
                    s += lr
                floor >>= 2
            s += (n - i - 1) * ' ' + '\n'
        return s

    def xpm(self):
        '''
        PBM image bitmap representation of the path.
        '''
        n = len(self.path)
        return """! XPM2
{} {} 2 1
. m black
  m None
""".format(2*n, n) + self.ascii_art(chars=('. ', ' .', '..', '  '))

    def bitmap(self):
        """
        Return a bitmap image representation of this strategy.

        Use `s.bitmap.save('strategy.png')` to save to a png file.
        """
        try:
            from PIL import Image
        except ImportError:
            raise NotImplementedError("Please install the PIL library")
        n = len(self.path)
        bytes = self.ascii_art(chars=('. ', ' .', '..', '  ')).replace('\n', '').encode()
        assert len(bytes) == 2*n**2
        im = Image.frombytes('P', (2*n, n), bytes)
        palette = 768*[255]
        palette[3*ord('.'):3*ord('.')+3] = [0]*3
        im.putpalette(palette)
        return im
    
    def __repr__(self):
        return "\n" + self.ascii_art()

    def __hash__(self):
        return hash(sum(f<<(2*i) for i,f in enumerate(self.path)))

    def __eq__(self, other):
        return self.path == other.path

    def height(self):
        """
        The height of a path is the number of levels, i.e. one plus
        the number of edges.
        """
        return len(self.path) + 1

    def count(self):
        """
        Count the number of left and right edges in the path and
        return it as a tuple.
        """
        left = 0
        right = 0
        for floor in self.path:
            while floor > 0:
                left += floor & 1
                floor >>= 1
                right += floor & 1
                floor >>= 1
        return left, right

    def crosses(self):
        """
        Returns False if the path is non-crossing (in the sense of the
        paper).
        """
        for floor in self.path:
            floor >>= 1
            while floor > 0:
                if floor & 3 == 3:
                    return True
                floor >>= 2
        return False

    def well_formed(self):
        """
        Returns True if the path is well-formed (in the sense of the
        paper).
        """
        reachable = 2**self.height() - 1
        for floor in reversed(self.path):
            i = 1
            newreachable = 0

            r = reachable & i
            if r != (floor & 1):
                return False
            reachable ^= r
            newreachable |= r
            floor >>= 1

            while floor > 0:
                i <<= 1
                edge = floor & 3
                r = reachable & i
                if r:
                    if edge == 1:
                        newreachable |= i >> 1
                    elif edge == 2:
                        newreachable |= i
                    else:
                        return False
                elif edge != 0:
                        return False
                reachable ^= r
                floor >>= 2
                
            if reachable > 0:
                return False
            reachable = newreachable
        return True
            

    def cat(self, floor):
        "Add one level to the bottom of the path. For internal use."
        if floor >= 4**self.height():
            raise ValueError("Malformed path: floor %s is too large" % floor)
        return Path(self.path + [floor])

    def __mul__(self, other):
        """
        Multiplication of two paths: construct the path with minimal
        number of edges having `self` as left sub-path and `other` as
        right sub-path.
        """
        sh = self.height()
        oh = other.height()
        newpath = []
        if sh >= oh:
            h, H = oh, sh
        else:
            h, H = sh, oh
            
        for i in range(h):
            newpath.append(1 | (1 << (2*i + 1)))
            
        lit = iter(self.path)
        rit = iter(other.path)

        for i in range(h, H):
            if sh > oh:
                newpath.append(next(lit) | (1 << (2*i + 1)))
            else:
                newpath.append(1 | (next(rit) << (2*sh)))

        for l, r in zip(lit, rit):
            newpath.append(l | (r << (2*sh)))
            
        return Path(newpath)
        


def paths(n):
    """
    Construct all paths of height `n`. Even non well-formed ones.
    """
    if n <= 0:
        raise ValueError
    elif n == 1:
        return [Path([])]
    else:
        subpaths = paths(n-1)
        return [p.cat(j) for p in subpaths for j in range(4**(n-1))]

def wf_paths(reachable):
    """
    Construct all well-formed paths satisfying a given condition.

    The condition is as follows: all the paths have height equal to
    the ceiling of log_2(`reachable` + 1). `reachable` is interpreted
    as a bitfield, with 1 meaning that the corresponding leaf on the
    floor of the path should be reachable from the root, 0 meaning the
    opposite.

    This function has been used to count well-formed paths and guess
    the link with Gelfand-Zetlin polytopes.
    """
    if reachable <= 0:
        raise ValueError
    elif reachable == 1:
        return [Path([])]
    else:
        floors = [reachable & 1]
        reachable >>= 1

        left = 2; right = 4
        while reachable > 1:
            if reachable & 1:
                floors = [f | left for f in floors] + [f | right for f in floors]
            left <<= 2; right <<= 2
            reachable >>= 1

        floors = [f | left for f in floors]
        
        paths = []
        for f in floors:
            paths.extend([p.cat(f) for p in wf_paths(_h4(f))])

        return paths


def _h4(n):
    """
    Utility function on bitfields.
    """
    res = 0
    i = 1
    while n > 0:
        if n & 3 != 0:
            res |= i
        i <<= 1
        n >>= 2

    return res


def optimal_paths(n, l, r, construct=True):
    """
    Compute optimal paths (in the sense of the paper) of height up to
    `n`, with the cost of a left (resp. right) edge being given by `l`
    (resp. `r`).
    
    If `construct` is True, the output is a list of pairs cost-path,
    with `output[i]` holding the optimal pair of height `i`
    (`output[0]` contains nothing interesting).

    If `construct` is False, the output is a list of pairs of
    integers. The first integer is the same cost as for the
    `construct` case, the second integer is the height of the topmost
    branching in the left branch of the optimal path (so that the
    sub-path starting at that point is given by the optimal path of
    that height).
    """
    if n <= 0:
        raise ValueError
    else:
        if l > r:
            l, r = r, l
            swapped = True
        else:
            swapped = False

        opaths = [(0,0), (0,1)]
        for i in range(2, n+1):
            opaths.append(((i + 1)**2*(l+r), 0))
            for j in range(1, i//2 + 1):
                lscore = opaths[j][0]
                rscore = opaths[i-j][0]
                score = lscore + rscore + (i-j)*l + j*r
                if score <= opaths[i][0]:
                    opaths[i] = (score, j)

        if swapped:
            opaths = [(x,i-j) for (i,(x,j)) in enumerate(opaths)]
                    
        if construct:
            opaths[1] = (0, Path([]))
            for i in range(2,len(opaths)):
                j = opaths[i][1]
                opaths[i] = (opaths[i][0], opaths[j][1] * opaths[i-j][1])

        return opaths
