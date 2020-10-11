import numpy

def calc_dihedral(c1, c2, c3, c4):
    
    v1 = Vector(c1[0], c1[1], c1[2])
    v2 = Vector(c2[0], c2[1], c2[2])
    v3 = Vector(c3[0], c3[1], c3[2])
    v4 = Vector(c4[0], c4[1], c4[2])
    
    ab = v1 - v2
    cb = v3 - v2
    db = v4 - v3
    u = ab ** cb
    v = db ** cb
    w = u ** v
    angle = u.angle(v)
    # Determine sign of angle
    try:
        if cb.angle(w) > 0.001:
            angle = -angle
    except ZeroDivisionError:
        # dihedral=pi
        pass
    return round(angle,2)

class Vector(object):
    "3D vector"

    def __init__(self, x, y=None, z=None):
        if y is None and z is None:
            # Array, list, tuple...
            if len(x) != 3:
                raise ValueError("Vector: x is not a "
                                "list/tuple/array of 3 numbers")
            self._ar = numpy.array(x, 'd')
        else:
            # Three numbers
            self._ar = numpy.array((x, y, z), 'd')

    def __repr__(self):
        x, y, z = self._ar
        return "<Vector %.2f, %.2f, %.2f>" % (x, y, z)

    def __neg__(self):
        "Return Vector(-x, -y, -z)"
        a = -self._ar
        return Vector(a)

    def __add__(self, other):
        "Return Vector+other Vector or scalar"
        if isinstance(other, Vector):
            a = self._ar + other._ar
        else:
            a = self._ar + numpy.array(other)
        return Vector(a)

    def __sub__(self, other):
        "Return Vector-other Vector or scalar"
        if isinstance(other, Vector):
            a = self._ar - other._ar
        else:
            a = self._ar - numpy.array(other)
        return Vector(a)

    def __mul__(self, other):
        "Return Vector.Vector (dot product)"
        return sum(self._ar * other._ar)

    def __div__(self, x):
        "Return Vector(coords/a)"
        a = self._ar / numpy.array(x)
        return Vector(a)

    def __pow__(self, other):
        "Return VectorxVector (cross product) or Vectorxscalar"
        if isinstance(other, Vector):
            a, b, c = self._ar
            d, e, f = other._ar
            c1 = numpy.linalg.det(numpy.array(((b, c), (e, f))))
            c2 = -numpy.linalg.det(numpy.array(((a, c), (d, f))))
            c3 = numpy.linalg.det(numpy.array(((a, b), (d, e))))
            return Vector(c1, c2, c3)
        else:
            a = self._ar * numpy.array(other)
            return Vector(a)

    def __getitem__(self, i):
        return self._ar[i]

    def __setitem__(self, i, value):
        self._ar[i] = value

    def __contains__(self, i):
        return (i in self._ar)

    def norm(self): 
        "Return vector norm" 
        return numpy.sqrt(sum(self._ar * self._ar)) 

    def left_multiply(self, matrix): 
        "Return Vector=Matrix x Vector" 
        a = numpy.dot(matrix, self._ar) 
        return Vector(a) 

    def angle(self, other):
        "Return angle between two vectors"
        n1 = self.norm()
        n2 = other.norm()
        c = (self * other) / (n1 * n2)
        # Take care of roundoff errors
        c = min(c, 1)
        c = max(-1, c)
        return numpy.arccos(c)/numpy.pi*180
    def normalize(self): 
        "Normalize the Vector" 
        self._ar = self._ar / self.norm() 

    def copy(self): 
        "Return a deep copy of the Vector" 
        return Vector(self._ar) 
    def get_array(self): 
        "Return (a copy of) the array of coordinates" 
        return numpy.array(self._ar)         
        
