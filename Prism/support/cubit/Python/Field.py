#
# Field class
#
# $Id$

import _cubit, sys

class Field:

    def __init__(self, mesh):
	self._field = _cubit.build_field(mesh._mesh)
	self.mesh   = mesh
	self.type   = 'u'
	return

    # Computational methods
    def parse(self,expr):
	self._field.parse(expr)
	return self

    # Scalar reduction
    def min(self):
	return self._field.min()
    def max(self):
	return self._field.max()
    def amax(self):
	return self._field.amax()
    def integral(self):
	return self._field.integral()
    def nrm2(self):
	return self._field.nrm2()
    def L2(self):
	return self._field.L2()
    def H1(self):
	return self._field.H1()

    # Computational methods (all return new Field objects)
    def dx(self):
	return self.gradient(0)
    def dy(self):
	return self.gradient(1)
    def dz(self):
	return self.gradient(2)
    def gradient(self,dir):
	res = Field(self.mesh)
	res._field = self._field.gradient(dir)
	return res
    def FFT(self, dir):
	res = Field(self.mesh)
	res._field = self._field.FFT(dir)
	return res
    def copy(self):
	res = Field(self.mesh)
	res._field = self._field.copy()
	return res
    def scal(self,d):
	res = Field(self.mesh)
	res._field = self._field.scal(d)
	return res
    def shift(self,d):
	res = Field(self.mesh)
	res._field = self._field.shift(d)
	return res

    def __repr__(self):
	return repr(self._field)
    def __add__(self,other):
	res = Field(self.mesh)
	res._field = self._field + other._field
	return res
    def __sub__(self,other):
	res = Field(self.mesh)
	res._field = self._field - other._field
	return res
    def __mul__(self,other):
	res = Field(self.mesh)
	res._field = self._field * other._field
	return res
    def __div__(self,other):
	res = Field(self.mesh)
	res._field = self._field / other._field
	return res
    def __pos__(self):
	res = Field(self.mesh)
	res._field = self._field
	return res
    def __neg__(self):
	res = Field(self.mesh)
	res._field = -self._field
	return res
    def __abs__(self):
	res = Field(self.mesh)
	res._field = abs(self._field)
	return res
    def __pow__(self,other):
	res = Field(self.mesh)
	res._field = pow(self._field,other)
	return res

# -----------------------------------------------------------------------------

def test_reduction(u):
    print '\nScalar reduction:'
    print 'u.min      = ', u.min()
    print 'u.max      = ', u.max()
    print 'u.amax     = ', u.amax()
    print 'u.integral = ', u.integral()
    print 'u.nrm2     = ', u.nrm2()
    print 'u.L2       = ', u.L2()
    print 'u.H1       = ', u.H1()
    return

def test_binary(u,v):
    print '\nBinary operations:'
    print 'u =', u.max()
    print 'v =', v.max()
    w = u + v
    print 'w = u + v =', w.max()
    w = u - v
    print 'w = u - v =', w.max()
    w = u * v
    print 'w = u * v =', w.max()
    w = u / v
    print 'w = u / v =', w.max()
    return

def test_misc(u):
    print '\nMisc operations:'

    print 'v = u.copy()'
    v = u.copy()
    print v

    print 'v = 0.5*u + v [v = u.scal(0.5)+v]'
    v = u.scal(0.5) + v
    print v

    print 'v = 4.5 + u [v = u.shift(4.5)]'
    v = u.shift(4.5)
    print v

    print 'v = 0.5*u + 1 [v = u.scal(0.5).shift(1.)]'
    v = u.scal(0.5).shift(1.)
    print v
    
    print 'w = v.dx()-u.dy() [where v = x and u = -y]'
    w = v.parse('x').dx()-u.parse('-y').dy()
    print w

    return
    
if __name__=='__main__':

    from Mesh import Mesh

    mesh = Mesh('q1')

    u = Field(mesh)
    v = Field(mesh)

    test_reduction (u.parse('x'))
    test_binary    (u.parse('1'), v.parse('2'))
    test_misc      (u)
