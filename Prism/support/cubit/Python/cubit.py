#
# Just a stub to initialize the cubit (adaptive spectral element) library
#
# $Id$

import _cubit
import Mesh
import Field


Mesh  = Mesh.Mesh
Field = Field.Field


# -----------------------------------------------------------------------------

if __name__=='__main__':
    
    from Mesh  import Mesh
    from Field import Field

    # Create a mesh by reading its description in a file
    mesh = Mesh('q1')
    print mesh

    # Create a field using this mesh and evaluate the integral of several 
    # functions using the "integrate" method.
    u = Field(mesh)
    for expr in ('1', '2', 'sin(x)', 'x', 'x^2'):
	u.set(expr)
	print "Integrate[%s] = %g" % (expr, u.integral())


    
