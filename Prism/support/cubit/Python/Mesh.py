#
# Mesh class
#
# $Id$

import _cubit, os

class Mesh:

    def __init__(self, fname=None):
	self._mesh = _cubit.build_mesh()
	if fname != None:
	    self.read(fname)
	return

    def read(self,fname):
	if os.path.exists(fname):
	    self._mesh.read(fname)
	elif os.path.exists(fname+'.rea'):
	    self._mesh.read(fname+'.rea')
	else:
	    raise IOError('cannot read mesh from %s or %s.rea' % (fname,fname))
	return

    def info(self):
	self._mesh.info()
	return

    def connect(self):
	self._mesh.connect()
	return

    def __repr__(self):
	return repr(self._mesh)

	
