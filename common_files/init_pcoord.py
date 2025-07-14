import mdtraj
import numpy

parent = mdtraj.load('basis.pdb', top='basis.pdb')
dist1_parent = mdtraj.compute_distances(parent, [[1965,4239]], periodic=True)
dist2_parent = mdtraj.compute_distances(parent, [[4260,681]], periodic=True)
angle1_parent = mdtraj.compute_angles(parent, [[4150, 4233, 4260]], periodic=True)

d2_arr = numpy.asarray(dist2_parent)
d2_arr = d2_arr*10

d1_arr = numpy.asarray(dist1_parent)
d1_arr = d1_arr*10

a1_arr = numpy.asarray(angle1_parent)
a1_arr = numpy.degrees(a1_arr)

numpy.savetxt("pcoord2.dat", d1_arr)
numpy.savetxt("pcoord1.dat", d2_arr)
numpy.savetxt("pcoord3.dat", a1_arr)
