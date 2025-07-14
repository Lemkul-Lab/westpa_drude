import mdtraj
import numpy

parent = mdtraj.load('abl1_i1_eq_drude.pdb', top='abl1_i1_eq_drude.pdb')
traj = mdtraj.load('seg.dcd', top='abl1_i1_eq_drude.pdb')

dist1_parent = mdtraj.compute_distances(parent, [[1965,4239]], periodic=True)
dist1_traj = mdtraj.compute_distances(traj, [[1965,4239]], periodic=True)

dist2_parent = mdtraj.compute_distances(parent, [[4260,681]], periodic=True)
dist2_traj = mdtraj.compute_distances(traj, [[4260,681]], periodic=True)

dih1_parent = mdtraj.compute_dihedrals(parent, [[4187, 4203, 4220, 4233]], periodic=True)
dih1_traj = mdtraj.compute_dihedrals(traj, [[4187, 4203, 4220, 4233]], periodic=True)

dih2_parent = mdtraj.compute_dihedrals(parent, [[4187, 4203, 4247, 4260]], periodic=True)
dih2_traj = mdtraj.compute_dihedrals(traj, [[4187, 4203, 4247, 4260]], periodic=True)

angle1_parent = mdtraj.compute_angles(parent, [[4150, 4233, 4260]], periodic=True)
angle1_traj = mdtraj.compute_angles(traj, [[4150, 4233, 4260]], periodic=True)


###

###

dih1 = numpy.append(dih1_parent,dih1_traj)
dih1_arr = numpy.asarray(dih1)
dih1_arr = numpy.degrees(dih1_arr)
dih1_arr = dih1_arr+180
dih1_arr = dih1_arr

dih2 = numpy.append(dih2_parent,dih2_traj)
dih2_arr = numpy.asarray(dih2)
dih2_arr = numpy.degrees(dih2_arr)
dih2_arr = dih2_arr+180
dih2_arr = dih2_arr

dist2 = numpy.append(dist2_parent,dist2_traj)
d2_arr = numpy.asarray(dist2)
d2_arr = d2_arr*10

dist1 = numpy.append(dist1_parent,dist1_traj)
d1_arr = numpy.asarray(dist1)
d1_arr = d1_arr*10


angle1 = numpy.append(angle1_parent,angle1_traj)
a1_arr = numpy.asarray(angle1)
a1_arr = numpy.degrees(a1_arr)

numpy.savetxt("dist.dat", d1_arr)
numpy.savetxt("dist2.dat", d2_arr)
numpy.savetxt("dih1.dat", dih1_arr)
numpy.savetxt("dih2.dat", dih2_arr)
numpy.savetxt("angle1.dat", a1_arr)
