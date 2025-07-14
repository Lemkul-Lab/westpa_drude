#!/usr/bin/python
###########################################
# Runs a MD simulation with Drude FF.
#
# Outputs DCD files, last frame (CRD) and RST files with positions and
# velocities. Requires a PSF, a CRD and RST from equilibration. It allows
# restarting a simulation at any given point (using both .chk and .rst files).
# It also makes a backup of restarted .log file in order to save previous saved
# energies.
#
# USAGE: python Production_drude.py -h
#
# mdpoleto@vt.edu -> version 2023
###########################################
import openmm as mm
from openmm import *
from openmm.app import *
from openmm.unit import *
import MDAnalysis as mda
import warnings
import parmed as pmd
import time

from sys import stdout, exit, stderr
import os, math, fnmatch
import argparse

warnings.filterwarnings("ignore", message="Found no information for attr:")
warnings.filterwarnings("ignore", message="Found missing chainIDs")
warnings.filterwarnings("ignore", message="Supplied AtomGroup was missing the following attributes")
warnings.filterwarnings("ignore", message="DCDReader currently makes independent timesteps")


# Parse user input and options
ap = argparse.ArgumentParser(description=__doc__)

# Mandatory
ap.add_argument('-crd', type=str, default=None, required=True,
                help='Input coordinate file (.crd)')
ap.add_argument('-psf', type=str, default=None, required=True,
                help='Topology file in XPLOR format (.psf)')
ap.add_argument('-toppar', type=str, default='toppar.str', required=True,
                help='Force field stream file (ex. "toppar.str").')
ap.add_argument('-state', type=str, required=True,
                help='XML file to read positions/velocities from (.rst).')

# Options
ap.add_argument('-outname', type=str, default="output",
                help='Default name for output files. Default is "output".')

ap.add_argument('-runtime', default=100, type=float,
                help='Simulation length of each stride (in ps). Default 100.')
ap.add_argument('-nstride', default=5, type=int,
                help='Number of strides/chunks in which the simulation will be splitted. Default 5.')
ap.add_argument('-dt', default=1, type=int,
                help='Integration step (in fs). Default 1.')

ap.add_argument('-savefreq', type=float, default=5,
                help='Frequency (in ps) to save coordinates, checkpoints and trajectory.')
ap.add_argument('-printfreq', type=float, default=5,
                help='Frequency (in ps) to print and write in .log file.')

ap.add_argument('-firststride', type=int, default=1,
                help='First stride number. Default 1.')


ap.add_argument('-temp', default=298, type=float,
                help='Target temperature, in Kelvin. Default is 298.')
ap.add_argument('-pressure', default=1.0, type=float,
                help='Target pressure, in bar. Default is 1.0.')

cmd = ap.parse_args()

with open('Production_drude.dat', 'w') as out:
	out.write(' '.join(sys.argv[1:]))

#############################################
# Although we use checkpoints to restart simulations, an unexpected crash may
# harm the dcd integrity beyond repair. It is rare, but that may happenself.
# Therefore, we use 10 simulation strides over a loop, creating 10 dcd files that
# are concatenated afterwards.

jobname = cmd.outname

stride_time = cmd.runtime*picosecond		# in ps
dt = cmd.dt*femtosecond						#fs
nstride = cmd.nstride

print_freq  = cmd.printfreq*picosecond
savcrd_freq = cmd.savefreq*picosecond

temperature = cmd.temp*kelvin
pressure	= cmd.pressure*bar

nsteps = int(stride_time.value_in_unit(picosecond)/dt.value_in_unit(picosecond))
nprint  = int(print_freq.value_in_unit(picosecond)/dt.value_in_unit(picosecond))
nsavcrd = int(savcrd_freq.value_in_unit(picosecond)/dt.value_in_unit(picosecond))

#############################################
# Defining functions to use below:
def backup_old_log(pattern, string):
	result = []
	for root, dirs, files in os.walk("./"):
		for name in files:
			if fnmatch.fnmatch(name, pattern):

				try:
					number = int(name[-2])
					avail = isinstance(number, int)
					#print(name,avail)
					if avail == True:
						result.append(number)
				except:
					pass

	if len(result) > 0:
		maxnumber = max(result)
	else:
		maxnumber = 0

	backup_file = "\#" + string + "." + str(maxnumber + 1) + "#"
	os.system("mv " + string + " " + backup_file)
	return backup_file

def get_cubic_box(psf, rstfile):

	f = open(rstfile, 'r')

	box = {}

	while True:
		line = f.readline()
		if not line: break

		if line.split()[0] == "<A":
			size = line.split()[1].strip('x="')
			box['A'] = float(size)
		elif line.split()[0] == "<B":
			size = line.split()[2].strip('y="')
			box['B'] = float(size)
		elif line.split()[0] == "<C":
			size = line.split()[3].strip('z="').strip('"/>')
			box['C'] = float(size)
		else:
			pass

	boxX = box['A']*nanometer
	boxY = box['B']*nanometer
	boxZ = box['C']*nanometer

	psf.setBox(boxX, boxY, boxZ)

	return psf

def read_toppar(filename):
	extlist = ['rtf', 'prm', 'str']

	parFiles = ()
	for line in open(filename, 'r'):
		if '!' in line: line = line.split('!')[0]
		parfile = line.strip()
		if len(parfile) != 0:
			ext = parfile.lower().split('.')[-1]
			if not ext in extlist: continue
			parFiles += ( parfile, )

	params = CharmmParameterSet( *parFiles )
	return params, parFiles
##############################################

#############################################
print("\n> Simulation details:\n")
print("\tJob name = " + jobname)
print("\tCRD file = " + str(cmd.crd))
print("\tPSF file = " + str(cmd.psf))
print("\tToppar stream file = " + str(read_toppar(cmd.toppar)[1]))

print("\n\tSimulation_time = " + str(stride_time*nstride))
print("\tIntegration timestep = " + str(dt))
print("\tTotal number of steps = " +  str(nsteps*nstride))
print("\tNumber of strides = " + str(cmd.nstride) + " (" + str(stride_time) + " in each stride)")

print("\n\tSave coordinates each " + str(savcrd_freq))
print("\tSave checkpoint each " + str(savcrd_freq))
print("\tPrint in log file each " + str(print_freq))

print("\n\tTemperature = " + str(temperature))
print("\tPressure = " + str(pressure))
#############################################

print("\n> Setting the system:\n")
print("\t- Reading force field directory...")
charmm_params = read_toppar(cmd.toppar)[0]

print("\t- Reading topology and structure file...")
psf = CharmmPsfFile(cmd.psf)
crd = PDBFile(cmd.crd) #CharmmCrdFile(cmd.crd)

print("\t- Setting box (using information on -state file)...")
psf = get_cubic_box(psf, cmd.state)

print("\t- Creating system and setting parameters...")
system = psf.createSystem(charmm_params, nonbondedMethod=PME, nonbondedCutoff=1.2*nanometer, switchDistance=1.0*nanometer, ewaldErrorTolerance = 0.0001, constraints=HBonds)

nbforce = [system.getForce(i) for i in range(system.getNumForces()) if isinstance(system.getForce(i), NonbondedForce)][0]
nbforce.setNonbondedMethod(NonbondedForce.PME)
nbforce.setEwaldErrorTolerance(0.0001)
nbforce.setCutoffDistance(1.2*nanometer)
nbforce.setUseSwitchingFunction(True)
nbforce.setSwitchingDistance(1.0*nanometer)

# not every system has NBFIX terms, so check
cstnb = [system.getForce(i) for i in range(system.getNumForces()) if isinstance(system.getForce(i), CustomNonbondedForce)]
if cstnb:
	nbfix = cstnb[0]
	nbfix.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
	nbfix.setCutoffDistance(1.2*nanometer)
	nbfix.setUseSwitchingFunction(True)
	nbfix.setSwitchingDistance(1.0*nanometer)

print("\t- Setting barostat...")
system.addForce(MonteCarloBarostat(pressure, temperature))

print("\t- Setting integrator and thermostat...")
integrator = DrudeLangevinIntegrator(temperature, 5/picosecond, 1*kelvin, 20/picosecond, 0.001*picoseconds)
integrator.setMaxDrudeDistance(0.02) # Drude Hardwall

##
integrator.setRandomNumberSeed(RAND)
platform = Platform.getPlatformByName('CUDA')
print(os.environ['WM_PROCESS_INDEX'])
process_id = os.environ['WM_PROCESS_INDEX']
properties = {'Precision': 'mixed'}
properties["DeviceIndex"] = process_id
##

print('\t- Setting simulation context...')
simulation = Simulation(psf.topology, system, integrator, platform, properties)
simulation.context.setPositions(crd.positions)
simulation.context.computeVirtualSites()

print('\t- Using platform:', simulation.context.getPlatform().getName())

# Opening a loop of extension NSTRIDE to simulate the entire STRIDE_TIME*NSTRIDE
print("\n\n>>> Simulating... <<<")

dcd_file = jobname + ".dcd"
chk_file = jobname + ".chk"
log_file = jobname + ".log"
rst_file = jobname + ".rst"
pdb_out_file = jobname + ".pdb"
crd_out_file = jobname + ".crd"


print("> Loading previous state from > " + cmd.state + " <")
with open(cmd.state, 'r') as f:
	simulation.context.setState(XmlSerializer.deserialize(f.read()))
	currstep = 0
	currtime = currstep*dt.in_units_of(picosecond)
	simulation.currentStep = currstep
	simulation.context.setTime(currtime)
	print("> Current time: " + str(currtime) + " (Step = " + str(currstep) + ")")


dcd = DCDReporter(dcd_file, nsavcrd)
firstdcdstep = (currstep) + nsavcrd
dcd._dcd = DCDFile(dcd._out, simulation.topology, simulation.integrator.getStepSize(), firstdcdstep, nsavcrd) # charmm doesn't like first step to be 0

simulation.reporters.append(dcd)
simulation.reporters.append(StateDataReporter(stdout, nprint, step=True, speed=True, progress=True, totalSteps=nsteps, remainingTime=True, separator='\t\t'))
simulation.reporters.append(StateDataReporter(log_file, nprint, step=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True))

print("\n> Simulating " + str(nsteps) + " steps... ")
simulation.step(nsteps)

simulation.reporters.clear() # remove all reporters so the next iteration don't trigger them.


##################################
# Writing last frame information of stride
print("\n> Writing stride state file (" + str(rst_file) + ")...")
state = simulation.context.getState( getPositions=True, getVelocities=True )
with open(rst_file, 'w') as f:
	f.write(XmlSerializer.serialize(state))


print("> Writing last coordinate (" + str(pdb_out_file) + ")...")
u = mda.Universe(cmd.psf, dcd_file)
mda_system = u.select_atoms('all')
for ts in u.trajectory[-1]:
	with mda.Writer(str(pdb_out_file), system.n_atoms) as W:
		W.write(mda_system)

	#with mda.Writer(str(crd_out_file), system.n_atoms, extended=True) as W2: #uncomment to write a .crd file too
	#	W2.write(mda_system)
	
print("\n> Finished!\n")
