from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
#from dcdreporter import DCDReporter
import numpy as np


def OPLS_LJ(system):
    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    lorentz = CustomNonbondedForce('epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)*4.0')
    lorentz.setNonbondedMethod(nonbonded_force.CutoffPeriodic)
    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')
    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
    system.addForce(lorentz)
    LJset = {}
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        LJset[index] = (sigma, epsilon)
        lorentz.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(
            index, charge, sigma, epsilon * 0)
    for i in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
        # ALL THE 12,13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED
        # FORCE
        lorentz.addExclusion(p1, p2)
        if eps._value != 0.0:
            # print p1,p2,sig,eps
            sig14 = sqrt(LJset[p1][0] * LJset[p2][0])
            eps14 = sqrt(LJset[p1][1] * LJset[p2][1])
            nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
    return system
'''
def Minimize(simulation,iters=0):
    simulation.minimizeEnergy(maxIterations=iters)
    position = simulation.context.getState(getPositions=True).getPositions()
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    PDBFile.writeFile(simulation.topology, position,open('gasmin.pdb', 'w'))
    print('Energy at Minima is %3.3f kcal/mol' % (energy._value * KcalPerKJ))
    return simulation
'''
temperature = 640 * kelvin
pdb = PDBFile('in.pdb')
modeller = Modeller(pdb.topology, pdb.positions)
forcefield = ForceField('ptz.xml')
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=None)
#system = forcefield.createSystem(modeller.topology,constraints=None)
#system.addForce(AndersenThermostat(temperature*kelvin, 1/picosecond))
system.addForce(MonteCarloBarostat(1*bar, temperature*kelvin))
system = OPLS_LJ(system)
#integrator=VerletIntegrator(0.001*picoseconds)
integrator = LangevinIntegrator(temperature, 1 / picosecond,  0.001 * picoseconds)
platform = Platform.getPlatformByName('CUDA')
properties = {'Precision': 'double'}
simulation = Simulation(modeller.topology, system, integrator,platform,properties)
#simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
#print('Initial Potential Energy %3.3f kcal/mol' % (energy._value * KcalPerKJ))
#simulation.minimizeEnergy(tolerance=0.1*kilojoule/mole)
#energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
#print('Energy at Minima is %3.3f kcal/mol' % (energy._value * KcalPerKJ))
#position = simulation.context.getState(getPositions=True).getPositions()
#PDBFile.writeFile(simulation.topology, position,open('minised.pdb', 'w'))

simulation.reporters.append(DCDReporter('output.dcd', 10000))
simulation.reporters.append(StateDataReporter(stdout, 1000, time=True, volume=True,temperature=True,kineticEnergy=True, potentialEnergy=True,totalEnergy=True,density=True,speed=True))

simulation.step(5000000)
simulation.saveState(f"final1.xml")
simulation.step(5000000)
simulation.saveState(f"final2.xml")
#for i in range(1,21):
#    integrator.setTemperature((450+15*i)*kelvin)
#    simulation.step(1000000)
#position = simulation.context.getState(getPositions=True).getPositions()
#PDBFile.writeFile(simulation.topology, position,open("final.pdb", 'w'))
#for i in range(5):
#    simulation.step(1000000)
#    position = simulation.context.getState(getPositions=True).getPositions()
#    PDBFile.writeFile(simulation.topology, position,open("output"+str(i+1)+".pdb", 'w'))



