# /usr/bin/env hoomd
import hoomd
from hoomd import md
from hoomd import deprecated
hoomd.context.initialize()

import os
import shutil
import re
import sys
import xml.etree.ElementTree as ET
import math as m
import numpy as np


######################
##### Init params #####
######################

# free parameters, units
# of particle diameter(1.0)
N_drop = 1 #number of droplets, only 0 or 1
r_drop = 6.0 #radius of droplet
lx = 60.0
ly = 40.0
lz = 40.0
mass = 1.0

# dumps
dump_period = 1e2
rsz_steps = 1e3
eql_steps = 1e3

# geometry
wall_dens = 61.35 #from Millan et al., JCP (2007)
solv_dens = 4.0

#Lx_fnl = Ly_fnl = Lz_fnl = 10.0
#v_box = Lx_fnl*Ly_fnl*Lz_fnl
v_box = lx*ly*lz
a = (4/wall_dens)**(1./3.) #fcc number density
N_solv = int(solv_dens*v_box) #solv number density

# get the job index from PBS_ARRAYID, or return 0 if it is not specified (for test jobs)
def get_array_id():
    pbs_arrayid = os.getenv('SLURM_ARRAY_TASK_ID');
    if pbs_arrayid is None:
        return 0
    else:
        return int(pbs_arrayid) - 1;

# set up the initial and final box dimensions
# beads are all diameter = 1.0
id = get_array_id();
fx_arr = [0.001,0.003,0.005,0.007,0.009]
fx_fnl = fx_arr[id]
fName = '-eql-' + str(fx_fnl)

# random solvent **placer**, no system init via hoomd
def random_solv_plcr():
    # from lx,ly,lz randomly pluck positions
    for i in range(int(N_solv)):
        x = lx*np.random.rand()-lx/2.
        y = ly*np.random.rand()-ly/2.
        z = lz*np.random.rand()-lz/2.
        grid_solv_pos.append([x,y,z])
    return

# an alternative wall builder
def fcc_wall_builder():
    # create an array of bead positions
    xs = np.linspace(-lx/2.,lx/2.,m.ceil(lx/a))
    ys = np.linspace(-ly/2.,ly/2.,m.ceil(ly/a))
    zs = np.array([a/2.,0,-a/2.])#don't want overlaps
    for z in zs:
       for y in ys:
            for x in xs:
                # place the solvent particles
                if z == 0:
                    x_new = x+a/2.
                    y_new = y+a/2.
                    if x_new > lx/2.:
                        x_new = x_new - lx
                    if y_new > ly/2.:
                        y_new = y_new - ly
                    grid_wall_pos.append([x_new,y_new,z])
                else:
                    grid_wall_pos.append([x,y,z])
    return

# declare pos array(s)
star_pos = []
grid_solv_pos = []
grid_wall_pos = []
masses = []
diameters = []
bodies = []
types = []
bonds = []

# call parser, which will fill pos and map_pos
#parse = xml_parser()

# create the particle grid
solv_build = random_solv_plcr() #solvent
wall_build = fcc_wall_builder()

# ***TESTING***
print('Length of solv_pos: ' + str(len(grid_solv_pos)))
print('Length of wall_pos: ' + str(len(grid_wall_pos)))

# aggreate all arrays
full_pos = grid_solv_pos + grid_wall_pos

# build a single list for each particle field
masses.extend(int(N_solv)*[str(mass)]) #tack on the solvent
masses.extend(len(grid_wall_pos)*[str(mass)]) #tack on the walls
diameters.extend(int(N_solv)*['1.0']) #tack on the solvent
diameters.extend(len(grid_wall_pos)*['1.0']) #tack on the walls
bodies.extend(int(N_solv)*['-1']) #tack on the solvent
bodies.extend(len(grid_wall_pos)*['-1']) #tack on the walls
types_solv = int(N_solv)*['S1']
types_wall = len(grid_wall_pos)*['wall']
types = types_solv + types_wall

#print('Print for testing')
print(len(masses))
print(len(diameters))
print(len(bodies))
print(len(types))
print(len(bonds))
print(len(full_pos))

# catch improper configurations
if len(full_pos) != len(masses):
    print('The position and mass arrays are not the same length!')
    sys.exit(-1)

# write out the file
with open('Init' + fName + '.xml', 'w') as inpFile:
    inpFile.write('\n'.join(['<?xml version="1.1" encoding="UTF-8"?>',
                             '<hoomd_xml version="1.5">',
                             '<configuration time_step = "0" dimensions = "3">']))
    inpFile.write('<box lx="{size1}" ly="{size2}" lz="{size3}" />\n'.format(size1=lx,size2=ly,size3=lz))
    inpFile.write('<position>\n' + '\n'.join('{} {} {}'.format(x, y, z) for (x, y, z) in full_pos)
        + '</position>\n')
    inpFile.write('<mass>\n' + '\n'.join(str(mass) for mass in masses) + '</mass>\n')
    inpFile.write('<diameter>\n' + '\n'.join(str(diameter) for diameter in diameters) + '</diameter>\n')
    inpFile.write('<body>\n' + '\n'.join(str(body) for body in bodies) + '</body>\n')
    inpFile.write('<type>\n' + '\n'.join(str(type) for type in types) + '</type>\n')
    inpFile.write('\n</configuration>\n</hoomd_xml>\n')



######################
###### Main Run ######
######################

# read in the Init file
system = deprecated.init.read_xml(filename='Init' + fName + '.xml', wrap_coordinates = True)

# NOTE: in these simulations epsilon has ben rescaled so that the ODT for the given volume fraction
# is closer to that of the cores.
nl = md.nlist.cell()
dpd = md.pair.dpd(r_cut=1.0, nlist=nl, kT=1.0, seed=0)

dpd.pair_coeff.set(system.particles.types, system.particles.types, A=20.0, gamma= 1.0)
dpd.pair_coeff.set('wall', system.particles.types, A=0.0, gamma= 1.0)

md.integrate.mode_standard(dt=0.01)

# set up groups
all = hoomd.group.all()
groupWALL = hoomd.group.type(name='groupWALL', type='wall')
notWALL=hoomd.group.difference(name="particles-not-typeWALL", a=all, b=groupWALL)
md.integrate.nve(group=notWALL)
#md.integrate.nve(group=group.all())

#make sure that the neighbor lists do not include particles that are bonded or are in the same body do not interact with one another
nl.reset_exclusions(exclusions = ['bond', 'body'])

# run to set up geometry
hoomd.run(rsz_steps)

######################
#### Drop Init ####
######################

# add the types via snapshot
snap = system.take_snapshot(bonds=True)
snap.bonds.types = ['S2']
system.restore_snapshot(snap)

# manipulate particle type(s)
if N_drop != 0:
    system.particles.types.add('S2')
    for p in system.particles:
        # create two half-spheres of particles, at box top & bottom
        # NOTE: origin is @ box center, as well as wall
        r_diff_pos = (p.position[0]**2 + p.position[1]**2 + \
                      (p.position[2]-8)**2)**(1./2.)
        if r_diff_pos <= r_drop:
            # do not allow wall types to change
            if p.type=='wall':
                print('Error: radius too large for size.')
                sys.exit(-1)
            p.type='S2'

# dump the final xml
deprecated.dump.xml(filename='Init' +  fName  + '.xml', all=True, group=all, image=False)