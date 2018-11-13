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

# dumps
dump_period = 1e4
flw_steps = 1e7

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

# force params
id = get_array_id();
fx_arr = [0.001,0.003,0.005,0.007,0.009]
fx_fnl = fx_arr[id]
fName = '-eql-' + str(fx_fnl)

# read in the Init file
system = deprecated.init.read_xml(filename='Init' + fName + '.xml')


# pair forces
nl = md.nlist.cell()
dpd = md.pair.dpd(r_cut=1.0, nlist=nl, kT=0.1, seed=id)

dpd.pair_coeff.set('wall', 'S1', A=3.0, gamma=4.5)
dpd.pair_coeff.set('wall', 'wall', A=0.0, gamma=4.5)
dpd.pair_coeff.set('S2','S1',A=60.0,gamma=4.5)
dpd.pair_coeff.set('S1','S1',A=25.0,gamma=4.5)
dpd.pair_coeff.set('S2','S2',A=25.0,gamma=4.5)
dpd.pair_coeff.set('S2','wall',A=10.0,gamma=4.5)

md.integrate.mode_standard(dt=0.01)

# set up groups
all = hoomd.group.all()
groupWALL = hoomd.group.type(name='groupWALL', type='wall')
notWALL=hoomd.group.difference(name="particles-not-typeWALL", a=all, b=groupWALL)
md.integrate.nve(group=notWALL)

#make sure that the neighbor lists do not include particles that are bonded or are in the same body do not interact with one another
nl.reset_exclusions(exclusions = ['bond', 'body'])

#start the logs, including restart file
hoomd.dump.dcd('traj-flow' + fName + '.dcd', period = dump_period*10, overwrite = True)
#deprecated.dump.xml(group = all, filename='Rstrt'+fName,  period = dump_period*10)

# dump the system data - position, velocity
# NOTE: this is compressed system data that
# will be post processed later with gtar
zipfile = hoomd.dump.getar('dump' + fName + '.zip',
                static=['type'],
                dynamic={'position': dump_period, 'velocity': dump_period})

# logs
hoomd.analyze.log('energies-flow' + fName + '.txt', quantities=['temperature', 'potential_energy', 'kinetic_energy'], header_prefix='#', period = dump_period, overwrite = True)

# start up the pos writer
pos = deprecated.dump.pos(filename="dump-flow" + fName + ".pos", period=dump_period)
pos.set_def('S1', 'sphere 1 CC0000')
pos.set_def('wall', 'sphere 1 336600')
pos.set_def('S2', 'sphere 1 0000FF')

# main run
# apply the constant force
const = md.force.constant(fx=fx_fnl, fy=0.0, fz=0.0)
hoomd.run(flw_steps)

# dump the final xml
deprecated.dump.xml(filename='Final-flow' + str(fx_fnl) + '.xml', all=True, group=all)