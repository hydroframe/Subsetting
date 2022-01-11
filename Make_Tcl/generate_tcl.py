#!/usr/bin/env python3

import argparse
import os
import sys
import pfio


def read_infile(infile):

    # critical information
    crit_info = ['runname',
                 'Process.Topology.P', 'Process.Topology.Q',
                 'Process.Topology.R', 'file copy -force',
                 'ComputationalGrid.NX', 'ComputationalGrid.NY',
                 'ComputationalGrid.NZ', 'ComputationalGrid.DX',
                 'ComputationalGrid.DY', 'ComputationalGrid.DZ',
                 'GeomInput.domaininput.FileName', 'Geom.domain.Patches',
                 'dzScale.nzListNumber', 'Cell.0.dzScale.Value',
                 'Cell.1.dzScale.Value', 'Cell.2.dzScale.Value',
                 'Cell.3.dzScale.Value', 'Cell.4.dzScale.Value',
                 'Geom.domain.Perm.Value', 'TimingInfo.BaseUnit',
                 'TimingInfo.StartTime', 'TimingInfo.StopTime',
                 'TimingInfo.DumpInterval', 'TimeStep.Value',
                 'Geom.domain.Porosity.Value', 'BCPressure.PatchNames',
                 'Patch.top.BCPressure.Type', 'Patch.top.BCPressure.Cycle',
                 'Patch.top.BCPressure.rain.Value',
                 'Patch.top.BCPressure.rec.Value',
                 'Patch.top.BCPressure.alltime.Value',
                 'Solver.EvapTransFile',
                 'Solver.EvapTrans.FileName', 'TopoSlopesX.FileName',
                 'TopoSlopesY.FileName', 'pfdist',
                 'Geom.domain.ICPressure.Value',
                 'Geom.domain.ICPressure.RefPatch']

    with open(infile, 'r') as fi:
        content = fi.read()

    results = {}

    content = content.split('\n')

    # filter out commented, empty lines
    for ii, line in enumerate(content):
        if line:
            if line[0] != '#':
                for info in crit_info:
                    if info in line:
                        if info not in results.keys():
                            results[info] = {}
                            results[info]['locs'] = []
                            results[info]['vals'] = []
                        line = [x.strip() for x in line.split(' ') if x]
                        results[info]['locs'].append(ii)
                        results[info]['vals'].append(line)

    return results, content


parser = argparse.ArgumentParser(description='Generate .tcl script from '
                                             'parkinglot template')

# CORES
parser.add_argument('-P', type=int, default=2,
                    help='Processor in P(optional). Default is 2',
                    required=False)
parser.add_argument('-Q', type=int, default=1,
                    help='Processor in Q(optional). Default is 1',
                    required=False)
parser.add_argument('-R', type=int, default=1,
                    help='Processor in R(optional). Default is 1',
                    required=False)

# OUT FILE
parser.add_argument('-o', '--out_file', type=str,
                    help='output file', required=True)

# INPUT FILES
parser.add_argument('-i', '--temp_file', type=str,
                    help='template file for modifying', required=True)
parser.add_argument('--runname', type=str,
                    help='run name of the simulation', required=True)
parser.add_argument('-sl', '--slope_file', type=str,
                    help='name of the slope file, either in x or y direction',
                    required=True)
parser.add_argument('-so', '--solid_file', type=str,
                    help='name of the solid file', required=True)
parser.add_argument('-evap', choices=[0, 1], default=0,
                    help='PME file for simulation (optional).Default is 0')
parser.add_argument('--evap_file', type=str, help='name of the PME file')

# PARAMETERS
parser.add_argument('-K', '--perm', type=float, default=0.02849,
                    help='Permeability (optional). Default is 0.02849 (i.e. '
                         'average K of the Upper Colorado basin)',
                         required=False)
parser.add_argument('--porosity', type=float, default=0.39738,
                    help='Porosity (optional). Default is 0.39738 (i.e. '
                         'average porosity of the Upper Colorado basin)',
                         required=False)
parser.add_argument('--rain', type=float, default=-0.05,
                    help='Rain value (optional). Default is -0.05',
                    required=False)
parser.add_argument('--rec', type=float, default=0.0,
                    help='Recession value (optional). Default is 0.0',
                    required=False)
parser.add_argument('--flow', type=str, default='OverlandFlow',
                    help='Flow type for top layer (optional). Default is '
                    'OverlandFlow', required=False)
parser.add_argument('--constant', type=int, default=0,
                    help='Constant boundary condition for top layer '
                         '(optional). Default is No (0)', required=False)
parser.add_argument('--initw', choices=['top', 'bottom'], default='bottom',
                    help='Where to set initial pressure (optional). '
                         'Default is bottom', required=False)
parser.add_argument('--init', type=float, default=0.0,
                    help='Initial pressure for bottom layer (optional). '
                         'Default is 0.0', required=False)

# TIMING INFO
parser.add_argument('-s', '--start', type=float, default=0.0,
                    help='Start time(optional). Default is 0.0',
                    required=False)
parser.add_argument('-e', '--end', type=float,
                    help='End time (required)', required=True)
parser.add_argument('--baseu', type=float, default=0.1,
                    help='Base unit(optional). Default is 0.1',
                    required=False)
parser.add_argument('-ts', '--timestep', type=float, default=0.1,
                    help='Constant time step (optional). Default is 0.1',
                    required=False)
parser.add_argument('--dump', type=float, default=-1.0,
                    help='Dump interval (optional). Default is -1',
                    required=False)

# DOMAIN RESOLUTIONS
parser.add_argument('-dx', type=float, default=1000.,
                    help='dx resolution (optional). Default is 1000.')
parser.add_argument('-dy', type=float, default=1000.,
                    help='dy resolution (optional). Default is 1000.')
parser.add_argument('-dz', type=float, default=1000.,
                    help='dz resolution (optional). Default is 1000.')
parser.add_argument('-nz', type=int, default=1,
                    help='Number of z layers (optional). Default is 1.')
parser.add_argument('--dz_scales', nargs='+',
                    help='dz scale to be multiplied with dz (optional). '
                         'Default is 0.5', required=False)
parser.add_argument('--batches', nargs='+',
                    help='batches in domain (required)', required=True)


# parsing arguments
args = parser.parse_args()

out_file = args.out_file

# parsing processor info
P = args.P
Q = args.Q
R = args.R

# parsing input files
temp_file = args.temp_file
runname = args.runname
slope_file = args.slope_file
solid_file = args.solid_file
evap_choice = args.evap
evap_file = args.evap_file

if evap_choice == 0 and evap_file:
    parser.error('--evap_file can only be set when -evap=1.')
elif evap_choice == 1 and evap_file is None:
    parser.error('--evap_file must be set when -evap=1')

# parsing parameters
K = args.perm
poros = args.porosity
rain = args.rain
rec = args.rec
constant = args.constant
init = args.init
initw = args.initw
flow = args.flow

if constant == 1 and flow == 'OverlandFlow':
    parser.error('--flow should be set to FluxConst when --constant==1.')

# parsing timing info
start_time = args.start
end_time = args.end
baseu = args.baseu
timestep = args.timestep
dump = args.dump

# parsing domain info
dx = args.dx
dy = args.dy
dz = args.dz
nz = args.nz
batches = sorted([int(x) for x in args.batches])

if args.dz_scales is None:
    dz_scales = [0.5]
else:
    dz_scales = args.dz_scales

# read input file
results, content = read_infile(temp_file)

# change the content step by step
# change runname
results['runname']['vals'][0][-1] = runname

# change processors
results['Process.Topology.P']['vals'][0][-1] = str(P)
results['Process.Topology.Q']['vals'][0][-1] = str(Q)
results['Process.Topology.R']['vals'][0][-1] = str(R)

# change input files
if 'slopex' in slope_file:
    slope_file_y = slope_file.replace('slopex','slopey')
    slope_file_x = slope_file
elif 'slopey' in slope_file:
    slope_file_x = slope_file.replace('slopey','slopex')
    slope_file_y = slope_file
else:
    print('please check slopefile names')
    sys.exit()

results['file copy -force']['vals'][0][-2] = slope_file_x
results['file copy -force']['vals'][1][-2] = slope_file_y
results['file copy -force']['vals'][2][-2] = solid_file
if evap_choice == 0:
    results['file copy -force']['vals'][3][0] = '#'+results['file copy -force']['vals'][3][0]
else:
    results['file copy -force']['vals'][3][-2] = evap_file

# change parameters
results['Geom.domain.Perm.Value']['vals'][0][-1] = str(K)
results['Geom.domain.Porosity.Value']['vals'][0][-1] = str(poros)

# change boundary pressure of top layer
if evap_choice == 0:
    # comment out evap file
    results['Solver.EvapTransFile']['vals'][0][-1] = 'False'
    results['Solver.EvapTrans.FileName']['vals'][0][0] = '#'+results['Solver.EvapTrans.FileName']['vals'][0][0]
    if constant == 1:
        # comment out rain rec
        results['Patch.top.BCPressure.Cycle']['vals'][0][0] = '#'+results['Patch.top.BCPressure.Cycle']['vals'][0][0]
        results['Patch.top.BCPressure.rain.Value']['vals'][0][0] = '#'+results['Patch.top.BCPressure.rain.Value']['vals'][0][0]
        results['Patch.top.BCPressure.rec.Value']['vals'][0][0] = '#'+results['Patch.top.BCPressure.rec.Value']['vals'][0][0]
        # set constant input
        results['Patch.top.BCPressure.Cycle']['vals'][1][-1] = '\"constant\"'
        results['Patch.top.BCPressure.alltime.Value']['vals'][0][-1] = str(0.0)
        # set flow type for top layer
        results['Patch.top.BCPressure.Type']['vals'][0][-1] = flow
    else:
        # change rain and rec values
        results['Patch.top.BCPressure.rain.Value']['vals'][0][-1] = str(rain)
        results['Patch.top.BCPressure.rec.Value']['vals'][0][-1] = str(rec)
        # comment out constant input
        results['Patch.top.BCPressure.Cycle']['vals'][1][0] = '#'+results['Patch.top.BCPressure.Cycle']['vals'][1][0]
        results['Patch.top.BCPressure.alltime.Value']['vals'][0][0] = '#'+results['Patch.top.BCPressure.alltime.Value']['vals'][0][0]
else:
    # Change evap files
    results['Patch.top.BCPressure.Cycle']['vals'][1][-1] = '\"constant\"'
    results['Patch.top.BCPressure.alltime.Value']['vals'][0][-1] = str(0.0)
    results['Solver.EvapTransFile']['vals'][0][-1] = 'True'
    results['Solver.EvapTrans.FileName']['vals'][0][-1] = os.path.basename(evap_file)

# change initial pressure height
results['Geom.domain.ICPressure.Value']['vals'][0][-1] = str(init)
results['Geom.domain.ICPressure.RefPatch']['vals'][0][-1] = initw

# change timing info
results['TimingInfo.BaseUnit']['vals'][0][-1] = str(baseu)
results['TimingInfo.StartTime']['vals'][0][-1] = str(start_time)
results['TimingInfo.StopTime']['vals'][0][-1] = str(end_time)
results['TimingInfo.DumpInterval']['vals'][0][-1] = str(dump)
results['TimeStep.Value']['vals'][0][-1] = str(timestep)

# finally change domain info
results['GeomInput.domaininput.FileName']['vals'][0][-1] = os.path.basename(solid_file)
results['TopoSlopesX.FileName']['vals'][0][-1] = os.path.basename(slope_file_x)
results['TopoSlopesY.FileName']['vals'][0][-1] = os.path.basename(slope_file_y)

results['pfdist']['vals'][0][-1] = os.path.basename(slope_file_x)
results['pfdist']['vals'][1][-1] = os.path.basename(slope_file_y)

if evap_choice == 1:
    results['pfdist']['vals'][2][-1] = os.path.basename(evap_file)
else:
    results['pfdist']['vals'][2][0] = '#'+results['pfdist']['vals'][2][0]

# read slope file x into array
slope_x = pfio.pfread(slope_file_x)
nz0, ny0, nx0 = slope_x.shape

results['ComputationalGrid.NX']['vals'][0][-1] = str(nx0)
results['ComputationalGrid.NY']['vals'][0][-1] = str(ny0)
results['ComputationalGrid.NZ']['vals'][0][-1] = str(nz)

results['ComputationalGrid.DX']['vals'][0][-1] = str(dx)
results['ComputationalGrid.DY']['vals'][0][-1] = str(dy)
results['ComputationalGrid.DZ']['vals'][0][-1] = str(dz)

results['dzScale.nzListNumber']['vals'][0][-1] = str(nz)

for ni in range(nz):
    results['Cell.'+str(ni)+'.dzScale.Value']['vals'][0][-1] = str(dz_scales[ni])

# get batch name
batch_dict = {0: 'land', 1: 'ocean', 2: 'land', 3: 'top',
              4: 'lake', 5: 'sink', 6: 'bottom'}

patch_str = '\"'
for batch in batches:
    patch_str += batch_dict[batch]+' '

patch_str += '\"'
results['Geom.domain.Patches']['vals'][0] = results['Geom.domain.Patches']['vals'][0][:2]+[patch_str]

results['BCPressure.PatchNames']['vals'][0] = results['BCPressure.PatchNames']['vals'][0][:2]+[patch_str]

# write new file
for k in results.keys():
    locs = results[k]['locs']
    vals = results[k]['vals']
    for ix, loci in enumerate(locs):
        content[loci] = ' '.join(vals[ix])

with open(out_file, 'w') as fo:
    for line in content:
        fo.write(line+'\n')

