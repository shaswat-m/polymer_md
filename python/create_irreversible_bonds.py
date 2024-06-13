
import sys, os
import numpy as np
import re, random
from util import verletlist, find_quartic_bonds, uniform_bond_indices, excluded_volume_indices, uniform_spacing_indices
import argparse

######## INPUT PARSER ########
parser = argparse.ArgumentParser(description="Process some integers.")

# Add the arguments
parser.add_argument('--n_chain', type=int, default=500, help='Number of chains')
parser.add_argument('--n_bead', type=int, default=500, help='Number of beads')
parser.add_argument('--n_quartic', type=int, default=13600, help='Number of quartic')
parser.add_argument('--assignment_type', type=str, default='random', help='Type of assignment')
parser.add_argument('--system_type', type=str, default='bonded', help='Type of system')
parser.add_argument('--lfrac', type=float, default=0.8, help='Fraction')
parser.add_argument('--verbose', action='store_true', default = False, help='Verbose')
parser.add_argument('--infile', type=str, default='polymer_melt.data', help='Input file')
parser.add_argument('--outfile', type=str, default='polymer_irreversible.data', help='Output file')

# Parse the arguments
args = parser.parse_args()

# Assign arguments to variables
n_chain = args.n_chain
n_bead = args.n_bead
n_quartic = args.n_quartic
assignment_type = args.assignment_type
system_type = args.system_type
lfrac = args.lfrac
infile = args.infile
outfile = args.outfile

# Set other variables based on the system_type
if system_type == 'bonded':
    at_start, pos_start, ncs, frc = 2, 16, 9, 1.25
    angular = False
else:
    at_start, pos_start, ncs, frc = 1, 8, 6, 2000
    angular = True

if args.verbose:
# Now you can use the variables as needed
    print(f'number of chains in input file: {n_chain}')
    print(f'number of beads per chain: {n_bead}')
    print(f'number of irreversible bonds to assign: {n_quartic}')
    print(f'Bonds will be chosen from candidate sites using: {assignment_type}')
    if assignment_type == 'hetero':
        print(f'fraction of crosslinks in the left half: {lfrac}')
    print(f'Polymer melt data file that is being read: {infile}')
    print(f'Cross-linked network will be written into: {outfile}')
##################################


fid=open(infile,'r')
for i in range(at_start):
    li=fid.readline()
li=fid.readline()
a=re.findall("[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", li.strip())
Atoms=int(a[0])
li=fid.readline()
a=re.findall("[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", li.strip())
atom_types=int(a[0])
li=fid.readline()
a=re.findall("[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", li.strip())
Bonds=int(a[0])
li=fid.readline()
a=re.findall("[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", li.strip())
bond_types=int(a[0])
if angular:
    li=fid.readline()
    a=re.findall("[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", li.strip())
    Angles=int(a[0])
    li=fid.readline()
    a=re.findall("[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", li.strip())
    angle_types=int(a[0])

if not angular:
    li=fid.readline()
li=fid.readline()
a=re.findall("[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", li.strip())
box_x=np.array(a,dtype=float)
li=fid.readline()
a=re.findall("[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", li.strip())
box_y=np.array(a,dtype=float)
li=fid.readline()
a=re.findall("[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", li.strip())
box_z=np.array(a,dtype=float)
for i in range(pos_start):
    li=fid.readline()
    
DATA_ATOMs = np.zeros((int(Atoms),ncs))
DATA_ATOM = np.zeros((int(Atoms),ncs))
DATA_BONDS = np.zeros((int(Bonds)+n_quartic,4))
DATA_VELOCITY = np.zeros((int(Atoms),4))
for i in range(int(Atoms)):
    li=fid.readline()
    a=re.findall("[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", li.strip())
    DATA_ATOMs[i,:]=np.array(a,dtype=float)


DATA_ATOM = DATA_ATOMs.copy()
DATA_ATOM=DATA_ATOM[np.argsort(DATA_ATOM[:,0])]
if not angular:
    for i in range(3):
        li=fid.readline()
    for i in range(int(Atoms)):
        li=fid.readline()
        a=re.findall("[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", li.strip())
        DATA_VELOCITY[i,:]=np.array(a,dtype=float)
        
    DATA_VELOCITY[:,0]=DATA_VELOCITY[:,0].astype(int)
for i in range(3):
    li=fid.readline()
for i in range(int(Bonds)):
    li=fid.readline()
    a=re.findall("[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", li.strip())
    DATA_BONDS[i,:]=np.array(a,dtype=int)

if angular:
    DATA_ANGLES = np.zeros((int(Angles),5))
    for i in range(3):
        li=fid.readline()
    for i in range(int(Angles)):
        li=fid.readline()
        a=re.findall("[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", li.strip())
        DATA_ANGLES[i,:]=np.array(a,dtype=int)

pos = DATA_ATOM[:,3:6].copy()
h = np.diag([box_x[1]-box_x[0],box_y[1]-box_y[0],box_z[1]-box_z[0]])
if os.path.exists('bond_info.npz'):
    cand_bonds = np.load('bond_info.npz')['b']
    scaled_bond_pos = np.load('bond_info.npz')['sbp']
else:
    cand_bonds = find_quartic_bonds(pos,h,frc,n_bead)
    b_pos = 0.5*(pos[cand_bonds[:,0]-1,:]+pos[cand_bonds[:,1]-1,:])
    s_pos = np.dot(b_pos,np.linalg.inv(h))
    s_pos = s_pos - np.round(s_pos)
    scaled_bond_pos = s_pos.copy()
    np.savez('bond_info.npz',b = cand_bonds, sbp = scaled_bond_pos)

if assignment_type == 'random':
    n_cand = cand_bonds.shape[0]
    pick_range = np.linspace(0,n_cand-1,n_cand,dtype=int)
    pick_ind = np.random.choice(pick_range,n_quartic,replace=False)
elif assignment_type == 'kmeans':
    pick_ind = uniform_bond_indices(scaled_bond_pos, n_quartic)
elif assignment_type == 'exclusion':
    selected_indices = excluded_volume_indices(scaled_bond_pos, h, n_quartic)
    selected_indices = np.array(selected_indices)
    n_cand = len(selected_indices)
    pick_range = np.linspace(0, n_cand - 1, n_cand, dtype=int)
    choose_ind = np.random.choice(pick_range, n_quartic, replace=False)
    pick_ind = selected_indices[choose_ind.astype(int)]
elif assignment_type == 'spacing':
    selected_indices = uniform_spacing_indices(cand_bonds, n_quartic, n_bead)
    selected_indices = np.array(selected_indices)
    n_cand = len(selected_indices)
    pick_range = np.linspace(0, n_cand - 1, n_cand, dtype=int)
    choose_ind = np.random.choice(pick_range, n_quartic, replace=False)
    pick_ind = selected_indices[choose_ind.astype(int)]
elif assignment_type == 'hetero':
    l_quartic = int(n_quartic*lfrac)
    r_quartic = n_quartic - l_quartic
    l_cand = np.where(scaled_bond_pos[:,0]<0)[0]
    r_cand = np.where(scaled_bond_pos[:,0]>=0)[0]
    pick_ind = np.append(np.random.choice(l_cand,l_quartic,replace=False),np.random.choice(r_cand,r_quartic,replace=False))

for i in range(n_quartic):
    DATA_BONDS[int(Bonds)+i,:] = np.array([int(Bonds)+i+1,2,cand_bonds[pick_ind[i],0],cand_bonds[pick_ind[i],1]])
        
fid.close()
fid=open(outfile,'w')
fid.write('LAMMPS data file by inserting crosslink bond using Python Script\n\n');
fid.write('%d atoms\n'%Atoms);
fid.write('%d bonds\n\n'%(Bonds+n_quartic));
atom_types = 1;
fid.write('%d atom types\n'%atom_types);
if angular:
    fid.write('%d angles\n'%Angles);
    fid.write('%d angle types\n'%angle_types);
fid.write('%d bond types\n\n'%(2));
fid.write('%.16e %.16e xlo xhi\n'%(box_x[0],box_x[1]));
fid.write('%.16e %.16e ylo yhi\n'%(box_y[0],box_y[1]));
fid.write('%.16e %.16e zlo zhi\n\n'%(box_z[0],box_z[1]));
fid.write('Masses\n\n1 1.0\n\n');
fid.write('Atoms\n\n');
if not angular:
    for i in range(int(Atoms)):
        fid.write('%d %d %d %.16e %.16e %.16e %d %d %d \n'%(DATA_ATOM[i,0],DATA_ATOM[i,1],DATA_ATOM[i,2],DATA_ATOM[i,3],DATA_ATOM[i,4],DATA_ATOM[i,5],DATA_ATOM[i,6],DATA_ATOM[i,7],DATA_ATOM[i,8]))
    fid.write('\nVelocities \n\n');
    for i in range(int(Atoms)):
        fid.write('%d %.16e %.16e %.16e \n'%(DATA_VELOCITY[i,0],DATA_VELOCITY[i,1],DATA_VELOCITY[i,2],DATA_VELOCITY[i,3]))
else:
    for i in range(int(Atoms)):
        fid.write('%d %d %d %.16e %.16e %.16e \n'%(DATA_ATOM[i,0],DATA_ATOM[i,1],DATA_ATOM[i,2],DATA_ATOM[i,3],DATA_ATOM[i,4],DATA_ATOM[i,5]))
fid.write('\nBonds \n\n');
for i in range(int(Bonds)+n_quartic):
    fid.write('%d %d %d %d \n'%(DATA_BONDS[i,0],DATA_BONDS[i,1],DATA_BONDS[i,2],DATA_BONDS[i,3]))
if angular:
    fid.write('\nAngles \n\n');
    for i in range(int(Angles)):
        fid.write('%d %d %d %d %d \n'%(DATA_ANGLES[i,0],DATA_ANGLES[i,1],DATA_ANGLES[i,2],DATA_ANGLES[i,3],DATA_ANGLES[i,4]))

fid.close()
