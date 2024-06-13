
import sys, os
import numpy as np
import re, random, argparse
######## INPUT PARSER ########
parser = argparse.ArgumentParser(description="Process some integers.")

# Add the arguments
parser.add_argument('--n_chain', type=int, default=500, help='Number of chains')
parser.add_argument('--n_bead', type=int, default=500, help='Number of beads')
parser.add_argument('--n_crosslink', type=int, default=8000, help='Number of quartic')
parser.add_argument('--n_id', type=int, default=1, help='ID of assignment')
parser.add_argument('--n_ratio', type=int, default=2, help='Ratio of ligands to metal ion for reversible cross-link')
parser.add_argument('--spacing', type=str, default='equal', help='Type of assignment')
parser.add_argument('--verbose', action='store_true', default = False, help='Verbose')
parser.add_argument('--infile', type=str, default='polymer_melt.data', help='Input file')
parser.add_argument('--outfile', type=str, default='linear_equal.data', help='Output file')

# Parse the arguments
args = parser.parse_args()

# Assign arguments to variables
n_chain = args.n_chain
n_bead = args.n_bead
n_crosslink = args.n_crosslink
spacing = args.spacing
n_id = args.n_id
n_ratio = args.n_ratio
infile = args.infile
outfile = args.outfile

# n_id is a bookkeeping variable
if n_id == 1:
	add_atom_id, at_start, pos_start, rm_bonds = 0, 2, 11, 0
else:
	add_atom_id, add_line, rm_bonds = 2, 4, 2

if args.verbose:
# Now you can use the variables as needed
    print(f'number of chains in input file: {n_chain}')
    print(f'number of beads per chain: {n_bead}')
    print(f'number of reversible bonds to assign: {n_crosslink}')
    print(f'Bonds will be placed with this spacing: {spacing}')
    print(f'Number of ligands to metal ion ratio: {n_ratio}')
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
    
DATA_ATOMs = np.zeros((int(Atoms),9)) if n_id ==1 else np.zeros((int(Atoms),10))
DATA_ATOM = np.zeros((int(Atoms),10))
DATA_BONDS = np.zeros((int(Bonds),4))
DATA_VELOCITY = np.zeros((int(Atoms),4))
for i in range(int(Atoms)):
    li=fid.readline()
    a=re.findall("[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", li.strip())
    DATA_ATOMs[i,:]=np.array(a,dtype=float)

if n_id == 1:
	DATA_ATOM[:,:3]=DATA_ATOMs[:,:3].astype(int)
	DATA_ATOM[:,3]=0.0
	DATA_ATOM[:,4:7]=DATA_ATOMs[:,3:6]
	DATA_ATOM[:,7:]=DATA_ATOMs[:,6:].astype(int)
else:
	DATA_ATOM = DATA_ATOMs.copy()
DATA_ATOM=DATA_ATOM[np.argsort(DATA_ATOM[:,0])]
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

nonbond_chain = int(n_crosslink/n_bead); # debond 16 chains
nonbond_atom = n_crosslink
if n_id != 1:
	bond_ref =int(Bonds)-(nonbond_atom-nonbond_chain) - rm_bonds
	aa = np.where(DATA_BONDS[:,2]==1)[0]	
	bb = np.where(DATA_BONDS[:,3]==250001)[0]	
	cc = np.where(DATA_BONDS[:,2]==n_bead+1)[0]
	DATA_BONDS = np.delete(DATA_BONDS,[np.intersect1d(aa,bb),np.intersect1d(bb,cc)],0)	
DATA_BONDS=DATA_BONDS[np.argsort(DATA_BONDS[:,2])]
for i in range(int(Bonds)-rm_bonds):
    DATA_BONDS[i,0]=i+1

DATA_BOND=DATA_BONDS[:int(Bonds)-(nonbond_atom-nonbond_chain)-rm_bonds,:]
if n_id!=1:
	DATA_BOND = np.append(DATA_BOND,[[bond_ref+1,2,1,250001],[bond_ref+2,2,n_bead+1,250001]],axis=0)
if n_crosslink>0:    
    if n_id ==1:	
        DATA_ATOM[int(Atoms)-nonbond_atom:,2]=2+add_atom_id
    else:
        DATA_ATOM[int(Atoms)-nonbond_atom-1:-1,2]=2+add_atom_id

    new_atom_per_chain = int(np.floor(n_ratio*nonbond_atom/(n_chain-nonbond_chain)));
    if new_atom_per_chain == 0:
        new_atom_per_chain = 1
    space_new_atom = int(np.floor(n_bead/new_atom_per_chain));
    
    tail = int(np.ceil(space_new_atom/2));
    if spacing == 'equal':
        for i in range(n_chain-nonbond_chain):
    	    for j in range(new_atom_per_chain):
    	        DATA_ATOM[n_bead*i+space_new_atom*(j+1)-tail, 2] = 3+add_atom_id;
    elif spacing == 'random':
        for i in range(n_chain-nonbond_chain):
    	    dumct=random.sample(range(0,n_bead-1),new_atom_per_chain)
    	    for j in range(new_atom_per_chain):
    	        DATA_ATOM[n_bead*i+dumct[j], 2] = 3+add_atom_id;
        
fid.close()
fid=open(outfile,'w')
fid.write('LAMMPS data file by inserting crosslink bond using Python Script\n\n');
fid.write('%d atoms\n'%Atoms);
fid.write('%d bonds\n\n'%(Bonds-(nonbond_atom-nonbond_chain)));
atom_types = 3 + add_atom_id;
fid.write('%d atom types\n'%atom_types);
fid.write('%d bond types\n\n'%(1+add_atom_id//2));
fid.write('%.16e %.16e xlo xhi\n'%(box_x[0],box_x[1]));
fid.write('%.16e %.16e ylo yhi\n'%(box_y[0],box_y[1]));
fid.write('%.16e %.16e zlo zhi\n\n'%(box_z[0],box_z[1]));
if n_id == 1:
	fid.write('Masses\n\n1 1\n2 1\n3 1\n\n');
else:
	fid.write('Masses\n\n1 1\n2 1\n3 1\n4 1\n5 1\n\n');
fid.write('Atoms\n\n');
for i in range(int(Atoms)):
    fid.write('%d %d %d %.16e %.16e %.16e %.16e %d %d %d \n'%(DATA_ATOM[i,0],DATA_ATOM[i,1],DATA_ATOM[i,2],DATA_ATOM[i,3],DATA_ATOM[i,4],DATA_ATOM[i,5],DATA_ATOM[i,6],DATA_ATOM[i,7],DATA_ATOM[i,8],DATA_ATOM[i,9]))
fid.write('\nVelocities \n\n');
for i in range(int(Atoms)):
    fid.write('%d %.16e %.16e %.16e \n'%(DATA_VELOCITY[i,0],DATA_VELOCITY[i,1],DATA_VELOCITY[i,2],DATA_VELOCITY[i,3]))
fid.write('\nBonds \n\n');
for i in range(int(Bonds-(nonbond_atom-nonbond_chain))):
    fid.write('%d %d %d %d \n'%(DATA_BOND[i,0],DATA_BOND[i,1],DATA_BOND[i,2],DATA_BOND[i,3]))
fid.close()
