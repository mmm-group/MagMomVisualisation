from pymatgen.io.vasp.outputs import Outcar
from pymatgen.io.vasp import Poscar
from pymatgen.core import Structure
from sklearn.preprocessing import normalize
import numpy as np 
import re

#SETTINGS:
scalelength=0.4 #scale length of vectors to better visualise
minmom=0.1 #minimum magnetic moment to show
vecradius=0.5 #vector radius
veccolour=[255,0,0] #vector RGB colour

#Read in lattice vectors and work out matrix to transform magnetization vectors
struct = Structure.from_file(filename="./CONTCAR")
matrix = struct.lattice.matrix.copy()
normed_matrix = normalize(matrix, axis=1, norm='l2')
T=np.linalg.inv(normed_matrix)

#Read in x,y,z compnents of total magnetization and transform to a,b,c basis
vaspout=Outcar('OUTCAR')
magvecs=[]
for entry in vaspout.magnetization:
  mag=[entry['tot'][0],entry['tot'][1],entry['tot'][2]]
  magvec = np.array(mag)
  transmagvec= np.matmul(magvec,T)
  magvecs.append(transmagvec)

#Prepare strings to insert into VESTA file
VECTR_str=r"\1 "
VECTT_str=r"\1 "
VECTS_str="VECTS {0}".format(scalelength)
for i,atommag in enumerate(magvecs):
  if np.linalg.norm(atommag) > minmom: # only create vector if more than cutoff modulus
    VECTR_str += "{0} {1} {2} {3} 0\n {0} 0 0 0 0\n 0 0 0 0 0\n".format(i+1,atommag[0],atommag[1],atommag[2])
    VECTT_str += "{0} {1} {2} {3} {4} 1\n".format(i+1,vecradius,veccolour[0],veccolour[1],veccolour[2]) # set vector radius and colour   
#Read in VESTA file, insert strings and write file
data = open('struct.vesta','r').read()
newdata = re.sub(r'(VECTR\n)',VECTR_str,data)   
newdata = re.sub(r'(VECTT\n)',VECTT_str,newdata)
newdata=newdata.replace('VECTS 1.000000',VECTS_str)
file_out = open('struct-with-mag.vesta','w+')
file_out.write(newdata)
file_out.close()

