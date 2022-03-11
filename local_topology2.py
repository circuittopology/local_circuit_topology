import os
import numpy as np

f = open('./structure_path.txt',"r")
pdbids = f.readlines()
f.close()

p_num_out = []
ip_num_out = []
x_num_out = []
s_num_out = []
residues = []
coi = []

for n in range(len(pdbids)):
  print n
  pdb = pdbids[n].split()
  if len(pdb)>0:

    numrelations = np.genfromtxt('./out/len_cmap_'+pdb[0]+pdb[1]+'.txt',delimiter='\n')
    numresidues = np.genfromtxt('./out/nres_cmap_'+pdb[0]+pdb[1]+'.txt',delimiter='\n')
    prot_co = np.genfromtxt('./out/co_cmap_'+pdb[0]+pdb[1]+'.txt',delimiter='\n')
    numbering = np.genfromtxt('./out/n_'+pdb[0]+pdb[1]+'.txt',delimiter='\n')

    p_num = np.genfromtxt('./out/p_cmap_'+pdb[0]+pdb[1]+'.txt', delimiter=',')
    ip_num = np.genfromtxt('./out/ip_cmap_'+pdb[0]+pdb[1]+'.txt', delimiter=',')
    x_num = np.genfromtxt('./out/x_cmap_'+pdb[0]+pdb[1]+'.txt', delimiter=',')
    s_num = np.genfromtxt('./out/s_cmap_'+pdb[0]+pdb[1]+'.txt', delimiter=',')
    coi_all = np.genfromtxt('./out/coi_cmap_'+pdb[0]+pdb[1]+'.txt', delimiter=',')

    residues.append(numresidues)
    p_num_out.append(p_num[np.where(((int(pdb[2]))==numbering)>0)[0][0]])
    ip_num_out.append(ip_num[np.where(((int(pdb[2]))==numbering)>0)[0][0]])
    x_num_out.append(x_num[np.where(((int(pdb[2]))==numbering)>0)[0][0]])
    s_num_out.append(s_num[np.where(((int(pdb[2]))==numbering)>0)[0][0]])
    coi.append(coi_all[np.where(((int(pdb[2]))==numbering)>0)[0][0]])

  else:
    residues.append(' ')
    p_num_out.append(' ')
    ip_num_out.append(' ')
    x_num_out.append(' ')
    s_num_out.append(' ')
    coi.append(' ')

f = open('residues.txt',"a")
for i in range(len(residues)):
  f.write(str(residues[i])+'\n')
f.close()

f = open('p_number.txt',"a")
for i in range(len(p_num_out)):
  f.write(str(p_num_out[i])+'\n')
f.close()

f = open('ip_number.txt',"a")
for i in range(len(ip_num_out)):
  f.write(str(ip_num_out[i])+'\n')
f.close()

f = open('x_number.txt',"a")
for i in range(len(x_num_out)):
  f.write(str(x_num_out[i])+'\n')
f.close()

f = open('s_number.txt',"a")
for i in range(len(s_num_out)):
  f.write(str(s_num_out[i])+'\n')
f.close()
