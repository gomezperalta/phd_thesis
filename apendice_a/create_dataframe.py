# -*- coding: utf-8 -*-
"""
@author: iG
"""

import pandas as pd
import numpy as np
import pymatgen
import Wyckoff_finder as wf
import time

directorio = input('Provide the directory where the CIF-files are:'+ \
              '\n')
cif_list = input('Provide the txt-file where the CIF-file list is (no-extension)' +\
                 '\n')

df = pd.read_csv(cif_list + '.txt', header=None)
df = df.rename(columns={0:'cif'})

formulas = list()
sgnums = list()
element_list = list()
wyck_list = list()
site_list = list()
atom_list = list()
not_added_cifs = list()

diccio = np.load('WyckoffSG_dict.npy').item()['wyckmul']

samples = df.shape[0]

start = time.time()
for row in range(df.shape[0]):
    
    try:
        estructura = pymatgen.Structure.from_file(directorio + \
                                                  str(df['cif'][row]) + \
                                                  '.cif')
        
        sgnum = pymatgen.symmetry.analyzer.SpacegroupAnalyzer(estructura).get_space_group_number()
        formula = estructura.composition.reduced_formula
        elements = len(estructura.composition.elements)
        wyckdict = wf.wyckoff_occupation(ruta = directorio,
                                     archivo = str(df['cif'][row]))
        sites = len(wyckdict)
        
        sg_diccio = diccio[str(sgnum).zfill(3)]  
        
        atoms = 0      
        for item in range(len(wyckdict)):
            
            label = wyckdict[item].keys()
            label = list(label)[0]
            fracsum = np.sum(list(wyckdict[item][label].values()))
            multiplicity = int(sg_diccio[label[0]])
            atoms += multiplicity*fracsum
       
    except:
        formula = None
        sgnum = None
        elements = None
        wyckdict = None
        sites = None
        atoms = None
        not_added_cifs += [df['cif'][row]]
   
    formulas += [formula]
    sgnums += [sgnum]
    element_list += [elements]
    wyck_list += [wyckdict]    
    site_list += [sites]
    atom_list += [atoms]
    
df['formula'] = formulas
df['sgnum'] = sgnums
df['elements'] = element_list
df['WyckOcc'] = wyck_list
df['sites'] = site_list
df['atoms'] = atom_list

df = df[df['WyckOcc'] != None].reset_index(drop=True)

if len(not_added_cifs) != 0:
    with open('not_added_cif.txt','w') as f:
        for item in not_added_cifs:
            f.write(str(item)+'\n')
        f.close()
        
df.to_csv('cod_dataframe.csv',index=None)
df.to_pickle('cod_dataframe.pkl')

print('Process lasted',np.round(time.time()-start,2),
      's to treat',samples,'samples') 
