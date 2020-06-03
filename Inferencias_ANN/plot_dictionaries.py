#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 11:38:46 2019

@author: ig
"""

import numpy as np
import re
import matplotlib
import matplotlib.pyplot as plt



'''lattice parameter dictionaries are activated in next lines'''

oxides_latpar = np.load('pure_aristotype_joined_oxides_latpar.npy').item()
sulfides_latpar = np.load('pure_aristotype_joined_sulfides_latpar.npy').item()
selenides_latpar = np.load('pure_aristotype_joined_selenides_latpar.npy').item()
tellurides_latpar = np.load('pure_aristotype_joined_tellurides_latpar.npy').item()

fluorides_latpar = np.load('pure_aristotype_joined_fluorides_latpar.npy').item()
chlorides_latpar = np.load('pure_aristotype_joined_chlorides_latpar.npy').item()
bromides_latpar = np.load('pure_aristotype_joined_bromides_latpar.npy').item()
iodides_latpar = np.load('pure_aristotype_joined_iodides_latpar.npy').item()

'''probability dictionaries are activated in next lines'''

oxides_prob = np.load('pure_aristotype_joined_oxides_prob.npy').item()
sulfides_prob = np.load('pure_aristotype_joined_sulfides_prob.npy').item()
selenides_prob = np.load('pure_aristotype_joined_selenides_prob.npy').item()
tellurides_prob = np.load('pure_aristotype_joined_tellurides_prob.npy').item()

fluorides_prob = np.load('pure_aristotype_joined_fluorides_prob.npy').item()
chlorides_prob = np.load('pure_aristotype_joined_chlorides_prob.npy').item()
bromides_prob = np.load('pure_aristotype_joined_bromides_prob.npy').item()
iodides_prob = np.load('pure_aristotype_joined_iodides_prob.npy').item()

names = {'O':'oxides', 'S':'sulfides', 'Se': 'selenides', 'Te':'tellurides',
         'F':'fluorides', 'Cl':'chlorides', 'Br':'bromides', 'I':'iodides'}

def prob(formula = ''):
    cubocta, octa, anion = re.findall('[A-Z][^A-Z]*', formula)
    anion = anion[:-1]
    return globals() [names.get(anion,None) + '_prob'][formula][:,0]

def latpar(formula = ''):
    cubocta, octa, anion = re.findall('[A-Z][^A-Z]*', formula)
    anion = anion[:-1]
    return globals() [names.get(anion,None) + '_latpar'][formula]

formulas = ['LiBeF3', 'LiMgF3', 'LiCaF3','LiSrF3','LiBaF3']

def compare_compounds(formulas = list()):
    fig, axes = plt.subplots(nrows=len(formulas), ncols=1)

    cmap = plt.cm.get_cmap('RdYlBu')
    norm = matplotlib.colors.BoundaryNorm(np.arange(0,1.1,0.1), cmap.N)
    
    for ax, idx, formula in zip(axes.flat, np.arange(0, len(axes.flat),1), formulas):
        nombre = formula[:-1]
        todos = ax.imshow(prob(formula)[:,np.newaxis].T,cmap=cmap,norm = norm)
        
        parameter = latpar(formula)
        #print(parameter)
        parameter = ["%.4f" % parameter[i] for i in range(len(parameter)) if i%5 == 0]
        #print(parameter)
        ax.set_title(nombre+'$_{3}$', fontsize= 16)
        ax.set_xticks(np.arange(0,41,5))
        ax.set_xticklabels(parameter, fontsize = 14)
        ax.set_xlabel('Lattice parameter in ' + '$\AA$', fontsize= 14)
        ax.set_xticks(np.arange(-0.5,41,1),minor=True)
        ax.set_yticks([])
        ax.grid(which='minor', color='black', linewidth = '0.50')
    cb = plt.colorbar(todos, cmap = cmap, norm = norm, ax = list(axes.ravel()), shrink = 0.75, 
                      orientation = 'vertical')
    cb.ax.tick_params(labelsize = 14) #, extend='min')
    cb.set_label('Probability to crystallize as perovskite', fontsize = 14)
    plt.show()
    return 
