#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 17:45:29 2017

@author: iG
"""
import pandas as pd
import pymatgen as pg
import numpy as np
import os
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

def drop_characters(in_str):
    '''
    Elimina los caracteres numericos y los signos + y -.
    Esta funcion se tomo de stackoverflow.com
    Parametros: string
    Regresa: string sin digitos y signos +-
    '''
    char_list = "1234567890+-"
    for char in char_list:
        in_str = in_str.replace(char, "")

    return in_str

def move_element(odict, thekey, newpos):
    '''
    Esta funcion se tomo de stackoverflow.com
    Tiene como objetivo cambiar de posicion la clave de un diccionario.
    Parametros:
        odict: diccionario
        thekey: la clave que se quiere mover
        newpos: posicion a la que se quiere mover
    Regresa:
        diccionario con la clave en la posicion requerida        
    '''
    
    odict[thekey] = odict.pop(thekey)
    i = 0
    for key, value in odict.items():
        if key != thekey and i >= newpos:
            odict[key] = odict.pop(key)
        i += 1
    return odict

def wyckoff_occupation(ruta='/home/ivan/Alles/',archivo='1010902'):
    
    '''
    Esta funcion permite encontrar tanto todas las multiplicidad y las etiquetas de los sitios
    de Wyckoff. Se vale de la libreria de pymatgen para conseguir este fin
    Parametros:
        ruta(string): Direccion donde se encuentra el archivo cif.
        archivo(string): Nombre del archivo cif.
    Regresa:
        Diccionario de los elementos de la formula con un array de todas las multiplicidades 
        y sitios de Wyckoff
    '''
    
    estructura=pg.Structure.from_file(str(ruta)+str(archivo)+'.cif')
    eqpos=SpacegroupAnalyzer(estructura).get_symmetry_dataset()

    df=pd.DataFrame({'eq_atoms': eqpos['equivalent_atoms'],'Wyckoff': eqpos['wyckoffs']})
    
    ocupacion=[]
    clave={}    
    for i in range(len(estructura)):
        ocupacion.append(estructura[i].species_string.split(','))
        #print(ocupacion)
        for occ in range(len(estructura[i].species_string.split(','))):
            especie=estructura[i].species_string.split(',')[occ].split(':')[0].lstrip(' ')
            clave.setdefault(i,[]).append(especie)
      
    df['elemento']=drop_characters(df['eq_atoms'].map(clave))
    
    df=df.groupby(['Wyckoff','eq_atoms']).size().reset_index(name='multiplicidad')
    
    df=df[['multiplicidad','Wyckoff','eq_atoms']]
    df=df.drop_duplicates()
    df=df[['eq_atoms','Wyckoff']]
    df=df.reset_index(drop=True)
    df=df.drop_duplicates()
    
    dicc={}
    for item in range(len(ocupacion)):
        
            frac={}
            if len(ocupacion[item]) == 1:
                try:
                    frac[ocupacion[item][0].split(':')[0].lstrip(' ')] = np.round(float(ocupacion[item][0].split(':')[1].lstrip(' ')),5)
                except:
                    frac[ocupacion[item][0]] = 1
            
            else:
                for Z in ocupacion[item]:
                    frac[Z.split(':')[0].lstrip(' ')] = np.round(float(Z.split(':')[1].lstrip(' ')),5)
                
            dicc[item]=frac
       
    df=df[['Wyckoff', 'eq_atoms']]
    df['eq_atoms']=df['eq_atoms'].map(dicc)
    
    diccionario={}
    for row in range(len(df)):
        diccionario[row]={df['Wyckoff'][row]:df['eq_atoms'][row]}

    return diccionario


def wyckoff_positions(ruta='/home/ivan/Alles/',archivo='1010902'):
    
    '''
    Esta funcion permite encontrar tanto todas las multiplicidad y las etiquetas de los sitios
    de Wyckoff. Se vale de la libreria de pymatgen para conseguir este fin
    Parametros:
        ruta(string): Direccion donde se encuentra el archivo cif.
        archivo(string): Nombre del archivo cif.
    Regresa:
        Diccionario de los elementos de la formula con un array de todas las multiplicidades 
        y sitios de Wyckoff
    '''
    
    estructura=pg.Structure.from_file(str(ruta)+str(archivo)+'.cif')
    eqpos=SpacegroupAnalyzer(estructura).get_symmetry_dataset()
    archivo=SpacegroupAnalyzer(estructura).get_conventional_standard_structure()
    text=str(archivo)
    
    sitios=int(text.split('\n')[4].split('(')[1].split(')')[0])
            
    abc=[float(item) for item in list(filter(None,text.split('\n')[2].split(':')[1].split(' ')))]
    angles=[float(item) for item in list(filter(None,text.split('\n')[3].split(':')[1].split(' ')))]
    lista=text.split('\n')[-sitios:]
        
    newlist=[list(filter(None,line.split(' '))) for line in lista]
    newlist = [[item[0]] + [str(item[1:-3])] +  item[-3:] for item in newlist]
    newlist=np.asarray(newlist)
    #print(archivo)    
    motif=pd.DataFrame(newlist)[[1,2,3,4]]
    motif[2]=[float(i) for i in motif[2].values]
    motif[3]=[float(i) for i in motif[3].values]
    motif[4]=[float(i) for i in motif[4].values]

    volumen=abc[0]*abc[1]*abc[2]*np.sqrt(1-(np.cos(np.deg2rad(angles[0])))**2-(np.cos(np.deg2rad(angles[1])))**2-(np.cos(np.deg2rad(angles[2])))**2+2*np.cos(np.deg2rad(angles[0]))*np.cos(np.deg2rad(angles[1]))*np.cos(np.deg2rad(angles[2])))

    
    df=pd.DataFrame({'eq_atoms': eqpos['equivalent_atoms'],'Wyckoff': eqpos['wyckoffs']})
    
    df=df.groupby(['Wyckoff','eq_atoms']).size().reset_index(name='multiplicidad')
    df=df[['multiplicidad','Wyckoff','eq_atoms']]
    df=df.drop_duplicates()
    df=df[['eq_atoms','Wyckoff']]
    df=df.reset_index(drop=True)
    df=df.drop_duplicates()
    #print(df)
    #motif=pd.read_csv('motif.csv', header=None, delim_whitespace=True)[[1,2,3,4]]
    motif.columns = np.arange(len(motif.columns))
    df=pd.concat([df,(motif.loc[list(df['eq_atoms'].values),1:]).reset_index(drop=True)], axis=1, join='inner')
    #print(df)
    idx = sorted(df['eq_atoms'].values) + [len(motif)]
    #print(idx)
    diccionario={}
    for row in range(len(df)):
        #vector=np.around(df.iloc[row,2:].values.astype(np.double),decimals=4)
        init = df['eq_atoms'][row].item()
        finit = idx[idx.index(init) + 1]
        #print(init,finit)
        diccionario[row]={df['Wyckoff'][row] : motif.iloc[init:finit,-3:].values}
    
    return diccionario,motif, angles, abc

def wyckoff_finder(ruta='/home/ivan/Alles/',archivo='1010902'):
    
    '''
    Esta funcion permite encontrar tanto todas las multiplicidad y las etiquetas de los sitios
    de Wyckoff. Se vale de la libreria de pymatgen para conseguir este fin
    Parametros:
        ruta(string): Direccion donde se encuentra el archivo cif.
        archivo(string): Nombre del archivo cif.
    Regresa:
        Diccionario de los elementos de la formula con un array de todas las multiplicidades 
        y sitios de Wyckoff
    '''
    
    estructura=pg.Structure.from_file(str(ruta)+str(archivo)+'.cif')
    eqpos=SpacegroupAnalyzer(estructura).get_symmetry_dataset()

    df=pd.DataFrame({'eq_atoms': eqpos['equivalent_atoms'],'Wyckoff': eqpos['wyckoffs']})
    
    clave={k:v for k,v in enumerate([[occ.split(':')[0].lstrip(' ') for occ in estructura[i].species_string.split(',')] for i in range(len(estructura))])}
        
    df['elemento']=drop_characters(df['eq_atoms'].map(clave))
    df=df.groupby(['Wyckoff','eq_atoms']).size().reset_index(name='multiplicidad')
    df=df[['multiplicidad','Wyckoff','eq_atoms']]
    df=df.drop_duplicates()
    df['mult_wyc']=df['multiplicidad'].map(str)+df['Wyckoff']
    df=df[['eq_atoms','mult_wyc']]
    df=df.reset_index(drop=True)
    df=df.drop_duplicates()
    
    dicc1={k:[v] for k,v in zip(df['eq_atoms'],df['mult_wyc'])}
       
    dicc={}
    for keys in dicc1.keys():
        for el in range(len(clave[keys])):
            dicc.setdefault(clave[keys][el].rstrip('0123456789+-.'),[]).append(dicc1[keys])
    return dicc

