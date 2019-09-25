# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

import pickle
import os, sys
import pandas as pd
import numpy as np
import scipy.sparse as sp
from numpy import *

# FILTRATION GIVEN BY THE USER, MUST BE A DICTIONARY WITH KEYS ALL THE SIMPLICES IN THE SIMPLICIAL COMPLEX AND AS VALUE THE STEP OF THE FILTRATION WHEN THEY ARE CREATED (AND A WEIGHT NOT NECESSARY)
#test_filtration=sys.argv[1];
#fil=pickle.load(open(test_filtration,'r'));

# THE TEST FILTRATION IN THE C-ELEGANS
f=open("./rand1_single_step_fil.pck",'rb');
bin_data=f.read()
fil=pickle.loads(bin_data);

# FUNCTION CODIMENSION: GIVEN 2 SIMPLICES RETURNS 1 AND -1 IF THAT IS THAT IS THEIR CODIMENSION AND 0 OTHERWISE
def codimension(simplex1,simplex2,verbose=False):

    if verbose==True:
         ( a, len(a), b, len(b))
    if (simplex1 < simplex2) and (len(simplex2-simplex1)==1):
        return 1;
    if (simplex2 < simplex1) and (len(simplex1-simplex2)==1):
        return -1;
    return 0;

# FUNCTION: INVERTS THE DICTIONARY GIVEN SO THAT THE KEYS ARE COUPLES WHERE THE FIRST NUMBER IS THE FILTRATION STEP AND THE SECOND IS THE DIMENSION+1 OF THE SIMPLICES IN THE CORRESPONDING VALUE

def invert_filtration_dictionary(fil):
    import ast
    inv_dict={};
    
    for c in fil: # TAKES A SIMPLEX IN THE KEYS OF FILTRATION
        try:
            inv_dict[(int(fil[c][0]),len(set(ast.literal_eval(c))))].append(set(ast.literal_eval(c))); 
        except: # IF THE KEY DOESN'T ALREADY EXIST
            inv_dict[(int(fil[c][0]),len(set(ast.literal_eval(c))))]=[]; 
            inv_dict[(int(fil[c][0]),len(set(ast.literal_eval(c))))].append(set(ast.literal_eval(c))); 
    for l in inv_dict:
        inv_dict[l]=sorted(inv_dict[l], key=lambda x: len(x));
        
    return inv_dict;

# FUNCTION SPARSE_BOUNDARY_MATRIX_ZERO: RETURNS THE kth-BOUNDARY MATRIX FOR THE SIMPLICIAL COMPLEX AT STEP c IN THE FILTRATION WHERE THE FIRST k-1-SIMPLEX IS ADDED TO THE FILTRATION (YOU DON'T NEED TO ADD THE PREVIOUS kth-BOUNDARY MATRIX)
def sparse_boundary_matrix_zero(inv_dict,c,k,verbose=False): 
    if k==0:
	return zeros((1,len(inv_dict[(c,k+1)])))
    ordered_simplex_list=[];
    try:
        ordered_simplex_list.extend(inv_dict[(c,k)]); 
        R=len(ordered_simplex_list);
    except KeyError:
        return []# IF THERE ARE NO (k-1)-SIMPLICES TO ADD, SINCE THERE WHERE NONE IN THE PREVIOUS STEPS EITHER IT RETURN AN EMPTY MATRIX AND AN EMPTY LIST
    try:
        ordered_simplex_list.extend(inv_dict[(c,k+1)]); 
    except KeyError:
        return []# IF THERE ARE NO k-SIMPLICES TO ADD, THE MATRIX HAS NO COLUMNS SO IT RETURNS AN EMPTY MATRIX AND THE LIST OF (k-1)-SIMPLICES ADDED AT THIS STEP
    C=len(ordered_simplex_list)-R;
    ordered_simplex_series=pd.Series(ordered_simplex_list,index=range(len(ordered_simplex_list))); 
    del ordered_simplex_list;
    bm=zeros((R,C)); 
    for i in ordered_simplex_series.index:
        if len(ordered_simplex_series[i])==k:
            for j in ordered_simplex_series.index[i:]: 
                cod=codimension(ordered_simplex_series[i], ordered_simplex_series[j]);
                if cod==1:
		    piu=list(ordered_simplex_series[j]-ordered_simplex_series[i])
		    s=sorted(ordered_simplex_series[j])
		    esp=s.index(piu[0])
		    if esp%2==0:
		    	bm[i,j-R]=1;
		    else:
			bm[i,j-R]=-1;
                    if verbose==True:
                        print (i,j, cod,ordered_simplex_series[i], ordered_simplex_series[j]);
    bm= matrix(bm);
    return bm

inv_fil=invert_filtration_dictionary(fil)
print inv_fil.keys();
for i in range(sorted(inv_fil.keys())[-1][1]):# 6 IS CHOSEN AS THE MAX k for H_k
    print '\n', i,'\n';
    delta=sparse_boundary_matrix_zero(inv_fil,1,i);
    delta=sp.dok_matrix(delta); 
    pickle.dump(delta, open('boundary%d.pck' %i,'w')) 



#inv_fil={(0,1):[{1},{2}],(0,2):[{1,2}],(1,1):[{3}],(1,2):[{2,3}],(2,2):[{1,3},{5,6},{5,3},{6,3}],(2,1):[{4},{5},{6}],(2,3):[{1,2,3}], (3,1):[{7}],(3,2):[{7,5},{7,3},{7,6}],(3,3):[{7,5,3},{7,5,6},{7,3,6},{5,3,6}],(4,4):[{7,5,6,3}]}
