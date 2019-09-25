# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import pickle
import os, sys
import networkx as nx
import pandas as pd

# <codecell>

sys.path.append('./');
#import Holes as ho

# <codecell>

notebook_mode=False;
if notebook_mode==True:
    test_filtration='../data/celegans/celegans_weighted_clique_filtration.pck'
    output_bm_file='/Users/lordgrilo/Desktop/celegans_weighted_clique_boundary_hdm5';
else:
    if len(sys.argv)<=3:
        print 'Inputs required: \n1) filtration file \n2) tag name for output files \n';
    test_filtration=sys.argv[1];
    output_bm_file=sys.argv[2];

fil=pickle.load(open(test_filtration,'r'));

# <codecell>


# <headingcell level=2>

# Sparse representation of boundary matrix

# <markdowncell>

# The problem here is that for large filtrations, the boundary matrix easily becomes extremely large. 
# It is therefore necessary to have a viable sparse representation. 
# For some reason the pandas rep. has problems with the pickling part.   
# I see two solutions: 
# 
# 1. use HDF5 for storing  
# 2. use pure numpy-ness and see whether it works.  

# <codecell>

def codimension(simplex1,simplex2,verbose=False):

    if verbose==True:
        print a, len(a), b, len(b)
    if (simplex1 < simplex2) and (len(simplex2-simplex1)==1):
        return 1;
    if (simplex2 < simplex1) and (len(simplex1-simplex2)==1):
        return -1;
    return 0;
                
def invert_filtration_dictionary(holes_filtration):
    import ast
    inv_dict={};
    
    for c in fil:
        try:
            inv_dict[int(fil[c][0])].append(set(ast.literal_eval(c)));
        except:
            inv_dict[int(fil[c][0])]=[];
            inv_dict[int(fil[c][0])].append(set(ast.literal_eval(c)));
    for l in inv_dict:
        inv_dict[l]=sorted(inv_dict[l], key=lambda x: len(x));
        
    return inv_dict;

def boundary_matrix(inv_dict,verbose=False):

    ordered_simplex_list=[];
    count=0;
    n={};
    for k in sorted(inv_dict):
        ordered_simplex_list.extend(inv_dict[k]);
        n[k]=len(ordered_simplex_list);
        print k, len(ordered_simplex_list);
    ordered_simplex_series=pd.Series(ordered_simplex_list,index=range(len(ordered_simplex_list)));    
    
    del ordered_simplex_list;
    print 'Simplex ordering complete. Building boundary matrix now.';
    bm=pd.DataFrame(index=ordered_simplex_series.index, columns=ordered_simplex_series.index);
    for i in ordered_simplex_series.index:
        if i%100==0:
            print i;
        for j in ordered_simplex_series.index[i:]:
            cod=codimension(ordered_simplex_series[i], ordered_simplex_series[j]);
            if cod==1:
                bm[j][i]=1;
                if verbose==True:
                    print i,j, cod,ordered_simplex_series[i], ordered_simplex_series[j];

            elif cod==-1:
                bm[i][j]=1;
                if verbose==True:
                    print i,j, cod,ordered_simplex_series[i], ordered_simplex_series[j];
    #print 'Sparsifying boundary matrix representation...'
    bm=bm.fillna(0)#.to_sparse(fill_value=0)
    print 'Done.';
    return bm,n, ordered_simplex_series;

def sparse_boundary_matrix(inv_dict,verbose=False):
    import scipy
    from scipy.sparse import csc_matrix, lil_matrix;
    import numpy as np
    
    ordered_simplex_list=[];
    count=0;
    n={};
    for k in sorted(inv_dict):
        ordered_simplex_list.extend(inv_dict[k]);
        n[k]=len(ordered_simplex_list);
        print k, len(ordered_simplex_list);
    ordered_simplex_series=pd.Series(ordered_simplex_list,index=range(len(ordered_simplex_list)));    
    L=len(ordered_simplex_list);
    del ordered_simplex_list;
    print 'Simplex ordering complete. Building boundary matrix now.';
    bm=lil_matrix((L,L));
    for i in ordered_simplex_series.index:
        if i%100==0:
            print i;
        for j in ordered_simplex_series.index[i:]:
            cod=codimension(ordered_simplex_series[i], ordered_simplex_series[j]);
            if cod==1:
                bm[i,j]=1;
                if verbose==True:
                    print i,j, cod,ordered_simplex_series[i], ordered_simplex_series[j];

            elif cod==-1:
                bm[j,i]=1;
                if verbose==True:
                    print i,j, cod,ordered_simplex_series[i], ordered_simplex_series[j];
    #print 'Sparsifying boundary matrix representation...'
    bm=csc_matrix(bm);
    print 'Done.';
    return bm,n;
            

# <codecell>

print 'Creating inverse filtration dictionary.'
invf=invert_filtration_dictionary(fil)
#print 'Creating boundary matrix. ';
#bound, n=sparse_boundary_matrix(invf,verbose=True)

# <markdowncell>

# Let's see if we can fix this problem of saving.  
# Now the output of the boundary matrix calculation is a sparse pd.DataFrame, filled and zipped with zeros.   
# This does not solve the problem completely but it is a beginning.  

# <codecell>

print 'Saving filtration dictionary'+output_bm_file+'.pck';
pickle.dump(invf, open(output_bm_file+'_fil_numbering.pck','wb'))
#print 'Saving boundary matrix dictionary to '+output_bm_file+'_fil_numbering.pck';
#pickle.dump(n, open(output_bm_file+'_fil_numbering.pck','wb'))

