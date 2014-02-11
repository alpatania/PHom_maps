# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import pickle
import os, sys
import pandas as pd
from numpy import *
# <codecell>

#test_filtration='code/celegans_resolved_weighted_clique_filtration_100.pck';#sys.argv[1];
#output_bm_file=sys.argv[2];
#fil=pickle.load("/Users/alicepatania/Dropbox/quiver_resolution/data/celegans_resolved_weighted_clique_filtration_100.pck");

f=open("/Users/alicepatania/Desktop/ISI/celegans_weighted_clique_filtration.pck",'rb');
bin_data=f.read()
fil=pickle.loads(bin_data);



# <codecell>

def codimension(simplex1,simplex2,verbose=False):

    if verbose==True:
         ( a, len(a), b, len(b))
    if (simplex1 < simplex2) and (len(simplex2-simplex1)==1):
        return 1;
    if (simplex2 < simplex1) and (len(simplex1-simplex2)==1):
        return -1;
    return 0;

def invert_filtration_dictionary(fil):
    import ast
    inv_dict={};
    
    for c in fil:
        try:
            inv_dict[(int(fil[c][0]),len(set(ast.literal_eval(c))))].append(set(ast.literal_eval(c)));
        except:
            inv_dict[(int(fil[c][0]),len(set(ast.literal_eval(c))))]=[];
            inv_dict[(int(fil[c][0]),len(set(ast.literal_eval(c))))].append(set(ast.literal_eval(c)));
    for l in inv_dict:
        inv_dict[l]=sorted(inv_dict[l], key=lambda x: len(x));
        
    return inv_dict;

def sparse_boundary_matrix_zero(inv_dict,c,k,verbose=False):#c is the filtration step needed
    
    ordered_simplex_list=[];
    Ord=[]
    count=0;
    try:
        ordered_simplex_list.extend(inv_dict[(c,k)]);
        R=len(ordered_simplex_list);
        Ord=list(ordered_simplex_list)#
	#print 'here i am born %d'% c,Ord
    except KeyError:
        return Matrix(ZZ,[]),[]
    try:
        ordered_simplex_list.extend(inv_dict[(c,k+1)]);
	#print 'here i am %d'% c,Ord
    except KeyError:
        return Matrix(ZZ,[]),Ord
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
    bm= Matrix(ZZ,bm);
    return Matrix(ZZ,bm),Ord

def sparse_boundary_matrix(inv_dict,c,k,deltak=Matrix(ZZ,[]),ordered_ksimplex_list=[]):#where deltak is the kth boundary matrix of c-1
    if c==0:
        return sparse_boundary_matrix_zero(inv_dict,c,k,verbose=False)
    elif not deltak:
	print 'I have found that delta %d is empty' %k, deltak;
	if not ordered_ksimplex_list:
	   print 'I have found that there are no %d simplices'%(k-1),ordered_ksimplex_list;
	   return sparse_boundary_matrix_zero(inv_dict,c,k,verbose=False)
    ordered_simplex_list=[];
    new=False;
    Ord=[]
    Ord=list(ordered_ksimplex_list);
    try:
        ordered_simplex_list.extend(inv_dict[(c,k)]);
        R=len(ordered_simplex_list);
        Ord.extend(ordered_simplex_list)
        if deltak:
            D=Matrix(ZZ,vstack([deltak, Matrix(ZZ,zeros((R,shape(deltak)[1])))]))
	    row_sum_D=[sum([x]) for x in D.rows()]
	else:
	    #print 'uso il numero a casaccio'
	    row_sum_D=[0]*(len(Ord))
    except KeyError:
        R=0;
        print 'no new %d simplex'%(k-1) 
        D=deltak;
	row_sum_D=[sum([x]) for x in D.rows()]

    try:
        ordered_simplex_list.extend(inv_dict[(c,k+1)]);
        if not deltak:
            new=True;
    except KeyError:
        if not deltak:
            return Matrix(ZZ,[]),Ord
        else:
            return Matrix(ZZ,D),Ord
    C=len(ordered_simplex_list)-R;
    ordered_ksimplex_list.extend(ordered_simplex_list)
    ordered_simplex_series=pd.Series(ordered_ksimplex_list,index=range(len(ordered_ksimplex_list)));
    del ordered_ksimplex_list;
    del ordered_simplex_list;
    bm=zeros((len(Ord),C));
    for i in ordered_simplex_series.index:
        if len(ordered_simplex_series[i])==k:
            for j in ordered_simplex_series.index[i:]:
                cod=codimension(ordered_simplex_series[i], ordered_simplex_series[j]);
                if cod==1:
		    piu=list(ordered_simplex_series[j]-ordered_simplex_series[i])
		    s=sorted(ordered_simplex_series[j])
		    esp=s.index(piu[0])
		    if esp%2==0:
		    	bm[i,j-len(Ord)]=1;
		    else:
			bm[i,j-len(Ord)]=-1;
		
    if new:
	bm= Matrix(ZZ,bm);
        BM=bm;
        del bm
    else:
        BM=hstack([D,bm]);
	BM= Matrix(ZZ,BM);
        del bm,D
    return Matrix(ZZ,BM),Ord;
        
    
def Laplacian(inv_fil,c,k,deltak=Matrix(ZZ,[]),Ord_k=[],deltak1=Matrix(ZZ,[]),Ord_k1=[],save_boundary=True):
    if k==0:
	Dk1,Ordk1=sparse_boundary_matrix(inv_fil,c,k+1,deltak1,Ord_k1)
	print (c,k+1)
	if Dk1:
	   Dk1=Matrix(ZZ,Dk1); 
	   L=Dk1*(Dk1.transpose())
	   return L,Matrix(ZZ,[]),Dk1,[],Ordk1
	else:
	   return Matrix(ZZ,[]),Matrix(ZZ,[]),Matrix(ZZ,[]),[],[]
    Dk,Ordk=sparse_boundary_matrix(inv_fil,c,k,deltak,Ord_k)
    print (c,k)
    Dk1,Ordk1=sparse_boundary_matrix(inv_fil,c,k+1,deltak1,Ord_k1)
    print (c,k+1) 
    if not Dk1:
        if not Dk:
            return Matrix(ZZ,[]),Matrix(ZZ,[]),Matrix(ZZ,[]),[],[]
        else:
	    Dk=Matrix(ZZ,Dk)
            L=(Dk.transpose())*Dk
    else:
	#print 'cosa fai? ecco il Laplaciano passo passo:\n D1\n',Dk,'\n la sua trasposta\n',Dk.transpose(),'\n la parte di D1\n',(Dk.transpose())*Dk,'\n la parte di D2:\n',Dk1*(Dk1.transpose());
	Dk=Matrix(ZZ,Dk)
	Dk1=Matrix(ZZ,Dk1)
        L=(Dk.transpose())*Dk+Dk1*(Dk1.transpose())
    if save_boundary==True:
        return L,Dk,Dk1,Ordk,Ordk1
    else:
        return L

def Proiettore(L):
    LT=L.transpose()
    PP=LT*L
    invPP=(PP).inverse()
    P=((L*invPP)*LT)
    return P

# <codecell>

inv_fil=invert_filtration_dictionary(fil)

# <codecell>


def right_kernel_space(L):
    if L==0:
        return Matrix(ZZ,[])
    u, s, vh = scipy.linalg.svd(L)
    null_mask = (s <= eps)
    null_space = scipy.compress(null_mask, vh, axis=0)
    if null_space.any():
        return  Matrix(ZZ,null_space.transpose())
    else:
        return Matrix(ZZ,[])

def column_space(L):
    #print('COLUMN SPACE - I am using toll:', eps)
    if L==0:
        return Matrix(ZZ,[])
    u, s, vh = scipy.linalg.svd(L)
    column_mask = (s >= eps)
    column = scipy.compress(column_mask, u, axis=1)
    return  Matrix(ZZ,column)

def H(LambdaD, BD, BBD,verbose=False):
    k=[1,1,1]
    if verbose==True:
	print parent(LambdaD).dims(), parent(BD).dims(), parent(BBD).dims();
    if not LambdaD:
        k[0]=0
    if not BD:
        k[1]=0
    if not BBD:
        k[2]=0
    if verbose==True:
	print k;
    if k==[1,1,1]:
        M=hstack([hstack([LambdaD,BD]),BBD])
    elif k==[1,0,1]:
        M=hstack([LambdaD,BBD])
    elif k==[1,1,0]:
        M=hstack([LambdaD,BD])
    elif k==[1,0,0]:
        return LambdaD
    elif k==[0,1,1]:
        M=hstack([BD,BBD])
    elif k==[0,1,0]:
        return BD
    elif k==[0,0,1]:
        return BBD
    return M

# <codecell>

def Persistent_Homology_maps(k):
    from numpy import zeros
    n=sorted(inv_fil.keys())[-1][0]
    homCD={}
    LapC,D1C,D2C,O1C,O2C=Laplacian(inv_fil,0,k)
    LambdaC=( Matrix(ZZ,LapC )).kernel()
    LambdaC=LambdaC.basis_matrix()
    del LapC
    for c in range(1,n+1):
        print c 
        LapD,D1D,D2D,O1D,O2D=Laplacian(inv_fil,c,k,D1C,O1C,D2C,O2C)
        if [LapD,D1D,D2D]!=[0,0,0]:
            if [D1C,D2C]==[0,0]:
                LambdaD=(Matrix(ZZ,LapD)).kernel()
		LambdaD=LambdaD.basis_matrix()
                LambdaC=(Matrix(ZZ,LambdaD)).transpose()
                D1C=D1D
                if D2D!=0:
                    D2C=D2D
            else:
                HD= Matrix(ZZ,[[],[]])
                while HD.is_square()==False:
                    LD=(Matrix(ZZ,LapD)).column_space()
		    LD=LD.basis_matrix();
                    LambdaD=(Matrix(ZZ,LapD)).kernel()
		    LambdaD=(Matrix(ZZ,LambdaD.basis_matrix())).transpose()
		    PD=Proiettore(LD.transpose())

		    BD=(Matrix(ZZ,D2D)).column_space()
		    BD=(BD.basis_matrix()).transpose()

                    fuffa=(D1D).transpose()
                    BBD=(Matrix(ZZ,(fuffa))).column_space()
		    BBD=(BBD.basis_matrix()).transpose()

                    IDC=eye(shape(D1C)[1]);
                    if shape(D1D)[1]>shape(D1C)[1]:
			r=shape(D1D)[1]-shape(D1C)[1]
                        ZERO=zeros((r,shape(D1C)[1]));
                        F1C=Matrix(ZZ,vstack([IDC,ZERO]));
                    else:
                        F1C=Matrix(ZZ,IDC);
		    
                    HD=H(LambdaD,BD,BBD)
		    HD=Matrix(ZZ,HD)
		    if  HD.is_square():
			print 'HD is square';
			HD=Matrix(ZZ,HD).inverse()
		    else:
                        print 'HD=LambdaD|BD|BBD',(LambdaD.nrows(),LambdaD.ncols()),(BD.nrows(),BD.ncols()),(D1D.nrows(),D1D.ncols()),(BBD.nrows(),BBD.ncols());
                        raise ValueError ('ERROR HD NOT SQUARE')

		if not LambdaC:
		    homCD[c-1]=Matrix(ZZ,[])
		    LambdaC=LambdaD
                    D1C=D1D
               	    if D2D!=0:
                    	D2C=D2D
		else:
		    #print 'HD,PD,F1C,LambdaC',(HD.nrows(),HD.ncols()),(PD.nrows(),PD.ncols()), (F1C.nrows(),F1C.ncols()),(LambdaC.nrows(),LambdaC.ncols());
                    HOM=HD*(PD*(F1C*LambdaC))
		    homCD[c-1]=HOM[:HOM.ncols()][:]
                    LambdaC=LambdaD
                    D1C=D1D
                    if D2D!=0:
                  	  D2C=D2D
        O1C=O1D#
        O2C=O2D#
    IDD=Matrix.identity(shape(D1D)[1]);
    HOM=HD*(PD*(IDD*LambdaD))
    if not HOM:
	homCD[c]=HOM[:HOM.ncols()][:]
    #print 'THIS IS THE LAST STEP\n HD\n',HD,'\n PD\n',PD,'\n IDD\n',IDD,'\n LambdaD \n',LambdaD
    del D1C,D2C,O1C,O2C,LapD,D1D,D2D,O1D,O2D,LambdaC
    print 'Done.' 
    return homCD

# <codecell>

#inv_fil={(0,1):[{1},{2}],(0,2):[{1,2}],(1,1):[{3}],(1,2):[{2,3}],(2,2):[{1,3},{5,6},{5,3},{6,3}],(3,1):[{4},{5},{6}],(2,3):[{1,2,3}]}
#print inv_fil;
#D=Persistent_Homology_maps(4)