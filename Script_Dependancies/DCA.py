# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 13:19:49 2020

@author: glham
"""

import numpy
import scipy
import scipy.spatial
# =============================================================================
# Direct Coupling Analysis (DCA)
# 
# function dca(infile , outputfile)
# 
# INPUTS: 
#   infile  - file containing the FASTA alignment
#   outputfile - file for dca results. The file is composed by N(N-1)/2 
#                (N = length of the sequences) rows and 4 columns: 
#                residue i (column 1), residue j (column 2),
#                MI(i,j) (Mutual Information between i and j), and 
#                DI(i,j) (Direct Information between i and j).
#                Note: all insert columns are removed from the alignment.
# 
# SOME RELEVANT VARIABLES:
#   N        number of residues in each sequence (no insert)
#   M        number of sequences in the alignment
#   Meff     effective number of sequences after reweighting
#   q        equal to 21 (20 aminoacids + 1 gap)
#   align    M x N matrix containing the alignmnent
#   Pij_true N x N x q x q matrix containing the reweigthed frequency
#            counts.
#   Pij      N x N x q x q matrix containing the reweighted frequency 
#            counts with pseudo counts.
#   C        N(q-1) x N(q-1) matrix containing the covariance matrix.
# 
# 
# Copyright for this implementation: 
#             2011/12 - Andrea Pagnani and Martin Weigt
#                       andrea.pagnani@gmail.com 
#                       martin.weigt@upmc.fr
# 
# Permission is granted for anyone to copy, use, or modify this
# software and accompanying documents for any uncommercial
# purposes, provided this copyright notice is retained, and note is
# made of any changes that have been made. This software and
# documents are distributed without any warranty, express or
# implied. All use is entirely at the user's own risk.
# 
# Any publication resulting from applications of DCA should cite:
# 
#     F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander, 
#     R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), Direct-coupling
#     analysis of residue co-evolution captures native contacts across 
#     many protein families, Proc. Natl. Acad. Sci. 108:E1293-1301.
# =============================================================================
# =============================================================================
# 
# class chisq_functions:
#     def meaner(k,vec): #function written as if to find the best while having split k
#         return len(vec[0:k])*numpy.var(vec[0:k])+len(vec[k:])*numpy.var(vec[k:])
#     def mean(arrays): #works with multiple arrays, easier to generalize if the splitting is done outside the function
#         return numpy.sum([len(array)*numpy.var(array) for array in arrays])
#     string={}
#     string['meaner']=meaner
#     string['mean']=mean
# =============================================================================

def dca(infile,outfile,pseudocount_weight=.5,theta=.2,writemode="w+"):
    #pseudocount weight: relative weight of pseudocount
    #theta: threshold for sequence id in reweighting
    
    N,M,q,align=return_alignment(infile)
    Pij_true,Pi_true,Meff=Compute_True_Frequencies(align,M,N,q,theta)

    Pij,Pi=with_pc(Pij_true,Pi_true,pseudocount_weight,N,q)

    C=Compute_C(Pij,Pi,N,q)
    invC=numpy.linalg.inv(C)
    #good through here
    file=open(outfile,writemode)
    print('### N = {0} M = {1} Meff = {2:.2f} q = {3}\n'.format(N,M,Meff,q))
    Compute_Results(Pij,Pi,Pij_true,Pi_true,invC,N,q,file)
    file.close()

def return_alignment(infile):
    #reads alignment from inputfile, removes inserts and converts into numbers
    align_full=fastaread(infile)
    M=len(align_full)
    ind=[(x!='.')&(x==x.upper()) for x in align_full[0][1]]
    N=sum(ind)
    Z=numpy.zeros((M,N))
    for i in range(1,M+1):
        counter=0
        for j in range(1,len(ind)+1):
            if ind[j-1]:
                counter=counter+1
                #print(i,j)
                Z[i-1,counter-1]=letter2number(align_full[i-1][1][j-1])
    q=numpy.max(Z)
    return [N,M,int(q),Z.astype(int)]

class fasta:
    def __init__(self,header,sequence):
        self.Header=header
        self.Sequence=sequence
    

def fastaread(infile):
    fasta=open(infile,"r")
    a=list(fasta.readlines())
    b=[]
    for i in range(len(a[::2])):
        header=a[::2][i][1:-1]
        sequence=a[1::2][i][:-1]
        temp=[header,sequence]
        b.append(temp)
    fasta.close()
    return b
    

def Compute_Results(Pij,Pi,Pij_true,Pi_true,invC,N,q,file):
    #computes and prints mutual and direct infromations
    N=int(N)
    q=int(q)
    for i in range(1,N):
        for j in range(i+1,N+1):
            MI_true,si_true,sj_true=calculate_mi(i,j,Pij_true,Pi_true,q)
            W_mf=ReturnW(invC,i,j,q)
            DI_mf_pc=bp_link(i,j,W_mf,Pi,q)
            file.write('{0} {1} {2} {3}\n'.format(i,j,MI_true,DI_mf_pc))

def Compute_True_Frequencies(align,M,N,q,theta):
    #computes reweighted frequency counts
    W=numpy.ones(M)
    if theta>0.0:
        W=(1./(1+sum(scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(align,'hamming')<theta))))
    Meff=sum(W)
    N,q=int(N),int(q)
    Pij_true=numpy.zeros((N,N,q,q))
    Pi_true=numpy.zeros((N,q))
    
    for j in range(1,M+1):
        for i in range(1,N+1):
            Pi_true[i-1,align[j-1,i-1].astype(int)-1]=Pi_true[i-1,align[j-1,i-1].astype(int)-1]+W[j-1]
    Pi_true=Pi_true/Meff
    
    for l in range(1,M+1):
        for i in range(1,N):
            for j in range(i+1,N+1):
                Pij_true[i-1,j-1,align[l-1,i-1].astype(int)-1,align[l-1,j-1].astype(int)-1]=Pij_true[i-1,j-1,align[l-1,i-1].astype(int)-1,align[l-1,j-1].astype(int)-1]+W[l-1]
                Pij_true[j-1,i-1,align[l-1,j-1].astype(int)-1,align[l-1,i-1].astype(int)-1]=Pij_true[i-1,j-1,align[l-1,i-1].astype(int)-1,align[l-1,j-1].astype(int)-1]
    Pij_true=Pij_true/Meff
    
    scra=numpy.identity(q)
    for i in range(1,N+1):
        for alpha in range(1,q+1):
            for beta in range(1,q+1):
                Pij_true[i-1,i-1,alpha-1,beta-1]=Pi_true[i-1,alpha-1]*scra[alpha-1,beta-1]
    return [Pij_true,Pi_true,int(Meff)]
                
def letter2number(a):
    switch={'-':1,
        'A':2,
        'C':3,
        'D':4,
        'E':5,
        'F':6,
        'G':7,
        'H':8,
        'I':9,
        'K':10,
        'L':11,
        'M':12,
        'N':13,
        'P':14,
        'Q':15,
        'R':16,
        'S':17,
        'T':18,
        'V':19,
        'W':20,
        'Y':21,}
    try:
        return switch[a]
    except KeyError:
        return 1

def with_pc(Pij_true, Pi_true, pseudocount_weight,N,q):
    #adds pseudocount
    N=int(N)
    q=int(q)
    Pij=(1.-pseudocount_weight)*Pij_true + pseudocount_weight/q/q*numpy.ones((N,N,q,q))
    Pi=(1.-pseudocount_weight)*Pi_true + pseudocount_weight/q*numpy.ones((N,q))

    scra=numpy.identity(q)
    
    for i in range(1,N+1):
        for alpha in range(1,q+1):
            for beta in range(1,q+1):
                Pij[i-1,i-1,alpha-1,beta-1]=(1.-pseudocount_weight)*Pij_true[i-1,i-1,alpha-1,beta-1]+pseudocount_weight/q*scra[alpha-1,beta-1]
    return [Pij,Pi]

def mapkey(i,alpha,q):
    q=int(q)
    return (q-1)*(i-1)+numpy.array(alpha)

def Compute_C(Pij,Pi,N,q):
    #computes correlation matrix
    N=int(N)
    q=int(q)
    C=numpy.zeros((N*(q-1),N*(q-1)))
    
    for i in range(1,N+1):
        for j in range(1,N+1):
            for alpha in range(1,q):
                for beta in range(1,q):
                    C[mapkey(i,alpha,q)-1,mapkey(j,beta,q)-1]=Pij[i-1,j-1,alpha-1,beta-1]-Pi[i-1,alpha-1]*Pi[j-1,beta-1]                    
    return C

def calculate_mi(i,j,P2,P1,q):
    #computes mutual information between columns i and j
    M=0.
    for alpha in range(1,q+1):
        for beta in range(1,q+1):
            if P2[i-1,j-1,alpha-1,beta-1]>0:
                M=M+P2[i-1,j-1,alpha-1,beta-1]*numpy.log(P2[i-1,j-1,alpha-1,beta-1]/P1[i-1,alpha-1]/P1[j-1,beta-1])
    
    s1=0.
    s2=0.
    for alpha in range(1,q+1):
        if P1[i-1,alpha-1]>0:
            s1=s1-P1[i-1,alpha-1]*numpy.log(P1[i-1,alpha-1])
        if P1[j-1,alpha-1]>0:
            s2=s2-P1[j-1,alpha-1]*numpy.log(P1[j-1,alpha-1])
    return [M,s1,s2]

def ReturnW(C,i,j,q): #has problem
    #extracts coupling matrix for columns i and j
    q=int(q)
    W=numpy.ones((q,q))
    W[0:q-1,0:q-1]=numpy.exp(-C[mapkey(i,range(1,q),q).min()-1:mapkey(i,range(1,q),q).max(),mapkey(j,range(1,q),q).min()-1:mapkey(j,range(1,q),q).max()])
    
    return W

def bp_link(i,j,W,P1,q):
    #computes direct information
    [mu1,mu2]=compute_mu(i,j,W,P1,q)
    DI=compute_di(i,j,W,mu1,mu2,P1)
    #continue
    return DI

def compute_mu(i,j,W,P1,q):
    epsilon=float(.0001)
    diff=float(1.0)
    mu1=numpy.ones(q)/q
    mu2=numpy.ones(q)/q
    pi=P1[i-1,:]
    pj=P1[j-1,:]
    while diff>epsilon:
        
        scra1=numpy.dot(mu2,numpy.conjugate(W.T))
        scra2=numpy.dot(mu1,W)
        new1=numpy.divide(pi,scra1)
        new1=new1/numpy.sum(new1)
        
        new2=numpy.divide(pj,scra2)
        new2=new2/numpy.sum(new2)
        
        diff=numpy.max(numpy.maximum(abs(new1-mu1),abs(new2-mu2)),0)

        mu1=new1
        mu2=new2
        
    return [mu1,mu2]
        
def compute_di(i,j,W, mu1,mu2, Pia):
    #computes direct information
    tiny=1.0e-100
    
    Pdir=numpy.multiply(W,numpy.outer(mu1,mu2))
    Pdir=Pdir/sum(sum(Pdir))
    
    Pfac=numpy.outer(Pia[i-1,:],Pia[j-1,:])
    
    DI=numpy.trace(numpy.matmul(Pdir.T,numpy.log(numpy.divide((Pdir+tiny),(Pfac+tiny)))))
        
    return DI
    

    