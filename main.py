# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 22:46:10 2021

@author: Axel
"""

from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import math
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


'''
=======
Hluti 1
=======
'''

def get_skilyrdi(i,L1,L2,h):
    #byrja að telja í 1
    #skilar skilyrði á punkti i
    N=L1/h
    M=L2/h

    dxv=[1+(N+1)*j for j in range(1,int(M))]
    dxh=[(N+1)*k for k in range(2,int(M+1))]
    w=[1*k for k in range(1,int(N+2))]
    v=[1+(N+1)*(M)+i for i in range(int(N+1))]
    
    if i in w:
        return 'w'
    if  i in v:
        return 'v'
    if i in dxv:
        return 'dxv'
    if i in dxh:
        return 'dxh'
    elif i<(N+1)*(M+1): return 'miðja'

def get_netpunktur(j,k,L1,h):
    #tekur inn (j,k)
    #skilar i
    N=L1/h
    return j+(k-1)*(N+1)
        

def get_jk(i,L1,h):
    #tekur inn i
    #skilar (j,k)
    N=L1/h
    if i%(N+1)==0:
        return ((N+1),i/(N+1))
    else:
        return (i%(N+1)),(1+i//(N+1))
    
def get_hnit(i,L1,h):
    #tekur inn i
    #skilar (x_j,y_k)
    j=get_jk(i,L1,h)[0]
    k=get_jk(i,L1,h)[1]
    return (((j-1)*h,(k-1)*h))


    
    

    
#w =lambda x: -10*x/1*(x-1)**2*(1+x)
#v =lambda x:  x*(1-x)*(1+x)**2    
    
w =lambda x: 1
v =lambda x: 0   

#def non_zero(l,m):
#    return 2*l+4*2*(m-2)+(m*l-2*l-2*(m-2))*5

def construct_matrix(L1,L2,h,q,v,w):
     N=L1/h
     M=L2/h
     l=int(N+1)
     m=int(M+1)
     
     A = sp.lil_matrix((l*m,l*m))
     b = np.zeros((l*m,1))
     
     h2=1/(h**2)
     
     for i in range(1,int(l*m)+1):
         if get_skilyrdi(i,L1,L2,h)=='miðja':
            A[i-1,i-1]=4*h2-q**2
            A[i-1,i-2]=-h2
            A[i-1,i]=-h2
            A[i-1,i-1-l]=-h2
            A[i-1,i-1+l]=-h2
            
         if get_skilyrdi(i,L1,L2,h)=='dxh':
            A[i-1,i-1]=4*h2-q**2
            A[i-1,i-2]=-2*h2
            A[i-1,i-1+l]=-h2
            A[i-1,i-1-l]=-h2
            
         if get_skilyrdi(i,L1,L2,h)=='dxv':
            A[i-1,i-1]=4*h2-q**2
            A[i-1,i]=-2*h2
            A[i-1,i-1+l]=-h2
            A[i-1,i-1-l]=-h2
            
         if get_skilyrdi(i,L1,L2,h)=='v':
            A[i-1,i-1]=1
            b[i-1]= v(get_hnit(i,L1,h)[0])
         
         if get_skilyrdi(i,L1,L2,h)=='w':
             A[i-1,i-1]=1
             b[i-1]=w(get_hnit(i,L1,h)[0])
         
     A = A.tocsr()
     return A,b




def helmholtzeq(L1,L2,h,q,v,w):
     N=round(L1/h)
     M=round(L2/h)
     l=N+1
     m=M+1
     
     A,b=construct_matrix(L1,L2,h,q,v,w)
    
     c = spsolve(A,b).tolist()
    
     HZ=np.ones((m,l))
     k=0
     for i in range(m):
         for j in range(l):
             HZ[i][j]=c[k]
             k=k+1
     return HZ
             
#fræðileg lausn
def func(h,q,L1,L2):
    l = round(L1 / h) + 1
    m = round(L2 / h) + 1
    x = np.linspace(0, L1, l)
    y = np.linspace(0, L2, m)
    X,Y=np.meshgrid(x,y)
    return np.sin(q*(L2-Y))/(np.sin(q*L2))
    
'''
skilgreinið fyrst fylki til að plotta t.d. V=Varmajafnvaegi(2,1,1/50,1,2) 
og passið því inn í plot_grid, það þarf að hafa sömu L1 og L2 og h.
'''

def plot_grid(L1,L2,h,Z):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    fig1, ax1 = plt.subplots()
    l = round(L1 / h) + 1
    m = round(L2 / h) + 1
    x = np.linspace(0, L1, l)
    y = np.linspace(0, L2, m)
    X,Y=np.meshgrid(x,y)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    ax.set_title('')
    
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    C=ax1.contour(X,Y, Z)
    ax1.clabel(C, inline=True, fontsize=10)
    ax1.set_title('')
    plt.show()
    
    
'''
=======
Hluti 2
=======
'''

#fourier
def b(n):
    return 2/((2*n+1)*np.pi)
#fourier
def S(x):
    s=0
    for i in range(0,150):
        s=s+b(i)*np.sin((2*i+1)*x)
    return 4*s

#psi1 = lambda x,b,a: b*(np.sin(2*np.pi/a*(x-1/2))+1)
#psi2 = lambda x,b : b*(np.cos(np.pi/4*(x-1))-1/(np.sqrt(2)))

psi1 = lambda x,b,a: np.cos(5*x+np.pi/2)**3*(-4)-1
#psi2 = lambda x,b : S(x*8)

psi2= lambda x,b: np.sin(x*3-np.pi/2)*4-1

def construct_matrix2(L1,L2,h,beta1,beta2):
    N=L1/h
    M=L2/h
    l=int(N+1)
    m=int(M+1)
    
    A = sp.lil_matrix((l*m,l*m))
    b = np.zeros((l*m,1))
    for i in range(1,int(l*m)+1):
        j=i-1
        if get_skilyrdi(i,L1,L2,h)=='miðja':
            A[j,j]=4
            A[j,j-1]=-1
            A[j,j+1]=-1
            A[j,j-l]=-1
            A[j,j+l]=-1
        if get_skilyrdi(i,L1,L2,h)=='v':
            A[j,j]=1
            b[j]=psi2(get_hnit(i,L1,h)[0],beta2)
        if get_skilyrdi(i,L1,L2,h)=='w':
            A[j,j]=1
            b[j]=psi1(get_hnit(i,L1,h)[0],beta1,L1)
            
        if get_skilyrdi(i,L1,L2,h)=='dxh':
            A[j,j]=2
            A[j,j-1]=-1
            A[j,j+l]=-1/2
            A[j,j-l]=-1/2
        if get_skilyrdi(i,L1,L2,h)=='dxv':
            A[j,j]=2
            A[j,j+1]=-1
            A[j,j-l]=-1/2
            A[j,j+l]=-1/2
    A=A.tocsr()
    return A,b


def Varmajafnvaegi(L1,L2,h,beta1,beta2):
    N=round(L1/h)
    M=round(L2/h)
    l=N+1
    m=M+1
    
    A,b=construct_matrix2(L1,L2,h,beta1,beta2)
    
    c = spsolve(A,b).tolist()
    V=np.ones((m,l))
    k=0
    
    for i in range(m):
        for j in range(l):
            V[i][j]=c[k]
            k=k+1
    return V