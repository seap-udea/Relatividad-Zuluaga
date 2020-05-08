#!/usr/bin/env python
# coding: utf-8

# ########################################
#  .//RelatividadGeneral.ipynb
# ########################################

# ########################################
#  .//RelatividadGeneral.PostuladosBasicos.ipynb
# ########################################

# ########################################
#  .//RelatividadGeneral.TransporteParalelo.ipynb
# ########################################

# ########################################
#  .//RelatividadGeneral.SimbolosChristoffel.ipynb
# ########################################

def Gamma(xmu,gfun,gargs=(),N=4):
    """
    Calcula todos los símbolos de Christoffel
    gfun: función métrica
    xmu: evento
    """
    from scipy.misc import derivative
    from numpy import where,arange
    from numpy import zeros
    #Indices
    index=arange(N)

    #Gamma
    G=zeros((N,N,N))
    
    for pi in range(N):
        for nu in range(N):
            #Inversa
            gpipi=1/gfun(xmu,pi,*gargs) #g^pipi
            #Coeficientes diagonales
            xd=xmu[pi] #Punto en el que estoy derivando
            dx=max(0.01,0.1*abs(xd))
            gnunu_pi=derivative(lambda x:gfun(where(index==pi,x,xmu),nu,*gargs),xd,dx)
            G[pi,nu,nu]=-0.5*gpipi*gnunu_pi
            #Coeficientes mixtos
            if nu==pi:continue
            xd=xmu[nu] #Punto en el que estoy derivando
            dx=max(0.01,0.1*abs(xd))
            gpipi_nu=derivative(lambda x:gfun(where(index==nu,x,xmu),pi,*gargs),xd,dx)
            G[pi,pi,nu]=0.5*gpipi*gpipi_nu
            G[pi,nu,pi]=G[pi,pi,nu]
    return G


from numpy import array
def g_cilindricas_4d(xmu,mu,p=1):
    """
    Coeficiente métrico g_mumu calculados en el evento xmu 
    para espacio-tiempo plano con coordenadas cilíndricas.
    
    g_munu=diag(1,-1,-r^2,-1)
    """
    from numpy import sin
    t,r,teta,z=xmu
    if mu==0:
        g=1
    elif mu==1:
        g=-1
    elif mu==2:
        g=-r**2
    elif mu==3:
        g=-1
    return g


# ########################################
#  .//RelatividadGeneral.Geodesicas.ipynb
# ########################################

def A_parallel(A,u,xfun,dxdufun,gfun,gargs=(),N=4):
    """
    Calcula la derivada de las componentes de un vector A 
    respecto al parámetro u de una función xfun
    
    Parametros:
        A: Arreglo con valores del vector
        u: Valor del parámetro
    
    Opciones:
        xfun: función que da la posición sobre la trayectoria
        dxdufun: función que da la derivadad de la trayectoria
        gfun: función que da la métrica
        N: Número de dimensiones
    """
    from export import Gamma
    from numpy import zeros
    dAdu=zeros(N)
    xmu=xfun(u)
    dxmudu=dxdufun(u)
    G=Gamma(xmu,gfun,gargs,N)
    for pi in range(N):
        for mu in range(N):
            for nu in range(N):
                dAdu[pi]+=-G[pi,mu,nu]*A[mu]*dxmudu[nu]
    return dAdu


from numpy import array
def g_cilindricas_2d(xmu,mu):
    """
    Coeficiente métrico g_mumu calculados en el evento xmu 
    para espacio-tiempo plano con coordenadas cilíndricas.
    
    g_munu=diag(1,r^2)
    """
    r,teta=xmu
    if mu==0:
        g=1
    elif mu==1:
        g=r**2
    return g


def ecuacion_geodesica(Y,s,gfun,gargs,N=4):
    """
    Opciones:
        gfun: función que da la métrica
        N: Número de dimensiones
    """
    from export import Gamma
    from numpy import zeros
    dYdu=zeros(2*N)
    x=Y[:N]
    dxds=Y[N:]

    dYdu[:N]=dxds
    G=Gamma(x,gfun,gargs,N)
    for pi in range(N):
        for mu in range(N):
            for nu in range(N):
                dYdu[N+pi]+=-G[pi,mu,nu]*dxds[mu]*dxds[nu]
    return dYdu

