#!/usr/bin/env python
# coding: utf-8

# ########################################
#  .//Prefacio.ipynb
# ########################################

# ########################################
#  .//Agradecimientos.ipynb
# ########################################

# ########################################
#  .//Introduccion.ipynb
# ########################################

def calcula_discriminante(a,b,c):
    disc=b**2-4*a*c
    return disc


# ########################################
#  .//RelatividadEspecial.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.ConceptosBasicos.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.TransformacionesLorentzEinstein.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.TransformacionesLorentzEinstein.Propiedades.ipynb
# ########################################

def Lambda_TLE(u):
    from numpy import zeros
    Lambda=zeros((4,4))
    
    #Factor de Lorentz
    umag=(u[0]**2+u[1]**2+u[2]**2)**0.5
    gamma=(1-umag**2)**(-0.5)
    
    #Lambda
    Lambda[0,0]=gamma
    Lambda[0,1:]=-u*gamma
    Lambda[1:,0]=-u*gamma
    for i in range(1,4):
        for j in range(1,4):
            dij=0
            if i==j:dij=1
            Lambda[i,j]=dij+(gamma-1)*u[i-1]*u[j-1]/umag**2
    return Lambda


def mapa_TLE(ux=0.0,uy=0.0,uz=0.0,
             rmax=10,ngrid=10,nticks=10,
             interact=False):
    from numpy import array
    u=array([ux,uy,uz])
    Lambda=Lambda_TLE(-u)

    #Escoge valores de x:
    from numpy import linspace
    xs=linspace(0,rmax,ngrid+1,endpoint=True)
    ts=linspace(0,rmax,ngrid+1,endpoint=True)

    #Calcula valores de t' y x' usando la matriz:
    from numpy import zeros_like
    tps=zeros_like(xs)
    xps=zeros_like(xs)

    from numpy import matmul

    import matplotlib.pyplot as plt
    fig=plt.figure(figsize=(5,5))
    ax=fig.gca()

    for t in xs:
        for i,x in enumerate(xs):
            tps[i],xps[i],yp,zp=matmul(Lambda,[t,x,0,0])
        ax.plot(tps,xps,'r-',alpha=0.5)

    for x in xs:
        for i,t in enumerate(ts):
            tps[i],xps[i],yp,zp=matmul(Lambda,[t,x,0,0])
        ax.plot(tps,xps,'r-',alpha=0.5)

    #Decoración
    ax.set_xticks(linspace(0,rmax,nticks+1,endpoint=True))
    ax.set_yticks(linspace(0,rmax,nticks+1,endpoint=True))
    ax.set_xlabel("$t$")
    ax.set_ylabel("$x$")
    ax.set_xlim((0,rmax))
    ax.set_ylim((0,rmax))
    ax.grid()
    fig.tight_layout()
    if not interact:
        return fig


# ########################################
#  .//RelatividadEspecial.TransformacionesLorentzEinstein.Consecuencias.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.Minkowski.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.OpticaRelativista.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.CinematicaRelativista.Definiciones.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.CinematicaRelativista.AceleracionPropiaConstante.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.DinamicaRelativista.EnergiaMomentum.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.DinamicaRelativista.ColisionesRelativistas.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.DinamicaRelativista.Cuadrifuerza.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.DinamicaRelativista.FuerzaDeLorentz.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.SintesisMecanicaRelativista.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.Electrodinamica.EcuacionesClasicas.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.Electrodinamica.CuadriCorrientePotencial.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.Electrodinamica.TensorFaraday.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.Electrodinamica.Ejemplos.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.ProblemasSeleccionados.ipynb
# ########################################

# ########################################
#  ./build/probs/RelatividadEspecial.Problemas.ipynb
# ########################################

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
def g_cilindricas_4d(xmu,mu):
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


# ########################################
#  .//RelatividadGeneral.InerciaYGeodesicas.ipynb
# ########################################

from numpy import array
def g_newtoniana_4d(xmu,mu,R=1):
    """
    Coeficiente métrico g_mumu calculados en el evento xmu 
    para espacio-tiempo plano con coordenadas cilíndricas.
    
    g_munu=diag(A,-1,-r^2,-r^2 sin^2 teta)
    """
    from numpy import sin
    t,r,fi,teta=xmu
    A=(1-R/r)
    if mu==0:
        g=A
    elif mu==1:
        g=-1
    elif mu==2:
        g=-r**2
    elif mu==3:
        g=-r**2*sin(teta)**2
    return g

