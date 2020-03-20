#!/usr/bin/env python
# coding: utf-8

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


def mapa_TLE(ux=0.0,uy=0.0,uz=0.0):
    from numpy import array
    u=array([ux,uy,uz])
    Lambda=Lambda_TLE(-u)

    #Escoge valores de x:
    from numpy import linspace
    rmax=10
    xs=linspace(0,rmax,rmax+1,endpoint=True)
    ts=linspace(0,rmax,rmax+1,endpoint=True)

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
        ax.plot(tps,xps,'r-')

    for x in xs:
        for i,t in enumerate(ts):
            tps[i],xps[i],yp,zp=matmul(Lambda,[t,x,0,0])
        ax.plot(tps,xps,'r-')

    #Decoración
    ax.set_xticks(xs)
    ax.set_yticks(xs)
    ax.set_xlabel("$t$")
    ax.set_ylabel("$x$")
    ax.set_xlim((0,rmax))
    ax.set_ylim((0,rmax))
    ax.grid()
    fig.tight_layout()
    return fig


# ########################################
#  .//RelatividadEspecial.TransformacionesLorentzEinstein.Consecuencias.ipynb
# ########################################

# ########################################
#  .//RelatividadEspecial.ProblemasSeleccionados.ipynb
# ########################################

# ########################################
#  ./build/probs/RelatividadEspecial.Problemas.ipynb
# ########################################
