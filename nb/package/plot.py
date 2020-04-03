def fija_ejes_proporcionales(ax,values=(),margin=0,xcm=None,ycm=None,xmin=None,ymin=None):
    """Ajusta los ejes para hacerlos proporcionales de acuerdo a un
    conjunto de valores.

    Normalmente esta tarea es realizada ax.set_aspect('equal','box')
    pero este comando solo se puede ejecutar después de que se han
    graficado los datos.  Esta rutina se puede ejecutar antes, si se
    pasan (como una tupla), todos los datos que van en el gráfico.

    Args:
      ax (matplotlib.axes): axes de matplotlib.

    Keyword Args:
      values (tuple): tupla de datos.
          Los datos deben corresponder a objetos de que puedan
          convertirse en arreglos de numpy. Si no se pasa nada se
          usan los valores en axes.

      margin (float): margen alrededor del gráfico.
          En unidades del (ancho o alto del mismo)

    Returns:
      (xlims,ylims) (tuple,tuple): Límites en x e y.

    """
    import numpy as np
    
    #values
    vals=np.array([])
    for value in values:
        vals=np.append(vals,np.array(value).flatten())
    #Center of values
    rcm=vals.mean()
    vals=vals-rcm

    if xcm is None:
        xcm=rcm
    if ycm is None:
        ycm=rcm
    
    fig=ax.figure
    bbox=ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width,height=bbox.width,bbox.height
    fx=width/height
    fy=1
    if fx<1:
        factor=fy
        fy=(1+margin)*1/fx
        fx=(1+margin)*factor
    else:
        fx*=(1+margin)
        fy*=(1+margin)

    max_value=np.abs(vals).max()
    ax.set_xlim((xcm-fx*max_value,xcm+fx*max_value))
    ax.set_ylim((ycm-fy*max_value,ycm+fy*max_value))

    if xmin is not None:
        xinf,xsup=ax.get_xlim()
        dx=xsup-xinf
        ax.set_xlim((xmin,xmin+dx))

    if ymin is not None:
        yinf,ysup=ax.get_ylim()
        dy=ysup-yinf
        ax.set_ylim((ymin,ymin+dy))

    return ax.get_xlim(),ax.get_ylim()

def fija_ejes3d_proporcionales(ax):
    """Ajusta los ejes en 3d para hacelos proporcionales.

    Hace que los ejes de un gráfico en 3d tengan la misma escala, de
    modo que las esferas aparezcan como esferas, los cubos como cubos
    y así sucesivamente.

    Esta es una de las soluciones alternativas para los comandos de
    matplotlib ax.set_aspect('equal') and ax.axis('equal') que no
    funcionan en 3D.

    Args:
      ax (matplotlib.axes): axis de matplotlib.
          Este debe ser el axis donde esta la figura.
    
    References: 
      tomado originalmente de https://stackoverflow.com/a/31364297

    """
    import numpy as np
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

    return ax.get_xlim3d(),ax.get_ylim3d(),ax.get_zlim3d()
