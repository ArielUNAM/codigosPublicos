"""Método de Newton-Rapshon
    Ariel Cerón G
    Computació  I
    Licenciatura en Matemáticas
    21/10/2019"""

def Newton(x0, fun, Dfun):
#def Newton(x0, fun, Dfun,error=0.01):
    """A la entrada se necesita un valor inicial, una función 
    y la derivada de la función, a la salida se obtiene la raíz
    de la función cercana al valor x0"""
    #e = error
    #efun = 10
    #x = x0
    #while(efun <= e):
    vfun = fun(x0)
    if type(Dfun) == int:
        vDfun = Dfun
    else:
        vDfun = Dfun(x0)
    x = x0 - (vfun/vDfun)
    return(x)

def funx(x):
    return(x)

def funx2(x):
    return(x**2)

def Dfunx2(x):
    return(2*x)