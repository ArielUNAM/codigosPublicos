###############################
import numpy as np
import matplotlib.pyplot as plt
from sympy import *

##############Definiciones y funciones
def tempH(x:object, n:int):
    if n == 0:
        return(1)
    elif n == 1:
        return(2*x)
    else:
        return(2*x*tempH(x,n-1)-2*(n-1)*tempH(x,n-2))
    
def IHermite(xl,yl,sl,xr,yr,sr,X:list)->list:
    """
        Devuelve una lista de coeficientes
    """
    x = Symbol('x')
    dif = xr - xl
    yp = (yr - yl)/dif
    c = (yp - sl)/dif
    d = (sl + sr - 2*yp )/(dif**2)

    p = yl + sl*(x - X[0]) + c*(x-X[0])**2 + d*(x-X[-1])*(x-X[0])**2
    p = Poly(p,x)
    print(expand(p))
    return([xl,sl,c,d])

def IHermite2(B:list,X:list)->list:
    """
        Devuelve una lista de coeficientes
    """
    x = Symbol('x')
    p = B[0] + B[1]*(x - X[0]) + B[2]*(x-X[0])**2 + B[3]*(x-X[-1])*(x-X[0])**2
    p = Poly(p,x)
    print(expand(p))
    return(p.all_coeffs())

def INewton(X:list,Y:object)->list:
    """
        Devuelve una lista de coeficientes

        Devuelve una lista de coeficientes con el método de diferencias divididas de Newton  que siguen la ecuación

            P_N(x) = sum^N f[x_1,...,x_i] mul^N (x - x_j)

        Entradas:
        ---------
        X: list
            Lista de puntos
        Y: list
            Lista de imágenes de puntos

        Salidas
        a: list
            Lista de coeficientes del polinómio
    """
    B = DDivididas(X,Y)
    x = Symbol('x')
    p = 0
    for i in range(len(X)):
        if i == 0:
            p += B[i]
        else:
            aux = 1
            for j in range(i):
                aux *= (x-X[j])
            p += B[i]*aux
    a = Poly(p, x) 
    print("Polinomio: ",expand(p))
    return(a.all_coeffs())

def DDivididas(X:list,Y:object)->list:
    """
        Devuelve la diagonal de una matriz de diferencias divididas
    """
    l = []
    for i in range(len(X)):
        for j in range(len(X)-i):
            if j == 0:
                if type(Y) == list:
                    print("i:{},j:{}".format(i,j))
                    l.append(DDaux(X[:i+1],i+1,Y[:i+1]))
                else:
                    l.append(DDaux(X[:i+1],i+1,Y))
    return(l)
            
def DDaux(X:list,n:int,Y:object):
    """
        Devuelve el valor de la diferencia dividida
    """
    if type(Y) == list:
        if n == 1:
            return(Y[0])
        elif n == 2:
            #print("1")
            #print(X)
            #print(Y)
            num = Y[-1] - Y[0]
            den = X[-1] - X[0]
            #print("1) {} - {}/({} - {})".format(Y[-1], Y[0],X[-1], X[0]))
            return(num/den)
        else:
            #print("2")
            #print(X)
            #print(Y)
            izq = DDaux(X[1:n],n-1,Y[1:n])
            der = DDaux(X[0:n-1],n-1,Y[0:n-1])
            #print("2) {} - {}/({} - {})".format(izq, der,X[-1], X[0]))
            num = izq-der
            den = X[-1] - X[0]
            return(num/den)
    else:
        if n == 1:
            return(Y(X[0]))
        elif n == 2:
            num = Y(X[-1]) - Y(X[0])
            den = X[-1] - X[0]
            #print("1) {}/{}".format(num,den))
            return(num/den)
        else:
            izq = DDaux(X[1:n],n-1,f)
            der = DDaux(X[0:n-1],n-1,f)
            num = izq-der
            den = X[-1] - X[0]
            #print("2) {}/{}".format(num,den))
            return(num/den)

def ILegendre(X:list, Y:list)->list:
    """
        Devuelve un vector de coeficientes 

        Devuelve un vector de coeficientes relacionados a un polinómio de grado N, obtenido con la formula
            P_N(x) = sum_N Y[i]*l[i]
        Con l[i] el polinómio de Legendre definido
            l[i] = mul_N (x - X[i])/(X[j] - X[i])
        con i != j

        Entradas:
        ---------
        X: list
            Lista de puntos
        Y: list
            Lista de valores de la función en los puntos X
        
        Salidas:
        --------
        B:list
            Lista de coeficientes para el polinómio de orden n
                P_N = b_0x**N + b_1x^{N-1} + ... + b_N
    """
    #Definimos los coeficientes de la combinación lineal
    tmp = 0
    x = Symbol('x')
    for i in range(len(X)):
        aux = 1
        p = 1
        for j in range(len(X)):
            if i != j:
                aux *= (X[i] - X[j])
                p *= (x - X[j])
        tmp += (Y[i]/aux)*p
    a = Poly(tmp, x) 
    print("Polinomio: ",expand(tmp))
    return(a.all_coeffs())

def MHorner(B:list, x0:float) -> float:
    """
    Devuelve el valor resultante al evaluar un polinómio p(x) con un valor x_o

    El método de Horner busca evaluar de forma iterativa los elementos del polinómio de forma monómica haciendo una simple sustitución de valores, en general el algoritmo sigue la siguiente forma
        r0 = a_0
        r1 = a_1 + r0*x0

    Entradas:
    --------
    a: list
        Coeficientes del polinómio
    x0: float
        Punto a evaluar del polinómio

    Salidas:
    --------
    r: floar
        Valor del polinómio evaluado
    
    """
    r = B[0]
    for i in range(1,len(B)):
        r = B[i] + x0*r
    return(r)

######### Evalución

X = [-0.1,0.8,1.4]
Y = [-0.198669,0.999574,0.334988]

B = [4.64373783264746, -13.2354993470508, 9.77209654211249, -0.620204453772291, 0.527106346677640, -0.128614236939369]

X2 = np.linspace(-1,2,20)
Y2 = [MHorner(B,x) for x in np.linspace(-1,2,20)]

plt.figure(figsize=(10,10))
plt.plot(X,Y,'o',label="Dispersión")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Diagrama de dispersión")
plt.savefig("./F6.png")
#plt.show()

plt.figure(figsize=(10,10))
plt.plot(X2,Y2,'o-',label="Polinómio obtenido")
plt.plot(X,Y,'o',label="Valores obtenidos")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Polinómio de Newton")
plt.savefig("./F7.png")
#plt.show()


f = lambda x: np.sin(2*x)

X3 = [-0.2,0.2,0.4,1.0,1.6]
Y3 = [f(x) for x in X3]
Y4 = [MHorner(B,x) for x in X3]

e = list()
for i in range(len(Y3)):
    e.append(abs(Y3[i] - Y4[i]))

er = list()
for j in range(len(e)):
    er.append(abs(e[i]/Y3[i]))

plt.figure(figsize=(10,10))
plt.plot(X3,Y4,'o-',label="Polinómio obtenido")
plt.plot(X3,Y3,'o-',label="Función")
plt.plot([X3[0],X3[0]],[Y4[0],Y3[0]],'o-',color="red",label="Error")
for i in range(1,len(X3)):
    plt.plot([X3[i],X3[i]],[Y4[i],Y3[i]],'o-',color="red")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Error grafico")
plt.savefig("./F8.png")
#plt.show()

plt.figure(figsize=(10,10))
plt.plot(e,'o-')
plt.grid()
plt.xlabel("x")
plt.ylabel("e(x)")
plt.title("Error absoluto")
plt.savefig("./F9.png")
#plt.show()

plt.figure(figsize=(10,10))
plt.plot(er,'o-')
plt.grid()
plt.xlabel("x")
plt.ylabel("e(x)")
plt.title("Error relativo")
plt.savefig("./F10.png")
#plt.show()



"""
>>> x0 = -0.1
>>> x1 = 0.8
>>> x2 = 1.4
>>> y0 = -0.198669
>>> y1 = 0.999574
>>> y2 = 0.334988
>>> s0 = 0.995004
>>> s1 = 0.696707
>>> s2 = 0.169967
>>> dd1 = 1.327121
>>> dd2 = -1.101253
>>> dd3 = -1.107643
>>> x0x1 = dd1
>>> x1x2 = dd2
>>> x1x2 = dd3
>>> (x0x1 - s0)/0.9
0.3690188888888889
>>> x0x0x1 = (x0x1 - s0)/0.9
>>> x0x1x1 = (s1 - x0x1)/0.9
>>> x0x1x1
-0.70046
>>> x1x1x2 = (x1x2 - x1x1)/(x2 - x1)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'x1x1' is not defined
>>> x1x1x2 = (x1x2 - s1)/(x2 - x1)
>>> x1x1x2
-3.0072500000000004
>>> x1x2x2 = (s2 - x1x2)/(x2 - x1)
>>> x1x2x2
2.1293500000000005
>>> x0x0x1x1 = (x0x1x1 - x0x0x1)/(x1 - x0)
>>> x0x0x1x1
-1.1883098765432099
>>> x0x1x1x2 = (x1x1x2 - x0x1x1)/(x2 - x0)
>>> x0x1x1x2
-1.5378600000000002
>>> x1x1x2x2 = (x1x2x2 - x1x1x2)/(x2 - x1)
>>> x1x1x2x2
8.561000000000003
>>> x0x0x1x1x2 = (x0x1x1x2 - x0x0x1x1)/(x2 - x0)
>>> x0x0x1x1x2
-0.23303341563786026
>>> x0x1x1x2x2 = (x1x1x2x2 - x0x1x1x2)/(x2 - x0)
>>> x0x1x1x2x2
6.7325733333333355
>>> x0x0x1x1x2x2 = (x0x1x1x2x2 - x0x0x1x1x2)/(x2 - x0)
>>> x0x0x1x1x2x2
4.643737832647464
>>> from sympy import *
>>> x = Symbol('x')
>>> B = [y0,y1, x0x0x1,x0x0x1x1, x0x0x1x1x2,x0x0x1x1x2x2]
>>> p = B[0] + B[1]*(x-x0) + B[2]*(x-x0)*(x-x0) + B[3]*(x-x0)*(x-x0)*(x-x1) + B[4]*(x-x0)*(x-x0)*(x-x1)*(x-x1) + B[5]*(x-x0)(x-x0)*(x-x1)*(x-x1)*(x-x2) 
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: 'Add' object is not callable
>>> p = B[0] + B[1]*(x-x0) + B[2]*(x-x0)*(x-x0) 
>>> p = B[0] + B[1]*(x-x0) + B[2]*(x-x0)*(x-x0) + B[3]*(x-x0)*(x-x0)*(x-x1) 
>>> p = B[0] + B[1]*(x-x0) + B[2]*(x-x0)*(x-x0) + B[3]*(x-x0)*(x-x0)*(x-x1) + B[4]*(x-x0)*(x-x0)*(x-x1)*(x-x1) 
>>> p = B[0] + B[1]*(x-x0) + B[2]*(x-x0)*(x-x0) + B[3]*(x-x0)*(x-x0)*(x-x1) + B[4]*(x-x0)*(x-x0)*(x-x1)*(x-x1) + B[5]*(x-x0)*(x-x0)*(x-x1)*(x-x1)*(x-x2) 
>>> p
0.999574*x + (-1.18830987654321*x - 0.118830987654321)*(x - 0.8)*(x + 0.1) + (-0.23303341563786*x - 0.023303341563786)*(x - 0.8)**2*(x + 0.1) + (0.369018888888889*x + 0.0369018888888889)*(x + 0.1) + (x - 1.4)*(x - 0.8)**2*(x + 0.1)*(4.64373783264746*x + 0.464373783264746) - 0.0987116
>>> a = Poly(p,x)
>>> C = a.all_coeffs()
>>> C
[4.64373783264746, -13.2354993470508, 9.77209654211249, -0.620204453772291, 0.527106346677640, -0.128614236939369]
>>> 
KeyboardInterrupt
>>> B
[-0.198669, 0.999574, 0.3690188888888889, -1.1883098765432099, -0.23303341563786026, 4.643737832647464]
>>> p
0.999574*x + (-1.18830987654321*x - 0.118830987654321)*(x - 0.8)*(x + 0.1) + (-0.23303341563786*x - 0.023303341563786)*(x - 0.8)**2*(x + 0.1) + (0.369018888888889*x + 0.0369018888888889)*(x + 0.1) + (x - 1.4)*(x - 0.8)**2*(x + 0.1)*(4.64373783264746*x + 0.464373783264746) - 0.0987116
>>> expand(p)


X = [-0.1,0.8,1.4]
Y = [-0.198669,0.999574,0.334988]

print(DDivididas(X,Yf))
B = INewton(X,Y)

X2 = np.linspace(-1,2,20)
Y2 = [MHorner(B,x) for x in np.linspace(-1,2,20)]

plt.figure(figsize=(10,10))
plt.plot(X,Y,'o',label="Dispersión")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Diagrama de dispersión")
plt.savefig("./F6.png")
#plt.show()

plt.figure(figsize=(10,10))
plt.plot(X2,Y2,'o-',label="Polinómio obtenido")
plt.plot(X,Y,'o',label="Valores obtenidos")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Polinómio de Newton")
plt.savefig("./F7.png")
#plt.show()


f = lambda x: np.sin(2*x)

X3 = [-0.2,0.2,0.4,1.0,1.6]
Y3 = [f(x) for x in X3]
Y4 = [MHorner(B,x) for x in X3]

e = list()
for i in range(len(Y3)):
    e.append(abs(Y3[i] - Y4[i]))

er = list()
for j in range(len(e)):
    er.append(abs(e[i]/Y3[i]))

plt.figure(figsize=(10,10))
plt.plot(X3,Y4,'o-',label="Polinómio obtenido")
plt.plot(X3,Y3,'o-',label="Función")
plt.plot([X3[0],X3[0]],[Y4[0],Y3[0]],'o-',color="red",label="Error")
for i in range(1,len(X3)):
    plt.plot([X3[i],X3[i]],[Y4[i],Y3[i]],'o-',color="red")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Error grafico")
plt.savefig("./F8.png")
#plt.show()

plt.figure(figsize=(10,10))
plt.plot(e,'o-')
plt.grid()
plt.xlabel("x")
plt.ylabel("e(x)")
plt.title("Error absoluto")
plt.savefig("./F9.png")
#plt.show()

plt.figure(figsize=(10,10))
plt.plot(er,'o-')
plt.grid()
plt.xlabel("x")
plt.ylabel("e(x)")
plt.title("Error relativo")
plt.savefig("./F10.png")
#plt.show()




X = [1.1,1.6]
Y = [9.0250,24.5325]

B = INewton(X,Y)

X2 = np.linspace(0,2,20)
Y2 = [MHorner(B,x) for x in np.linspace(0,2,20)]

plt.figure(figsize=(10,10))
plt.plot(X,Y,'o',label="Dispersión")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Diagrama de dispersión")
plt.savefig("./F1.png")
#plt.show()

plt.figure(figsize=(10,10))
plt.plot(X2,Y2,'o-',label="Polinómio obtenido")
plt.plot(X,Y,'o',label="Valores obtenidos")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Polinómio de Newton")
plt.savefig("./F2.png")
#plt.show()


f = lambda x: np.exp(2*x)

X3 = [1, 1.2, 1.4, 1.6, 1.8]
Y3 = [f(x) for x in X3]
Y4 = [MHorner(B,x) for x in X3]
e = list()
for i in range(len(Y3)):
    e.append(abs(Y3[i] - Y4[i]))

er = list()
for j in range(len(e)):
    er.append(e[i]/Y3[i])

plt.figure(figsize=(10,10))
plt.plot(X3,Y4,'o-',label="Polinómio obtenido")
plt.plot(X3,Y3,'o-',label="Función")
plt.plot([X3[0],X3[0]],[Y4[0],Y3[0]],'o-',color="red",label="Error")
for i in range(1,len(X3)):
    plt.plot([X3[i],X3[i]],[Y4[i],Y3[i]],'o-',color="red")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Error grafico")
plt.savefig("./F3.png")
#plt.show()

plt.figure(figsize=(10,10))
plt.plot(e,'o-')
plt.grid()
plt.xlabel("x")
plt.ylabel("e(x)")
plt.title("Error absoluto")
plt.savefig("./F4.png")
#plt.show()

plt.figure(figsize=(10,10))
plt.plot(er,'o-')
plt.grid()
plt.xlabel("x")
plt.ylabel("e(x)")
plt.title("Error relativo")
plt.savefig("./F5.png")
#plt.show()

X = [-0.1,0.8,1.4]
Y = [-0.198669,0.999574,0.334988]

B = ILegendre(X,Y)

X2 = np.linspace(-1,2,20)
Y2 = [MHorner(B,x) for x in np.linspace(-1,2,20)]

plt.figure(figsize=(10,10))
plt.plot(X,Y,'o',label="Dispersión")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.ylabel("Y(x)")
plt.title("Diagrama de dispersión")
plt.savefig("./F6.png")
#plt.show()

plt.figure(figsize=(10,10))
plt.plot(X2,Y2,'o-',label="Polinómio obtenido")
plt.plot(X,Y,'o',label="Valores obtenidos")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Polinómio de Legendre")
plt.savefig("./F7.png")
#plt.show()


f = lambda x: np.sin(2*x)

X3 = [-0.2,0.2,0.4,1.0,1.6]
Y3 = [f(x) for x in X3]
Y4 = [MHorner(B,x) for x in X3]

e = list()
for i in range(len(Y3)):
    e.append(abs(Y3[i] - Y4[i]))

er = list()
for j in range(len(e)):
    er.append(abs(e[i]/Y3[i]))

plt.figure(figsize=(10,10))
plt.plot(X3,Y4,'o-',label="Polinómio obtenido")
plt.plot(X3,Y3,'o-',label="Función")
plt.plot([X3[0],X3[0]],[Y4[0],Y3[0]],'o-',color="red",label="Error")
for i in range(1,len(X3)):
    plt.plot([X3[i],X3[i]],[Y4[i],Y3[i]],'o-',color="red")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Error grafico")
plt.savefig("./F8.png")
#plt.show()

plt.figure(figsize=(10,10))
plt.plot(e,'o-')
plt.grid()
plt.xlabel("x")
plt.ylabel("e(x)")
plt.title("Error absoluto")
plt.savefig("./F9.png")
#plt.show()

plt.figure(figsize=(10,10))
plt.plot(er,'o-')
plt.grid()
plt.xlabel("x")
plt.ylabel("e(x)")
plt.title("Error relativo")
plt.savefig("./F10.png")
#plt.show()




X = [1.1,1.6]
Y = [9.0250,24.5325]

B = ILegendre(X,Y)

X2 = np.linspace(0,2,20)
Y2 = [MHorner(B,x) for x in np.linspace(0,2,20)]

plt.figure(figsize=(10,10))
plt.plot(X,Y,'o',label="Dispersión")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Diagrama de dispersión")
plt.savefig("./F1.png")
#plt.show()

plt.figure(figsize=(10,10))
plt.plot(X2,Y2,'o-',label="Polinómio obtenido")
plt.plot(X,Y,'o',label="Valores obtenidos")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Polinómio de Legendre")
plt.savefig("./F2.png")
#plt.show()


f = lambda x: np.exp(2*x)

X3 = [1, 1.2, 1.4, 1.6, 1.8]
Y3 = [f(x) for x in X3]
Y4 = [MHorner(B,x) for x in X3]
e = list()
for i in range(len(Y3)):
    e.append(abs(Y3[i] - Y4[i]))

er = list()
for j in range(len(e)):
    er.append(e[i]/Y3[i])

plt.figure(figsize=(10,10))
plt.plot(X3,Y4,'o-',label="Polinómio obtenido")
plt.plot(X3,Y3,'o-',label="Función")
plt.plot([X3[0],X3[0]],[Y4[0],Y3[0]],'o-',color="red",label="Error")
for i in range(1,len(X3)):
    plt.plot([X3[i],X3[i]],[Y4[i],Y3[i]],'o-',color="red")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Error grafico")
plt.savefig("./F3.png")
#plt.show()

plt.figure(figsize=(10,10))
plt.plot(e,'o-')
plt.grid()
plt.xlabel("x")
plt.ylabel("e(x)")
plt.title("Error absoluto")
plt.savefig("./F4.png")
#plt.show()

plt.figure(figsize=(10,10))
plt.plot(er,'o-')
plt.grid()
plt.xlabel("x")
plt.ylabel("e(x)")
plt.title("Error relativo")
plt.savefig("./F5.png")
#plt.show()
"""