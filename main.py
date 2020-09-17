import numpy as np
import matplotlib.pyplot as plt

def TGauss(a:float,b:float,x:float):
    num = (b-a)*x + (b+a)
    den =  2
    return(num/den)

t = [1.1, 1.2, 1.3, 1.4]
I = [0.9849, 2.2368, 3.7011, 5.3154]
It = [8.71, 8.31, 7.95, 7.65]
Is = [8.71, 8.31, 7.95, 7.65]

plt.figure(figsize=(10,10))
plt.plot(t,I,'o-',label="Valores del ejercicio")
plt.plot(t,It,'o-',label="Valores por el método del trapecio")
plt.plot(t,Is,'o-',label="Valores por el método de Simpson")
plt.plot([t[0],t[0],t[0]],[I[0],It[0],Is[0]],'o-',color="red",label="Error")
for i in range(1,len(t)):
    plt.plot([t[i],t[i],t[i]],[I[i],It[i],Is[i]],'o-',color="red")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Error grafico")
plt.savefig("/home/aceron/Onedrive/UnAD/Quinto Semestre/Análisis Numérico II/Unidad III/EA/f4.png")
#plt.show()
"""
#w = [5/9,-5/9,8/9]
#x = [np.sqrt(3/5), -np.sqrt(3/5), 0]
w = [1,1]
x = [np.sqrt(3)/3,-np.sqrt(3)/3]
f = lambda x: 10*np.exp((-x/10)*(1/10*np.sin(2*x) - 2*np.cos(2*x)))

a = 1
b = 1.4

S = [((b - a)/2)*w[i]*f(TGauss(a,b,x[i])) for i in range(len(x))]
print(S)
print(sum(S))

h1 = 0.1
aux = 1
x1 = []
while aux < 1.5:
    x1.append(aux)
    aux += 0.1

f = lambda x: 10*np.exp((-x/10)*(1/10*np.sin(2*x) - 2*np.cos(2*x)))

aux = 0
print("\\begin{table}[H]")
print("\t\\centering")
print("\t\\begin{tabular}{|c|c|c|}")
print("\t\\hline")
print("i & $x_i$ & f_i  \\\\ \\hline")
for i in range(len(x1)-1):
    print("{:.2f} & {:.2f} & {:.2f}\\\\ \\hline".format(x1[i],f(x1[i]),h1/2*(f(x1[i]) + f(x1[i+1]))))
    aux += h1/2*(f(x1[i]) + f(x1[i+1]))
    #print("{:.2f} &{:.2f} & {:.2f} \\\\ \\hline".format(i, x1[i],f(x1[i])))

print("{:.2f} & {:.2f} & {:.2f}\\\\ \\hline".format(x1[-1],f(x1[-1]),h1/2*(f(x1[-2]) + f(x1[-1]))))
aux += h1/2*(f(x1[i]) + f(x1[i+1]))
print("\t\\end{tabular}")
print("\\end{table}")

print(aux)

w1 = 5/9
w2 = 8/9

x1 = np.sqrt(3/5)
x2 = 0

f = lambda x: 4*(np.sqrt(1 - x*x))

print(w1*f(x1) + w1*f(-x1) + w2*f(x2))



h1 = 0.1
s = -1.1
x1 = []
while s <= 1:
    x1.append(s)
    s += 0.1

#f = lambda x: 1 + np.exp(-x)
f = lambda x: 4*(np.sqrt(1 - x*x))
aux = 0

fa = f(x1[1])
fb = f(x1[-1])
P = 0
I = 0
for i in range(1,len(x1)-1):
    print("x{}",i)
    print("f(x) = ",f(x1[i]))
    if i%2 == 0:
        P += f(x1[i])
    else:
        I += f(x1[i])

print("\\begin{table}[H]")
print("\t\\centering")
print("\t\\begin{tabular}{|c|c|c|}")
print("\t\\hline")
print("i & $x_i$ & f_i  \\\\ \\hline")
for i in range(1,len(x1)):
    #print("{} & {} & {}\\\\ \\hline".format(x1[i],f(x1[i]),h1/2*(f(x1[i]) + f(x1[i+1]))))
    print("{:.2f} &{:.2f} & {:.2f} \\\\ \\hline".format(i, x1[i],f(x1[i])))
print("\t\\end{tabular}")
print("\\end{table}")

print(fa)
print(fb)
print(4*I)
print(2*P)
print(h1/3 * (fa+fb+4*I+2*P))


h1 = 0.1
s = -1.1
x1 = []
while s <= 1:
    x1.append(s)
    s += 0.1

f = lambda x: 4*(np.sqrt(1 - x*x))

aux = 0

print("\\begin{table}[H]")
print("\t\\centering")
print("\t\\begin{tabular}{|c|c|c|}")
print("\t\\hline")
print("$x_i$ & $f_i$ & $h/2(f_i + f_{i+1} $  \\\\ \\hline")
for i in range(1,len(x1)-1):
    print("{:.2f} & {:.2f} & {:.2f}\\\\ \\hline".format(x1[i],f(x1[i]),h1/2*(f(x1[i]) + f(x1[i+1]))))
    aux += h1/2*(f(x1[i]) + f(x1[i+1]))
    
#    print("{:.2f} &{:.2f} & {:.2f} \\\\ \\hline".format(i, x1[i],f(x1[i])))
print("{:.2f} & {:.2f} & {:.2f}\\\\ \\hline".format(x1[-1],f(x1[-1]),h1/2*(f(x1[-2]) + f(x1[-1]))))
aux += h1/2*(f(x1[i]) + f(x1[i+1]))
print("\t\\end{tabular}")
print("\\end{table}")

print(aux)

x = np.linspace(-1.5,1.5,20)
f = lambda x: np.cos(x)
y = [f(i) for i in x]

plt.figure(figsize=(10,10))
plt.plot(x,y,'o-')
plt.title("Gráfica de la función")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.grid()
plt.savefig("./f1.png")

#h1 = 1
#x1 = [0,1,2]

h1 = 0.5
x1 = [0,0.5,1,1.5,2,2.5,3,3.5]

#f = lambda x: 1 + np.exp(-x)
f = lambda x: 1 + np.exp(x)*np.cos(4*x)
aux = 0

fa = f(x1[0])
fb = f(x1[-1])
P = 0
I = 0
for i in range(1,len(x1)-1):
    print("x{}",i)
    print("f(x) = ",f(x1[i]))
    if i%2 == 0:
        P += f(x1[i])
    else:
        I += f(x1[i])

print("\\begin{table}[H]")
print("\t\\centering")
print("\t\\begin{tabular}{|c|c|c|}")
print("\t\\hline")
print("i & $x_i$ & f_i  \\\\ \\hline")
for i in range(len(x1)):
    #print("{} & {} & {}\\\\ \\hline".format(x1[i],f(x1[i]),h1/2*(f(x1[i]) + f(x1[i+1]))))
    print("{} &{} & {} \\\\ \\hline".format(i, x1[i],f(x1[i])))
print("\t\\end{tabular}")
print("\\end{table}")

print(fa)
print(fb)
print(4*I)
print(2*P)
print(h1/3 * (fa+fb+4*I+2*P))
"""
