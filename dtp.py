#   Ariel Cerón Gonzáles
#   14/10/2019
#   Estructuras de datos en python para la clase de Computación I de la 
#   licenciatura en matemáticas en la Universidad Abierta y a Distancia de 
#   de México.
#
#   Se presentan cuatro estructuras:
#       1. Listas ligadas
#       2. Pilas
#       3. Colas
#       4. Árboles

#Primera estructura Listas Ligadas

class Nodo(object):

    def __init__(self, dato=None, siguiente=None):
        self.n = dato
        self.nx = siguiente
    
    def __str__(self):
        return(str(self.n))

        
class Llist(object):

    def __init__(self):
        """La estrucutra contiene un unico elemento, el nodo
        de la cabeza todos los demas elementos estan ligados
        a este por una liga"""
        self.head = None
        self.len = 0

    def isEmpty(self):
        """True si la lista esta vacia, False en caso contrario"""
        return(True if not self.head else False)

    def add(self, elemento):
        """Agrega un elemento al final"""     
        if self.head == None:
            self.head = Nodo(elemento)
            self.len += 1
        else:
            tmp = self.head
            while(tmp.nx != None):
                tmp = tmp.nx
            tmp.nx = Nodo(elemento)
            self.len += 1

    def addF(self, elemento):
        """Agrega un elemento al inicio de la lista"""
        if self.head == None:
            self.head = Nodo(elemento)
            self.len += 1
        else:
            aux = self.head
            self.head = Nodo(elemento,aux)#Cambia la cabeza
            self.len += 1

    def remove(self,n):
        """Elimina el elemento de la posición n"""
        tmp = self.len
        self.len = 0
        if (n > tmp):
            print("Se sale de la lista")
        else:
            for elemento in range(n-1):
                aux = self.get()
                self.add(aux)
                self.len += 1
            print(self.get())
            for elemento in range(n+1,tmp):
                aux = self.get()
                self.addF(aux)
                self.len += 1

    def get(self):
        """Obtiene el primer elemento de la lista y lo elimina"""
        if self.isEmpty():
            return False
        else:
            aux = self.head
            self.head = self.head.nx
            return(aux)

    def getAll(self):
        """Imprime todos los elementos ligados al nodo principal"""
        if self.isEmpty():
            return False
        else:
            lista = self.head
            while(lista != None):
                print(lista, end =" => ")
                lista = lista.nx
            print(None)

class p(object):

    def __init__(self):
        """Crea una estructura pila"""
        self.head = None

    def push(self, elemento):
        """Añadir elemento a la pila"""
        if self.head == None:
            self.head = Nodo(elemento)
        else:
            aux = self.head
            self.head = Nodo(elemento,aux)

    def pop(self):
        """Regresa el último elemento que fue ingresado
        y lo elimina de la pila"""
        if self.isEmpty():
            return False
        else:
            aux = self.head
            self.head = self.head.nx
            return(aux)

    def isEmpty(self):
        """True si la lista esta vacia, False en caso contrario"""
        return(True if not self.head else False)
        
