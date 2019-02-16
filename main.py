import tkinter as tk  # for python 3
import pygubu, cv2
import random
import numpy as np
import math
from tkinter import messagebox
import scipy.stats as stats

class Application:
    def __init__(self, master):

        #1: Create a builder
        self.builder = builder = pygubu.Builder()

        #2: Load an ui file
        builder.add_from_file('AG.ui')

        #3: Create the widget using a master as parent
        self.mainwindow = builder.get_object('AG', master)

        builder.connect_callbacks(self)

    def empezar(self):
        # Get all my constants
        TAMANIO_GENOTIPOS = 32
        NUMERO_DE_CROMOSOMA = int(self.builder.get_variable('numCrom').get())
        # NUMERO_DE_GENERACIONES = int(self.builder.get_variable('numeroGeneraciones').get())
        INTERVAL_X = self.builder.get_variable('interval_X').get()
        INTERVAL_Y = self.builder.get_variable('interval_Y').get() 
        
        #Inicializar poblacion
        LISTA_BINARIOS, LISTA_ENTEROS = self.crear_poblacion(NUMERO_DE_CROMOSOMA, TAMANIO_GENOTIPOS)
        matrix_xy = self.obtener_matrix(NUMERO_DE_CROMOSOMA, TAMANIO_GENOTIPOS, LISTA_BINARIOS)
        matrix_xy = self.mapear_valores(NUMERO_DE_CROMOSOMA, matrix_xy, TAMANIO_GENOTIPOS, INTERVAL_X, INTERVAL_Y)
        
        #Fitness
        lista_valores_z = self.obtener_valores_Z(matrix_xy)
        lista_probabilidades = self.obtener_probabilidades(lista_valores_z)

        #Apareamiento
        self.apareamiento(lista_probabilidades, LISTA_BINARIOS)

        print("Valores Z")
        print(lista_valores_z)
        print("Probabilidades:")
        print(lista_probabilidades)

    def crear_poblacion(self, _NUMERO_DE_CROMOSOMA, _TAMANIO_GENOTIPOS):
        NUMERO_DE_CROMOSOMA, TAMANIO_GENOTIPOS = _NUMERO_DE_CROMOSOMA, 2**_TAMANIO_GENOTIPOS
        lista_binarios, lista_enteros = [], []

        for i in range(NUMERO_DE_CROMOSOMA):
            numero_random = random.randint(0, TAMANIO_GENOTIPOS-1)
            binario = '{0:032b}'.format(numero_random)
            lista_enteros.append(numero_random)
            lista_binarios.append(binario)

        return lista_binarios, lista_enteros
    
    def calculo_media(self, _lista_valores_z):
        lista = _lista_valores_z
        TAMANIO_LISTA = len(lista)
        media = 0
        for i in range(TAMANIO_LISTA):
            media = media + lista[i]
        
        media = media / TAMANIO_LISTA
        return media

    def calculo_desviazion_estandar(self, _lista_valores_z, media):
        lista = _lista_valores_z
        TAMANIO_GENOTIPOS = len(lista)
        sumatoria = 0
        
        for i in range(TAMANIO_GENOTIPOS):
            valor = (lista[i] - media)**2
            sumatoria = sumatoria + valor
        
        sumatoria = math.sqrt(sumatoria / TAMANIO_GENOTIPOS)

        return sumatoria

    def calculo_desvizacion_normal_estandar(self, _lista_valores_z, media, varianza):
        TAMANIO_LISTA = len(_lista_valores_z)
        lista, auxLista = [], _lista_valores_z

        for i in range(TAMANIO_LISTA):
            valor = (auxLista[i] - media) / varianza
            lista.append(valor)

        return lista

    def mapear_valores(self, _NUMERO_DE_CROMOSOMA, _MATRIX, _TAMANIO_GENOTIPOS, _INTERVALOS_X, _INTERVALOS_Y):
        NUMERO_DE_CROMOSOMA, MATRIX, TAMANIO_GENOTIPOS = _NUMERO_DE_CROMOSOMA, _MATRIX, _TAMANIO_GENOTIPOS / 2        
        TAMANIO_MATRIX_X = len(MATRIX[0])
        nueva_matrix = np.zeros_like(MATRIX)
        INTERVALOS_X = list(map(int, _INTERVALOS_X.split(":"))) #X[A:B]
        INTERVALOS_Y = list(map(int, _INTERVALOS_Y.split(":"))) #Y[C:D]

        DIFERENCIA_Y = (INTERVALOS_Y[1] - INTERVALOS_Y[0])/(2**TAMANIO_GENOTIPOS) #Y[D:C]
        DIFERENCIA_X = (INTERVALOS_X[1] - INTERVALOS_X[0])/(2**TAMANIO_GENOTIPOS) #X[B:A]

        for i in range(TAMANIO_MATRIX_X):
            numero_x = int(MATRIX[0][i], 2)
            numero_y = int(MATRIX[1][i], 2)
            nueva_matrix[0][i] = INTERVALOS_X[0] + numero_x * DIFERENCIA_X
            nueva_matrix[1][i] = INTERVALOS_Y[0] + numero_y * DIFERENCIA_Y

        return nueva_matrix

    def obtener_valores_Z(self, _MATRIX):
        MATRIX = _MATRIX
        TAMANIO_X = len(MATRIX[0])
        lista_valores_z = []
        for i in range(TAMANIO_X):
            valor_x, valor_y = float(MATRIX[0][i]), float(MATRIX[1][i])
            aux_valor_x, valor_y = valor_x/math.pi,  valor_y/math.pi
            valor_seno, valor_cos = math.sin(math.radians(valor_y)), math.cos(math.radians(aux_valor_x))
            valor_z = (valor_x**2)*(valor_cos + valor_seno)
            lista_valores_z.append(valor_z)
        
        return lista_valores_z

    def obtener_matrix(self, _NUMERO_DE_CROMOSOMA, _TAMANIO_GENOTIPOS, _LISTA_BINARIOS):
        NUMERO_DE_CROMOSOMA, TAMANIO_GENOTIPOS = _NUMERO_DE_CROMOSOMA, int(_TAMANIO_GENOTIPOS / 2)
        lista_binarios_x, lista_binarios_y, LISTA_BINARIOS = [], [], _LISTA_BINARIOS

        for i in range(NUMERO_DE_CROMOSOMA):
            binario = LISTA_BINARIOS[i]
            value_x, value_y = binario[0:TAMANIO_GENOTIPOS], binario[TAMANIO_GENOTIPOS:]
            lista_binarios_x.append(value_x)
            lista_binarios_y.append(value_y)

        #matriz
        #[x, x1, xn]
        #[y, y1, yn]  

        return [lista_binarios_x, lista_binarios_y]

    def obtener_desviacion_acumulada(self, lista_valores_z, desviacion_estandar):
        TAMAMIO_LISTA = len(lista_valores_z)
        lista = []

        for i in range(TAMAMIO_LISTA):
            valor = stats.norm.cdf(lista_valores_z[i], loc=0, scale=desviacion_estandar)
            lista.append(valor)

        return lista

    def obtener_probabilidades(self, _lista_valores_z):
        lista_valores_z = np.copy(_lista_valores_z)
        # print("Valores de z: ")
        # print(lista_valores_z)
        media = self.calculo_media(lista_valores_z)
        deviasion_estandar = self.calculo_desviazion_estandar(lista_valores_z, media)
        # print('media: {}, deviasion_estandar: {}'.format(media, deviasion_estandar))
        lista_valores_normal = self.calculo_desvizacion_normal_estandar(lista_valores_z, media, deviasion_estandar)
        # print("desviazion normal: ")
        # print(lista_valores_z)
        lista_probabilidades = self.obtener_desviacion_acumulada(lista_valores_normal, deviasion_estandar)
        return lista_probabilidades

    def apareamiento(self, _lista_probabilidades, _lista_valores_binarios):
        lista_probabilidades, lista_valores_binarios =  np.copy(_lista_probabilidades), np.copy(_lista_valores_binarios)
        TAMANIO_LISTA, j, lista_apareamiento, numero_apareamientos = len(_lista_valores_binarios), 0, [], 0
        buscar_maximo = True

        for i in range (TAMANIO_LISTA):
            maxima_probabilidad = lista_probabilidades[i]

            if not(buscar_maximo):
                maxima_probabilidad = (maxima_probabilidad-1)

            j = i + 1

            while j < TAMANIO_LISTA:
                probabilidad_random_value = random.random()
                #print('Pm: {}, Pr: {}'.format(maxima_probabilidad, probabilidad_random_value))
                if not(buscar_maximo):
                    maxima_probabilidad = (maxima_probabilidad-1) 

                if(maxima_probabilidad > probabilidad_random_value):
                    numero_apareamientos = numero_apareamientos + 1
                    print("Con------------------------------------------------------------------>")
                    print('P: {}, M: {}'.format(lista_valores_binarios[i], lista_valores_binarios[j]))
                    # hijo_uno, hijo_dos = self, aparear(lista_valores_binarios[i], lista_valores_binarios[j])

                j = j + 1

            if(numero_apareamientos == 0):
                #lista_apareamiento.append(lista_valores_z[i])
                print("Sin------------------------------------------------------------------------->")
                print('PS: {}'.format(lista_valores_binarios[i]))

    # def aparear(self, _padre, _madre):
        # padre, madre = _padre, _madre


    def imprimir_lista(self, LISTA):
        for i in range(len(LISTA)):
            print('i:{}, v:{}'.format(i, LISTA[i]))       

    def imprimir_matrix(self, MATRIX):
        TAMANIO_Y, TAMANIO_X =  len(MATRIX), len(MATRIX[0])
        for j in range(TAMANIO_X):
            print('x:{}, y:{}'.format(MATRIX[0][j], MATRIX[1][j]))

if __name__ == '__main__':
    root = tk.Tk()
    app = Application(root)
    root.mainloop()