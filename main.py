import tkinter as tk  # for python 3
import pygubu, cv2
import random
import numpy as np
import math
from tkinter import messagebox
import scipy.stats as stats
import statistics as _stats

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
        matrix_xy = self.obtener_componentes_xy(NUMERO_DE_CROMOSOMA, TAMANIO_GENOTIPOS, LISTA_BINARIOS)
        matrix_xy = self.mapear_componente_xy(NUMERO_DE_CROMOSOMA, matrix_xy, TAMANIO_GENOTIPOS, INTERVAL_X, INTERVAL_Y)

        #Fitness
        lista_valores_z = self.obtener_valores_Z(matrix_xy)
        lista_probabilidades = self.obtener_probabilidades(lista_valores_z)

        #Apareamiento
        hijos_binarios = self.apareamiento(lista_probabilidades, LISTA_BINARIOS)

        #Mutacion
        hijos_binarios_mutados = self.mutacion(hijos_binarios, 0.28, 0.23)


        #Obtener posiciones de los hijos


        print("\n")
        print('Hijos binarios [{}]'.format(len( hijos_binarios )))
        print('Hijos binarios mutados [{}]'.format(len( hijos_binarios )))
        print('Lista binarios [{}]'.format(len(LISTA_BINARIOS)))
        print("\n\n")

    def crear_poblacion(self, _NUMERO_DE_CROMOSOMA, _TAMANIO_GENOTIPOS):
        NUMERO_DE_CROMOSOMA, TAMANIO_GENOTIPOS = _NUMERO_DE_CROMOSOMA, 2**_TAMANIO_GENOTIPOS
        lista_binarios, lista_enteros = [], []

        for i in range(NUMERO_DE_CROMOSOMA):
            numero_random = random.randint(0, TAMANIO_GENOTIPOS-1)
            binario = '{0:032b}'.format(numero_random)
            lista_enteros.append(numero_random)
            lista_binarios.append(binario)

        return lista_binarios, lista_enteros

    def mapear_componente_xy(self, _NUMERO_DE_CROMOSOMA, _MATRIX, _TAMANIO_GENOTIPOS, _INTERVALOS_X, _INTERVALOS_Y):
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

    def obtener_componentes_xy(self, _NUMERO_DE_CROMOSOMA, _TAMANIO_GENOTIPOS, _LISTA_BINARIOS):
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

    def calculo_desvizacion_normal_estandar(self, lista_valores_z, media, varianza):
        lista = []
        lista = list(map(lambda x:(x - media) / varianza, lista_valores_z))
        return lista

    def obtener_desviacion_acumulada(self, lista_valores_z, desviacion_estandar):
        lista = []
        lista = list(map(lambda x:stats.norm.cdf(x, loc=0, scale=desviacion_estandar) , lista_valores_z))
        return lista

    def obtener_probabilidades(self, _lista_valores_z):
        lista_valores_z = np.copy(_lista_valores_z)
        media = _stats.mean(lista_valores_z)
        deviasion_estandar = _stats.pstdev(lista_valores_z)
        print('media: {}, deviasion_estandar: {}'.format(media, deviasion_estandar))
        lista_valores_normal = self.calculo_desvizacion_normal_estandar(lista_valores_z, media, deviasion_estandar) 
        lista_probabilidades = self.obtener_desviacion_acumulada(lista_valores_normal, deviasion_estandar)
        return lista_probabilidades

    def apareamiento(self, _lista_probabilidades, _lista_valores_binarios):
        lista_probabilidades, lista_valores_binarios, hijos =  np.copy(_lista_probabilidades), np.copy(_lista_valores_binarios), []
        TAMANIO_LISTA, buscar_maximo = len(_lista_valores_binarios), True

        for i in range (TAMANIO_LISTA-1):
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
                    print('P: {}, M: {} '.format(lista_valores_binarios[i], lista_valores_binarios[j]))
                    hijo_uno, hijo_dos = self.aparear(lista_valores_binarios[i], lista_valores_binarios[j])
                    hijos.append(hijo_uno); hijos.append(hijo_dos)
                    print("\n")

                j = j + 1

        return hijos

    def aparear(self, _padre, _madre):
        lista_padre, lista_madre = _padre, _madre
        hijo_uno, hijo_dos, acumulador, i = "", "", 0, 0
        lista_point_crossover = self.obtener_point_crossover(len(lista_padre))
        TAMANIO_LISTA_PC = len(lista_point_crossover) + 1

        while i < TAMANIO_LISTA_PC:
            j = 0
            if(i == ( TAMANIO_LISTA_PC - 1 )):
                part_padre, part_madre = lista_padre[acumulador:], lista_madre[acumulador:]
                j = 1
            else:
                salto = acumulador + lista_point_crossover[i]
                part_padre, part_madre = lista_padre[acumulador:salto], lista_madre[acumulador:salto]
            
            acumulador = acumulador + lista_point_crossover[i - j]

            #print('c1: {}, c2: {}'.format(part_padre, part_madre))

            if(((i+1) % 2) != 0):
                hijo_uno, hijo_dos = hijo_uno + part_padre, hijo_dos + part_madre 
            if(((i+1) % 2) == 0):
                hijo_uno, hijo_dos = hijo_uno + part_madre, hijo_dos + part_padre 

            i = i + 1

        print("Hijos:")
        print('uno: {}, dos: {} '.format(hijo_uno, hijo_dos))

        return hijo_uno, hijo_dos

    def mutacion(self, _lista_hijos_binarios, probabilidad_mutar_individuo, probabilidad_mutar_gen):
        lista_hijos_binarios, TAMANIO_LISTA_HIJOS, lista_hijos_binarios_mutados = np.copy(_lista_hijos_binarios), len(_lista_hijos_binarios), []

        for i in range(TAMANIO_LISTA_HIJOS):
            random_probabilidad_mutacion = random.random()
            if(probabilidad_mutar_individuo > random_probabilidad_mutacion):
                list_individuo = list(map(lambda individuo: self.mutar_gen(individuo, probabilidad_mutar_gen), lista_hijos_binarios[i]))
                lista_hijos_binarios_mutados.append(''.join(str(individuo) for individuo in list_individuo))
            else:
                lista_hijos_binarios_mutados.append(lista_hijos_binarios[i])
        return lista_hijos_binarios_mutados
    
    def mutar_gen(self, individuo, probabilidad_mutar_gen):
        random_probabilidad = random.random()
        if(probabilidad_mutar_gen > random_probabilidad):
            if individuo == 0:
                return 1
            else: 
                return 0
        return individuo

    def obtener_point_crossover(self, TAMANIO_CROMOSOMA):
        acumulador_cut, i, cantidad_point_crossover = 0, 0, random.randint(1, 5)
        lista_point_crossover = []

        while i < cantidad_point_crossover:
            longitud_cut = random.randint(3, 29)
            acumulador_cut = acumulador_cut + longitud_cut
            diferencia = TAMANIO_CROMOSOMA - acumulador_cut

            if(diferencia > 0 and i < cantidad_point_crossover):
                lista_point_crossover.append(longitud_cut)
            else:
                i = cantidad_point_crossover
            
            i = i + 1
        print('Cantidad de cruza: {}'.format(lista_point_crossover))

        return lista_point_crossover
        
if __name__ == '__main__':
    root = tk.Tk()
    app = Application(root)
    root.mainloop()