import tkinter as tk  # for python 3
import pygubu
import random
import numpy as np
import math
from tkinter import messagebox
import scipy.stats as stats
import statistics as _stats
import matplotlib.pyplot as plt
import os
from os.path import isfile, join
import cv2

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
        TAMANIO_GENOTIPOS, TAMANIO_POBLACION, MAXIMOS, GENERACIONES_REALES, i = 32, 0.60, True, 1, 0
        NUMERO_DE_CROMOSOMA = int(self.builder.get_variable('numCrom').get())
        NUMERO_DE_GENERACIONES = int(self.builder.get_variable('numeroGeneraciones').get())
        INTERVAL_X = self.builder.get_variable('interval_X').get()
        INTERVAL_Y = self.builder.get_variable('interval_Y').get()
        PORCENTAJE_INDIVIDUO = self.builder.get_variable('porcentaje_individuo').get() 
        PORCENTAJE_GEN = self.builder.get_variable('porcentaje_gen').get() 

        media_mejor, media_peor, media_normal = [], [], []

        print('PI: {}, PG: {}'.format( PORCENTAJE_INDIVIDUO, PORCENTAJE_GEN ))
        print('IX: {}, IY: {}'.format(INTERVAL_X, INTERVAL_Y))

        #Inicializar poblacion
        hijos, lista_enteros = self.crear_poblacion(NUMERO_DE_CROMOSOMA, TAMANIO_GENOTIPOS)
        
        #Comenzar algoritmo
        while i < NUMERO_DE_GENERACIONES:
            print("h")
            try:
                #Obtener componente XY mapeados
                matrix_xy = self.crear_componentes_xy(hijos, TAMANIO_GENOTIPOS, INTERVAL_X, INTERVAL_Y)

                #Fitness
                lista_valores_z = self.obtener_valores_Z(matrix_xy)
                lista_probabilidades = self.obtener_probabilidades(lista_valores_z)
                hijos_probabilidad = self.agregar_probabilidad(lista_probabilidades, hijos)

                #Apareamiento
                hijos_binarios = self.apareamiento(hijos_probabilidad, len(hijos_probabilidad), MAXIMOS, media_mejor, media_peor, media_normal, lista_valores_z )

                #Mutacion
                hijos = self.mutacion(hijos_binarios, PORCENTAJE_INDIVIDUO, PORCENTAJE_GEN)
                
                coincidencia, matrix_xy = self.buscar_coincidencia(hijos, TAMANIO_GENOTIPOS, INTERVAL_X, INTERVAL_Y)

                #Graficar
                self.create_image_hijos(matrix_xy, 'imagen_{}'.format(i), INTERVAL_X, INTERVAL_Y)

                if(coincidencia):
                    i = NUMERO_DE_GENERACIONES
                else:
                    GENERACIONES_REALES = GENERACIONES_REALES + 1
                
                i = i + 1
            except ex:
                print(ex) 
                i = NUMERO_DE_GENERACIONES
        
        self.crear_chart_media(media_normal, media_mejor, media_peor, GENERACIONES_REALES)

        self.crear_video()

    def crear_poblacion(self, _NUMERO_DE_CROMOSOMA, _TAMANIO_GENOTIPOS):
        NUMERO_DE_CROMOSOMA, TAMANIO_GENOTIPOS = _NUMERO_DE_CROMOSOMA, 2**_TAMANIO_GENOTIPOS
        lista_binarios, lista_enteros = [], []

        for i in range(NUMERO_DE_CROMOSOMA):
            numero_random = random.randint(0, TAMANIO_GENOTIPOS)
            binario = '{0:032b}'.format(numero_random)
            lista_enteros.append(numero_random)
            lista_binarios.append(binario)

        return lista_binarios, lista_enteros

    def mapear_componente_xy(self, _NUMERO_DE_CROMOSOMA, _MATRIX, _TAMANIO_GENOTIPOS, _INTERVALOS_X, _INTERVALOS_Y):
        NUMERO_DE_CROMOSOMA, MATRIX, TAMANIO_GENOTIPOS = _NUMERO_DE_CROMOSOMA, np.copy(_MATRIX), _TAMANIO_GENOTIPOS / 2        
        TAMANIO_MATRIX_X = len(MATRIX[0])
        nueva_matrix = np.zeros_like(MATRIX)
        INTERVALOS_X = list(map(int, _INTERVALOS_X.split(":"))) #X[A:B]
        INTERVALOS_Y = list(map(int, _INTERVALOS_Y.split(":"))) #Y[C:D]

        DIFERENCIA_Y = (INTERVALOS_Y[1] - INTERVALOS_Y[0])/(2**TAMANIO_GENOTIPOS) #Y[D:C]
        DIFERENCIA_X = (INTERVALOS_X[1] - INTERVALOS_X[0])/(2**TAMANIO_GENOTIPOS) #X[B:A]

        #print('D_Y: {}, D_X: {}'.format(DIFERENCIA_X, DIFERENCIA_Y))

        for i in range(TAMANIO_MATRIX_X):
            numero_x = int(MATRIX[0][i], 2)
            numero_y = int(MATRIX[1][i], 2)
            #print('D_nx:{}, D_ny: {}'.format(numero_x, numero_y))
            nueva_matrix[0][i] = INTERVALOS_X[0] + numero_x * DIFERENCIA_X
            nueva_matrix[1][i] = INTERVALOS_Y[0] + numero_y * DIFERENCIA_Y
            #print('M_nx:{}, M_ny: {}'.format(nueva_matrix[0][i], nueva_matrix[1][i]))

        return nueva_matrix

    def obtener_valores_Z(self, _MATRIX):
        MATRIX = _MATRIX
        TAMANIO_X = len(MATRIX[0])
        lista_valores_z = []
        for i in range(TAMANIO_X):
            valor_x, valor_y = float(MATRIX[0][i]), float(MATRIX[1][i])
            valor_seno, valor_cos = math.sin(math.radians(valor_x)), math.cos(math.radians(valor_y))
            valor_z = (valor_x**2)*(valor_cos + valor_seno)
            lista_valores_z.append(valor_z)

        return lista_valores_z

    def obtener_componentes_xy(self, _NUMERO_DE_CROMOSOMA, _TAMANIO_GENOTIPOS, _lista_binarios):
        NUMERO_DE_CROMOSOMA, TAMANIO_GENOTIPOS = _NUMERO_DE_CROMOSOMA, int(_TAMANIO_GENOTIPOS / 2)
        lista_binarios_x, lista_binarios_y, lista_binarios = [], [], np.copy(_lista_binarios)
        print('TG: {}'.format(TAMANIO_GENOTIPOS))
        print('Cl: {}'.format(len(_lista_binarios)))
        for i in range(NUMERO_DE_CROMOSOMA):
            binario = lista_binarios[i]
            value_x, value_y = binario[0:TAMANIO_GENOTIPOS], binario[TAMANIO_GENOTIPOS:]
            lista_binarios_x.append(value_x)
            lista_binarios_y.append(value_y)

        #print('Bx: {}, By: {}'.format(lista_binarios_x, lista_binarios_y))

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
        #print('media: {}, deviasion_estandar: {}'.format(media, deviasion_estandar))
        lista_valores_normal = self.calculo_desvizacion_normal_estandar(lista_valores_z, media, deviasion_estandar) 
        lista_probabilidades = self.obtener_desviacion_acumulada(lista_valores_normal, deviasion_estandar)
        return lista_probabilidades

    def agregar_probabilidad(self, _lista_probabilidad, _hijos):
        lista_probabilidad, lista_hijos_probabilidad =  _lista_probabilidad, []
        TAMANIO_PROBABILIDAD = len(_hijos)

        for i in range (TAMANIO_PROBABILIDAD):
            tupla_hijos = (lista_probabilidad[i], _hijos[i])
            lista_hijos_probabilidad.append(tupla_hijos)

        return lista_hijos_probabilidad

    def apareamiento(self, _lista_binarios_probabilidad, _TAMANIO_LISTA, _MAXIMOS, list_media_mejor, list_media_peor, list_media, lista_valores_z):
        lista_binarios_probabilidad, hijos =  _lista_binarios_probabilidad, []
        TAMANIO_LISTA, BUSCAR_MAXIMO, exiteCruze = _TAMANIO_LISTA, _MAXIMOS, False
        media_mejor, cantidad_mejor, media_peor, cantidad_peor, media, cantidad_media = 0, 0, 0, 0, 0, 0

        for i in range (TAMANIO_LISTA-1):
            maxima_probabilidad,_ = _lista_binarios_probabilidad[i]
            exiteCruze = False

            if not(BUSCAR_MAXIMO):
                maxima_probabilidad = (1-maxima_probabilidad)

            j = i + 1

            while j < TAMANIO_LISTA:
                probabilidad_random_value = random.random()
                #print('Pm: {}, Pr: {}'.format(maxima_probabilidad, probabilidad_random_value))
                if(maxima_probabilidad > probabilidad_random_value):
                    _, binario_uno = lista_binarios_probabilidad[i]
                    _, binario_dos = lista_binarios_probabilidad[j]
                    #print('P: {}, M: {} '.format(binario_uno, binario_dos))
                    hijo_uno, hijo_dos = self.aparear(binario_uno, binario_dos)
                    hijos.append(hijo_uno); hijos.append(hijo_dos)
                    exiteCruze = True

                j = j + 1
            
            if( exiteCruze ):
                media_mejor = media_mejor + lista_valores_z[i]
                cantidad_mejor = cantidad_mejor + 1
            else:
                media_peor = media_peor + lista_valores_z[i]
                cantidad_peor = cantidad_peor + 1

            media = media + lista_valores_z[i]
            cantidad_media = cantidad_media + 1

        media_mejor, media_peor, media = self.obtener_promedios( media_mejor, cantidad_mejor, media_peor, cantidad_peor, media, cantidad_media )
        list_media_mejor.append(media_mejor); list_media_peor.append(media_peor); list_media.append(media)

        return  self.meteoro(hijos, 20)

    def obtener_promedios(self, mejor, cantidad_mejor, peor, cantidad_peor, media, cantidad_media):
        if(cantidad_mejor != 0):
            mejor =  (mejor/cantidad_mejor)
        else:
            mejor = 0

        if(cantidad_peor != 0):
            peor = (peor/cantidad_peor)
        else:
            peor = 0

        if(cantidad_media != 0):
            media =  (media/cantidad_media)
        else:
            media = 0

        return mejor, peor, media

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

        #print("Hijos:")
        #print('uno: {}, dos: {} '.format(hijo_uno, hijo_dos))

        return hijo_uno, hijo_dos

    def mutacion(self, _lista_hijos_binarios, probabilidad_mutar_individuo, probabilidad_mutar_gen):
        lista_hijos_binarios, TAMANIO_LISTA_HIJOS, lista_hijos_binarios_mutados = np.copy(_lista_hijos_binarios), len(_lista_hijos_binarios), []

        for i in range(TAMANIO_LISTA_HIJOS):
            random_probabilidad_mutacion = random.random()
            if(probabilidad_mutar_individuo > random_probabilidad_mutacion):
                #print("Hm: {}".format(lista_hijos_binarios[i]))
                list_individuo = list(map(lambda individuo: self.mutar_gen(individuo, probabilidad_mutar_gen), lista_hijos_binarios[i]))
                #print('LHC: {}'.format(list_individuo))
                lista_hijos_binarios_mutados.append(''.join(str(individuo) for individuo in list_individuo))
                #print('LHC: {}'.format(lista_hijos_binarios_mutados))
            else:
                lista_hijos_binarios_mutados.append(lista_hijos_binarios[i])

        #print("Mutados: {}".format(lista_hijos_binarios_mutados))

        return lista_hijos_binarios_mutados
    
    def mutar_gen(self, individuo, probabilidad_mutar_gen):
        random_probabilidad = random.random()
        #print('Individuo:{}'.format(individuo))
        #print("PG:{}, RP:{}".format(probabilidad_mutar_gen, random_probabilidad))
        individuo = int(individuo)
        if(probabilidad_mutar_gen > random_probabilidad):
            if individuo == 0:
                #print("Ahora es: 1")
                return 1
            else: 
                #print("Ahora es: 0")
                return 0
        #print("Lo mismo: {}".fomat(individuo))
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
        #print('Cantidad de cruza: {}'.format(lista_point_crossover))

        return lista_point_crossover

    def meteoro(self, hijos, MAXIMOS):

        sobrevivientes = hijos[:MAXIMOS]

        return sobrevivientes

    def crear_componentes_xy(self, _hijos, CANTIDAD_GENOTIPOS, INTERVAL_X, INTERVAL_Y):
        CANTIDAD_HIJOS = len(_hijos)
        matrix_xy = self.obtener_componentes_xy(CANTIDAD_HIJOS, CANTIDAD_GENOTIPOS, _hijos)
        matrix_xy = self.mapear_componente_xy(CANTIDAD_HIJOS, matrix_xy, CANTIDAD_GENOTIPOS, INTERVAL_X, INTERVAL_Y)
        return matrix_xy

    def buscar_coincidencia(self, _hijos, CANTIDAD_GENOTIPOS, INTERVAL_X, INTERVAL_Y):
        CANTIDAD_HIJOS, cantidad_coincidencias = len(_hijos), 0

        matrix_xy = self.crear_componentes_xy(_hijos, CANTIDAD_GENOTIPOS, INTERVAL_X, INTERVAL_Y)

        VALUE_INICIAL_X = round(float(matrix_xy[0][0]),2)
        VALUE_INICIAL_Y = round(float(matrix_xy[1][0]),2)

        for i in range( CANTIDAD_HIJOS ):
            coordenada_x = round(float(matrix_xy[0][i]),2)
            coordenada_y = round(float(matrix_xy[1][i]),2)
            #print('Coor_x: {}, Coor_y:{}'.format(coordenada_x, coordenada_y))
            if((VALUE_INICIAL_Y == coordenada_y) and (VALUE_INICIAL_X == coordenada_x )):
                cantidad_coincidencias = cantidad_coincidencias + 1

        print('cantidad_coincidencias: {}, CH: {}'.format( cantidad_coincidencias, CANTIDAD_HIJOS ))

        if( cantidad_coincidencias == CANTIDAD_HIJOS ):
            return True, matrix_xy
        else:
            return False, matrix_xy

    def create_image_hijos(self, _matrix_xy, name_image, _INTERVALOS_X, _INTERVALOS_Y):
        CANTIDAD_HIJOS = len(_matrix_xy[0])
        INTERVALOS_X = list(map(int, _INTERVALOS_X.split(":"))) #X[A:B]
        INTERVALOS_Y = list(map(int, _INTERVALOS_Y.split(":"))) #Y[C:D]

        print("IX: {}, IY: {}".format( INTERVALOS_X, INTERVALOS_Y ))

        # Create plot
        fig = plt.figure()
        ax, area = fig.add_subplot(111), 8**2
        plt.title('Generacion: image_{}'.format(name_image))
        plt.ylim( int(INTERVALOS_Y[0]) , int(INTERVALOS_Y[1]) )
        plt.xlim( int(INTERVALOS_X[0]) , int(INTERVALOS_X[1]) )

        for i in range( CANTIDAD_HIJOS ):
            coordenada_x = round(float(_matrix_xy[0][i]),2)
            coordenada_y = round(float(_matrix_xy[1][i]),2)

            #print("Cx: {}, Cy:{}".format(coordenada_x, coordenada_y))
            ax.scatter(coordenada_x, coordenada_y, s=area, c="green", alpha=0.5)

        plt.savefig('image/{}'.format(name_image))
        plt.close(fig)

    def crear_chart_media(self, list_media, list_media_mejor, list_media_peor, CANTIDAD_GENERACIONES):
        CANTIDAD_MEDIA = len(list_media)
        lista_generaciones = []

        print('Cantidad:{}, CANTIDAD_MEDIA: {} '.format(CANTIDAD_GENERACIONES, CANTIDAD_MEDIA))

        for i in range(CANTIDAD_GENERACIONES):
            lista_generaciones.append(i+1)

        print('Min: {}'.format(min(list_media_mejor + list_media_peor + list_media)))
        print('Max: {}'.format(max(list_media_mejor + list_media_peor + list_media)))

        valor_minimo = min(list_media_mejor + list_media_peor + list_media)
        valor_maximo = max(list_media_mejor + list_media_peor + list_media)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.title('Valor de la media')
        plt.xlabel('Generaciones')
        plt.ylabel('Fitness')
        ax.set_ylim(bottom=valor_minimo, top=valor_maximo)
        plt.plot(lista_generaciones, list_media_peor, 'ro-')
        plt.plot(lista_generaciones, list_media_mejor, 'go-')
        plt.plot(lista_generaciones, list_media,'b')
        plt.show()


    def imprimir_matrix(self, matrix, _TAMANIO_FILAS):
        for i in range(_TAMANIO_FILAS):
            print("x:{}, y:{}".format(matrix[0][i], matrix[1][i]))

    def crear_video(self):
        pathIn= 'image/'
        pathOut = 'video/video.avi'
        fps = 1.0
        frame_array = []
        files = [f for f in os.listdir(pathIn) if isfile(join(pathIn, f))]
        #for sorting the file names properly
        files.sort(key = lambda x: int(x[7:-4]))
        for i in range(len(files)):
            filename=pathIn + files[i]
            #reading each files
            img = cv2.imread(filename)
            height, width, layers = img.shape
            size = (width,height)
            #inserting the frames into an image array
            frame_array.append(img)
    
        out = cv2.VideoWriter(pathOut,cv2.VideoWriter_fourcc(*'DIVX'), fps, size)
    
        for i in range(len(frame_array)):
            # writing to a image array
            out.write(frame_array[i])
        out.release()

if __name__ == '__main__':
    root = tk.Tk()
    app = Application(root)
    root.mainloop()