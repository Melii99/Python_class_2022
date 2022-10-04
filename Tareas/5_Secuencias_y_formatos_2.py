"""
NAME
       5_Secuencias_y_formatos_2.py
VERSION
        [1.0]
AUTHOR
        Melissa Mayén Quiroz
DESCRIPTION
        El programa contiene una función 'umbral' que toma como parámetros una ruta
        a un archivo con formato FASTQ y un umbral de calidad deseada, y con ellos 
        selecciona las secuencias cuyos valores (cada nucleótido) superen el umbral 
        y regresa el número de secuencias que pasan el umbral y sus respectivos IDs.
        En el programa se le pide primero al usuario ingresar una ruta a un archivo 
        con formato FASTQ y un umbral de calidad deseada y usando la función umbral 
        obtenemos obtenemos el numero de secuencias y sus IDs los cuales guargaremos
        en un archivo para lo que pediremos una ruta al usuario.
        En la segunda parte se obtiene la información solicitada del archivo virus.gb
        y se imprime a pantalla.
        
        *Para el correcto funcionamiento del programa es necesario tener instalado Biopython
        *Dentro del programa importamos algunas librerías:
        - from Bio.Seq import Seq
        *Dentro de la función umral se realiza importación de librerías
        *Para el archivo virus.gb la ruta escrita en este programa no funcionará
        en otros equipos por lo que deberá modificarse  
INPUT
        1) Una ruta al archivo FASTQ que se desea evaluar 
        2) Un umbral (numérico) de calidad deseado 
        3) Una ruta para crear el archivo donde guardar el numero de secuencias
          que pasan el umbral de calidad y sus IDs
OUTPUT
        1) Archivo con el numero de secuencias que pasan el umbral de calidad y sus IDs
        2) Se imprimen a pantalla la información solicitada del archivo virus.gb
EXAMPLES
        Input
        1) C:\\Users\\Melissa\\Downloads\\sample.fastq
        2) 29
        3) C:\\Users\\Melissa\\Downloads\\IDs_fastq
        Output
        1) Archivo IDs_fastq
        2) Fha:  13-AUG-2018
           Organismo del que proviene:  Isfahan virus
           País de orígen: ['Iran:Isfahan province']
GITHUB 
        https://github.com/Melii99/Python_class_2022/tree/main/Tareas/5_Secuencias_y_formatos_2.py
       
"""
#Importación de librerías (necesaria para la 2ª parte)
from Bio import SeqIO

'''
Primera parte: La función umbral (parámetros: archivo y calidad) selecciona
las secuencias cuyos valores (cada nucleótido) superen el umbral y regresa el 
número de secuencias que pasan el umbral y sus IDs
'''
def umbral(archivo, calidad):
  #Importación de librerías
  from Bio import SeqIO
  # Contador de secuencias
  n = 0
  # Lista de ids seleccionados 
  ids = []
   
  # Se recorren todos los records del archivo fastq
  for record in SeqIO.parse(archivo, "fastq"):
    # Se obtiene el score de cada record
    score = record.letter_annotations["phred_quality"]

    # Lista de bases individuales que no pasan la calidad requerida
    no_calidad = []
    #Se recorre la calidad de cada base
    for i in score:
    # Si la calidad en algúna base es menor a la requerida se agrega a la lista
      if i < calidad:
        no_calidad.append(i)
    
    # Si la calidad requerida de todas las bases es igual o mayor se agregan a la lista los IDs    
    if len(no_calidad) == 0:
      ids.append(record.id)
      # Se incrementa el contador
      n += 1

  return(n, ids)


'''
Pedimos los parámetros y usamos la funcion umbral. Obtenemos y guardamos el 
número de secuencias que pasan el umbral y sus IDs en un archivo 
para el cual pedirémos una ruta al usuario
'''
# Le pedimos la ruta al archivo FASTQ al usuario 
archivo = str(input('Ingrese la ruta del archivo FASTQ a procesar: '))
# Pedimos la calidad deseada de las secuencias
quality = int(input('Ingrese la calidad deseada para sus secuencias: '))

# Usamos la función con los datos proporcionados por el usuario
cantidad, ids = umbral(archivo,quality)

# Pedimos una ruta para guardar el archivo que contenga los IDs
path = str(input('Ingrese una ruta para el nuevo archivo que contenga los IDs de la calidad deseada: '))
# Creamos un archivo para escribir los IDs
archivo_calidad = open(path, "w")
# Escribimos la cantidad de secuencias que pasaron el umbral de calidad
archivo_calidad.write("Numero de secuencias con la calidad requerida: " + str(cantidad))
# Recorremos la lista con los IDs y los escribimos en el archivo
for id in ids:
  archivo_calidad.write("\n ID: " + id)

archivo_calidad.close()

'''
Segunda parte: Obtener información de las anotaciones Utilizando el archivo 
virus.gb obtengan de qué organismo viene y la fecha Incluir algún valor de features
'''
print('Segunda parte: Información de archivo virus.gb \n')

# Procesamos y parseamos el archivo virus.gb
for gb_record in SeqIO.parse("C:\\Users\\Melissa\\Downloads\\virus.gb", "genbank"):
    pass
# Imprimimos la fecha
print('Fha: ', gb_record.annotations['date'])
# Imprimimos esl organismo de donde viene
print('Organismo del que proviene: ', gb_record.annotations['organism'])

#Imprimimos el país de orígen desde features
print('País de orígen:', gb_record.features[0].qualifiers['country'])
