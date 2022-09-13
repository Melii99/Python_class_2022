"""
NAME
       1_ATrich.py
VERSION
        [1.0]
AUTHOR
        Melissa Mayén Quiroz
DESCRIPTION
        La función ATrich nos pide como parámetros el path a la ubicación del archivo (raw) con la secuencia
        a analizar y el numero que indica en tamaño de la región mínima a buscar como región rica en ATs.
        Si existen caracteres raros en la secuencia imprime el caracter y la posición de ocurrencia, de no
        haber caracteres raros en caso de haber regiones ricas se imprimen las regiones y posición de regiones
        ricas en A o T encontradas (en caso de no encontrar se indicará también).
INPUT
        Path a la ubicación del archivo con la secuencia a analizar
        Numero que indica en tamaño de la región mínima a buscar como región rica en ATs
OUTPUT
        Se imprime en pantalla la posición de ocurrencia de caracteres raros en la secuencia y la posición
        de regiones ricas en A o T encontradas
EXAMPLES
        Input
          ATrich('/content/dna.seq', 5)
        Output
          n found at position 5
          r found at position 19
          y found at position 29
GITHUB
        https://github.com/Melii99/Python_class_2022/blob/main/Tareas/1_ATrich.py

"""
# Importar librerias
import argparse
import re

# Paso de argumentos mediante argparse
arg_parser = argparse.ArgumentParser(description="Search AT rich regions")
arg_parser.add_argument("-f", "--file",
                    metavar="path/to/file",
                    help="Archivo con secuencia de ADN",
                    required=True)

arg_parser.add_argument("-s", "--search",
                    help="Tamaño minimo de las regiones de AT",
                    type=int,
                    required=False)

args = arg_parser.parse_args()

def ATrich(path, min):
  """
  La función ATrich nos pide como parámetros el path a la ubicación del archivo (raw) con la secuencia
  a analizar y el numero que indica en tamaño de la región mínima a buscar como región rica en ATs.
  Si existen caracteres raros en la secuencia imprime el caracter y la posición de ocurrencia, de no
  haber caracteres raros en caso de haber regiones ricas se imprimen las regiones y posición de regiones
  ricas en A o T encontradas (en caso de no encontrar se indicará también).

  """
  #Se obtiene la secuencia del archivo de texto
  archivo = open(path, "r")
  lineas = archivo.readlines()
  dnaseq = lineas[0].rstrip('\n').upper()

  #Se obtienen las ocurrencias de caracteres raros en la secuencia
  m = re.finditer("[^ATCG]", dnaseq, re.IGNORECASE)
  counts = len([*re.finditer("[^ATGC]+", dnaseq)])
  #condicion si m tiene caracteres raros
  if counts > 0:
    #Impresión a pantalla de caracteres raos encontrados
    for match in m:
      base = match.group()
      posicion = match.start()
      print(base + " found at position " + str(posicion))

  else:
    #Busqueda de regiones ricas en A o T en caso de no existir caracteres raros
    n = re.finditer("[A|T]+", dnaseq)
    matches = len([*re.finditer("A+|T+", dnaseq)])
    if matches:
    #Impresión a pantalla de regiones ricas en AT encontradas
      for match in n:
        base = match.group()
        posicion = match.start()
        print(base + " found at position " + str(posicion))
    #En caso de no encontrar regiones ricas en AT se indica
    else:
        print("No se encontraron regiones ricas en AT mayores a "+ str(min) + " bases")

  #Se cierra el archivo de texto
  archivo.close()

ATrich(args.file, args.search)
