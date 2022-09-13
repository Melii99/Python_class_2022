"""
NAME
       2_Gencode.py
VERSION
        [1.0]
AUTHOR
        Melissa Mayén Quiroz
DESCRIPTION
        El programa busca resolver el problema de Rosalind relacionada al código genético:
        The 20 commonly occurring amino acids are abbreviated by using 20 letters from the English alphabet 
        (all letters except for B, J, O, U, X, and Z). 
        Protein strings are constructed from these 20 symbols. Henceforth, the term genetic string will incorporate protein 
        strings along with DNA strings and RNA strings.
        The RNA codon table dictates the details regarding the encoding of specific codons into the amino acid alphabet.
        Given: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp).
        Return: The protein string encoded by s.
INPUT
        Una cadena de texto que contenga la secuencia a traducir
OUTPUT
        Se imprime en pantalla la posición de ocurrencia de caracteres raros en la secuencia y la posición 
        de regiones ricas en A o T encontradas
EXAMPLES
        Input
        AGUACCCGUAUUAACGGGUGA 
        Output
         ['M', 'A', 'M', 'A', 'P', 'R', 'T', 'E', 'I', 'N', 'S', 'T', 'R', 'I', 'N', 'G', '_']      
GITHUB 
        https://github.com/Melii99/Python_class_2022/tree/main/Tareas/2_Gencode.py

"""

#Importamos el módulo re
import re

#Se genera el diccionario del código genético
gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T',
    'ACC':'T', 'ACG':'T', 'ACT':'T', 'AAC':'N', 'AAT':'N',
    'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R',
    'AGG':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H',
    'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R',
    'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V',
    'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G',
    'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S',
    'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L',
    'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

#Se usa la cadena en mayusculas y se reemplazan las Us por Ts en el transcrito (para coincidir con los codones del diccionario)
s = input('Introduzca la secuencia a traducir:  ')
s = s.upper()
s = s.replace('U','T')

#Inicialización de contadores y creación de listas vacía
a = 0
b = 3
codones = []
transcrito = []

#Condición para tener codones completos (secuencia divisible entre 3)
if len(s)/3 is float:
  print('La secuencia contiene un codón incompleto')
  
else: 
  #Condición para encontrar e imprimir caracteres raros en la secuencia
  m = re.finditer("[^ATGC]", s)
  if len(list(m)) > 0:  
    #Impresión a pantalla de caracteres raros encontrados
    for match in m:
      base = match.group()
      posicion = match.start()
      print(base + " found at position " + str(posicion))
#Si hay codones completos y no hay caracteres raros se divide la cadena de texto proporcionada en codones (3 letras)
  else: 
    for i in range(int(len(s)/3)):
      codones.append(s[a:b])
      a += 3
      b += 3

#Se traduce cada codon de la lista de codones para generar el trancrito
for codon in codones:
  if codon in gencode.keys():
    transcrito.append(gencode[codon])
  else: transcrito.append('El codón encontrado no codifica un aminoácido')

#Se imprime a pantalla el transcrito
print("".join(transcrito))
