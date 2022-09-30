"""
NAME
       4_Secuencias_y_formatos.py
VERSION
        [1.0]
AUTHOR
        Melissa Mayén Quiroz
DESCRIPTION
        El programa contiene una función 'cadena' que toma como parámetros una secuencia
        de nucleótidos y un codón de inicio, y con ellos encuentra los posibles ORFs de la
        secuencia dada, los traduce a cadenas protéicas y regresa la cadena protéica más larga.
        En el programa se le pide primero al usuario ingresar una secuencia y un codón de 
        inicio y usando la función cadena se obtiene la cadena protéica más larga y se imprime
        a pantalla, posteriormente se pide una ruta a un archivo en formato FASTA y un codón de 
        inicio y con la misma función se obtiene la cadena protéica más larga de cada entrada
        y se imprimen a pantalla.
        *Para el correcto funcionamiento del programa es necesario tener instalado Biopython
        *Dentro del programa importamos algunas librerías:
        - from Bio.Seq import Seq
        - from Bio.SeqUtils import nt_search
        - from Bio import SeqIO
INPUT
        1) Una secuencia de nucleótidos 
        2) El codón de inicio para dicha secuencia de nucleótidos
        3) Un 'path' a un archivo FASTA que contenga las secuencias que deseamos procesar
        4) El codón de inicio para las secuencias del archivo FASTA 
OUTPUT
        1) Se imprime en pantalla la cadena protéica más larga encontrada en la secuencia 
        2) Se imprimen a pantalla las cadenas protéicas más largas de cada secuencia
        del archivo FASTA
EXAMPLES
        Input
        1) AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG
        2) ATG
        3) /content/prueba.fasta
        4) ATG
        Output
        1) MGMTPRLGLESLLE
        2) MKRHQFKSWIFELREIKNYHYFLDSWIKFDSVGSFTHIFFHQERFMKLFDPRIWSILLSRDSQGATSNRYFTIKGVVLL ...
           MKRHQFKSWIFELREIKNYHYFLDSWIKFDSVGSFTHIFFHQERFMKLFDPRILSILLSRDSQGATSNRYFMIMIKGVV ...
           MYVNGKILRIFKLKKDGFRQQNFLYPLLFQEYIYSLAHDHNFNSLIFYEPVEIIGYDNKSSLVLVKRLITQMYQQNFFI ...
GITHUB 
        https://github.com/Melii99/Python_class_2022/tree/main/Tareas/4_Secuencias_y_formatos.py
       
"""

# Importación de librerías 
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
from Bio import SeqIO

''' Función cadena
cadena es una fución donde a partir de una secuencia de nucleótidos y un codón 
de inicio obtenemos las cadenas proteicas de todos sus ORFs, y nos regresa
la cadena protéica de mayor longitud encontrada 
'''
def cadena(secuencia, codon_inicio):
  # Buscamos las posiciones donde aparezca  el codón de inicio y las guardamos
  pattern = Seq(codon_inicio)
  positions = nt_search(str(secuencia), pattern)

  # Se crea una lista para guardar las cadenas protéicas encontradas 
  cadena_proteica = []
  # Se recorre positions para obtener la posición del codón de inicio
  for i in range(1,len(positions)):
    sub_sec = secuencia[positions[i]:]
    # Condiciones para que la cadena (nucleótidos) sean múltiplos de 3
    if len(sub_sec) % 3 != 0:
      sub_sec= sub_sec + 'n'
      if len(sub_sec) % 3 != 0:
        sub_sec= sub_sec + 'n'
    # Cada secuencia es traducida y se agrega a la lista 'cadena_proteica' como string   
    cadena_proteica.append(str(sub_sec.translate(to_stop=True)))

  # Se crea una lista para guardar la longitud de cada cadena protéica
  len_cadena = []
  # Se guarda la longitud de cada cadena
  for cadena in cadena_proteica:
    len_cadena.append(len(cadena))

  # Se obtiene la cadena protéica más larga encontrada  
  long_prot = cadena_proteica[len_cadena.index(max(len_cadena))]
  return long_prot


''' Parte 1
Se pide al usuario ingresar la secuencia y el codón de inicio y se imprime a pantalla 
la cadena protéica más larga encontrada en la secuencia
'''
# Pedimos la secuencia a procesar al usuario
sequence = str(input('Ingrese la secuencia nucleotídica a procesar: '))
# Creamos el objeto secuencia a partir de la cadena que se nos proporcione
sequence = Seq(sequence)
# pedimos el codón de inicio a usar
codon = str(input('Ingrese el codón de inicio de su secuencia:  '))
# Usamos la función para obtener la cadena protéica más larga con los datos que
# nos proporcionó el usuario 
prot_chain = cadena(sequence, codon)
print(prot_chain)

''' Parte 2
Se pide al usuario ingresar la ruta a un archivo con formato FASTA y un codón de 
inicio y se imprime a pantalla la cadena protéica más larga encontrada en cada
record del archivo 
'''
# Pedimos una dirección de archivo FASTA 
filename = str(input('Ingrese un path a un archivo FASTA: '))
# Pedimos el codón de inicio a usar
codon = str(input('Ingrese el codón de inicio de su archivo FASTA:  '))
# Por cada record (seq) usaremos la funcion cadena para obtener la cadena protéica más larga
for seq_record in SeqIO.parse(filename, "fasta"):
  prot_chain = cadena(seq_record.seq, codon)
  print(prot_chain)