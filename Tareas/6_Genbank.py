"""
NAME
       6_Genbank.py
VERSION
        [1.0]
AUTHOR
        Melissa Mayén Quiroz
DESCRIPTION
        El programa contiene una función 'resumen' que recibe como parámetros la ubicación
        de un archivo Genbank determiado y una lista de genes y nos imprime la información
        solicitada.
        En la primera parte usando el archivo virus.gb imprimimos  la secuencia, transcripción
        y traducción del gen L.
        En la segunda parte usamos la función resumen con los datos que le pedimos al usuario
        para imprimir la información solicitada.

        *Para el correcto funcionamiento del programa es necesario tener instalado Biopython
        *Dentro del programa importamos algunas librerías:
        - from Bio import SeqIO
        *Dentro de la función resumen se realiza importación de librerías
        *Para el archivo virus.gb la ruta escrita en este programa no funcionará
        en otros equipos por lo que deberá modificarse  
INPUT
        1) Una ruta al archivo Genbank deseado
        2) Una lista de genes a buscar (escritos en mayúsculas y separados por una coma) 
OUTPUT
        1) Se imprime la secuencia, transcripción y traducción del gen L del archivo virus.gb
        2) Se imprime a pantalla la información solicitada de cada gen buscado en el archivo dado
EXAMPLES
        Input
        1) C:\\Users\\Melissa\\Downloads\\virus.gb
        2) M,N,G
        Output
        1) Primera parte: 
        Gen  L
        MDEYSEEKWGDSDEESFGTGKYSDESRIRGLNSVDYNLNSPLIQDDLYYLMERVRGRPVPPIWKAKNWTETIHL ...
        2) Segunda parte:
        Organismo:  Isfahan virus
        Pais:  ['Iran:Isfahan province']
        Aislado:  ['Phlebotomus papatasi']
        Gen:  ['N']
        N gene
        nucleocapsid protein
        ADN:  ATGACTTCTGTAGTA
        ARN:  AUGACUUCUGUAGUA
        Proteína:  MTSVV 
        Gen:  ['M']
        M gene
        matrix protein
        ADN:  ATGAAGAGCTTAAAG
        ARN:  AUGAAGAGCUUAAAG
        Proteína:  MKSLK 
        Gen:  ['G']
        G gene
        glycoprotein
        ADN:  ATGACTTCAGTCTTA
        ARN:  AUGACUUCAGUCUUA
        Proteína:  MTSVL 
GITHUB 
        https://github.com/Melii99/Python_class_2022/tree/main/Tareas/6_Genbank.py
       
"""
# Importación de librerías
from Bio import SeqIO

'''
Resumen es una función que recibe como parámetros la ubicación de un archivo Genbank
determiado y una lista de genes y nos imprime el organismo del que proviene, el país de 
la muestra y la fuente del aislado (a falta de número de aislado), y con la lista de 
genes obtenemos el nombre del gen, el nombre de la proteína que produce, los primeros 
15 nucleótidos de ADN, su respectivo ARN y su respectiva proteína.
'''
def resumen(archivo, genes):
    # Importación de librerías 
    from Bio import SeqIO
    # Parseamos el archivo 
    for gb_record in SeqIO.parse(archivo, "genbank"):

        # Se imprime el organismo, país y aislado
        print('Organismo: ', gb_record.annotations['organism'])
        print('Pais: ', gb_record.features[0].qualifiers['country'])
        # El archivo 'virus.gb' no contiene numero de aislado (ni en annotations ni en features)
        # Por ello se remplazó por la fuente del aislado*
        print('Aislado: ', gb_record.features[0].qualifiers['isolation_source'])

        # Con un for recorremos los features (desde la posición 1 y saltando los CDS)
        for i in range(1, len(gb_record.features),2):

          # Recorremos la lista de genes a buscar
          for gene in genes:

            # Si el gen de la lista y el gen del archivo coinciden, se imprime la información deseada
            if gb_record.features[i].qualifiers['gene'][0] == gene:
              
              # Se imprime el gen que se está buscando desde qualifiers
              print('Gen: ', gb_record.features[i].qualifiers['gene'])
              # Se busca e imprime la el gen y la Proteína que produce en las keywords
              for j in gb_record.annotations['keywords']:
                if j == str(gene + ' gene'):
                  # Gen
                  print(j)
                  #La indice siguiente al gen dentro de keywords corresponde a su proteína
                  x = gb_record.annotations['keywords'].index(j)
                  # Proteína
                  print(gb_record.annotations['keywords'][x+1])
            
              # Se crean la variable a y b para indicar los primeros 15 nucleótidos 
              a = gb_record.features[i].location.nofuzzy_start
              b = a + 15
              # Se imprime la secuencia de DNA trunca
              print('ADN: ', gb_record.seq[a:b])
              # Se imprime la secuencia de RNA trunca
              print('ARN: ', gb_record.seq[a:b].transcribe())
              # Se imprime la secuencia de la proteína trunca
              print('Proteína: ', gb_record.seq[a:b].translate(), '\n') 

'''
Primera parte: Utilizando el archivo virus.gb obtengan la secuencia, 
transcripción y traducción del gen L
'''
print('Primera parte: \n')

# Parseamos el archivo 
for gb_record in SeqIO.parse("C:\\Users\\Melissa\\Downloads\\virus.gb", "genbank"):

  # Con un for recorremos los features (desde la posición 1 y saltando los CDS)
  for i in range(1, len(gb_record.features),2):
    # Encontramos el gen 'L'
    if gb_record.features[i].qualifiers['gene'][0] == 'L':
      print('Gen ',gb_record.features[i].qualifiers['gene'][0])
      # Se obtiene el inicio y final de la secuencia
      a = gb_record.features[i].location.nofuzzy_start
      b = gb_record.features[i].location.nofuzzy_end
      # Se imprime el transcrito de la secuencia
      print(gb_record.seq[a:b].translate())

'''
Segunda parte:

Generalizar una función resumen
Input:
- Nombre del archivo y lista con los nombres de los genes que buscamos 
Output:
- Organismo
- País de la muestra
- Número del aislado
- Nombre del gen y de la proteína que produce
- Los primeros 15 nucleótidos de ADN
- De estos 15 nucleótidos su respectivo ARN
- Y respectiva proteína
'''

print('Segunda parte: \n')

# Le pedimos la ubicación del archivo Genbank al usuario
path = str(input('Ingrese la ubicación del archivo Genbank a procesar: '))
# Le pedimos al usuario los genes de los cuales se desea imprimir la información
genes =  str(input('Ingrese los genes (en mayúsculas) a buscar separados por una coma: '))
genes = genes.split(',')

# Usamos la función resumen con los parámetros indicados por el usuario
resumen(path, genes)