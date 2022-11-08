"""
NAME
       7_NCBI_1.py
VERSION
        [1.0]
AUTHOR
        Melissa Mayén Quiroz
DESCRIPTION
        Dentro de este programa se encuentran 2 funciones: crear_termino y buscar_db.
        La primera recibe un diccionario con el organismo y los genes que se desean buscar
        y genera los términos necesarios para la búsqueda con Entrez (egquerry u esearch), 
        devolviendo una lista con ellos.
        La segunda recibe la lista de términos generada por la primera función y realiza la
        búsqueda de estos regresando un diccionario con el organismo (keys) y los IDs 
        encontrados en la base de datos correspondiente (values).
        Usando estas funciones hacemos una búsqueda como ejemplo.
        *Para el correcto funcionamiento del programa es necesario tener instalado Biopython
        *Dentro de la función buscar_db se realiza importación de Entrez 
        (el email es modificable) 
INPUT
        1) Un diccionario que contenga el organismo (key) y los genes a buscar (values)
        2) Output de la 1º función: Una lista con los términos a buscar
        OUTPUT
        1) Una lista con los términos a buscar
        2) Un diccionario con el organismo (key) y una lista de bases de datos y sus 
            IDs correspondientes (values)
EXAMPLES
        Input
        1) {'Drosophila melanogaster':'dpp,Ubx,Adh',
        'Caenorhabditis elegans': 'daf-16,daf-2,daf-12,hsf-1'}
        2) Output de la 1º parte (lista de términos a buscar)
        Output
        1) Primera parte: 
        ['(Drosophila melanogaster[Orgn] AND dpp[Gene] AND Ubx[Gene] AND Adh[Gene] )',
        '(Caenorhabditis elegans[Orgn] AND daf-16[Gene] AND daf-2[Gene] AND daf-12[Gene] ...
        2) Segunda parte:
        {'Drosophila melanogaster': ['pmc',
        ['8343187', '7463283', '6366919', '4949521', '2965233', '2813253', '2712966', ...
        'books',
        ['3312398'],
        'genome',
        ['47']],
        'Caenorhabditis elegans': ['pmc',
        ['9465036', '9441740', '9470895', '9304931', '9217932', '9168093', '8885179', ...
        'books',
        ['4684026', '4078908', '3069336'],
        'genome',
        ['41']]}
GITHUB 
        https://github.com/Melii99/Python_class_2022/tree/main/Tareas/7_NCBI_1.py
       
"""

'''
Primera parte: La función crear_termino recibe un diccionario con el organismo y los genes 
que se desean buscar y genera los términos de búsqueda devolviendolos en una lista.
* Se incluyen todos los genes correspondientes al organismo en el término a buscar
** No sólamente 2 (como en el ejemplo)
'''
def crear_termino(informacion):
  terminos = []
  # Iteramos sobre cada organismo
  for i in informacion:
    # Obtenemos los genes correspondientes a cada organismo
    genes = str(informacion[i])
    genes = genes.split(',')
    # Escribimos el organismo en el término
    cadena = '(' + i + '[Orgn] '
    # Escribimos los genes correspondientes a dicho organismo en el término
    for gen in genes:
      cadena = cadena + 'AND ' + gen + '[Gene] '
    cadena = cadena + ')'
    # Agregamos cada termino completo a la lista de términos
    terminos.append(cadena)

  return(terminos)

'''
Segunda parte: La función buscar_db recibe la lista de términos generada por la primera función 
y realiza la búsqueda de cada uno regresando un diccionario con el organismo (key) y los IDs 
encontrados en la base de datos correspondiente (value).
'''
def buscar_db(terminos):
  #Realizamos la importación de librerías
  from Bio import Entrez
  #Proporcionamos un email
  Entrez.email = "mmayen@lcg.unam.mx"  # El e-mail es modificable 

  # Creamos un diccionario vacío para añadir los organismos, la base de datos y sus ids 
  IDs = {}
  # Se recorre cada término a buscar
  for termino in terminos:
    # Se realiza la búsqueda para saber en qué bases de datos se encuentran resultados
    handle = Entrez.egquery(term=termino)
    record = Entrez.read(handle)
    handle.close()
    # A partir del término obtenemos a qué organismo pertenecen los genes buscados
    organismo = termino.split('[Orgn]')[0]
    organismo = organismo.replace('(','')

    # Se crea una lista vacía para añadir los ids de cada base de datos de cada organismo
    ids = []
    # Recorremos los resultados que nos arrojó la búsqueda anterior
    for resultados in record['eGQueryResult']:
      # Si existen resultados en una base de datos buscamos los IDs
      if (resultados['Count'] != '0') and (resultados['Status'] == 'Ok'):
        handle2 = Entrez.esearch(db=resultados['DbName'], term=termino)
        record2 = Entrez.read(handle2)
        handle2.close()
         # Se agrega cada base de datos junto con los IDs encontrados a la lista ids
        ids.append(resultados['DbName'])
        ids.append(record2['IdList'])

    # Añadimos el organismo como llave a el diccionario de IDs
    # Añadimos como valor los ids correspondientes a cada base de dato 
    IDs[organismo] = ids

  return(IDs)


'''
Ejémplo de búsqueda
'''
# Objeto a leer en la función crear_termino (diccionario)
info = {'Drosophila melanogaster':'dpp,Ubx,Adh',
        'Caenorhabditis elegans': 'daf-16,daf-2,daf-12,hsf-1'}

# usamos la función crear_termino para crear los términos del objeto info
terminos = crear_termino(info)

# Usamos la funcion buscar_db para encontrar los IDs en diferentes bases de datos de cada organismo
IDs = buscar_db(terminos)