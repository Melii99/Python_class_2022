"""
NAME
       8_NCBI_2.py
VERSION
        [1.0]
AUTHOR
        Melissa Mayén Quiroz
DESCRIPTION
        Dentro del programa se encuentra la función extraer_y_guardar que recibe un diccionario:
        el output de la función 2 (buscar_db) de la tarea anterior (7_NCBI_1).
        A partir de este se realiza una búsqueda sobre cada organismo de cada ID en bases de datos 
        que contienen artículos. Se extraen y agregan los títulos de los artículos encontrados y
        los IDs de las citas. Con la lista de IDs se realiza una nueva búsqueda de la cual se 
        extraen los títulos y abstracts de las citas y se agregan al archivo. En caso de no 
        encontrarse citas, también se escribe en el archivo.
        Al final del programa se encuentra como comentario lo necesario para realizar una búsqueda 
        de ejemplo, usada para generar el archivo archivo_articulos encontrado en 
        https://github.com/Melii99/Python_class_2022/tree/main/Data/archivo_articulos y los ejemplos 
        de input y output.
        *Para el correcto funcionamiento del programa es necesario tener instalado Biopython
        *Dentro de la función buscar_db se realiza importación de Entrez y Re
        (el email es modificable) 
        *Para el archivo archivo_articulos la ruta escrita en este programa no funcionará
        en otros equipos por lo que deberá modificarse 
INPUT
        Un diccionario con el organismo (key) y una lista de bases de datos y sus 
        IDs correspondientes (values)
OUTPUT
        Un archivo que contiene el organismo, título de los artículos encontrados con los IDs 
        proporcionados y título y abstracts de los artículos correspondientes a las citas de cada
        artículo.
EXAMPLES
        Input
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

        Output
        Organismo: Drosophila melanogaster

        Titulo principal: TI  - Developmental Robustness: The Haltere Case in Drosophila.

        Citas: 

        Título: TI  - Association of HLA-DR1, HLA-DR13, and HLA-DR16 Polymorphisms with ...
        Abstract:
        1:  Association of HLA-DR1, HLA-DR13, and HLA-DR16 Polymorphisms with Systemic ...
        Tingrui Wang, Hong Wang, Lijuan Qiu, Lingling Wu, Huayun Ling, Yu Xue, Ying Zho ...
        J Immunol Res. 2022; 2022: 8140982.  Published online 2022 Apr 16. doi: 10 ...
        PMCID: PMC9034954


        Título: TI  - TNFB Gene Polymorphism in Patients with Systemic Lupus Erythemato ...
        Abstract:
        1:  TNFB Gene Polymorphism in Patients with Systemic Lupus Erythematosus in Kore ...
        Ho-Youn Kim, Sang-Heon Lee, Hyung-In Yang, Sung-Hwan Park, Chul-Soo Cho, Tai-Gyu ...
        Korean J Intern Med. 1995 Jul; 10(2): 130–136.  doi: 10.3904/kjim.1995.10.2 ...
        PMCID: PMC4532044

        ...
GITHUB 
        https://github.com/Melii99/Python_class_2022/tree/main/Tareas/8_NCBI_2.py
       
"""

'''
La función extraer_y_guardar recibe un diccionario con el nombre del organismo y una lista
de bases de datos y IDs correspondientes de los que se realiza una búsqueda sobre cada organismo 
de cada ID en bases de datos que contienen artículos. Se extraen y agregan los títulos de los 
artículos encontrados y los IDs de las citas. Con la lista de IDs se realiza una nueva búsqueda 
de la cual se extraen los títulos y abstracts de las citas y se agregan al archivo. En caso de no 
encontrarse citas, también se escribe en el archivo.
'''
def extraer_y_guardar(diccionario):
  #Realizamos la importación de librerías
  from Bio import Entrez
  import re
  #Proporcionamos un email
  Entrez.email = "mmayen@lcg.unam.mx"
  #Abrimos un archivo en modo escritura 
  path = "C:\\Users\\Melissa\\Downloads\\archivo_articulos"
  archivo = open(path, "w")

  # Recorremos cada organismo del diccionario
  for i in diccionario:
    # ids_list contiene los ids de las bases de datos de cada organismo
    ids_list = diccionario[i]
    # Agregamos el organismo 
    archivo.write('Organismo: ' + i + '\n\n')
    # Recorremos la lista con los ids de diferentes bases de datos
    for db in ids_list:
      # Si encontramos una base de datos que contengan artículos (pubmed,pmc) se obtienen los ids 
      if db == ('pmc' or 'pubmed'):
        ids_articulos = ids_list[ids_list.index(db)+1]
        # Iteramos sobre cada id a buscar
        for id in ids_articulos:
          #Realizamos la búsqueda de los ids (papers) en la base de datos correspondiente
          handle = Entrez.efetch(db=db, id=id, rettype="medline", retmode="text")
          title = handle.read()
          #Realizamos la búsqueda de los links (IDs de los artículos que citan al artículo en la misma base de datos)
          citas = Entrez.read(Entrez.elink(dbfrom=db, db=db, LinkName="pubmed_pmc_refs", from_uid=id))
          # Cerramos el handle
          handle.close()
          # Buscamos el titulo separando las lineas e iterando sobre ellas
          lines = str(title).split('\n')
          for line in lines:
            m = re.search('TI', line)
            if m:
              # Se imprime el título
              #print(m.string)
              archivo.write('Titulo principal: ' + m.string + '\n\n')
          # Si se encuentran citas se guardan los IDs
          if len(citas[0]['LinkSetDb']) > 0:
            cite_ids = [link["Id"] for link in citas[0]['LinkSetDb'][0]['Link']]
            # Imprimimos los IDs de las citas 
            #print(cite_ids)

            archivo.write('Citas: \n\n')
            ### Búsqueda de los títulos y abstracts de las referencias ###
            for c_id in cite_ids:
              #Realizamos la búsqueda de los ids (papers) en la base de datos correspondiente
              handle2 = Entrez.efetch(db=db, id=c_id, rettype="medline", retmode="text")
              #Realizamos la búsqueda de los abstracts de las citas
              handle3 = Entrez.efetch(db=db, id=c_id, rettype="abstract", retmode="text")
              # Guerdamos el handle del que recuperaremos el título de las citas
              title2 = handle2.read()
              # Guardamos el abstract de cada referencia
              abstract = handle3.read()
              # Cerramos los handles
              handle2.close()
              handle3.close()
              # Buscamos el título separando las lineas e iterando sobre ellas
              lines2 = str(title2).split('\n')
              # Se agregan los títulos de las citas 
              for line2 in lines2:
                m2 = re.search('TI', line2)
                if m2:
                  archivo.write('Título: ' + m2.string + '\n')
              # Se agregan los abstracts de las citas
              archivo.write('Abstract:\n' + abstract + '\n')
          else:
            # En caso de no encontrar citas lo especificamos
            archivo.write('No se encontraron citas\n\n')
      archivo.write('\n\n\n')

  return(archivo)        

'''
Ejemplo de búsqueda: se agregan los objetos necesarios para correr una búsqueda

# Objeto que se le pasa a la función
dicc = {'Drosophila melanogaster': ['pmc',
  ['8343187', '7463283', '6366919', '4949521', '2965233', '2813253', '2712966', '2047340', '1852585', '1440717', '1472285', '1461466', '1460734', '1204120'],
  'books',
  ['3312398'],
  'genome',
  ['47']],
 'Caenorhabditis elegans': ['pmc',
  ['9465036', '9441740', '9470895', '9304931', '9217932', '9168093', '8885179', '8869312', '8753122', '8521664', '8487828', '8465798', '8314386', '8143881', '8001149', '7728442', '9057432', '7768255', '7504583', '7786404'],
  'books',
  ['4684026', '4078908', '3069336'],
  'genome',
  ['41']]}

extraer_y_guardar(dicc)

'''
