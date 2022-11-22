
"""
NAME
       9_Numpy.py
VERSION
        [1.0]
AUTHOR
        Melissa Mayén Quiroz
DESCRIPTION
        Este programa consta de ejercicios y ejemplos del uso de Numpy. Se imprimirán algunas
        operaciones y variables a pantalla.

        *Para el correcto funcionamiento del programa es necesario tener instalado Biopython
        *Dentro del programa importamos algunas librerías:
        - Import Numpy as np
        *Para los archivos Biorreactor.csv, prod_gml.csv, ind_cost.csv y Produccion_ej.csv la 
        ruta escrita en este programa no funcionará en otros equipos por lo que deberá modificarse 
INPUT
       
OUTPUT
       
EXAMPLES
       
GITHUB 
        https://github.com/Melii99/Python_class_2022/tree/main/Tareas/9_Numpy.py
"""

##### Ejercicios Clase Numpy #####
import numpy as np

### Array de 1 dimensión ###
array_1D = np.array([1,2,3])
array_1D
# Shape, len y ndim de array_1D
print(array_1D.shape)
print(len(array_1D))
print(array_1D.ndim)

### Array de 2 dimensiones ###
array_2D = np.array([[1,2,3],
                     (2,3,4)
                    ])
array_2D
# Shape, len y ndim de array_2D
print(array_2D.shape)
print(len(array_2D))
print(array_2D.ndim)

### Otro array de 2 dimensiones ###
array_2D = np.array([ ['F','C'],
                     ('i','o'),
                     ('l','l'),
                     ('a','m'),
                     ('s','s')])
array_2D
# Shape, len y ndim del nuevo array_2D
print(array_2D.shape)
print(len(array_2D))
print(array_2D.ndim)

### Array de 3 dimensiones ###
array_3D = np.array([  
                     [ [1, 2, 2],
                       [3, 4, 4] ]
                    ])
array_3D
# Shape, len y ndim del nuevo array_3D
print(array_3D.shape)
print(len(array_3D))
print(array_3D.ndim)

### Podemos crear arrays a partir de archivos (csv) ###
ecoli_m_b = np.genfromtxt('Biorreactor.csv', delimiter=',')
# Para realizar operaciones no es necesario iterar sobre cada elemento del array!!! #

# Multiplicación 
ecoli_m_b * 0.39
# Multiplicación y asignación
ecoli_m_b *= 0.39
# División
originales = ecoli_m_b/0.39


##### Ejercicio 1 #####
'''
Parte 1
Al inducir 4 genes de producción a diferentes temperaturas se obtuvieron las siguientes producciones del metabolito de interés en g/mL:

30 °C	35 °C
Gen 1	0.005	0.003
Gen 2	0.011	0.007
Gen 3	0.004	0.009
Gen 4	0.002	0.006

Pero necesitamos las cantidades en g/L
xg1mL×1000mL1L
xg1×10001L
1000xg1L

Conviertan los valores a g/L
'''
# Creamos el array (de genes y temperatura) a partir del archivo prod_gml.csv
prod_g_mL = np.genfromtxt('prod_gml.csv', delimiter=',')

# Multiplicamos por 1000 el array prod_g_mL para cambiar de g/ml a g/L
prod_g_L = prod_g_mL * 1000

'''
Parte 2
Cada gen tiene un inductor diferente. Estos inductores tienen distintos costos en el mercado

Costo de inductor (g)
Ind Gen 1	3.5
Ind Gen 2	5
Ind Gen 3	7
Ind Gen 4	4.3

Obtener los costos de inductor en cada temperatura
- 30°C necesito agregar 1.75 unidades del inductor
- 35°C necesito agregar 0.8 unidades del inductor
'''
# Creamos el array de los costos de los inductore a partir del archivo idx_cost 
idx_cost = np.genfromtxt('ind_cost.csv', delimiter=',')
# Obtenemos los costos del inductor a 30º
idx_cost_30 = idx_cost * 1.75
# Obtenemos los costos del inductor a 35º
idx_cost_35 = idx_cost * 0.8


### Ejemplo producción (suma, resta) ###
'''
Tengo 3 bacterias produciendo 2 metabolitos de interés biotecnológico.

Si al final de mi producción en un biorreactor con 50 L tengo las siguientes cantidades de metabolito en g/L:
Metabolito A	Metabolito B
Bacteria 1	16	14
Bacteria 2	12	9
Bacteria 3	8	15

- ¿Cuál sería el total de mi producción si tengo 2 biorreactores?
- Del total del metabolito, necesito 5 unidades para hacer pruebas de pureza. ¿Qué total me queda?

'''
# Creamos un array a partir del archivo produccion_ej.csv
produccion = np.genfromtxt('Produccion_ej.csv', delimiter=',')
# Suma dos arrays
produccion + produccion
# También se puede multiplica el array x 2
produccion*2

# Por metabolito
# Numpy sum y elección de axis 
np.sum(produccion*2, axis=0 )

# Por baceria
# Numpy sum y elección de axis 
np.sum(produccion*2, axis=1)

# Numpy sum y elección de axis 
produ_met = np.sum(produccion*2, axis=0)
# Resta 
prod_raw = produ_met-5

### Es posible accesar por índices ###
extra = produccion[:,0] - produccion[:,1]

##### Ejercicio 2 #####
'''
Obtener los costos de producción por cada g/L
Ej. El gen 1 a 30°, utiliza 1.75 g de inductor. Esto tiene un costo Y.
Si el gen 1 produce 5 g/L con un costo de Y, ¿Cuánto cuesta producir 1 g/L?
'''
# Obtenemos los costos de cada gen por cada g/L a 30º
cost_gL_30 = idx_cost * 1.75 /prod_g_L[:,0]
# Obtenemos los costos de cada gen por cada g/L a 30º
cost_gL_35 = idx_cost * 0.8 / prod_g_L[:,1]

### Es posible usar operadores booleanos ###
bool_np = np.array([True, False, True]) 
bool_np.dtype
# Acceso con booleanos
extra[extra > 0]

##### Ejercicio 3 #####
'''
¿En algún caso es más barato producir en 30° en lugar de 35°?

Costo por g/L en 30°C	Costo por g/L en 35°C	Diferencia de costos
Gen 1	x1	y1	x1-y1 = +z1
Gen 2	x2	y2	x2-y2 = -z2
Gen 3	x3	y3	x3-y3 = +z3
Gen 4	x4	y4	x4-y4 = +z4

Acceder con un booleano donde x-y < 0
'''
# Creamos un vector booleano donde comprobamos en qué caso es más barato producir a 30º
booleano = idx_cost * 1.75 /prod_g_L[:,0] - idx_cost * 0.8 / prod_g_L[:,1] < 0
# Podemos hacer lo mismo con los arrays cost_gL_30 y cost_gL_35
booleano = cost_gL_30 - cost_gL_35 < 0
# Ingresamos al aray idx_cost con el vector booleano
idx_cost[booleano]