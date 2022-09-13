"""
NAME
       3_Clases.py
VERSION
        [1.0]
AUTHOR
        Melissa Mayén Quiroz
DESCRIPTION
        Este programa crea varias clases (Alumnos, Alumnos_2, Transporte,
        Terrestre, Acuatico y Aereo) con diferentes métodos y atributos para 
        explicar dentro de la documentación su uso y funcionamiento, además de
        añadir algunos ejemplos específicos.
        Para las clases Terrestre, Acuatico y Aereo usamos herencia de la clase 
        Transpote, existiendo polimorfismos entre las mismas.
GITHUB 
        https://github.com/Melii99/python_class/tree/branch1/Tareas/3_Clases.py
       
"""

# Creamos una clase llamada alumno que guardará su calificación #
class Alumnos():
  def __init__ (self,nombre,calif):
    self.nombre = nombre
    self.calif = calif

  def aprobado(self): 
    if(self.calif > 5):
      print("El alumno aprobó") 
    else:
      print("El alumno reprobó")

# Creamos un ejemplo para una alumna llamada Ana cuya calificación es 8
ana = Alumnos('Ana',8)
# Aplicamos el método 'aprobado' 
ana.aprobado()


# Podemos seguir la misma lógica para calificaciones de varias materias #
class Alumnos_2():
  def __init__(self,nombre,calif1,calif2,calif3):
    self.nombre = nombre
    self.calif1 = calif1
    self.calif2 = calif2
    self.calif3 =calif3

  def promedio(self):
    self.promedio = (self.calif1+self.calif2+self.calif3)/3

# Creamos un ejemplo para un alumno llamado pedro con calificaciónes 10, 9 y 7
pedro = Alumnos_2('Pedro',10,9,7)
# Aplicamos el método 'promedio'
pedro.promedio()
# Ingresamos al atributo 'promedio'
pedro.promedio


# Podemos hacer una clasificación jerárquica de medios de transporte #

# Todos son transportes pero la forma en que se trasladadan es distinta
class Transporte():
  viajar = True
  def __init__(self,modelo,color,tamaño,uso):
    self.modelo = modelo
    self.color = color
    self.tamaño = tamaño
    self.uso = uso

# A partir de la clase transporte creamos otras (herencia)
# Podemos cambiar los valores de los atributos (overriding)
# Así tenemos polimorfismos entre las nuevas clases que creamos

class Terrestre(Transporte):
  terrestre = True  
  acuatico = False
  aereo = False

class Acuatico(Transporte):
  terrestre = False
  acuatico = True
  aereo = False

class Aereo(Transporte):
  terrestre = False
  acuatico = False
  aereo = True

# Creamos un ejemplo para un avión
avion = Aereo('boeing','blanco','grande','carga')
# Ingresamos al atributo 'modelo'
avion.modelo
# Creamos un ejemplo para un tren
tren = Terrestre('suburbano','gris','grande','pasajeros')
# Ingresamos al atributo 'uso'
tren.uso