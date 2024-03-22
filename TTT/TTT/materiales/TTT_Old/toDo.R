#Test de bondad de ajuste basado en sizer map 
#Para evaluar el método elegimos un modelo concreto  Phi(t) es la curva TTT teórica (TRUE-TTT)

#Podemos empezar con una distribución Weibull (forma=3,escala=1)



#H0: Phi es cóncava (TRUE MODEL ES IFR) 
#H1: ---------------------------------
  
# 1. Generamos una muestra de tamaño n=100, 500, 1000 (probamos varios casos distintos)
# 2. Construimos el estimador local cuadrático de la curva TTT basada en los datos: Tenemos estimaciones de Phi, Phi' y Phi''
# 3. Producimos el SiZer map para Phi''(segunda derivada, esto es lo que has hecho hasta ahora). Este mapa le vamos a denotar "EmpiMap" (mapa empírirco. luego simplificamos la notación.)





# 4. Si los datos son generados de un modelo IFR, como es el caso de la W(3,1), el EmpiMap debería ser "significativamente" azul. 

# Bajo la H0 el mapa debería ser entero azul, este es nuestro "TrueMap". 

# Podemos traducir el contraste de arriba al siguiente en lenguaje de SiZer map:

# H0: EmpiMap=TrueMap  
## H1: ---------------

# 4. Error tipo I: Prob{Rechazar H0 |H0}
# Podemos estimar esta cantidad mediante la proporción de pixels no azules en el EmpiMap.

# Repetimos el experimento digamos R=1000 veces. El test es 'bueno' si el error tipo I se estima muy pequeño en gran cantidad de ocasiones entre las R réplicas, de modo que podemos hacer un summary de los errores tipo I obtenidos mediante el punto 4 a lo largo de las R repeticiones.

# 5. Potencia del test: Esto es 1- Prob{Aceptar H0| H1} Este problema es más difícil que el anterior porque la hipótesis H1 consiste en negar H0, es decir que bajo H1 el TrueMap es un mapa que no es azul entero. Y esto se puede dar de muy diferentes maneras. Tenemos que darle alguna vuelta más a este punto. No lo hacemos de momento.


# 6. p-valor: Para una muestra concreta, se trata de calcular la probabilidad de que bajo la hipótesis nula obtengamos una muestra peor que la que tenemos.  Este problema lo tenemos que (podemos) implementar basándonos en unos datos reales

# 6.1 Creamos el EmpiMap basado en los datos y lo comparamos con el TrueMap. Calculamos la proporción de pixels que no coinciden con lo que queremos.
# 6.2 Generamos una muestra bootstrap con reemplazamiento a partir de los datos
# 6.3 Construimos el EmpiMap* bootstrap y lo comparamos con el TrueMap igual que en el punto 6.1
# 6.4 Repetimos los puntos 6.2-6-3 R=1000 veces
# El p-valor lo podemos definir como la proporción de EmpiMap* (bootstrap) que dieron un resultado peor que el de EmpiMap basado en los datos.

# Esto sería más o menos, espero que funcione bien en la práctica. El miércoles vemos la comparación con los test clásicos que te comenté ayer y también los test sizer para el análisis de otras propiedades de envejecimiento. 
