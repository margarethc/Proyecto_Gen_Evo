# Proyecto_Gen_Evo
Proyecto del curso de Genómica Evolutiva
# Análisis de conservación evolutiva para secuencias proteicas
## Objetivo
El objetivo de este trabajo es estimar la conservación evolutiva de cada residuo de una proteínas y mapear dichos valores en una estructura 3D. 
Esto permitirá identificar residuos funcionales importantes, localizar sitios activos, diferenciar regiones conservadas vs. variables y comparar modelos estructurales.

## Abstract

## Variables
### Inputs
- Secuencia de aminoácidos
- Estructura 3D de la proteína
- MSA (a partir de la secuencia)
- Árbol filogenético (a partir del MSA)
### Outputs
- Estructura 3D con anotaciones de conservación
- Secuencia con anotaciones de conservación

## Flujo
### 1. Búsqueda de secuencia homólogas
Las secuencias homólogas a la secuencia input se obtienen a partir de un hmmsearch en la base de datos de Pfam. 
### 2. Alineamiento de secuencia múltiple
A partir de la secuencia y los mejores homólogos (se podría insertar un propio archivo de homología) se realiza el alineamiento con T-coffee o mafft. 
### 3. Obtención de árbol filogenético
A partir del MSA generado
### 4. Cálculo de conservación utilizando Rate4Site
### 5. Visualización en 3D 
Las estructuras 3D se anotan con los resultados de Rate4Site utilizando Jmol. 

