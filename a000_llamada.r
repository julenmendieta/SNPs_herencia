## Programa para llamar a los otros programas de filtrado
## a000_llamada.r
## 2016-03-04 julenmendieta92@gmail.com
## Modificado: 

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.2.2 (2015-08-14)"
try (source (".job.r")); try (.job)

############################################################################################################
## Primero, que es lo que quieres buscar? "novo", "novo2", "dominante", "recesivo"
## Instrucciones al final de a010_seleccion_3_herencias.r
busqueda = c("novo")
############################################################################################################

############################################################################################################################
## Tienes que indicar por que letra/letras empiezan los nombres de las familias, para que el programa
## no se confunda con otros ficheros. Ej. Si empezara por S (en mayuscula), seria "^S".
patron_familia = "^S"
############################################################################################################################

# Leemos los ficheros de filtrado
filtrado1 <- dget("a010_seleccion_3_herencias.r")
filtrado2 <- dget("a020_tabla_presencia_genes.r")
filtrado3 <- dget("a030_filtrado_MAF.r")
filtrado4 <- dget("a040_contador_variantes.r")

# Lanzamos cada uno de los pasos
filtrado1(busqueda)
filtrado2(patron_familia)
filtrado3(patron_familia)
filtrado4(patron_familia)

##EXIT
warnings ()
sessionInfo ()
q ("no")
