## Programa para filtrar el fichero de familias (global y global con panel) segun el MAF de poblacion española
## (BIER) y X1000G < 0.001
## a030_filtrado_MAF.r
## 2016-02-02 julenmendieta92@gmail.com
## Modificado: 2016-02-02
# 
# date ()
# Sys.info ()[c("nodename", "user")]
# commandArgs ()
# rm (list = ls ())
# R.version.string ##"R version 3.2.2 (2015-08-14)"
# 
# try (source (".job.r")); try (.job)
function(patron_familia) {
  setwd (file.path(.job$dir$proces))
  
  # ############################################################################################################################
  # ## Tienes que indicar por que letra/letras empiezan los nombres de las familias, para que el programa
  # ## no se confunda con otros ficheros. Ej. Si empezara por S (en mayuscula), seria "^S".
  # patron_familia = "^S"
  # ############################################################################################################################
  
  filtro_maf = c("X1000G", "BIER")
  final = c("_global", "_panel_global")
  
  ficheros = list.files()
  herencias = grep("dominante|novo|recesivo", ficheros, value = TRUE)
  
  for (her in herencias) {
    dir_herencia = paste(her, "/", sep = "")
    familias = list.files(dir_herencia)
    familias = familias[grep(patron_familia, familias)]
    
    for (fa in familias) {
      for (fin in final) {
        fich = paste(fa, fin, sep = "")
        dir = paste(dir_herencia, fa, "/", fich, sep = "")
        tabla = read.table(dir, header = TRUE, sep = ";")
        for (fil in filtro_maf) {
          # Primero guardamos las frecuencias, les quitamos el alelo y las guardamos como numerico (las missing pasan a NA)
          frec_temp = as.character(tabla[,fil])
          head(frec_temp)
          frec_temp = as.numeric(gsub("[ ([:alpha:])]", "", unlist(frec_temp)))
          head(frec_temp)
          # Guardamos las missing y las que tienen una frecuencia menor que 0.01
          pos = is.na(frec_temp) | frec_temp < 0.01
          table(pos)
          # pos hace referencia a las filas de la tabla con missing un MAF < 0.01
          tabla2 = tabla[pos,]
    #       # Si antes de guardar queremos reordenar la tabla para que los missing esten en la parte de abajo
    #       # Para ello hay cambiar los NA a un valor imposible, porque 0 ya hay
    #       frec_temp = as.character(tabla2[,fil])
    #       frec_temp = as.numeric(gsub("[ ([:alpha:])]", "", unlist(frec_temp)))
    #       frec_temp[is.na(frec_temp)] <- 999999
    #       tabla2 = tabla2[ order(frec_temp), ]
          
          dir2 = paste(dir_herencia, fa, "/", fich, "_MAF_", fil, sep = "")
          write.csv2(tabla2, dir2, row.names = FALSE, quote = FALSE)
        }
        
      }
    }
  }
}
###EXIT
# warnings ()
# sessionInfo ()
# q ("no")