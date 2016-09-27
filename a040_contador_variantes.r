## Programa para modificar el fichero que lleva las variantes en cada paso, y ampliarlo
## a040_contador_variantes.r
## 2016-02-02 julenmendieta92@gmail.com
## Modificado: 2016-02-02

# date ()
# Sys.info ()[c("nodename", "user")]
# commandArgs ()
# rm (list = ls ())
# R.version.string ##"R version 3.2.2 (2015-08-14)"
# 
# try (source (".job.r")); try (.job)
function(patron_familia) {
  setwd (file.path(.job$dir$proces))
  
#   ############################################################################################################################
#   ## Tienes que indicar por que letra/letras empiezan los nombres de las familias, para que el programa
#   ## no se confunda con otros ficheros. Ej. Si empezara por S (en mayuscula), seria "^S".
#   patron_familia = "^S"
#   ############################################################################################################################
  
  ficheros = list.files()
  herencias = grep("dominante|novo|recesivo", ficheros, value = TRUE)
  
  for (her in herencias) {
    dir_her = paste(her, "/nvariantes", sep = "")
    fich_variant = list.files(dir_her)
    fich_variant = fich_variant[grep(patron_familia, fich_variant)]
    
    familias = list.files(her)
    familias = familias[grep(patron_familia, familias)]
    
    # Primero tomamos los ficheros con el conteo de variantes generado, y los agrupamos
    tabla_variantes <- data.frame()
    for (fi in fich_variant) {
      dir_her = paste(her, "/nvariantes/", fi, sep = "")
      tabla = read.table(dir_her, header = TRUE, sep = ";")
      nrow(tabla)
      if (fi == fich_variant[1]) {
        tabla_variantes = tabla
      }
      else {
        tabla_variantes = rbind(tabla_variantes, tabla)
      }
    }
    row.names(tabla_variantes) = tabla_variantes[,1]
    
    # Ahora añadimos las comparaciones que faltan, con los MAF
    MAF = c("_MAF_BIER", "_MAF_X1000G")
    final = c("_global", "_panel_global")
    
    for (fa in familias) {
      for (fin in final) {
        for (m in MAF) {
          fich = paste(fa, fin, m, sep = "")
          dir = paste(her, "/", fa, "/", fich, sep = "")
          tabla = read.table(dir, header = TRUE, sep = ";")
          nomb = paste("variantes", fin, gsub("_MAF", x=m, ""), sep = "")
          tabla_variantes[fa, nomb] = nrow(tabla)
        }
      }
    }
    
    colnames(tabla_variantes)   
    # Ponemos todo en un orden mas logico
    tabla_variantes = tabla_variantes[,c(1,2,3,5,6,4,7,8)]
    
    # Y guardamos
    salida = paste(her, "/nvariantes/nvariantes", sep = "")
    write.csv2(tabla_variantes, salida, row.names = FALSE, quote = FALSE)
  }
}

###EXIT
# warnings ()
# sessionInfo ()
# q ("no")
