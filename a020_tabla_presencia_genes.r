## Programa para generar un archivo con la presencia de cada gen en las familias
## a020_tabla_presencia_genes.r
## 2015-03-26 julenmendieta92@gmail.com
## Modificado: 2015-04-23

# date ()
# Sys.info ()[c("nodename", "user")]
# commandArgs ()
# rm (list = ls ())
# R.version.string ##"R version 3.2.2 (2015-08-14)"
# 
# try (source (".job.r")); try (.job)

function(patron_familia) {
#   ############################################################################################################################
#   ## Tienes que indicar por que letra/letras empiezan los nombres de las familias, para que el programa
#   ## no se confunda con otros ficheros. Ej. Si empezara por S (en mayuscula), seria "^S".
#   patron_familia = "^S"
#   ############################################################################################################################
  
  # Primero cargamos la lista de genes de interes
  setwd (file.path (.job$dir$raw))
  ficheros = list.files()
  genes = ficheros[grep("^genes.txt", ficheros)]
  if (length(genes) > 0) {
    genes_filtrado = read.table(genes, header = FALSE, sep = "\t")
    genes_filtrado = as.character(genes_filtrado[,1])
    genes = genes_filtrado
  }
  setwd (file.path(.job$dir$proces))
  
  
  ficheros = list.files()
  herencias = grep("dominante|novo|recesivo", ficheros, value = TRUE)
  for (her in herencias) {
    dir_herencia = paste(her, "/", sep = "")
    familias = list.files(dir_herencia)
    familias = familias[grep(patron_familia, familias)]
    tabla_genes <- data.frame(S1010 = numeric(0), S1315 = numeric(0), S1349 = numeric(0), S1651 = numeric(0), 
                              nfam = numeric(0), panel = numeric(0))
    for (fa in familias) {
      fich = paste(fa, "_global", sep = "")
      dir = paste(dir_herencia, fa, "/", fich, sep = "")
      tabla = read.table(dir, header = TRUE, sep = ";")
      genes_temp = tabla[,"Gene"]
      head(genes_temp)
      # Como de momento no se que hacer con los missing genes, los elimino de la lista, y también los genes duplicados
      table(genes_temp == ".")
      genes_temp = genes_temp[genes_temp != "."]
      genes_temp = as.character(genes_temp)
      # Separamos los genes que estan solapados y guardamos una lista con todos los genes que aparecen sin repetir
      length(genes_temp)
      genes_temp = strsplit(genes_temp, ",")
      genes_temp = unlist(genes_temp, recursive = TRUE, use.names = TRUE)
      length(genes_temp)
      genes_temp = unique(genes_temp)
      # Ahora lo guardamos en la tabla
      # Puede ser que no tengamos genes en alguna de las herencias, así que para que el programa pueda continuar habra que dejar alternativas
      if (length(genes_temp) > 0) {
        tabla_genes[genes_temp, fa] = 1
      }
      else {
        print_salida = paste ("Con la herencia ", her, " y para la familia ", fa, ", no tenemos variantes", sep = "")
        print (print_salida)
      }
    }
    
    # Ahora Eliminamos todo lo que sea NA y lo pasamos a 0
    tabla_genes[is.na(tabla_genes)] = 0
    # Esto estaria bien ahcerlo con un apply cuando sepa usarlos
    for (i in 1:dim(tabla_genes)[1]) {
      tabla_genes[i, "nfam"] = rowSums(tabla_genes[i,])
    }
    # a=lapply(tabla_genes, function(x) tabla_genes[x, "nfam"] = rowSums(tabla_genes[x,]))
    # tabla_genes = tabla_genes[order(tabla_genes[,"nfam"], decreasing=TRUE),]
    
    # Añadimos la info de si el gen está presente en el panel
    pos = rownames(tabla_genes) %in% genes
    table(pos)
    tabla_genes[pos, "panel"] = 1
    tabla_genes = tabla_genes[ order(-tabla_genes[,"nfam"], -tabla_genes[,"panel"]), ]
    
    # Ahora ordenamos por la ultima columna y guardamos la salida en un fichero
    salida = paste(dir_herencia, "Tabla_presencia_genes", sep = "")
    write.csv2(tabla_genes, salida, row.names = TRUE, quote = FALSE)
    
    ## Ahora aprovechamos esto para generar los ficheros con genes que tienen variantes de interés en todas las familias
    
    # Primero guardamos los genes
    pos = tabla_genes[, "nfam"] == length(familias)
    genes = rownames(tabla_genes[pos, ])
    
    # Ahora habra que buscarlos familia por familia
    final = c("_global", "_panel_global")
    suppressWarnings(dir.create(file.path(.job$dir$proces, her, "comunes")))
    
    for (fa in familias) {
      for (fin in final) {
        fich = paste(fa, fin, sep = "")
        dir = paste(dir_herencia, fa, "/", fich, sep = "")
        dir2 = paste(dir_herencia, "comunes", "/", fich, "_comun", sep = "")
        #fichero = paste(fa, fin)
        tabla = read.table(dir, header = TRUE, sep = ";")
        # Seleccionamos las posiciones de interes
        matches <- unique(grep(paste("\\b",genes, "\\b", collapse="|", sep=""), 
                               tabla[,"Gene"], value=TRUE))
        pos = tabla[,"Gene"] %in% matches
        table(pos)
        write.csv2(tabla[pos,], dir2, row.names = FALSE, quote = FALSE)
        
      }
    }
  }
}

###EXIT
# warnings ()
# sessionInfo ()
# q ("no")