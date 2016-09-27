## Programa para filtrar en busca de un genotipo dominante que lleve a la enfermedad.
## a010_seleccion_3_herencias.r
## 2016-02-02 julenmendieta92@gmail.com
## Modificado: 2016-02-02

## OJO Mira a ver si los nombres de individuos se diferencia por "_" o ".". Yo he tenido en cuenta que en los ficheros se 
## diferencian por "_", y dentro del csv por "."


# date ()
# Sys.info ()[c("nodename", "user")]
# commandArgs ()
# rm (list = ls ())
# R.version.string ##"R version 3.2.2 (2015-08-14)"

function(busqueda) {
  try (source (".job.r")); try (.job)
  setwd (file.path (.job$dir$raw))
  
#   ############################################################################################################
#   ## Primero, que es lo que quieres buscar? "novo", "novo2", "dominante", "recesivo"
#   ## Instrucciones al final de a010_seleccion_3_herencias.r
#   busqueda = c("novo")
#   ############################################################################################################
  
  
  
  # No es que no me fie de ti
  opciones = c("novo", "novo2", "dominante", "recesivo")
  pos = busqueda %in% opciones
  if (length(pos) != sum(pos)) {
    print("Alguna de las opciones no es correcta")
    stop("¡¡ERROR!!")
  }
  # Necesitamos los csv y los ficheros .ped en la carpeta de data raw, y en caso de tener genes candidatos, un
  # fichero de texto con el texto genes, y los genes en una columna 
  ficheros = list.files()
  pedigris = ficheros[grep("*.ped$", ficheros)]
  genes = ficheros[grep("^genes.txt", ficheros)]
  ficheros_inic = ficheros[grep("*.csv$", ficheros)]
  
  # Comprobamos que al menos el numero de ficheros este bien
  if (length(pedigris) != length(ficheros_inic)) {
    stop("El numero de ficheros no coincide")
  }
  
  # Guardamos los genes que se van a utilizar para filtrar en una variable
  # El fichero tiene que ser un txt o csv con los genes en la primera columna
  genes_filtrado = c()
  if (length(genes) > 0) {
    genes_filtrado = read.table(genes, header = FALSE, sep = "\t")
    genes_filtrado = as.character(genes_filtrado[,1])
  }
  
  
  for (fich_inic in ficheros_inic) {
    tabla = read.table(fich_inic, header = TRUE, sep = "\t")
    colnames(tabla)
    
    # Antes de nada guardamos el número de filas que hay para el fichero con numero de variantes
    variantes_csv = nrow(tabla)
    familia = gsub("_bierapp.csv$", "",fich_inic)
    # Tenemos que seleccionar presencia en heterocigosis en enfermos y no presencia en sanos
    # Para ello, primero hay que saber cuales son enfermos y cuales sanos
    ped_familia = paste(familia, ".ped", sep = "")
    pedigri = read.table(ped_familia, header = FALSE, sep = "\t", colClasses = "character")
    
    # Guardamos quien es que, y mas adelante corregiremos diferencias de formato ("-" y ".")
    individuos = as.character(pedigri[,2])
    sanos = as.character(pedigri[pedigri[,6] == 1, 2])
    enfermos = as.character(pedigri[pedigri[,6] == 2, 2])
    # Guardamos los individuos con dos progenitores que sean caso
    # Primero se cogen los que tienen info de los dos progenitores
    pos = pedigri[,3] != 0 & pedigri[,4] != 0
    posp = pedigri[pos, 3] %in% enfermos
    posm = pedigri[pos, 4] %in% enfermos
    progenitores_caso = as.character(pedigri[pos, 2])
    # Despues borramos los que no tienen los dos progenitores enfermos
    wi = 1
    i = 1
    while (wi <= length(posp)) {
      if (!(posp[i] == TRUE & posm[i] == TRUE)) {
        progenitores_caso = progenitores_caso[-i]
      }
      else {
        i = 1 + 1
      }
      wi = wi + 1
      }
    remove(i)
    remove(wi)
    # Nuestra estructura de datos diferencia familias por guiones, pero BierApp lo hace con puntos
    progenitores_caso = gsub("-", ".", progenitores_caso)
    individuos = gsub("-", ".", individuos)
    sanos = gsub("-", ".", sanos)
    enfermos = gsub("-", ".", enfermos)
    # Guardamos listas de individuos con datos de progenitores para novo y novo2
    if ("novo" %in% busqueda | "novo2" %in% busqueda) {
      pos = pedigri[,3] != 0 & pedigri[,4] != 0 & pedigri[,6] == 2
      # Para novo2, con que tengamos datos de los dos progenitores y el individuo este enfermo sirve
      individuos_novo2 = as.character(pedigri[pos, 2])
      # Para novo tenemos que mirar si los progenitores estan sanos, si no no sirven
      individuos_novo = individuos_novo2
      for (i in individuos_novo2) {
        posi = pedigri[,2] == i
        padre = pedigri[posi,3] 
        madre = pedigri[posi,4]
        # Si ninguno de los progenitores esta enfermo podemos plantearla como una enfermedad de novo, y en el caso de q esto
        # no se cumpla, eliminamos de la lista a este individuo
        if (!(pedigri[pedigri[,2] == padre, 6] == 1 & pedigri[pedigri[,2] == madre, 6] == 1)) {
          individuos_novo = individuos_novo[individuos_novo != i]
        }
      }
    }
   
    
    
    ############# Primer filtrado, con controles solo a 0/0, y casos con 1/0 o 0/1
    ## DOMINANTE ##
    if ("dominante" %in% busqueda) {
      tabla_dominante = tabla
      for (enf in enfermos) {
        if (enf %in% progenitores_caso) {
          pos = tabla_dominante[, enf] == "0/1" | tabla_dominante[, enf] == "1/0" | tabla_dominante[, enf] == "1/1"
        }
        else {
          pos = tabla_dominante[, enf] == "0/1" | tabla_dominante[, enf] == "1/0"
        }
        tabla_dominante = tabla_dominante[pos,]
      }
      
      for (san in sanos) {
        pos = tabla_dominante[, san] == "0/0"
        tabla_dominante = tabla_dominante[pos,]
      }
      
  #     unique(tabla_dominante[, enfermos])
  #     unique(tabla_dominante[, sanos])
    }
  
    ## RECESIVO ##
    if ("recesivo" %in% busqueda) {
      tabla_recesivo = tabla
      for (en in enfermos) {
        pos = tabla_recesivo[,en] == "1/1"
        tabla_recesivo = tabla_recesivo[pos,]
      }
      for (san in sanos) {
        pos = tabla_recesivo[,san] == "0/1" | tabla_recesivo[,san] == "1/0"
        tabla_recesivo = tabla_recesivo[pos,]
      }  
      #     unique(tabla_recesivo[, enfermos])
      #     unique(tabla_recesivo[, sanos])
    }
  
    ## DE NOVO ##
    if ("novo" %in% busqueda) {
      tabla_novo = tabla
        
      tabla_novo3 = data.frame()
      tabla_novo4 = data.frame()
      for (indv in individuos_novo) {
        indv1 = gsub("-", "\\.", indv)
        posi = pedigri[,2] == indv
        padre = gsub("-", ".", pedigri[posi,3])
        madre = gsub("-", ".", pedigri[posi,4])
        
        pos = tabla_novo[,indv1] == "0/1" | tabla_novo[,indv1] == "1/0" | tabla_novo[,indv1] == "1/1"
        tabla_novo = tabla_novo[pos,]
        pos = tabla_novo[,padre] == "0/0"
        tabla_novo = tabla_novo[pos,]
        pos = tabla_novo[,madre] == "0/0"
        tabla_novo = tabla_novo[pos,]
        tabla_novo3 = rbind(tabla_novo3, tabla_novo)
        
        # Para el caso de mutación de novo que confiere homocigosidad
        # En este caso el descendiente seria homocigoto
        # Uno de los progenitores tendria la variante de interés
        # El otro no tendria la variante
        tabla_novo2 = tabla  
        pos1 = (tabla_novo2[,padre] == "0/1" | tabla_novo2[,padre] == "1/0") & 
          (tabla_novo2[,madre] == "0/0")
        pos2 = (tabla_novo2[,madre] == "0/1" | tabla_novo2[,madre] == "1/0") & 
          (tabla_novo2[,padre] == "0/0")
        # Guardamos las posiciones que se cumplen en los dos casos
        pos = pos1 + pos2
        pos = pos >= 1
        tabla_novo2 = tabla_novo2[pos,]
        pos = tabla_novo2[,indv1] == "1/1"
        tabla_novo2 = tabla_novo2[pos,]
        tabla_novo4 = rbind(tabla_novo4, tabla_novo2)
      }
      tabla_novo = rbind(tabla_novo4, tabla_novo3)
      tabla_novo = unique(tabla_novo)
      remove(tabla_novo2)
      remove(tabla_novo3)
      remove(tabla_novo4)
    }
    
    ## DE NOVO 2 ##
    if ("novo2" %in% busqueda) {
      tabla_novo2_1 = tabla
      
      tabla_novo2_3 = data.frame()
      tabla_novo2_4 = data.frame()
      for (indv in individuos_novo2) {
        indv1 = gsub("-", "\\.", indv)
        posi = pedigri[,2] == indv
        padre = gsub("-", ".", pedigri[posi,3])
        madre = gsub("-", ".", pedigri[posi,4])
        
        pos = tabla_novo2_1[,indv1] == "0/1" | tabla_novo2_1[,indv1] == "1/0" | tabla_novo2_1[,indv1] == "1/1"
        tabla_novo2_1 = tabla_novo2_1[pos,]
        pos = tabla_novo2_1[,padre] == "0/0"
        tabla_novo2_1 = tabla_novo2_1[pos,]
        pos = tabla_novo2_1[,madre] == "0/0"
        tabla_novo2_1 = tabla_novo2_1[pos,]
        tabla_novo2_3 = rbind(tabla_novo2_3, tabla_novo2_1)
        
        # Para el caso de mutación de novo que confiere homocigosidad
        # En este caso el descendiente seria homocigoto
        # Uno de los progenitores tendria la variante de interés
        # El otro no tendria la variante
        tabla_novo2_2 = tabla  
        pos1 = (tabla_novo2_2[,padre] == "0/1" | tabla_novo2_2[,padre] == "1/0") & 
          (tabla_novo2_2[,madre] == "0/0")
        pos2 = (tabla_novo2_2[,madre] == "0/1" | tabla_novo2_2[,madre] == "1/0") & 
          (tabla_novo2_2[,padre] == "0/0")
        # Guardamos las posiciones que se cumplen en los dos casos
        pos = pos1 + pos2
        pos = pos >= 1
        tabla_novo2_2 = tabla_novo2_2[pos,]
        pos = tabla_novo2_2[,indv1] == "1/1"
        tabla_novo2_2 = tabla_novo2_2[pos,]
        tabla_novo2_4 = rbind(tabla_novo2_4, tabla_novo2_2)
      }
      tabla_novo2_1 = rbind(tabla_novo2_4, tabla_novo2_3)
      tabla_novo2_1 = unique(tabla_novo2_1)
      remove(tabla_novo2_2)
      remove(tabla_novo2_3)
      remove(tabla_novo2_4)
    }
  
    for ( bus in busqueda) {
      if (bus == "dominante") {
        tabla2 = tabla_dominante
      }
      if (bus == "recesivo") {
        tabla2 = tabla_recesivo
      }
      if (bus == "novo") {
        tabla2 = tabla_novo
      }
      if (bus == "novo2") {
        tabla2 = tabla_novo2_1
      }
    
      # guardamos el fichero con esto, y si ya esta creada la carpeta evitamos la aparicion de warnings
      suppressWarnings(dir.create(file.path(.job$dir$proces, bus)))
      suppressWarnings(dir.create(file.path(.job$dir$proces, bus, familia)))
      
      salida = paste(bus, "/", familia, "/", familia, sep = "")
      write.csv2(tabla2, file.path (.job$dir$proces, salida), row.names = FALSE, quote = FALSE)
      
      #Buscamos los genes de interes
      matches <- unique (grep(paste("\\b",genes_filtrado, "\\b", collapse="|", sep=""), 
                              tabla2[,"Gene"], value=TRUE))
      pos = tabla2[,"Gene"] %in% matches
      table(pos)
      tabla2 = tabla2[pos,]
      salida = paste(bus, "/", familia, "/", familia, "_panel", sep = "")
      write.csv2(tabla2, file.path (.job$dir$proces, salida), row.names = FALSE, quote = FALSE)
    }
    
    
    ########### Parte dos. Filtramos por diferentes genotipos.
    # Lo primero, vemos que hay
    table(tabla[, individuos])
    
    # Primero guardamos una lista con todas las variantes
    variantes = c()
    for (indv in individuos) {
      variantes = append(variantes, names(table(tabla[,indv])))
    }
    variantes = unique(variantes)
    variantes
    
    # Vamos a generar todas las posibles combinaciones de interes
    variantes = unique(unlist(strsplit(variantes, split = "/")))
    genotipos = as.numeric(variantes[variantes != "."])
    
    #Generamos todas las posibles combinaciones con este genotipo
    combinaciones_tot = c()
    for (e in genotipos) {
      for (a in genotipos) {
        genot = paste(e, "/", a, sep = "")
        combinaciones_tot = append(combinaciones_tot, genot)
      }
    }
                
    #########################
    tabla_global <- as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla)))
    tabla_global_panel <- as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla)))
  
    if ("dominante" %in% busqueda) {
      tabla_global_panel_dominante = as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla)))
      tabla_global_dominante = as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla)))
    }
    if ("recesivo" %in% busqueda) {
      tabla_global_panel_recesivo = as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla)))
      tabla_global_recesivo = as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla)))
    }
    if ("novo" %in% busqueda) {
      tabla_global_panel_novo = as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla)))
      tabla_global_novo = as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla))) 
    }
    if ("novo2" %in% busqueda) {
      tabla_global_panel_novo2 = as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla)))
      tabla_global_novo2 = as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla)))
    }
    ###########################
    for (i in genotipos) {
  #     combinaciones_san = c()
  #     combinaciones_enf = c()
  #     tabla2 = tabla
      # Suponemos que la variante 0 no va a ser la causante de la enfermedad
      if (i != 0) {
  #       # generamos las combinaciones de sanos
  #       for (e in genotipos) {
  #         if (e != i) {
  #           for (a in genotipos) {
  #             if (a != i) {
  #               # En los controles, ante herencia dominante no se contempla que el alelo de la "enfermedad"(i) este presente
  #               genot = paste(e, "/", a, sep = "")
  #               combinaciones_san = append(combinaciones_san, genot)
  #             }
  #           }
  #           # Generamos las combinaciones de enfermos
  #           # En enfermos, el alelo de la enfermedad tiene que estar al menos presente con menos de 2 parentales enfermos
  #           genot = paste(i, "/", e, sep = "")
  #           combinaciones_enf = append(combinaciones_enf, genot)
  #           genot = paste(e, "/", i, sep = "")
  #           combinaciones_enf = append(combinaciones_enf, genot)
  #         }
  #       }
  #       combinaciones_san = append(combinaciones_san, "./.")
  #       combinaciones_enf = append(combinaciones_enf, "./.")
  #       # Cabe la posibilidad de tener dos parentales enfermos, así que vamos a ñadir las posibles combinaciones
  #       combinaciones_doble_enf = paste(i, "/", i, sep = "")
  #       combinaciones_doble_enf = append(combinaciones_doble_enf, combinaciones_enf)
  #       #Si se diera el caso de que se pueda plantear una mutación de novo
  #       combinaciones_novo = paste(i, "/", i, sep = "")
  #       combinaciones_novo = append(combinaciones_novo, combinaciones_enf)
          
  #       print("El alelo enfermo es: ")
  #       print(i)
  #       print(combinaciones_san)
  #       print(combinaciones_enf)
        
        
        # Guardamos en variable los que cumplan las combinaciones de caso y control, y luego generamos el documento
        
        ## DE NOVO ##
        if ("novo" %in% busqueda) {
          tabla_novo = tabla
          pos_tot1 = grep(i, combinaciones_tot)
          combinaciones_san_novo =  combinaciones_tot[-pos_tot1]
          combinaciones_enf_novo = combinaciones_tot[pos_tot1]
          
          # Aqui hay q añadir la muy remota posibilidad de un progenitor heterocigoto y el hijo homocigoto
          # En este caso el descendiente seria homocigoto
          combinaciones_enf_novo2 = paste(i, "/", i, sep = "")
          # Uno de los progenitores tendria la variante de interés
          combinaciones_san_novo_1 = combinaciones_enf_novo
          # El otro no tendria la variante
          combinaciones_san_novo_2 = combinaciones_san_novo
          
          tabla_novo3 = data.frame()
          tabla_novo4 = data.frame()
          for (indv in individuos_novo) {
            indv1 = gsub("-", "\\.", indv)
            posi = pedigri[,2] == indv
            padre = gsub("-", ".", pedigri[posi,3])
            madre = gsub("-", ".", pedigri[posi,4])
            
            pos = tabla_novo[,indv1] %in% combinaciones_enf_novo
            tabla_novo = tabla_novo[pos,]
            pos = tabla_novo[,padre] %in% combinaciones_san_novo
            tabla_novo = tabla_novo[pos,]
            pos = tabla_novo[,madre] %in% combinaciones_san_novo
            tabla_novo = tabla_novo[pos,]
            tabla_novo3 = rbind(tabla_novo3, tabla_novo)
            
            # Para el caso de mutación de novo que confiere homocigosidad
            tabla_novo2 = tabla  
            pos1 = tabla_novo2[,padre] %in% combinaciones_san_novo_1 & tabla_novo2[,madre] %in% combinaciones_san_novo_2
            pos2 = tabla_novo2[,madre] %in% combinaciones_san_novo_1 & tabla_novo2[,padre] %in% combinaciones_san_novo_2
            # Guardamos las posiciones que se cumplen en los dos casos
            pos = pos1 + pos2
            pos = pos >= 1
            tabla_novo2 = tabla_novo2[pos,]
            pos = tabla_novo2[,indv1] %in% combinaciones_enf_novo2
            tabla_novo2 = tabla_novo2[pos,]
            tabla_novo4 = rbind(tabla_novo4, tabla_novo2)
          }
          tabla_novo = rbind(tabla_novo4, tabla_novo3)
          tabla_novo = unique(tabla_novo)
          remove(tabla_novo2)
          remove(tabla_novo3)
          remove(tabla_novo4)
        }
        
        ## DE NOVO 2##
        if ("novo2" %in% busqueda) {
          tabla_novo2_1 = tabla
          pos_tot1 = grep(i, combinaciones_tot)
          combinaciones_san_novo =  combinaciones_tot[-pos_tot1]
          combinaciones_enf_novo = combinaciones_tot[pos_tot1]
          
          # Aqui hay q añadir la muy remota posibilidad de un progenitor heterocigoto y el hijo homocigoto
          # En este caso el descendiente seria homocigoto
          combinaciones_enf_novo2 = paste(i, "/", i, sep = "")
          # Uno de los progenitores tendria la variante de interés
          combinaciones_san_novo_1 = combinaciones_enf_novo
          # El otro no tendria la variante
          combinaciones_san_novo_2 = combinaciones_san_novo
          
          tabla_novo2_3 = data.frame()
          tabla_novo2_4 = data.frame()
          for (indv in individuos_novo2) {
            indv1 = gsub("-", "\\.", indv)
            posi = pedigri[,2] == indv
            padre = gsub("-", ".", pedigri[posi,3])
            madre = gsub("-", ".", pedigri[posi,4])
            
            pos = tabla_novo2_1[,indv1] %in% combinaciones_enf_novo
            tabla_novo2_1 = tabla_novo2_1[pos,]
            pos = tabla_novo2_1[,padre] %in% combinaciones_san_novo
            tabla_novo2_1 = tabla_novo2_1[pos,]
            pos = tabla_novo2_1[,madre] %in% combinaciones_san_novo
            tabla_novo2_1 = tabla_novo2_1[pos,]
            tabla_novo2_3 = rbind(tabla_novo2_3, tabla_novo2_1)
            
            # Para el caso de mutación de novo que confiere homocigosidad
            tabla_novo2_2 = tabla  
            pos1 = tabla_novo2_2[,padre] %in% combinaciones_san_novo_1 & tabla_novo2_2[,madre] %in% combinaciones_san_novo_2
            pos2 = tabla_novo2_2[,madre] %in% combinaciones_san_novo_1 & tabla_novo2_2[,padre] %in% combinaciones_san_novo_2
            # Guardamos las posiciones que se cumplen en los dos casos
            pos = pos1 + pos2
            pos = pos >= 1
            tabla_novo2_2 = tabla_novo2_2[pos,]
            pos = tabla_novo2_2[,indv1] %in% combinaciones_enf_novo2
            tabla_novo2_2 = tabla_novo2_2[pos,]
            tabla_novo2_4 = rbind(tabla_novo2_4, tabla_novo2_2)
          }
          tabla_novo2_1 = rbind(tabla_novo2_4, tabla_novo2_3)
          tabla_novo2_1 = unique(tabla_novo2_1)
          remove(tabla_novo2_2)
          remove(tabla_novo2_3)
          remove(tabla_novo2_4)
        }
        
        ## RECESIVO ##
        if ("recesivo" %in% busqueda) {
          tabla_recesivo = tabla
          pos_tot1 = grep(i, combinaciones_tot)
          combinaciones_enf_rec = paste(i, "/", i, sep = "")
          combinaciones_san_rec =  combinaciones_tot[combinaciones_tot != combinaciones_enf_rec]
          combinaciones_enf_rec = append(combinaciones_enf_rec, "./.")
          combinaciones_san_rec = append(combinaciones_san_rec, "./.")
          
          for (en in enfermos) {
            pos = tabla_recesivo[,en] %in% combinaciones_enf_rec
            tabla_recesivo = tabla_recesivo[pos,]
          }
          for (san in sanos) {
            pos = tabla_recesivo[,san] %in% combinaciones_san_rec
            tabla_recesivo = tabla_recesivo[pos,]
          }
        }
        
        ## DOMINANTE ##
        if ("dominante" %in% busqueda) {
          tabla_dominante = tabla
          pos_tot1 = grep(i, combinaciones_tot)
          combinaciones_san_dom =  combinaciones_tot[-pos_tot1]
          combinaciones_doble_enf_dom = combinaciones_tot[pos_tot1]
          combinaciones_enf_dom = combinaciones_doble_enf_dom[combinaciones_doble_enf_dom != paste(i, "/", i, sep = "")]
          combinaciones_san_dom = append(combinaciones_san_dom, "./.")
          combinaciones_doble_enf_dom = append(combinaciones_doble_enf_dom, "./.")
          combinaciones_enf_dom = append(combinaciones_enf_dom, "./.")
          
          for (en in enfermos) {
            # Si un individuo tiene dos progenitores enfermos puede tener homocigosis para el mutado
            if (en %in% progenitores_caso) {
              pos = tabla_dominante[,en] %in% combinaciones_doble_enf_dom
            }
            else {
              pos = tabla_dominante[,en] %in% combinaciones_enf_dom
            }
            tabla_dominante = tabla_dominante[pos,]
          }
          for (san in sanos) {
            pos = tabla_dominante[,san] %in% combinaciones_san_dom
            tabla_dominante = tabla_dominante[pos,]
          }
        }
              
        print(unique(tabla2[,enfermos]))
        print(unique(tabla2[,sanos]))
        
        # Solo permitimos un alelo missing en los enfermos por cada variante, y los que sea en sanos
        # Vamos a crear un vector con las posiciones donde hay un missing, si esa posicion se repite es que hay que quitarla
        for ( bus in busqueda) {
          if (bus == "dominante") {
            tabla2 = tabla_dominante
          }
          if (bus == "recesivo") {
            tabla2 = tabla_recesivo
          }
          if (bus == "novo") {
            tabla2 = tabla_novo
          }
          if (bus == "novo2") {
            tabla2 = tabla_novo2_1
          }
          
          pos_mis = c()
          for (en in enfermos) {
            pos = grep("./.", tabla2[, en], fixed = TRUE)
            pos_mis = append(pos_mis, pos)
          }
          pos_del = unique(pos_mis[duplicated(pos_mis)])
          print(pos_del)
          # Guardamos en pos_missing los casos que solo tienen un missing
          pos_missing = pos_mis[!pos_mis %in% pos_del]
          # Con controles lo mismo, pero permitiendo 2
          pos_mis2 = c()
          for (san in sanos) {
            pos = grep("./.", tabla2[, san], fixed = TRUE)
            pos_mis2 = append(pos_mis2, pos)
          }
          pos_mis2 = pos_mis2[!pos_mis2 %in% pos_del]
          print(pos_mis2)
          # Añadimos a pos_missing los controles con missing
          pos_missing = append(pos_missing, unique(pos_mis2))
          # guardamos el fichero con esto
          nom_fichero = paste(bus, "/", familia, "/", familia, "_genot_causante_", i, "_missing", sep = "")
          write.csv2(tabla2[pos_missing,], file.path (.job$dir$proces, nom_fichero), row.names = FALSE, 
                     quote = FALSE)
          
          # y la añadimas a la tabla global
          tabla_global <- rbind(tabla_global,tabla2[pos_missing,])
          
          #Guardamos el fichero tras filtrar con los genes del panel
          if (length(genes) > 0) {
            matches <- unique (grep(paste("\\b",genes_filtrado, "\\b", collapse="|", sep=""), 
                                    tabla2[pos_missing,"Gene"], value=TRUE))
            pos = tabla2[pos_missing,"Gene"] %in% matches
            print(table(pos))
            tabla3 = tabla2[pos,]
            nom_fichero = paste(bus, "/", familia, "/", familia, "_genot_causante_", i, "_missing_panel", sep = "")
            write.csv2(tabla3, file.path (.job$dir$proces, nom_fichero), row.names = FALSE, quote = FALSE)
          }
          
          ####################################################
          # y la añadimas a la tabla global
          tabla_global_panel <- rbind(tabla_global_panel,tabla3)
          
          # Eliminamos los que tienen missing y guardamos normal
          pos = unique(sort(append(pos_mis, pos_mis2), decreasing = TRUE))
          if (length(pos) > 0) {
            tabla2 = tabla2[-c(pos),]
          }
          
          # guardamos el fichero sin missings
          nom_fichero = paste(bus, "/", familia, "/", familia, "_genot_causante_", i, sep = "")
          write.csv2(tabla2, file.path (.job$dir$proces, nom_fichero), row.names = FALSE, quote = FALSE)
          
          # y la añadimas a la tabla global
          tabla_global <- rbind(tabla_global, tabla2)
          
          #Guardamos el fichero sin missings tras filtrar con los genes del panel
          if (length(genes) > 0) {
            matches <- unique (grep(paste("\\b",genes_filtrado, "\\b", collapse="|", sep=""), 
                                    tabla2[,"Gene"], value=TRUE))
            pos = tabla2[,"Gene"] %in% matches
            print(table(pos))
            tabla2 = tabla2[pos,]
            nom_fichero = paste(bus, "/", familia, "/", familia, "_genot_causante_", i, "_panel", sep = "")
            write.csv2(tabla2, file.path (.job$dir$proces, nom_fichero), row.names = FALSE, quote = FALSE)
            
            # y la añadimas a la tabla global
            tabla_global_panel <- rbind(tabla_global_panel,tabla2)
          }
          if (bus == "dominante") {
            tabla_global_panel_dominante = rbind(tabla_global_panel_dominante, tabla_global_panel) 
            tabla_global_dominante = rbind(tabla_global_dominante, tabla_global)  
            tabla_global <- as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla)))
            tabla_global_panel <- as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla)))
          }
          if (bus == "recesivo") {
            tabla_global_panel_recesivo = rbind(tabla_global_panel_recesivo, tabla_global_panel) 
            tabla_global_recesivo = rbind(tabla_global_recesivo, tabla_global)
            tabla_global <- as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla)))
            tabla_global_panel <- as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla)))
          }
          if (bus == "novo") {
            tabla_global_panel_novo = rbind(tabla_global_panel_novo, tabla_global_panel) 
            tabla_global_novo = rbind(tabla_global_novo, tabla_global)  
            tabla_global <- as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla)))
            tabla_global_panel <- as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla)))
          }
          if (bus == "novo2") {
            tabla_global_panel_novo2 = rbind(tabla_global_panel_novo2, tabla_global_panel) 
            tabla_global_novo2 = rbind(tabla_global_novo2, tabla_global)  
            tabla_global <- as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla)))
            tabla_global_panel <- as.data.frame(setNames(replicate(length(colnames(tabla)),character(0), simplify = F), colnames(tabla)))
          }
        }
      }
    }
  
    for (bus in busqueda) {
      if (bus == "dominante") {
        tabla_global_panel = tabla_global_panel_dominante
        tabla_global = tabla_global_dominante
      }
      if (bus == "recesivo") {
        tabla_global_panel = tabla_global_panel_recesivo
        tabla_global = tabla_global_recesivo
      }
      if (bus == "novo") {
        tabla_global_panel = tabla_global_panel_novo
        tabla_global = tabla_global_novo
      }
      if (bus == "novo2") {
        tabla_global_panel = tabla_global_panel_novo2
        tabla_global = tabla_global_novo2
      }  
      dim(tabla_global)
      variantes_global = nrow(tabla_global)
      dim(tabla_global_panel)
      variantes_global_panel = nrow(tabla_global_panel)
      
      # Ordenamos por cromosoma y posicion
      chrOrder <-c((1:22),"X","Y")
      tabla_global[,"Chromosome"] <- factor(tabla_global[,"Chromosome"], chrOrder, ordered=TRUE)
      tabla_global = tabla_global[do.call(order, tabla_global[, c("Chromosome", "Position")]), ]
      
      tabla_global_panel[,"Chromosome"] <- factor(tabla_global_panel[,"Chromosome"], chrOrder, ordered=TRUE)
      tabla_global_panel = tabla_global_panel[do.call(order, tabla_global_panel[, c("Chromosome", "Position")]), ]
      
      
      # Guardamos los ficheros a entregar a los investigadores, con todas las variantes agrupadas
      nom_fichero = paste(bus, "/", familia, "/", familia, "_global", sep = "")
      write.csv2(tabla_global, file.path (.job$dir$proces, nom_fichero), row.names = FALSE, quote = FALSE)
      nom_fichero = paste(bus, "/", familia, "/", familia, "_panel_global", sep = "")
      write.csv2(tabla_global_panel, file.path (.job$dir$proces, nom_fichero), row.names = FALSE, quote = FALSE)
      
      # Guardamos el fichero con el numero de variantes
      # OJO, a la hora de automatizar habŕa que mirar si esta creado ya o no
      nom_fichero = paste(bus, "/", "nvariantes", sep = "")
      suppressWarnings(dir.create(file.path(.job$dir$proces, nom_fichero)))
      nom_fichero = paste(bus, "/", "nvariantes/", familia, "_nvariantes", sep = "")
      tabla_genes <- data.frame(familia = numeric(0), variantes_csv = numeric(0), variantes_global = numeric(0), 
                                variantes_panel_global = numeric(0))
      tabla_genes[1, "familia"] = familia
      tabla_genes[1, "variantes_csv"] = variantes_csv
      tabla_genes[1, "variantes_global"] = variantes_global
      tabla_genes[1, "variantes_panel_global"] = variantes_global_panel
      
      write.csv2(tabla_genes, file.path (.job$dir$proces, nom_fichero), row.names = FALSE, quote = FALSE)
    }
  }
}
###EXIT
#warnings ()
#sessionInfo ()
#q ("no")





## Busqueda dominante
# Se buscan variantes que esten presentes en todos los casos y ausentes en los controles. Solo se contempla la homocigosis 
# cuando hay información de los dos progenitores y ambos estan enfermos.

## Busqueda recesiva
# Se buscan variantes que esten en homocigosis en todos los casos y en heterocigosis o ausentes en los controles.

## Busqueda novo
# Se buscan variantes que se presenten en homocigosis para descendientes, y solo esten presentes en un alelo de los 
# progenitores e.j. (progenitor1 = 1/0, progenitor2 = 0/0, descendiente = 1/1)
# Se buscan variantes presentes en descendientes y ausentes en los progenitores
# SOLO SE ANALIZAN DESCENDIENTES CON LOS DOS PROGENITORES SANOS
# No tiene sentido utilizarlo si no se tienen datos de dos progenitores y un descendeinte directo

## Busqueda novo2
# Se buscan variantes que se presenten en homocigosis para descendientes, y solo esten presentes en un alelo de los 
# progenitores e.j. (progenitor1 = 1/0, progenitor2 = 0/0, descendiente = 1/1)
# Se buscan variantes presentes en descendientes y ausentes en los progenitores
# SE ANALIZAN TODOS LOS INDIVIDUOS ENFERMOS DE LOS QUE TENGAMOS DATOS DE LOS PROGENITORES
# No tiene sentido utilizarlo si no se tienen datos de dos progenitores y un descendeinte directo