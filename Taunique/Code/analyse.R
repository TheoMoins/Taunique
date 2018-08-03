Marches_coord <- function(resultat){
  #
  # Calculs de différentes caractéristiques de segmentations obtenues par FusionALL ou Calcul_SegmentsALL :
  # nombres de graphes finalement analysés (à 2 ou 3 marches), tailles des marches, abscisse des points de rupture
  #
  # resultat : liste obtenue avec les fonctions FusionALL ou Calcul_SegmentsALL
  #
  
  nbgraphes_i = length(resultat$nbmarches)
  nbgraphes = 0
  
  absgraphes = NULL
  
  new_nbmarches = rep(0,nbgraphes_i)
  Premiere_marche = rep(0,nbgraphes_i)
  Deuxieme_marche = rep(0,nbgraphes_i)
  Premiere_abs = rep(0,nbgraphes_i)
  Deuxieme_abs = rep(0,nbgraphes_i)
  
  for(i in 1:nbgraphes_i){
    if (resultat$nbmarches[i,1]<4 & resultat$nbmarches[i,1] > 1){
      if(all(sort(resultat$moymarches[i,], decreasing = TRUE) == resultat$moymarches[i,])){
        nbgraphes = nbgraphes + 1
        absgraphes = c(absgraphes,i)
        new_nbmarches[i] = resultat$nbmarches[i,1]-1
        Premiere_marche[i] = resultat$moymarches[i,1]-resultat$moymarches[i,2]
        Premiere_abs[i] = resultat$absmarches[i,1]
        if(resultat$nbmarches[i,1] == 3){
          Deuxieme_marche[i] = resultat$moymarches[i,2]-resultat$moymarches[i,3]
          Deuxieme_abs[i] = resultat$absmarches[i,2]
        }
      }
    }
  }
  list(nbgraphes = nbgraphes, Premiere_marche = Premiere_marche, Deuxieme_marche = Deuxieme_marche, Premiere_abs = Premiere_abs, Deuxieme_abs = Deuxieme_abs, absgraphes = absgraphes)
}


f1 <- function(x){
  if(x>0){
    return(0.5*(tanh(x-1.5)+1))
  } else {
    return(0)
  }
}

f2 <- function(x){
  if(x<0.5){
    return(1)
  } else {
    return(2-2*x)
  }
}

Indice_confiance <- function(marches,resultat){
  #
  # Calcul de l'indice de confiance d'une segmentation.
  #
  # marches : liste obtenue avec la fonction Marches_Coord
  # resultat : liste obtenue avec les fonctions FusionALL ou Calcul_SegmentsALL
  #
  ratio1 = marches$Premiere_marche/resultat$sigma
  ratio2 = marches$Deuxieme_marche/resultat$sigma

  p1 = marches$Premiere_marche
  p1[marches$Deuxieme_marche == 0] = sapply(ratio1[marches$Deuxieme_marche == 0], f1)
  if(length(ratio2[marches$Deuxieme_marche != 0])>0){
    p1[marches$Deuxieme_marche != 0] = sqrt(sapply(2*ratio1[marches$Deuxieme_marche != 0], f1)*sapply(2*ratio2[marches$Deuxieme_marche != 0], f1))
  }
  
  x = marches$Premiere_abs
  x[marches$Deuxieme_abs != 0] = marches$Deuxieme_abs[marches$Deuxieme_abs != 0] 
  x = x/resultat$size
  p2 = sapply(x,f2)
  
  borne_inf = 0.4
  p3 = 1-borne_inf*(resultat$L-min(resultat$L))/(max(resultat$L)-min(resultat$L))
  
  alpha1 = 0.5
  alpha2 = 0.2
  alpha3 = 0.3
  
  return((p1**alpha1)*(p2**alpha2)*(p3**alpha3))
}
  

affiche_marche<-function(resultat, data, marches, i, IC){
  #
  # Fonction permettant d'afficher la segmentation numéro i
  #
  # resultat : liste obtenue avec les fonctions FusionALL ou Calcul_SegmentsALL
  # data : matrice contenant l'ensemble des intensités en fonction du temps
  # marches : liste obtenue avec la fonction Marches_Coord
  # i : index du signal
  # IC : Liste contenant les indices de confiance du signal (obtenue avec la fonction Indice_confiance)
  #
  
  plot(data[,i], xlab = "Time", ylab = "Intensity", type = "l")
  mtext(paste("IC : ", round(IC[i],3), sep = ""), side = 3, line = -2, col = "blue")
  title(main = paste("Graphe n°",i,sep=""))
  
  if(resultat$nbmarches[i] == 1){
    
    segments(x0 = 0, y0 = resultat$moymarches[i,1], 
             x1 = resultat$size[i], y1 = resultat$moymarches[i,1], 
             col="red")
    
  } else if(resultat$nbmarches[i] == 2){
    
    segments(x0 = 0, y0 = resultat$moymarches[i,1], 
             x1 = marches$Premiere_abs[i], y1 = resultat$moymarches[i,1], 
             col="red")    
    segments(x0 = marches$Premiere_abs[i], y0 = resultat$moymarches[i,1], 
             x1 = marches$Premiere_abs[i], y1 = resultat$moymarches[i,2], 
             col="red")
    segments(x0 = marches$Premiere_abs[i], y0 = resultat$moymarches[i,2], 
             x1 = resultat$size[i], y1 = resultat$moymarches[i,2], 
             col="red")
    
  } else if(resultat$nbmarches[i] == 3){
    
    segments(x0 = 0, y0 = resultat$moymarches[i,1], 
             x1 = marches$Premiere_abs[i], y1 = resultat$moymarches[i,1], 
             col="red")
    segments(x0 = marches$Premiere_abs[i], y0 = resultat$moymarches[i,1], 
             x1 = marches$Premiere_abs[i], y1 = resultat$moymarches[i,2], 
             col="red")
    segments(x0 = marches$Premiere_abs[i], y0 = resultat$moymarches[i,2], 
             x1 = marches$Deuxieme_abs[i], y1 = resultat$moymarches[i,2], 
             col="red")
    segments(x0 = marches$Deuxieme_abs[i], y0 = resultat$moymarches[i,2], 
             x1 = marches$Deuxieme_abs[i], y1 = resultat$moymarches[i,3],
             col="red")
    segments(x0 = marches$Deuxieme_abs[i], y0 = resultat$moymarches[i,3],
             x1 = resultat$size[i], y1 = resultat$moymarches[i,3], 
             col="red")
  }
}

