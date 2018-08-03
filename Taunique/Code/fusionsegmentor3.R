FusionALL<-function(file, kmax=10, ratio=8, affichage = 0){
  #
  # Calcul des segmentations avec la pénalité de Lebarbier
  #
  # file : fichier contenant les valeurs des signaux en fonction du temps
  # kmax : nombre maximale de segments possibles
  # ratio : rapport minimum entre l'amplitude de signal (max - min) et la taille des marches calculées
  # Affichage : 0 pour ne pas afficher l'ensemble des graphes analysés, 1 sinon
  #
  
  # model = 2 correspond au modèle Gaussien dans Segmentor3IsBack
  model = 2 
  
  # Importation des intensités en enlevant la première colonne qui contient un label :
  datamatrix=read.table(file, header=TRUE, sep = ",")
  datamatrix = datamatrix[,-1]
  N=dim(datamatrix)[2]
  
  Kfinal = matrix(0, nrow = N, ncol = 1)
  marchesfin = matrix(0, nrow = N, ncol = 5)
  breaksfinal = matrix(0, nrow = N, ncol = 5)
  pen = NULL
  ecart_type = matrix(0, nrow = N, ncol = 1)
  len = matrix(0, nrow = N, ncol = 1)
  likelihood = matrix(0, nrow = N, ncol = 1)
  
  # Pour chaque signal annoté :
  for(i in 1:N){
    
    Intensity = datamatrix[,i]
    Intensity = Intensity[!is.na(Intensity)]

    if(!is.numeric(Intensity)){
      print("ATTENTION : Données non numériques")
    } else {
      et = sd(Intensity)
      ecart_type[i,] = et
      len[i,] = length(Intensity)
      tmp = Fusion(Intensity, kmax, ratio, affichage,i)
      Kfinal[i,] = tmp$Kfinal
      marchesfin[i,1:Kfinal[i,]] = tmp$marchesfin
      breaksfinal[i,1:Kfinal[i,]] = tmp$breaksfinal
      likelihood[i] = tmp$L
      pen = rbind(pen, c(i, log(2*log(length(Intensity)/tmp$Kfinal)+5)))
    }
  }
  list(nbmarches = Kfinal, moymarches = marchesfin, absmarches= breaksfinal, sigma = ecart_type, size = len, L = likelihood)
}


Fusion<-function(Intensity, kmax, ratio, affichage, ii = 0){
  #
  # Calcul de la segmentation d'un signal avec la pénalité de Lebarbier.
  #
  # Intensity : vecteur contenant l'intensité d'un signal au cours du temps
  # kmax : nombre maximale de segments possibles
  # ratio : rapport minimum entre l'amplitude de signal (max - min) et la taille des marches calculées
  # affichage : 0 pour ne pas afficher l'ensemble des graphes analysés, 1 sinon
  # ii : index du signal
  #

  # model = 2 correspond au modèle Gaussien dans Segmentor3IsBack
  model = 2
  
  # Normalisation de l'intensité (nécessaire pour utiliser le modèle Gaussien de Segmentor3IsBack)
  et = sd(Intensity)
  Intensity = Intensity/et
  
  # segmentations de K=1 à K=kmax
  segIntensity=Segmentor(Intensity,model,Kmax=kmax)
  
  # Selection de la segmentation optimale
  Ks=SelectModel(segIntensity)
  
  # Fusion des fausses marches et des breaks 
  # vecteur des marches (moyennes des lois segments) Ks marches ou segments
  marches=getParameters(segIntensity)[Ks,1:Ks]
  breaks=getBreaks(segIntensity)[Ks,1:Ks]
  marchesuivantes=getParameters(segIntensity)[Ks,2:(Ks+1)]
  breaksuivants=getBreaks(segIntensity)[Ks,2:(Ks+1)]
  
  hauteurs=abs((marches-marchesuivantes)[1:(Ks-1)])

  # calcul du seuil pour éliminer les "petites contre"fausses" marches: max/ratio
  seuilhauteur = (max(Intensity)-min(Intensity))/ratio
  
  # indices des contre-marches sous le seuil
  indices=seq(1:(Ks-1))[hauteurs<seuilhauteur]
  
  # nb de contre marches inférieures au seuil
  nbindices=sum(hauteurs<seuilhauteur)
  
  while(nbindices>0){
    newbreaks = NULL
    newmarches = NULL
    for(i in 1:Ks){
      if(i < Ks && hauteurs[i]<seuilhauteur){
        s1 = breaks[i]
        s2 = breaks[i+1]-breaks[i]
        if(i>1){
          s1 = s1 - breaks[i-1]
        }
        newmarches = c(newmarches, (marches[i]*s1+marches[i+1]*s2)/(s1+s2))
        newbreaks = c(newbreaks,breaks[i+1])
      } else if(i == 1 || hauteurs[i-1]>seuilhauteur){
        newmarches = c(newmarches, marches[i])
        newbreaks = c(newbreaks, breaks[i])
      }
    }
    newlong=length(newmarches)
    newlongbreaks=length(newbreaks)
    marches=newmarches
    breaks=newbreaks
    marchesuivantes=c(newmarches[2:newlong],0)
    Ks=newlong
    hauteurs=(marches-marchesuivantes)[1:(Ks-1)]
    # hauteur max: contre marche la plus grande
    maxhauteur=max(hauteurs)
    # calcul du seuil pour éliminer les petites contre marches: max/ratio
    seuilhauteur=maxhauteur/ratio
    nbindices=sum(hauteurs<seuilhauteur)
    if(is.na(nbindices)){
      nbindices = 0
    }
  }
  
  # Affichage des courbes :
  marches = marches * et
  if(affichage != 0){
    ts.plot(Intensity*et, xlab = "Time", ylab =ii)
    for (i in 1:Ks){
      segments(x0 = c(0,breaks)[i], y0 = marches[i], x1 = c(0,breaks)[i+1], y1 = marches[i], col="red")
      if(i<Ks){
        segments(x0 = c(0,breaks)[i+1], y0 = marches[i], x1 = c(0,breaks)[i+1], y1 = marches[i+1], col="red")
      }
    }
  }
  list(marchesfin=marches, Kfinal=Ks, breaksfinal=breaks, L = getLikelihood(segIntensity)[Ks])
} 

