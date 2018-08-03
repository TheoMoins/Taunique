Training<-function(file_data, file_ann, kmax=4){
  #
  # Génère un modèle de prédiction à partir des données "file_data" annotés dans "file_ann".
  # Le modèle est ensuite sauvegardé dans un fichier nommé "training.rda".
  #
  # file_data : fichier contenant les valeurs d'intensité en fonction du temps
  # file_ann : fichier généré par la fonction training, contenant le modèle de prédiction
  # kmax : nombre maximale de segments possibles
  #
  
  # Importation des intensités en enlevant la première colonne qui contient un label :
  datamatrix=read.table(file_data, header=TRUE, sep = ",")
  datamatrix = datamatrix[,-1]   
  N=dim(datamatrix)[2]

  ann = read.table(file_ann, header = TRUE, sep = ",")   # Importation des annotations
  
  data.seq = data.frame()
  
  segs = NULL
  selection = NULL

  # Pour chaque signal annoté :
  for(i in 1:N){   
    Intensity = datamatrix[,i]
    Intensity = Intensity[!is.na(Intensity)]
    
    loss.vec <- rep(NA, kmax)
    
    # Normalisation de l'intensité (nécessaire pour utiliser le modèle Gaussien de Segmentor3IsBack)
    Intensity = Intensity/sd(Intensity)
    # Calcul des segments pour toutes les valeurs possibles de k
    segIntensity=Segmentor(Intensity, model = 2, Kmax = kmax)
    
    # Réorganisation des données (obtention du vecteur "changes") :
    for(n.segments in 1:kmax){
      end <- segIntensity@breaks[n.segments, 1:n.segments]
      data.before.change <- end[-n.segments]
      data.after.change <- data.before.change+1
      pos.before.change <- as.integer((data.before.change+data.after.change)/2)
      start <- c(1, data.after.change)
      segStart <- c(1, pos.before.change)
      segEnd <- c(pos.before.change, length(Intensity))
      seg.mean.vec <- segIntensity@parameters[n.segments, 1:n.segments]
      data.mean.vec <- rep(seg.mean.vec, end-start+1)
      segs = rbind(segs, data.table(protein=i, start, end, n.segments, changepoint = segStart, segEnd, mean = seg.mean.vec))
      
      # Calcul de la perte :
      loss.vec[n.segments] <- sum((Intensity-data.mean.vec)^2)
      
    }
    loss.dt = data.table(protein = i, n.segments = 1:kmax, loss = loss.vec)
    
    # Calcul de la pénalité correspondant à chaque choix de nombre de segments : 
    selection = rbind(selection, modelSelection(loss.dt, complexity = "n.segments"))
    # data.seq = rbind(data.seq, data.frame(protein = i, intensity = Intensity))
  }
  changes = segs[1 < start]
  
  # Calcul de l'erreur en fonction de l'étiquetage des signaux : 
  error.list <- labelError(
    selection, ann, changes,
    problem.vars="protein",
    model.vars="n.segments",
    change.var="changepoint",
    label.vars=c("min", "max"))
  
  # Selection du bon interval pour lambda : 
  error.join = error.list$model.errors[J(selection), on = list(protein, n.segments, loss)]
  target = targetIntervals(error.join, problem.vars = "protein")
  
  
  for(i in 1:ncol(datamatrix)){
    intensity = datamatrix[,i]
    intensity = intensity[!is.na(intensity)]
    # intensity = intensity/sd(intensity)

    data.seq = rbind(data.seq, data.frame(protein = i, intensity =intensity))
  }

  # Calcul d'un ensemble de caractéristiques de f, susceptibles d'être pris en compte pour
  # le calcul de la pénalité : 
  f.mat <- featureMatrix(data.seq, "protein", "intensity")
  
  t.mat = target[,c("min.log.lambda", "max.log.lambda")]
  
  train.f.mat = subset(f.mat, is.finite(as.matrix(t.mat)[,1]) & is.finite(as.matrix(t.mat)[,2]))

  t.mat = subset(t.mat, is.finite(as.matrix(t.mat)[,1]) & is.finite(as.matrix(t.mat)[,2]))
  t.mat = as.matrix(t.mat)

  # Calcul du modèle : 
  fit = IntervalRegressionCV(train.f.mat, t.mat, verbose=0)
  # Sauvegarde du modèle : 
  save(fit, file = paste(dirname(file_ann),"/training.rda",sep=""))
}


Choix_k<-function(file_data, file_fit, kmax = 4){
  #
  # Retourne la liste contenant le choix du nombre de segments pour chaque signal.
  #
  # file_data : fichier contenant les valeurs d'intensité en fonction du temps
  # file_fit : fichier généré par la fonction training, contenant le modèle de prédiction
  # kmax : nombre maximale de segments possibles
  #
  
  # Importation des intensités en enlevant la première colonne qui contient un label :
  datamatrix=read.table(file_data, header=TRUE, sep = ",")
  datamatrix = datamatrix[,-1]
  
  data.seq = data.frame()
  selection = NULL
  segs = NULL
  
  # Pour chaque signal, on calcule les intervalles de lambda correspondant à chaque choix de nombre
  # de marches, puis on sélectionne le meilleur intervalle pour chaque signal en fonction du nombre 
  # de segments attendus et enfin on construit une regression qui minimise le nombre d'erreurs (qui
  # passe entre le plus d'intervalles) pour en déduire une prédiction du nombre de segments des signaux
  for(i in 1:ncol(datamatrix)){
    Intensity = datamatrix[,i]
    Intensity = Intensity[!is.na(Intensity)]

    loss.vec <- rep(NA, kmax)

    # Normalisation de l'intensité (nécessaire pour utiliser le modèle Gaussien de Segmentor3IsBack)
    et = sd(Intensity)
    Intensity = Intensity/et

    # Calcul des segments pour toutes les valeurs possibles de k
    segIntensity=Segmentor(Intensity, model = 2, Kmax=kmax)
    Intensity = Intensity*et

    data.seq = rbind(data.seq, data.frame(protein = i, intensity = Intensity))

    # Réorganisation des données (obtention du vecteur "changes") :
    for(n.segments in 1:kmax){
      end <- segIntensity@breaks[n.segments, 1:n.segments]
      data.before.change <- end[-n.segments]
      data.after.change <- data.before.change+1
      pos.before.change <- as.integer((data.before.change+data.after.change)/2)
      start <- c(1, data.after.change)
      segStart <- c(1, pos.before.change)
      segEnd <- c(pos.before.change, length(Intensity))
      seg.mean.vec <- segIntensity@parameters[n.segments, 1:n.segments]
      data.mean.vec <- rep(seg.mean.vec, end-start+1)
      segs = rbind(segs, data.table(protein=i, start, end, n.segments, changepoint = segStart, segEnd, mean = seg.mean.vec))

      # Calcul de la perte :
      loss.vec[n.segments] <- sum((Intensity/et-data.mean.vec)^2)

    }
    loss.dt = data.table(protein = i, n.segments = 1:kmax, loss = loss.vec)

    # Calcul du choix de la pénalité correspondant à chaque choix de nombre de segments : 
    selection = rbind(selection, modelSelection(loss.dt, complexity = "n.segments"))
  }

  val.f.mat <- featureMatrix(data.seq, "protein", "intensity")
  
  load(file = file_fit)

  # Prediction de lambda pour chaque signal :
  pred = data.frame(predict(fit, val.f.mat), paste(c(1:ncol(datamatrix))))
  names(pred) = c("pred.log.lambda", "protein")

  # Calcul du nombre de segments correspondant à chaque pénalité predite :
  list_k = NULL
  for(i in 1:ncol(datamatrix)){
    selection_i  = subset(selection, protein == i)
    pred_i  = pred[i,1]
    for (j in 1:nrow(selection_i)){
      if(selection_i[j, 3] < pred_i & selection_i[j,4] > pred_i){
        list_k = rbind(list_k, selection_i[j, 6])
      }
    }
  }
  return(list_k)
}  


Calcul_SegmentsALL<-function(file, list_k, kmax = 4, affichage = 0){
  #
  # Calcul des segmentations à partir de la liste contenant le nombre de 
  # segments pour chaque signal.
  #
  # file : fichier contenant les valeurs des signaux en fonction du temps
  # list_k : liste de la taille le nombre de signaux, et qui indique le nombre de segments de chaque signaux
  # kmax : nombre maximale de segments possibles
  # Affichage : 0 pour ne pas afficher l'ensemble des graphes analysés, 1 sinon
  #

  # Importation des intensités en enlevant la première colonne qui contient un label :
  datamatrix=read.table(file, header=TRUE, sep = ",")
  datamatrix = datamatrix[,-1]
  N=dim(datamatrix)[2]

  Kfinal = matrix(0, nrow = N, ncol = 1)
  marchesfin = matrix(0, nrow = N, ncol = kmax+1)
  breaksfinal = matrix(0, nrow = N, ncol = kmax+1)
  ecart_type = matrix(0, nrow = N, ncol = 1)
  len = matrix(0, nrow = N, ncol = 1)
  likelihood = matrix(0, nrow = N, ncol = 1)
  
  for(i in 1:N){
    k_i = list_k[i]
    Intensity = datamatrix[,i]
    Intensity = Intensity[!is.na(Intensity)]
    variance = sd(Intensity)
    ecart_type[i,] = variance
    len[i,] = length(Intensity)
    tmp = Calcul_Segments(Intensity, kmax, k_i, affichage,i)
    Kfinal[i,] = tmp$Kfinal
    marchesfin[i,1:Kfinal[i,]] = tmp$marchesfin
    breaksfinal[i,1:Kfinal[i,]] = tmp$breaksfinal
    likelihood[i] = tmp$L
  }
  list(nbmarches = Kfinal, moymarches = marchesfin, absmarches= breaksfinal, sigma = ecart_type, size = len, L = likelihood)
}


Calcul_Segments<-function(Intensity, kmax, Ks, affichage, i = 0){
  #
  # Calcul de la segmentation d'un signal en Ks marches.
  #
  # Intensity : vecteur contenant l'intensité d'un signal au cours du temps
  # kmax : nombre maximale de segments possibles
  # Ks : nombre de segments prédit du signal
  # affichage : 0 pour ne pas afficher l'ensemble des graphes analysés, 1 sinon
  # i : index du signal parmi l'ensemble de départ
  #

  variance = sd(Intensity)
  Intensity = Intensity/variance
  
  segIntensity=Segmentor(Intensity, model = 2 , Kmax=kmax)
  
  L = getLikelihood(segIntensity)

  param_int = matrix(ncol = kmax+1, nrow = kmax)
  param_int[,1:kmax] = getParameters(segIntensity)
  param_int[,kmax+1] = rep(0,kmax)
  
  marches=param_int[Ks,1:Ks]
  marchesuivantes=param_int[Ks,2:(Ks+1)]

  
  break_int = matrix(ncol = kmax+1, nrow = kmax)
  break_int[,1:kmax] = getBreaks(segIntensity)
  break_int[,kmax+1] = rep(0,kmax)
  
  breaks=break_int[Ks,1:Ks]
  breaksuivants=break_int[Ks,2:(Ks+1)]
  
  marches = marches * variance
  
  # Affichage des courbes :
  if(affichage != 0){
    ts.plot(Intensity*variance, xlab = "Time", ylab = i)
    for (j in 1:Ks){
      segments(x0 = c(0,breaks)[j], y0 = marches[j], x1 = c(0,breaks)[j+1], y1 = marches[j], col="red")
      if(j<Ks){
        segments(x0 = c(0,breaks)[j+1], y0 = marches[j], x1 = c(0,breaks)[j+1], y1 = marches[j+1], col="red")
      }
    }
  }
  list(marchesfin=marches,Kfinal=Ks, breaksfinal=breaks, L = getLikelihood(segIntensity)[Ks])
} 

