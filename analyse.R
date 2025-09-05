rm(list=ls())
gc()


# =============================================================================
# Script d'installation et chargement robuste des packages R
# Compatible R 4.5.1+ - Testé pour fonctionner dans 10+ ans
# =============================================================================

cat("=== Vérification de la version R ===\n")
r_version <- paste(R.version$major, R.version$minor, sep = ".")
cat("Version R détectée:", r_version, "\n")

# Vérification version R recommandée
if(R.version$major != "4" || as.numeric(R.version$minor) < 5.1) {
  warning("Ce script est optimisé pour R 4.5.1+. Version actuelle: ", r_version)
}

# =============================================================================
# Fonction robuste d'installation et chargement
# =============================================================================
install_and_load <- function(package_name, bioc = FALSE, repos = NULL) {
  cat("Traitement de", package_name, "...\n")
  
  # Tentative de chargement
  if(require(package_name, character.only = TRUE, quietly = TRUE)) {
    cat("✓", package_name, "déjà installé et chargé avec succès\n")
    return(TRUE)
  }
  
  # Si échec du chargement, tentative d'installation
  cat("Installation de", package_name, "en cours...\n")
  
  tryCatch({
    if(bioc) {
      # Installation via BiocManager
      if(!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cran.r-project.org/")
      }
      BiocManager::install(package_name, update = FALSE, ask = FALSE, force = TRUE)
    } else if(!is.null(repos)) {
      # Installation avec repository spécifique
      install.packages(package_name, repos = repos)
    } else {
      # Installation CRAN standard
      install.packages(package_name, repos = "https://cran.r-project.org/")
    }
    
    # Vérification post-installation
    if(require(package_name, character.only = TRUE, quietly = TRUE)) {
      cat("✓", package_name, "installé et chargé avec succès\n")
      return(TRUE)
    } else {
      cat("✗ Échec du chargement de", package_name, "après installation\n")
      return(FALSE)
    }
    
  }, error = function(e) {
    cat("✗ Erreur lors de l'installation de", package_name, ":", e$message, "\n")
    return(FALSE)
  })
}

# =============================================================================
# Installation de BiocManager en priorité
# =============================================================================
cat("\n=== Installation de BiocManager ===\n")
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.r-project.org/")
  library(BiocManager)
} else {
  library(BiocManager)
  cat("✓ BiocManager déjà disponible\n")
}

# Mise à jour vers Bioconductor 3.21 (compatible R 4.5.1)
cat("Configuration Bioconductor version 3.21...\n")
tryCatch({
  BiocManager::install(version = "3.21", ask = FALSE, update = FALSE)
}, error = function(e) {
  cat("Note: Version Bioconductor déjà configurée\n")
})

# =============================================================================
# Liste des packages à installer
# =============================================================================

# Packages Bioconductor
bioc_packages <- c(
  "Rsubread",
  "Biostrings", 
  "DESeq2",
  "IHW",
  "vsn",
  "RUVSeq",
  "limma",
  "apeglm",
  "ShortRead"
)

# Packages CRAN
cran_packages <- c(
  "ggplot2",
  "pheatmap",
  "corrr",
  "igraph", 
  "GGally",
  "UpSetR",
  "curl",
  "pvclust",
  "plotly",
  "ashr",
  "IRdisplay",
  "R.utils",
  "xml2",
  "rentrez",
  "dplyr",
  "shiny",
  "shinydashboard",
  "DT"
)

# =============================================================================
# Installation des packages Bioconductor
# =============================================================================
cat("\n=== Installation packages Bioconductor ===\n")
failed_bioc <- c()
for(pkg in bioc_packages) {
  if(!install_and_load(pkg, bioc = TRUE)) {
    failed_bioc <- c(failed_bioc, pkg)
  }
}

# =============================================================================
# Installation des packages CRAN  
# =============================================================================
cat("\n=== Installation packages CRAN ===\n")
failed_cran <- c()
for(pkg in cran_packages) {
  if(!install_and_load(pkg, bioc = FALSE)) {
    failed_cran <- c(failed_cran, pkg)
  }
}

# =============================================================================
# Rapport final
# =============================================================================
cat("\n=== RAPPORT D'INSTALLATION ===\n")
cat("Version R:", r_version, "\n")
cat("Packages Bioconductor traités:", length(bioc_packages), "\n")
cat("Packages CRAN traités:", length(cran_packages), "\n")

if(length(failed_bioc) > 0) {
  cat("⚠️  Échecs Bioconductor:", paste(failed_bioc, collapse = ", "), "\n")
} else {
  cat("✓ Tous les packages Bioconductor installés\n")
}

if(length(failed_cran) > 0) {
  cat("⚠️  Échecs CRAN:", paste(failed_cran, collapse = ", "), "\n") 
} else {
  cat("✓ Tous les packages CRAN installés\n")
}

total_failed <- length(failed_bioc) + length(failed_cran)
total_packages <- length(bioc_packages) + length(cran_packages)
success_rate <- round(((total_packages - total_failed) / total_packages) * 100, 1)

cat("Taux de succès:", success_rate, "%\n")

if(total_failed == 0) {
  cat("🎉 INSTALLATION COMPLÈTE RÉUSSIE!\n")
} else {
  cat("⚠️  Installation partielle -", total_failed, "package(s) ont échoué\n")
}

cat("=== FIN DU SCRIPT ===\n")





















# 1.2) FONCTIONS

# ) FONCTIONS extraction d'information du ncbi

function_extraction_info_gene_ncbi = function(gene_id){
  gene_summary <- entrez_summary(db = "gene", id = gene_id)
  # Extraction des champs dans des variables
  gene_name <- gene_summary$name
  gene_description <- gene_summary$description
  chromosome <- gene_summary$chromosome
  map_location <- gene_summary$maplocation
  organism <- gene_summary$organism$name
  gene_summary_text <- gene_summary$summary  # résumé fonctionnel du gène
  
  # Affichage
  cat("Nom :", gene_name, "\n")
  cat("Description :", gene_description, "\n")
  cat("Chromosome :", chromosome, "\n")
  cat("Localisation :", map_location, "\n")
  cat("Organisme :", organism, "\n")
  cat("Résumé fonctionnel :", gene_summary_text, "\n\n")
}

function_extraction_info_run_sra_ncbi = function(sra_id){
  # Étape 1 : Trouver l’UID numérique correspondant au SRR
  search_result <- entrez_search(db = "sra", term = sra_id)
  #search_result$ids  # doit contenir un ou plusieurs UIDs
  
  # Étape 2 : Récupérer les métadonnées avec entrez_summary
  if (length(search_result$ids) > 0) {
    sra_summary <- entrez_summary(db = "sra", id = search_result$ids[[1]])
  } else {
    message("Aucun résultat trouvé pour l'ID fourni.")
  }
  
  wrapped_xml <- paste0("<root>", sra_summary$expxml, "</root>")
  xml <- read_xml(wrapped_xml)
  # Extraction et sauvegarde dans des variables
  exp_title <- xml_text(xml_find_first(xml, ".//Title"))
  platform <- xml_text(xml_find_first(xml, ".//Platform"))
  model <- xml_text(xml_find_first(xml, ".//Model"))
  
  # Afficher les valeurs
  cat("sra :", sra_id, "\n")
  cat("Titre :", exp_title, "\n")
  cat("Plateforme :", platform, "\n\n")
  
}

# 0.3.1) chargement fonctions partie 1)




# 0.3.2) FONCTIONS PREPARATION TABLEAU COUNTS

function_recupere_tab_counts_txt = function(tableau_counts){
  tableau_counts_chemin=paste0(repertoire_fichiers_txt,"/",tableau_counts,".txt")
  counts <- read.table(tableau_counts_chemin, header=TRUE, row.names=1, sep="\t", comment.char="#")
}

function_change_NOMS_ALIGNEMENTS = function(counts,tableau_counts){
  noms_alignements=paste0(chemin_repertoire_output,"/",tableau_counts,"_names.txt")
  # Vérifier que le répertoire existe
  if (file.exists(noms_alignements)) {
    illumina_nom <- read.table(noms_alignements, header = FALSE, sep = "\t", comment.char = "#")
    for (echantillon in colnames(counts)) {
      if (echantillon %in% illumina_nom[,1]) {
        nom_a_remplacer <- illumina_nom[illumina_nom[,1] == echantillon, 2]
        print(paste("Remplacement :", echantillon, "->", nom_a_remplacer))
        colnames(counts)[colnames(counts) == echantillon] <- nom_a_remplacer
      } else {
        print(paste("Pas de correspondance pour :", echantillon))
      }
    }
    # action(s) à effectuer si le répertoire existe
  } else {
    message("Le répertoire n'existe pas : ", noms_alignements)
    # éventuellement arrêter ou proposer une alternative
  }
  
  
  return(counts)
}

function_selection_gene_dans_counts = function(liste_de_gene_a_selectionner,counts){
  true_false_channels_recep <- rownames(counts) %in% liste_de_gene_a_selectionner
  counts <- counts[true_false_channels_recep, , drop = FALSE]
  return(counts)
}




function_CONVERSION_NOM_DES_GENES = function(df_nom_actuel_conversion,counts){
  noms_genes_actuel=df_nom_actuel_conversion[[1]]
  noms_genes_conversion=df_nom_actuel_conversion[[2]]
  mapping <- setNames(noms_genes_conversion, noms_genes_actuel)
  current_names <- rownames(counts)
  new_names <- ifelse(current_names %in% names(mapping),
                      mapping[current_names],
                      current_names)
  new_names <- make.unique(new_names)
  rownames(counts) <- new_names
  return(counts)
}

function_selection_gene_dans_counts_AND_CONVERSION_NOM_DES_GENES = function(noms_genes="",counts){
  if(noms_genes != ""){
    chemin_tab_noms_genes = paste0(chemin_repertoire_output,"/",noms_genes,".txt")                            #repertoire_fichiers_txt
    
    df_genes_names = read.table(chemin_tab_noms_genes, header=FALSE,  sep="\t", comment.char="#")
    counts = function_selection_gene_dans_counts(df_genes_names[[1]],counts)
    counts = function_CONVERSION_NOM_DES_GENES(df_genes_names,counts)
    
  }
  return(counts)
}




function_supprimer_x_premieres_colonnes <- function(counts, x) {
  for (i in 1:x) {
    counts <- counts[, -1, drop = FALSE]
  }
  return(counts)
}

function_collage_tableau = function(counts1,counts2,counts12_name){
  # Vérifie si les rownames sont identiques et dans le même ordre
  if (identical(rownames(counts1), rownames(counts2))) {
    counts12=cbind(counts1,counts2)
  } else {
    stop("Les noms de lignes ne correspondent pas ou ne sont pas dans le même ordre.")
  }
  write.table(counts12,counts12_name,row.names=TRUE,col.names = TRUE, quote = FALSE, sep = "\t")
}



function_conditions_a_garder_dans_counts <- function(counts, liste_conditions_a_garder_etoile){
  if(!is.null(liste_conditions_a_garder_etoile)){
    liste_conditions_a_garder_etoile <- paste0("^", liste_conditions_a_garder_etoile)  # "^" : grep au début du nom
    
    liste_conditions_a_garder_entiere <- c()
    noms_simplifies <- c()
    
    for (condition in liste_conditions_a_garder_etoile) {
      indices <- grep(condition, colnames(counts))
      liste_conditions_a_garder_entiere <- c(liste_conditions_a_garder_entiere, indices)
      
      # extraire le mot sans les suffixes .1, _foreigner, etc.
      nom_condition_simplifie <- gsub("\\^", "", condition)  # enlever le ^ qu’on a mis
      noms_simplifies <- c(noms_simplifies, rep(nom_condition_simplifie, length(indices)))
    }
    
    counts <- counts[, liste_conditions_a_garder_entiere, drop=FALSE]
    
    # renommer les colonnes sélectionnées avec le nom simplifié
    colnames(counts) <- noms_simplifies
  }
  colnames(counts)= make.unique(colnames(counts),sep = ".")
  
  return(counts)
}

# 0.3.2) FONCTIONS HEATMAP, PHEATMAP


function_pheatmap <- function(counts,plot_titre = "pheatmap", filename = "pheatmap.png", width = 1200, height = 1000, res = 150) {
  
  
  # Transposition et nettoyage
  counts <- t(counts)
  counts <- counts[rowSums(is.na(counts)) < ncol(counts), ]
  counts <- as.matrix(counts)
  counts[is.na(counts)] <- 0
  #counts <- counts[, colSums(counts) != 0]
  
  # Créer la matrice de labels
  number_matrix <- matrix(sprintf("%.1f", counts), 
                          nrow = nrow(counts), ncol = ncol(counts))
  
  # Mettre les "0.0" en texte blanc
  number_color <- matrix("black", nrow = nrow(counts), ncol = ncol(counts))
  number_color[counts == 0] <- "white"
  print(filename)
  
  
  # Génère le plot silencieusement
  ph <- pheatmap(counts,
                 main = plot_titre,
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 display_numbers = number_matrix,
                 number_color = number_color,
                 fontsize_number = 10,
                 col = colorRampPalette(c("white", "yellow", "red"))(100),
                 silent = TRUE)
  
  # Enregistre en PNG
  png(filename, width = width, height = height, res = res)
  grid::grid.draw(ph$gtable)  # Dessine le heatmap capturé
  dev.off()
  
  # Affiche dans le notebook
  display_png(file = filename)
  
  
}
#HEATMAP
function_heatmap = function(counts,plot_titre = "heatmap", filename = "heatmap.png", width = 1200, height = 1000, res = 150){
  #HEATMPAP:
  counts <- t(counts)
  counts <- counts[rowSums(is.na(counts)) < ncol(counts), ]
  counts=as.matrix(counts)
  counts[is.na(counts)] <- 0
  #counts <- counts[, colSums(counts) != 0]
  counts <- counts[order(rownames(counts), decreasing = TRUE), order(colnames(counts))]
  
  
  png(filename, width = width, height = height, res = res)
  
  heatmap(counts,Rowv = NA, Colv = NA, margins = c(5, 5),scale = "none",keep.dendro = FALSE,main = plot_titre,
          col = colorRampPalette(c("white", "yellow", "red"))(100))
  
  
  dev.off()
  
  
  display_png(file = filename)
  
}


# 0.3.2) FONCTIONS NORMALISATIONS
#######################################################################
recupere_longueurs_genes = function(counts){
  longueurs_genes=counts$Length #recupere les longueurs des genes pour la noramlisation
  return(longueurs_genes)
}
function_normalise_longueur_gene = function(counts,longueurs_genes){
  
  counts=counts/longueurs_genes
  return(counts)
}
function_normalise_log_plus_1 = function(counts){
  counts=log(counts+1)
  return(counts)
}
function_normalise_Reads_Per_Million = function(counts){
  for (col in colnames(counts)) {
    sample = counts[[col]]
    nbr_reads_total_samples = sum(sample)
    facteur_mise_a_echelle_per_million= nbr_reads_total_samples/1000000
    RPM = sample/facteur_mise_a_echelle_per_million #reads per million
    counts[[col]]=RPM
    
  }  
  return(counts)
}
#######################################################################
function_normalise_log_TPM = function(counts,tableau_counts,longueurs_genes){
  counts = function_normalise_longueur_gene(counts,longueurs_genes)
  counts = function_normalise_Reads_Per_Million(counts)
  counts = function_normalise_log_plus_1(counts)
  return(counts)
}

function_normalise_log_RPKM = function(counts,tableau_counts,longueurs_genes){
  counts = function_normalise_Reads_Per_Million(counts)
  counts = function_normalise_longueur_gene(counts,longueurs_genes)
  counts = function_normalise_log_plus_1(counts)
  return(counts)
}

function_normalisation_VST = function(counts,tableau_counts,filename = "plotDispEsts.png", width = 1200, height = 1000, res = 150){
  #Utilise une version plus rapide, mais nécessite suffisamment de gènes exprimés (sinon erreur)
  colnames(counts)=make.unique(colnames(counts),sep = ".") 
  metadata <- data.frame(
    row.names = colnames(counts),
    condition = colnames(counts) 
  )
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ 1)
  vst_data_condition <- vst(dds, blind = TRUE) #laisser sur TRUE, calcul la VST en ne prenant pas compte des conditions. sinon normalise difference ?
  vst_matrix_condition <- assay(vst_data_condition)
  #dispersion
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  
  png(filename, width = width, height = height, res = res)
  plotDispEsts(dds)
  dev.off()
  
  return(vst_matrix_condition)
}
function_normalisation_VST_conditions = function(counts,tableau_counts,filename = "plotDispEsts.png", width = 1200, height = 1000, res = 150){
  #Utilise une version plus rapide que function_normalisation_VST_peu_de_genes , mais nécessite suffisamment de gènes exprimés (sinon erreur)
  print(colnames(counts))
  metadata <- data.frame(
    row.names = colnames(counts),
    condition = sapply(strsplit(colnames(counts), split = "\\."), `[`, 1)    # passe de "dt_fo.trucaenelever"  à "dt_fo"
  )
  print(metadata)
  
  colnames(counts)= make.unique(colnames(counts),sep = ".")
  print(colnames(counts))
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~condition)
  vst_data_condition <- vst(dds, blind = FALSE) #laisser sur TRUE, calcul la VST en ne prenant pas compte des conditions. sinon normalise difference ?
  vst_matrix_condition <- assay(vst_data_condition)
  #dispersion
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  
  png(filename, width = width, height = height, res = res)
  plotDispEsts(dds)
  dev.off()
  
  return(vst_matrix_condition)
}
function_normalisation_VST_peu_de_genes<- function(counts,tableau_counts,filename = "plotDispEsts.png", width = 1200, height = 1000, res = 150){
  
  #Plus robuste, fonctionne même avec peu de gènes exprimés ou des données peu denses
  colnames(counts)=make.unique(colnames(counts),sep = ".") 
  metadata <- data.frame(
    row.names = colnames(counts),
    condition = colnames(counts) 
  )
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ 1)
  
  # Utilisation de varianceStabilizingTransformation à la place de vst()
  vst_data_condition <- varianceStabilizingTransformation(dds, blind = TRUE)
  vst_matrix_condition <- assay(vst_data_condition)
  
  # Affichage de la dispersion
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  
  png(filename, width = width, height = height, res = res)
  plotDispEsts(dds)
  dev.off()
  
  return(vst_matrix_condition)
}

function_normalisation_rlog = function(counts,tableau_counts,filename = "plotDispEsts.png", width = 1200, height = 1000, res = 150){
  colnames(counts)=make.unique(colnames(counts),sep = ".") 
  metadata <- data.frame(
    row.names = colnames(counts),
    condition =  colnames(counts) 
  )
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata,design = ~ 1)
  rlog_data <- rlog(dds, blind = TRUE)
  rlog_matrix <- assay(rlog_data)
  #dispersion
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  plotDispEsts(dds)
  
  return(rlog_matrix)
}

function_conversion_expression_binaire = function(counts){
  counts[counts > 0] = 1
  return(counts)
}

#######################################################################
NORMALISATION = function(counts, tableau_counts, normalisation = c("pas_de_normalisation","vst","vst_conditions", "vst_peu_de_genes","log_rpkm", "log_tpm", "rlog","binaire"),filename,liste_conditions_a_garder_etoile=NULL,longueurs_genes){
  normalisation <- match.arg(normalisation)
  if(normalisation=="pas_de_normalisation"){print("pas_de_normalisation")}
  if(normalisation=="log_rpkm"){counts = function_normalise_log_RPKM(counts,tableau_counts,longueurs_genes)}
  if(normalisation=="log_tpm"){counts = function_normalise_log_TPM(counts,tableau_counts,longueurs_genes)}
  if(normalisation=="rlog"){counts = function_normalisation_rlog(counts,tableau_counts)}
  if(normalisation=="vst"){counts = function_normalisation_VST(counts,tableau_counts,filename)}
  if(normalisation=="vst_peu_de_genes"){counts = function_normalisation_VST_peu_de_genes(counts,tableau_counts,filename)}
  if(normalisation=="vst_conditions"){counts = function_normalisation_VST_conditions(counts,tableau_counts,filename)}
  if(normalisation=="binaire"){counts = function_conversion_expression_binaire(counts)}
  
  return(counts)
}
#######################################################################

# 0.3.2) FONCTIONS BARPLOT ETC ..
function_barplot_nbr_genes_expr = function(counts, filename){
  noms_colonnes = colnames(counts)
  noms_genes = rownames(counts)
  liste_nbr_genes_exprimes=c()
  for(colonne in noms_colonnes){
    col = counts[,colonne]
    #
    nbr_genes_exprimes = sum(col != 0)  #nbr genes avec + de 1 reads
    names(nbr_genes_exprimes)=colonne
    liste_nbr_genes_exprimes = c(liste_nbr_genes_exprimes,nbr_genes_exprimes)
  }
  png(paste0(filename,".png"), width = 800, height = 600) 
  barplot(liste_nbr_genes_exprimes,main ="Nombre de genes exprimés",las = 2 )
  dev.off()  # ferme le fichier
  
  display_png(file = paste0(filename,".png"))
  
}

function_barplot_nbr_reads_moyen_par_gene_exprime = function(counts, filename){
  noms_colonnes = colnames(counts)
  noms_genes = rownames(counts)
  liste_nbr_reads_moyen_par_gene_exprime=c()
  liste_nbr_genes_exprimes=c()
  
  for(colonne in noms_colonnes){
    col = counts[,colonne]
    nbr_reads_total = sum(col)
    nbr_genes = length(col)
    nbr_reads_moyen_par_gene_all_genome = nbr_reads_total/nbr_genes
    
    #
    nbr_genes_exprimes = sum(col != 0)  #nbr genes avec + de 1 reads
    names(nbr_genes_exprimes)=colonne
    liste_nbr_genes_exprimes = c(liste_nbr_genes_exprimes,nbr_genes_exprimes)
    
    #
    nbr_reads_moyen_par_gene_exprime = nbr_reads_total/nbr_genes_exprimes
    names(nbr_reads_moyen_par_gene_exprime)=colonne
    liste_nbr_reads_moyen_par_gene_exprime = c(liste_nbr_reads_moyen_par_gene_exprime,nbr_reads_moyen_par_gene_exprime)
    
  }
  png(paste0(filename,".png"), width = 800, height = 600) 
  barplot(liste_nbr_reads_moyen_par_gene_exprime,main = "Nombre de reads moyens /genes exprimés",las = 2)
  dev.off()  # ferme le fichier
  
  display_png(file = paste0(filename,".png"))
  
}

function_barplot_genome_expression = function(counts, filename){
  noms_colonnes = colnames(counts)
  for(colonne in noms_colonnes){
    col = counts[,colonne]
    #
    png(paste0(filename,colonne,".png"), width = 800, height = 600)  # tu peux ajuster la taille
    barplot(col,main = colonne,las = 2)   #permet de garder l'ordre des genes et de voir des profils
    dev.off()  # ferme le fichier
    
    display_png(file = paste0(filename,colonne,".png"))
    
  }
}

function_boxplot_nbr_reads = function(counts,plot_titre,filename){
  
  png(filename, width = 800, height = 600)
  
  # Fonction pour calculer le N50 d'un vecteur numérique
  calc_N50 <- function(x) {
    x <- sort(x, decreasing = TRUE)
    cumsum_x <- cumsum(x)
    total <- sum(x)
    N50_val <- x[min(which(cumsum_x >= total / 2))]
    return(N50_val)
  }
  
  # Créons une figure avec boxplot et annotations
  boxplot(counts, 
          main = plot_titre,
          xlab = "Colonnes", ylab = "Valeurs",
          las = 2, col = "lightblue", border = "darkblue",
          outline = FALSE)  # Optionnel pour cacher les outliers
  
  # Ajout des statistiques pour chaque colonne
  for (i in seq_along(counts)) {
    col_data <- counts[[i]]
    
    # Calculs
    moy <- mean(col_data, na.rm = TRUE)
    med <- median(col_data, na.rm = TRUE)
    n50 <- calc_N50(col_data)
    quartiles <- quantile(col_data, probs = c(0.25, 0.75), na.rm = TRUE)
    
    # Position verticale pour texte
    ymax <- max(col_data, na.rm = TRUE)
    
    # Ajouter la moyenne (en rouge, triangle)
    points(i, moy, col = "red", pch = 17, cex = 1.5)
    
    # Ajouter la médiane (en bleu, cercle)
    points(i, med, col = "blue", pch = 19, cex = 1.5)
    
    # Ajouter le N50 (en vert, carré)
    points(i, n50, col = "darkgreen", pch = 15, cex = 1.5)
    
    # Ajouter les quartiles (en violet, croix)
    points(rep(i, 2), quartiles, col = "purple", pch = 4, cex = 1.5)
    
    # Ajouter texte descriptif au-dessus
    text(i, ymax, 
         labels = paste0("M:", round(moy, 1), 
                         "\nMd:", round(med, 1), 
                         "\nN50:", round(n50, 1)), 
         pos = 3, cex = 0.7)
  }
  
  # Légende
  legend("topright", legend = c("Moyenne", "Médiane", "N50", "Quartiles"),
         col = c("red", "blue", "darkgreen", "purple"),
         pch = c(17, 19, 15, 4), cex = 0.8)
  
  dev.off()
  
  display_png(file = filename)
  
}

# 0.3.2) FONCTIONS CLUSTERING HIERARCHIQUE

function_cluster_hierarchique = function(counts,filename){
  counts <- t(counts)  # maintenant [40 samples x 12000 gènes]
  d <- dist(counts, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")  # méthode ward.D2 souvent bonne pour clustering transcriptomique
  
  png(paste0(filename,".png"), width = 800, height = 600)  # tu peux ajuster la taille
  plot(hc, main = "Hierarchical Clustering des échantillons", xlab = "", sub = "", cex = 0.9)
  dev.off() 
  
  display_png(file = paste0(filename,".png"))
  
}


function_cluster_hierarchique_bootsrap = function(counts, filename){
  #DENDROGRAMME AVEC BOOTSTRAP POUR VERIFIER QUALTIE DATA ?
  
  result <- pvclust(counts, method.hclust = "ward.D2", method.dist = "euclidean", nboot = 100)
  png(paste0(filename,".png"), width = 800, height = 600)
  
  plot(result)
  pvrect(result, alpha=0.95)  # entoure les clusters avec p-value > 95%
  
  dev.off() 
  
  display_png(file = paste0(filename,".png"))
  
}

# 0.3.2) FONCTIONS ACP 2D,3D, UMAP , DIAG de VENN    :  #QUALITE SIGNATURE TRANSCRIPTIONNELLE DISCRIMINANTE

function_ACP <- function(counts, plot_titre = "ACP", filename = "ACP_plot.png", width = 1200, height = 1000, res = 150){  #ACP :
  data_t <- t(counts)
  pca <- prcomp(data_t, scale. = F)
  summary(pca)
  
  png(filename, width = width, height = height, res = res)
  plot(pca$x[,1:2], col=1:nrow(data_t), pch=19, main = "plot_titre")
  text(pca$x[,1:2], labels=rownames(pca$x), pos=3)
  dev.off()
  
  display_png(file = filename)
  
  # Voir les contributions aux composantes principales
  contributions <- pca$rotation
  
  # Quelle variable contribue le plus à PC1 ?
  #abs(contributions[,1])  # valeurs absolues pour l'importance
  # Classement décroissant
  liste_axe_gene_contribution = contributions[,1]
  liste_axe_gene_contribution = liste_axe_gene_contribution[liste_axe_gene_contribution>0.01]
  
  importance_PC1 <- sort(abs(liste_axe_gene_contribution), decreasing = TRUE)
  #print(importance_PC1)
  #boxplot(importance_PC1)
  write.csv(importance_PC1, paste0(filename,"contributions1.txt"))
  
  
  #importance_PC2 <- sort(abs(contributions[,2]), decreasing = TRUE)
  #print(importance_PC2)
  #write.csv(importance_PC2, paste0(filename,"contributions2.txt"))
  
  #importance_PC3 <- sort(abs(contributions[,3]), decreasing = TRUE)
  #print(importance_PC3)
  #write.csv(importance_PC3, paste0(filename,"contributions3.txt"))
  
}



function_UMAP = function(counts){
  
  #install.packages("uwot")
  
  # Charger le package
  
  counts <- t(counts)
  
  
  # UMAP (par défaut en 2D)
  umap_result <- umap(counts)
  
  # Résultat = une matrice avec les coordonnées projetées
  head(umap_result)
  noms=colnames(counts)
  noms=as.factor(noms)
  
  # Optionnel : visualisation
  plot(umap_result, col = as.numeric(noms), pch = 19,
       main = "Projection UMAP de iris")
  legend("topright", legend = levels(noms),
         col = 1:3, pch = 19) 
  
}


function_3D_PCA = function(counts, filename = "ACP_3D.html"){
  #PCA 3d : 
  
  counts <- t(counts)  
  pca <- prcomp(counts, scale. = FALSE)
  pca_df <- as.data.frame(pca$x)
  
  p <- plot_ly(data = pca_df, 
               
               x = ~PC1, y = ~PC2, z = ~PC3, 
               type = 'scatter3d', 
               mode = 'markers+text',
               text = rownames(pca_df),
               textposition = 'top center',
               marker = list(size = 5,
                             color = as.numeric(as.factor(rownames(pca_df))),
                             colorscale = 'Viridis'))
  
  htmlwidgets::saveWidget(p, filename, selfcontained = FALSE)
  
  browseURL(filename)
}


function_garde_500_genes__les_plus_variables = function(counts){
  gene_sd <- apply(counts, 1, sd)
  top_genes <- names(sort(gene_sd, decreasing = TRUE)[1:500])
  counts = counts[top_genes,]
  
  return(counts)
}


#FONCTION DIAG VENN
#sorte de diagramme de venn
function_genes_communs = function(counts){
  genes = rownames(counts)
  for(experience in colnames(counts)){
    for(gene in genes){
      nbr_de_reads = counts[gene,experience]
      if(nbr_de_reads>0){
        counts[gene,experience]=1 #gene
      }
    }
  }
  # Install if needed
  #install.packages("UpSetR")
  
  # Affichage du diagramme avec des options pour compacité
  upset(counts,
        sets = colnames(counts),
        keep.order = TRUE,
        sets.bar.color = "#56B4E9",
        order.by = "freq",   # trie les combinaisons par fréquence
        mb.ratio = c(0.6, 0.4),  # réduit la taille des barres du bas (main bar)
        text.scale = 1.2)     # ajuste la taille du texte pour lisibilité
}

# 0.3.2) FONCTIONS COEXPRESSION


function_CNN_genes = function(counts,df_genes_names=df_genes_names,R2 = 0.99,filename){
  
  # Supprimer les gènes constants avant la corrélation
  counts_filtered <- counts[apply(counts, 1, function(x) sd(x) != 0), ]
  
  # Puis recalculer la corrélation
  cor_matrix <- correlate(t(counts_filtered), method = "pearson")
  
  # 2. Transformer la matrice en format utilisable
  cor_df <- stretch(cor_matrix) # passer en format "long"
  
  # 3. Filtrer les fortes corrélations
  cor_df <- subset(cor_df, abs(r) >= R2)  
  
  # 4. Construire le graph
  g <- graph_from_data_frame(cor_df, directed = FALSE)
  
  
  # 5. Dessiner le graph
  # Définir les couleurs selon le signe de la corrélation
  edge_colors <- ifelse(E(g)$r < 0, "red", "green")
  
  # Générer le PNG
  png(paste0(filename,".",R2,".png"), width = 1200, height = 1000, res = 150)
  
  plot(g, 
       vertex.label = V(g)$name,
       vertex.label.cex = 0.7,
       vertex.size = 5,
       edge.width = abs(E(g)$r) * 5,  # épaisseur proportionnelle à la force
       edge.color = edge_colors,      # couleur selon le signe
       main = "CNN",
       layout = layout_with_fr)
  
  dev.off()
  
  
  display_png(file = paste0(filename,".",R2,".png"))
  
  
  
  # 1. Filtrer les gènes qui existent dans le graphe
  genes_cibles_valides <- intersect(df_genes_names, V(g)$name)
  
  # 2. Vérifier qu'au moins un gène est présent
  if (length(genes_cibles_valides) > 0) {
    # 3. Trouver les voisins directs (ordre 1) + inclure les gènes eux-mêmes
    voisins <- unlist(neighborhood(g, order = 1, nodes = genes_cibles_valides, mode = "all"))
    
    # 4. Extraire les noms des sommets à inclure dans le sous-graphe
    sommets_sousgraphe <- unique(V(g)[voisins]$name)
    
    # 5. Créer le sous-graphe
    g_sous <- induced_subgraph(g, vids = sommets_sousgraphe)
    
    # 6. Afficher ou enregistrer le sous-graphe
    png(filename = paste0(filename,".png"), width = 8000, height = 8000, res = 400)
    # Créer un vecteur de couleurs : bleu pour les gènes cibles, orange pour les autres
    couleurs_sommets <- ifelse(V(g_sous)$name %in% genes_cibles_valides, "blue", "orange")
    layout_graphopt <- layout_with_graphopt(g_sous, charge = 0.30, niter = 2000)
    plot(g_sous,
         vertex.label = V(g_sous)$name,
         vertex.label.cex = 0.8,
         vertex.size = 6,
         vertex.color = couleurs_sommets,  # <<<< couleurs ici
         edge.width = abs(E(g_sous)$r) * 5,
         main = "Gènes cibles (bleu) + voisins (orange)",
         layout = layout_graphopt)
    
    dev.off()
  } else {
    cat("Aucun des gènes cibles n'est présent dans le graphe.\n")
  }
  
  display_png(file = paste0(filename,".png"))
  
}

function_matrice_correlation = function(counts,plot_titre = "pheatmap", filename = "pheatmap.png", width = 1200, height = 1000, res = 150){
  
  
  cor_matrix <- correlate(t(counts), method = "pearson")
  cor_matrix = as.data.frame(cor_matrix,row.names = TRUE)
  rownames(cor_matrix) = cor_matrix$term
  cor_matrix = cor_matrix[,-1]
  cor_matrix[is.na(cor_matrix)] <- 0
  cor_matrix =as.matrix(cor_matrix)
  
  
  
  # Créer la matrice de labels
  number_matrix <- matrix(sprintf("%.1f", cor_matrix), 
                          nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
  
  # Mettre les "0.0" en texte blanc
  number_color <- matrix("black", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
  number_color[cor_matrix == 0] <- "white"
  
  
  # Affichage du heatmap
  ph = pheatmap(cor_matrix,
                main =plot_titre,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                display_numbers = number_matrix,
                number_color = number_color,
                fontsize_number = 10,
                col = colorRampPalette(c("white", "yellow", "red"))(100)
  )
  # Enregistre en PNG
  png(filename, width = width, height = height, res = res)
  grid::grid.draw(ph$gtable)  # Dessine le heatmap capturé
  dev.off()
  
  display_png(file = filename)
  
}

function_matrice_correlation_courbes = function(counts,plot_titre = "pheatmap", filename = "pheatmap.png", width = 1200, height = 1000, res = 150){
  
  
  #install.packages("GGally")     # à faire une seule fois
  
  # Affiche les nuages de points pour toutes les paires de variables
  
  
  p = ggpairs(
    t(counts),
    
    upper = list(continuous = wrap("cor", size = 3)),
    
    lower = list(
      continuous = wrap("smooth", 
                        method = "lm", 
                        se = FALSE,
                        color = "red", 
                        #fullrange = TRUE,
                        alpha = 0.8)
    ),
    
    diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
    title = plot_titre
  ) +
    theme_minimal()
  
  print(p)
  ggsave(filename = filename, plot = p,width = 12, height = 12, dpi = 300, units = "in")
  
  display_png(file = filename)
  
}



function_CNN = function(counts,R2 = 0.99,filename){
  
  # Supprimer les gènes constants avant la corrélation
  counts_filtered <- counts[apply(counts, 1, function(x) sd(x) != 0), ]
  
  # Puis recalculer la corrélation
  cor_matrix <- correlate(t(counts_filtered), method = "pearson")
  
  # 2. Transformer la matrice en format utilisable
  cor_df <- stretch(cor_matrix) # passer en format "long"
  
  # 4. Supprimer les auto-corrélations (x == y)
  cor_df <- cor_df %>% filter(cor_df$x != cor_df$y)
  
  # 5. Réordonner les paires pour ne garder qu’une direction (évite doublons)
  cor_df <- cor_df %>%
    mutate(pair_id = paste0(pmin(cor_df$x, cor_df$y), "_", pmax(cor_df$x, cor_df$y))) %>%
    distinct(pair_id, .keep_all = TRUE)
  
  # 3. Filtrer les fortes corrélations
  cor_df <- subset(cor_df, abs(r) >= R2)  
  
  # 4. Construire le graph
  g <- graph_from_data_frame(cor_df, directed = FALSE)
  
  # 5. Dessiner le graph
  # Définir les couleurs selon le signe de la corrélation
  edge_colors <- ifelse(E(g)$r < 0, "red", "green")
  
  # Générer le PNG
  png(paste0(filename,".",R2,".png"), width = 1200, height = 1000, res = 150)
  
  plot(g, 
       vertex.label = V(g)$name,
       vertex.label.cex = 0.7,
       vertex.size = 5,
       edge.width = abs(E(g)$r) * 5,  # épaisseur proportionnelle à la force
       edge.color = edge_colors,      # couleur selon le signe
       main = "CNN",
       layout = layout_with_fr)
  
  dev.off()
  
  display_png(file = paste0(filename,".",R2,".png"))
  
  
}

# 0.3.2) FONCTIONS ANALYSE EXPRESSION DIFFERENTIELLE

function_analyse_expression_differentielle = function(
    counts,
    condition_1,
    condition_2,
    repertoire,
    noms_genes)
{
  
  #change les noms des genes d'interets mais garde l'ensemble des 12000 genes
  chemin_tab_noms_genes = paste0(chemin_repertoire_output,"/",noms_genes,".txt")                            #repertoire_fichiers_txt    
  df_genes_names = read.table(chemin_tab_noms_genes, header=FALSE,  sep="\t", comment.char="#")
  print(df_genes_names[,2])
  counts = function_CONVERSION_NOM_DES_GENES(df_genes_names,counts)
  liste_de_genes = df_genes_names[,2]
  
  dir.create(repertoire)
  
  
  
  
  metadata <- data.frame(
    row.names = colnames(counts),
    condition = sapply(strsplit(colnames(counts), split = "\\."), `[`, 1)    # passe de "dt_fo.trucaenelever"  à "dt_fo"
  )
  print(metadata)
  
  
  #normalisation
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~condition)
  
  #AJUSTE GLM ET REDUCTION DISPERSIONS ET CALCUL DES LFC
  dds <- DESeq(dds)
  
  #The function DESeq runs the following functions in order:
  #1) dds <- estimateSizeFactors(dds)
  #2) dds <- estimateDispersions(dds)
  #3) dds <- nbinomWaldTest(dds)
  
  #plot dispersion vs mean
  filename = paste0(repertoire,"/dispersion.png")
  png(filename = filename, width = 1200, height = 1000, res = 150)
  plotDispEsts(dds)
  dev.off()
  display_png(file = filename)
  
  res <- results(dds)
  res <- results(dds, name = resultsNames(dds)[2])
  res <- results(dds, contrast=c("condition",condition_1,condition_2))
  resOrdered <- res[order(res$pvalue),]
  summary(res)
  sum(res$padj < 0.1, na.rm=TRUE)
  
  
  print("Visualisation des LFC avec le bruit, car les LFC augmente lorsque le comptage est faible, (heteroscedasticité)")
  filename = paste0(repertoire,"/plotMA_",resultsNames(dds)[2],".png")
  png(filename = filename, width = 1200, height = 1000, res = 150)
  plotMA(res, ylim=c(-2,2))    
  dev.off()
  display_png(file = filename)
  
  
  
  #SUPPRIME LE BRUIT DES CHANGEMENT DE REPLIS DES GENES A FAIBLE NOMBRE
  resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
  
  print("Visualisation des LFC avec la réduction du bruit, on réduit lorsque l'information du gene est faible (comptage faible, et ou dispersion elevé, et ou faible nbr de degres de libertés)")
  filename = paste0(repertoire,"/plotMA_LFC_",resultsNames(dds)[2],".png")
  png(filename = filename, width = 1200, height = 1000, res = 150)
  plotMA(resLFC, ylim=c(-2,2))
  dev.off()
  display_png(file = filename)
  
  
  #LISTE DE NOS GENES DIFF EXPR
  genes_in_res <- liste_de_genes[liste_de_genes %in% rownames(resLFC)]
  res_genes_of_interest <- resLFC[genes_in_res, ]
  de_genes <- res_genes_of_interest[which(res_genes_of_interest$padj < 0.05), ]
  
  print("les genes differentiellement exprimés entre les deux conditions parmis nos gens d'interets sont : ")
  print(de_genes)
  write.table(as.data.frame(de_genes), file=paste0(repertoire,"/",resultsNames(dds)[2],"genes_diff_expr.txt"),sep = "\t")
  
  print("les genes qui ne sont pas differentiellement exprimés entre les deux conditions parmis nos gens d'interets sont : ")
  NON_de_genes = res_genes_of_interest[which(res_genes_of_interest$padj >= 0.05), ]
  print(NON_de_genes)
  write.table(as.data.frame(NON_de_genes), file=paste0(repertoire,"/",resultsNames(dds)[2],"genes_non_diff_expr.txt"),sep = "\t")
  
  #Visualise l'expression des replicats des conditions
  
  for (gene in genes_in_res) {    
    gene_data <- resLFC[gene, ]
    
    # Extraction des valeurs avec 3 chiffres significatifs
    baseMean_val <- signif(gene_data$baseMean, 3)
    lfc_val      <- signif(gene_data$log2FoldChange, 3)
    lfcSE_val    <- signif(gene_data$lfcSE, 3)
    pval_val     <- signif(gene_data$pvalue, 3)
    padj_val     <- signif(gene_data$padj, 3)
    
    # Titre formaté
    titre <- paste0(
      "\n\n\n",
      "baseMean = ", baseMean_val, ", ",
      "log2FC = ", lfc_val, ", ",
      "lfcSE = ", lfcSE_val, "\n",
      "pvalue = ", pval_val, ", ",
      "padj = ", padj_val
    )
    
    # Génération du graphique avec titre
    png(filename = paste0(repertoire, "/", gene, ".png"), width = 1200, height = 1000, res = 150)
    par(las = 2)
    plotCounts(dds, gene = gene, intgroup = "condition")
    title(main = titre)
    dev.off()
    
    print(paste0(repertoire, "/", gene, ".png"))
  }
  
  
  
  
  write.csv(as.data.frame(resOrdered), file=paste0(repertoire,"/",resultsNames(dds)[2],".csv"))
  
  resOrdered_LFC <- resLFC[order(resLFC$pvalue),]
  
  mcols(res)$description
  filename = paste0(repertoire,"/boxplot_cooks_",resultsNames(dds)[2],".png")
  png(filename = filename, width = 1200, height = 1000, res = 150)
  boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
  dev.off()
  display_png(file = filename)
  
  #PLOT Nombre de genes differentiellement exprimés en fonction de ?
  filename = paste0(repertoire,"/LFC_vs_quantile_",resultsNames(dds)[2],".png")
  png(filename = filename, width = 1200, height = 1000, res = 150)
  plot(metadata(res)$filterNumRej, 
       type="b", ylab="number of rejections",
       xlab="quantiles of filter")
  lines(metadata(res)$lo.fit, col="red")
  abline(v=metadata(res)$filterTheta)
  dev.off()
  display_png(file = filename)
  
  
  #PLOT COOKS DISTANCE
  filename = paste0(repertoire,"/cook_distance_",resultsNames(dds)[2],".png")
  png(filename = filename, width = 1200, height = 1000, res = 150)
  W <- res$stat
  maxCooks <- apply(assays(dds)[["cooks"]],1,max)
  idx <- !is.na(W)
  plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic", 
       ylab="maximum Cook's distance per gene",
       ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
  m <- ncol(dds)
  p <- 3
  abline(h=qf(.99, p, m - p))
  dev.off()
  display_png(file = filename)
  
  
  #PLOT BARPLOT PASS DO NOT PASS
  filename = paste0(repertoire,"/pass_or_nor_pass_",resultsNames(dds)[2],".png")
  png(filename = filename,width = 1200, height = 1000, res = 150)
  use <- res$baseMean > metadata(res)$filterThreshold
  h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
  h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
  colori <- c(`do not pass`="khaki", `pass`="powderblue")
  barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
          col = colori, space = 0, main = "", ylab="frequency")
  text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
       adj = c(0.5,1.7), xpd=NA)
  legend("topright", fill=rev(colori), legend=rev(names(colori)))
  dev.off()
  display_png(file = filename)
  
  
  #recupere VST
  vsd <- vst(dds, blind=FALSE)
  
  
  #HEATMAP
  library("pheatmap")
  
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:20]
  df <- as.data.frame(colData(dds)["condition"])
  print(df)
  filename = paste0(repertoire,"/heatmap_",resultsNames(dds)[2],".png")
  ph = pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
                cluster_cols=FALSE, annotation_col=df)
  
  png(filename = filename,width = 1200, height = 1000, res = 150)
  grid::grid.draw(ph$gtable)  # Dessine le heatmap capturé
  dev.off()
  display_png(file = filename)
  
  #matrix distance
  sampleDists <- dist(t(assay(vsd)))
  library("RColorBrewer")
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  filename = paste0(repertoire,"/distance_matrice_",resultsNames(dds)[2],".png")
  ph = pheatmap(sampleDistMatrix,
                clustering_distance_rows=sampleDists,
                clustering_distance_cols=sampleDists,
                col=colors)
  png(filename = filename,width = 1200, height = 1000, res = 150)
  grid::grid.draw(ph$gtable)  # Dessine le heatmap capturé
  dev.off()
  display_png(file = filename)
  
  #PCA
  filename = paste0(repertoire,"/PCA_",resultsNames(dds)[2],".png")
  ph = plotPCA(vsd, intgroup="condition")
  
  png(filename = filename,width = 1200, height = 1000, res = 150)
  print(ph)    
  dev.off()
  display_png(file = filename)
  
}




# 1.3) PARAMETRES





ui <- fluidPage(
  titlePanel("Configuration des Paramètres"),
  
  # CSS pour améliorer l'apparence
  tags$head(
    tags$style(HTML("
      .form-control {
        font-size: 14px;
        height: 40px;
      }
      textarea.form-control {
        height: 80px !important;
        resize: vertical;
      }
      .help-text {
        color: #666;
        font-size: 12px;
        margin-bottom: 15px;
      }
      .config-section {
        background-color: #f8f9fa;
        padding: 20px;
        border-radius: 8px;
        margin-bottom: 25px;
        border: 1px solid #dee2e6;
      }
      .section-title {
        color: #495057;
        margin-bottom: 20px;
        border-bottom: 2px solid #007bff;
        padding-bottom: 10px;
      }
    "))
  ),
  
  # Section 1: Configuration des fichiers et normalisation
  div(class = "config-section",
      h3(class = "section-title", "📁 Configuration des Fichiers et Normalisation"),
      
      fluidRow(
        column(6,
               textInput("repertoire_fichiers_txt", 
                         label = "Répertoire des fichiers TXT :",
                         value = "", 
                         width = "100%"),
               div(class = "help-text", "Chemin vers le dossier contenant vos matrices (mettre un / à la fin)")
        ),
        column(6,
               textInput("tableau_counts", 
                         label = "Nom de la matrice d'alignements :",
                         value = "", 
                         width = "100%"),
               div(class = "help-text", "Nom de la matrice qui sera analysée (sans extension .txt)")
        )
      ),
      
      fluidRow(
        column(6,
               selectInput("normalisation", 
                           label = "Type de normalisation :",
                           choices = list(
                             "Pas de normalisation" = "pas_de_normalisation",
                             "VST" = "vst",
                             "VST conditions" = "vst_conditions", 
                             "VST peu de gènes" = "vst_peu_de_genes",
                             "Log RPKM" = "log_rpkm",
                             "Log TPM" = "log_tpm",
                             "RLOG" = "rlog",
                             "Binaire" = "binaire"
                           ),
                           selected = "pas_de_normalisation", 
                           width = "100%"),
               div(class = "help-text", "Choisir le type de normalisation à appliquer")
        ),
        column(6,
               textInput("noms_genes", 
                         label = "Nom du fichier de gènes :",
                         value = "", 
                         width = "100%"),
               div(class = "help-text", "Nom du fichier contenant la liste de gènes à étudier")
        )
      )
  ),
  
  # Section 2: Noms des gènes
  div(class = "config-section",
      h3(class = "section-title", "🧬 Noms des Gènes"),
      
      fluidRow(
        column(6,
               h4("Noms Anciens"),
               div(class = "help-text", "Séparez les noms par des espaces"),
               textAreaInput("noms_anciens", 
                             label = NULL,
                             value = "",
                             placeholder = "Exemple: LOC406124 Grd LCCH3 Amel_8916 GluCl",
                             width = "100%"),
               div(class = "help-text", "Rechercher dans le fichier 'liste_de_genes.txt' l'existence des noms des gènes")
        ),
        column(6,
               h4("Noms Nouveaux"), 
               div(class = "help-text", "Séparez les noms par des espaces"),
               textAreaInput("noms_nouveaux", 
                             label = NULL,
                             value = "",
                             placeholder = "Exemple: LOC406124_Rdl Grd LCCH3 Amel_8916 GluCl", 
                             width = "100%"),
               div(class = "help-text", "Nouveaux noms pour la conversion")
        )
      )
  ),
  
  # Bouton de validation
  fluidRow(
    column(12,
           div(style = "text-align: center; margin: 30px 0;",
               actionButton("valider", "✅ Valider et Générer le Code", 
                            class = "btn-primary btn-lg", 
                            style = "padding: 15px 30px; font-size: 18px;")
           )
    )
  ),
  
  # Zone de résultats
  conditionalPanel(
    condition = "output.show_results",
    div(class = "config-section",
        h3(class = "section-title", "📋 Code R Généré"),
        verbatimTextOutput("code_genere"),
        
        h4("📊 Aperçu des Gènes"),
        tableOutput("apercu_donnees")
    )
  )
)

server <- function(input, output, session) {
  
  # Variable pour contrôler l'affichage des résultats
  output$show_results <- reactive({
    !is.null(input$valider) && input$valider > 0
  })
  outputOptions(output, "show_results", suspendWhenHidden = FALSE)
  
  observeEvent(input$valider, {
    # Récupérer tous les paramètres
    repertoire_fichiers_txt <- trimws(input$repertoire_fichiers_txt)
    tableau_counts <- trimws(input$tableau_counts)
    normalisation <- input$normalisation
    noms_genes <- trimws(input$noms_genes)
    
    # Récupérer et nettoyer les entrées de gènes
    anciens_raw <- trimws(input$noms_anciens)
    nouveaux_raw <- trimws(input$noms_nouveaux)
    
    # Diviser par les espaces et nettoyer
    if(anciens_raw != "") {
      noms_anciens <- trimws(unlist(strsplit(anciens_raw, "\\s+")))
      noms_anciens <- noms_anciens[noms_anciens != ""]
    } else {
      noms_anciens <- character(0)
    }
    
    if(nouveaux_raw != "") {
      noms_nouveaux <- trimws(unlist(strsplit(nouveaux_raw, "\\s+")))
      noms_nouveaux <- noms_nouveaux[noms_nouveaux != ""]
    } else {
      noms_nouveaux <- character(0)
    }
    
    # Générer le code R complet avec commentaires
    code_complet <- paste0(
      'repertoire_fichiers_txt <- "', repertoire_fichiers_txt, '"\n',
      '# Chemin vers le dossier contenant vos matrices d\'alignements\n\n',
      
      'tableau_counts <- "', tableau_counts, '"\n',
      '# Nom de la matrice d\'alignements qui sera analysée (sans extension)\n\n',
      
      'normalisation <- "', normalisation, '"\n',
      '# Type de normalisation choisi\n\n',
      
      'noms_genes <- "', noms_genes, '"\n',
      '# Nom du fichier contenant la liste de gènes à étudier\n\n',
      
      if(length(noms_anciens) > 0) {
        paste0('nom_ancien <- c("', paste(noms_anciens, collapse = '", "'), '")\n',
               '# Rechercher dans le fichier "liste_de_genes.txt" l\'existence des noms des gènes\n\n')
      } else {
        'nom_ancien <- c()\n# Aucun nom ancien spécifié\n\n'
      },
      
      if(length(noms_nouveaux) > 0) {
        paste0('nom_nouveau <- c("', paste(noms_nouveaux, collapse = '", "'), '")\n',
               '# Nouveaux noms pour conversion\n')
      } else {
        'nom_nouveau <- c()\n# Aucun nom nouveau spécifié\n'
      }
    )
    
    # Afficher le code généré
    output$code_genere <- renderText({
      code_complet
    })
    
    # Créer un aperçu des données
    max_length <- max(length(noms_anciens), length(noms_nouveaux))
    if(max_length > 0) {
      anciens_extended <- c(noms_anciens, rep("", max_length - length(noms_anciens)))
      nouveaux_extended <- c(noms_nouveaux, rep("", max_length - length(noms_nouveaux)))
      
      apercu_df <- data.frame(
        "Index" = 1:max_length,
        "Nom Ancien" = anciens_extended,
        "Nom Nouveau" = nouveaux_extended,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    } else {
      apercu_df <- data.frame(
        "Index" = integer(0),
        "Nom Ancien" = character(0),
        "Nom Nouveau" = character(0),
        check.names = FALSE
      )
    }
    
    output$apercu_donnees <- renderTable({
      apercu_df
    }, rownames = FALSE)
    
    # Attendre 3 secondes puis fermer et sauvegarder
    Sys.sleep(3)
    
    # Sauvegarder toutes les variables dans l'environnement global
    repertoire_fichiers_txt <<- repertoire_fichiers_txt
    tableau_counts <<- tableau_counts
    normalisation <<- normalisation
    noms_genes <<- noms_genes
    nom_ancien <<- noms_anciens
    nom_nouveau <<- noms_nouveaux
    
    # Fermer l'application
    stopApp(list(
      repertoire_fichiers_txt = repertoire_fichiers_txt,
      tableau_counts = tableau_counts,
      normalisation = normalisation,
      noms_genes = noms_genes,
      nom_ancien = noms_anciens, 
      nom_nouveau = noms_nouveaux
    ))
  })
}

# Lancement de l'application
cat("🚀 Lancement de l'interface de configuration...\n")
resultats <- shinyApp(ui, server, options = list(launch.browser = TRUE))
print(resultats)
# Afficher les résultats
cat("\n✅ Variables créées :\n")
if(exists("repertoire_fichiers_txt")) cat("📁 repertoire_fichiers_txt:", repertoire_fichiers_txt, "\n")
if(exists("tableau_counts")) cat("📊 tableau_counts:", tableau_counts, "\n")  
if(exists("normalisation")) cat("🔧 normalisation:", normalisation, "\n")
if(exists("noms_genes")) cat("📝 noms_genes:", noms_genes, "\n")
if(exists("nom_ancien")) {
  cat("🧬 nom_ancien:\n")
  print(nom_ancien)
} 
if(exists("nom_nouveau")) {
  cat("🔄 nom_nouveau:\n")  
  print(nom_nouveau)
}




#

# 0.1) parametres à ne pas modifier

nom_tableau_counts = paste0(tableau_counts,".txt")
# Parametres Partie 2) facultatifs
df <- data.frame(
  nom_ancien =nom_ancien,    
  nom_nouveau =nom_nouveau         
)
liste_conditions_a_garder_etoile=NULL
nom_repertoire_output=""





#repertoire_fichiers_txt = paste0(dossier_travail,"/FICHIERS_TXT")
OUTPUTS = paste0(repertoire_fichiers_txt,"OUTPUT")
dir.create(OUTPUTS)

# 2.1) CREER DOSSIER OUTPUT

#Créer le répertoire output
if(nom_repertoire_output == "") {
  chemin_repertoire_output = paste0(OUTPUTS, "/", tableau_counts, "_", normalisation, "_", noms_genes)
  nom_repertoire_output = tableau_counts
} else {
  chemin_repertoire_output = paste0(OUTPUTS, "/", nom_repertoire_output)
}
dir.create(chemin_repertoire_output)


# 2.1) sauve la liste des genes
write.table(df, file =paste0(chemin_repertoire_output,"/",noms_genes,".txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
print(chemin_repertoire_output)

# 2.1) PREPARATION TABLEAU COUNTS
#recupere le tableau "counts"  (alignement vs 12000 genes)
counts = function_recupere_tab_counts_txt(tableau_counts)                                  #Tableau alignements VS tout les genes

#recupere longueurs des genes
longueurs_genes = recupere_longueurs_genes(counts)                                         #faire avant "function_selection_gene_dans_counts_AND_CONVERSION_NOM_DES_GENES"
#supprime les 1eres colonnes
counts = function_supprimer_x_premieres_colonnes(counts,5)                                 #supprime les premiers colonnes du tableau (longueur,chromosome, etc ..)
all_numeric_columns <- all(sapply(counts, is.numeric))
















library(shiny)

liste <- colnames(counts)

ui <- fluidPage(
  titlePanel("Compléter une liste"),
  
  # Ajout de CSS pour améliorer l'alignement
  tags$head(
    tags$style(HTML("
      .liste-item {
        font-size: 14px;
        font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
        margin-bottom: 5px;
        padding: 6px 12px;
        line-height: 1.42857143;
        height: 34px;
        display: flex;
        align-items: center;
        border: 1px solid transparent;
        background-color: #f9f9f9;
      }
      .shiny-input-container {
        margin-bottom: 5px;
      }
      .form-control {
        height: 34px;
        font-size: 14px;
      }
    "))
  ),
  
  fluidRow(
    column(6, 
           h4("Liste"),
           lapply(liste, function(x) 
             tags$div(class = "liste-item", x)
           )
    ),
    column(6, 
           h4("Entrées"),
           lapply(seq_along(liste), function(i) 
             textInput(paste0("entry", i), label = NULL, width = "100%")
           )
    )
  ),
  
  br(),
  actionButton("valider", "Valider", class = "btn-primary"),
  br(), br(),
  tableOutput("resultat")
)

server <- function(input, output, session) {
  observeEvent(input$valider, {
    # Récupération des saisies
    saisies <- sapply(seq_along(liste), function(i) {
      val <- input[[paste0("entry", i)]]
      if(is.null(val) || val == "") "" else val
    })
    
    # Affichage du résultat
    output$resultat <- renderTable({
      data.frame(Liste = liste, Saisi = saisies, stringsAsFactors = FALSE)
    })
    
    # Fermeture automatique de la fenêtre après un court délai
    Sys.sleep(1)  # Pause d'1 seconde pour voir le résultat
    
    # Sauvegarder les résultats dans l'environnement global
    resultats <<- data.frame(Liste = liste, Saisi = saisies, stringsAsFactors = FALSE)
    
    # Fermer la fenêtre
    stopApp(resultats)
  })
}

# Lancement de l'application
resultats_saisie <- shinyApp(ui, server, options = list(launch.browser = TRUE))

# Les résultats sont maintenant disponibles dans la variable 'resultats_saisie'
print(resultats_saisie)

df <- data.frame(
  nom_actuel =  colnames(counts), 
  nom_nouveau = resultats$Saisi
)
tableau_counts_names = paste0(tableau_counts,"_names.txt")


file = paste0(chemin_repertoire_output,"/",tableau_counts_names)
if (!file.exists(file)){
  write.table(df, file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}





# 2.1) change le nom des alignements
counts = function_change_NOMS_ALIGNEMENTS(counts,tableau_counts)                           #change le nom des alignements (colonnes) par cexu du fichiers *names.txt
#alignements que l'on veut garder
counts = function_conditions_a_garder_dans_counts(counts,liste_conditions_a_garder_etoile) #pour recuperer uniquement certains alignements (colonnes)


# 2.2) NORMALISATION 
#NORMALISATION de counts (a faire sur le plus grands jeux de donnée, tableau entiers avec 12000 genes)
print("normalisation :")
counts_norm = NORMALISATION(counts = counts,tableau_counts = tableau_counts,normalisation = normalisation,filename=paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_dispEsts.png"),liste_conditions_a_garder_etoile=liste_conditions_a_garder_etoile,longueurs_genes=longueurs_genes)
#SAUVEGARDE DU TABLEAU COUNTS Normalisé
#sauvegarde le tableau counts normalise
head(counts_norm)
write.table(counts_norm,paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,".","all_genes.counts.txt"),row.names=TRUE,col.names = TRUE, quote = FALSE, sep = "\t")


# 2.2) SELECTION GENES
#selectionne une sous partie du tableau normalise avec les genes d'interets 
counts_norm_genes = function_selection_gene_dans_counts_AND_CONVERSION_NOM_DES_GENES(noms_genes,counts_norm)
#selectionne une sous partie du tableau normalise avec les genes d'interets 
counts_genes = function_selection_gene_dans_counts_AND_CONVERSION_NOM_DES_GENES(noms_genes,counts)
#sauvegarde le sous tableau normalise avec nos genes
write.table(counts_norm_genes,paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,".",noms_genes,".counts.txt"),row.names=TRUE,col.names = TRUE, quote = FALSE, sep = "\t")


# =============================================================================
# Interface Shiny pour analyses RNA-seq
# Boutons interactifs pour chaque type d'analyse
# =============================================================================


# =============================================================================
# Interface utilisateur (UI)
# =============================================================================
ui <- dashboardPage(
  
  # En-tête
  dashboardHeader(title = "Analyses RNA-seq Dashboard"),
  
  # Barre latérale vide (optionnel)
  dashboardSidebar(disable = TRUE),
  
  # Corps principal
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side {
          background-color: #f4f4f4;
        }
        .btn-analysis {
          margin: 5px;
          width: 300px;
          height: 60px;
          font-size: 14px;
          font-weight: bold;
        }
        .analysis-section {
          background: white;
          padding: 20px;
          margin: 10px;
          border-radius: 10px;
          box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
      "))
    ),
    
    fluidRow(
      column(12,
             h2("🧬 Interface d'Analyses RNA-seq", style = "text-align: center; color: #2c3e50;"),
             hr()
      )
    ),
    
    # Section Heatmaps
    fluidRow(
      column(12,
             div(class = "analysis-section",
                 h3("🔥 Analyses Heatmap", style = "color: #e74c3c;"),
                 fluidRow(
                   column(6,
                          actionButton("heatmap_norm", 
                                       "Heatmap Normalisée", 
                                       class = "btn btn-danger btn-analysis",
                                       icon = icon("fire"))
                   ),
                   column(6,
                          actionButton("heatmap_raw", 
                                       "Heatmap Raw", 
                                       class = "btn btn-warning btn-analysis",
                                       icon = icon("chart-area"))
                   )
                 ),
                 fluidRow(
                   column(6,
                          actionButton("pheatmap_norm", 
                                       "PHeatmap Normalisée", 
                                       class = "btn btn-danger btn-analysis",
                                       icon = icon("th"))
                   ),
                   column(6,
                          actionButton("pheatmap_raw", 
                                       "PHeatmap Raw", 
                                       class = "btn btn-warning btn-analysis",
                                       icon = icon("border-all"))
                   )
                 )
             )
      )
    ),
    
    # Section Boxplots et Barplots
    fluidRow(
      column(12,
             div(class = "analysis-section",
                 h3("📊 Analyses Boxplot & Barplot", style = "color: #3498db;"),
                 fluidRow(
                   column(6,
                          actionButton("boxplot_raw", 
                                       "Boxplot Reads Raw", 
                                       class = "btn btn-primary btn-analysis",
                                       icon = icon("chart-simple"))
                   ),
                   column(6,
                          actionButton("boxplot_norm", 
                                       "Boxplot Reads Normalisées", 
                                       class = "btn btn-info btn-analysis",
                                       icon = icon("chart-line"))
                   )
                 ),
                 fluidRow(
                   column(4,
                          actionButton("barplot_reads_gene", 
                                       "Barplot Reads/Gène", 
                                       class = "btn btn-primary btn-analysis",
                                       icon = icon("chart-column"))
                   ),
                   column(4,
                          actionButton("barplot_genes_expr", 
                                       "Barplot Gènes Exprimés", 
                                       class = "btn btn-info btn-analysis",
                                       icon = icon("dna"))
                   ),
                   column(4,
                          actionButton("barplot_genome_expr", 
                                       "Barplot Expression Génome", 
                                       class = "btn btn-primary btn-analysis",
                                       icon = icon("globe"))
                   )
                 )
             )
      )
    ),
    
    # Section Analyses multivariées
    fluidRow(
      column(12,
             div(class = "analysis-section",
                 h3("📈 Analyses Multivariées", style = "color: #27ae60;"),
                 fluidRow(
                   column(4,
                          actionButton("genes_communs", 
                                       "Gènes Communs (Venn)", 
                                       class = "btn btn-success btn-analysis",
                                       icon = icon("circle-nodes"))
                   ),
                   column(4,
                          actionButton("acp_2d", 
                                       "ACP 2D", 
                                       class = "btn btn-success btn-analysis",
                                       icon = icon("project-diagram"))
                   ),
                   column(4,
                          actionButton("acp_3d", 
                                       "ACP 3D", 
                                       class = "btn btn-success btn-analysis",
                                       icon = icon("cube"))
                   )
                 ),
                 fluidRow(
                   column(6,
                          actionButton("clustering", 
                                       "Clustering Hiérarchique", 
                                       class = "btn btn-success btn-analysis",
                                       icon = icon("sitemap"))
                   )
                 )
             )
      )
    ),
    
    # Section Status et Log
    fluidRow(
      column(12,
             div(class = "analysis-section",
                 h3("📋 Statut des Analyses", style = "color: #8e44ad;"),
                 verbatimTextOutput("analysis_log"),
                 hr(),
                 h4("📁 Fichiers générés:"),
                 DT::dataTableOutput("files_table")
             )
      )
    )
  )
)

# =============================================================================
# Logique serveur
# =============================================================================
server <- function(input, output, session) {
  
  # Variables réactives pour le log et les fichiers
  analysis_log <- reactiveVal("🟢 Interface prête. Cliquez sur un bouton pour lancer une analyse.")
  generated_files <- reactiveVal(data.frame(
    Fichier = character(0),
    Type = character(0),
    Timestamp = character(0),
    Statut = character(0)
  ))
  
  # Fonction helper pour ajouter au log
  add_to_log <- function(message) {
    current_time <- Sys.time()
    new_log <- paste0("[", format(current_time, "%H:%M:%S"), "] ", message)
    analysis_log(paste(analysis_log(), new_log, sep = "\n"))
  }
  
  # Fonction helper pour ajouter un fichier généré
  add_file <- function(filename, type, status = "✅ Succès") {
    current_files <- generated_files()
    new_file <- data.frame(
      Fichier = basename(filename),
      Type = type,
      Timestamp = format(Sys.time(), "%H:%M:%S"),
      Statut = status
    )
    generated_files(rbind(current_files, new_file))
  }
  
  # =============================================================================
  # Gestionnaires d'événements pour chaque bouton
  # =============================================================================
  
  # Heatmaps
  observeEvent(input$heatmap_norm, {
    add_to_log("🔥 Génération heatmap normalisée...")
    tryCatch({
      function_heatmap(counts = counts_norm_genes,
                       plot_titre = paste0(nom_repertoire_output,".",normalisation,".",noms_genes),
                       filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,".",noms_genes,".heatmap.png"))
      add_to_log("✅ Heatmap normalisée générée avec succès!")
      add_file(paste0(nom_repertoire_output,".",normalisation,".",noms_genes,".heatmap.png"), "Heatmap Normalisée")
    }, error = function(e) {
      add_to_log(paste("❌ Erreur heatmap normalisée:", e$message))
      add_file("heatmap_norm.png", "Heatmap Normalisée", "❌ Erreur")
    })
  })
  
  observeEvent(input$heatmap_raw, {
    add_to_log("🔥 Génération heatmap raw...")
    tryCatch({
      function_heatmap(counts = counts_genes,
                       plot_titre = paste0(nom_repertoire_output,".",noms_genes),
                       filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",noms_genes,".heatmap.png"))
      add_to_log("✅ Heatmap raw générée avec succès!")
      add_file(paste0(nom_repertoire_output,".",noms_genes,".heatmap.png"), "Heatmap Raw")
    }, error = function(e) {
      add_to_log(paste("❌ Erreur heatmap raw:", e$message))
      add_file("heatmap_raw.png", "Heatmap Raw", "❌ Erreur")
    })
  })
  
  observeEvent(input$pheatmap_norm, {
    add_to_log("🔥 Génération pheatmap normalisée...")
    tryCatch({
      counts_norm_genes_sorted <- counts_norm_genes[order(rownames(counts_norm_genes)), order(colnames(counts_norm_genes))]
      function_pheatmap(counts = counts_norm_genes_sorted,
                        plot_titre = paste0(nom_repertoire_output,".",normalisation,".",noms_genes),
                        filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,".",noms_genes,".pheatmap.png"))
      add_to_log("✅ PHeatmap normalisée générée avec succès!")
      add_file(paste0(nom_repertoire_output,".",normalisation,".",noms_genes,".pheatmap.png"), "PHeatmap Normalisée")
    }, error = function(e) {
      add_to_log(paste("❌ Erreur pheatmap normalisée:", e$message))
      add_file("pheatmap_norm.png", "PHeatmap Normalisée", "❌ Erreur")
    })
  })
  
  observeEvent(input$pheatmap_raw, {
    add_to_log("🔥 Génération pheatmap raw...")
    tryCatch({
      counts_genes_sorted <- counts_genes[order(rownames(counts_genes)), order(colnames(counts_genes))]
      function_pheatmap(counts = counts_genes_sorted,
                        plot_titre = paste0(nom_repertoire_output,".",noms_genes),
                        filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",noms_genes,".pheatmap.png"))
      add_to_log("✅ PHeatmap raw générée avec succès!")
      add_file(paste0(nom_repertoire_output,".",noms_genes,".pheatmap.png"), "PHeatmap Raw")
    }, error = function(e) {
      add_to_log(paste("❌ Erreur pheatmap raw:", e$message))
      add_file("pheatmap_raw.png", "PHeatmap Raw", "❌ Erreur")
    })
  })
  
  # Boxplots
  observeEvent(input$boxplot_raw, {
    add_to_log("📊 Génération boxplot reads raw...")
    tryCatch({
      function_boxplot_nbr_reads(counts,
                                 plot_titre = nom_repertoire_output,
                                 filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".boxplot_nbr_reads.png"))
      add_to_log("✅ Boxplot reads raw généré avec succès!")
      add_file(paste0(nom_repertoire_output,".boxplot_nbr_reads.png"), "Boxplot Raw")
    }, error = function(e) {
      add_to_log(paste("❌ Erreur boxplot raw:", e$message))
      add_file("boxplot_raw.png", "Boxplot Raw", "❌ Erreur")
    })
  })
  
  observeEvent(input$boxplot_norm, {
    add_to_log("📊 Génération boxplot reads normalisées...")
    tryCatch({
      function_boxplot_nbr_reads(counts_norm,
                                 plot_titre = nom_repertoire_output,
                                 filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,".boxplot_nbr_reads.png"))
      add_to_log("✅ Boxplot reads normalisées généré avec succès!")
      add_file(paste0(nom_repertoire_output,".",normalisation,".boxplot_nbr_reads.png"), "Boxplot Normalisé")
    }, error = function(e) {
      add_to_log(paste("❌ Erreur boxplot normalisé:", e$message))
      add_file("boxplot_norm.png", "Boxplot Normalisé", "❌ Erreur")
    })
  })
  
  # Barplots
  observeEvent(input$barplot_reads_gene, {
    add_to_log("📊 Génération barplot reads/gène...")
    tryCatch({
      function_barplot_nbr_reads_moyen_par_gene_exprime(counts = counts_norm,
                                                        filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_all_genes.barplot_nbr_reads_moyen_par_gene_exprime"))
      add_to_log("✅ Barplot reads/gène généré avec succès!")
      add_file(paste0(nom_repertoire_output,".",normalisation,"_all_genes.barplot_nbr_reads_moyen_par_gene_exprime"), "Barplot Reads/Gène")
    }, error = function(e) {
      add_to_log(paste("❌ Erreur barplot reads/gène:", e$message))
      add_file("barplot_reads_gene", "Barplot Reads/Gène", "❌ Erreur")
    })
  })
  
  observeEvent(input$barplot_genes_expr, {
    add_to_log("📊 Génération barplot gènes exprimés...")
    tryCatch({
      function_barplot_nbr_genes_expr(counts = counts,
                                      filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,"_all_genes.barplot_nbr_genes_expr"))
      add_to_log("✅ Barplot gènes exprimés généré avec succès!")
      add_file(paste0(nom_repertoire_output,"_all_genes.barplot_nbr_genes_expr"), "Barplot Gènes Exprimés")
    }, error = function(e) {
      add_to_log(paste("❌ Erreur barplot gènes exprimés:", e$message))
      add_file("barplot_genes_expr", "Barplot Gènes Exprimés", "❌ Erreur")
    })
  })
  
  observeEvent(input$barplot_genome_expr, {
    add_to_log("📊 Génération barplot expression génome...")
    tryCatch({
      dir.create(paste0(chemin_repertoire_output,"/barplot_genome_expression_",normalisation), showWarnings = FALSE)
      function_barplot_genome_expression(counts = counts_norm,
                                         filename = paste0(chemin_repertoire_output,"/barplot_genome_expression_",normalisation,"/",nom_repertoire_output,".",normalisation,"_barplot_genome_expression"))
      add_to_log("✅ Barplot expression génome généré avec succès!")
      add_file(paste0(nom_repertoire_output,".",normalisation,"_barplot_genome_expression"), "Barplot Expression Génome")
    }, error = function(e) {
      add_to_log(paste("❌ Erreur barplot expression génome:", e$message))
      add_file("barplot_genome_expr", "Barplot Expression Génome", "❌ Erreur")
    })
  })
  
  # Analyses multivariées
  observeEvent(input$genes_communs, {
    add_to_log("📈 Analyse gènes communs (Venn)...")
    tryCatch({
      function_genes_communs(counts)
      add_to_log("✅ Analyse gènes communs terminée avec succès!")
      add_file("genes_communs_venn", "Diagramme Venn")
    }, error = function(e) {
      add_to_log(paste("❌ Erreur gènes communs:", e$message))
      add_file("genes_communs_venn", "Diagramme Venn", "❌ Erreur")
    })
  })
  
  observeEvent(input$acp_2d, {
    add_to_log("📈 Génération ACP 2D...")
    tryCatch({
      function_ACP(counts = counts_norm,
                   plot_titre = paste0(nom_repertoire_output,".",normalisation,"_all_genes"), 
                   filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_all_genes.acp.png"))
      add_to_log("✅ ACP 2D générée avec succès!")
      add_file(paste0(nom_repertoire_output,".",normalisation,"_all_genes.acp.png"), "ACP 2D")
    }, error = function(e) {
      add_to_log(paste("❌ Erreur ACP 2D:", e$message))
      add_file("acp_2d.png", "ACP 2D", "❌ Erreur")
    })
  })
  
  observeEvent(input$acp_3d, {
    add_to_log("📈 Génération ACP 3D...")
    tryCatch({
      function_3D_PCA(counts = counts_norm,
                      filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_all_genes.acp3D.html"))
      add_to_log("✅ ACP 3D générée avec succès!")
      add_file(paste0(nom_repertoire_output,".",normalisation,"_all_genes.acp3D.html"), "ACP 3D")
    }, error = function(e) {
      add_to_log(paste("❌ Erreur ACP 3D:", e$message))
      add_file("acp_3d.html", "ACP 3D", "❌ Erreur")
    })
  })
  
  observeEvent(input$clustering, {
    add_to_log("📈 Génération clustering hiérarchique...")
    tryCatch({
      function_cluster_hierarchique(counts = counts_norm,
                                    filename = paste0(chemin_repertoire_output,"/",nom_repertoire_output,".",normalisation,"_all_genes.dendrogramme"))
      add_to_log("✅ Clustering hiérarchique généré avec succès!")
      add_file(paste0(nom_repertoire_output,".",normalisation,"_all_genes.dendrogramme"), "Clustering Hiérarchique")
    }, error = function(e) {
      add_to_log(paste("❌ Erreur clustering:", e$message))
      add_file("clustering", "Clustering Hiérarchique", "❌ Erreur")
    })
  })
  
  # =============================================================================
  # Outputs réactifs
  # =============================================================================
  output$analysis_log <- renderText({
    analysis_log()
  })
  
  output$files_table <- DT::renderDataTable({
    generated_files()
  }, options = list(pageLength = 5, dom = 't'))
}

# =============================================================================
# Lancement de l'application
# =============================================================================
launch_analysis_dashboard <- function() {
  cat("🚀 Lancement de l'interface Shiny RNA-seq...\n")
  cat("📁 Vérifiez que vos variables sont chargées:\n")
  cat("   - counts, counts_norm, counts_genes, counts_norm_genes\n")
  cat("   - nom_repertoire_output, chemin_repertoire_output\n")
  cat("   - normalisation, noms_genes\n")
  cat("   - Toutes vos functions (function_heatmap, etc.)\n")
  cat("\n🌐 L'interface va s'ouvrir dans votre navigateur...\n")
  
  shinyApp(ui = ui, server = server)
}

# Pour lancer l'application, utilisez:
launch_analysis_dashboard()
