# Daniel Lauer
# 21 June 2019 - 8 June 2023

# Paper title: Disruption of trait-environment relationships in African megafauna occurred in the middle Pleistocene

###########################################################################################################################
###########################################################################################################################

# PART I: INTEGRATING BODY MASS DATA INTO FAITH ET AL. DATASETS, AND DATA CLEANING

# Set the working directory to the folder containing "Phylacine" and "PanTHERIA" (and other data sources):

setwd('./Data')

# Read in the subsections of "Phylacine" and "PanTHERIA" containing taxonomic information and body masses of species:

Phylacine <- read.csv('Phylacine_Traits.csv')
Phylacine <- Phylacine[, c(which(colnames(Phylacine) == 'Binomial.1.2'):which(colnames(Phylacine) == 
  'Species.1.2'), which(colnames(Phylacine) == 'Mass.g'))]

PanTHERIA <- read.table('PanTHERIA.txt', header = TRUE, sep = '\t')
PanTHERIA <- PanTHERIA[, c(1:5, 7)]

# Process the two dataframes to prepare them for integration into the Faith et al. datasets:

  # Rename the columns to make the names simpler:
  
  colnames(PanTHERIA) <- c('Order', 'Family', 'Genus', 'Species', 'Binomial', 'Mass_Kg_PanTHERIA')
  colnames(Phylacine) <- c('Binomial', 'Order', 'Family', 'Genus', 'Species', 'Mass_Kg_Phylacine')
  
  # Format "Phylacine" to match the format of "PanTHERIA":
  
  Phylacine$Binomial <- gsub('_', ' ', Phylacine$Binomial)
  Phylacine <- Phylacine[, c(2:5, 1, ncol(Phylacine))]

  # Remove all rows from "PanTHERIA" with a missing body mass value:
  
  PanTHERIA <- PanTHERIA[!PanTHERIA$Mass_Kg_PanTHERIA == -999, ]
  
  # Convert body masses from grams to kilograms (Faith et al. data uses kilograms):
  
  PanTHERIA[, ncol(PanTHERIA)] <- PanTHERIA[, ncol(PanTHERIA)] / 1000
  Phylacine[, ncol(Phylacine)] <- Phylacine[, ncol(Phylacine)] / 1000

# Read in the Faith et al. datasets containing fossil and modern herbivore information:
  
Modern_herbivores <- read.csv('Modern_Herbivores.csv', skip = 1) # Skip first row, which contains irrelevant info.

Fossil_herbivores <- read.csv('Fossil_Herbivores_ml_ras_Updated.csv', skip = 1) # Skip first row here too for same reason.
Fossil_herbivores <- Fossil_herbivores[1:349, ] # Remove empty rows.
  
# Process both the "Fossil_herbivores" and "Modern_herbivores" datasets:

  # Update the levels of the "Size.Class" column in "Fossil_herbivores":

  Fossil_herbivores$Size.Class <- droplevels(Fossil_herbivores$Size.Class, '')

  # Remove reference columns from "Fossil_herbivores":

  Fossil_herbivores <- Fossil_herbivores[, -grep('Reference', colnames(Fossil_herbivores))]
  
  # Rename the "Species" and "Mass" columns from both datasets to names that make data processing below easier:
  
  Datasets_newnames <- lapply(list(Fossil_herbivores, Modern_herbivores), function(Dataset) { # For each dataset...
    colnames(Dataset)[grep('Species', colnames(Dataset))] <- 'Binomial' # Change one column's name.
    colnames(Dataset)[grep('Mass', colnames(Dataset))] <- 'Mass_Kg_Herbivores' # Change another column's name.
    return(Dataset) }) # Return the dataset with the updated column names.
  
  Fossil_herbivores <- Datasets_newnames[[1]] # Assign "Fossil_herbivores" to equal the first modified dataset.
  Modern_herbivores <- Datasets_newnames[[2]] # Assign "Modern_herbivores" to equal the second modified dataset.
  
  # Handle taxonomic mismatches and uncertainties within the datasets, to avoid erroneous merging with "Phylacine"/
  # "PanTHERIA" and to ensure that all subsequent analyses are accurate, clean, and consistent with Faith et al.:
  
    # Count how many words are in each "Binomial" entry (should ideally be 2 words - genus, species) in each dataset. This 
    # requires converting the "Binomial" column from a factor to a string/character:
  
    library(stringr) # For using the functions "str_count" and "str_length" below.
    Binomial_wordcounts <- list() # List to hold the word counts of the "Binomial" column entries in each dataset.
    Counter <- 1 # Counter to refer to the datasets in the "for" loop below.
    
    for (Dataset in list(Fossil_herbivores, Modern_herbivores, PanTHERIA, Phylacine)) { # For each dataset...
      Binomial_wordcounts[[Counter]] <- str_count(as.character(Dataset$Binomial), ' ') + 1 # No. of spaces + 1 = no. words.
      print(all(Binomial_wordcounts[[Counter]] == 2)) # Print out if all "Binomial" entries have exactly two words.
      Counter <- Counter + 1 } # Update the counter, such that it can refer to the next dataset in the next loop iteration.
    
    # Remove variables that are no longer necessary to keep (every time "rm" is used hereafter, it is for this reason):
    
    rm(Dataset, Counter)
  
    # The "Fossil_herbivores" dataset contains a number of taxonomic "Binomial" entries with more than two words. These
    # extra words pertain to important paleontological code. Thus, reduce the dataset's entries to two words by removing
    # small words and phrases (like "cf", "aff", etc.) from the "Binomial" column, BUT while keeping a record of the
    # original entries in a column called "Binomial_Original":
    
    library(tibble) # For using the function "add_column" below.
    Fossil_herbivores <- add_column(Fossil_herbivores, Binomial_Original=Fossil_herbivores$Binomial, .before = 'Binomial')
    
    Removed_phrases <- c('\\?', 'cf ', 'aff ', 'et ', ' nov', ' indet', ' KH2', ' 2', ' 3', ' 4', ' A', ' B', ' D')
    for (Phrase in Removed_phrases) { Fossil_herbivores$Binomial <- gsub(Phrase, '', Fossil_herbivores$Binomial) }
    rm(Phrase)
    
    # Remove all rows from "Fossil_herbivores" whose "Binomial_Original" column indicates that at the genus level, the
    # taxonomic identity of the species in a given row is uncertain. To do this, convert "Binomial_Original" from a factor
    # to a character:
    
    Fossil_herbivores$Binomial_Original <- as.character(Fossil_herbivores$Binomial_Original)
    
    Removed_phrases_genuslevel <- c('cf', 'Gen', '?', 'Alcelaphini') # Phrases denoting genus-level uncertainty.
    for (Phrase in Removed_phrases_genuslevel) { # For each phrase...
      Fossil_herbivores <- Fossil_herbivores[startsWith(Fossil_herbivores$Binomial_Original, Phrase) == FALSE, ] }
    rm(Phrase)

    # To avoid double-counting certain "Fossil_herbivores" taxa in the analyses that follow, merge rows that contain the
    # same taxonomic name in the "Binomial" column, UNLESS an original name contains " aff " (which Faith et al. mentions
    # refers to a distinct taxonomic group). This double-counting occurs in cases where one "Binomial_Original" entry 
    # contains " cf " and the other doesn't, so handle this accordingly:
      
      # Reset the row names of "Fossil_herbivores" to better visualize/understand the manipulations below:
    
      rownames(Fossil_herbivores) <- seq(nrow(Fossil_herbivores))
    
      # Find dataset row indices that contain " cf " in their "Binomial_Original" entries:
        
      Doublecount_rows_1 <- grep(' cf ', Fossil_herbivores$Binomial_Original)
      
      # Find the indices within "Doublecount_rows_1" that refer to each row in "Fossil_herbivores" that differs from the
      # row immediately below it in the "Binomial_Original" column by only " cf ", and create a variable called
      # "Doublecount_rows_below" with these subsetted indices accordingly:
        
      Doublecount_indices_below <- which(Fossil_herbivores$Binomial[Doublecount_rows_1] == Fossil_herbivores$Binomial[
        Doublecount_rows_1 + 1])
      Doublecount_rows_below <- Doublecount_rows_1[Doublecount_indices_below]
      
      # For each row in "Doublecount_rows_below", update its occurrence data to equal "1" if the row immediately below it
      # (its double-count analog) is equal to "1":
      
      Fossil_herbivores[Doublecount_rows_below, ][Fossil_herbivores[Doublecount_rows_below + 1, ] == 1] <- 1
      
      # Now that occurrence data has been updated in each row in "Doublecount_rows_below", remove each double-count analog
      # row below those rows:
      
      Fossil_herbivores <- Fossil_herbivores[-(Doublecount_rows_below + 1), ]
      
      # Repeat the above procedure, this time addressing rows "above" instead of "below". Make sure to subset
      # "Doublecount_rows_2" such that indexing a row index of "0" is avoided when calculating ("Doublecount_rows_2" - 1),
      # and make sure to exclude cases that involve the phrase " aff ":
      
      rownames(Fossil_herbivores) <- seq(nrow(Fossil_herbivores))

      Doublecount_rows_2 <- grep(' cf ', Fossil_herbivores$Binomial_Original)
      Doublecount_rows_2 <- Doublecount_rows_2[2:length(Doublecount_rows_2)] # Perform subsetting here.
      
      Doublecount_indices_above <- which(Fossil_herbivores$Binomial[Doublecount_rows_2] == Fossil_herbivores$Binomial[
        Doublecount_rows_2 - 1] & grepl(' aff ', Fossil_herbivores$Binomial_Original[Doublecount_rows_2 - 1]) == FALSE)
      Doublecount_rows_above <- Doublecount_rows_2[Doublecount_indices_above]
      
      Fossil_herbivores[Doublecount_rows_above, ][Fossil_herbivores[Doublecount_rows_above - 1, ] == 1] <- 1
      
      Fossil_herbivores <- Fossil_herbivores[-(Doublecount_rows_above - 1), ]
        
    # To avoid double-counting of species that are only recognized as "sp" in the "Binomial_Original" column of
    # "Fossil_herbivores", remove instances of their occurrences at sites where another species of their same genus are
    # also known to be found there. Otherwise, keep their occurrences and rows as-is:
      
    for (Row in which(endsWith(Fossil_herbivores$Binomial_Original, 'sp'))) { # For each such species...
      for (Col in which(colnames(Fossil_herbivores) == 'Amboseli_Pleis'):ncol(Fossil_herbivores)) { # For each site...
        if (Fossil_herbivores[Row,Col] == 1 & # If the species occurs at the site and...
          sum(Fossil_herbivores[grep(strsplit(Fossil_herbivores$Binomial[Row], ' ')[[1]][1],
            Fossil_herbivores$Binomial), Col]) > 1) { # ...other species of the same genus also occur at the site...
              Fossil_herbivores[Row,Col] <- 0 }}}; rm(Row, Col) # Remove the occurrence of that species from that site.
          
# Perform merges of "Phylacine"/"PanTHERIA" to both "Fossil_herbivores" and "Modern_herbivores" datasets, in order to
# integrate the body masses from "Phylacine"/"PanTHERIA" into them. Prioritize "Phylacine", since it is newer:

Datasets_merged_bodymass <- lapply(list(Fossil_herbivores, Modern_herbivores), function(Dataset) { # For each dataset...
  
  # Merge on the "Binomial" column of the datasets, first for "Phylacine" and then for "PanTHERIA":
  
  Dataset <- merge(x = Dataset, y = Phylacine[, c('Binomial', 'Mass_Kg_Phylacine')], by = 'Binomial', all.x = TRUE)
  Dataset <- merge(x = Dataset, y = PanTHERIA[, c('Binomial', 'Mass_Kg_PanTHERIA')], by = 'Binomial', all.x = TRUE)
  
  # Create a column of finalized body masses that contain the mass reported by "Phylacine"/"PanTHERIA". If the mass from
  # "Phylacine" is missing, use the "PanTHERIA" mass, and if that is missing, keep the original mass value:
  
  Dataset <- add_column(Dataset, Mass_Kg_Final = ifelse(!is.na(Dataset$Mass_Kg_Phylacine), Dataset$Mass_Kg_Phylacine,
    ifelse(!is.na(Dataset$Mass_Kg_PanTHERIA), Dataset$Mass_Kg_PanTHERIA, Dataset$Mass_Kg_Herbivores)),
    .before = 'Size.Class')
  
  # Move the "Binomial" column immediately before the finalized body masses column:
  
  Dataset <- add_column(Dataset, Binomial_Final = Dataset$Binomial, .before = 'Mass_Kg_Final')

  # Remove the "Mass_Kg_PanTHERIA", "Mass_Kg_Phylacine", "Mass_Kg_Herbivores", and "Binomial" columns, and return result:
  
  Dataset <- subset(Dataset, select = -c(Mass_Kg_PanTHERIA, Mass_Kg_Phylacine, Mass_Kg_Herbivores, Binomial))
  return(Dataset) })
    
Fossil_herbivores <- Datasets_merged_bodymass[[1]] # Assign "Fossil_herbivores" to equal the first modified dataset.
Modern_herbivores <- Datasets_merged_bodymass[[2]] # Assign "Modern_herbivores" to equal the second modified dataset.

# For rows in "Fossil_herbivores" with missing body mass data, obtain values from the literature where possible, or from
# the values of other species of the same genus:

  # Create a marker column that indicates which rows have these approximations, and which do not:

  Fossil_herbivores <- add_column(Fossil_herbivores, Mass_approximated = ifelse(is.na(Fossil_herbivores$Mass_Kg_Final),
    TRUE, FALSE), .before = 'Size.Class')
    
  # Body masses derived from the literature:  
  
    # Alcelaphini sp. (NA if "Alcelaphini" genus removed from data, as "Alcelaphini" is a tribe):
    
    Fossil_herbivores$Mass_Kg_Final[startsWith(Fossil_herbivores$Binomial_Final, 'Alcelaphini')] <- mean(c(
      Fossil_herbivores$Mass_Kg_Final[startsWith(Fossil_herbivores$Binomial_Final, 'Alcelaphus')], 158))
    
    # Anancus sp. (including that Hautier et al. claims that kenyensis = osiris):
    
    Fossil_herbivores$Mass_Kg_Final[startsWith(Fossil_herbivores$Binomial_Final, 'Anancus') & is.na(
      Fossil_herbivores$Mass_Kg_Final)] <- c(Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 
      'Anancus kenyensis'], 6803.89, 9071.85)
  
    # Ancylotherium cheboitense:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Ancylotherium cheboitense'] <- 3000
  
    # Antidorcas sp.:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Antidorcas sp'] <- mean(c(42.1, 36.8, 45.5, 
      43.6, 46.8, 45.0))
  
    # Archaeopotamus harvardi and lothagamensis:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Archaeopotamus harvardi'] <- 
      Modern_herbivores$Mass_Kg_Final[Modern_herbivores$Binomial_Final == 'Hippopotamus amphibius']
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Archaeopotamus lothagamensis'] <- 1000
  
    # Bos sp.:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Bos sp'] <- mean(c(350, 800))
  
    # Cainochoerus africanus:
    
    Fossil_herbivores$Mass_Kg_Final[startsWith(Fossil_herbivores$Binomial_Final, 'Cainochoerus')] <- c(5, 7.5, 10)

    # Camelus sp.:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Camelus sp'] <- mean(c(300, 320))
  
    # Cephalophus sp.:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Cephalophus sp'] <- 36
  
    # Ceratotherium efficax:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Ceratotherium efficax'] <- 
      Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Ceratotherium simum']

    # Chemositia tugenensis:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Chemositia tugenensis'] <- mean(c(
      Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Ancylotherium cheboitense'], 
      Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Ancylotherium hennigi']))
    
    # Diceros sp.:
    
    Fossil_herbivores$Mass_Kg_Final[startsWith(Fossil_herbivores$Binomial_Final, 'Diceros') & is.na(
      Fossil_herbivores$Mass_Kg_Final)] <- c(700, 1400, mean(c(700, 1400)))

    # Dytikodorcas sp.:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Dytikodorcas sp'] <- mean(
      Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Aepyceros premelampus'])
    
    # Equus capensis and stenonis:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final %in% c('Equus capensis', 'Equus stenonis')] <- 500

    # Eurygnathohippus sp.:
    
    Fossil_herbivores$Mass_Kg_Final[startsWith(Fossil_herbivores$Binomial_Final, 'Eurygnathohippus') & is.na(
      Fossil_herbivores$Mass_Kg_Final)] <- Modern_herbivores$Mass_Kg_Final[Modern_herbivores$Binomial_Final == 
      'Equus zebra']
    
    # Hexaprotodon bruneti:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Hexaprotodon bruneti'] <- 1000

    # Hippopotamus and Hippotragus sp.:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final %in% c('Hippopotamus coryndonae',
      'Hippopotamus protoamphibius')] <- 1200
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Hippopotamus dulu'] <- 1300
    
    for (Pair in list(c('Hippopotamus afarensis', 'Hippopotamus amphibius'), c('Hippopotamus kaisensis',
      'Hippopotamus gorgops'), c('Hippopotamus karumensis', 'Hippopotamus protoamphibius'), c('Hippotragus cookei',
      'Hippotragus gigas'))) { # For each of these pairs of species...
        Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == Pair[1]] <- # Set the first species' mass...
          Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == Pair[2]] }; rm(Pair) # ...equal to second.
    
    # Kobus and Kolpochoerus sp.:	
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final %in% c(
      'Kobus ancystrocera', 'Kobus subdolus', 'Kolpochoerus afarensis', 'Kolpochoerus deheinzelini')] <- 80
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Kolpochoerus olduvaiensis'] <- 1000
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Kolpochoerus phillipi'] <- 
      Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Kolpochoerus majus']

    # Loxodonta sp.:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final %in% c('Loxodonta cookei', 'Loxodonta exoptata')] <-
      Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Mammuthus subplanifrons']

    # Makapania sp.:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Makapania sp'] <- 263
    
    # Neotragus sp.:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Neotragus sp'] <- mean(
      Modern_herbivores$Mass_Kg_Final[startsWith(as.character(Modern_herbivores$Binomial_Final), 'Neotragus')])

    # Notochoerus capensis and clarkei:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Notochoerus capensis'] <- mean(
      Fossil_herbivores$Mass_Kg_Final[startsWith(as.character(Fossil_herbivores$Binomial_Final), 'Notochoerus')],
      na.rm = TRUE)
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Notochoerus clarkei'] <- 250
    
    # Nyanzachoerus sp.:
    
    for (Sp in which(grepl('Nyanzachoerus', Fossil_herbivores$Binomial_Final) & is.na(Fossil_herbivores$Mass_Kg_Final))) {
      Fossil_herbivores$Mass_Kg_Final[Sp] <- mean(Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Family == 
      'Suidae' & Fossil_herbivores$Size.Class == Fossil_herbivores$Size.Class[Sp]], na.rm = TRUE) }; rm(Sp)
    
    # Paleotragus sp.:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Paleotragus sp'] <- 1000
    
    # Parantidorcas latifrons:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Parantidorcas latifrons'] <- 30
    
    # Parmularius pachyceras:

    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Parmularius pachyceras'] <- 175
    
    # Phacochoerus antiquus:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Phacochoerus antiquus'] <- mean(
      Fossil_herbivores$Mass_Kg_Final[grepl('Phacochoerus|Potamochoerus', Fossil_herbivores$Binomial_Final)], na.rm = TRUE)
    
    # Primelephas korotorensis:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Primelephas korotorensis'] <-
      Fossil_herbivores$Mass_Kg_Final[grepl('Stegotetrabelodon', Fossil_herbivores$Binomial_Final)]
    
    # Rabaticeras arambourgi and lemutai:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final %in% c('Rabaticeras arambourgi',
      'Rabaticeras lemutai')] <- mean(Fossil_herbivores$Mass_Kg_Final[grepl('Damaliscus', 
      Fossil_herbivores$Binomial_Final)], na.rm = TRUE)
    
    # Sivatherium hendeyi:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Sivatherium hendeyi'] <- 1200

    # Ugandax demissum and gautieri:
    
    Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final %in% c('Ugandax demissum', 'Ugandax gautieri')] <-
      Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Syncerus caffer']
      
  # Body masses derived from known means:
  
    # Size class-dependent:
  
    for (Class in c('s2', 's3')) { # For taxa with missing values in these two size classes, and with knowns in both...
      Fossil_herbivores$Mass_Kg_Final[startsWith(Fossil_herbivores$Binomial_Final, 'Aepyceros') & 
        Fossil_herbivores$Size.Class == Class & is.na(Fossil_herbivores$Mass_Kg_Final)] <- mean(
        Fossil_herbivores$Mass_Kg_Final[startsWith(Fossil_herbivores$Binomial_Final, 'Aepyceros') & 
        Fossil_herbivores$Size.Class == Class], na.rm = TRUE) # Only take the mean of species from the same size class.
      if (Class == 's3') { # For taxa with missing vals only in s3, but with knowns in both s3 and other classes...
        for (Sp in c('Equus', 'Kobus', 'Kolpochoerus', 'Metridiochoerus', 'Nyanzachoerus', 'Parmularius', 'Tragelaphus')) {
            Fossil_herbivores$Mass_Kg_Final[startsWith(Fossil_herbivores$Binomial_Final, Sp) & 
              Fossil_herbivores$Size.Class == Class & is.na(Fossil_herbivores$Mass_Kg_Final)] <- mean(
              Fossil_herbivores$Mass_Kg_Final[startsWith(Fossil_herbivores$Binomial_Final, Sp) & 
              Fossil_herbivores$Size.Class == Class], na.rm = TRUE) }}}; rm(Class, Sp) # Only take mean of s3 species.
    
    # Not size class-dependent:
    
      # Genus-level or higher:
    
      Fossil_herbivores$Mass_Kg_Final[startsWith(Fossil_herbivores$Binomial_Final, 'Beatragus') & is.na(
        Fossil_herbivores$Mass_Kg_Final)] <- mean(Fossil_herbivores$Mass_Kg_Final[startsWith(
        Fossil_herbivores$Binomial_Final, 'Beatragus') & !is.na(Fossil_herbivores$Mass_Kg_Final)]) # Mean of same genus.
      Fossil_herbivores$Mass_Kg_Final[startsWith(Fossil_herbivores$Binomial_Final, 'Damalacra')] <- 
        mean(Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Tribe == 'Alcelaphini'], na.rm = TRUE) # Of same tribe.
    
      # Species-level (for taxa with missing values and only one size class):
      
      for (Gen in c('Ceratotherium', 'Connochaetes', 'Damalborea', 'Damaliscus', 'Gazella', 'Giraffa', 
        'Hippopotamus', 'Hippotragus', 'Loxodonta', 'Madoqua', 'Megalotragus', 'Menelikia', 'Notochoerus', 'Oryx',
        'Pelorovis', 'Phacochoerus', 'Potamochoerus', 'Praedamalis', 'Rabaticeras', 'Raphicerus', 'Redunca',
        'Simatherium', 'Sivatherium', 'Syncerus', 'Taurotragus', 'Tragoportax', 'Ugandax')) { # For each genus...
          Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == paste(Gen, 'sp', sep = ' ')] <- mean(
            Fossil_herbivores$Mass_Kg_Final[startsWith(Fossil_herbivores$Binomial_Final, Gen)], na.rm = TRUE) }; rm(Gen)
      
  # Body masses corrected from above:
      
  Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Gazella granti'] <- 55.46446
  Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores$Binomial_Final == 'Tragelaphus spekei'] <- 77.99923

###########################################################################################################################
###########################################################################################################################

# PART II: INTEGRATING DENTAL DATA INTO FAITH ET AL. DATASETS, AND FURTHER DATA CLEANING

# Read in the dental data for both fossil and modern herbivores:

Fossil_herbivores_dental <- read.csv('Fossil_Herbivores_Dental.csv')
    
Modern_herbivores_dental <- read.csv('Modern_Herbivores_Dental.csv', nrows = 650)
colnames(Modern_herbivores_dental)[1] <- 'TAXON'

# Perform a merge of the "Modern_herbivores_dental" and "Modern_herbivores" datasets, in order to integrate dental traits
# into the latter. This can be done simply, since there is relatively little taxonomic mismatch here:
  
  # Merge on the respective "TAXON" and "Binomial_Final" columns, which correspond to each other:

  Modern_herbivores <- merge(x = Modern_herbivores, y = Modern_herbivores_dental[, c('TAXON', 'HYP', 'LOP', 'AL')],
    by.x = 'Binomial_Final', by.y = 'TAXON', all.x = TRUE)
  
  # Re-order and re-name specific columns in the resultant "Modern_herbivores" dataset accordingly:
  
  Modern_herbivores <- Modern_herbivores[, c('Order', 'Family', 'Binomial_Final', colnames(Modern_herbivores)[
    which(colnames(Modern_herbivores) == 'Mass_Kg_Final'):which(colnames(Modern_herbivores) == 'Diet')],
    colnames(Modern_herbivores)[which(colnames(Modern_herbivores) == 'HYP'):which(colnames(Modern_herbivores) == 'AL')],
    colnames(Modern_herbivores)[which(colnames(Modern_herbivores) == 'ABER'):which(colnames(Modern_herbivores) == 'ZIN')])]
                                             
  colnames(Modern_herbivores)[colnames(Modern_herbivores) == 'LOP.y'] <- 'LOP'
  
  # Examine and fill in the resultant entries with missing dental data on a case-by-case basis, in case of taxonomic
  # mismatches or other factors that prevented them from being filled in properly from the merge. Also modify the
  # "Mass_Kg_Final" column of these entries (and one more entry with a mismatch), drawing their mass from "PanTHERIA":
  
  Modern_fillin_cols <- c('HYP', 'LOP', 'AL', 'Mass_Kg_Final')
  Modern_fillin_taxa <- c('Cephalophus monticola', 'Damaliscus hunteri', 'Gazella dama', 'Gazella granti', 
    'Gazella rufifrons', 'Gazella soemmerringii', 'Gazella thomsonii', 'Hippopotamus amphibius')
  Modern_fillin_vals <- list(c(1, 2, 0, 4.89605), c(3, 2, 0, 79.13217), c(3, 2, 0, 71.42481), c(3, 2, 0, 55.46446),
    c(3, 2, 1, 26.99977), c(3, 2, 0, 41.58293), c(3, 2, 0, 22.90743), c(2, 1, 0, 1417.49))
  
  for (Row in 1:length(Modern_fillin_taxa)) { # For each taxon which requires data filling...
    Modern_herbivores[Modern_herbivores$Binomial_Final == Modern_fillin_taxa[Row], Modern_fillin_cols] <-
      Modern_fillin_vals[[Row]] }; rm(Row) # Perform the filling in.

  Modern_herbivores[Modern_herbivores$Binomial_Final == 'Cephalophus maxwellii', 'Mass_Kg_Final'] <- 8.55769

# Resolve taxonomic differences in the "TAXON" columns of both the datasets "Fossil_herbivores_dental" and 
# "Modern_herbivores_dental", and the "Binomial_Final" column of "Fossil_herbivores". Then, perform a merge for
# "Fossil_herbivores" to integrate dental traits into it:

  # Convert the "TAXON" columns from factor to character:
  
  Fossil_herbivores_dental$TAXON <- as.character(Fossil_herbivores_dental$TAXON)
  Modern_herbivores_dental$TAXON <- as.character(Modern_herbivores_dental$TAXON)
  
  # For each row in "Fossil_herbivores", determine the row index in either of "Modern_herbivores_dental" or
  # "Fossil_herbivores_dental" for which there is the highest degree of string match, using the Levenshtein edit
  # distance and Levenshtein similarity formula to quantify this match. Classify a match as one with at least 60% match:
  
  String_match_modern <- sapply(Fossil_herbivores$Binomial_Final, function(Row1) { # Fossil compared to modern...
    sapply(Modern_herbivores_dental$TAXON, function(Row2) { 1 - (adist(Row2, Row1) / max(str_length(c(Row1, Row2)))) }) })
  String_match_fossil <- sapply(Fossil_herbivores$Binomial_Final, function(Row1) { # Fossil compared to fossil...
    sapply(Fossil_herbivores_dental$TAXON, function(Row2) { 1 - (adist(Row2, Row1) / max(str_length(c(Row1, Row2)))) }) })
  
  String_match <- rbind(String_match_modern, String_match_fossil) # Combine above results.
  rm(String_match_modern, String_match_fossil)
  
  String_match_indices <- apply(String_match, 2, function(Column) { ifelse(max(Column) >= 0.6, which.max(Column), NA) })
  
  # Create empty columns in "Fossil_herbivores" for the dental traits to be integrated:
  
  Fossil_herbivores <- add_column(Fossil_herbivores, HYP = NA, LOP = NA, AL = NA, .after = 'Diet')
  
  # For each row in "Fossil_herbivores", fill in dental traits from the dental datasets based on string match:
  
  for (Row in seq(nrow(Fossil_herbivores))) { # For each row...
    
    # Set all rows without matching dental data to have missing dental data values, or else...:
    
    if (is.na(String_match_indices[[Row]])) { Fossil_herbivores[Row, c('HYP', 'LOP', 'AL')] <- NA } else {
      
      # Keep in mind that string indices in the "String_match_indices" variable will be offset from those indices in the
      # "Fossil_herbivores_dental" dataset by a number equal to the length of "Modern_herbivores_dental", since both
      # dental datasets were concatenated row-wise before calculating "String_match_indices" above. Account for this
      # difference here, and assign dental data accordingly. Or else...:
      
      if (String_match_indices[[Row]] > nrow(Modern_herbivores_dental)) {
        Fossil_herbivores[Row, c('HYP', 'LOP', 'AL')] <- Fossil_herbivores_dental[String_match_indices[[Row]] - 
          nrow(Modern_herbivores_dental), c('HYP', 'LOP', 'AL')] } else {
      
        # Otherwise, the string indices in "String_match_indices" referring to data in "Modern_herbivores_dental" is not
        # offset, since this dental dataset was put on top of the other in the concatenation process. Account for this
        # here, and assign dental data accordingly:
      
        Fossil_herbivores[Row, c('HYP', 'LOP', 'AL')] <- Modern_herbivores_dental[String_match_indices[[Row]], 
          c('HYP', 'LOP', 'AL')] }}}
  
  # Resolve case-by-case instances in which dental trait data was entered incorrectly or was missing when it should have
  # in fact not been missing:
  
    # Create lists, much like those created for the "Modern_herbivores" dataset above, with information on specific
    # entries in "Fossil_herbivores" to be addressed:
  
    Fossil_fillin_cols <- c('HYP', 'LOP', 'AL')
    Fossil_fillin_genera <- c('Aepyceros', 'Anancus', 'Antidorcas', 'Beatragus', 'Bos', 'Budorcas', 'Cainochoerus',
      'Camelus', 'Diceros', 'Elephas', 'Equus', 'Gazella', 'Giraffa', 'Hippopotamus', 'Kobus', 'Loxodonta', 'Madoqua',
      'Mammuthus', 'Notochoerus', 'Parmularius', 'Phacochoerus', 'Redunca', 'Syncerus')
    Fossil_fillin_generavals <- list(c(3, 2, 0), c(1, 0, 0), c(3, 2, 0), c(3, 2, 0), c(3, 2, 0), c(3, 2, 0), c(3, 0, 0),
      c(3, 2, 0), c(2, 1, 1), c(3, 0, 0), c(3, 2, 0), c(3, 2, 0), c(1, 2, 1), c(2, 1, 0), c(3, 2, 0), c(3, 0, 0),
      c(2, 2, 0), c(3, 0, 0), c(NA, NA, NA), c(NA, NA, NA), c(3, 0, 0), c(3, 2, 0), c(3, 2, 0))
    Fossil_fillin_taxa <- c('Kobus sp', 'Antilope subtorta', 'Kolpochoerus majus', 'Prostrepsiceros vinayaki',
      'Tragelaphus gaudryi', 'Tragelaphus kyaloi', 'Tragelaphus lockwoodi', 'Tragelaphus moroitu', 'Tragelaphus nakuae',
      'Tragelaphus oryx', 'Tragelaphus pricei', 'Tragelaphus rastafari', 'Tragelaphus saraitu')
    Fossil_fillin_taxavals <- list(c(3, 2, 0), c(3, 2, 0), c(NA, NA, NA), c(NA, NA, NA), c(NA, NA, NA), c(NA, NA, NA),
      c(NA, NA, NA), c(NA, NA, NA), c(NA, NA, NA), c(NA, NA, NA), c(NA, NA, NA), c(NA, NA, NA), c(NA, NA, NA))
    
    # Handle/fill in all relevant "Binomial_Final" dental entries involving a common genus:
    
    for (Row in 1:length(Fossil_fillin_genera)) { # For each genus that requires data filling...
      Fossil_herbivores[startsWith(Fossil_herbivores$Binomial_Final, Fossil_fillin_genera[Row]) == TRUE, 
        Fossil_fillin_cols] <- aperm(array(Fossil_fillin_generavals[[Row]], dim = c(length(Fossil_fillin_cols), 
        nrow(Fossil_herbivores[startsWith(Fossil_herbivores$Binomial_Final, Fossil_fillin_genera[Row]) == TRUE, 
        Fossil_fillin_cols])))) }; rm(Row)
    
    # Handle/fill in all relevant "Binomial_Final" dental entries involving a specific species:
    
    for (Row in 1:length(Fossil_fillin_taxa)) { # For each species that requires data filling...
      Fossil_herbivores[Fossil_herbivores$Binomial_Final == Fossil_fillin_taxa[Row], Fossil_fillin_cols] <-
        aperm(array(Fossil_fillin_taxavals[[Row]], dim = c(length(Fossil_fillin_cols), nrow(Fossil_herbivores[
        Fossil_herbivores$Binomial_Final == Fossil_fillin_taxa[Row], Fossil_fillin_cols])))) }; rm(Row)

# For rows in "Fossil_herbivores" with missing dental data, obtain values from the literature where possible, or from the
# values of other species of the same genus:

  # Create a marker column that indicates which rows have these approximations, and which do not:
    
  Fossil_herbivores <- add_column(Fossil_herbivores, Dental_approximated = ifelse(is.na(Fossil_herbivores$HYP),
    TRUE, FALSE), .after = 'AL')

  # Insert the values obtained from the literature:
  
  Fossil_herbivores[grep('Ancylotherium|Chemositia', Fossil_herbivores$Binomial_Final), c('HYP', 'LOP')] <- list(1, 1)
  Fossil_herbivores[grep('Archaeopotamus|Deinotherium|Nyanzachoerus|Stegotetrabelodon',
    Fossil_herbivores$Binomial_Final), c('HYP', 'LOP')] <- list(1, 0)
  Fossil_herbivores[grep(paste('Alcelaphini|Damalacra|Damalborea|Eurygnathohippus|Hippotherium|Makapania|Megalotragus|',
    'Menelikia|Nitidarcus|Notochoerus|Numidocapra|Parantidorcas|Parmularius|Pelorovis|Rabaticeras|',
    'Rusingoryx|Simatherium|Thaleroceros|Ugandax|Zephyreduncinus', sep = ''), 
    Fossil_herbivores$Binomial_Final), c('HYP', 'LOP')] <- list(3, 2)
  Fossil_herbivores[grep('Bouria|Brabovus|Dytikodorcas|Praedamalis|Prostrepsiceros|Sivatherium|Tragelaphus|Tragoportax',
    Fossil_herbivores$Binomial_Final), c('HYP', 'LOP')] <- list(2, 2)
  Fossil_herbivores[grep('Paleotragus|Palaeotragus|Tragelaphus buxtoni|Tragoportax abyssinicus',
    Fossil_herbivores$Binomial_Final), c('HYP', 'LOP')] <- list(1, 2)
  Fossil_herbivores[grep('Tragelaphus oryx', Fossil_herbivores$Binomial_Final), c('HYP', 'LOP')] <- list(3, 2)
  Fossil_herbivores[grep('Kolpochoerus|Notochoerus|Primelephas|Stegodon', Fossil_herbivores$Binomial_Final),
    c('HYP', 'LOP')] <- list(2, 0)
  Fossil_herbivores[grep('Metridiochoerus|Notochoerus scotti', Fossil_herbivores$Binomial_Final),
    c('HYP', 'LOP')] <- list(3, 0)
  Fossil_herbivores[grep('Brachypotherium|Ugandax demissum|Gazella sp', Fossil_herbivores$Binomial_Final),
    c('HYP', 'LOP')] <- list(2, 1)

###########################################################################################################################
###########################################################################################################################
  
# PART III: INTEGRATING TEMPORAL DATA INTO FAITH ET AL. DATASETS, AND FINAL DATA CLEANING
  
# From this point forward, certain sites may be removed from analysis based on the degree to which we can rely on their
# measurements of woody/grassy cover (see below for details). Indicate here whether or not those sites should be removed:
  
Fossil_sites_choose <- 'Paleosol' # "All" or "Paleosol" - include all sites or only those with paleosol cover data?
Fossil_sites_term <- ifelse(Fossil_sites_choose == 'Paleosol', 'Soil', 'Soil|Literature') # Relevant key term for below.
  
# Read in the temporal data for fossil herbivore sites, and subset it to remove two blank columns:
  
Fossil_herbivores_temporal <- read.csv('Fossil_Herbivores_SiteAges.csv')
Fossil_herbivores_temporal <- Fossil_herbivores_temporal[, -c(9, 10)]

# Convert the "Min..Age" and "Max..Age" columns in the dataset from "factor" to "numeric", and set dataset entries with a
# "-" to being a missing value:

Fossil_herbivores_temporal[Fossil_herbivores_temporal == '-'] <- NA

Fossil_herbivores_temporal$Min..Age <- as.numeric(as.character(Fossil_herbivores_temporal$Min..Age))
Fossil_herbivores_temporal$Max..Age <- as.numeric(as.character(Fossil_herbivores_temporal$Max..Age))

# For each row/species in the "Fossil_herbivores" dataset, return the minimum, maximum, and mean ages during which each
# species lived, and save the results as columns in the dataset:

Fossil_herbivores <- add_column(Fossil_herbivores, Min_Age = NA, Max_Age = NA, Mean_Age = NA,
  .after = 'Dental_approximated') # Create the baseline columns.

Occurrence_sites <- list() # List to hold the names of sites that each species occurred in.
Occurrence_sites_indices <- list() # List to hold the row indices of those sites in "Fossil_herbivores_temporal".
Occurrence_sites_remove_orig <- read.csv('Fossil_Herbivores_SiteWoodyCover_V2.csv') # Df with potential sites to remove.
Occurrence_sites_remove <- as.character(Occurrence_sites_remove_orig$Site)[grepl('Literature',
  Occurrence_sites_remove_orig$Derived_From)] # Names of potential sites to remove.
Occurrence_sites_remove_indices <- match(Occurrence_sites_remove, Fossil_herbivores_temporal$Site) # Indices of sites.
Occurrence_sites_index_calc <- which(colnames(Fossil_herbivores) == 'Amboseli_Pleis') - 1 # Column index.

for (Row in seq(nrow(Fossil_herbivores))) { # For each species...
  
  # Return the names of the sites that the species occurred in:
  
  Occurrence_sites[[Row]] <- colnames(Fossil_herbivores[(Occurrence_sites_index_calc + 1):ncol(Fossil_herbivores)])[
    Fossil_herbivores[Row, (Occurrence_sites_index_calc + 1):ncol(Fossil_herbivores)] == 1]
  
  # Determine the row indices that those site names pertain to in the "Fossil_herbivores_temporal" dataset:
  
  Occurrence_sites_indices[[Row]] <- which(colnames(Fossil_herbivores) %in% Occurrence_sites[[Row]]) - 
    Occurrence_sites_index_calc
  
  # Potentially remove certain indices from "Occurrence_sites_indices" based on "Occurrence_sites_remove_indices" above:
  
  if (Fossil_sites_choose == 'Paleosol') { # If only including sites with paleosol data...
    Occurrence_sites_indices[[Row]] <- # Remove on next line.
      Occurrence_sites_indices[[Row]][!Occurrence_sites_indices[[Row]] %in% Occurrence_sites_remove_indices] }
  
  # Determine the minimum, maximum, and mean ages in "Fossil_herbivores_temporal" pertaining to those row indices, and
  # enter those values into the new columns in "Fossil_herbivores":
  
  suppressWarnings({ # R throws a long list of warnings here about setting calculations that only have "NA" in them to
    # "Inf". Because "Inf" will be replaced by "NA", the warnings are inconsequential.
  
  Fossil_herbivores$Min_Age[[Row]] <- min(Fossil_herbivores_temporal$Min..Age[Occurrence_sites_indices[[Row]]], na.rm = T)
  Fossil_herbivores$Max_Age[[Row]] <- max(Fossil_herbivores_temporal$Max..Age[Occurrence_sites_indices[[Row]]], na.rm = T)
  Fossil_herbivores$Mean_Age[[Row]] <- mean(Fossil_herbivores_temporal$Mean.Age[Occurrence_sites_indices[[Row]]],
    na.rm = T) }) }; rm(Row)

# In the three columns created above, replace all "Inf" and "NaN" with "NA", as they are missing values:

Fossil_herbivores$Min_Age[Fossil_herbivores$Min_Age == Inf] <- NA
Fossil_herbivores$Max_Age[Fossil_herbivores$Max_Age == -Inf] <- NA
Fossil_herbivores$Mean_Age[is.nan(Fossil_herbivores$Mean_Age)] <- NA

# Subset the "Fossil_herbivores" dataset to only include those rows with "Min_Age" and "Max_age" information:

Fossil_herbivores <- Fossil_herbivores[!is.na(Fossil_herbivores$Min_Age) & !is.na(Fossil_herbivores$Max_Age), ]

# Subset the "Modern_herbivores" dataset to only include rows of species that occur in the select sites with available
# grassy/woody cover data:

Occurrence_sites_include <- read.csv('Modern_Herbivores_SiteWoodyCover.csv') # Dataframe with the sites to include.
Occurrence_sites_include_indices <- match(Occurrence_sites_include$Code, colnames(Modern_herbivores)) # Indices of sites.
Modern_herbivores <- Modern_herbivores[sapply(1:nrow(Modern_herbivores), function(Row) { # For each row/species...
  sum(Modern_herbivores[Row, Occurrence_sites_include_indices], na.rm = TRUE) > 0 }), ] # Subset species.

# Visualize the distributions of body masses (the only continuous trait) to determine if any log transformations are
# necessary with it in "Fossil_herbivores" and "Modern_herbivores":

par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1,2)) # Set up plot layout settings.

lapply(list(Fossil_herbivores, Modern_herbivores), function(Dataset) { # For each dataset...
  hist(Dataset$Mass_Kg_Final, col = 'darkolivegreen3', xlab = '', ylab = '', cex.axis = 1.5, cex.main = 2.5, 
    main = ifelse(nrow(Dataset) > nrow(Modern_herbivores), 'Fossil', 'Modern')) # Histogram of masses in the dataset.
  mtext(side = 1, line = 3, 'Body Mass (Kg)', font = 2, cex = 1.5) # Format x-axis label.
  mtext(side = 2, line = 3, 'Frequency', font = 2, cex = 1.5) } ) # Format y-axis label.

# Because body mass distributions are skewed and not normal, log transform them in both datasets:

Fossil_herbivores$Mass_Kg_Final <- log(Fossil_herbivores$Mass_Kg_Final)
Modern_herbivores$Mass_Kg_Final <- log(Modern_herbivores$Mass_Kg_Final)

# Subset the datasets to only include species >= X kg in mass, based upon which "megafauna" will be analyzed (comment this
# out if all species are to be analyzed):

Megafauna_threshold <- 44

Fossil_herbivores <- Fossil_herbivores[which(Fossil_herbivores$Mass_Kg_Final >= log(Megafauna_threshold)), ]
Modern_herbivores <- Modern_herbivores[which(Modern_herbivores$Mass_Kg_Final >= log(Megafauna_threshold)), ]

# Determine the correlation between the age ranges and number of taxa at sites in "Fossil_herbivores":

Subset <- which(!is.na(Fossil_herbivores_temporal$Min..Age) & grepl(Fossil_sites_term, # Which sites relevant?
  Occurrence_sites_remove_orig$Derived_From[match(Fossil_herbivores_temporal$Site, Occurrence_sites_remove_orig$Site)]))
Fossil_herbivores_sitesonly <- Fossil_herbivores[, colnames(Fossil_herbivores)[which(colnames(Fossil_herbivores) == 
  'Amboseli_Pleis'):ncol(Fossil_herbivores)]] # Site occurrence data from "Fossil_herbivores".
Fossil_herbivores_cor <- cor(sapply(Fossil_herbivores_sitesonly[, Subset], sum),
  Fossil_herbivores_temporal$Max..Age[Subset] - Fossil_herbivores_temporal$Min..Age[Subset], use = 'complete.obs')
rm(Subset, Fossil_herbivores_sitesonly) # Remove variables created above that are no longer needed.

###########################################################################################################################
###########################################################################################################################

# PART IV: OBTAINING, VISUALIZING, AND ANALYZING WOODY/GRASSY COVER DATA

# Compile woody and/or grassy cover data for fossil herbivore sites, either by identifying it directly from literature
# sources and/or deriving it from paleosol soil carbonate records. Much of this data entry was done separately/manually.

# library(readxl) # For using the "read_xlsx" function.
# Fossil_herbivores_paleosol <- read_xlsx('Fossil_Herbivores_Paleosol.xlsx', skip = 3) # Paleosol data from Levin, 2016.
# Fossil_herbivores_paleosol <- Fossil_herbivores_paleosol[Fossil_herbivores_paleosol$`Data type` %in% c(
#   'soil carbonate - d13C, d18O', 'soil carbonate - d13C'), ] # Subset data to relevant rows.

Paleosol_collector <- function(Df_subset, Min_age, Max_age) { # Collects paleosol measurements for a given site.
  Values <- Df_subset$`δ13CVPDB (‰)`[Df_subset$`Age (Ma)` <= Max_age & Df_subset$`Age (Ma)` >= Min_age]
  return(Values[is.finite(Values)]) } # Return non-NA paleosol values.
Paleosol_to_woodycover <- function(Carbonate_measurements) { # Converts soil carbonate to fraction woody cover.
  return(sin(-1.06688 - 0.08538 * (Carbonate_measurements - 14))^2) } # Subtract values by 14 to account for Suess effect.
Paleosol_to_grassycover <- function(Carbonate_measurements) { # Converts soil carbonate to fraction grassy cover.
  return((((Carbonate_measurements - -12) * (1 - 0)) / (2 - -12)) + 0) } # From carbonate range to fraction grass range.

# Read in and process the final compiled dataframe(s):

Cover_type <- 'Woody' # "Woody" or "Grassy", depending on the cover type to be analyzed moving forward.

Fossil_herbivores_cover_sitenames <- colnames(Fossil_herbivores)[which(grepl('Amboseli',
  colnames(Fossil_herbivores))):ncol(Fossil_herbivores)] # Site names that are more compatible with our cover data.
Fossil_herbivores_cover_sitenames <- Fossil_herbivores_cover_sitenames[!is.na(Fossil_herbivores_temporal$Min..Age)]

Fossil_herbivores_cover <- read.csv(paste0('Fossil_Herbivores_Site', Cover_type, 'Cover_V2.csv')) # Read in final df.
Fossil_herbivores_cover$Site <- Fossil_herbivores_cover_sitenames # Put the site names above in a new df column.
Fossil_herbivores_cover$Min_Age <- as.numeric(as.character(Fossil_herbivores_cover$Min_Age)) # Convert to numeric.
Fossil_herbivores_cover$Max_Age <- as.numeric(as.character(Fossil_herbivores_cover$Max_Age)) # Convert to numeric.

Cover_SD <- 'Mean' # "Mean", "Above", "Below" - should cover values be analyzed at their mean, one SD above, or one below?
if (Cover_SD == 'Above') { Cover_shift <- Fossil_herbivores_cover[, grepl('_SD', colnames(Fossil_herbivores_cover))] }
if (Cover_SD == 'Below') { Cover_shift <- -Fossil_herbivores_cover[, grepl('_SD', colnames(Fossil_herbivores_cover))] }
if (Cover_SD == 'Mean') { Cover_shift <- 0 } # No addition or subtraction of SD for this option.
Fossil_herbivores_cover[, grepl('_Mean', colnames(Fossil_herbivores_cover))] <- # Potentially modify the "Mean" column.
  Fossil_herbivores_cover[, grepl('_Mean', colnames(Fossil_herbivores_cover))] + Cover_shift # Apply shift in vals.

# Determine the correlation between sample size and standard deviation of vegetation cover in "Fossil_herbivores_cover":

Subset <- grep('Soil', Fossil_herbivores_cover$Derived_From) # Which sites have SD?
Fossil_herbivores_cover_cor <- cor(Fossil_herbivores_cover[Subset, grepl('_SD', colnames(Fossil_herbivores_cover))],
  Fossil_herbivores_cover[Subset, grepl('Sample_Size', colnames(Fossil_herbivores_cover))], use = 'complete.obs')
rm(Subset) # Remove the subset created above.

# Read in woody cover data for modern herbivore sites (at 2.5-minute resolution), taken directly from Barr et al., 2020.
# Site codes for that data were manually entered, such that sites could be matched to sites in the Faith et al. datasets:

if (Cover_type == 'Woody') { Modern_herbivores_cover <- read.csv('Modern_Herbivores_SiteWoodyCover.csv') }

###########################################################################################################################
###########################################################################################################################
 
# PART V: ANALYSIS OF HERBIVORE BODY MASS AND DENTAL TRAIT COMPOSITION THROUGH TIME (FUNCTIONAL DIVERSITY)

# Define a variable that will represent the number of sampling iterations to perform for all bootstrap methods done:

Reps <- 1000
set.seed(500) # Set a random seed to ensure reproducibility of randomized results below.  

# Output and analyze plots depicting, across all sites, the regional fraction of woody cover and the compositions of each
# separate herbivore trait of interest over time:

  # Set up the single lists that will be used and processed in the "Trait" "for" loop below:

  Fossil_herbivores_traits <- list()
  Modern_traitmetrics <- list()
  Modern_traitmetrics_sites_sd <- list()
  Modern_traitmetrics_samples_sd <- list()
  Overall_traitmetrics <- list()
  Overall_traitmetrics_sites_sd <- list()
  Overall_traitmetrics_samples_sd <- list()

  # Set up the trait and metric names that will be used and processed in the "Trait" "for" loop:
  
  Trait_names <- list('Body Mass', 'Hypsodonty', 'Loph Count')
  Traitmetric_names <- list('Mean', 'Median', 'SD')
  
  # Set up the nested lists that will be used and processed in the "Trait" "for" loop:
  
  Fossil_traitmetrics <- list()
  Fossil_traitmetrics_sites_sd <- list()
  Fossil_traitmetrics_samples_sd <- list()
  Fossil_traitsubset <- list()
  Fossil_traitsubset_size <- list()
  Plot_df <- list()
  Plot_df_breakpoints <- list()
  Plot_df_breakpoints_sig <- list()
  Trait_plot <- list()
  Y_ticks_traits <- list()
  
  for (Index in seq(length(Trait_names))) { Fossil_traitmetrics[[Index]] <- list(); 
    Fossil_traitmetrics_sites_sd[[Index]] <- list(); Fossil_traitmetrics_samples_sd[[Index]] <- list();
    Fossil_traitsubset[[Index]] <- list(); Fossil_traitsubset_size[[Index]] <- list(); Plot_df[[Index]] <- list();
    Plot_df_breakpoints[[Index]] <- list(); Plot_df_breakpoints_sig[[Index]] <- list(); Trait_plot[[Index]] <- list();
    Y_ticks_traits[[Index]] <- list() }; rm(Index)
  
  # Set up the bounds of a set of time bins spanning the past 7.5 million years. Organize the bins such that they each
  # cover a portion of time and are either sequential or overlapping, the latter of which creates a moving average:

  Fossil_timebins <- list() # List to hold the bounds of bins of time.
  Counter <- 1 # Counter for adding elements to "Fossil_timebins".
  
  Timebins_shift <- 0.25 # Set up the time-shift interval from the previous to next bin (0.1 OR 0.25 OR 0.5).
  Timebins_max <- ifelse(Timebins_shift == 0.1, 7.4, 7.5) # Set maximum age past which bins will no longer be established.

  Timebin_current_min <- 0 # Set the minimum age of first bin.
  Timebin_current_max <- Timebin_current_min + 0.25 # Set the maximum age of first bin (0.1 OR 0.25 OR 0.5 above minimum).
  
  while (Timebin_current_max <= Timebins_max) { # Ensure that time bins do not exceed 7.5 million years.
    Fossil_timebins[[Counter]] <- c(Timebin_current_min, Timebin_current_max) # Append the current bin.
    
    Counter <- Counter + 1 # Update the "Fossil_timebins" counter.
    Timebin_current_min <- Timebin_current_min + Timebins_shift # Update the minimum age of the next bin.
    Timebin_current_max <- Timebin_current_max + Timebins_shift } # Update the maximum age of the next bin.
  rm(Counter, Timebins_max, Timebin_current_min, Timebin_current_max)
  
  Fossil_timebins <- rev(Fossil_timebins) # Reverse the order of "Fossil_timebins" so it goes from oldest to newest.
  
  # Set up the information about climate and hominin events that will be used and processed in the analyses that follow.
  # This includes the creation of rectangles that will be placed in the plots below to signify when the events occurred:

  library(ggplot2) # For general plotting.
  library(plyr) # For using the "round_any" function to format plot axis labels.
  
  Events <- c(7.5, 7, 5, 3.3, 3.0, 3.0, 1.9, 1.9, 1.7, 1.7, 1, 1, 0.0) # Assign time points/periods to each event.
  Events_breaks <- sapply(Events, function(Time) { # For each time point, find the index of the time bin representing it...
    Tmp <- sapply(append(Fossil_timebins, 0), max) - Time; return(which.min(Tmp[Tmp>=0])) }) # Find index.
  Events_IDs <- c('Clim', 'Hom', 'Clim', 'Clim2', 'Clim2', 'Hom', 'Clim3', 'Hom', 'Clim3', 'Hom', 'Clim4', 'Hom', 'Clim4')
  
  Breaks_width <- 0.25 # Width of single-event rectangles.
  
  Breaks_clim_rects <- seq(Events_breaks[Events_IDs == 'Clim'][1], Events_breaks[Events_IDs == 'Clim'][2],
    length.out = 50) # Set up a sequence of rectangles that will give the grassland expansion rectangle a gradient fill.
  Breaks_clim_cols <- colorRampPalette(c('deepskyblue3', 'white'))(length(Breaks_clim_rects)-1) # Gradient colors.
  Breaks_clim <- list() # List to hold the sequence of rectangles and their respective colors.
  for (Rect in 1:(length(Breaks_clim_rects)-1)) { # For each rectangle in the sequence...
    Breaks_clim[[Rect]] <- annotate('rect', xmin = Breaks_clim_rects[Rect], xmax = Breaks_clim_rects[Rect+1], ymin = -Inf,
      ymax = Inf, alpha = 0.5, fill = Breaks_clim_cols[Rect]) }; rm(Rect) # Set it up with its associated color.
  
  for (Event in unique(Events_IDs[grepl('Clim.', Events_IDs)])) { # For each additional climate event...
    Tmp <- annotate('rect', xmin = Events_breaks[Events_IDs == Event][1], xmax = Events_breaks[Events_IDs == Event][2],
      ymin = -Inf, ymax = Inf, alpha = 0.3, fill = 'deepskyblue3') # Create a rectangle for it.
    assign(paste('Breaks_', tolower(Event), sep = ''), Tmp) }; rm(Event, Tmp) # Assign the rectangle to a variable.
  
  Breaks_hom <- sapply(Events_breaks[Events_IDs == 'Hom'], function(Pt) { # For each hominin event...
    annotate('rect', xmin = Pt - Breaks_width, xmax = Pt + Breaks_width, ymin = -Inf, ymax = Inf, alpha = 0.3,
    fill = 'orange2') }) # Hominin rectangle surrounding single time point.
  
  # Using "Fossil_timebins" and the information about events above, create a version of the plot of woody or grassy cover
  # over time that is compatible with the trait/biodiversity plots that will be created below:
  
    # Create a dataframe for the plot. In so doing, perform bootstrap resampling (described in more detail in other areas
    # of this script below) to account for temporal sampling bias:
  
    Bin_choice <- ifelse(Timebins_shift < 0.5, 2, 1) # Which is the oldest time bin with sample size >1?
    Sample_size_cover <- nrow(Fossil_herbivores_cover[!is.na(Fossil_herbivores_cover$Max_Age) & 
      Fossil_herbivores_cover$Max_Age >= Fossil_timebins[[Bin_choice]][1] & Fossil_herbivores_cover$Min_Age <=
      Fossil_timebins[[Bin_choice]][2], ]) # Sample size for bootstrap resampling (oldest time bin with sample size > 1).
    rm(Bin_choice) # No longer needed.
  
    Plot_df_cover_bins <- sapply(Fossil_timebins, function(Bin) { # For each time bin...
      Column <- grep(paste0(Cover_type, 'Cover_Mean'), colnames(Fossil_herbivores_cover)) # Which column has cover data?
      Samples <- Fossil_herbivores_cover[!is.na(Fossil_herbivores_cover$Max_Age) & Fossil_herbivores_cover$Max_Age >=
        Bin[1] & Fossil_herbivores_cover$Min_Age <= Bin[2] & grepl(Fossil_sites_term,
        Fossil_herbivores_cover$Derived_From), Column] # Cover values for sites in bin.
      Replications <- replicate(Reps, mean(sample(Samples, Sample_size_cover, replace = TRUE), na.rm = TRUE))
      return(list(mean(Replications, na.rm = TRUE), sd(Replications, na.rm = TRUE))) }) # Bootstrapped mean and SD for bin.
    
    Replications_modern <- replicate(Reps, mean(sample(Modern_herbivores_cover[, grep(Cover_type, colnames(
      Modern_herbivores_cover))], Sample_size_cover, replace = TRUE), na.rm = TRUE)) # Bootstrap replications for modern.
    Plot_df_cover_bins <- cbind(Plot_df_cover_bins, list(mean(Replications_modern, na.rm = TRUE), sd(
      Replications_modern, na.rm = TRUE))); rm(Replications_modern) # Add in modern cover data.
    Plot_df_cover_bins <- data.frame(Age = 1:(length(Fossil_timebins) + 1), FractionCover = unlist(Plot_df_cover_bins[1,]),
      FractionCover_SD = unlist(Plot_df_cover_bins[2,])) # Store all data compiled above in dataframe.
  
    # Record the primary breakpoint that is associated with that dataframe from breakpoint analysis:
    
    library(segmented) # For using the "segmented" function to perform breakpoint analysis.
    
    if (exists('FractionCover')) { rm(FractionCover) }; if (exists('Age')) { rm(Age) } # Remove these variables.
    Plot_df_breakpoints_cover <- segmented(lm(FractionCover ~ Age, data = Plot_df_cover_bins), seg.Z = ~Age,
      control = seg.control(n.boot = Reps)) # Model for breakpoint.
    Plot_df_breakpoints_cover <- round(Plot_df_breakpoints_cover$psi[, 2]) # Actual point.
    Plot_df_breakpoints_cover_sig <- davies.test(lm(FractionCover ~ Age, data = Plot_df_cover_bins))$p.value # Needed?
  
    # Create the plot. Do not actually plot it here, as it will be plotted later on with other plots:
    
    Plot_cover <- ggplot(Plot_df_cover_bins, aes(Age, FractionCover)) + # Select dataframe and axes for plot.
      geom_point(size = 3) + # Set the plot to be a scatter plot.
      geom_smooth(color = 'black', span = 0.75) + # Set the plot to have a curve with LOESS smoothing.
      geom_errorbar(aes(ymin = FractionCover - FractionCover_SD, ymax = FractionCover + FractionCover_SD), width = 0.5,
        size = 1) + # Bootstrap samples error.
      geom_vline(xintercept = ifelse(is.na(Plot_df_breakpoints_cover), min(Events_breaks),
        Plot_df_breakpoints_cover), linetype = 'solid', size = 3, alpha = 0.5) + # Plot the breakpoint.
      labs(y = paste(Cover_type, 'Cover\n(Fraction)')) + # Set only Y label.
      ylim(0, 0.75) + # max(Plot_df_cover_bins$FractionCover + Plot_df_cover_bins$FractionCover_SD, na.rm = T)) + # Y bound.
      scale_x_continuous(breaks = Events_breaks, labels = rep('', length(Events)), expand = c(0.0075,0)) + # X ticks.
      # scale_y_continuous(breaks = c(0.25, 0.5, 0.75), labels = function(Lab) { sprintf('%.2f', Lab) }) + # Y ticks.
      theme(panel.background = element_rect(fill = 'ivory2'), panel.grid = element_blank(), plot.background = element_rect(
        fill = 'white'), axis.line = element_line(color = 'black'), axis.title.x = element_blank(), axis.text.x = 
        element_text(size = 2, angle = 60, hjust = 1), axis.title.y = element_text(size = 20, margin = margin(0,15,0,0)),
        axis.text.y = element_text(size = 20, color = 'black'), plot.margin = unit(c(0.5, 0.5, 0, 0.5), 'lines')) # Format.
    
  # Define a function that will calculate each trait's mean, median, and standard deviation values over time, keeping in
  # mind that the metric values of no numbers is "NA", and that the SD of a single number is also "NA":
  
  Metrics <- function(Trait) { c(ifelse(length(Trait) == 0, NA, mean(Trait, na.rm = TRUE)), ifelse(length(Trait) == 0,
    NA, median(Trait, na.rm = TRUE)), ifelse(length(Trait) <= 1, NA, sd(Trait, na.rm = TRUE))) }
  
  # Define the sample size used to calculate metrics (using the "Metrics" function) from data subsets in the "Trait" "for"
  # loop, as well as metrics in the 3D plots below. Specifically set the sample size to the number of species occurring in
  # the oldest time period, as older time periods have smaller absolute sample sizes due to temporal sampling bias:
  
  Sample_size_func <- nrow(Fossil_herbivores[Fossil_herbivores$Max_Age >= Fossil_timebins[[1]][1] &
    Fossil_herbivores$Min_Age <= Fossil_timebins[[1]][2], ])
    
  # Load libraries that will be used in the "Trait" "for" loop and after, create settings for the loop, and initiate it:
  
  library(gridExtra) # For making and arranging plots on the same panel.
  library(grid) # For making and arranging plots on the same panel.
  library(cowplot) # For making and arranging plots on the same panel.
  Counter <- 1 # Counter for referring to each trait through the loop.
  Samplesize_choice <- 'Species' # "Sites" or "Species", depending upon the sample size of interest.

  for (Trait in list('Mass_Kg_Final', 'HYP', 'LOP')) { # For each trait...
    
    # Subset "Fossil_herbivores" to only include those rows with information about the trait of interest across all sites:
    
    Fossil_herbivores_traits[[Counter]] <- Fossil_herbivores[!is.na(Fossil_herbivores[Trait]), ]
    
    # Bin the rows in "Fossil_herbivores_traits" in the time bins specified in "Fossil_timebins" above, according to each
    # row's "Min_Age" and "Max_Age". For each bin, calculate the mean, median, and standard deviation of the trait of
    # interest using the "Metrics" function. Also calculate the standard deviations of each "Metrics" calculation, based
    # on how they vary from each site to site, as well as from bootstrap sample to sample (to account for sampling bias):
  
    for (Bin in 1:length(Fossil_timebins)) { # For each time bin...
      
      # Create a "Fossil_herbivores_traits" data subset containing all of the species occurring in the time bin:
      
      Fossil_traitsubset[[Counter]][[Bin]] <- Fossil_herbivores_traits[[Counter]][
        Fossil_herbivores_traits[[Counter]]$Max_Age >= Fossil_timebins[[Bin]][1] & 
        Fossil_herbivores_traits[[Counter]]$Min_Age <= Fossil_timebins[[Bin]][2], ]
      
      # Remove any rows from "Fossil_traitsubset" that are recognized only as "sp" and are not the only members of their
      # genus present in the dataframe, as such rows refer to non-distinct species:
      
      Rows <- which(endsWith(Fossil_traitsubset[[Counter]][[Bin]]$Binomial_Original, 'sp')) # Rows denoted only as "sp".
      if (length(Rows) > 1) { # If such rows exist...
        Rows_genera <- sapply(strsplit(Fossil_traitsubset[[Counter]][[Bin]]$Binomial_Final[Rows], ' '), '[[', 1) # Gen.
        Rows_genera_all <- sapply(strsplit(Fossil_traitsubset[[Counter]][[Bin]]$Binomial_Final, ' '), '[[', 1) # All gen.
        Rows <- Rows[unname(which(table(Rows_genera_all)[match(Rows_genera, names(table(Rows_genera_all)))] > 1))] # Sub.
        if (length(Rows) > 1) { # If such rows still exist at this phase...
          Fossil_traitsubset[[Counter]][[Bin]] <- Fossil_traitsubset[[Counter]][[Bin]][-Rows, ] }} # Remove rows.

      # Determine the sample size of interest - either the number of rows (Samplesize_choice = "Species") or the number of
      # site columns with occurrences + occurring within the present "Fossil_timebins" bin (Samplesize_choice = "Sites"),
      # per element of "Fossil_traitsubset":
      
      if (Samplesize_choice == 'Species') { # If the number of rows is the sample size of interest...
        Fossil_traitsubset_size[[Counter]][[Bin]] <- nrow(Fossil_traitsubset[[Counter]][[Bin]]) }
      if (Samplesize_choice == 'Sites') { # If the number of site columns is the sample size of interest...
        Fossil_traitsubset_size[[Counter]][[Bin]] <- sum(sapply(Fossil_traitsubset[[Counter]][[Bin]][, (
          Occurrence_sites_index_calc + 1):ncol(Fossil_herbivores)], sum) > 0 & Fossil_herbivores_temporal$Max..Age >=
          Fossil_timebins[[Bin]][1] & Fossil_herbivores_temporal$Min..Age <= Fossil_timebins[[Bin]][2], na.rm = TRUE) }
      
      # Calculate the metrics for each trait in each time bin. Because the sample sizes in each data subset (where each
      # subset represents a given trait and time bin) are very different from each other (see sample size plots below),
      # calculate these metrics by randomly sampling the number of rows defined in "Sample_size_func" above from each data
      # subset "Reps" times, and taking the mean of those sampled calculations. Then, calculate the standard deviation of
      # those sampled calculations to get a sense of bootstrap confidence. Also, separately calculate the site-to-site
      # standard deviations for each metric, by first calculating each metric at each individual site per data subset,
      # and then determining the SDs from those site-specific values:
      
      Tmp <- replicate(Reps, { Matrix <- Fossil_traitsubset[[Counter]][[Bin]][sample(nrow( # Bootstrap resampling.
        Fossil_traitsubset[[Counter]][[Bin]][Trait]), Sample_size_func, replace = TRUE), Trait]; Metrics(Matrix) })
      
      Fossil_traitmetrics[[Counter]][[Bin]] <- apply(Tmp, 1, mean) # Calculate the samples' mean.
      Fossil_traitmetrics_samples_sd[[Counter]][[Bin]] <- apply(Tmp, 1, sd) # Calculate the samples' sd.
      
      Fossil_traitmetrics_sites_sd[[Counter]][[Bin]] <- apply( # Across all sites...
        
        sapply(colnames(Fossil_traitsubset[[Counter]][[Bin]][which(colnames(Fossil_traitsubset[[Counter]][[
        Bin]]) == 'Amboseli_Pleis'):ncol(Fossil_traitsubset[[Counter]][[Bin]])]), function(Column) { # For each site...
          Metrics(Fossil_traitsubset[[Counter]][[Bin]][Fossil_traitsubset[[Counter]][[Bin]][Column] == 1,
            Trait]) }), # Calculate "Metrics" for the trait of interest across all species occurring at that site.
        
        1, function(Row) { sd(Row, na.rm = TRUE) }) } # Calculate the sites' sd.
    
    # Create lists that mirror "Fossil_traitmetrics" as well as the two "sd"-based variables, derived from the 
    # "Modern_herbivores" dataset, and concatenate those to create "Overall" variables:
    
    Tmp2 <- replicate(Reps, { Matrix <- Modern_herbivores[sample(nrow(Modern_herbivores[Trait]), Sample_size_func,
      replace = TRUE), Trait]; Metrics(Matrix) })
    
    Modern_traitmetrics[[Counter]] <- apply(Tmp2, 1, mean)
    Overall_traitmetrics[[Counter]] <- append(Fossil_traitmetrics[[Counter]], list(Modern_traitmetrics[[Counter]]))
    
    Modern_traitmetrics_samples_sd[[Counter]] <- apply(Tmp2, 1, sd)
    Overall_traitmetrics_samples_sd[[Counter]] <- append(Fossil_traitmetrics_samples_sd[[Counter]],
      list(Modern_traitmetrics_samples_sd[[Counter]]))
    
    Modern_traitmetrics_sites_sd[[Counter]] <- apply(sapply(colnames(Modern_herbivores[which(colnames(
      Modern_herbivores) == 'ABER'):ncol(Modern_herbivores)]), function(Column) { Metrics(
      Modern_herbivores[Modern_herbivores[Column] == 1, Trait]) }), 1, function(Row) { sd(Row, na.rm = TRUE) })
    Overall_traitmetrics_sites_sd[[Counter]] <- append(Fossil_traitmetrics_sites_sd[[Counter]],
      list(Modern_traitmetrics_sites_sd[[Counter]]))
    
    # Plot the mean, median, and standard deviation for each trait in "Overall_traitmetrics" over time as data points and
    # LOESS regression curves. Plot standard deviations of those calculations as error bars:
    
    for (Metric in 1:length(Overall_traitmetrics[[Counter]][[Bin]])) { # For each of mean, median, and sd...

        # Create a dataframe that contains the X, Y, and error/SD values for each plot, spanning all dataset sites:
      
        Plot_df[[Counter]][[Metric]] <- data.frame(X = 1:(length(Fossil_timebins) + 1), # Account for a modern bin too.
          Y = sapply(Overall_traitmetrics[[Counter]], function(M) { M[Metric] }), # Extract the metric per time bin.
          SD_samples = sapply(Overall_traitmetrics_samples_sd[[Counter]], function(M) { M[Metric] }), # Bootstrap sd.
          SD_sites = sapply(Overall_traitmetrics_sites_sd[[Counter]], function(M) { M[Metric] })) # Site-to-site sd.
        
        # Record the primary breakpoint that is associated with that dataframe from breakpoint analysis:
        
        if (exists('X')) { rm(X) }; if (exists('Y')) { rm(Y) } # Remove these so they do not interfere with below.
        Plot_df_breakpoints[[Counter]][[Metric]] <- segmented(lm(Y ~ X, data = Plot_df[[Counter]][[Metric]]),
          seg.Z = ~X, control = seg.control(n.boot = Reps)) # Model for breakpoint.
        Plot_df_breakpoints[[Counter]][[Metric]] <- round(Plot_df_breakpoints[[Counter]][[Metric]]$psi[, 2]) # Actual pt.
        Plot_df_breakpoints_sig[[Counter]][[Metric]] <- davies.test(lm(Y ~ X, data = Plot_df[[Counter]][[Metric]]))$p.value
         
        # Determine the Y tick labels that will be used in the plot. Make them clean, round numbers that are evenly spaced:
  
        Y_ticks_traits[[Counter]][[Metric]] <- seq(min(Plot_df[[Counter]][[Metric]]$Y - 
          Plot_df[[Counter]][[Metric]]$SD_samples, na.rm = TRUE), max(Plot_df[[Counter]][[Metric]]$Y +
          Plot_df[[Counter]][[Metric]]$SD_samples, na.rm = TRUE), length.out = 5)[2:4] # Tick points, based on range.
        Y_ticks_traits[[Counter]][[Metric]] <- round_any(Y_ticks_traits[[Counter]][[Metric]], 0.1, f = ceiling) # Round.
        Y_ticks_traits[[Counter]][[Metric]] <- c(Y_ticks_traits[[Counter]][[Metric]][1], mean(
          Y_ticks_traits[[Counter]][[Metric]][c(1,3)]), Y_ticks_traits[[Counter]][[Metric]][3]) # Evenly space.
                             
        # Create each plot, using the information above:
        
        Trait_plot[[Counter]][[Metric]] <- local({ # Wrap plot in local environment to ensure proper indexing.
          
          Metric <- Metric # Define "Metric" within the local environment.
          
          Plot <- ggplot(Plot_df[[Counter]][[Metric]], aes(X, Y)) + # Select the dataframe and axes for the plot.
            Breaks_clim + Breaks_clim2 + Breaks_clim3 + Breaks_clim4 + Breaks_hom + # Input the rectangles of past events.
            geom_point(size = 3) + # Set the plot to be a scatter plot.
            geom_smooth(color = 'black', span = 0.75) + # Set the plot to have a curve with LOESS smoothing.
            # geom_errorbar(aes(ymin = Y - SD_sites, ymax = Y + SD_sites), width = 0.5, size = 1) + # Site-site error.
            geom_errorbar(aes(ymin = Y - SD_samples, ymax = Y + SD_samples), width = 0.5, size = 1) + # Samples error.
            geom_vline(xintercept = ifelse(is.na(Plot_df_breakpoints[[Counter]][[Metric]]), min(Events_breaks), 
              Plot_df_breakpoints[[Counter]][[Metric]]), linetype = 'solid', size = 3, alpha = 0.5) +
              # Input and plot a vertical line corresponding to the breakpoint determined in breakpoint analysis above.
            labs(y = paste(Traitmetric_names[[Metric]], Trait_names[[Counter]])) + # Set only Y label for now.
            ylim(min(Plot_df[[Counter]][[Metric]]$Y - Plot_df[[Counter]][[Metric]]$SD_samples, na.rm = T),
              max(Plot_df[[Counter]][[Metric]]$Y + Plot_df[[Counter]][[Metric]]$SD_samples, na.rm = T)) + # Y-axis bounds.
            # scale_x_continuous(breaks = seq(length(Fossil_timebins) + 1), labels = c(sapply(Fossil_timebins, 
            #   function(Bin) { paste(as.character(Bin[1]), '-', as.character(Bin[2])) }), '0')) + 
              # Format full X tick labels designed for detail (comment out for publication purposes and use below).
            scale_x_continuous(breaks = Events_breaks, labels = if (Metric == 3) { format(Events, nsmall = 1) } else { 
              rep('', length(Events)) }, expand = c(0.0075,0)) + # (Use this for publication sans above).
            scale_y_continuous(breaks = Y_ticks_traits[[Counter]][[Metric]], labels = function(Lab) {
              sprintf('%.2f', Lab) }) + # Format Y tick label locations (by range of values) and decimals.
            theme(panel.background = element_blank(), axis.line = element_line(color = 'black'), axis.title.x = 
              element_blank(), # axis.text.x = element_text(angle = 45, size = 6, face = 'bold', hjust = 1), # Full labels.
              axis.text.x = element_text(size = 23, angle = 60, hjust = 1, color = 'black'), # (For publication sans above).
              axis.title.y = element_text(size = 31, margin = margin(0,20,0,0)), axis.text.y = element_text(size = 23,
              color = 'black'), plot.margin = unit(c(ifelse(Metric == 1, 0.5, -1.1), 0.5, 0, 0.5), 'lines')) }) } # Format.
    
    # For each trait, arrange the trait's plots (one plot per metric) all on the same panel:
    
    grid.arrange(plot_grid(Trait_plot[[Counter]][[1]], # Trait_plot[[Counter]][[2]], # (Comment to remove median plot).
      Trait_plot[[Counter]][[length(Traitmetric_names)]], ncol = 1, align = 'v'))
    
    # Update the counter to reference the next trait in the next loop iteration:
          
    Counter <- Counter + 1 }; rm(Bin, Metric, Trait, Counter, Rows, Rows_genera, Rows_genera_all, Tmp, Tmp2)

# Update the body mass-related element of "Trait_names":
  
Trait_names[grepl('Mass', Trait_names)] <- 'Body Mass (ln kg)'
  
# Output and analyze plots for the compositions of all herbivore traits of interest together over time:
  
  # Choose the measure of interest that will be used to analyze the plots:
  
  Functional_measure_choice <- 'Func_disp' # Or "Func_rich".
  
  # Set up the single lists that will be used and processed in the "Bin" "for" loop below:
  
  Functional_measure <- list()
  Functional_measure_sampled <- list()

  # Set up the nested list that will be used and processed in the "Bin" "for" loop, accounting for the inclusion of modern
  # data in addition to the time bins in "Fossil_timebins":
  
  Trait_subset_3D <- list()
  for (Index in 1:(length(Fossil_timebins) + 1)) { Trait_subset_3D[[Index]] <- list() }; rm(Index)
  
  # For each time bin, including modern data, create a 3D plot of body mass, hypsodonty index, and loph count values, and
  # calculate the measure assigned in "Functional_measure_choice" of each plot:
  
  library(plot3D) # For using the "scatter3D" function to plot 3D plots.
  library(geometry) # For using the "convhulln" function to calculate convex hull volumes.
  
  Euclidean_dist <- function(Point1, Point2) { # To calculate the Euclidean distance between two points in space...
    Dist <- sapply(1:length(Point1), function(Idx) { (Point1[Idx] - Point2[Idx])^2 }) # Squared differences of pt coords.
    Dist <- sqrt(sum(unlist(Dist))) # Add up those squared difference and take the square root of their sum.
    return(Dist) } # Return the result.

  par(mar = c(1,1,1,1), mfrow = c(3,3)) # Set up plot layout settings.
  
  for (Bin in 1:(length(Fossil_timebins) + 1)) { # For each bin of time in "Fossil_timebins"...
    
    # Consider the last bin to be modern data, and all previous bins to be fossil data from "Fossil_timebins". Subset the
    # data in the bin to only include rows that contain information about all traits of interest:
    
    if (Bin == length(Fossil_timebins) + 1) { Sub <- Modern_herbivores } else { Sub <- Fossil_traitsubset[[1]][[Bin]] }
    Trait_subset_3D[[Bin]] <- Sub[!is.na(Sub['Mass_Kg_Final']) & !is.na(Sub['HYP']) & !is.na(Sub['LOP']), ]
    
    # Make the associated 3D plot for the bin from its corresponding subset defined above. Only show select plots that
    # are separated from each other by a set temporal interval:
    
    if (Bin %in% c(which(sapply(Fossil_timebins, function(B) { round(B, 1)[length(B)] %in% 7:1 })),
      length(Fossil_timebins) + 1)) { # If the bin in question represents a million-year mark or modern times...
      scatter3D(Trait_subset_3D[[Bin]]$'HYP', Trait_subset_3D[[Bin]]$'LOP', Trait_subset_3D[[Bin]]$'Mass_Kg_Final', # Axes.
        xlim = c(min(Fossil_herbivores$'HYP', na.rm = TRUE), max(Fossil_herbivores$'HYP', na.rm = TRUE)), # Axis ranges...
        ylim = c(min(Fossil_herbivores$'LOP', na.rm = TRUE), max(Fossil_herbivores$'LOP', na.rm = TRUE)),
        zlim = c(min(Fossil_herbivores$'Mass_Kg_Final', na.rm = TRUE), max(Fossil_herbivores$'Mass_Kg_Final', na.rm = T)),
        colvar = NULL, # Set the color of the points based upon their depth into the 3D space of the plot...
          col = colorRampPalette(c('orangered3', 'deepskyblue3'))(length(levels(as.factor(Trait_subset_3D[[Bin]]$'HYP' -
            Trait_subset_3D[[Bin]]$'LOP'))))[as.factor(Trait_subset_3D[[Bin]]$'HYP' - Trait_subset_3D[[Bin]]$'LOP')],
        main = ifelse(Bin == length(Fossil_timebins) + 1, '0 Mya', paste(as.character(round(Fossil_timebins[[Bin]][2], 2)),
          '-', as.character(round(Fossil_timebins[[Bin]][1], 2)), 'Mya')), # Set the plot title, based on years it covers.
        xlab = 'Hypsodonty', ylab = 'Loph Count', zlab = 'Body Mass (ln kg)', # Set the axis labels.
        pch = 19, phi = 0) } # Set data-point style, as well as plot orientation.

    # Compute the measure assigned in "Functional_measure_choice" for each 3D plot. Because the sample sizes in each 3D
    # plot are very different from each other, calculate the measure both according to the absolute data points in each
    # plot, as well as by randomly sampling "Sample_size_func" points from each 3D plot 1,000 times and calculating the
    # measure from each of those samples:
    
    if (Functional_measure_choice == 'Func_rich') { # If measuring func richness, or vol of convex hull enwrapping pts...
      Tmp <- Trait_subset_3D[[Bin]][c('HYP', 'LOP', 'Mass_Kg_Final')] # Define 3D plot dataframe.
      Functional_measure[[Bin]] <- convhulln(Tmp, options = 'FA')$vol # Perform the measure on absolute points.
      Functional_measure_sampled[[Bin]] <- replicate(Reps, tryCatch({ convhulln(Tmp[sample(nrow(
        Tmp), Sample_size_func, replace = TRUE), colnames(Tmp)], options = 'FA')$vol }, # Perform the measure on samples.
        error = function(E) {})) } # Bypass errored hull samples, since some samples do not form a viable polygon.
    
    if (Functional_measure_choice == 'Func_disp') { # If measuring func dispersion, or mean distance of pts to centroid...
      Tmp <- Trait_subset_3D[[Bin]][c('HYP', 'LOP', 'Mass_Kg_Final')] # Define 3D plot dataframe.
      Functional_measure[[Bin]] <- mean(apply(Tmp, 1, function(Row) { Euclidean_dist(Row, sapply(Tmp, mean)) }))
      Functional_measure_sampled[[Bin]] <- replicate(Reps, { mean(apply(Tmp[sample(nrow(Tmp), Sample_size_func,
        replace = TRUE), colnames(Tmp)], 1, function(Row) { Euclidean_dist(Row, sapply(Tmp, mean)) })) }) }}
    rm(Sub, Bin, Tmp)
    
  # Plot the sampled measures over time with LOESS regression and sample error bars. Use similar settings for this plot as
  # for previous plots above:
  
    # Create a dataframe for the plot:
    
    Plot_df_functional <- data.frame(X = 1:(length(Fossil_timebins) + 1), Y = unlist(Functional_measure), 
      Y_samples_mean = sapply(Functional_measure_sampled, function(List) { mean(unlist(List)) }),
      Y_samples_sd = sapply(Functional_measure_sampled, function(List) { sd(unlist(List)) }))
    
    # Record the primary breakpoint that is associated with that dataframe from breakpoint analysis:
    
    if (exists('X')) { rm(X) }; if (exists('Y')) { rm(Y) } # Remove these variables so they do not interfere with below.
    Plot_df_breakpoints_functional <- segmented(lm(Y_samples_mean ~ X, data = Plot_df_functional), seg.Z = ~X,
      control = seg.control(n.boot = Reps)) # Model for breakpoint.
    Plot_df_breakpoints_functional <- round(Plot_df_breakpoints_functional$psi[, 2]) # Actual point.
    Plot_df_breakpoints_functional_sig <- davies.test(lm(Y_samples_mean ~ X, data = Plot_df_functional))$p.value
        
    # Determine the Y tick labels that will be used in the plot. Make them clean, round numbers that are evenly spaced:
  
    Y_ticks_functional <- seq(min(Plot_df_functional$Y_samples_mean - Plot_df_functional$Y_samples_sd, na.rm = TRUE),
      max(Plot_df_functional$Y_samples_mean + Plot_df_functional$Y_samples_sd, na.rm = TRUE), length.out = 5)[2:4]
    Y_ticks_functional <- round_any(Y_ticks_functional, 0.1, f = ceiling) # Round the tick labels to the nearest 0.1.
    Y_ticks_functional <- c(Y_ticks_functional[1], mean(Y_ticks_functional[c(1,3)]), Y_ticks_functional[3]) # Evenly space.
    
    # Create the plot. Do not actually plot it here, as it will be plotted later on with taxonomic and phylogenetic plots:

    Plot_functional <- ggplot(Plot_df_functional, aes(X, Y_samples_mean)) + # Select the dataframe and axes for the plot.
      Breaks_clim + Breaks_clim2 + Breaks_clim3 + Breaks_clim4 + Breaks_hom + # Input the rectangles of past events.
      geom_point(size = 3) + # Set the plot to be a scatter plot.
      geom_smooth(color = 'black', span = 0.75) + # Set the plot to have a curve with LOESS smoothing.
      geom_errorbar(aes(ymin = Y_samples_mean - Y_samples_sd, ymax = Y_samples_mean + Y_samples_sd), width = 0.5,
        size = 1) + # Plot error bars that reflect the error across all bootstrap samples.
      geom_vline(xintercept = ifelse(is.na(Plot_df_breakpoints_functional), min(Events_breaks), 
        Plot_df_breakpoints_functional), linetype = 'solid', size = 3, alpha = 0.5) + # Plot the breakpoint.
      labs(y = ifelse(Functional_measure_choice == 'Func_disp', 'Functional\n(Distance)',
        'Functional\n(Volume)')) + # Set only Y axis label for now.
      ylim(min(Plot_df_functional$Y_samples_mean - Plot_df_functional$Y_samples_sd, na.rm = TRUE),
        max(Plot_df_functional$Y_samples_mean + Plot_df_functional$Y_samples_sd, na.rm = TRUE)) + # Y-axis bounds.
      # scale_x_continuous(breaks = seq(length(Fossil_timebins) + 1), labels = c(sapply(Fossil_timebins, 
        # function(Bin) { paste(as.character(Bin[1]), '-', as.character(Bin[2])) }), '0')) + 
        # Format full X tick labels designed for detail (comment out for publication purposes and use below).
      scale_x_continuous(breaks = Events_breaks, labels = format(Events, nsmall = 1), expand = c(0.0075,0)) + # (For pub).
      scale_y_continuous(breaks = Y_ticks_functional, labels = function(Lab) { sprintf('%.2f', Lab) }) + # Format Y ticks.
      theme(panel.background = element_blank(), axis.line = element_line(color = 'black'), axis.title.x = 
        element_blank(), # axis.text.x = element_text(angle = 45, size = 6, face = 'bold', hjust = 1), # For full labels.
        axis.text.x = element_text(size = 20, angle = 60, hjust = 1, color = 'black'), # (For publication sans above).
        axis.title.y = element_text(size = 20, margin = margin(0,15,0,0)), axis.text.y = element_text(size = 20,
        color = 'black'), plot.margin = unit(c(-1.5, 0.5, 0, 0.5), 'lines')) # Format plot elements.
    
###########################################################################################################################
###########################################################################################################################
  
# PART VI: ANALYSIS OF HERBIVORE SPECIES RICHNESS THROUGH TIME (TAXONOMIC DIVERSITY)
    
# Determine the number of species occurring within each time bin in "Fossil_timebins". To account for sampling bias,
# implement coverage-based rarefaction:

  # Set up the dataframe that will hold richness values:
    
  Species_richness <- data.frame(matrix(nrow = length(Fossil_timebins) + 1, ncol = 1))
  colnames(Species_richness) <- 'Richness'

  # Perform richness calculations across all species and time bins. To do so, perform coverage-based rarefaction on the
  # species occurring in time bins. Doing so requires setting up an incidence list, where each element represents a time
  # bin. Within each element is a vector, where the vector's first value is the total number of sites that occur during
  # the time bin's duration, and the remaining values refer to the number of sites at which each species occurs, out of
  # the total number for the time bin (incidence frequencies). Put calculations of richness in "Species_richness":
  
  library(iNEXT) # For performing the rarefaction.
  
  Species_richness_incidences <- list() # List to hold species incidence frequency data across all time bins.

  for (Bin in 1:nrow(Species_richness)) { # For each time bin...
    
    if (Bin != nrow(Species_richness)) { # If the modern time bin is not being addressed...
      Cols <- which(colnames(Fossil_traitsubset[[1]][[Bin]]) == 'Amboseli_Pleis'):ncol(Fossil_traitsubset[[1]][[Bin]])
      Cols <- Cols[which(Fossil_herbivores_temporal$Max..Age >= Fossil_timebins[[Bin]][1] &
        Fossil_herbivores_temporal$Min..Age <= Fossil_timebins[[Bin]][2])] # Which sites occur during the time bin?
      Data <- as.matrix(Fossil_traitsubset[[1]][[Bin]][, Cols]) # Occurrence matrix of bin's component species and sites.
    } else { # If the modern time bin is being addressed, create its occurrence matrix...
      Data <- as.matrix(Modern_herbivores[, which(colnames(Modern_herbivores) == 'ABER'):ncol(Modern_herbivores)]) }

    if (ncol(Data) == 1) { # If the time bin is associated with only one site...
      Species_richness_incidences[[Bin]] <- rep(1, nrow(Data) + 1) # Turn its occurrence vector into a frequency vector.
      } else { # Or else, if the time bin is associated with more than one site...
      Species_richness_incidences[[Bin]] <- unname(as.incfreq(Data)) # Turn its occurrence matrix into frequency vector.
      Species_richness_incidences[[Bin]] <- Species_richness_incidences[[Bin]][Species_richness_incidences[[Bin]] != 0] }

    }; rm(Bin, Cols, Data)
  
  names(Species_richness_incidences) <- as.character(1:length(Species_richness_incidences)) # Give names to time bins.
  Species_richness_estimates <- iNEXT(Species_richness_incidences, datatype = 'incidence_freq') # Estimate richness.
    # NOTE: the default bootstrap resampling number for this is 50, but 1000 resamples produces the same result.
  Species_richness$Richness <- Species_richness_estimates$AsyEst$Estimator[
    Species_richness_estimates$AsyEst$Diversity == 'Species richness'] # Integrate estimates into "Species_richness".
  Species_richness$Richness_se <- Species_richness_estimates$AsyEst$s.e.[
    Species_richness_estimates$AsyEst$Diversity == 'Species richness'] # Integrate standard errors into "Species_richness".

# Create a plot of the values in "Species_richness" over time, using similar settings as in previous plots above. Do not
# actually plot it here, as it will be plotted later on with functional and phylogenetic diversity plots:

  # Record the primary breakpoint that is associated with "Species_richness" from breakpoint analysis:

  Species_richness$Richness_Index <- 1:nrow(Species_richness) # Create an index column.

  Plot_df_breakpoints_richness <- segmented(lm(Richness ~ Richness_Index, data = Species_richness), 
    seg.Z = ~Richness_Index, control = seg.control(n.boot = Reps)) # Model for breakpoint.
  Plot_df_breakpoints_richness <- round(Plot_df_breakpoints_richness$psi[, 2]) # Actual point.
  Plot_df_breakpoints_richness_sig <- davies.test(lm(Richness ~ Richness_Index, data = Species_richness))$p.value

  # Determine the Y tick labels that will be used in the plot. Make them clean, round numbers that are evenly spaced:
  
  Y_ticks_richness <- seq(min(Species_richness, na.rm = TRUE), max(Species_richness, na.rm = TRUE), length.out = 5)[2:4]
  Y_ticks_richness <- round_any(Y_ticks_richness, 10, f = ceiling) # Round the tick labels up to the nearest 10.
  Y_ticks_richness <- c(Y_ticks_richness[1], mean(Y_ticks_richness[c(1,3)]), Y_ticks_richness[3]) # Evenly space.
  
  # Make the plot, including the breakpoint:
  
  Plot_taxonomic <- ggplot(Species_richness, aes(Richness_Index, Richness)) + # Select the dataframe and axes for the plot.
    Breaks_clim + Breaks_clim2 + Breaks_clim3 + Breaks_clim4 + Breaks_hom + # Input the rectangles of past events.
    geom_point(size = 3) + # Set the plot to be a scatter plot.
    geom_smooth(color = 'black', span = 0.75) + # Set the plot to have a curve with LOESS smoothing.
    geom_errorbar(aes(ymin = Richness - Richness_se, ymax = Richness + Richness_se), width = 0.5,
      size = 1) + # Plot error bars that reflect the error across richness estimates.
    geom_vline(xintercept = ifelse(is.na(Plot_df_breakpoints_richness), min(Events_breaks),
      Plot_df_breakpoints_richness), linetype = 'solid', size = 3, alpha = 0.5) + # Plot breakpoint.
    labs(y = 'Taxonomic\n(Richness)') + # Set Y axis label.
    # scale_x_continuous(breaks = seq(length(Fossil_timebins) + 1), labels = c(sapply(Fossil_timebins, 
      # function(Bin) { paste(as.character(Bin[1]), '-', as.character(Bin[2])) }), '0')) + 
      # Format full X tick labels designed for detail (comment out for publication purposes and use below).
    scale_x_continuous(breaks = Events_breaks, labels = rep('', length(Events)), expand = c(0.0075,0)) + # (For pub).
    scale_y_continuous(breaks = Y_ticks_richness, labels = function(Lab) { sprintf('%.0f', Lab) }) + # Format Y tick labs.
    theme(panel.background = element_blank(), axis.line = element_line(color = 'black'), axis.title.x = 
      element_blank(), # axis.text.x = element_text(angle = 45, size = 6, face = 'bold', hjust = 1), # For full labels.
      axis.text.x = element_text(size = 20, angle = 60, hjust = 1), # (Use this for publication sans above).
      axis.title.y = element_text(size = 20, margin = margin(0,15,0,0)), axis.text.y = element_text(size = 20,
      color = 'black'), plot.margin = unit(c(0.5, 0.5, -0.5, 0.5), 'lines')) # Format plot elements.
  
###########################################################################################################################
###########################################################################################################################
  
# PART VII: ANALYSIS OF HERBIVORE SPECIES PHYLOGENETIC DIVERSITY THROUGH TIME (PHYLOGENETIC DIVERSITY)
  
# Identify all genera in both "Fossil_herbivores" and "Modern_herbivores" to be used below:

Overall_genera <- unique(sapply(strsplit(c(as.character(Fossil_herbivores$Binomial_Final), as.character(
  Modern_herbivores$Binomial_Final)), ' '), '[[', 1)) # All unique genera across all herbivores.

# For each full mammalian phylogenetic tree downloaded from Vertlife, subset it down to just the genera in 
# "Overall_genera", and in the process, convert it from being at the species to the genus level (perform this operation
# once, and comment it out thereafter, as it only needs to be performed once):

library(phytools) # For conducting phylogenetic analyses below.
  
# setwd('./Overall_Herbivores_Phylogeny_Vertlife') # Set the working directory to the folder with the trees.
# for (Tree in 1:length(list.files())) { # For each tree...
#   Tmp <- read.tree(list.files()[Tree]) # Read the tree in.
#   Tmp <- keep.tip(Tmp, Tmp$tip.label[gsub('_.*', '', Tmp$tip.label) %in% Overall_genera]) # Subset to "Overall_genera".
#   Tmp <- drop.tip(Tmp, Tmp$tip.label[duplicated(gsub('_.*', '', Tmp$tip.label))]) # Keep only one member of each genus.
#   Tmp$tip.label <- gsub('_.*', '', Tmp$tip.label) # Have all tree tips only show genus names.
#   write.tree(Tmp, file = paste('Tree_', as.character(Tree), '.tre', sep = '')) }; rm(Tree, Tmp) # Save the modified tree.
# setwd('../') # Restore the working directory to its previous location.

# Read in a sample phylogenetic tree that was downloaded and processed from Vertlife, and plot it:

Overall_phylo_tree_sample <- read.tree('Overall_Herbivores_Phylogeny_Vertlife/Tree_1.tre')
par(mar = c(0.1,0.1,0.1,0.1), mfrow = c(1,1)); plot(Overall_phylo_tree_sample)

# For each phylogenetic tree from Vertlife, impute the genera in "Overall_genera" that are not included within it,
# ensuring that the imputation occurs according to the taxonomic affiliations of each genus:

Overall_phylo_trees <- list() # List to hold the finalized, imputed trees.

Overall_herbivores_taxonomy <- rbind(Fossil_herbivores[, c('Order', 'Family', 'Binomial_Final')], Modern_herbivores[,
  c('Order', 'Family', 'Binomial_Final')]) # Taxonomies of all species, which will be used in the "for" loop below.
Overall_herbivores_taxonomy <- cbind(Overall_herbivores_taxonomy, c(as.character(Fossil_herbivores$Subfamily),
  rep(NA, nrow(Modern_herbivores))), c(as.character(Fossil_herbivores$Tribe), rep(NA, nrow(Modern_herbivores))))
colnames(Overall_herbivores_taxonomy) <- c('Order', 'Family', 'Binomial_Final', 'Subfam', 'Tribe') # Add in fossil cols.
Overall_herbivores_taxonomy[Overall_herbivores_taxonomy == ''] <- NA # Remove blanks from "Subfamily" and "Tribe".

Overall_genera_nonunique <- gsub(' .*', '', Overall_herbivores_taxonomy$Binomial_Final) # Genera of all taxa.
Overall_genera_missing <- Overall_genera[!(Overall_genera %in% Overall_phylo_tree_sample$tip.label)] # Genera to impute.
Overall_genera_known <- Overall_phylo_tree_sample$tip.label # Genera that were in the tree before any imputation.
Overall_genera_known_tax <- Overall_herbivores_taxonomy[match(Overall_genera_known, Overall_genera_nonunique),
  c('Order', 'Family', 'Subfam', 'Tribe')] # Taxonomies of the genera that were in the tree before any imputation.

setwd('./Overall_Herbivores_Phylogeny_Vertlife') # Set the working directory to the folder with the trees.

for (Tree in 1:length(list.files())) { # For each tree...
  
  # Read the tree in, and index it into the list where all trees will be held. Then, force the tree to be ultrametric,
  # such that all of its leaves are equidistant from the root (they already are, but R does not recognize that). This is
  # necessary for ensuring that branch lengths are properly calibrated when imputing new branches into the tree below:
  
  Overall_phylo_trees[[Tree]] <- read.tree(list.files()[Tree])
  Overall_phylo_trees[[Tree]] <- force.ultrametric(Overall_phylo_trees[[Tree]])
  
  # Impute each genus that is missing into the tree, basing the imputation upon the genus' taxonomic information: 
    
  for (Genus in Overall_genera_missing) { # For each genus that is missing and should be imputed...
    
    # Determine the genus' taxonomic hierarchy, and determine which known genus/genera in the tree share the same
    # taxonomy, starting with finer-resolution taxonomic levels and moving up to higher ones if necessary:
    
    Genus_tax <- Overall_herbivores_taxonomy[which(Genus == Overall_genera_nonunique)[1], c('Order', 'Family',
      'Subfam', 'Tribe')] # Taxonomic hierarchy of genus of interest.

    if (Genus == 'Alcelaphini') { Genera_common <- 'Alcelaphus' } else if ( # NA if "Alcelaphini" genus removed from data.
      !is.na(Genus_tax[, 'Tribe']) & sum(Genus_tax[, 'Tribe'] == Overall_genera_known_tax[, 'Tribe'], na.rm = T) > 0) { 
        Genera_common <- Overall_genera_known[Genus_tax[, 'Tribe'] == Overall_genera_known_tax[, 'Tribe']] } else if (
      !is.na(Genus_tax[, 'Subfam']) & sum(Genus_tax[, 'Subfam'] == Overall_genera_known_tax[, 'Subfam'], na.rm = T) > 0) { 
        Genera_common <- Overall_genera_known[Genus_tax[, 'Subfam'] == Overall_genera_known_tax[, 'Subfam']] } else if (  
      !is.na(Genus_tax[, 'Family']) & sum(Genus_tax[, 'Family'] == Overall_genera_known_tax[, 'Family'], na.rm = T) > 0) { 
        Genera_common <- Overall_genera_known[Genus_tax[, 'Family'] == Overall_genera_known_tax[, 'Family']] } else if (
      !is.na(Genus_tax[, 'Order']) & sum(Genus_tax[, 'Order'] == Overall_genera_known_tax[, 'Order'], na.rm = T) > 0) { 
        Genera_common <- Overall_genera_known[Genus_tax[, 'Order'] == Overall_genera_known_tax[, 'Order']] }
    
    Genera_common <- Genera_common[!(is.na(Genera_common))] # Remove and "NA" values from "Genera_common".
    
    # Perform the imputation by adding the genus of interest to the node representing the most recent common ancestor of
    # the genera in "Genera_common", or to the tip representing the singular genus in "Genera_common":
    
    if (length(Genera_common) > 1) { # If the genus of interest shares taxonomy with more than one other genus...
      Where <- getMRCA(Overall_phylo_trees[[Tree]], Genera_common) # Find the common ancestor node of those other genera.
      Position <- 0 } # Add the genus to the common ancestor node exactly, as opposed to upstream or downstream of it.
    if (length(Genera_common) == 1) { # If the genus of interest shares taxonomy with only one other genus...
      Where <- which(Overall_phylo_trees[[Tree]]$tip.label == Genera_common) # Find the leaf of that genus.
      Position <- 2 } # Add the genus with a branch length that is reasonable for a leaf.
    
    Overall_phylo_trees[[Tree]] <- bind.tip(Overall_phylo_trees[[Tree]], Genus, where = Where, position = Position) # Add.
    
    # The imputation step above generates polytomies (nodes with more than two descendents), so resolve those by randomly
    # collapsing them into multiple nodes, each with a pair of descendents:
    
    Overall_phylo_trees[[Tree]] <- multi2di(Overall_phylo_trees[[Tree]]) }}
rm(Tree, Genus, Genus_tax, Genera_common, Where, Position)

setwd('../') # Restore the working directory to its previous location.

# Determine the phylogenetic diversity of all species occurring in each time bin in "Fossil_timebins", as well as in
# modern times. Calculate the sum of the total phylogenetic branch lengths connecting all species (based on Faith, 1992)
# for each of a group of bootstrap resampled taxa per time bin, to account for sampling bias, and across a group of
# sampled phylogenetic trees from those collected and processed above, to account for tree and imputing uncertainty.
# Repeat this analysis for all taxa, as well as separately for each of other groups of taxa of interest for comparison.
# (perform this operation once and comment it out thereafter with its key outputs saved, in order to save time):

# Phylo_taxa <- 'All' # "All", "Ruminants" or "Non-Ruminants" - taxa to analyze?
# if (Phylo_taxa == 'All') { Phylo_genera <- Overall_genera } else if (Phylo_taxa == 'Ruminants') { # All genera.
#   Phylo_genera <- unique(Overall_genera_nonunique[Overall_herbivores_taxonomy$Family == 'Bovidae']) } else { # Ruminants.
#   Phylo_genera <- unique(Overall_genera_nonunique[Overall_herbivores_taxonomy$Family != 'Bovidae']) } # Non-Ruminants.
# 
# Sample_size_phylo_taxa <- sum(unique(sapply(strsplit(Fossil_traitsubset[[1]][[1]][, 'Binomial_Final'], ' '),
#   '[[', 1)) %in% Phylo_genera) # Set the sample size of genera equal to the number of genera in the oldest time bin.
# Sample_size_phylo_trees <- length(Overall_phylo_trees) # Set the sample size of trees to be considered.
# Fossil_herbivores_phylo <- list() # List to hold subsets of taxa for each time bin.
# Plot_df_phylo <- data.frame(X = 1:(length(Fossil_timebins) + 1), Y = NA, SD_samples = NA) # Dataframe for plotting.
# 
# for (Bin in 1:(length(Fossil_timebins) + 1)) { # For each age/time bin, including modern times...
# 
#   if (Bin != length(Fossil_timebins) + 1) { # If the age in question is not modern...
# 
#     # Subset "Fossil_herbivores_phylo" to include only those genera that occur in each time bin:
# 
#     Fossil_herbivores_phylo[[Bin]] <- Fossil_traitsubset[[1]][[Bin]][, 'Binomial_Final']
# 
#     # Modify "Fossil_herbivores_phylo" to only include names of genera:
# 
#     Fossil_herbivores_phylo[[Bin]] <- unique(sapply(strsplit(Fossil_herbivores_phylo[[Bin]], ' '), '[[', 1))
# 
#     # Reduce "Fossil_herbivores_phylo" down to only those genera found in "Phylo_genera" (ruminant vs. non-ruminant):
# 
#     Fossil_herbivores_phylo[[Bin]] <- Fossil_herbivores_phylo[[Bin]][Fossil_herbivores_phylo[[Bin]] %in% Phylo_genera]
# 
#     # Reduce "Fossil_herbivores_phylo" down to only those remaining genera found in "Overall_phylo_trees":
# 
#     Fossil_herbivores_phylo[[Bin]] <- Fossil_herbivores_phylo[[Bin]][Fossil_herbivores_phylo[[Bin]] %in%
#       Overall_phylo_trees[[1]]$tip.label]
# 
#     # Bootstrap resample the genera in "Fossil_herbivores_phylo", calculate the mean sum of branch lengths for each
#     # sample across a set of sampled trees, and then calculate the mean and standard deviation of those mean sum
#     # calculations. Finally, insert the resultant mean and standard deviation values into "Plot_df_phylo":
# 
#     Tmp <- replicate(Reps, { # Bootstrap resample genera X times, calculating a mean sum per sample.
#       Sampled_genera <- sample(Fossil_herbivores_phylo[[Bin]], Sample_size_phylo_taxa, replace = TRUE)
#       mean(sapply(sample(1:length(Overall_phylo_trees), Sample_size_phylo_trees, replace = FALSE), function(Tree) {
#         Sampled_tree_clipped <- keep.tip(Overall_phylo_trees[[Tree]], Sampled_genera) # Clip tree to the sampled genera.
#         return(sum(Sampled_tree_clipped$edge.length)) })) }) # Perform sum for each sampled tree, and find mean of sums.
# 
#     Plot_df_phylo[Bin, c(2:ncol(Plot_df_phylo))] <- c(mean(Tmp, na.rm = TRUE), sd(Tmp, na.rm = TRUE)) }
# 
#   else { # If the age in question is modern...
# 
#     # Perform the analogous bootstrap analysis for "Modern_herbivores" as was done with the fossil data:
# 
#     Modern_herbivores_phylo <- unique(sapply(strsplit(as.character(Modern_herbivores$Binomial_Final), ' '), '[[', 1))
#     Modern_herbivores_phylo <- Modern_herbivores_phylo[Modern_herbivores_phylo %in% Phylo_genera] # Subset genera.
#     Modern_herbivores_phylo <- Modern_herbivores_phylo[Modern_herbivores_phylo %in% Overall_phylo_trees[[1]]$tip.label]
# 
#     Tmp <- replicate(Reps, { # Bootstrap resample genera X times, calculating a mean sum per sample.
#       Sampled_genera <- sample(Modern_herbivores_phylo, Sample_size_phylo_taxa, replace = TRUE)
#       mean(sapply(sample(1:length(Overall_phylo_trees), Sample_size_phylo_trees, replace = FALSE), function(Tree) {
#         Sampled_tree_clipped <- keep.tip(Overall_phylo_trees[[Tree]], Sampled_genera) # Clip tree to the sampled genera.
#         return(sum(Sampled_tree_clipped$edge.length)) })) }) # Perform sum for each sampled tree, and find mean of sums.
# 
#     Plot_df_phylo[Bin, c(2:ncol(Plot_df_phylo))] <- c(mean(Tmp, na.rm = TRUE), sd(Tmp, na.rm = TRUE)) }}; rm(Bin, Tmp)

Plot_df_phylo_dir <- paste('Overall_Herbivores_Phylogeny_Output_Sites=', Fossil_sites_choose, '/', sep = '') # Directory.
# saveRDS(Plot_df_phylo, file = paste(Plot_df_phylo_dir, 'Test.rds', sep = '')) # Save the dataframe to that directory.
Plot_df_phylo <- readRDS(paste(Plot_df_phylo_dir, 'Overall_Herbivores_Phylo_Plotdf_250000-Yr-Non-Moving-Av_Since7-5.rds',
  sep = '')) # Read in the dataframe from that directory.

# Record the breakpoints that are associated with the "Plot_df_phylo" dataframe from breakpoint analysis:

if (exists('X')) { rm(X) }; if (exists('Y')) { rm(Y) } # Remove these variables so they do not interfere with below.
Plot_df_breakpoints_phylo <- segmented(lm(Y ~ X, data = Plot_df_phylo), seg.Z = ~X, control = seg.control(n.boot = Reps))
Plot_df_breakpoints_phylo <- round(Plot_df_breakpoints_phylo$psi[, 2])
Plot_df_breakpoints_phylo_sig <- davies.test(lm(Y ~ X, data = Plot_df_phylo))$p.value # Check if breakpoint needed.
  
# Determine the Y tick labels that will be used in the plot below. Make them clean, round numbers that are evenly spaced:

Y_ticks_phylo <- seq(min(Plot_df_phylo$Y - Plot_df_phylo$SD_samples, na.rm = TRUE), max(Plot_df_phylo$Y +
  Plot_df_phylo$SD_samples, na.rm = TRUE), length.out = 5)[2:4] # Make baseline labels.
Y_ticks_phylo <- round_any(Y_ticks_phylo, 10, f = ceiling) # Round the tick labels to the nearest ten.
Y_ticks_phylo <- c(Y_ticks_phylo[1], mean(Y_ticks_phylo[c(1,3)]), Y_ticks_phylo[3]) # Evenly space.

# Plot the values in "Plot_df_phylo", using similar settings as the various LOESS-regression plots made above:
  
Plot_phylogenetic <- ggplot(Plot_df_phylo, aes(X, Y)) + # Select the dataframe and axes for the plot.
  Breaks_clim + Breaks_clim2 + Breaks_clim3 + Breaks_clim4 + Breaks_hom + # Input the rectangles of past events.
  geom_point(size = 3) + # Set the plot to be a scatter plot.
  geom_smooth(color = 'black', span = 0.75) + # Set the plot to have a curve with LOESS smoothing.
  geom_errorbar(aes(ymin = Y - SD_samples, ymax = Y + SD_samples), width = 0.5, size = 1) + # Bootstrap error.
  geom_vline(xintercept = ifelse(is.na(Plot_df_breakpoints_phylo), min(Events_breaks), Plot_df_breakpoints_phylo), 
    linetype = 'solid', size = 3, alpha = 0.5) + # Plot the breakpoint.
  labs(y = 'Phylogenetic\n(Length)') + # Set just the Y label for now.
  # scale_x_continuous(breaks = seq(length(Fossil_timebins) + 1), labels = c(sapply(Fossil_timebins, 
    # function(Bin) { paste(as.character(Bin[1]), '-', as.character(Bin[2])) }), '0')) + 
    # Format full X tick labels designed for detail (comment out for publication purposes and use below).
  scale_x_continuous(breaks = Events_breaks, labels = rep('', length(Events)), expand = c(0.0075,0)) + # (For pub).
  scale_y_continuous(breaks = Y_ticks_phylo, labels = function(Lab) { sprintf('%.0f', Lab) }) + # Format Y tick labels.
  theme(panel.background = element_blank(), axis.line = element_line(color = 'black'), axis.title.x = 
    element_blank(), # axis.text.x = element_text(angle = 45, size = 6, face = 'bold', hjust = 1), # For full labels.
    axis.text.x = element_text(size = 20, angle = 60, hjust = 1), # (Use this for publication sans above).
    axis.title.y = element_text(size = 20, margin = margin(0,15,0,0)), axis.text.y = element_text(size = 20,
    color = 'black'), plot.margin = unit(c(-0.5, 0.5, 0.5, 0.5), 'lines')) # Format plot elements.

# Now that all three of "Plot_functional", "Plot_taxonomic", and "Plot_phylogenetic" have been made, plot them together.
# Also plot "Plot_cover" (grassy or woody) with these plots:

grid.arrange(plot_grid(Plot_cover, Plot_taxonomic, Plot_phylogenetic, Plot_functional, ncol = 1, align = 'v'))

# Make a plot like "Plot_functional", except divided into different groups of taxa:

setwd(paste('./', Plot_df_phylo_dir, sep = '')) # Set the working directory to the folder containing the plot data.

Plot_phylogenetic_groups <- lapply(c('_Ruminants', '_Non-Ruminants'), function(Group) { # For each group...
  
  Plot_data <- readRDS(list.files()[grep(Group, list.files())]) # Read in group's data.
  
  Plot_breakpoints <- segmented(lm(Y ~ X, data = Plot_data), seg.Z = ~X, control = seg.control(n.boot = Reps))
  Plot_breakpoints <- round(Plot_breakpoints$psi[, 2]) # Obtain breakpoints for the plot below.
  
  Y_ticks <- seq(min(Plot_data$Y - Plot_data$SD_samples, na.rm = TRUE), max(Plot_data$Y + Plot_data$SD_samples,
    na.rm = TRUE), length.out = 5)[2:4] # Baseline Y ticks.
  Y_ticks <- round_any(Y_ticks, 1, f = ceiling) # Round the tick labels.
  Y_ticks <- c(Y_ticks[1], mean(Y_ticks[c(1,3)]), Y_ticks[3]) # Evenly space.
  
  ggplot(Plot_data, aes(X, Y)) + # Select the dataframe and axes for the plot.
    Breaks_clim + Breaks_clim2 + Breaks_clim3 + Breaks_clim4 + Breaks_hom + # Input the rectangles of past events.
    geom_point(size = 3) + # Set the plot to be a scatter plot.
    geom_smooth(color = 'black', span = 0.75) + # Set the plot to have a curve with LOESS smoothing.
    geom_errorbar(aes(ymin = Y - SD_samples, ymax = Y + SD_samples), width = 0.5, size = 1) + # Bootstrap error.
    geom_vline(xintercept = ifelse(is.na(Plot_breakpoints), min(Events_breaks), Plot_breakpoints), 
      linetype = 'solid', size = 3, alpha = 0.5) + # Plot the breakpoint.
    labs(y = 'Phylogenetic\n(Length)') + # Set just the Y label for now.
    scale_x_continuous(breaks = Events_breaks, labels = if (Group == '_Non-Ruminants') { 
      format(Events, nsmall = 1) } else { rep('', length(Events)) }, expand = c(0.0075,0)) + # Format X tick labels.
    scale_y_continuous(breaks = Y_ticks, labels = function(Lab) { sprintf('%.0f', Lab) }) + # Format Y tick labels.
    theme(panel.background = element_blank(), axis.line = element_line(color = 'black'), axis.title.x = 
      element_blank(), axis.text.x = element_text(size = 20, angle = 60, hjust = 1, color = 'black'), axis.title.y =
      element_text(size = 20, margin = margin(0,15,0,0)), axis.text.y = element_text(size = 20, color = 'black'),
      plot.margin = unit(c(ifelse(Group == '_Non-Ruminants', -1.5, 0.5), 0.5, 0.5, 0.5), 'lines')) }) # Format.

grid.arrange(plot_grid(Plot_phylogenetic_groups[[1]], Plot_phylogenetic_groups[[2]], ncol = 1, align = 'v'))

setwd('../') # Reset the working directory.

###########################################################################################################################
###########################################################################################################################
  
# PART VIII: COLLECTION OF ENVIRONMENTAL CONDITIONS DATA FOR ECOMETRIC ANALYSIS

# NOTE: a previous version of our ecometric analyses used climate data, as opposed to vegetation cover data. Much of the
# code that follows was written based on the climate data used. However, we have since incorporated our vegetation data 
# into this code, such that our downstream analyses work properly and are valid when only considering vegetation data
# alone. Thus, the code that follows gives reference to climate data, which is necessary for the code to run appropriately,
# but these references to climate can be ignored, as all analyses based on vegetation cover data are valid.

# Read in the libraries that will be used for plotting and analyzing environmental and ecometric data:

library(raster); library(sp); library(rgdal); library(rgeos); library(colorspace) # For analyzing spatial data.
library(fields) # For adding legends below/for using the "image.plot" function below.

# Create spatial objects containing the properly projected coordinates of all fossil and modern sites analyzed. Then,
# plot the fossil and modern site locations:

Fossil_sitelocations <- Fossil_herbivores_cover[grepl(Fossil_sites_term, Fossil_herbivores_cover$Derived_From),
  c('Longitude', 'Latitude')] # Coordinates of fossil sites.
Fossil_sitelocations <- SpatialPoints(Fossil_sitelocations, proj4string = crs(
  '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')) # Coordinates as spatial points.

Modern_sitelocations <- Modern_herbivores_cover[, c('Longitude', 'Latitude')] # Coordinates of modern sites.
Modern_sitelocations <- SpatialPoints(Modern_sitelocations, proj4string = crs(
  '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')) # Coordinates as spatial points.

par(mfrow = c(1,1), mar = c(2.5,2.5,2.5,2.5)) # Set plot layout settings.
Africa_mask <- shapefile('Africa_Mask/Africa.shp') # Read in a shapefile of Africa.
Africa_mask <- Africa_mask[Africa_mask@data$CONTINENT == 'Africa', ] # Subset to just Africa.

plot(Africa_mask, col = 'wheat2', axes = TRUE, cex.axis = 2) # Plot the shapefile/map of Africa.
plot(Modern_sitelocations, add = TRUE, col = 'cornflowerblue', pch = 19, cex = 2) # Plot the locations of modern sites.

Tmp <- Africa_mask # Replicate "Africa_mask" temporarily.
Tmp <- gBuffer(Tmp, byid = TRUE, width = 0) # Put a buffer around "Tmp" geometries to avoid cropping error.
Tmp <- crop(Tmp, extent(32.58333, 41.91667, -4.791667, 12.41667)) # Crop to East Africa.
plot(Tmp, col = 'wheat2', axes = FALSE, cex.axis = 2) # Plot the cropped map of Africa.
lines(spPolygons(extent(Tmp), crs = crs(Tmp))) # Put a box around the plot.
Tmp2 <- Fossil_herbivores_cover[grepl(Fossil_sites_term, Fossil_herbivores_cover$Derived_From), ] # Temporary dataset.
plot(Fossil_sitelocations, add = TRUE, col = ifelse(Tmp2$Mean_Age >= 5.3, 'dodgerblue3', ifelse(Tmp2$Mean_Age >= 2.6,
  'aquamarine3', ifelse(Tmp2$Mean_Age >= 0.01,'deeppink3', 'darkorchid3'))), pch = ifelse(Tmp2$Mean_Age >= 5.3, 15, ifelse(
  Tmp2$Mean_Age >= 2.6, 17, ifelse(Tmp2$Mean_Age >= 0.01, 18, 19))), cex = 3) # Plot fossil sites, distinguished by period.
rm(Tmp, Tmp2) # No longer needed.

# Create a dataframe in which each site in "Fossil_herbivores" in each time period of interest is represented by its mean,
# median, standard deviation, skewness, and kurtosis of its trait values, as well as its associated number of species,
# MAT and AP values, and woody or grassy cover fractions:

Traitmetric_fullnames <- c('Mean Body Mass (ln kg)', 'Median Body Mass (ln kg)', 'SD Body Mass (ln kg)',
  'Body Mass Skewness', 'Body Mass Kurtosis', 'Mean Hypsodonty Index', 'Median Hypsodonty Index', 
  'SD Hypsodonty Index', 'Hypsodonty Index Skewness', 'Hypsodonty Index Kurtosis', 'Mean Loph Count', 
  'Median Loph Count', 'SD Loph Count', 'Loph Count Skewness', 'Loph Count Kurtosis') # For reference.
Traitmetric_abbrevnames <- c('Mass_mean', 'Mass_median', 'Mass_sd', 'Mass_sk', 'Mass_ku', 'HYP_mean', 'HYP_median',
  'HYP_sd', 'HYP_sk', 'HYP_ku', 'LOP_mean', 'LOP_median', 'LOP_sd', 'LOP_sk', 'LOP_ku') # For reference as well.
Traitmetric_number <- sum(grepl('Body Mass', Traitmetric_fullnames)) # Number of metrics.

Fossil_sites_metrics <- list() # List to hold metric values.
Fossil_sites_specnum <- list() # List to hold species sample size values.

Time <- c(3.3, 3.205, 0.787, 0.13, 0.021, 0.017, 0.015, 0.013, 0.012, 0.008) # Used in our previous climate analyses.
Time <- round(sort(unique(c(Time, seq(7.5, 0.008, -0.1))), decreasing = TRUE), 3) # Used in previous analyses.

Fossil_sites_odd <- colnames(Fossil_herbivores[which(colnames(Fossil_herbivores) == 'Amboseli_Pleis'):ncol(
  Fossil_herbivores)]) # Initiate a variable that will contain sites whose age ranges overlap with no age in "Time".
Fossil_sites_odd <- Fossil_sites_odd[sapply(1:nrow(Fossil_herbivores_temporal), function(Row) { # Find those odd sites...
  sum(Time >= Fossil_herbivores_temporal$Min..Age[Row] & Time <= Fossil_herbivores_temporal$Max..Age[Row]) == 0 })]
Fossil_sites_odd <- Fossil_sites_odd[!is.na(Fossil_sites_odd)] # Remove "NA" values.

Fossil_herbivores_traitenvironment <- data.frame(matrix(vector(), 0, 22)) # Create the baseline dataframe.

library(e1071) # For using the "skewness" and "kurtosis" functions below.

Metrics_long <- function(Trait) { c(ifelse(length(Trait) == 0, NA, mean(Trait, na.rm = TRUE)), ifelse(length(Trait) == 0,
  NA, median(Trait, na.rm = TRUE)), ifelse(length(Trait) <= 1, NA, sd(Trait, na.rm = TRUE)), ifelse(length(Trait) <= 1,
  NA, skewness(Trait, na.rm = TRUE)), ifelse(length(Trait) <= 1, NA, kurtosis(Trait, na.rm = TRUE))) } # Add sk and ku.

for (Age in seq(length(Time))) { # For each time period...
  Fossil_sites_metrics[[Age]] <- list() # Create a separate list of metrics for each time period.
  Fossil_sites_specnum[[Age]] <- list() # Create a separate list of sample sizes for each time period.
  Counter <- 1 # Create a site counter.
  
  for (Site_name in colnames(Fossil_herbivores[which(colnames(Fossil_herbivores) == 'Amboseli_Pleis'):
    ncol(Fossil_herbivores)])) { # For each site within the time period...
    
    # Define a variable, "Logical", that evaluates whether the site occurs during or sufficiently close to the time period:
    
    if (!(Site_name %in% Fossil_sites_odd)) { # If the site is not a member of "Fossil_sites_odd"...
      Logical <- !is.na(Fossil_herbivores_temporal$Min..Age[Counter]) & # Does the site have a temporal range AND...
        Fossil_herbivores_temporal$Min..Age[Counter] <= Time[Age] & # ...is the site's minimum age <= period AND...
        Fossil_herbivores_temporal$Max..Age[Counter] >= Time[Age] } # ...is the site's maximum age >= period?
    else { # If the site is a member of "Fossil_sites_odd"...
      Logical <- Age == which.min(abs(Time - Fossil_herbivores_temporal$Mean.Age[Counter])) } # Closest period to site age.
    
    # Calculate metrics of traits from the species occurring at the site, assuming the site has fulfilled the temporal
    # criteria of "Logical" above. Also calculate the sample sizes of species that went into those trait calculations,
    # keeping in mind that species with missing trait values should not be included:
    
    if (Logical) { # If the criteria of "Logical" have been fulfilled by the site...
      
      # Compute the metrics from the species occurring at the site:
      
      Fossil_sites_metrics[[Age]][[Counter]] <- c(
        Metrics_long(Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores[Site_name] == 1]),
        Metrics_long(Fossil_herbivores$HYP[Fossil_herbivores[Site_name] == 1]),
        Metrics_long(Fossil_herbivores$LOP[Fossil_herbivores[Site_name] == 1]))
      
      # Compute the associated sample sizes for those metric calculations:
      
      Fossil_sites_specnum[[Age]][[Counter]] <- c(
        length(which(!is.na(Fossil_herbivores$Mass_Kg_Final[Fossil_herbivores[Site_name] == 1]))),
        length(which(!is.na(Fossil_herbivores$HYP[Fossil_herbivores[Site_name] == 1]))),
        length(which(!is.na(Fossil_herbivores$LOP[Fossil_herbivores[Site_name] == 1])))) }
    
    else { # If the site's minimum and maximum age put its temporal range outside of the time period...
      Fossil_sites_metrics[[Age]][[Counter]] <- rep(NA, length(Traitmetric_fullnames))
      Fossil_sites_specnum[[Age]][[Counter]] <- rep(NA, length(Trait_names)) }
    
    # Integrate metrics and other associated information (site, age index, MAT, AP, sample sizes) into the dataframe:
    
    Fossil_herbivores_traitenvironment <- rbind(Fossil_herbivores_traitenvironment, data.frame('Site' = Site_name,
      'Age_Index' = Age, matrix(Fossil_sites_metrics[[Age]][[Counter]], nrow = 1, dimnames = list(c(),
      Traitmetric_abbrevnames)), 'MAT' = 0, 'AP' = 0, matrix(Fossil_sites_specnum[[Age]][[Counter]], nrow = 1,
      dimnames = list(c(), c('SampSize_mass', 'SampSize_HYP', 'SampSize_LOP')))))
    
    # Update the site index/counter:
    
    Counter <- Counter + 1 }}; rm(Age, Site_name, Counter, Logical)

# Subset "Fossil_herbivores_traitenvironment" down to only rows with proper trait metric values in them:

Fossil_herbivores_traitenvironment <- Fossil_herbivores_traitenvironment[!is.na(
  Fossil_herbivores_traitenvironment$SampSize_mass), ]

# Integrate woody or grassy cover data from "Fossil_herbivores_cover" into "Fossil_herbivores_traitenvironment":

Fossil_herbivores_traitenvironment$VC <- Fossil_herbivores_cover[match(Fossil_herbivores_traitenvironment$Site,
  Fossil_herbivores_cover$Site), grep('_Mean', colnames(Fossil_herbivores_cover))] # "VC" = vegetation cover.

# Create dataframes of modern data that mirror those created above for fossil data:
  
  # Create a dataframe in which each site in "Modern_herbivores" is represented by the mean, median, standard deviation
  # skewness, and kurtosis of all trait values present there. Also determine the species sample sizes for each site:
  
  Modern_sites_metrics <- list() # List to store metric values.
  Modern_sites_specnum <- list() # List to store sample sizes.
  Modern_herbivores_traitenvironment <- data.frame(matrix(vector(), 0, 21)) # Create the baseline dataframe.
  
  Counter <- 1 # Create a site counter.
    
  Traitmetric_colnames <- colnames(Fossil_herbivores_traitenvironment)[grep('Mass|HYP|LOP', colnames(
    Fossil_herbivores_traitenvironment))][1:length(Traitmetric_fullnames)] # For reference.
  
  for (Site in colnames(Modern_herbivores[which(colnames(Modern_herbivores) == 'ABER'):ncol(Modern_herbivores)])) {
      
    # Calculate metrics for each site:
      
    Modern_sites_metrics[[Counter]] <- c(Metrics_long(Modern_herbivores$Mass_Kg_Final[Modern_herbivores[Site] == 1]),
      Metrics_long(Modern_herbivores$HYP[Modern_herbivores[Site] == 1]),
      Metrics_long(Modern_herbivores$LOP[Modern_herbivores[Site] == 1]))
    
    # Calculate the species sample sizes for each site:
    
    Modern_sites_specnum[[Counter]] <- c(
      length(which(!is.na(Modern_herbivores$Mass_Kg_Final[Modern_herbivores[Site] == 1]))),
      length(which(!is.na(Modern_herbivores$HYP[Modern_herbivores[Site] == 1]))),
      length(which(!is.na(Modern_herbivores$LOP[Modern_herbivores[Site] == 1]))))
    
    # Integrate the metrics and other associated information (site, MAT, AP) in the dataframe:
    
    Modern_herbivores_traitenvironment <- rbind(Modern_herbivores_traitenvironment, data.frame('Site' = Site, matrix(
      Modern_sites_metrics[[Counter]], nrow = 1, dimnames = list(c(), Traitmetric_colnames)), 'MAT' = 0, 'AP' = 0,
      matrix(Modern_sites_specnum[[Counter]], nrow = 1, dimnames = list(c(), c('SampSize_mass', 'SampSize_HYP',
      'SampSize_LOP')))))
    
    # Update the site index/counter:
    
    Counter <- Counter + 1 }; rm(Site, Counter)
  
# Integrate woody or grassy cover data from "Modern_herbivores_cover" into "Modern_herbivores_traitenvironment":
  
Modern_herbivores_traitenvironment$VC <- NA # Initiate blank column, because not all sites have cover values.
Modern_herbivores_traitenvironment$VC[match(Modern_herbivores_cover$Code, Modern_herbivores_traitenvironment$Site)] <-
  Modern_herbivores_cover[, grepl(Cover_type, colnames(Modern_herbivores_cover))]

# Subset "Modern_herbivores_traitenvironment" down to only rows with proper trait metric values in them, and/or with
# vegetation cover data. Perform this subsetting differently depending on if the ecometric analyses that follow use MAT
# and AP, or alternatively vegetation cover, as the number of sites with each of these values may vary:

Maxlike_env <- 'VegCover' # "Climate" or "VegCover" - which will be used for ecometric analyses?
  
if (Maxlike_env == 'Climate') { # If MAT and AP will be used...
  Modern_herbivores_traitenvironment <- Modern_herbivores_traitenvironment[complete.cases(
    Modern_herbivores_traitenvironment[, match(Traitmetric_abbrevnames, colnames(Modern_herbivores_traitenvironment))]), ]
  } else { # If woody or grassy cover will be used...
  Modern_herbivores_traitenvironment <- Modern_herbivores_traitenvironment[which(!(is.na(
    Modern_herbivores_traitenvironment$VC)) & Modern_herbivores_traitenvironment$SampSize_mass > 1), ] }
          
###########################################################################################################################
###########################################################################################################################
  
# PART IX: MAXIMUM LIKELIHOOD ANALYSIS OF HERBIVORE ECOMETRIC RELATIONSHIPS THROUGH TIME 
    
# Perform analyses on both fossil and modern data together to determine how the relationship between traits and
# environmental variables has changed through time. Start by setting up various parameters, options, and variables:

library(animation) # For creating animations of plots below.

Color_number <- 1000 # Number of colors to include in plots below.
Traitmetric_choices <- 'mean|sd' # Or 'sk', 'ku' - determines which metrics will be used in ecometric analyses below.
setwd('../../Figures/V19') # Set the working directory to where visualizations will be stored.
ani.record(reset = TRUE) # Reset/clear animation memory.

Overall_herbivores_traitenvironment <- rbind(Fossil_herbivores_traitenvironment[, -grep('Age', colnames(
  Fossil_herbivores_traitenvironment))], Modern_herbivores_traitenvironment) # Combine "traitenvironment" datasets.

Max <- max(Fossil_herbivores_cover$Max_Age[grepl(Fossil_sites_term, Fossil_herbivores_cover$Derived_From)], na.rm = TRUE)
Min <- min(Fossil_herbivores_cover$Min_Age[grepl(Fossil_sites_term, Fossil_herbivores_cover$Derived_From)], na.rm = TRUE)
Fossil_ecometric_agebreaks <- c(Max, Events[3], mean(c(Events[4], Events[5])), mean(c(Events[8], Events[9])),
  Events[11], Min) # Time point breaks.

Fossil_ecometric_ageranges <- list() # List to hold ranges of time, based upon those breaks.
for (Break in 1:(length(Fossil_ecometric_agebreaks) - 1)) { # For each break, minus the last one...
  Fossil_ecometric_ageranges[[Break]] <- c(which.min(abs(Time - Fossil_ecometric_agebreaks[Break])), which.min(abs(Time -
    Fossil_ecometric_agebreaks[Break+1]))) }; rm(Break) # Turn the break into a range of time between it and the next break.

# Check for sampling bias in "Overall_herbivores_traitenvironment" by determining if a relationship exists between trait
# metric values and the sample sizes of species (at sites) from which those values were derived. If such a relationship
# exists for a given metric, consider converting metric values to the residuals of that relationship:

  # Pre-set conditions for this analysis:

  par(mar = c(9,7,1,1), mfrow = c(3,2)) # Set plot panel settings.
  Tmp <- Overall_herbivores_traitenvironment[!duplicated(Overall_herbivores_traitenvironment$Site), ] # Remove dup sites.
  Traitmetric_rsq <- list() # List to hold R-squared values of metric vs. sample size relationships.
  Traitmetric_resids <- list() # List to hold residuals of metric vs. sample size relationships.
  Traitmetric_method <- 'Raw' # "Residuals" or "Raw", depending upon if residual conversion will be done.
  Counter <- 1 # Counter for the loop below.

  # Check if a relationship exists, and potentially perform the residual conversion of traitmetric values:

  for (Met in grep(Traitmetric_choices, colnames(Tmp))) { # For each trait metric of interest...
    
    # Check if a relationship exists between the metric and sample size, through plotting and correlation analysis:
    
    plot(Tmp$SampSize_mass, Tmp[, Met], xlab = '', ylab = '', cex.axis = 2.5, cex = 3, pch = 19) # Metric vs. sample sizes.
    mtext('Number of Species', side = 1, line = 4, font = 2, cex = 2) # Format X-axis label.
    mtext(Traitmetric_fullnames[Met-1], side = 2, line = 4, font = 2, cex = 1) # Format Y-axis label.
    abline(lm(Tmp[, Met] ~ Tmp$SampSize_mass), lwd = 3, col = 'red') # Regression line of metric vs. sample sizes.
    
    # Perform the residual conversion of the metric, if it is chosen to be done:
    
    Traitmetric_rsq[[Counter]] <- cor(Tmp$SampSize_mass, Tmp[, Met])^2 # R-squared value of metric vs. sample sizes.
    Traitmetric_resids[[Counter]] <- lm(Tmp[, Met] ~ Tmp$SampSize_mass)$residuals # Residuals of metric vs. sample sizes.
    if (Traitmetric_method == 'Residuals') { # If this method is to be done...
      Overall_herbivores_traitenvironment[, Met] <- unname(Traitmetric_resids[[Counter]][match(
        Overall_herbivores_traitenvironment$Site, Tmp$Site)]) } # Apply residuals to "Overall_herbivores_traitenvironment".
    Counter <- Counter + 1 }; rm(Met, Tmp, Counter) # Update counter.
  
  # Apply residual conversion (or not) to "Fossil_herbivores_traitenvironment" and "Modern_herbivores_traitenvironment":
  
  for (Dat in c('Fossil_herbivores_traitenvironment', 'Modern_herbivores_traitenvironment')) { # Per dataset...
    Tmp <- eval(parse(text = Dat)) # Create a temporary variable equal to the dataset.
    Tmp[, colnames(Tmp)[grep(Traitmetric_choices, colnames(Tmp))]] <- Overall_herbivores_traitenvironment[match(Tmp$Site,
      Overall_herbivores_traitenvironment$Site), colnames(Tmp)[grep(Traitmetric_choices, colnames(Tmp))]] # Update values.
    assign(Dat, Tmp) }; rm(Dat, Tmp) # Assign the updated dataset to the original variable name.

# Create a function that builds and plots a raster grid, which will be used for visualizations below. Its inputs are:
# "Mat" = the matrix of values from the maximum likelihood analysis below to be visualized, "Trait" = the trait of
# interest, "Trait_met1max" = the maximum value of the first metric specified by "Traitmetric_choices" for the trait,
# and accordingly, "Trait_met2max", "Trait_met1min", "Trait_met2min"; and "Clim" = MAT, AP, or VC (the code in this fxn
# is based upon code written by Michelle and her lab in the file titled "Africa_Mammal_Ecometrics.R"):
  
RasterGrid <- function(Mat, Trait, Trait_met1max, Trait_met2max, Trait_met1min, Trait_met2min, Clim) {
  
  # Create a raster object, and set the values of that object to be those in "Mat":
  
  Ras <- raster(extent(0, ncol(Mat), 0, nrow(Mat)), resolution = 1)
  Ras <- setValues(Ras, Mat)
  
  # Create an empty grid with appropriate grid lines, within which the values of "Mat" will be plotted below:
  
  print(plot(1, type = 'n', xlim = c(0, ncol(Mat)), ylim = c(0, nrow(Mat)), xaxs = 'i', yaxs = 'i', asp = 1, axes = FALSE,
    xlab = '', ylab = '', main = '')) # Set up a blank plot.
  rect(0, 0, ncol(Mat), nrow(Mat), lwd = 5) # Plot a rectangle.
  lines(x = floor(c(ncol(Mat)/3, ncol(Mat)/3)), y = c(0, nrow(Mat))) # Add a line within the rectangle.
  lines(x = ceiling(c(ncol(Mat)/1.5, ncol(Mat)/1.5)), y = c(0, nrow(Mat))) # Add a line within the rectangle.
  lines(x = c(0, ncol(Mat)), y = floor(c(nrow(Mat)/3, nrow(Mat)/3))) # Add a line within the rectangle.
  lines(x = c(0, ncol(Mat)), y = ceiling(c(nrow(Mat)/1.5, nrow(Mat)/1.5))) # Add a line within the rectangle.
  
  # Plot the values of "Mat" in the empty grid. In so doing, set the colors used to portray the climate values in the
  # grid. Set the breaks of those colors to be based upon the distribution of the climate values in the "traitenvironment"
  # datasets, such that they are scaled to show as much differentiation between values as possible:
  
  Color_number_ras <- Color_number/10 # Change number of colors to display.
  if (Maxlike_env == 'Climate') { Ras_color <- colorRampPalette(c('steelblue4', 'steelblue3', 'yellow2',
    'darkgoldenrod1', 'red3'))(Color_number_ras) } else { Ras_color <- colorRampPalette(c('olivedrab1',
    'olivedrab4'))(Color_number_ras) } # Different colors for different environmental variables.
  
  Ras_column <- ifelse(Maxlike_env == 'VegCover', 'VC', ifelse(grepl('MAT', Clim), 'MAT', 'AP')) # For reference below.
  
  Ras_breaks <- unname(quantile(Overall_herbivores_traitenvironment[, Ras_column], seq(0, 1, length.out =
    Color_number_ras + 1), na.rm = TRUE)) # Color breaks, based upon quantiles of values.
  for (Val in which(duplicated(Ras_breaks))) { Ras_breaks[Val] <- Ras_breaks[Val-1] + 0.0000000000001 }; rm(Val) # Rm dups.

  Ras_range <- c(min(Overall_herbivores_traitenvironment[, Ras_column], na.rm = TRUE),
    max(Overall_herbivores_traitenvironment[, Ras_column], na.rm = TRUE)) # Range of values.

  print(plot(Ras, col = Ras_color, breaks = Ras_breaks, cex = 3, add = TRUE, legend = FALSE, main = '')) # Plot raster.
    
  # Add a legend/colorbar to the grid, using settings consistent with the colors in the grid itself:
    
  image.plot(legend.only = TRUE, zlim = Ras_range, col = Ras_color, breaks = Ras_breaks, axis.args = list(cex.axis = 1.5))
      
  # Format the tick marks and labels of both axes, as well as the label of the legend/colorbar, in the grid:
  
  axis(1, at = round(seq(0, ncol(Mat), length.out = ncol(Mat)/(ncol(Mat)/4))), labels = format(round(seq(Trait_met1min,
    Trait_met1max, length.out = ncol(Mat)/(ncol(Mat)/4)), 1), nsmall = 1), las = 1, cex.axis = 2)
  axis(2, at = round(seq(0, nrow(Mat), length.out = nrow(Mat)/(nrow(Mat)/4))), labels = format(round(seq(Trait_met2min,
    Trait_met2max, length.out = nrow(Mat)/(nrow(Mat)/4)), 1), nsmall = 1), las = 1, cex.axis = 2, pos = 0)
  title(ylab = paste(capitalize(strsplit(Traitmetric_choices, '\\|')[[1]][2]), Trait), line = 4, cex.lab = 2, font.lab = 1)
  title(xlab = paste(capitalize(strsplit(Traitmetric_choices, '\\|')[[1]][1]), Trait), line = 4, cex.lab = 2, font.lab = 1)
  mtext(Clim, side = 4, line = 2, cex = 1.5, font = 1) }

# Initiate an analysis of the maximum likelihoods of climate, as predicted by the "Traitmetric_choices" of sites in the
# "traitenvironment" datasets. Begin the analysis by pre-setting parameters that will be applied throughout:

library(grDevices) # For using the "nclass" function below.
library(Hmisc) # For using the "capitalize" function below.

Maxlike_sampsize_site <- 2 # Minimum number of sites that must be included in each maximum likelihood calculation.
Maxlike_sampsize_spec <- 3 # Minimum number of species that must occur at a site for the site to be considered.

Maxlike_fossilsiteages <- aggregate(Fossil_herbivores_traitenvironment[, 'Age_Index'], by =
  Fossil_herbivores_traitenvironment['Site'], FUN = mean) # Determine the mean age index of each fossil site.
Maxlike_fossilsiteages$x <- unlist(apply( # Determine per row which "Fossil_ecometrics_agerange" the age index falls into.
  sapply(Fossil_ecometric_ageranges, function(Age) { # Per age range...
    Maxlike_fossilsiteages$x <= Age[2] & Maxlike_fossilsiteages$x >= Age[1] }), # Determine if the age index falls into it.
  1, function(Row) { which(Row == TRUE)[1] })) # Convert that determination from a logical into the index of the age range.

# Set the method by which maximum likelihood climates will be calculated and analyzed below. The methods are as follows:
# "All" = all sites across all ages are used to calculate the maximum likelihood climates, and then anomalies from those
# calculations are determined per age range in "Fossil_ecometric_ageranges"; "ByAge" = only the sites in each age range
# are used to calculate maximum likelihood climates and associated anomalies; "Project" = the sites in a given age range
# are used to calculate maximum likelihood climates, and then the anomalies of sites from the next age range are
# determined from those calculations, such that the given age's calculations are projected to the next one's sites; OR
# "Modern" = the sites in the modern age are used to calculate maximum likehlihood climates, and then the anomalies from
# those calculations are determined per age range in "Fossil_ecometric_ageranges"

Maxlike_method <- 'All'
Maxlike_method_length <- ifelse(grepl('All', Maxlike_method), 1, length(Fossil_ecometric_ageranges) + 1)

# Set up the single lists, vector, and dataframe that will be used and processed in the "Age" "for" loop below:

Maxlike_agerange_dfs <- list() # List to hold processed subsets of the "traitenvironment" datasets by age range.
Maxlike_agerange_siteclimsd <- list() # List to hold SDs of climate for sites duplicated in "Maxlike_agerange_dfs".
Maxlike_cols <- list() # List to hold the trait metric columns of the "traitenvironment" datasets of interest.
Maxlike_breakpointn <- list() # List to hold the number of bins in which to place trait metric values per column.
Maxlike_bincodes <- list() # List to hold bin ID integers for each value in each column.

Maxlike_sites_which <- Fossil_sites_choose # "All" or "Paleosol" - include all sites or only those with paleosol data?
if (Maxlike_sites_which == 'Paleosol') { # If only sites with paleosol soil carbonate data will be analyzed...
  Maxlike_sites_remove <- Fossil_herbivores_cover$Site[grepl('Literature', Fossil_herbivores_cover$Derived_From)] } else {
  Maxlike_sites_remove <- c() } # Designate the sites to remove, or otherwise indicate that no sites should be removed. 

Maxlike_evaluate_colnames <- c('Site', 'MAT', 'AP', 'VC', 'Maxlike_BodyMass', 'Maxlike_HYP', 'Maxlike_LOP')
Maxlike_evaluate_colnames_actual <- c('MAT', 'AP', 'VC') # Columns referring to actual environment values.
Maxlike_evaluate_colnames_pred <- Maxlike_evaluate_colnames[grepl('Maxlike', Maxlike_evaluate_colnames)] # Predict columns.
Maxlike_evaluate <- data.frame(matrix(ncol = length(Maxlike_evaluate_colnames), nrow = 0, dimnames = list(c(),
  Maxlike_evaluate_colnames))) # Dataframe to hold actual versus trait-predicted values of climate.

# Set up the nested lists that will be used and processed in the "Age" "for" loop:

Maxlike_values <- list() # List to hold the matrices with maximum likelihood MAT, AP, or VC values.
Maxlike_values_sampsize <- list() # List to hold the associated sample sizes of sites in those matrices.
for (I in 1:Maxlike_method_length) { Maxlike_values[[I]] <- list(); Maxlike_values_sampsize[[I]] <- list() }; rm(I) # Nest.

# Set up parameters for creating a GIF from all plots made below:

saveGIF({ # Save all plots created below as a GIF animation.
ani.options(interval = 0.8) # Set the time interval between plot frames.

# Begin iterating through all age ranges (or the one age range, depending upon the "Maxlike_method") for which maximum
# likelihood climate values will be calculated, and perform operations to make those calculations:
  
for (Age in 1:Maxlike_method_length) { # For the number of age ranges for which to calculate maximum likelihood climates...
  
  # Set "Maxlike_agerange_dfs" to the portion(s) of the "traitenvironment" datasets of interest. If maximum likelihood
  # climates will be determined separately per age, set "Maxlike_agerange_dfs" to the subset of the "traitenvironment"
  # datasets that encompasses the sites occurring at the current age. Alternatively, set it to all rows in the datasets:

  if (!(grepl('All', Maxlike_method))) { # If maximum likelihood climates will be determined separately per age...
    if (Age != length(Fossil_ecometric_ageranges) + 1) { # If the age in question is not modern...
      Maxlike_agerange_dfs[[Age]] <- Fossil_herbivores_traitenvironment[Fossil_herbivores_traitenvironment$Age_Index <=
        Fossil_ecometric_ageranges[[Age]][2] & Fossil_herbivores_traitenvironment$Age_Index >=
        Fossil_ecometric_ageranges[[Age]][1], ] } # Subset to the sites occurring at the current age range.
    else { Maxlike_agerange_dfs[[Age]] <- Modern_herbivores_traitenvironment }} # If the age in question is modern.
      
  if (grepl('All', Maxlike_method)) { # If maximum likelihood climates will be determined from all ages at once...
    Maxlike_agerange_dfs[[Age]] <- Overall_herbivores_traitenvironment }
      
  # Based upon downstream analyses, remove specific sites from "Maxlike_agerange_dfs" which exhibit more uncertainty:
  
  if (Maxlike_sites_which != 'All') { # If not all sites will be analyzed...
    Maxlike_agerange_dfs[[Age]] <- Maxlike_agerange_dfs[[Age]][-which(Maxlike_agerange_dfs[[Age]]$Site %in%
      Maxlike_sites_remove), ] } # Remove the sites that were designated to remove above.

  # Before processing "Maxlike_agerange_dfs" further, calculate the standard deviations of MAT and AP for sites that
  # appear more than once within it, to check if climate is relatively stable across the age range at those sites:
  
  Maxlike_agerange_siteclimsd[[Age]] <- aggregate(Maxlike_agerange_dfs[[Age]][, c('MAT', 'AP')], by = 
    Maxlike_agerange_dfs[[Age]]['Site'], FUN = sd) # Calculate the standard deviations for all sites.
  Maxlike_agerange_siteclimsd[[Age]] <- Maxlike_agerange_siteclimsd[[Age]][complete.cases(
    Maxlike_agerange_siteclimsd[[Age]]), ] # Keep only the standard deviations for sites appearing more than once.
  
  # Aggregate the duplicate sites in "Maxlike_agerange_dfs" by taking the means of their values. In certain cases, 
  # update the "Maxlike_evaluate" dataframe with some of the resulting information:
  
  Maxlike_agerange_dfs[[Age]] <- aggregate(Maxlike_agerange_dfs[[Age]][, -which(colnames(
    Maxlike_agerange_dfs[[Age]]) == 'Site')], by = Maxlike_agerange_dfs[[Age]]['Site'], FUN = mean)
  
  if (grepl('All', Maxlike_method) | (grepl('Modern', Maxlike_method) & Age == Maxlike_method_length)) { # If true...
    Cols <- c('Site', Maxlike_evaluate_colnames_actual) # Columns of "Maxlike_evaluate" in which to fill values below.
    Values <- Maxlike_agerange_dfs[[Age]][, Cols] # Values to fill into those columns.
    Maxlike_evaluate <- rbind(Maxlike_evaluate, data.frame(matrix(ncol = ncol(Maxlike_evaluate), nrow = nrow(Values),
      dimnames = list(c(), Maxlike_evaluate_colnames)))) # Insert new rows in which values will be placed.
    Maxlike_evaluate[, Cols] <- Values } # Insert values into new rows.
  
  # Per column in "Maxlike_agerange_dfs" that refers to a metric depicted by "Traitmetric_choices", perform analyses to
  # determine, visualize, and employ the optimal number of bins that its values should be placed into:

    # Determine the optimal number of bins per column using a statistical test:
    
    Maxlike_cols[[Age]] <- grep(c(Traitmetric_choices), colnames(Maxlike_agerange_dfs[[Age]])) # Col numbers of interest.
    Maxlike_breakpointn[[Age]] <- sapply(Maxlike_cols[[Age]], function(Col) { # For each column...
      nclass.scott(Maxlike_agerange_dfs[[Age]][, Col]) }) # Statistically determine that column's optimal number of bins.
  
    # Visualize those bin numbers by plotting a histogram per column:
  
    par(mfrow = c(3,2), mar = c(7,7,1,1)) # Set up subplot layout.
    sapply(Maxlike_cols[[Age]], function(Col) { # For each column...
      hist(Maxlike_agerange_dfs[[Age]][, Col], breaks = Maxlike_breakpointn[[Age]][match(Col, Maxlike_cols[[Age]])],
        xlab = '', ylab = '', main = '', col = 'deepskyblue4', cex.axis = 2) # Use "Maxlike_breakpointn" to set hist bins.
      mtext(side = 1, line = 4, Traitmetric_fullnames[Col - 1], font = 2, cex = 2) # Format X-axis label.
      mtext(side = 2, line = 4, 'Frequency', font = 2, cex = 2) }) # Format Y-axis label.
  
    # Employ those bin numbers by binning each column's values into the number of bins specified in "Maxlike_breakpointn".
    # Encode the binning by assigning an integer to each value in the column, where the integer is a number between 1 and
    # the number of bins, and where the integer designates the identity of the bin in which the value falls:
  
    Maxlike_bincodes[[Age]] <- sapply(Maxlike_cols[[Age]], function (C) { # For each column...
      cut(Maxlike_agerange_dfs[[Age]][, C], Maxlike_breakpointn[[Age]][match(C, Maxlike_cols[[Age]])], labels = FALSE) })
      
  # Concatenate "Maxlike_bincodes" to their respective rows in "Maxlike_agerange_dfs":
  
  Maxlike_agerange_dfs[[Age]] <- cbind(Maxlike_agerange_dfs[[Age]], as.data.frame(Maxlike_bincodes[[Age]]))
  
  # Calculate the maximum likelihood of MAT for each binned combination of the "Traitmetric_choices" metrics for body
  # mass, and correspondingly, AP for HYP and LOP. Plot the resultant calculated values:
  
  par(mfrow = c(2,2), mar = c(8,8,8,8.2)) # Initiate raster grid layout and margin settings.
  
  for (Trait in seq(length(Trait_names))) { # For each trait...
    
    # Create a counter for denoting columns of "Maxlike_bincodes" and values of "Maxlike_breakpointn" to refer to below.
    # Also denote which dataframe column should be referred to below:
    
    Counter <- ifelse(Trait == 1, 1, ifelse(Trait == 2, 3, 5))
    Maxlike_clim_col <- ifelse(Maxlike_env == 'VegCover', 'VC', ifelse(Trait == 1, 'MAT', 'AP'))
  
    # Create an empty matrix to hold maximum likelihood values for binned combinations. Create a separate matrix to
    # hold the sample sizes of sites that went into determining those values:
    
    Maxlike_values[[Age]][[Trait]] <- array(NA, dim = c(Maxlike_breakpointn[[Age]][Counter+1],
      Maxlike_breakpointn[[Age]][Counter])) # Set the dimensions of the matrix based upon "Maxlike_breakpointn".
    Maxlike_values_sampsize[[Age]][[Trait]] <- array(NA, dim = c(Maxlike_breakpointn[[Age]][Counter+1],
      Maxlike_breakpointn[[Age]][Counter])) # Have this matrix mirror the "Maxlike_values" matrix.
    
    for (Col in 1:Maxlike_breakpointn[[Age]][Counter]) { # For each matrix column, referring to one trait metric...
      for (Row in 1:Maxlike_breakpointn[[Age]][Counter+1]) { # For each matrix row, referring to the second metric...
    
        # Select the MAT,AP, or VC values that pertain to the sites that are binned in the present column/row combination.
        # Also identify the sites and ages (the index of the age range in "Fossil_ecometric_ageranges" or a separate value
        # indicating the modern age) pertaining to each such value. Store all of this information in a data frame:
        
        Maxlike_clim <- Maxlike_agerange_dfs[[Age]][which( # Select the dataframe from which to get climate values.
          Maxlike_bincodes[[Age]][, Counter + 1] == Row & # Select the bin code referring to the row-based trait metric.
          Maxlike_bincodes[[Age]][, Counter] == Col & # Select the bin code referring to the column-based trait metric.
          Maxlike_agerange_dfs[[Age]]$SampSize_mass >= Maxlike_sampsize_spec & # Sample size constraint.
          Maxlike_agerange_dfs[[Age]]$SampSize_HYP >= Maxlike_sampsize_spec & # Sample size constraint.
          Maxlike_agerange_dfs[[Age]]$SampSize_LOP >= Maxlike_sampsize_spec), # Sample size constraint.
          which(colnames(Maxlike_agerange_dfs[[Age]]) %in% c('Site', Maxlike_clim_col))] # Select columns.
        
        if (nrow(Maxlike_clim) >= Maxlike_sampsize_site) { # If "Maxlike_clim" has enough MAT, AP, or VC values...
          Maxlike_clim$AgeRange_Index <- Maxlike_fossilsiteages$x[match(Maxlike_clim$Site, Maxlike_fossilsiteages$Site)]
          Maxlike_clim$AgeRange_Index[is.na(Maxlike_clim$AgeRange_Index)] <- length(Fossil_ecometric_ageranges) + 1 }
        
        # Make a density plot from the MAT, AP, or VC values in "Maxlike_clim", and calculate their maximum likelihood
        # values from that plot. To ensure that each age range in "Fossil_ecometric_ageranges" contributes equally to the
        # maximum likelihood calculation, bootstrap resample a set number of climate values from "Maxlike_clim" "Reps"
        # times, perform the calculation each time, and take the mean of those calculations. Ensure that the probabilities
        # of resampling a given value from a given age range are spread out evenly across the age ranges:
        
        if (nrow(Maxlike_clim) >= Maxlike_sampsize_site) { # If "Maxlike_clim" has enough MAT, AP, or VC values...
          Maxlike_clim$Resample_ProbPerAge <- 1 / length(unique(Maxlike_clim$AgeRange_Index)) # Prob given to each age.
          Maxlike_clim$Resample_NumPerAge <- vapply(1:nrow(Maxlike_clim), function(Row) { # Per row in "Maxlike_clim"...
            sum(Maxlike_clim$AgeRange_Index == Maxlike_clim$AgeRange_Index[Row]) }, numeric(1)) # Number of sites per age.
          Maxlike_clim$Resample_Probs <- Maxlike_clim$Resample_ProbPerAge / Maxlike_clim$Resample_NumPerAge } # P per site.
        
        if (nrow(Maxlike_clim) < Maxlike_sampsize_site) { Maxlike <- NaN } else { # If not enough values, NaN, or else...
          Maxlike <- mean(replicate(Reps, { # Take the mean of "Reps" replications of the process below.
            Sample <- Maxlike_clim[sample(1:nrow(Maxlike_clim), Maxlike_sampsize_site, replace = TRUE, 
              prob = Maxlike_clim$Resample_Probs), Maxlike_clim_col] # Choose sample.
            Den <- density(Sample, bw = 1) # Create a density plot from that sample of values.
            return(Den$x[which.max(Den$y)]) }), na.rm = TRUE) } # Calculate the maximum likelihood from that density plot.
        
        # Insert the maximum likelihood MAT, AP, or VC value of the row/column combination into the matrix, and record the
        # sample size of sites used to obtain that value in a separate matrix. Further, insert the maximum likelihood
        # value into the appropriate rows (based on site identities) of "Maxlike_evaluate":
          
        Maxlike_values[[Age]][[Trait]][(Maxlike_breakpointn[[Age]][Counter+1] + 1) - Row, Col] <- Maxlike
        Maxlike_values_sampsize[[Age]][[Trait]][(Maxlike_breakpointn[[Age]][Counter+1]+1)-Row, Col] <- nrow(Maxlike_clim)
        
        if (exists('Values')) { # If "Maxlike_evaluate" was set up priorly for this iteration of the loop...
          Maxlike_evaluate[match(Maxlike_clim$Site, Maxlike_evaluate$Site), Maxlike_evaluate_colnames_pred[Trait]] <-
            Maxlike }}}
      
    # Plot a grid of the values in "Maxlike_values", using the "RasterGrid" function made above:
    
    RasterGrid( # Initiate the function created above.
      Maxlike_values[[Age]][[Trait]], # Choose the raster grid to act upon.
      Trait_names[[Trait]], # Specify the trait of interest.
      max(Overall_herbivores_traitenvironment[, grep(c(Traitmetric_choices), colnames(
        Overall_herbivores_traitenvironment))[Counter]], na.rm = TRUE), # "Traitmetric_choice" metric 1 max.
      max(Overall_herbivores_traitenvironment[, grep(c(Traitmetric_choices), colnames(
        Overall_herbivores_traitenvironment))[Counter+1]], na.rm = TRUE), # "Traitmetric_choice" metric 2 max.
      min(Overall_herbivores_traitenvironment[, grep(c(Traitmetric_choices), colnames(
        Overall_herbivores_traitenvironment))[Counter]], na.rm = TRUE), # "Traitmetric_choice" metric 1 min.
      min(Overall_herbivores_traitenvironment[, grep(c(Traitmetric_choices), colnames(
        Overall_herbivores_traitenvironment))[Counter+1]], na.rm = TRUE), # "Traitmetric_choice" metric 2 min.
      ifelse(Maxlike_env == 'VegCover', paste(Cover_type, 'Cover (Fraction)'), ifelse(Trait == 1, 'MAT (degrees C)',
        'AP (mm/year)'))) } # Label.
  
  # Add in text to each plot panel indicating the age that the panel refers to:
  
  if (!(grepl('All', Maxlike_method))) { # If maximum likelihood climates were determined separately per age...
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') # Create an empty plot for text.
    text(0.5, 0.5, ifelse(Age == length(Fossil_ecometric_ageranges) + 1, '0.00 Mya', # If age is modern, add this text.
      paste(as.character(formatC(Time[Fossil_ecometric_ageranges[[Age]][1]], digits = 3, format = 'f')), '-\n',
      as.character(formatC(Time[Fossil_ecometric_ageranges[[Age]][2]], digits = 3, format = 'f')), 'Mya')), cex = 4,
      font = 2) } # Or else, add text identifying the current age range.

  # Add a progress monitor for this process:

  print(paste('Frame number', as.character(Age), 'complete')) }}
  
  # Set the name of the GIF file being saved (change this name depending upon the output settings determined above),
  # as well as its dimensions:
  
  , movie.name = paste('Ecometrics_Animation_', capitalize(gsub('\\|', '', Traitmetric_choices)), '_Method=', 
  Maxlike_method, '.gif', sep = ''), ani.height = 800, ani.width = 800)
rm(Age, Trait, Row, Col, Counter, Maxlike_clim, Maxlike)

# Per trait of interest, evaluate the degree to which its above ecometric model predicts actual climates:

if (Maxlike_method == 'All') { # Constraint this evaluation to this method (for now, at least).
  par(mar = c(8,5,1,1), mfrow = c(1,1)) # Adjust plot layout settings.
  Maxlike_evaluate_data <- sapply(Maxlike_evaluate_colnames_pred, function(Trait) { # Per trait...
    Column <- ifelse(Maxlike_env == 'VegCover', 'VC', ifelse(grepl('Mass', Trait), 'MAT', 'AP')) # Column to refer to.
    Predicted <- Maxlike_evaluate[, Trait]; Actual <- Maxlike_evaluate[, Column] # Get predicted and actual values.
    hist(Actual - Predicted, xlab = '', ylab = '', main = '', col = 'deepskyblue4', cex.axis = 2.75, breaks = ifelse(
      Maxlike_env == 'VegCover', 10, 20), mgp = c(3,2,0)) # Make histogram.
    mtext(side = 1, ifelse(Maxlike_env == 'VegCover', 'Woody Cover (Fraction)', ifelse(grepl('Mass', Trait),
      'MAT (\u00B0C)', 'AP (mm/year)')), line = 6, cex = 3) # Format x-axis label.
    return(Actual - Predicted) }) } # Return the data underlying the histogram.

# Per trait of interest, calculate the percentage of sites in the above analysis that fell into ecometric bins with a
# number of sites greater than or equal to "Maxlike_sampsize_site" (i.e., were included in the analysis):

if (grepl('All', Maxlike_method)) { # Constrain this calculation to this method (for now, at least).
  Maxlike_sitepcts <- sapply(Maxlike_values_sampsize[[1]], function(Grid) { # Per ecometric grid, one per trait...
    return((1 - (sum(Grid[which(Grid < Maxlike_sampsize_site & Grid != 0)]) / sum(Grid))) * 100) }) }

# Generate a backup of "Maxlike_bincodes", as it may be modified below:

Maxlike_bincodes_backup <- Maxlike_bincodes

# Initiate an analysis of the anomalies between the maximum likelihood estimates of climate at each site and each site's
# actual climate. Begin by creating lists to hold anomaly values, the site index which each of those values is associated
# with, and the bootstrapped means of those values for each trait in each age:

Maxlike_anomalies <- list(); Maxlike_anomalies_sites <- list(); Maxlike_anomalies_means <- list()

# Determine the sample size to use for bootstraps below, based upon which age has the least number of sites:

Maxlike_anomalies_sampsize <- min(table(Maxlike_fossilsiteages$x), na.rm = T)

# For each age range, calculate the difference between maximum likelihood values of climate and actual values of climate
# (anomalies) for each site based upon each trait, and quantify the deviation of those calculations from 0:

for (Age in 1:(length(Fossil_ecometric_ageranges) + 1)) { # For each age range, including modern times...
  
  # If the "Project" method of "Maxlike_method" is being used, then the first age range is not relevant for anomalies,
  # since anomalies will be calculated per age range based upon the prior's maximum likelihood climates. Skip over it:
  
  if (grepl('Project', Maxlike_method) & Age == 1) { next }
  
  # Initiate an index for determining which list elements for different lists to refer to below, depending upon
  # "Maxlike_method". After all, if the "All" method is used, then lists only have one element, and if the "Project"
  # method is used, then the first element of lists has been skipped over above:
  
  Idx <- ifelse(grepl('Age', Maxlike_method), Age, ifelse(grepl('Project', Maxlike_method), Age - 1, ifelse(grepl(
    'Modern', Maxlike_method), length(Maxlike_agerange_dfs), 1)))
  
  # If the "Project" method of "Maxlike_method" is being used, modify the values in "Maxlike_bincodes" such that those
  # values reflect the binning structure used for the prior age range, as the prior age is projecting upon this one:
    
  if (grepl('Project', Maxlike_method)) { # If anomalies will be calculated based upon a prior age...
    Maxlike_bincodes[[Age]] <- sapply(Maxlike_cols[[Age]], function (C) { # For each df column used for bin coding...
      
      # Get the break point values that were used to delineate between bin codes for the prior age:
      
      Col <- ifelse(Age == length(Fossil_ecometric_ageranges) + 1, C + 1, C) # Distinguish bw fossil and modern columns.
      Tmp <- cut(Maxlike_agerange_dfs[[Idx]][, Col], Maxlike_breakpointn[[Idx]][match(Col, Maxlike_cols[[Idx]])])
      Tmp <- gsub('.*,', '', levels(Tmp)); Tmp <- as.numeric(gsub(']', '', Tmp))[-length(Tmp)] # Convert break to numeric.
      Tmp <- c(-Inf, Tmp, Inf) # Add "Inf" buffers to either side of the break point values.
      
      # Use those break point values to create new bin codes for the current column in the current age:
      
      return(cut(Maxlike_agerange_dfs[[Age]][, C], Tmp, labels = FALSE)) }) }
  
  # If the "Modern" method of "Maxlike_method" is being used, modify the values in "Maxlike_bincodes" such that those
  # values reflect the binning structure used for the modern age, as the modern age is projecting upon all age ranges:
  
  if (grepl('Modern', Maxlike_method)) { # If anomalies will be calculated based upon the modern age...
    Maxlike_bincodes[[Age]] <- sapply(Maxlike_cols[[Age]], function (C) { # For each df column used for bin coding...
      Col <- ifelse(Age == length(Fossil_ecometric_ageranges) + 1, C, C - 1) # Distinguish bw fossil and modern columns.
      Tmp <- cut(Maxlike_agerange_dfs[[length(Maxlike_agerange_dfs)]][, Col], Maxlike_breakpointn[[length(
        Maxlike_breakpointn)]][match(C, Maxlike_cols[[Age]])]) # Get the breakpoints of the column from the modern age.
      Tmp <- gsub('.*,', '', levels(Tmp)); Tmp <- as.numeric(gsub(']', '', Tmp))[-length(Tmp)] # Convert break to numeric.
      Tmp <- c(-Inf, Tmp, Inf) # Add "Inf" buffers to either side of the break point values.
      return(cut(Maxlike_agerange_dfs[[Age]][, C], Tmp, labels = FALSE)) }) } # Return updated bin codes.
  
  # Create sublists of the above lists to accomodate each trait, differentiating this based upon "Maxlike_method":
  
  if (grepl('Project', Maxlike_method)) { # If anomalies will be calculated based upon a prior age...
    Maxlike_anomalies[[Idx]] <- list(); Maxlike_anomalies_sites[[Idx]] <- list();
    Maxlike_anomalies_means[[Idx]] <- list() } # Use "Idx" for indexing, or else, use "Age"...
  
  Maxlike_anomalies[[Age]] <- list(); Maxlike_anomalies_sites[[Age]] <- list(); Maxlike_anomalies_means[[Age]] <- list()

  for (Trait in seq(length(Trait_names))) { # For each trait...
    
    # Initiate a counter for determining the columns of "Maxlike_bincodes" to refer to below:
    
    Counter <- ifelse(Trait == 1, 1, ifelse(Trait == 2, 3, 5))
    
    # Determine the rows of "Maxlike_bincodes" that correspond to a non-NaN value in "Maxlike_values", in order to know
    # which rows are eligible to be used for bootstrap resampling below:
    
    if (!(grepl('Project|Modern', Maxlike_method))) { # If anomalies will be calculated in a non-projected style...
      Rows <- c() # Vector to hold the rows that should be used for bootstrap resampling.
      for (R in 1:nrow(Maxlike_bincodes[[Idx]])) { # For each row in "Maxlike_bincodes"...
        Rows[R] <- !is.nan(Maxlike_values[[Idx]][[Trait]][(Maxlike_breakpointn[[Idx]][Counter+1] + 1) -
          Maxlike_bincodes[[Idx]][R, Counter+1], Maxlike_bincodes[[Idx]][R, Counter]]) }; rm(R) # Is row with non-NaN?
      Rows <- (1:nrow(Maxlike_bincodes[[Idx]]))[Rows] } # Keep only the rows that correspond to non-NaN.
    if (grepl('Project|Modern', Maxlike_method)) { # If anomalies will be calculated based upon a projection...
      Rows <- c() # Vector to hold the rows that should be used for bootstrap resampling.
      for (R in 1:nrow(Maxlike_bincodes[[Age]])) { # For each row in "Maxlike_bincodes"...
        Rows[R] <- !is.nan(Maxlike_values[[Idx]][[Trait]][(Maxlike_breakpointn[[Idx]][Counter+1] + 1) -
          Maxlike_bincodes[[Age]][R, Counter+1], Maxlike_bincodes[[Age]][R, Counter]]) }; rm(R) # Is row with non-NaN?
      Rows <- (1:nrow(Maxlike_bincodes[[Age]]))[Rows] } # Keep only the rows that correspond to non-NaN.
    
    if (grepl('All', Maxlike_method)) { # If maximum likelihood climates were determined using all sites at once...
      if (Age != length(Fossil_ecometric_ageranges) + 1) { # If the age in question is not modern...
        Tmp <- Fossil_herbivores_cover[match(Maxlike_agerange_dfs[[Idx]]$Site[Rows], Fossil_herbivores_cover$Site), ]
        # Rows <- Rows[Tmp$Mean_Age <= Fossil_ecometric_agebreaks[Age] & Tmp$Mean_Age >= Fossil_ecometric_agebreaks[Age+1]]
        Rows <- Rows[Tmp$Max_Age >= Fossil_ecometric_agebreaks[Age+1] & Tmp$Min_Age <= Fossil_ecometric_agebreaks[Age]]
        Rows <- Rows[!is.na(Rows)]; rm(Tmp) } # Subset "Rows" to just those referring to sites occurring within the age.
      if (Age == length(Fossil_ecometric_ageranges) + 1) { # If the age in question is modern...
        Rows <- Rows[Rows >= which(Maxlike_agerange_dfs[[Idx]]$Site == 'ABER')] }} # Subset "Rows" to just modern sites.
    
    # Calculate the difference calculations/anomalies. Implement bootstrap resampling to account for temporal sampling
    # bias. Do this by resampling a subset of sites (the eligible pool being in "Rows" above), calculating the subset's
    # anomalies, and repeating "Reps" times. Assign all anomaly results to "Maxlike_anomalies":

      # If the "Project" method of "Maxlike_method" is being used, update "Age" and/or "Idx" to ensure proper indexing:
    
      if (grepl('Project', Maxlike_method)) { Age <- Age - 1 }
    
      # Choose the samples of sites based upon the "Rows" vector:
    
      Maxlike_anomalies_sites[[Age]][[Trait]] <- replicate(Reps, { sample(Rows, Maxlike_anomalies_sampsize, replace = T) })
    
      # Using the chosen samples, perform the anomaly calculation for each sample and store the results:
    
      Maxlike_anomalies[[Age]][[Trait]] <- vapply(Maxlike_anomalies_sites[[Age]][[Trait]], function(Samp) { # Per site...
      
        # Assign a value to be used for proper indexing below, depending upon "Maxlike_method":
        
        Val <- ifelse(grepl('Project', Maxlike_method), 1, ifelse(grepl('Modern', Maxlike_method), -(Idx-Age), 0))
        
        # Assign "Tmp" to the trait-environment dataframe that will be used for the calculations:
          
        Tmp <- Maxlike_agerange_dfs[[Idx+Val]]
          
        # Use "Tmp" to perform the calculations of anomalies:
          
        if (all(!is.na(Maxlike_bincodes[[Idx+Val]][Samp, c(Counter, Counter + 1)]))) { # If the site has bin codes...
            
          # Subtract the site's maximum likelihood climate value from its actual climate value:
            
          return(Tmp[Samp, which(colnames(Tmp) == ifelse(Maxlike_env == 'VegCover', 'VC', ifelse(
            Trait == 1, 'MAT', 'AP')))] - # Actual climate value.
          
          Maxlike_values[[Idx]][[Trait]][(Maxlike_breakpointn[[Idx]][Counter+1] + 1) - Maxlike_bincodes[[Idx+Val]][Samp,
            Counter+1], Maxlike_bincodes[[Idx+Val]][Samp, Counter]]) } else { # Maximum likelihood climate value.
        
        # If the site does not have bin codes, return an "NaN" as a filler. Also, specify that each anomaly calculation
        # outputted per site should be in numeric form:
              
        return(NaN) }}, numeric(1))
        
      # Restructure "Maxlike_anomalies" to have the same dimensions as "Maxlike_anomalies_sites", in which each column
      # refers to a different set of anomalies calculated for a different set of bootstrap resampled sites:
      
      Maxlike_anomalies[[Age]][[Trait]] <- structure(Maxlike_anomalies[[Age]][[Trait]], dim = dim(
        Maxlike_anomalies_sites[[Age]][[Trait]]))
      
    # Calculate the mean anomaly of each bootstrapped site subset in "Maxlike_anomalies", per trait and age:
    
    Maxlike_anomalies_means[[Age]][[Trait]] <- apply(Maxlike_anomalies[[Age]][[Trait]], 2, function(Col) {
      mean(Col, na.rm = TRUE) })
  
    # If the "Project" method of "Maxlike_method" is being used, update "Age" and/or "Idx" to ensure proper indexing:
    
    if (grepl('Project', Maxlike_method)) { Age <- Age + 1 } }}; rm(Age, Trait, Rows, Counter, Idx)
          
# Analyze and plot the values in "Maxlike_anomalies_means" to visualize how anomalies have changed over time:

Maxlike_anomalies_means <- Maxlike_anomalies_means[lapply(Maxlike_anomalies_means, length) > 0] # Remove empty elements.

# load('EcometricAnomalies_Bins=Uniform.rds') # Load ecometric result of interest (saved from prior code runs).

Maxlike_anomalies_plots_method <- 'Absolute' # Either plot "Absolute" anomaly values, or those "Relative" to another age.
Maxlike_anomalies_plots_method_firstage <- ifelse(grepl('Modern', Maxlike_method), 1, ifelse(
  Maxlike_anomalies_plots_method == 'Relative', 2, 1)) # If the first age range applicable or not?
Maxlike_anomalies_plots <- list(); Plot_df_breakpoints_ecom <- list(); Y_ticks_ecom <- list() # Lists to hold features.

for (Trait in seq(length(Trait_names))) { # For each trait...
  
  # Process the data that will be used for making the boxplot/violin plot, and convert it to a dataframe while doing so:
  
  Tmp <- as.data.frame(sapply(Maxlike_anomalies_means, '[[', Trait)) # Get mean anomaly values for trait across age(s).
  
  if (Maxlike_anomalies_plots_method == 'Relative') { # If the anomaly values to be plotted are relative to another age...
    if (!grepl('Project|Modern', Maxlike_method)) { # If anomaly values were not calculated via a projection...
      Tmp <- sapply(2:ncol(Tmp), function(Col) { Tmp[, Col] - Tmp[, 1] }) } # Difference in anomaly values from first age.
    if (grepl('Modern', Maxlike_method)) { # If anomaly values were calculated via a projection from the modern age...
      Tmp <- sapply(1:(ncol(Tmp)-1), function(Col) { Tmp[, Col] - Tmp[, ncol(Tmp)] }) }} # Difference from modern age.
  
  Tmp <- stack(as.data.frame(Tmp)) # Reshape "Tmp" into two columns, one with values and the other with factor ages.
  # Tmp$values <- Tmp$values - median(Tmp$values)

  # Derive, from "Tmp", the data that will be used for making the line plot:
  
  Tmp2 <- data.frame(values = sapply(unique(Tmp$ind), function(I) { median(Tmp$values[Tmp$ind == I], na.rm = TRUE) }),
    ind = unique(Tmp$ind)) # One column with medians of each boxplot category, and the other with the categories.
  
  # Record the breakpoints that are associated with the "Tmp2" dataframe from breakpoint analysis:
    
  if (exists('X')) { rm(X) }; if (exists('Y')) { rm(Y) } # Remove these variables so they do not interfere with below.
  Tmp2$X <- 1:nrow(Tmp2); Tmp2$Y <- Tmp2$values # Make an index column and a column of values.
  Plot_df_breakpoints_ecom[[Trait]] <- segmented(lm(Y ~ X, data = Tmp2), seg.Z = ~X, control = seg.control(n.boot = 100))
  Plot_df_breakpoints_ecom[[Trait]] <- ifelse(is.null(Plot_df_breakpoints_ecom[[Trait]]$psi), 1, round(
    Plot_df_breakpoints_ecom[[Trait]]$psi[, 2])) # If no breakpoint was estimated, make the breakpoint the first point.
  
  # Determine, from "Tmp", whether the distribution in each boxplot/violin plot differs significantly from a value. Insert
  # a column into "Tmp" indicating which rows refer to distributions that do or do not differ significantly as such:
  
  Val <- ifelse(Maxlike_env == 'Climate' & Maxlike_anomalies_plots_method == 'Absolute', median(sapply(
    Maxlike_anomalies_means, '[[', Trait)[, Maxlike_anomalies_plots_method_firstage], na.rm = TRUE), 0) # Which value?
  
  Tmp3 <- sapply(levels(Tmp$ind), function(Age) { # For each age in "Tmp", each representing a plot distribution...
    
    # If the median of the distribution is >= "Val", return the percentage of distribution points that are < "Val",
    # and vice versa:
    
    if (Tmp2$values[Tmp2$ind == Age] >= Val) { return(sum(Tmp$values[Tmp$ind == Age] < Val) / sum(Tmp$ind == Age)) }
    if (Tmp2$values[Tmp2$ind == Age] <= Val) { return(sum(Tmp$values[Tmp$ind == Age] > Val) / sum(Tmp$ind == Age)) } })
  
  Tmp$sig <- unname(Tmp3[match(Tmp$ind, names(Tmp3))]); Tmp$sig <- as.factor(ifelse(Tmp$sig < 0.05, 'Y', 'N'))
  
  # Insert a column into "Tmp" indicating whether or not each row is associated with the first age range or not:
  
  Tmp$first <- as.factor(ifelse(Tmp$ind == 'V1', 'Y', 'N'))

  # Determine the Y tick labels that will be used in the plot. Make them clean, round numbers that are evenly spaced:

  Y_ticks_ecom[[Trait]] <- seq(min(Tmp$values, na.rm = TRUE), max(Tmp$values, na.rm = TRUE), length.out = 5)[2:4]
  Y_ticks_ecom[[Trait]] <- round_any(Y_ticks_ecom[[Trait]], ifelse(Maxlike_env == 'VegCover', 0.05, ifelse(
    Trait == 1, 0.5, 10)), f = ceiling) # Round labels.
  Y_ticks_ecom[[Trait]] <- c(Y_ticks_ecom[[Trait]][1], mean(Y_ticks_ecom[[Trait]][c(1,3)]), Y_ticks_ecom[[Trait]][3])
  
  # Make the plot, using the above information:
  
  Maxlike_anomalies_plots[[Trait]] <- ggplot(data = Tmp, aes(x = ind, y = values)) + # Select data for plot.
    geom_violin(data = Tmp, aes_string(fill = ifelse(Maxlike_anomalies_plots_method == 'Absolute', 'sig', 'sig')),
      lwd = 1.5) + # Make violin plots, colored by whether referring to the first age OR by if different from 0.
    scale_fill_manual(values = c('deepskyblue3', 'orange2')) + # Determine the specific fill colors of the plots.
    geom_line(data = Tmp2, aes(ind, values, group = 1), color = 'black', lwd = 2) + # Draw line through plot medians.
    geom_boxplot(width = 0.1) + # Establish boxplots on top of the line and violin plots drawn previously.
    geom_hline(yintercept = Val, lty = '11', size = 1, color = 'black') + # Put line at y = "Val".
    # geom_vline(xintercept = ifelse(is.na(Plot_df_breakpoints_ecom[[Trait]]), 1, Plot_df_breakpoints_ecom[[Trait]]),
    #   linetype = 'solid', size = 3, alpha = 1) + # Plot the breakpoint.
    scale_x_discrete(labels = if (Trait == 3) { c(sapply( # First age range is NA if "Relative", as specified below.
      Maxlike_anomalies_plots_method_firstage:(length(Fossil_ecometric_agebreaks)-1), function(Idx) { # Per age range...
      paste(as.character(formatC(Fossil_ecometric_agebreaks[Idx], digits = 2, format = 'f')), '-', as.character(formatC(
      Fossil_ecometric_agebreaks[Idx+1], digits = 2, format = 'f')), sep = '') }), '0.00') } else { rep('', length(
      Maxlike_anomalies_plots_method_firstage:length(Fossil_ecometric_agebreaks))) }) +
      # Format X tick labels as age ranges, based upon those in "Fossil_ecometric_ageranges", and add modern "0.000" age.
    scale_y_continuous(breaks = Y_ticks_ecom[[Trait]], labels = ifelse(Maxlike_env == 'VegCover', function(Lab) { Lab },
      ifelse(Trait == 1, function(Lab) { sprintf('%.1f', Lab) }, function(Lab) { sprintf('%.0f', Lab) }))) + # Y ticks.
    # labs(y = ifelse(Maxlike_env == 'VegCover' & Trait == 2, paste(Cover_type, 'Cover (Fraction)'), ifelse(
    #   Maxlike_env == 'VegCover' & Trait != 2, '', ifelse(Trait==1, 'MAT (\u00B0C)', 'AP (mm/year)')))) + # Y-axis labels.
    theme(panel.background = element_blank(), axis.line = element_line(color = 'black'), axis.title = element_blank(),
      axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 0, color = 'black'), axis.text.y = element_text(size = 20,
      color = 'black'), plot.margin = unit(c(ifelse(Trait == 1, 0, ifelse(Trait == 2, -1, -2)), 0.5, ifelse(Trait == 1,
      -0.5, ifelse(Trait == 2, 0.5, 1)), 0.5), 'lines'), legend.position = 'none') }; rm(Trait, Tmp, Tmp2, Tmp3) # Format.

grid.arrange(plot_grid(Maxlike_anomalies_plots[[1]], Maxlike_anomalies_plots[[2]], Maxlike_anomalies_plots[[
  length(Traitmetric_names)]], ncol = 1, align = 'v')) # Plot.

setwd('../../Code_and_Data/Data') # Reset the working directory.
