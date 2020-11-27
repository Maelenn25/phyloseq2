Dada2 tutorial
================

  - [Préparation de l’environnement](#préparation-de-lenvironnement)
      - [Chargement des librairies Phyloseq et Dada2 afin de pouvoir
        utiliser ses
        packages](#chargement-des-librairies-phyloseq-et-dada2-afin-de-pouvoir-utiliser-ses-packages)
      - [Création d’une variable ‘path’ dans laquelle on met l’objet
        “Miseq\_SOP”. Ensuite on les
        liste](#création-dune-variable-path-dans-laquelle-on-met-lobjet-miseq_sop.-ensuite-on-les-liste)
      - [Création de deux listes.](#création-de-deux-listes.)
  - [Inspect read quality profiles](#inspect-read-quality-profiles)
  - [Learn the Error Rates](#learn-the-error-rates)
  - [Sample Inference](#sample-inference)
  - [Merge paired reads](#merge-paired-reads)
  - [Construct sequence table](#construct-sequence-table)
  - [Remove chimeras](#remove-chimeras)
  - [Track reads through the
    pipeline](#track-reads-through-the-pipeline)
  - [Assign taxonomy](#assign-taxonomy)
  - [Evaluate accuracy](#evaluate-accuracy)
  - [Sauvegarde de l’environnement pour pouvoir jouer le code de
    phyloseq](#sauvegarde-de-lenvironnement-pour-pouvoir-jouer-le-code-de-phyloseq)

# Préparation de l’environnement

## Chargement des librairies Phyloseq et Dada2 afin de pouvoir utiliser ses packages

``` r
library ("phyloseq")
library("dada2")
```

    ## Loading required package: Rcpp

## Création d’une variable ‘path’ dans laquelle on met l’objet “Miseq\_SOP”. Ensuite on les liste

``` r
path <- "~/MiSeq_SOP" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
```

    ##  [1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
    ##  [3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
    ##  [5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
    ##  [7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
    ##  [9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
    ## [11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
    ## [13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
    ## [15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
    ## [17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
    ## [19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
    ## [21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
    ## [23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
    ## [25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
    ## [27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
    ## [29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
    ## [31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
    ## [33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
    ## [35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
    ## [37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
    ## [39] "filtered"                      "HMP_MOCK.v35.fasta"           
    ## [41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
    ## [43] "mouse.dpw.metadata"            "mouse.time.design"            
    ## [45] "stability.batch"               "stability.files"

## Création de deux listes.

Une avec les Forwards et une avec les Reverses qui correspondent
respectivement au Read1 et Read2. La commande avec ‘sample.names’ permet
de mettre toutes les séquences sous un même format pour harmoniser leurs
noms.

``` r
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

# Inspect read quality profiles

Permet d’inspecter la qualité des différents reads et nous permet de
voir à partir de quel nucléotide nous allons avoir une baisse du score
de qualité pour pouvoir par la suite retirer les nucléotides de mauvaise
qualité. Sur les forwards nous allons retirer les 10 derniers
nucléotides avec une commande retrouvée plus bas dans le script.

``` r
plotQualityProfile(fnRs[1:4])
```

![](03_Dada2_tutorial_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
Ici nous faisons la même manipulation mais avec les Reverses qui sont
d’un peu moins bonne qualité que les Forwards, dû au sequençage avec
Illumina.Ici, nous retirerons tous les nucléotides à partir du 160e

``` r
plotQualityProfile(fnRs[1:4])
```

![](03_Dada2_tutorial_files/figure-gfm/unnamed-chunk-5-1.png)<!-- --> \#
Filter and trim Ici nous créeons une nouvelle variable qui recevra tous
les fichiers filtrés, que ce soit pour les Forwards ou les Reverses.
Nous allons en même temps indiquer à la machine que les noms utilisés
pour ranger les séquences dans les dossiers seront les mêmes que nous
avons standardisé plus haut.

``` r
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

Nous allons indiquer à la machine quels paramètres nous allons utiliser
pour filtrer les séquences avant de les ranger. On indique ainsi dans la
première fonction les deux fichiers d’où les séquences viendront (fnFs
et fnRs), puis le fichier où elles seront rangées (filtFs et filtRs
crées juste au dessus). Enuite nous allons indiquer où couper pour les
deux sortes de séquences.

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
```

    ##                               reads.in reads.out
    ## F3D0_S188_L001_R1_001.fastq       7793      7113
    ## F3D1_S189_L001_R1_001.fastq       5869      5299
    ## F3D141_S207_L001_R1_001.fastq     5958      5463
    ## F3D142_S208_L001_R1_001.fastq     3183      2914
    ## F3D143_S209_L001_R1_001.fastq     3178      2941
    ## F3D144_S210_L001_R1_001.fastq     4827      4312

# Learn the Error Rates

Nous allons ici utiliser des lignes de commandes qui vont permettre
d’apprendre à la machine les différents profils d’erreurs générées
lors du séquençage. L’opération est faite sur les deux types de
séquence.

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

Cette ligne nous permet de visualiser les erreurs qu’on a fait apprendre
à la machine

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](03_Dada2_tutorial_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

# Sample Inference

Ici nous créons une autre variable “dadaFs” dans laquelle nous mettons
les fichiers obtenus après avoir filtré et appliqué le profil d’erreur à
nos séquences. Nous allons faire la même chose avec dadaRS.

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1979 unique sequences.
    ## Sample 2 - 5299 reads in 1639 unique sequences.
    ## Sample 3 - 5463 reads in 1477 unique sequences.
    ## Sample 4 - 2914 reads in 904 unique sequences.
    ## Sample 5 - 2941 reads in 939 unique sequences.
    ## Sample 6 - 4312 reads in 1267 unique sequences.
    ## Sample 7 - 6741 reads in 1756 unique sequences.
    ## Sample 8 - 4560 reads in 1438 unique sequences.
    ## Sample 9 - 15637 reads in 3590 unique sequences.
    ## Sample 10 - 11413 reads in 2762 unique sequences.
    ## Sample 11 - 12017 reads in 3021 unique sequences.
    ## Sample 12 - 5032 reads in 1566 unique sequences.
    ## Sample 13 - 18075 reads in 3707 unique sequences.
    ## Sample 14 - 6250 reads in 1479 unique sequences.
    ## Sample 15 - 4052 reads in 1195 unique sequences.
    ## Sample 16 - 7369 reads in 1832 unique sequences.
    ## Sample 17 - 4765 reads in 1183 unique sequences.
    ## Sample 18 - 4871 reads in 1382 unique sequences.
    ## Sample 19 - 6504 reads in 1709 unique sequences.
    ## Sample 20 - 4314 reads in 897 unique sequences.

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1660 unique sequences.
    ## Sample 2 - 5299 reads in 1349 unique sequences.
    ## Sample 3 - 5463 reads in 1335 unique sequences.
    ## Sample 4 - 2914 reads in 853 unique sequences.
    ## Sample 5 - 2941 reads in 880 unique sequences.
    ## Sample 6 - 4312 reads in 1286 unique sequences.
    ## Sample 7 - 6741 reads in 1803 unique sequences.
    ## Sample 8 - 4560 reads in 1265 unique sequences.
    ## Sample 9 - 15637 reads in 3414 unique sequences.
    ## Sample 10 - 11413 reads in 2522 unique sequences.
    ## Sample 11 - 12017 reads in 2771 unique sequences.
    ## Sample 12 - 5032 reads in 1415 unique sequences.
    ## Sample 13 - 18075 reads in 3290 unique sequences.
    ## Sample 14 - 6250 reads in 1390 unique sequences.
    ## Sample 15 - 4052 reads in 1134 unique sequences.
    ## Sample 16 - 7369 reads in 1635 unique sequences.
    ## Sample 17 - 4765 reads in 1084 unique sequences.
    ## Sample 18 - 4871 reads in 1161 unique sequences.
    ## Sample 19 - 6504 reads in 1502 unique sequences.
    ## Sample 20 - 4314 reads in 732 unique sequences.

Cette commande nous permet de visualiser le résultat global qu’on
retrouve classé dans la liste dadaFs. Ils nous indiquent que sur les
séquences on retrouve 128 séquences qui correspondent aux vrais
variants, par rapport aux 1979 séquences. Ils nous indiquent aussi les
diagnostiques de qualité.

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 128 sequence variants were inferred from 1979 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# Merge paired reads

Ici nous voulons mettre en une seule séquence les Forwards et les
Reverses.Nous pouvons faire cette opération grâce aux overlaps de 12
paires de base. Cela se fait grâce à un alignement entre les forwards et
les reverses qui vont permettre de contruire les contigs.

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 6540 paired-reads (in 107 unique pairings) successfully merged out of 6891 (in 197 pairings) input.

    ## 5028 paired-reads (in 101 unique pairings) successfully merged out of 5190 (in 157 pairings) input.

    ## 4986 paired-reads (in 81 unique pairings) successfully merged out of 5267 (in 166 pairings) input.

    ## 2595 paired-reads (in 52 unique pairings) successfully merged out of 2754 (in 108 pairings) input.

    ## 2553 paired-reads (in 60 unique pairings) successfully merged out of 2785 (in 119 pairings) input.

    ## 3646 paired-reads (in 55 unique pairings) successfully merged out of 4109 (in 157 pairings) input.

    ## 6079 paired-reads (in 81 unique pairings) successfully merged out of 6514 (in 198 pairings) input.

    ## 3968 paired-reads (in 91 unique pairings) successfully merged out of 4388 (in 187 pairings) input.

    ## 14233 paired-reads (in 143 unique pairings) successfully merged out of 15355 (in 352 pairings) input.

    ## 10528 paired-reads (in 120 unique pairings) successfully merged out of 11165 (in 278 pairings) input.

    ## 11154 paired-reads (in 137 unique pairings) successfully merged out of 11797 (in 298 pairings) input.

    ## 4349 paired-reads (in 85 unique pairings) successfully merged out of 4802 (in 179 pairings) input.

    ## 17431 paired-reads (in 153 unique pairings) successfully merged out of 17812 (in 272 pairings) input.

    ## 5850 paired-reads (in 81 unique pairings) successfully merged out of 6095 (in 159 pairings) input.

    ## 3716 paired-reads (in 86 unique pairings) successfully merged out of 3894 (in 147 pairings) input.

    ## 6865 paired-reads (in 99 unique pairings) successfully merged out of 7191 (in 187 pairings) input.

    ## 4426 paired-reads (in 67 unique pairings) successfully merged out of 4603 (in 127 pairings) input.

    ## 4576 paired-reads (in 101 unique pairings) successfully merged out of 4739 (in 174 pairings) input.

    ## 6092 paired-reads (in 109 unique pairings) successfully merged out of 6315 (in 173 pairings) input.

    ## 4269 paired-reads (in 20 unique pairings) successfully merged out of 4281 (in 28 pairings) input.

``` r
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                       sequence
    ## 1 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACAGG
    ## 2 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAGG
    ## 3 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTGTTAAGTCAGCGGTCAAATGTCGGGGCTCAACCCCGGCCTGCCGTTGAAACTGGCGGCCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ## 4 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTTTTAAGTCAGCGGTAAAAATTCGGGGCTCAACCCCGTCCGGCCGTTGAAACTGGGGGCCTTGAGTGGGCGAGAAGAAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCCTTCCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCGAACAGG
    ## 5 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGACTCTCAAGTCAGCGGTCAAATCGCGGGGCTCAACCCCGTTCCGCCGTTGAAACTGGGAGCCTTGAGTGCGCGAGAAGTAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCCTACCGGCGCGCAACTGACGCTCATGCACGAAAGCGTGGGTATCGAACAGG
    ## 6 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGCCAAGTCAGCGGTAAAAAAGCGGTGCTCAACGCCGTCGAGCCGTTGAAACTGGCGTTCTTGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1       579       1       1    148         0      0      1   TRUE
    ## 2       470       2       2    148         0      0      2   TRUE
    ## 3       449       3       4    148         0      0      1   TRUE
    ## 4       430       4       3    148         0      0      2   TRUE
    ## 5       345       5       6    148         0      0      1   TRUE
    ## 6       282       6       5    148         0      0      2   TRUE

# Construct sequence table

Nous allons construire une table des variations de séquence dans les
amplicons (ASV) qui permet une meilleure résolution que les tables OTUs
97%

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]  20 293

``` r
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```

    ## 
    ## 251 252 253 254 255 
    ##   1  88 196   6   2

# Remove chimeras

Malgré qu’on ait pu appliquer les modèles d’erreurs aux séquences, il
reste des chimères. Ces chimères sont facilement reconnaissables par la
machine et peuvent etre réparées en y rajoutant les parties droites et
gauche des 2 séquences les plus abondantes.

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 61 bimeras out of 293 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]  20 232

Ici on peut voir qu’on à 3% de séquences chimériques dans notre jeu de
donnée.

``` r
1-sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.03596257

# Track reads through the pipeline

Ce code nous permet de visualiser le nombre de séquences obtenues à la
suite de toutes nos manipulations de filtrage. Ici nous pouvons voir
qu’on a pu récupérer la plupart de nos séquences brutes, ce qui est
signe d’une bonne qualité de séquençage.

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##        input filtered denoisedF denoisedR merged nonchim
    ## F3D0    7793     7113      6976      6979   6540    6528
    ## F3D1    5869     5299      5227      5239   5028    5017
    ## F3D141  5958     5463      5331      5357   4986    4863
    ## F3D142  3183     2914      2799      2830   2595    2521
    ## F3D143  3178     2941      2822      2868   2553    2519
    ## F3D144  4827     4312      4151      4228   3646    3507

# Assign taxonomy

Ici nous avons du récupérer silva afin d’analyser et d’assigner les
taxonomies.

``` bash
wget https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz
```

    ## --2020-11-27 16:31:46--  https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz
    ## Resolving zenodo.org (zenodo.org)... 137.138.76.77
    ## Connecting to zenodo.org (zenodo.org)|137.138.76.77|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 81840166 (78M) [application/octet-stream]
    ## Saving to: ‘silva_species_assignment_v138.fa.gz.5’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0% 7.42M 11s
    ##     50K .......... .......... .......... .......... ..........  0% 10.1M 9s
    ##    100K .......... .......... .......... .......... ..........  0% 15.3M 8s
    ##    150K .......... .......... .......... .......... ..........  0% 36.5M 6s
    ##    200K .......... .......... .......... .......... ..........  0% 15.5M 6s
    ##    250K .......... .......... .......... .......... ..........  0% 51.7M 5s
    ##    300K .......... .......... .......... .......... ..........  0% 14.3M 5s
    ##    350K .......... .......... .......... .......... ..........  0% 58.3M 5s
    ##    400K .......... .......... .......... .......... ..........  0% 15.5M 5s
    ##    450K .......... .......... .......... .......... ..........  0% 57.4M 4s
    ##    500K .......... .......... .......... .......... ..........  0% 18.8M 4s
    ##    550K .......... .......... .......... .......... ..........  0% 54.6M 4s
    ##    600K .......... .......... .......... .......... ..........  0% 14.1M 4s
    ##    650K .......... .......... .......... .......... ..........  0% 18.2M 4s
    ##    700K .......... .......... .......... .......... ..........  0% 29.1M 4s
    ##    750K .......... .......... .......... .......... ..........  1% 34.4M 4s
    ##    800K .......... .......... .......... .......... ..........  1% 27.4M 4s
    ##    850K .......... .......... .......... .......... ..........  1% 16.6M 4s
    ##    900K .......... .......... .......... .......... ..........  1% 38.5M 4s
    ##    950K .......... .......... .......... .......... ..........  1% 14.5M 4s
    ##   1000K .......... .......... .......... .......... ..........  1% 34.2M 4s
    ##   1050K .......... .......... .......... .......... ..........  1% 16.0M 4s
    ##   1100K .......... .......... .......... .......... ..........  1% 45.7M 4s
    ##   1150K .......... .......... .......... .......... ..........  1% 33.0M 4s
    ##   1200K .......... .......... .......... .......... ..........  1% 22.2M 4s
    ##   1250K .......... .......... .......... .......... ..........  1% 22.1M 4s
    ##   1300K .......... .......... .......... .......... ..........  1% 42.0M 4s
    ##   1350K .......... .......... .......... .......... ..........  1% 18.2M 4s
    ##   1400K .......... .......... .......... .......... ..........  1% 6.41M 4s
    ##   1450K .......... .......... .......... .......... ..........  1% 39.2M 4s
    ##   1500K .......... .......... .......... .......... ..........  1% 48.4M 4s
    ##   1550K .......... .......... .......... .......... ..........  2% 15.0M 4s
    ##   1600K .......... .......... .......... .......... ..........  2% 41.4M 4s
    ##   1650K .......... .......... .......... .......... ..........  2% 32.9M 4s
    ##   1700K .......... .......... .......... .......... ..........  2% 17.7M 4s
    ##   1750K .......... .......... .......... .......... ..........  2% 46.5M 4s
    ##   1800K .......... .......... .......... .......... ..........  2% 16.3M 4s
    ##   1850K .......... .......... .......... .......... ..........  2% 49.0M 4s
    ##   1900K .......... .......... .......... .......... ..........  2% 26.6M 4s
    ##   1950K .......... .......... .......... .......... ..........  2% 31.2M 4s
    ##   2000K .......... .......... .......... .......... ..........  2% 40.4M 4s
    ##   2050K .......... .......... .......... .......... ..........  2% 7.38M 4s
    ##   2100K .......... .......... .......... .......... ..........  2% 45.8M 4s
    ##   2150K .......... .......... .......... .......... ..........  2% 49.3M 4s
    ##   2200K .......... .......... .......... .......... ..........  2% 22.2M 4s
    ##   2250K .......... .......... .......... .......... ..........  2% 43.7M 4s
    ##   2300K .......... .......... .......... .......... ..........  2% 40.9M 4s
    ##   2350K .......... .......... .......... .......... ..........  3% 21.3M 4s
    ##   2400K .......... .......... .......... .......... ..........  3% 39.5M 3s
    ##   2450K .......... .......... .......... .......... ..........  3% 37.6M 3s
    ##   2500K .......... .......... .......... .......... ..........  3% 35.9M 3s
    ##   2550K .......... .......... .......... .......... ..........  3% 13.6M 3s
    ##   2600K .......... .......... .......... .......... ..........  3% 32.2M 3s
    ##   2650K .......... .......... .......... .......... ..........  3% 49.2M 3s
    ##   2700K .......... .......... .......... .......... ..........  3% 30.8M 3s
    ##   2750K .......... .......... .......... .......... ..........  3% 44.9M 3s
    ##   2800K .......... .......... .......... .......... ..........  3% 33.7M 3s
    ##   2850K .......... .......... .......... .......... ..........  3% 39.5M 3s
    ##   2900K .......... .......... .......... .......... ..........  3% 44.3M 3s
    ##   2950K .......... .......... .......... .......... ..........  3% 19.0M 3s
    ##   3000K .......... .......... .......... .......... ..........  3% 38.1M 3s
    ##   3050K .......... .......... .......... .......... ..........  3% 15.5M 3s
    ##   3100K .......... .......... .......... .......... ..........  3% 32.6M 3s
    ##   3150K .......... .......... .......... .......... ..........  4% 46.8M 3s
    ##   3200K .......... .......... .......... .......... ..........  4% 30.8M 3s
    ##   3250K .......... .......... .......... .......... ..........  4% 44.3M 3s
    ##   3300K .......... .......... .......... .......... ..........  4% 42.2M 3s
    ##   3350K .......... .......... .......... .......... ..........  4% 36.1M 3s
    ##   3400K .......... .......... .......... .......... ..........  4% 39.6M 3s
    ##   3450K .......... .......... .......... .......... ..........  4% 26.5M 3s
    ##   3500K .......... .......... .......... .......... ..........  4% 34.4M 3s
    ##   3550K .......... .......... .......... .......... ..........  4% 51.0M 3s
    ##   3600K .......... .......... .......... .......... ..........  4% 32.5M 3s
    ##   3650K .......... .......... .......... .......... ..........  4% 42.7M 3s
    ##   3700K .......... .......... .......... .......... ..........  4% 21.3M 3s
    ##   3750K .......... .......... .......... .......... ..........  4% 42.9M 3s
    ##   3800K .......... .......... .......... .......... ..........  4% 44.4M 3s
    ##   3850K .......... .......... .......... .......... ..........  4% 28.6M 3s
    ##   3900K .......... .......... .......... .......... ..........  4% 35.8M 3s
    ##   3950K .......... .......... .......... .......... ..........  5% 50.8M 3s
    ##   4000K .......... .......... .......... .......... ..........  5% 35.9M 3s
    ##   4050K .......... .......... .......... .......... ..........  5% 31.8M 3s
    ##   4100K .......... .......... .......... .......... ..........  5% 39.7M 3s
    ##   4150K .......... .......... .......... .......... ..........  5% 39.3M 3s
    ##   4200K .......... .......... .......... .......... ..........  5% 31.7M 3s
    ##   4250K .......... .......... .......... .......... ..........  5% 50.5M 3s
    ##   4300K .......... .......... .......... .......... ..........  5% 33.9M 3s
    ##   4350K .......... .......... .......... .......... ..........  5% 37.6M 3s
    ##   4400K .......... .......... .......... .......... ..........  5% 31.2M 3s
    ##   4450K .......... .......... .......... .......... ..........  5% 48.9M 3s
    ##   4500K .......... .......... .......... .......... ..........  5% 44.1M 3s
    ##   4550K .......... .......... .......... .......... ..........  5% 31.0M 3s
    ##   4600K .......... .......... .......... .......... ..........  5% 41.1M 3s
    ##   4650K .......... .......... .......... .......... ..........  5% 22.3M 3s
    ##   4700K .......... .......... .......... .......... ..........  5% 38.1M 3s
    ##   4750K .......... .......... .......... .......... ..........  6% 53.8M 3s
    ##   4800K .......... .......... .......... .......... ..........  6% 43.4M 3s
    ##   4850K .......... .......... .......... .......... ..........  6% 54.2M 3s
    ##   4900K .......... .......... .......... .......... ..........  6% 43.2M 3s
    ##   4950K .......... .......... .......... .......... ..........  6% 39.1M 3s
    ##   5000K .......... .......... .......... .......... ..........  6% 41.2M 3s
    ##   5050K .......... .......... .......... .......... ..........  6% 51.3M 3s
    ##   5100K .......... .......... .......... .......... ..........  6% 52.5M 3s
    ##   5150K .......... .......... .......... .......... ..........  6% 45.1M 3s
    ##   5200K .......... .......... .......... .......... ..........  6% 48.0M 3s
    ##   5250K .......... .......... .......... .......... ..........  6% 54.4M 3s
    ##   5300K .......... .......... .......... .......... ..........  6% 48.7M 3s
    ##   5350K .......... .......... .......... .......... ..........  6% 50.5M 3s
    ##   5400K .......... .......... .......... .......... ..........  6% 51.3M 3s
    ##   5450K .......... .......... .......... .......... ..........  6% 57.0M 3s
    ##   5500K .......... .......... .......... .......... ..........  6% 46.6M 3s
    ##   5550K .......... .......... .......... .......... ..........  7% 54.9M 3s
    ##   5600K .......... .......... .......... .......... ..........  7% 48.1M 3s
    ##   5650K .......... .......... .......... .......... ..........  7% 55.2M 3s
    ##   5700K .......... .......... .......... .......... ..........  7% 52.4M 3s
    ##   5750K .......... .......... .......... .......... ..........  7% 55.9M 3s
    ##   5800K .......... .......... .......... .......... ..........  7% 13.1M 3s
    ##   5850K .......... .......... .......... .......... ..........  7% 38.0M 3s
    ##   5900K .......... .......... .......... .......... ..........  7% 35.2M 3s
    ##   5950K .......... .......... .......... .......... ..........  7% 46.1M 3s
    ##   6000K .......... .......... .......... .......... ..........  7% 32.5M 3s
    ##   6050K .......... .......... .......... .......... ..........  7% 41.9M 3s
    ##   6100K .......... .......... .......... .......... ..........  7% 36.0M 3s
    ##   6150K .......... .......... .......... .......... ..........  7% 33.8M 3s
    ##   6200K .......... .......... .......... .......... ..........  7% 37.1M 3s
    ##   6250K .......... .......... .......... .......... ..........  7% 48.6M 2s
    ##   6300K .......... .......... .......... .......... ..........  7% 41.2M 2s
    ##   6350K .......... .......... .......... .......... ..........  8% 46.2M 2s
    ##   6400K .......... .......... .......... .......... ..........  8% 43.7M 2s
    ##   6450K .......... .......... .......... .......... ..........  8% 47.1M 2s
    ##   6500K .......... .......... .......... .......... ..........  8% 45.7M 2s
    ##   6550K .......... .......... .......... .......... ..........  8% 54.2M 2s
    ##   6600K .......... .......... .......... .......... ..........  8% 43.4M 2s
    ##   6650K .......... .......... .......... .......... ..........  8% 37.4M 2s
    ##   6700K .......... .......... .......... .......... ..........  8% 45.3M 2s
    ##   6750K .......... .......... .......... .......... ..........  8% 44.7M 2s
    ##   6800K .......... .......... .......... .......... ..........  8% 45.1M 2s
    ##   6850K .......... .......... .......... .......... ..........  8% 46.0M 2s
    ##   6900K .......... .......... .......... .......... ..........  8% 41.8M 2s
    ##   6950K .......... .......... .......... .......... ..........  8% 52.4M 2s
    ##   7000K .......... .......... .......... .......... ..........  8% 46.5M 2s
    ##   7050K .......... .......... .......... .......... ..........  8% 43.6M 2s
    ##   7100K .......... .......... .......... .......... ..........  8% 35.6M 2s
    ##   7150K .......... .......... .......... .......... ..........  9% 45.4M 2s
    ##   7200K .......... .......... .......... .......... ..........  9% 39.8M 2s
    ##   7250K .......... .......... .......... .......... ..........  9% 50.5M 2s
    ##   7300K .......... .......... .......... .......... ..........  9% 45.9M 2s
    ##   7350K .......... .......... .......... .......... ..........  9% 49.9M 2s
    ##   7400K .......... .......... .......... .......... ..........  9% 46.1M 2s
    ##   7450K .......... .......... .......... .......... ..........  9% 51.7M 2s
    ##   7500K .......... .......... .......... .......... ..........  9% 49.6M 2s
    ##   7550K .......... .......... .......... .......... ..........  9% 44.3M 2s
    ##   7600K .......... .......... .......... .......... ..........  9% 42.3M 2s
    ##   7650K .......... .......... .......... .......... ..........  9% 53.1M 2s
    ##   7700K .......... .......... .......... .......... ..........  9% 47.5M 2s
    ##   7750K .......... .......... .......... .......... ..........  9% 56.5M 2s
    ##   7800K .......... .......... .......... .......... ..........  9% 50.5M 2s
    ##   7850K .......... .......... .......... .......... ..........  9% 46.2M 2s
    ##   7900K .......... .......... .......... .......... ..........  9% 55.6M 2s
    ##   7950K .......... .......... .......... .......... .......... 10% 47.5M 2s
    ##   8000K .......... .......... .......... .......... .......... 10% 52.0M 2s
    ##   8050K .......... .......... .......... .......... .......... 10% 58.1M 2s
    ##   8100K .......... .......... .......... .......... .......... 10% 49.9M 2s
    ##   8150K .......... .......... .......... .......... .......... 10% 38.6M 2s
    ##   8200K .......... .......... .......... .......... .......... 10% 55.3M 2s
    ##   8250K .......... .......... .......... .......... .......... 10% 58.0M 2s
    ##   8300K .......... .......... .......... .......... .......... 10% 42.0M 2s
    ##   8350K .......... .......... .......... .......... .......... 10% 55.4M 2s
    ##   8400K .......... .......... .......... .......... .......... 10% 51.2M 2s
    ##   8450K .......... .......... .......... .......... .......... 10% 49.3M 2s
    ##   8500K .......... .......... .......... .......... .......... 10% 58.6M 2s
    ##   8550K .......... .......... .......... .......... .......... 10% 59.7M 2s
    ##   8600K .......... .......... .......... .......... .......... 10% 52.2M 2s
    ##   8650K .......... .......... .......... .......... .......... 10% 42.7M 2s
    ##   8700K .......... .......... .......... .......... .......... 10% 61.4M 2s
    ##   8750K .......... .......... .......... .......... .......... 11% 55.4M 2s
    ##   8800K .......... .......... .......... .......... .......... 11% 61.0M 2s
    ##   8850K .......... .......... .......... .......... .......... 11% 67.0M 2s
    ##   8900K .......... .......... .......... .......... .......... 11% 58.5M 2s
    ##   8950K .......... .......... .......... .......... .......... 11% 55.4M 2s
    ##   9000K .......... .......... .......... .......... .......... 11% 60.4M 2s
    ##   9050K .......... .......... .......... .......... .......... 11% 60.5M 2s
    ##   9100K .......... .......... .......... .......... .......... 11% 47.9M 2s
    ##   9150K .......... .......... .......... .......... .......... 11% 44.1M 2s
    ##   9200K .......... .......... .......... .......... .......... 11% 44.0M 2s
    ##   9250K .......... .......... .......... .......... .......... 11% 60.9M 2s
    ##   9300K .......... .......... .......... .......... .......... 11% 60.1M 2s
    ##   9350K .......... .......... .......... .......... .......... 11% 64.6M 2s
    ##   9400K .......... .......... .......... .......... .......... 11% 60.9M 2s
    ##   9450K .......... .......... .......... .......... .......... 11% 74.9M 2s
    ##   9500K .......... .......... .......... .......... .......... 11% 55.5M 2s
    ##   9550K .......... .......... .......... .......... .......... 12% 65.1M 2s
    ##   9600K .......... .......... .......... .......... .......... 12% 74.2M 2s
    ##   9650K .......... .......... .......... .......... .......... 12% 47.3M 2s
    ##   9700K .......... .......... .......... .......... .......... 12% 46.3M 2s
    ##   9750K .......... .......... .......... .......... .......... 12% 61.5M 2s
    ##   9800K .......... .......... .......... .......... .......... 12% 51.0M 2s
    ##   9850K .......... .......... .......... .......... .......... 12% 64.1M 2s
    ##   9900K .......... .......... .......... .......... .......... 12% 68.4M 2s
    ##   9950K .......... .......... .......... .......... .......... 12% 67.9M 2s
    ##  10000K .......... .......... .......... .......... .......... 12% 70.4M 2s
    ##  10050K .......... .......... .......... .......... .......... 12% 73.1M 2s
    ##  10100K .......... .......... .......... .......... .......... 12% 31.4M 2s
    ##  10150K .......... .......... .......... .......... .......... 12% 59.8M 2s
    ##  10200K .......... .......... .......... .......... .......... 12% 71.7M 2s
    ##  10250K .......... .......... .......... .......... .......... 12% 73.0M 2s
    ##  10300K .......... .......... .......... .......... .......... 12% 24.4M 2s
    ##  10350K .......... .......... .......... .......... .......... 13% 56.0M 2s
    ##  10400K .......... .......... .......... .......... .......... 13% 69.9M 2s
    ##  10450K .......... .......... .......... .......... .......... 13% 35.3M 2s
    ##  10500K .......... .......... .......... .......... .......... 13% 40.0M 2s
    ##  10550K .......... .......... .......... .......... .......... 13% 67.0M 2s
    ##  10600K .......... .......... .......... .......... .......... 13% 59.9M 2s
    ##  10650K .......... .......... .......... .......... .......... 13% 55.3M 2s
    ##  10700K .......... .......... .......... .......... .......... 13% 46.2M 2s
    ##  10750K .......... .......... .......... .......... .......... 13% 35.3M 2s
    ##  10800K .......... .......... .......... .......... .......... 13% 45.3M 2s
    ##  10850K .......... .......... .......... .......... .......... 13% 74.1M 2s
    ##  10900K .......... .......... .......... .......... .......... 13% 61.4M 2s
    ##  10950K .......... .......... .......... .......... .......... 13% 52.4M 2s
    ##  11000K .......... .......... .......... .......... .......... 13% 51.4M 2s
    ##  11050K .......... .......... .......... .......... .......... 13% 70.6M 2s
    ##  11100K .......... .......... .......... .......... .......... 13% 50.5M 2s
    ##  11150K .......... .......... .......... .......... .......... 14% 30.8M 2s
    ##  11200K .......... .......... .......... .......... .......... 14% 66.4M 2s
    ##  11250K .......... .......... .......... .......... .......... 14% 64.8M 2s
    ##  11300K .......... .......... .......... .......... .......... 14% 34.0M 2s
    ##  11350K .......... .......... .......... .......... .......... 14% 34.2M 2s
    ##  11400K .......... .......... .......... .......... .......... 14% 73.3M 2s
    ##  11450K .......... .......... .......... .......... .......... 14% 60.2M 2s
    ##  11500K .......... .......... .......... .......... .......... 14% 67.5M 2s
    ##  11550K .......... .......... .......... .......... .......... 14% 61.3M 2s
    ##  11600K .......... .......... .......... .......... .......... 14% 62.1M 2s
    ##  11650K .......... .......... .......... .......... .......... 14% 63.3M 2s
    ##  11700K .......... .......... .......... .......... .......... 14% 66.0M 2s
    ##  11750K .......... .......... .......... .......... .......... 14% 48.0M 2s
    ##  11800K .......... .......... .......... .......... .......... 14% 58.4M 2s
    ##  11850K .......... .......... .......... .......... .......... 14% 23.2M 2s
    ##  11900K .......... .......... .......... .......... .......... 14% 81.6M 2s
    ##  11950K .......... .......... .......... .......... .......... 15% 59.2M 2s
    ##  12000K .......... .......... .......... .......... .......... 15% 53.0M 2s
    ##  12050K .......... .......... .......... .......... .......... 15% 42.0M 2s
    ##  12100K .......... .......... .......... .......... .......... 15% 68.2M 2s
    ##  12150K .......... .......... .......... .......... .......... 15% 74.8M 2s
    ##  12200K .......... .......... .......... .......... .......... 15% 55.9M 2s
    ##  12250K .......... .......... .......... .......... .......... 15% 90.7M 2s
    ##  12300K .......... .......... .......... .......... .......... 15% 30.8M 2s
    ##  12350K .......... .......... .......... .......... .......... 15% 70.0M 2s
    ##  12400K .......... .......... .......... .......... .......... 15% 53.3M 2s
    ##  12450K .......... .......... .......... .......... .......... 15% 42.7M 2s
    ##  12500K .......... .......... .......... .......... .......... 15% 64.0M 2s
    ##  12550K .......... .......... .......... .......... .......... 15% 87.5M 2s
    ##  12600K .......... .......... .......... .......... .......... 15% 61.2M 2s
    ##  12650K .......... .......... .......... .......... .......... 15% 31.7M 2s
    ##  12700K .......... .......... .......... .......... .......... 15% 64.1M 2s
    ##  12750K .......... .......... .......... .......... .......... 16% 79.1M 2s
    ##  12800K .......... .......... .......... .......... .......... 16% 66.1M 2s
    ##  12850K .......... .......... .......... .......... .......... 16% 33.4M 2s
    ##  12900K .......... .......... .......... .......... .......... 16% 67.3M 2s
    ##  12950K .......... .......... .......... .......... .......... 16% 56.2M 2s
    ##  13000K .......... .......... .......... .......... .......... 16% 67.1M 2s
    ##  13050K .......... .......... .......... .......... .......... 16% 23.1M 2s
    ##  13100K .......... .......... .......... .......... .......... 16%  106M 2s
    ##  13150K .......... .......... .......... .......... .......... 16% 77.9M 2s
    ##  13200K .......... .......... .......... .......... .......... 16% 77.3M 2s
    ##  13250K .......... .......... .......... .......... .......... 16% 39.3M 2s
    ##  13300K .......... .......... .......... .......... .......... 16% 64.0M 2s
    ##  13350K .......... .......... .......... .......... .......... 16% 58.6M 2s
    ##  13400K .......... .......... .......... .......... .......... 16% 49.8M 2s
    ##  13450K .......... .......... .......... .......... .......... 16% 60.6M 2s
    ##  13500K .......... .......... .......... .......... .......... 16% 38.6M 2s
    ##  13550K .......... .......... .......... .......... .......... 17% 47.9M 2s
    ##  13600K .......... .......... .......... .......... .......... 17% 43.3M 2s
    ##  13650K .......... .......... .......... .......... .......... 17% 62.1M 2s
    ##  13700K .......... .......... .......... .......... .......... 17% 91.9M 2s
    ##  13750K .......... .......... .......... .......... .......... 17% 43.9M 2s
    ##  13800K .......... .......... .......... .......... .......... 17% 71.1M 2s
    ##  13850K .......... .......... .......... .......... .......... 17% 37.8M 2s
    ##  13900K .......... .......... .......... .......... .......... 17% 73.5M 2s
    ##  13950K .......... .......... .......... .......... .......... 17% 54.8M 2s
    ##  14000K .......... .......... .......... .......... .......... 17% 99.3M 2s
    ##  14050K .......... .......... .......... .......... .......... 17% 30.1M 2s
    ##  14100K .......... .......... .......... .......... .......... 17% 44.3M 2s
    ##  14150K .......... .......... .......... .......... .......... 17% 71.5M 2s
    ##  14200K .......... .......... .......... .......... .......... 17% 63.8M 2s
    ##  14250K .......... .......... .......... .......... .......... 17% 30.5M 2s
    ##  14300K .......... .......... .......... .......... .......... 17% 10.6M 2s
    ##  14350K .......... .......... .......... .......... .......... 18% 55.0M 2s
    ##  14400K .......... .......... .......... .......... .......... 18% 54.6M 2s
    ##  14450K .......... .......... .......... .......... .......... 18% 91.5M 2s
    ##  14500K .......... .......... .......... .......... .......... 18% 59.3M 2s
    ##  14550K .......... .......... .......... .......... .......... 18% 34.6M 2s
    ##  14600K .......... .......... .......... .......... .......... 18% 71.3M 2s
    ##  14650K .......... .......... .......... .......... .......... 18% 91.8M 2s
    ##  14700K .......... .......... .......... .......... .......... 18% 45.6M 2s
    ##  14750K .......... .......... .......... .......... .......... 18% 36.9M 2s
    ##  14800K .......... .......... .......... .......... .......... 18% 42.3M 2s
    ##  14850K .......... .......... .......... .......... .......... 18%  122M 2s
    ##  14900K .......... .......... .......... .......... .......... 18% 42.3M 2s
    ##  14950K .......... .......... .......... .......... .......... 18% 25.5M 2s
    ##  15000K .......... .......... .......... .......... .......... 18% 20.5M 2s
    ##  15050K .......... .......... .......... .......... .......... 18% 37.6M 2s
    ##  15100K .......... .......... .......... .......... .......... 18% 22.5M 2s
    ##  15150K .......... .......... .......... .......... .......... 19% 39.3M 2s
    ##  15200K .......... .......... .......... .......... .......... 19% 22.1M 2s
    ##  15250K .......... .......... .......... .......... .......... 19%  143M 2s
    ##  15300K .......... .......... .......... .......... .......... 19% 21.8M 2s
    ##  15350K .......... .......... .......... .......... .......... 19% 39.7M 2s
    ##  15400K .......... .......... .......... .......... .......... 19% 37.5M 2s
    ##  15450K .......... .......... .......... .......... .......... 19% 22.9M 2s
    ##  15500K .......... .......... .......... .......... .......... 19%  123M 2s
    ##  15550K .......... .......... .......... .......... .......... 19% 38.5M 2s
    ##  15600K .......... .......... .......... .......... .......... 19%  125M 2s
    ##  15650K .......... .......... .......... .......... .......... 19%  135M 2s
    ##  15700K .......... .......... .......... .......... .......... 19% 38.5M 2s
    ##  15750K .......... .......... .......... .......... .......... 19%  151M 2s
    ##  15800K .......... .......... .......... .......... .......... 19% 38.8M 2s
    ##  15850K .......... .......... .......... .......... .......... 19%  126M 2s
    ##  15900K .......... .......... .......... .......... .......... 19%  131M 2s
    ##  15950K .......... .......... .......... .......... .......... 20% 17.8M 2s
    ##  16000K .......... .......... .......... .......... .......... 20%  110M 2s
    ##  16050K .......... .......... .......... .......... .......... 20%  143M 2s
    ##  16100K .......... .......... .......... .......... .......... 20%  132M 2s
    ##  16150K .......... .......... .......... .......... .......... 20%  118M 2s
    ##  16200K .......... .......... .......... .......... .......... 20%  124M 2s
    ##  16250K .......... .......... .......... .......... .......... 20% 13.7M 2s
    ##  16300K .......... .......... .......... .......... .......... 20%  114M 2s
    ##  16350K .......... .......... .......... .......... .......... 20%  163M 2s
    ##  16400K .......... .......... .......... .......... .......... 20%  143M 2s
    ##  16450K .......... .......... .......... .......... .......... 20%  151M 2s
    ##  16500K .......... .......... .......... .......... .......... 20%  137M 2s
    ##  16550K .......... .......... .......... .......... .......... 20% 17.3M 2s
    ##  16600K .......... .......... .......... .......... .......... 20%  110M 2s
    ##  16650K .......... .......... .......... .......... .......... 20%  142M 2s
    ##  16700K .......... .......... .......... .......... .......... 20% 99.5M 2s
    ##  16750K .......... .......... .......... .......... .......... 21%  127M 2s
    ##  16800K .......... .......... .......... .......... .......... 21% 29.7M 2s
    ##  16850K .......... .......... .......... .......... .......... 21% 61.6M 2s
    ##  16900K .......... .......... .......... .......... .......... 21%  115M 2s
    ##  16950K .......... .......... .......... .......... .......... 21%  140M 2s
    ##  17000K .......... .......... .......... .......... .......... 21% 46.3M 2s
    ##  17050K .......... .......... .......... .......... .......... 21% 47.5M 2s
    ##  17100K .......... .......... .......... .......... .......... 21% 68.1M 2s
    ##  17150K .......... .......... .......... .......... .......... 21% 46.6M 2s
    ##  17200K .......... .......... .......... .......... .......... 21% 50.7M 2s
    ##  17250K .......... .......... .......... .......... .......... 21%  103M 2s
    ##  17300K .......... .......... .......... .......... .......... 21% 52.7M 2s
    ##  17350K .......... .......... .......... .......... .......... 21% 38.0M 2s
    ##  17400K .......... .......... .......... .......... .......... 21% 82.7M 2s
    ##  17450K .......... .......... .......... .......... .......... 21% 72.6M 2s
    ##  17500K .......... .......... .......... .......... .......... 21% 47.4M 2s
    ##  17550K .......... .......... .......... .......... .......... 22% 35.5M 2s
    ##  17600K .......... .......... .......... .......... .......... 22% 22.8M 2s
    ##  17650K .......... .......... .......... .......... .......... 22%  133M 2s
    ##  17700K .......... .......... .......... .......... .......... 22%  116M 2s
    ##  17750K .......... .......... .......... .......... .......... 22% 13.1M 2s
    ##  17800K .......... .......... .......... .......... .......... 22% 41.7M 2s
    ##  17850K .......... .......... .......... .......... .......... 22%  139M 2s
    ##  17900K .......... .......... .......... .......... .......... 22%  116M 2s
    ##  17950K .......... .......... .......... .......... .......... 22% 61.2M 2s
    ##  18000K .......... .......... .......... .......... .......... 22% 55.3M 2s
    ##  18050K .......... .......... .......... .......... .......... 22% 58.8M 2s
    ##  18100K .......... .......... .......... .......... .......... 22% 73.0M 2s
    ##  18150K .......... .......... .......... .......... .......... 22%  101M 2s
    ##  18200K .......... .......... .......... .......... .......... 22% 62.5M 2s
    ##  18250K .......... .......... .......... .......... .......... 22%  108M 2s
    ##  18300K .......... .......... .......... .......... .......... 22% 41.2M 2s
    ##  18350K .......... .......... .......... .......... .......... 23% 30.3M 2s
    ##  18400K .......... .......... .......... .......... .......... 23% 42.4M 2s
    ##  18450K .......... .......... .......... .......... .......... 23%  139M 2s
    ##  18500K .......... .......... .......... .......... .......... 23% 55.9M 2s
    ##  18550K .......... .......... .......... .......... .......... 23% 38.7M 2s
    ##  18600K .......... .......... .......... .......... .......... 23% 67.8M 1s
    ##  18650K .......... .......... .......... .......... .......... 23% 63.2M 1s
    ##  18700K .......... .......... .......... .......... .......... 23% 36.8M 1s
    ##  18750K .......... .......... .......... .......... .......... 23% 92.1M 1s
    ##  18800K .......... .......... .......... .......... .......... 23% 56.0M 1s
    ##  18850K .......... .......... .......... .......... .......... 23% 18.0M 1s
    ##  18900K .......... .......... .......... .......... .......... 23% 80.3M 1s
    ##  18950K .......... .......... .......... .......... .......... 23% 83.5M 1s
    ##  19000K .......... .......... .......... .......... .......... 23% 65.8M 1s
    ##  19050K .......... .......... .......... .......... .......... 23% 98.4M 1s
    ##  19100K .......... .......... .......... .......... .......... 23%  105M 1s
    ##  19150K .......... .......... .......... .......... .......... 24% 19.1M 1s
    ##  19200K .......... .......... .......... .......... .......... 24%  103M 1s
    ##  19250K .......... .......... .......... .......... .......... 24%  153M 1s
    ##  19300K .......... .......... .......... .......... .......... 24%  141M 1s
    ##  19350K .......... .......... .......... .......... .......... 24% 13.6M 1s
    ##  19400K .......... .......... .......... .......... .......... 24%  138M 1s
    ##  19450K .......... .......... .......... .......... .......... 24% 52.0M 1s
    ##  19500K .......... .......... .......... .......... .......... 24% 63.9M 1s
    ##  19550K .......... .......... .......... .......... .......... 24%  128M 1s
    ##  19600K .......... .......... .......... .......... .......... 24% 28.4M 1s
    ##  19650K .......... .......... .......... .......... .......... 24% 49.7M 1s
    ##  19700K .......... .......... .......... .......... .......... 24% 55.5M 1s
    ##  19750K .......... .......... .......... .......... .......... 24% 46.8M 1s
    ##  19800K .......... .......... .......... .......... .......... 24% 95.8M 1s
    ##  19850K .......... .......... .......... .......... .......... 24% 60.0M 1s
    ##  19900K .......... .......... .......... .......... .......... 24% 34.8M 1s
    ##  19950K .......... .......... .......... .......... .......... 25% 48.6M 1s
    ##  20000K .......... .......... .......... .......... .......... 25% 45.1M 1s
    ##  20050K .......... .......... .......... .......... .......... 25% 81.2M 1s
    ##  20100K .......... .......... .......... .......... .......... 25% 54.3M 1s
    ##  20150K .......... .......... .......... .......... .......... 25% 52.7M 1s
    ##  20200K .......... .......... .......... .......... .......... 25% 35.7M 1s
    ##  20250K .......... .......... .......... .......... .......... 25% 46.8M 1s
    ##  20300K .......... .......... .......... .......... .......... 25% 39.5M 1s
    ##  20350K .......... .......... .......... .......... .......... 25% 56.8M 1s
    ##  20400K .......... .......... .......... .......... .......... 25% 87.1M 1s
    ##  20450K .......... .......... .......... .......... .......... 25% 62.1M 1s
    ##  20500K .......... .......... .......... .......... .......... 25% 51.8M 1s
    ##  20550K .......... .......... .......... .......... .......... 25% 52.0M 1s
    ##  20600K .......... .......... .......... .......... .......... 25% 45.5M 1s
    ##  20650K .......... .......... .......... .......... .......... 25%  102M 1s
    ##  20700K .......... .......... .......... .......... .......... 25% 58.2M 1s
    ##  20750K .......... .......... .......... .......... .......... 26% 46.1M 1s
    ##  20800K .......... .......... .......... .......... .......... 26% 50.7M 1s
    ##  20850K .......... .......... .......... .......... .......... 26% 79.4M 1s
    ##  20900K .......... .......... .......... .......... .......... 26% 41.0M 1s
    ##  20950K .......... .......... .......... .......... .......... 26% 63.9M 1s
    ##  21000K .......... .......... .......... .......... .......... 26% 57.6M 1s
    ##  21050K .......... .......... .......... .......... .......... 26% 96.8M 1s
    ##  21100K .......... .......... .......... .......... .......... 26% 66.0M 1s
    ##  21150K .......... .......... .......... .......... .......... 26% 43.8M 1s
    ##  21200K .......... .......... .......... .......... .......... 26% 45.1M 1s
    ##  21250K .......... .......... .......... .......... .......... 26% 85.0M 1s
    ##  21300K .......... .......... .......... .......... .......... 26% 60.4M 1s
    ##  21350K .......... .......... .......... .......... .......... 26% 52.1M 1s
    ##  21400K .......... .......... .......... .......... .......... 26% 37.3M 1s
    ##  21450K .......... .......... .......... .......... .......... 26% 47.9M 1s
    ##  21500K .......... .......... .......... .......... .......... 26% 94.3M 1s
    ##  21550K .......... .......... .......... .......... .......... 27% 69.2M 1s
    ##  21600K .......... .......... .......... .......... .......... 27% 35.3M 1s
    ##  21650K .......... .......... .......... .......... .......... 27% 80.8M 1s
    ##  21700K .......... .......... .......... .......... .......... 27% 67.7M 1s
    ##  21750K .......... .......... .......... .......... .......... 27% 54.2M 1s
    ##  21800K .......... .......... .......... .......... .......... 27% 57.3M 1s
    ##  21850K .......... .......... .......... .......... .......... 27% 68.5M 1s
    ##  21900K .......... .......... .......... .......... .......... 27% 54.4M 1s
    ##  21950K .......... .......... .......... .......... .......... 27% 64.5M 1s
    ##  22000K .......... .......... .......... .......... .......... 27% 79.8M 1s
    ##  22050K .......... .......... .......... .......... .......... 27% 26.4M 1s
    ##  22100K .......... .......... .......... .......... .......... 27% 56.5M 1s
    ##  22150K .......... .......... .......... .......... .......... 27% 41.0M 1s
    ##  22200K .......... .......... .......... .......... .......... 27% 88.1M 1s
    ##  22250K .......... .......... .......... .......... .......... 27%  112M 1s
    ##  22300K .......... .......... .......... .......... .......... 27% 67.2M 1s
    ##  22350K .......... .......... .......... .......... .......... 28% 59.5M 1s
    ##  22400K .......... .......... .......... .......... .......... 28% 65.4M 1s
    ##  22450K .......... .......... .......... .......... .......... 28% 79.3M 1s
    ##  22500K .......... .......... .......... .......... .......... 28% 76.5M 1s
    ##  22550K .......... .......... .......... .......... .......... 28% 87.9M 1s
    ##  22600K .......... .......... .......... .......... .......... 28% 69.0M 1s
    ##  22650K .......... .......... .......... .......... .......... 28% 64.4M 1s
    ##  22700K .......... .......... .......... .......... .......... 28% 44.4M 1s
    ##  22750K .......... .......... .......... .......... .......... 28% 55.3M 1s
    ##  22800K .......... .......... .......... .......... .......... 28% 80.3M 1s
    ##  22850K .......... .......... .......... .......... .......... 28% 75.1M 1s
    ##  22900K .......... .......... .......... .......... .......... 28% 72.5M 1s
    ##  22950K .......... .......... .......... .......... .......... 28% 52.4M 1s
    ##  23000K .......... .......... .......... .......... .......... 28% 52.7M 1s
    ##  23050K .......... .......... .......... .......... .......... 28% 23.1M 1s
    ##  23100K .......... .......... .......... .......... .......... 28% 30.6M 1s
    ##  23150K .......... .......... .......... .......... .......... 29% 61.7M 1s
    ##  23200K .......... .......... .......... .......... .......... 29% 61.5M 1s
    ##  23250K .......... .......... .......... .......... .......... 29% 73.2M 1s
    ##  23300K .......... .......... .......... .......... .......... 29% 50.1M 1s
    ##  23350K .......... .......... .......... .......... .......... 29% 48.9M 1s
    ##  23400K .......... .......... .......... .......... .......... 29% 57.4M 1s
    ##  23450K .......... .......... .......... .......... .......... 29% 30.6M 1s
    ##  23500K .......... .......... .......... .......... .......... 29% 45.4M 1s
    ##  23550K .......... .......... .......... .......... .......... 29% 72.5M 1s
    ##  23600K .......... .......... .......... .......... .......... 29% 57.1M 1s
    ##  23650K .......... .......... .......... .......... .......... 29% 56.5M 1s
    ##  23700K .......... .......... .......... .......... .......... 29% 36.4M 1s
    ##  23750K .......... .......... .......... .......... .......... 29% 47.8M 1s
    ##  23800K .......... .......... .......... .......... .......... 29% 45.8M 1s
    ##  23850K .......... .......... .......... .......... .......... 29% 52.1M 1s
    ##  23900K .......... .......... .......... .......... .......... 29% 47.2M 1s
    ##  23950K .......... .......... .......... .......... .......... 30% 51.0M 1s
    ##  24000K .......... .......... .......... .......... .......... 30% 44.5M 1s
    ##  24050K .......... .......... .......... .......... .......... 30% 54.4M 1s
    ##  24100K .......... .......... .......... .......... .......... 30% 53.3M 1s
    ##  24150K .......... .......... .......... .......... .......... 30% 70.3M 1s
    ##  24200K .......... .......... .......... .......... .......... 30% 49.7M 1s
    ##  24250K .......... .......... .......... .......... .......... 30% 82.1M 1s
    ##  24300K .......... .......... .......... .......... .......... 30% 62.5M 1s
    ##  24350K .......... .......... .......... .......... .......... 30% 87.5M 1s
    ##  24400K .......... .......... .......... .......... .......... 30% 73.6M 1s
    ##  24450K .......... .......... .......... .......... .......... 30% 48.2M 1s
    ##  24500K .......... .......... .......... .......... .......... 30% 78.4M 1s
    ##  24550K .......... .......... .......... .......... .......... 30%  101M 1s
    ##  24600K .......... .......... .......... .......... .......... 30% 81.4M 1s
    ##  24650K .......... .......... .......... .......... .......... 30% 18.8M 1s
    ##  24700K .......... .......... .......... .......... .......... 30% 59.2M 1s
    ##  24750K .......... .......... .......... .......... .......... 31% 25.0M 1s
    ##  24800K .......... .......... .......... .......... .......... 31% 48.9M 1s
    ##  24850K .......... .......... .......... .......... .......... 31%  106M 1s
    ##  24900K .......... .......... .......... .......... .......... 31% 74.0M 1s
    ##  24950K .......... .......... .......... .......... .......... 31% 55.0M 1s
    ##  25000K .......... .......... .......... .......... .......... 31% 44.4M 1s
    ##  25050K .......... .......... .......... .......... .......... 31% 46.8M 1s
    ##  25100K .......... .......... .......... .......... .......... 31% 32.7M 1s
    ##  25150K .......... .......... .......... .......... .......... 31%  100M 1s
    ##  25200K .......... .......... .......... .......... .......... 31% 44.2M 1s
    ##  25250K .......... .......... .......... .......... .......... 31% 75.2M 1s
    ##  25300K .......... .......... .......... .......... .......... 31% 49.4M 1s
    ##  25350K .......... .......... .......... .......... .......... 31% 54.3M 1s
    ##  25400K .......... .......... .......... .......... .......... 31% 48.2M 1s
    ##  25450K .......... .......... .......... .......... .......... 31% 69.0M 1s
    ##  25500K .......... .......... .......... .......... .......... 31% 78.0M 1s
    ##  25550K .......... .......... .......... .......... .......... 32% 53.6M 1s
    ##  25600K .......... .......... .......... .......... .......... 32% 62.2M 1s
    ##  25650K .......... .......... .......... .......... .......... 32% 57.7M 1s
    ##  25700K .......... .......... .......... .......... .......... 32% 89.5M 1s
    ##  25750K .......... .......... .......... .......... .......... 32%  111M 1s
    ##  25800K .......... .......... .......... .......... .......... 32% 80.4M 1s
    ##  25850K .......... .......... .......... .......... .......... 32%  115M 1s
    ##  25900K .......... .......... .......... .......... .......... 32% 86.7M 1s
    ##  25950K .......... .......... .......... .......... .......... 32% 41.1M 1s
    ##  26000K .......... .......... .......... .......... .......... 32% 52.2M 1s
    ##  26050K .......... .......... .......... .......... .......... 32% 85.9M 1s
    ##  26100K .......... .......... .......... .......... .......... 32% 18.2M 1s
    ##  26150K .......... .......... .......... .......... .......... 32% 76.3M 1s
    ##  26200K .......... .......... .......... .......... .......... 32% 37.3M 1s
    ##  26250K .......... .......... .......... .......... .......... 32% 68.0M 1s
    ##  26300K .......... .......... .......... .......... .......... 32%  100M 1s
    ##  26350K .......... .......... .......... .......... .......... 33% 67.4M 1s
    ##  26400K .......... .......... .......... .......... .......... 33%  102M 1s
    ##  26450K .......... .......... .......... .......... .......... 33% 12.5M 1s
    ##  26500K .......... .......... .......... .......... .......... 33% 69.1M 1s
    ##  26550K .......... .......... .......... .......... .......... 33% 87.2M 1s
    ##  26600K .......... .......... .......... .......... .......... 33% 89.5M 1s
    ##  26650K .......... .......... .......... .......... .......... 33%  126M 1s
    ##  26700K .......... .......... .......... .......... .......... 33% 19.8M 1s
    ##  26750K .......... .......... .......... .......... .......... 33% 92.9M 1s
    ##  26800K .......... .......... .......... .......... .......... 33% 51.2M 1s
    ##  26850K .......... .......... .......... .......... .......... 33% 73.7M 1s
    ##  26900K .......... .......... .......... .......... .......... 33% 78.3M 1s
    ##  26950K .......... .......... .......... .......... .......... 33% 90.8M 1s
    ##  27000K .......... .......... .......... .......... .......... 33% 22.5M 1s
    ##  27050K .......... .......... .......... .......... .......... 33% 88.5M 1s
    ##  27100K .......... .......... .......... .......... .......... 33% 63.9M 1s
    ##  27150K .......... .......... .......... .......... .......... 34% 87.6M 1s
    ##  27200K .......... .......... .......... .......... .......... 34% 70.3M 1s
    ##  27250K .......... .......... .......... .......... .......... 34% 56.5M 1s
    ##  27300K .......... .......... .......... .......... .......... 34% 39.1M 1s
    ##  27350K .......... .......... .......... .......... .......... 34% 83.8M 1s
    ##  27400K .......... .......... .......... .......... .......... 34% 79.3M 1s
    ##  27450K .......... .......... .......... .......... .......... 34% 63.4M 1s
    ##  27500K .......... .......... .......... .......... .......... 34% 60.0M 1s
    ##  27550K .......... .......... .......... .......... .......... 34%  104M 1s
    ##  27600K .......... .......... .......... .......... .......... 34% 73.9M 1s
    ##  27650K .......... .......... .......... .......... .......... 34% 41.6M 1s
    ##  27700K .......... .......... .......... .......... .......... 34% 75.5M 1s
    ##  27750K .......... .......... .......... .......... .......... 34% 49.6M 1s
    ##  27800K .......... .......... .......... .......... .......... 34% 62.6M 1s
    ##  27850K .......... .......... .......... .......... .......... 34% 99.6M 1s
    ##  27900K .......... .......... .......... .......... .......... 34% 70.0M 1s
    ##  27950K .......... .......... .......... .......... .......... 35% 29.2M 1s
    ##  28000K .......... .......... .......... .......... .......... 35% 61.9M 1s
    ##  28050K .......... .......... .......... .......... .......... 35% 64.5M 1s
    ##  28100K .......... .......... .......... .......... .......... 35% 57.9M 1s
    ##  28150K .......... .......... .......... .......... .......... 35% 57.1M 1s
    ##  28200K .......... .......... .......... .......... .......... 35% 68.5M 1s
    ##  28250K .......... .......... .......... .......... .......... 35% 46.4M 1s
    ##  28300K .......... .......... .......... .......... .......... 35% 66.8M 1s
    ##  28350K .......... .......... .......... .......... .......... 35% 74.7M 1s
    ##  28400K .......... .......... .......... .......... .......... 35% 43.3M 1s
    ##  28450K .......... .......... .......... .......... .......... 35% 90.1M 1s
    ##  28500K .......... .......... .......... .......... .......... 35% 92.9M 1s
    ##  28550K .......... .......... .......... .......... .......... 35% 47.4M 1s
    ##  28600K .......... .......... .......... .......... .......... 35% 18.4M 1s
    ##  28650K .......... .......... .......... .......... .......... 35% 74.0M 1s
    ##  28700K .......... .......... .......... .......... .......... 35% 36.6M 1s
    ##  28750K .......... .......... .......... .......... .......... 36% 94.4M 1s
    ##  28800K .......... .......... .......... .......... .......... 36% 96.2M 1s
    ##  28850K .......... .......... .......... .......... .......... 36% 27.5M 1s
    ##  28900K .......... .......... .......... .......... .......... 36% 29.1M 1s
    ##  28950K .......... .......... .......... .......... .......... 36% 61.2M 1s
    ##  29000K .......... .......... .......... .......... .......... 36% 61.9M 1s
    ##  29050K .......... .......... .......... .......... .......... 36% 28.1M 1s
    ##  29100K .......... .......... .......... .......... .......... 36% 88.9M 1s
    ##  29150K .......... .......... .......... .......... .......... 36% 28.5M 1s
    ##  29200K .......... .......... .......... .......... .......... 36% 88.5M 1s
    ##  29250K .......... .......... .......... .......... .......... 36% 95.3M 1s
    ##  29300K .......... .......... .......... .......... .......... 36% 83.5M 1s
    ##  29350K .......... .......... .......... .......... .......... 36% 81.5M 1s
    ##  29400K .......... .......... .......... .......... .......... 36% 38.0M 1s
    ##  29450K .......... .......... .......... .......... .......... 36% 56.7M 1s
    ##  29500K .......... .......... .......... .......... .......... 36% 56.7M 1s
    ##  29550K .......... .......... .......... .......... .......... 37% 91.3M 1s
    ##  29600K .......... .......... .......... .......... .......... 37% 52.9M 1s
    ##  29650K .......... .......... .......... .......... .......... 37% 96.0M 1s
    ##  29700K .......... .......... .......... .......... .......... 37%  105M 1s
    ##  29750K .......... .......... .......... .......... .......... 37% 44.0M 1s
    ##  29800K .......... .......... .......... .......... .......... 37% 94.1M 1s
    ##  29850K .......... .......... .......... .......... .......... 37% 55.2M 1s
    ##  29900K .......... .......... .......... .......... .......... 37% 61.5M 1s
    ##  29950K .......... .......... .......... .......... .......... 37%  103M 1s
    ##  30000K .......... .......... .......... .......... .......... 37% 66.3M 1s
    ##  30050K .......... .......... .......... .......... .......... 37% 45.5M 1s
    ##  30100K .......... .......... .......... .......... .......... 37% 58.6M 1s
    ##  30150K .......... .......... .......... .......... .......... 37% 96.0M 1s
    ##  30200K .......... .......... .......... .......... .......... 37% 65.1M 1s
    ##  30250K .......... .......... .......... .......... .......... 37% 57.9M 1s
    ##  30300K .......... .......... .......... .......... .......... 37% 59.2M 1s
    ##  30350K .......... .......... .......... .......... .......... 38% 56.9M 1s
    ##  30400K .......... .......... .......... .......... .......... 38% 87.9M 1s
    ##  30450K .......... .......... .......... .......... .......... 38% 80.0M 1s
    ##  30500K .......... .......... .......... .......... .......... 38% 47.8M 1s
    ##  30550K .......... .......... .......... .......... .......... 38% 61.1M 1s
    ##  30600K .......... .......... .......... .......... .......... 38% 73.6M 1s
    ##  30650K .......... .......... .......... .......... .......... 38% 47.9M 1s
    ##  30700K .......... .......... .......... .......... .......... 38% 72.4M 1s
    ##  30750K .......... .......... .......... .......... .......... 38% 91.0M 1s
    ##  30800K .......... .......... .......... .......... .......... 38% 55.6M 1s
    ##  30850K .......... .......... .......... .......... .......... 38% 58.3M 1s
    ##  30900K .......... .......... .......... .......... .......... 38% 62.3M 1s
    ##  30950K .......... .......... .......... .......... .......... 38% 53.8M 1s
    ##  31000K .......... .......... .......... .......... .......... 38% 46.9M 1s
    ##  31050K .......... .......... .......... .......... .......... 38% 71.4M 1s
    ##  31100K .......... .......... .......... .......... .......... 38% 77.0M 1s
    ##  31150K .......... .......... .......... .......... .......... 39% 67.2M 1s
    ##  31200K .......... .......... .......... .......... .......... 39% 56.7M 1s
    ##  31250K .......... .......... .......... .......... .......... 39% 51.8M 1s
    ##  31300K .......... .......... .......... .......... .......... 39% 52.8M 1s
    ##  31350K .......... .......... .......... .......... .......... 39%  104M 1s
    ##  31400K .......... .......... .......... .......... .......... 39% 52.3M 1s
    ##  31450K .......... .......... .......... .......... .......... 39% 46.8M 1s
    ##  31500K .......... .......... .......... .......... .......... 39% 67.0M 1s
    ##  31550K .......... .......... .......... .......... .......... 39% 66.8M 1s
    ##  31600K .......... .......... .......... .......... .......... 39% 88.2M 1s
    ##  31650K .......... .......... .......... .......... .......... 39% 53.1M 1s
    ##  31700K .......... .......... .......... .......... .......... 39% 55.0M 1s
    ##  31750K .......... .......... .......... .......... .......... 39% 68.2M 1s
    ##  31800K .......... .......... .......... .......... .......... 39% 56.1M 1s
    ##  31850K .......... .......... .......... .......... .......... 39%  101M 1s
    ##  31900K .......... .......... .......... .......... .......... 39% 55.8M 1s
    ##  31950K .......... .......... .......... .......... .......... 40% 52.2M 1s
    ##  32000K .......... .......... .......... .......... .......... 40% 54.3M 1s
    ##  32050K .......... .......... .......... .......... .......... 40%  117M 1s
    ##  32100K .......... .......... .......... .......... .......... 40% 63.2M 1s
    ##  32150K .......... .......... .......... .......... .......... 40% 64.9M 1s
    ##  32200K .......... .......... .......... .......... .......... 40% 57.2M 1s
    ##  32250K .......... .......... .......... .......... .......... 40% 64.0M 1s
    ##  32300K .......... .......... .......... .......... .......... 40% 83.8M 1s
    ##  32350K .......... .......... .......... .......... .......... 40% 65.6M 1s
    ##  32400K .......... .......... .......... .......... .......... 40% 60.2M 1s
    ##  32450K .......... .......... .......... .......... .......... 40% 46.8M 1s
    ##  32500K .......... .......... .......... .......... .......... 40% 55.6M 1s
    ##  32550K .......... .......... .......... .......... .......... 40% 96.0M 1s
    ##  32600K .......... .......... .......... .......... .......... 40% 51.5M 1s
    ##  32650K .......... .......... .......... .......... .......... 40% 76.4M 1s
    ##  32700K .......... .......... .......... .......... .......... 40% 49.9M 1s
    ##  32750K .......... .......... .......... .......... .......... 41%  101M 1s
    ##  32800K .......... .......... .......... .......... .......... 41% 83.5M 1s
    ##  32850K .......... .......... .......... .......... .......... 41% 54.4M 1s
    ##  32900K .......... .......... .......... .......... .......... 41% 49.2M 1s
    ##  32950K .......... .......... .......... .......... .......... 41% 67.4M 1s
    ##  33000K .......... .......... .......... .......... .......... 41% 95.0M 1s
    ##  33050K .......... .......... .......... .......... .......... 41% 96.2M 1s
    ##  33100K .......... .......... .......... .......... .......... 41% 54.7M 1s
    ##  33150K .......... .......... .......... .......... .......... 41% 70.0M 1s
    ##  33200K .......... .......... .......... .......... .......... 41% 12.3M 1s
    ##  33250K .......... .......... .......... .......... .......... 41%  104M 1s
    ##  33300K .......... .......... .......... .......... .......... 41% 87.8M 1s
    ##  33350K .......... .......... .......... .......... .......... 41%  136M 1s
    ##  33400K .......... .......... .......... .......... .......... 41%  111M 1s
    ##  33450K .......... .......... .......... .......... .......... 41%  116M 1s
    ##  33500K .......... .......... .......... .......... .......... 41% 18.3M 1s
    ##  33550K .......... .......... .......... .......... .......... 42% 84.7M 1s
    ##  33600K .......... .......... .......... .......... .......... 42% 34.1M 1s
    ##  33650K .......... .......... .......... .......... .......... 42% 88.6M 1s
    ##  33700K .......... .......... .......... .......... .......... 42% 80.9M 1s
    ##  33750K .......... .......... .......... .......... .......... 42% 95.9M 1s
    ##  33800K .......... .......... .......... .......... .......... 42% 94.8M 1s
    ##  33850K .......... .......... .......... .......... .......... 42% 24.8M 1s
    ##  33900K .......... .......... .......... .......... .......... 42% 92.4M 1s
    ##  33950K .......... .......... .......... .......... .......... 42% 87.7M 1s
    ##  34000K .......... .......... .......... .......... .......... 42% 65.7M 1s
    ##  34050K .......... .......... .......... .......... .......... 42% 42.2M 1s
    ##  34100K .......... .......... .......... .......... .......... 42%  100M 1s
    ##  34150K .......... .......... .......... .......... .......... 42%  116M 1s
    ##  34200K .......... .......... .......... .......... .......... 42% 40.5M 1s
    ##  34250K .......... .......... .......... .......... .......... 42% 69.9M 1s
    ##  34300K .......... .......... .......... .......... .......... 42% 39.2M 1s
    ##  34350K .......... .......... .......... .......... .......... 43% 88.7M 1s
    ##  34400K .......... .......... .......... .......... .......... 43% 93.0M 1s
    ##  34450K .......... .......... .......... .......... .......... 43%  103M 1s
    ##  34500K .......... .......... .......... .......... .......... 43% 47.9M 1s
    ##  34550K .......... .......... .......... .......... .......... 43% 48.6M 1s
    ##  34600K .......... .......... .......... .......... .......... 43% 25.6M 1s
    ##  34650K .......... .......... .......... .......... .......... 43% 4.98M 1s
    ##  34700K .......... .......... .......... .......... .......... 43% 90.1M 1s
    ##  34750K .......... .......... .......... .......... .......... 43% 95.6M 1s
    ##  34800K .......... .......... .......... .......... .......... 43% 84.3M 1s
    ##  34850K .......... .......... .......... .......... .......... 43%  101M 1s
    ##  34900K .......... .......... .......... .......... .......... 43% 83.5M 1s
    ##  34950K .......... .......... .......... .......... .......... 43% 86.3M 1s
    ##  35000K .......... .......... .......... .......... .......... 43% 85.5M 1s
    ##  35050K .......... .......... .......... .......... .......... 43% 77.3M 1s
    ##  35100K .......... .......... .......... .......... .......... 43% 27.0M 1s
    ##  35150K .......... .......... .......... .......... .......... 44% 77.4M 1s
    ##  35200K .......... .......... .......... .......... .......... 44% 21.5M 1s
    ##  35250K .......... .......... .......... .......... .......... 44% 80.3M 1s
    ##  35300K .......... .......... .......... .......... .......... 44% 99.5M 1s
    ##  35350K .......... .......... .......... .......... .......... 44% 86.6M 1s
    ##  35400K .......... .......... .......... .......... .......... 44% 94.7M 1s
    ##  35450K .......... .......... .......... .......... .......... 44% 38.5M 1s
    ##  35500K .......... .......... .......... .......... .......... 44% 75.1M 1s
    ##  35550K .......... .......... .......... .......... .......... 44% 30.6M 1s
    ##  35600K .......... .......... .......... .......... .......... 44% 69.9M 1s
    ##  35650K .......... .......... .......... .......... .......... 44% 98.1M 1s
    ##  35700K .......... .......... .......... .......... .......... 44% 95.3M 1s
    ##  35750K .......... .......... .......... .......... .......... 44%  118M 1s
    ##  35800K .......... .......... .......... .......... .......... 44% 77.3M 1s
    ##  35850K .......... .......... .......... .......... .......... 44% 69.6M 1s
    ##  35900K .......... .......... .......... .......... .......... 44% 34.9M 1s
    ##  35950K .......... .......... .......... .......... .......... 45% 85.3M 1s
    ##  36000K .......... .......... .......... .......... .......... 45% 68.6M 1s
    ##  36050K .......... .......... .......... .......... .......... 45% 86.6M 1s
    ##  36100K .......... .......... .......... .......... .......... 45% 83.6M 1s
    ##  36150K .......... .......... .......... .......... .......... 45% 37.0M 1s
    ##  36200K .......... .......... .......... .......... .......... 45% 58.7M 1s
    ##  36250K .......... .......... .......... .......... .......... 45%  102M 1s
    ##  36300K .......... .......... .......... .......... .......... 45% 62.3M 1s
    ##  36350K .......... .......... .......... .......... .......... 45% 92.1M 1s
    ##  36400K .......... .......... .......... .......... .......... 45% 55.2M 1s
    ##  36450K .......... .......... .......... .......... .......... 45% 31.7M 1s
    ##  36500K .......... .......... .......... .......... .......... 45% 80.2M 1s
    ##  36550K .......... .......... .......... .......... .......... 45%  107M 1s
    ##  36600K .......... .......... .......... .......... .......... 45% 43.0M 1s
    ##  36650K .......... .......... .......... .......... .......... 45% 55.8M 1s
    ##  36700K .......... .......... .......... .......... .......... 45% 6.58M 1s
    ##  36750K .......... .......... .......... .......... .......... 46% 92.6M 1s
    ##  36800K .......... .......... .......... .......... .......... 46% 99.8M 1s
    ##  36850K .......... .......... .......... .......... .......... 46% 67.2M 1s
    ##  36900K .......... .......... .......... .......... .......... 46% 84.8M 1s
    ##  36950K .......... .......... .......... .......... .......... 46%  102M 1s
    ##  37000K .......... .......... .......... .......... .......... 46% 92.4M 1s
    ##  37050K .......... .......... .......... .......... .......... 46% 86.6M 1s
    ##  37100K .......... .......... .......... .......... .......... 46% 62.5M 1s
    ##  37150K .......... .......... .......... .......... .......... 46% 82.7M 1s
    ##  37200K .......... .......... .......... .......... .......... 46% 40.3M 1s
    ##  37250K .......... .......... .......... .......... .......... 46% 64.2M 1s
    ##  37300K .......... .......... .......... .......... .......... 46% 61.6M 1s
    ##  37350K .......... .......... .......... .......... .......... 46% 84.4M 1s
    ##  37400K .......... .......... .......... .......... .......... 46% 64.0M 1s
    ##  37450K .......... .......... .......... .......... .......... 46% 40.0M 1s
    ##  37500K .......... .......... .......... .......... .......... 46% 83.5M 1s
    ##  37550K .......... .......... .......... .......... .......... 47% 91.5M 1s
    ##  37600K .......... .......... .......... .......... .......... 47% 51.3M 1s
    ##  37650K .......... .......... .......... .......... .......... 47% 99.1M 1s
    ##  37700K .......... .......... .......... .......... .......... 47%  101M 1s
    ##  37750K .......... .......... .......... .......... .......... 47% 13.6M 1s
    ##  37800K .......... .......... .......... .......... .......... 47% 94.8M 1s
    ##  37850K .......... .......... .......... .......... .......... 47% 96.1M 1s
    ##  37900K .......... .......... .......... .......... .......... 47% 75.7M 1s
    ##  37950K .......... .......... .......... .......... .......... 47%  108M 1s
    ##  38000K .......... .......... .......... .......... .......... 47% 15.4M 1s
    ##  38050K .......... .......... .......... .......... .......... 47% 37.8M 1s
    ##  38100K .......... .......... .......... .......... .......... 47% 63.2M 1s
    ##  38150K .......... .......... .......... .......... .......... 47% 61.9M 1s
    ##  38200K .......... .......... .......... .......... .......... 47% 99.4M 1s
    ##  38250K .......... .......... .......... .......... .......... 47% 90.9M 1s
    ##  38300K .......... .......... .......... .......... .......... 47% 69.7M 1s
    ##  38350K .......... .......... .......... .......... .......... 48% 65.1M 1s
    ##  38400K .......... .......... .......... .......... .......... 48% 61.5M 1s
    ##  38450K .......... .......... .......... .......... .......... 48% 97.9M 1s
    ##  38500K .......... .......... .......... .......... .......... 48% 56.1M 1s
    ##  38550K .......... .......... .......... .......... .......... 48% 84.6M 1s
    ##  38600K .......... .......... .......... .......... .......... 48% 65.8M 1s
    ##  38650K .......... .......... .......... .......... .......... 48% 35.5M 1s
    ##  38700K .......... .......... .......... .......... .......... 48% 80.9M 1s
    ##  38750K .......... .......... .......... .......... .......... 48% 91.7M 1s
    ##  38800K .......... .......... .......... .......... .......... 48% 88.6M 1s
    ##  38850K .......... .......... .......... .......... .......... 48% 51.8M 1s
    ##  38900K .......... .......... .......... .......... .......... 48% 75.4M 1s
    ##  38950K .......... .......... .......... .......... .......... 48% 58.4M 1s
    ##  39000K .......... .......... .......... .......... .......... 48% 41.3M 1s
    ##  39050K .......... .......... .......... .......... .......... 48% 36.4M 1s
    ##  39100K .......... .......... .......... .......... .......... 48% 83.8M 1s
    ##  39150K .......... .......... .......... .......... .......... 49% 75.4M 1s
    ##  39200K .......... .......... .......... .......... .......... 49% 96.1M 1s
    ##  39250K .......... .......... .......... .......... .......... 49% 80.5M 1s
    ##  39300K .......... .......... .......... .......... .......... 49% 48.3M 1s
    ##  39350K .......... .......... .......... .......... .......... 49% 96.8M 1s
    ##  39400K .......... .......... .......... .......... .......... 49% 19.9M 1s
    ##  39450K .......... .......... .......... .......... .......... 49%  105M 1s
    ##  39500K .......... .......... .......... .......... .......... 49% 96.8M 1s
    ##  39550K .......... .......... .......... .......... .......... 49% 20.8M 1s
    ##  39600K .......... .......... .......... .......... .......... 49% 91.4M 1s
    ##  39650K .......... .......... .......... .......... .......... 49% 84.8M 1s
    ##  39700K .......... .......... .......... .......... .......... 49% 81.3M 1s
    ##  39750K .......... .......... .......... .......... .......... 49% 56.4M 1s
    ##  39800K .......... .......... .......... .......... .......... 49% 25.9M 1s
    ##  39850K .......... .......... .......... .......... .......... 49% 64.8M 1s
    ##  39900K .......... .......... .......... .......... .......... 49% 65.4M 1s
    ##  39950K .......... .......... .......... .......... .......... 50% 79.4M 1s
    ##  40000K .......... .......... .......... .......... .......... 50% 89.6M 1s
    ##  40050K .......... .......... .......... .......... .......... 50% 11.0M 1s
    ##  40100K .......... .......... .......... .......... .......... 50% 58.6M 1s
    ##  40150K .......... .......... .......... .......... .......... 50% 75.9M 1s
    ##  40200K .......... .......... .......... .......... .......... 50% 87.4M 1s
    ##  40250K .......... .......... .......... .......... .......... 50%  101M 1s
    ##  40300K .......... .......... .......... .......... .......... 50% 91.4M 1s
    ##  40350K .......... .......... .......... .......... .......... 50% 49.9M 1s
    ##  40400K .......... .......... .......... .......... .......... 50% 77.8M 1s
    ##  40450K .......... .......... .......... .......... .......... 50% 49.5M 1s
    ##  40500K .......... .......... .......... .......... .......... 50% 38.6M 1s
    ##  40550K .......... .......... .......... .......... .......... 50% 50.3M 1s
    ##  40600K .......... .......... .......... .......... .......... 50% 85.7M 1s
    ##  40650K .......... .......... .......... .......... .......... 50% 82.5M 1s
    ##  40700K .......... .......... .......... .......... .......... 50% 78.0M 1s
    ##  40750K .......... .......... .......... .......... .......... 51% 77.5M 1s
    ##  40800K .......... .......... .......... .......... .......... 51% 53.5M 1s
    ##  40850K .......... .......... .......... .......... .......... 51% 46.6M 1s
    ##  40900K .......... .......... .......... .......... .......... 51% 36.4M 1s
    ##  40950K .......... .......... .......... .......... .......... 51% 57.9M 1s
    ##  41000K .......... .......... .......... .......... .......... 51% 89.2M 1s
    ##  41050K .......... .......... .......... .......... .......... 51% 79.5M 1s
    ##  41100K .......... .......... .......... .......... .......... 51% 84.1M 1s
    ##  41150K .......... .......... .......... .......... .......... 51% 73.3M 1s
    ##  41200K .......... .......... .......... .......... .......... 51% 76.1M 1s
    ##  41250K .......... .......... .......... .......... .......... 51% 54.7M 1s
    ##  41300K .......... .......... .......... .......... .......... 51% 51.8M 1s
    ##  41350K .......... .......... .......... .......... .......... 51% 69.8M 1s
    ##  41400K .......... .......... .......... .......... .......... 51% 69.0M 1s
    ##  41450K .......... .......... .......... .......... .......... 51% 74.9M 1s
    ##  41500K .......... .......... .......... .......... .......... 51% 63.7M 1s
    ##  41550K .......... .......... .......... .......... .......... 52%  100M 1s
    ##  41600K .......... .......... .......... .......... .......... 52% 87.4M 1s
    ##  41650K .......... .......... .......... .......... .......... 52% 70.4M 1s
    ##  41700K .......... .......... .......... .......... .......... 52% 84.2M 1s
    ##  41750K .......... .......... .......... .......... .......... 52% 72.5M 1s
    ##  41800K .......... .......... .......... .......... .......... 52% 34.0M 1s
    ##  41850K .......... .......... .......... .......... .......... 52% 65.4M 1s
    ##  41900K .......... .......... .......... .......... .......... 52% 80.9M 1s
    ##  41950K .......... .......... .......... .......... .......... 52% 71.0M 1s
    ##  42000K .......... .......... .......... .......... .......... 52% 72.1M 1s
    ##  42050K .......... .......... .......... .......... .......... 52% 53.1M 1s
    ##  42100K .......... .......... .......... .......... .......... 52% 56.8M 1s
    ##  42150K .......... .......... .......... .......... .......... 52% 69.8M 1s
    ##  42200K .......... .......... .......... .......... .......... 52% 69.4M 1s
    ##  42250K .......... .......... .......... .......... .......... 52% 70.2M 1s
    ##  42300K .......... .......... .......... .......... .......... 52% 55.3M 1s
    ##  42350K .......... .......... .......... .......... .......... 53% 46.3M 1s
    ##  42400K .......... .......... .......... .......... .......... 53% 79.5M 1s
    ##  42450K .......... .......... .......... .......... .......... 53% 81.1M 1s
    ##  42500K .......... .......... .......... .......... .......... 53% 68.4M 1s
    ##  42550K .......... .......... .......... .......... .......... 53% 63.4M 1s
    ##  42600K .......... .......... .......... .......... .......... 53% 61.2M 1s
    ##  42650K .......... .......... .......... .......... .......... 53% 87.4M 1s
    ##  42700K .......... .......... .......... .......... .......... 53% 64.6M 1s
    ##  42750K .......... .......... .......... .......... .......... 53% 69.4M 1s
    ##  42800K .......... .......... .......... .......... .......... 53% 79.1M 1s
    ##  42850K .......... .......... .......... .......... .......... 53% 52.4M 1s
    ##  42900K .......... .......... .......... .......... .......... 53% 81.3M 1s
    ##  42950K .......... .......... .......... .......... .......... 53% 72.4M 1s
    ##  43000K .......... .......... .......... .......... .......... 53% 57.5M 1s
    ##  43050K .......... .......... .......... .......... .......... 53% 72.1M 1s
    ##  43100K .......... .......... .......... .......... .......... 53% 65.2M 1s
    ##  43150K .......... .......... .......... .......... .......... 54% 85.9M 1s
    ##  43200K .......... .......... .......... .......... .......... 54% 61.9M 1s
    ##  43250K .......... .......... .......... .......... .......... 54% 68.2M 1s
    ##  43300K .......... .......... .......... .......... .......... 54% 61.6M 1s
    ##  43350K .......... .......... .......... .......... .......... 54% 71.2M 1s
    ##  43400K .......... .......... .......... .......... .......... 54% 88.5M 1s
    ##  43450K .......... .......... .......... .......... .......... 54% 80.0M 1s
    ##  43500K .......... .......... .......... .......... .......... 54% 38.0M 1s
    ##  43550K .......... .......... .......... .......... .......... 54% 82.7M 1s
    ##  43600K .......... .......... .......... .......... .......... 54% 95.0M 1s
    ##  43650K .......... .......... .......... .......... .......... 54% 86.3M 1s
    ##  43700K .......... .......... .......... .......... .......... 54% 56.2M 1s
    ##  43750K .......... .......... .......... .......... .......... 54% 46.3M 1s
    ##  43800K .......... .......... .......... .......... .......... 54% 63.0M 1s
    ##  43850K .......... .......... .......... .......... .......... 54%  115M 1s
    ##  43900K .......... .......... .......... .......... .......... 54% 48.1M 1s
    ##  43950K .......... .......... .......... .......... .......... 55% 82.4M 1s
    ##  44000K .......... .......... .......... .......... .......... 55% 51.6M 1s
    ##  44050K .......... .......... .......... .......... .......... 55% 66.4M 1s
    ##  44100K .......... .......... .......... .......... .......... 55% 82.3M 1s
    ##  44150K .......... .......... .......... .......... .......... 55% 40.1M 1s
    ##  44200K .......... .......... .......... .......... .......... 55% 44.1M 1s
    ##  44250K .......... .......... .......... .......... .......... 55% 70.4M 1s
    ##  44300K .......... .......... .......... .......... .......... 55% 80.6M 1s
    ##  44350K .......... .......... .......... .......... .......... 55% 96.7M 1s
    ##  44400K .......... .......... .......... .......... .......... 55% 49.8M 1s
    ##  44450K .......... .......... .......... .......... .......... 55% 63.9M 1s
    ##  44500K .......... .......... .......... .......... .......... 55% 56.4M 1s
    ##  44550K .......... .......... .......... .......... .......... 55% 44.1M 1s
    ##  44600K .......... .......... .......... .......... .......... 55% 58.1M 1s
    ##  44650K .......... .......... .......... .......... .......... 55% 83.2M 1s
    ##  44700K .......... .......... .......... .......... .......... 55% 59.7M 1s
    ##  44750K .......... .......... .......... .......... .......... 56%  107M 1s
    ##  44800K .......... .......... .......... .......... .......... 56% 82.7M 1s
    ##  44850K .......... .......... .......... .......... .......... 56% 56.8M 1s
    ##  44900K .......... .......... .......... .......... .......... 56% 50.6M 1s
    ##  44950K .......... .......... .......... .......... .......... 56% 66.2M 1s
    ##  45000K .......... .......... .......... .......... .......... 56% 82.0M 1s
    ##  45050K .......... .......... .......... .......... .......... 56% 41.5M 1s
    ##  45100K .......... .......... .......... .......... .......... 56% 66.9M 1s
    ##  45150K .......... .......... .......... .......... .......... 56% 86.3M 1s
    ##  45200K .......... .......... .......... .......... .......... 56% 60.3M 1s
    ##  45250K .......... .......... .......... .......... .......... 56% 81.1M 1s
    ##  45300K .......... .......... .......... .......... .......... 56% 83.3M 1s
    ##  45350K .......... .......... .......... .......... .......... 56% 62.1M 1s
    ##  45400K .......... .......... .......... .......... .......... 56% 53.3M 1s
    ##  45450K .......... .......... .......... .......... .......... 56% 96.9M 1s
    ##  45500K .......... .......... .......... .......... .......... 56% 61.3M 1s
    ##  45550K .......... .......... .......... .......... .......... 57% 50.0M 1s
    ##  45600K .......... .......... .......... .......... .......... 57% 87.4M 1s
    ##  45650K .......... .......... .......... .......... .......... 57% 89.1M 1s
    ##  45700K .......... .......... .......... .......... .......... 57% 55.4M 1s
    ##  45750K .......... .......... .......... .......... .......... 57%  108M 1s
    ##  45800K .......... .......... .......... .......... .......... 57% 55.3M 1s
    ##  45850K .......... .......... .......... .......... .......... 57% 29.9M 1s
    ##  45900K .......... .......... .......... .......... .......... 57% 87.4M 1s
    ##  45950K .......... .......... .......... .......... .......... 57%  135M 1s
    ##  46000K .......... .......... .......... .......... .......... 57% 62.6M 1s
    ##  46050K .......... .......... .......... .......... .......... 57% 98.5M 1s
    ##  46100K .......... .......... .......... .......... .......... 57% 86.8M 1s
    ##  46150K .......... .......... .......... .......... .......... 57% 76.4M 1s
    ##  46200K .......... .......... .......... .......... .......... 57% 86.5M 1s
    ##  46250K .......... .......... .......... .......... .......... 57% 52.5M 1s
    ##  46300K .......... .......... .......... .......... .......... 57% 44.1M 1s
    ##  46350K .......... .......... .......... .......... .......... 58% 75.3M 1s
    ##  46400K .......... .......... .......... .......... .......... 58% 35.4M 1s
    ##  46450K .......... .......... .......... .......... .......... 58%  100M 1s
    ##  46500K .......... .......... .......... .......... .......... 58% 94.4M 1s
    ##  46550K .......... .......... .......... .......... .......... 58%  112M 1s
    ##  46600K .......... .......... .......... .......... .......... 58% 69.6M 1s
    ##  46650K .......... .......... .......... .......... .......... 58% 53.0M 1s
    ##  46700K .......... .......... .......... .......... .......... 58% 82.5M 1s
    ##  46750K .......... .......... .......... .......... .......... 58% 67.3M 1s
    ##  46800K .......... .......... .......... .......... .......... 58% 48.6M 1s
    ##  46850K .......... .......... .......... .......... .......... 58% 39.3M 1s
    ##  46900K .......... .......... .......... .......... .......... 58% 56.6M 1s
    ##  46950K .......... .......... .......... .......... .......... 58% 47.3M 1s
    ##  47000K .......... .......... .......... .......... .......... 58% 83.3M 1s
    ##  47050K .......... .......... .......... .......... .......... 58%  105M 1s
    ##  47100K .......... .......... .......... .......... .......... 58%  106M 1s
    ##  47150K .......... .......... .......... .......... .......... 59%  119M 1s
    ##  47200K .......... .......... .......... .......... .......... 59% 91.0M 1s
    ##  47250K .......... .......... .......... .......... .......... 59% 98.5M 1s
    ##  47300K .......... .......... .......... .......... .......... 59% 39.0M 1s
    ##  47350K .......... .......... .......... .......... .......... 59% 58.8M 1s
    ##  47400K .......... .......... .......... .......... .......... 59% 13.1M 1s
    ##  47450K .......... .......... .......... .......... .......... 59% 53.3M 1s
    ##  47500K .......... .......... .......... .......... .......... 59%  113M 1s
    ##  47550K .......... .......... .......... .......... .......... 59%  113M 1s
    ##  47600K .......... .......... .......... .......... .......... 59%  115M 1s
    ##  47650K .......... .......... .......... .......... .......... 59%  111M 1s
    ##  47700K .......... .......... .......... .......... .......... 59%  112M 1s
    ##  47750K .......... .......... .......... .......... .......... 59% 99.4M 1s
    ##  47800K .......... .......... .......... .......... .......... 59% 90.1M 1s
    ##  47850K .......... .......... .......... .......... .......... 59% 32.1M 1s
    ##  47900K .......... .......... .......... .......... .......... 59% 38.2M 1s
    ##  47950K .......... .......... .......... .......... .......... 60% 57.9M 1s
    ##  48000K .......... .......... .......... .......... .......... 60% 75.2M 1s
    ##  48050K .......... .......... .......... .......... .......... 60% 82.1M 1s
    ##  48100K .......... .......... .......... .......... .......... 60% 72.8M 1s
    ##  48150K .......... .......... .......... .......... .......... 60% 72.6M 1s
    ##  48200K .......... .......... .......... .......... .......... 60% 35.5M 1s
    ##  48250K .......... .......... .......... .......... .......... 60%  115M 1s
    ##  48300K .......... .......... .......... .......... .......... 60% 35.3M 1s
    ##  48350K .......... .......... .......... .......... .......... 60% 83.6M 1s
    ##  48400K .......... .......... .......... .......... .......... 60% 71.6M 1s
    ##  48450K .......... .......... .......... .......... .......... 60% 52.6M 1s
    ##  48500K .......... .......... .......... .......... .......... 60% 67.0M 1s
    ##  48550K .......... .......... .......... .......... .......... 60%  117M 1s
    ##  48600K .......... .......... .......... .......... .......... 60% 67.8M 1s
    ##  48650K .......... .......... .......... .......... .......... 60% 43.1M 1s
    ##  48700K .......... .......... .......... .......... .......... 60% 76.5M 1s
    ##  48750K .......... .......... .......... .......... .......... 61% 50.1M 1s
    ##  48800K .......... .......... .......... .......... .......... 61% 56.1M 1s
    ##  48850K .......... .......... .......... .......... .......... 61% 98.4M 1s
    ##  48900K .......... .......... .......... .......... .......... 61% 58.9M 1s
    ##  48950K .......... .......... .......... .......... .......... 61% 56.4M 1s
    ##  49000K .......... .......... .......... .......... .......... 61% 45.9M 1s
    ##  49050K .......... .......... .......... .......... .......... 61%  107M 1s
    ##  49100K .......... .......... .......... .......... .......... 61% 53.5M 1s
    ##  49150K .......... .......... .......... .......... .......... 61%  111M 1s
    ##  49200K .......... .......... .......... .......... .......... 61% 27.7M 1s
    ##  49250K .......... .......... .......... .......... .......... 61%  125M 1s
    ##  49300K .......... .......... .......... .......... .......... 61% 16.7M 1s
    ##  49350K .......... .......... .......... .......... .......... 61%  100M 1s
    ##  49400K .......... .......... .......... .......... .......... 61%  119M 1s
    ##  49450K .......... .......... .......... .......... .......... 61%  137M 1s
    ##  49500K .......... .......... .......... .......... .......... 61% 88.3M 1s
    ##  49550K .......... .......... .......... .......... .......... 62%  113M 1s
    ##  49600K .......... .......... .......... .......... .......... 62% 45.4M 1s
    ##  49650K .......... .......... .......... .......... .......... 62% 44.2M 1s
    ##  49700K .......... .......... .......... .......... .......... 62% 55.8M 1s
    ##  49750K .......... .......... .......... .......... .......... 62% 37.1M 1s
    ##  49800K .......... .......... .......... .......... .......... 62% 62.7M 1s
    ##  49850K .......... .......... .......... .......... .......... 62%  119M 1s
    ##  49900K .......... .......... .......... .......... .......... 62%  111M 1s
    ##  49950K .......... .......... .......... .......... .......... 62%  101M 1s
    ##  50000K .......... .......... .......... .......... .......... 62%  130M 1s
    ##  50050K .......... .......... .......... .......... .......... 62% 94.0M 1s
    ##  50100K .......... .......... .......... .......... .......... 62%  102M 1s
    ##  50150K .......... .......... .......... .......... .......... 62% 69.0M 1s
    ##  50200K .......... .......... .......... .......... .......... 62% 46.1M 1s
    ##  50250K .......... .......... .......... .......... .......... 62%  104M 1s
    ##  50300K .......... .......... .......... .......... .......... 62% 88.9M 1s
    ##  50350K .......... .......... .......... .......... .......... 63% 78.7M 1s
    ##  50400K .......... .......... .......... .......... .......... 63%  112M 1s
    ##  50450K .......... .......... .......... .......... .......... 63%  127M 1s
    ##  50500K .......... .......... .......... .......... .......... 63% 34.8M 1s
    ##  50550K .......... .......... .......... .......... .......... 63% 31.8M 1s
    ##  50600K .......... .......... .......... .......... .......... 63%  118M 1s
    ##  50650K .......... .......... .......... .......... .......... 63% 86.8M 1s
    ##  50700K .......... .......... .......... .......... .......... 63%  134M 1s
    ##  50750K .......... .......... .......... .......... .......... 63%  111M 1s
    ##  50800K .......... .......... .......... .......... .......... 63% 81.8M 1s
    ##  50850K .......... .......... .......... .......... .......... 63% 29.9M 1s
    ##  50900K .......... .......... .......... .......... .......... 63% 79.5M 1s
    ##  50950K .......... .......... .......... .......... .......... 63%  104M 1s
    ##  51000K .......... .......... .......... .......... .......... 63% 46.9M 1s
    ##  51050K .......... .......... .......... .......... .......... 63% 84.7M 1s
    ##  51100K .......... .......... .......... .......... .......... 63% 79.3M 1s
    ##  51150K .......... .......... .......... .......... .......... 64%  112M 1s
    ##  51200K .......... .......... .......... .......... .......... 64%  121M 1s
    ##  51250K .......... .......... .......... .......... .......... 64% 93.3M 1s
    ##  51300K .......... .......... .......... .......... .......... 64% 39.9M 1s
    ##  51350K .......... .......... .......... .......... .......... 64% 59.3M 1s
    ##  51400K .......... .......... .......... .......... .......... 64% 85.1M 1s
    ##  51450K .......... .......... .......... .......... .......... 64% 83.0M 1s
    ##  51500K .......... .......... .......... .......... .......... 64% 77.0M 1s
    ##  51550K .......... .......... .......... .......... .......... 64% 72.4M 1s
    ##  51600K .......... .......... .......... .......... .......... 64% 74.4M 1s
    ##  51650K .......... .......... .......... .......... .......... 64%  120M 1s
    ##  51700K .......... .......... .......... .......... .......... 64% 73.0M 1s
    ##  51750K .......... .......... .......... .......... .......... 64%  119M 1s
    ##  51800K .......... .......... .......... .......... .......... 64% 85.2M 1s
    ##  51850K .......... .......... .......... .......... .......... 64% 63.0M 1s
    ##  51900K .......... .......... .......... .......... .......... 65% 93.4M 1s
    ##  51950K .......... .......... .......... .......... .......... 65% 73.5M 1s
    ##  52000K .......... .......... .......... .......... .......... 65% 68.6M 1s
    ##  52050K .......... .......... .......... .......... .......... 65% 91.4M 1s
    ##  52100K .......... .......... .......... .......... .......... 65%  118M 1s
    ##  52150K .......... .......... .......... .......... .......... 65% 69.7M 1s
    ##  52200K .......... .......... .......... .......... .......... 65% 66.7M 1s
    ##  52250K .......... .......... .......... .......... .......... 65% 89.3M 1s
    ##  52300K .......... .......... .......... .......... .......... 65% 83.5M 1s
    ##  52350K .......... .......... .......... .......... .......... 65% 75.6M 1s
    ##  52400K .......... .......... .......... .......... .......... 65% 27.4M 1s
    ##  52450K .......... .......... .......... .......... .......... 65% 38.8M 1s
    ##  52500K .......... .......... .......... .......... .......... 65% 30.9M 1s
    ##  52550K .......... .......... .......... .......... .......... 65% 30.4M 1s
    ##  52600K .......... .......... .......... .......... .......... 65% 46.9M 1s
    ##  52650K .......... .......... .......... .......... .......... 65% 69.2M 1s
    ##  52700K .......... .......... .......... .......... .......... 66% 64.6M 1s
    ##  52750K .......... .......... .......... .......... .......... 66% 88.6M 1s
    ##  52800K .......... .......... .......... .......... .......... 66% 57.9M 1s
    ##  52850K .......... .......... .......... .......... .......... 66% 71.4M 1s
    ##  52900K .......... .......... .......... .......... .......... 66% 51.8M 1s
    ##  52950K .......... .......... .......... .......... .......... 66% 69.8M 1s
    ##  53000K .......... .......... .......... .......... .......... 66% 60.5M 1s
    ##  53050K .......... .......... .......... .......... .......... 66% 80.7M 1s
    ##  53100K .......... .......... .......... .......... .......... 66% 68.3M 1s
    ##  53150K .......... .......... .......... .......... .......... 66% 60.9M 1s
    ##  53200K .......... .......... .......... .......... .......... 66% 82.0M 1s
    ##  53250K .......... .......... .......... .......... .......... 66% 77.3M 1s
    ##  53300K .......... .......... .......... .......... .......... 66% 46.5M 1s
    ##  53350K .......... .......... .......... .......... .......... 66% 73.6M 1s
    ##  53400K .......... .......... .......... .......... .......... 66% 81.4M 1s
    ##  53450K .......... .......... .......... .......... .......... 66% 45.9M 1s
    ##  53500K .......... .......... .......... .......... .......... 67% 97.2M 1s
    ##  53550K .......... .......... .......... .......... .......... 67% 44.8M 1s
    ##  53600K .......... .......... .......... .......... .......... 67% 62.9M 1s
    ##  53650K .......... .......... .......... .......... .......... 67% 76.9M 1s
    ##  53700K .......... .......... .......... .......... .......... 67% 47.0M 1s
    ##  53750K .......... .......... .......... .......... .......... 67% 63.5M 1s
    ##  53800K .......... .......... .......... .......... .......... 67% 69.2M 1s
    ##  53850K .......... .......... .......... .......... .......... 67% 44.8M 1s
    ##  53900K .......... .......... .......... .......... .......... 67% 68.3M 1s
    ##  53950K .......... .......... .......... .......... .......... 67% 71.5M 1s
    ##  54000K .......... .......... .......... .......... .......... 67% 95.0M 1s
    ##  54050K .......... .......... .......... .......... .......... 67% 5.93M 1s
    ##  54100K .......... .......... .......... .......... .......... 67%  101M 1s
    ##  54150K .......... .......... .......... .......... .......... 67%  106M 1s
    ##  54200K .......... .......... .......... .......... .......... 67% 88.4M 1s
    ##  54250K .......... .......... .......... .......... .......... 67%  113M 1s
    ##  54300K .......... .......... .......... .......... .......... 68% 89.8M 1s
    ##  54350K .......... .......... .......... .......... .......... 68% 82.7M 1s
    ##  54400K .......... .......... .......... .......... .......... 68% 59.2M 1s
    ##  54450K .......... .......... .......... .......... .......... 68% 64.8M 1s
    ##  54500K .......... .......... .......... .......... .......... 68% 41.3M 1s
    ##  54550K .......... .......... .......... .......... .......... 68% 87.3M 1s
    ##  54600K .......... .......... .......... .......... .......... 68% 94.5M 0s
    ##  54650K .......... .......... .......... .......... .......... 68% 8.27M 1s
    ##  54700K .......... .......... .......... .......... .......... 68% 77.0M 0s
    ##  54750K .......... .......... .......... .......... .......... 68% 29.1M 0s
    ##  54800K .......... .......... .......... .......... .......... 68% 35.7M 0s
    ##  54850K .......... .......... .......... .......... .......... 68% 75.2M 0s
    ##  54900K .......... .......... .......... .......... .......... 68% 88.2M 0s
    ##  54950K .......... .......... .......... .......... .......... 68% 46.2M 0s
    ##  55000K .......... .......... .......... .......... .......... 68%  109M 0s
    ##  55050K .......... .......... .......... .......... .......... 68% 99.5M 0s
    ##  55100K .......... .......... .......... .......... .......... 69% 30.8M 0s
    ##  55150K .......... .......... .......... .......... .......... 69% 31.1M 0s
    ##  55200K .......... .......... .......... .......... .......... 69% 32.1M 0s
    ##  55250K .......... .......... .......... .......... .......... 69% 71.3M 0s
    ##  55300K .......... .......... .......... .......... .......... 69% 66.2M 0s
    ##  55350K .......... .......... .......... .......... .......... 69% 85.5M 0s
    ##  55400K .......... .......... .......... .......... .......... 69%  104M 0s
    ##  55450K .......... .......... .......... .......... .......... 69% 81.0M 0s
    ##  55500K .......... .......... .......... .......... .......... 69%  114M 0s
    ##  55550K .......... .......... .......... .......... .......... 69% 96.6M 0s
    ##  55600K .......... .......... .......... .......... .......... 69% 40.5M 0s
    ##  55650K .......... .......... .......... .......... .......... 69% 87.6M 0s
    ##  55700K .......... .......... .......... .......... .......... 69% 95.4M 0s
    ##  55750K .......... .......... .......... .......... .......... 69% 35.9M 0s
    ##  55800K .......... .......... .......... .......... .......... 69%  109M 0s
    ##  55850K .......... .......... .......... .......... .......... 69% 57.8M 0s
    ##  55900K .......... .......... .......... .......... .......... 70% 78.8M 0s
    ##  55950K .......... .......... .......... .......... .......... 70% 37.1M 0s
    ##  56000K .......... .......... .......... .......... .......... 70% 78.3M 0s
    ##  56050K .......... .......... .......... .......... .......... 70%  119M 0s
    ##  56100K .......... .......... .......... .......... .......... 70% 32.9M 0s
    ##  56150K .......... .......... .......... .......... .......... 70% 50.3M 0s
    ##  56200K .......... .......... .......... .......... .......... 70% 37.6M 0s
    ##  56250K .......... .......... .......... .......... .......... 70% 72.5M 0s
    ##  56300K .......... .......... .......... .......... .......... 70% 54.5M 0s
    ##  56350K .......... .......... .......... .......... .......... 70% 96.0M 0s
    ##  56400K .......... .......... .......... .......... .......... 70% 64.9M 0s
    ##  56450K .......... .......... .......... .......... .......... 70% 88.0M 0s
    ##  56500K .......... .......... .......... .......... .......... 70% 85.0M 0s
    ##  56550K .......... .......... .......... .......... .......... 70% 38.1M 0s
    ##  56600K .......... .......... .......... .......... .......... 70% 67.1M 0s
    ##  56650K .......... .......... .......... .......... .......... 70% 23.9M 0s
    ##  56700K .......... .......... .......... .......... .......... 71% 68.8M 0s
    ##  56750K .......... .......... .......... .......... .......... 71%  109M 0s
    ##  56800K .......... .......... .......... .......... .......... 71% 65.2M 0s
    ##  56850K .......... .......... .......... .......... .......... 71% 98.7M 0s
    ##  56900K .......... .......... .......... .......... .......... 71% 11.6M 0s
    ##  56950K .......... .......... .......... .......... .......... 71% 77.4M 0s
    ##  57000K .......... .......... .......... .......... .......... 71% 37.2M 0s
    ##  57050K .......... .......... .......... .......... .......... 71% 79.1M 0s
    ##  57100K .......... .......... .......... .......... .......... 71% 69.5M 0s
    ##  57150K .......... .......... .......... .......... .......... 71%  100M 0s
    ##  57200K .......... .......... .......... .......... .......... 71% 19.4M 0s
    ##  57250K .......... .......... .......... .......... .......... 71% 77.9M 0s
    ##  57300K .......... .......... .......... .......... .......... 71% 42.1M 0s
    ##  57350K .......... .......... .......... .......... .......... 71% 77.5M 0s
    ##  57400K .......... .......... .......... .......... .......... 71% 70.8M 0s
    ##  57450K .......... .......... .......... .......... .......... 71%  116M 0s
    ##  57500K .......... .......... .......... .......... .......... 72% 51.0M 0s
    ##  57550K .......... .......... .......... .......... .......... 72% 61.8M 0s
    ##  57600K .......... .......... .......... .......... .......... 72% 92.1M 0s
    ##  57650K .......... .......... .......... .......... .......... 72% 26.3M 0s
    ##  57700K .......... .......... .......... .......... .......... 72% 59.8M 0s
    ##  57750K .......... .......... .......... .......... .......... 72% 85.3M 0s
    ##  57800K .......... .......... .......... .......... .......... 72% 69.5M 0s
    ##  57850K .......... .......... .......... .......... .......... 72%  109M 0s
    ##  57900K .......... .......... .......... .......... .......... 72% 74.5M 0s
    ##  57950K .......... .......... .......... .......... .......... 72%  111M 0s
    ##  58000K .......... .......... .......... .......... .......... 72% 46.9M 0s
    ##  58050K .......... .......... .......... .......... .......... 72% 82.4M 0s
    ##  58100K .......... .......... .......... .......... .......... 72% 32.0M 0s
    ##  58150K .......... .......... .......... .......... .......... 72% 68.8M 0s
    ##  58200K .......... .......... .......... .......... .......... 72% 58.8M 0s
    ##  58250K .......... .......... .......... .......... .......... 72% 78.8M 0s
    ##  58300K .......... .......... .......... .......... .......... 73% 97.0M 0s
    ##  58350K .......... .......... .......... .......... .......... 73%  112M 0s
    ##  58400K .......... .......... .......... .......... .......... 73% 97.6M 0s
    ##  58450K .......... .......... .......... .......... .......... 73% 87.3M 0s
    ##  58500K .......... .......... .......... .......... .......... 73% 30.5M 0s
    ##  58550K .......... .......... .......... .......... .......... 73% 93.3M 0s
    ##  58600K .......... .......... .......... .......... .......... 73% 52.3M 0s
    ##  58650K .......... .......... .......... .......... .......... 73% 91.2M 0s
    ##  58700K .......... .......... .......... .......... .......... 73% 68.3M 0s
    ##  58750K .......... .......... .......... .......... .......... 73% 81.2M 0s
    ##  58800K .......... .......... .......... .......... .......... 73% 75.4M 0s
    ##  58850K .......... .......... .......... .......... .......... 73% 76.1M 0s
    ##  58900K .......... .......... .......... .......... .......... 73% 28.1M 0s
    ##  58950K .......... .......... .......... .......... .......... 73% 64.2M 0s
    ##  59000K .......... .......... .......... .......... .......... 73% 27.7M 0s
    ##  59050K .......... .......... .......... .......... .......... 73% 28.3M 0s
    ##  59100K .......... .......... .......... .......... .......... 74% 82.3M 0s
    ##  59150K .......... .......... .......... .......... .......... 74% 99.1M 0s
    ##  59200K .......... .......... .......... .......... .......... 74% 76.2M 0s
    ##  59250K .......... .......... .......... .......... .......... 74% 66.9M 0s
    ##  59300K .......... .......... .......... .......... .......... 74% 47.9M 0s
    ##  59350K .......... .......... .......... .......... .......... 74% 78.1M 0s
    ##  59400K .......... .......... .......... .......... .......... 74% 75.7M 0s
    ##  59450K .......... .......... .......... .......... .......... 74% 92.2M 0s
    ##  59500K .......... .......... .......... .......... .......... 74% 84.3M 0s
    ##  59550K .......... .......... .......... .......... .......... 74% 88.2M 0s
    ##  59600K .......... .......... .......... .......... .......... 74% 15.3M 0s
    ##  59650K .......... .......... .......... .......... .......... 74% 74.6M 0s
    ##  59700K .......... .......... .......... .......... .......... 74% 76.2M 0s
    ##  59750K .......... .......... .......... .......... .......... 74% 84.5M 0s
    ##  59800K .......... .......... .......... .......... .......... 74% 89.6M 0s
    ##  59850K .......... .......... .......... .......... .......... 74%  103M 0s
    ##  59900K .......... .......... .......... .......... .......... 75% 79.2M 0s
    ##  59950K .......... .......... .......... .......... .......... 75% 40.0M 0s
    ##  60000K .......... .......... .......... .......... .......... 75% 72.0M 0s
    ##  60050K .......... .......... .......... .......... .......... 75% 89.0M 0s
    ##  60100K .......... .......... .......... .......... .......... 75% 63.3M 0s
    ##  60150K .......... .......... .......... .......... .......... 75% 77.8M 0s
    ##  60200K .......... .......... .......... .......... .......... 75% 93.3M 0s
    ##  60250K .......... .......... .......... .......... .......... 75% 86.2M 0s
    ##  60300K .......... .......... .......... .......... .......... 75% 90.4M 0s
    ##  60350K .......... .......... .......... .......... .......... 75% 80.4M 0s
    ##  60400K .......... .......... .......... .......... .......... 75% 76.0M 0s
    ##  60450K .......... .......... .......... .......... .......... 75% 33.6M 0s
    ##  60500K .......... .......... .......... .......... .......... 75% 90.4M 0s
    ##  60550K .......... .......... .......... .......... .......... 75% 61.8M 0s
    ##  60600K .......... .......... .......... .......... .......... 75% 84.1M 0s
    ##  60650K .......... .......... .......... .......... .......... 75% 78.3M 0s
    ##  60700K .......... .......... .......... .......... .......... 76% 43.9M 0s
    ##  60750K .......... .......... .......... .......... .......... 76%  101M 0s
    ##  60800K .......... .......... .......... .......... .......... 76% 61.4M 0s
    ##  60850K .......... .......... .......... .......... .......... 76% 84.3M 0s
    ##  60900K .......... .......... .......... .......... .......... 76% 12.3M 0s
    ##  60950K .......... .......... .......... .......... .......... 76%  102M 0s
    ##  61000K .......... .......... .......... .......... .......... 76% 88.8M 0s
    ##  61050K .......... .......... .......... .......... .......... 76%  108M 0s
    ##  61100K .......... .......... .......... .......... .......... 76% 94.4M 0s
    ##  61150K .......... .......... .......... .......... .......... 76% 51.1M 0s
    ##  61200K .......... .......... .......... .......... .......... 76% 27.3M 0s
    ##  61250K .......... .......... .......... .......... .......... 76% 93.5M 0s
    ##  61300K .......... .......... .......... .......... .......... 76% 40.7M 0s
    ##  61350K .......... .......... .......... .......... .......... 76% 84.2M 0s
    ##  61400K .......... .......... .......... .......... .......... 76% 87.9M 0s
    ##  61450K .......... .......... .......... .......... .......... 76% 83.7M 0s
    ##  61500K .......... .......... .......... .......... .......... 77% 90.8M 0s
    ##  61550K .......... .......... .......... .......... .......... 77%  102M 0s
    ##  61600K .......... .......... .......... .......... .......... 77% 52.8M 0s
    ##  61650K .......... .......... .......... .......... .......... 77% 69.3M 0s
    ##  61700K .......... .......... .......... .......... .......... 77% 93.5M 0s
    ##  61750K .......... .......... .......... .......... .......... 77% 97.6M 0s
    ##  61800K .......... .......... .......... .......... .......... 77% 58.8M 0s
    ##  61850K .......... .......... .......... .......... .......... 77% 96.3M 0s
    ##  61900K .......... .......... .......... .......... .......... 77% 34.6M 0s
    ##  61950K .......... .......... .......... .......... .......... 77% 68.3M 0s
    ##  62000K .......... .......... .......... .......... .......... 77% 79.2M 0s
    ##  62050K .......... .......... .......... .......... .......... 77% 86.3M 0s
    ##  62100K .......... .......... .......... .......... .......... 77% 95.8M 0s
    ##  62150K .......... .......... .......... .......... .......... 77%  105M 0s
    ##  62200K .......... .......... .......... .......... .......... 77% 82.0M 0s
    ##  62250K .......... .......... .......... .......... .......... 77% 76.2M 0s
    ##  62300K .......... .......... .......... .......... .......... 78% 76.4M 0s
    ##  62350K .......... .......... .......... .......... .......... 78% 65.4M 0s
    ##  62400K .......... .......... .......... .......... .......... 78% 85.2M 0s
    ##  62450K .......... .......... .......... .......... .......... 78% 99.1M 0s
    ##  62500K .......... .......... .......... .......... .......... 78% 84.7M 0s
    ##  62550K .......... .......... .......... .......... .......... 78% 55.8M 0s
    ##  62600K .......... .......... .......... .......... .......... 78%  100M 0s
    ##  62650K .......... .......... .......... .......... .......... 78% 87.7M 0s
    ##  62700K .......... .......... .......... .......... .......... 78% 72.0M 0s
    ##  62750K .......... .......... .......... .......... .......... 78%  101M 0s
    ##  62800K .......... .......... .......... .......... .......... 78% 94.2M 0s
    ##  62850K .......... .......... .......... .......... .......... 78% 62.6M 0s
    ##  62900K .......... .......... .......... .......... .......... 78% 87.3M 0s
    ##  62950K .......... .......... .......... .......... .......... 78% 91.7M 0s
    ##  63000K .......... .......... .......... .......... .......... 78% 56.0M 0s
    ##  63050K .......... .......... .......... .......... .......... 78% 57.7M 0s
    ##  63100K .......... .......... .......... .......... .......... 79% 91.1M 0s
    ##  63150K .......... .......... .......... .......... .......... 79% 90.7M 0s
    ##  63200K .......... .......... .......... .......... .......... 79% 80.7M 0s
    ##  63250K .......... .......... .......... .......... .......... 79% 94.9M 0s
    ##  63300K .......... .......... .......... .......... .......... 79% 86.5M 0s
    ##  63350K .......... .......... .......... .......... .......... 79% 68.1M 0s
    ##  63400K .......... .......... .......... .......... .......... 79% 76.5M 0s
    ##  63450K .......... .......... .......... .......... .......... 79%  105M 0s
    ##  63500K .......... .......... .......... .......... .......... 79% 74.3M 0s
    ##  63550K .......... .......... .......... .......... .......... 79% 84.8M 0s
    ##  63600K .......... .......... .......... .......... .......... 79% 89.5M 0s
    ##  63650K .......... .......... .......... .......... .......... 79% 94.4M 0s
    ##  63700K .......... .......... .......... .......... .......... 79% 76.1M 0s
    ##  63750K .......... .......... .......... .......... .......... 79% 56.5M 0s
    ##  63800K .......... .......... .......... .......... .......... 79% 82.1M 0s
    ##  63850K .......... .......... .......... .......... .......... 79% 90.4M 0s
    ##  63900K .......... .......... .......... .......... .......... 80% 99.4M 0s
    ##  63950K .......... .......... .......... .......... .......... 80% 46.1M 0s
    ##  64000K .......... .......... .......... .......... .......... 80% 69.7M 0s
    ##  64050K .......... .......... .......... .......... .......... 80% 77.1M 0s
    ##  64100K .......... .......... .......... .......... .......... 80% 86.1M 0s
    ##  64150K .......... .......... .......... .......... .......... 80% 92.8M 0s
    ##  64200K .......... .......... .......... .......... .......... 80% 58.9M 0s
    ##  64250K .......... .......... .......... .......... .......... 80% 74.7M 0s
    ##  64300K .......... .......... .......... .......... .......... 80% 87.0M 0s
    ##  64350K .......... .......... .......... .......... .......... 80%  111M 0s
    ##  64400K .......... .......... .......... .......... .......... 80% 16.0M 0s
    ##  64450K .......... .......... .......... .......... .......... 80% 66.0M 0s
    ##  64500K .......... .......... .......... .......... .......... 80% 91.1M 0s
    ##  64550K .......... .......... .......... .......... .......... 80%  104M 0s
    ##  64600K .......... .......... .......... .......... .......... 80% 78.8M 0s
    ##  64650K .......... .......... .......... .......... .......... 80% 87.1M 0s
    ##  64700K .......... .......... .......... .......... .......... 81% 87.1M 0s
    ##  64750K .......... .......... .......... .......... .......... 81% 37.9M 0s
    ##  64800K .......... .......... .......... .......... .......... 81% 83.4M 0s
    ##  64850K .......... .......... .......... .......... .......... 81% 35.3M 0s
    ##  64900K .......... .......... .......... .......... .......... 81% 64.3M 0s
    ##  64950K .......... .......... .......... .......... .......... 81% 85.2M 0s
    ##  65000K .......... .......... .......... .......... .......... 81% 76.2M 0s
    ##  65050K .......... .......... .......... .......... .......... 81%  101M 0s
    ##  65100K .......... .......... .......... .......... .......... 81% 77.2M 0s
    ##  65150K .......... .......... .......... .......... .......... 81% 52.3M 0s
    ##  65200K .......... .......... .......... .......... .......... 81% 82.1M 0s
    ##  65250K .......... .......... .......... .......... .......... 81% 62.3M 0s
    ##  65300K .......... .......... .......... .......... .......... 81% 56.9M 0s
    ##  65350K .......... .......... .......... .......... .......... 81% 86.7M 0s
    ##  65400K .......... .......... .......... .......... .......... 81% 14.4M 0s
    ##  65450K .......... .......... .......... .......... .......... 81% 23.2M 0s
    ##  65500K .......... .......... .......... .......... .......... 82% 46.6M 0s
    ##  65550K .......... .......... .......... .......... .......... 82% 68.5M 0s
    ##  65600K .......... .......... .......... .......... .......... 82% 70.3M 0s
    ##  65650K .......... .......... .......... .......... .......... 82% 84.7M 0s
    ##  65700K .......... .......... .......... .......... .......... 82% 88.0M 0s
    ##  65750K .......... .......... .......... .......... .......... 82% 85.0M 0s
    ##  65800K .......... .......... .......... .......... .......... 82% 63.0M 0s
    ##  65850K .......... .......... .......... .......... .......... 82% 61.9M 0s
    ##  65900K .......... .......... .......... .......... .......... 82% 67.2M 0s
    ##  65950K .......... .......... .......... .......... .......... 82% 72.7M 0s
    ##  66000K .......... .......... .......... .......... .......... 82% 57.0M 0s
    ##  66050K .......... .......... .......... .......... .......... 82% 76.2M 0s
    ##  66100K .......... .......... .......... .......... .......... 82% 62.2M 0s
    ##  66150K .......... .......... .......... .......... .......... 82% 64.4M 0s
    ##  66200K .......... .......... .......... .......... .......... 82% 74.3M 0s
    ##  66250K .......... .......... .......... .......... .......... 82% 84.9M 0s
    ##  66300K .......... .......... .......... .......... .......... 83% 64.5M 0s
    ##  66350K .......... .......... .......... .......... .......... 83% 65.4M 0s
    ##  66400K .......... .......... .......... .......... .......... 83% 59.2M 0s
    ##  66450K .......... .......... .......... .......... .......... 83% 86.3M 0s
    ##  66500K .......... .......... .......... .......... .......... 83% 80.6M 0s
    ##  66550K .......... .......... .......... .......... .......... 83%  104M 0s
    ##  66600K .......... .......... .......... .......... .......... 83% 79.1M 0s
    ##  66650K .......... .......... .......... .......... .......... 83%  100M 0s
    ##  66700K .......... .......... .......... .......... .......... 83% 68.2M 0s
    ##  66750K .......... .......... .......... .......... .......... 83%  104M 0s
    ##  66800K .......... .......... .......... .......... .......... 83% 80.0M 0s
    ##  66850K .......... .......... .......... .......... .......... 83% 52.8M 0s
    ##  66900K .......... .......... .......... .......... .......... 83% 76.2M 0s
    ##  66950K .......... .......... .......... .......... .......... 83% 57.2M 0s
    ##  67000K .......... .......... .......... .......... .......... 83% 65.8M 0s
    ##  67050K .......... .......... .......... .......... .......... 83% 84.6M 0s
    ##  67100K .......... .......... .......... .......... .......... 84% 86.0M 0s
    ##  67150K .......... .......... .......... .......... .......... 84% 93.9M 0s
    ##  67200K .......... .......... .......... .......... .......... 84% 78.8M 0s
    ##  67250K .......... .......... .......... .......... .......... 84% 87.1M 0s
    ##  67300K .......... .......... .......... .......... .......... 84% 46.5M 0s
    ##  67350K .......... .......... .......... .......... .......... 84% 53.3M 0s
    ##  67400K .......... .......... .......... .......... .......... 84% 18.4M 0s
    ##  67450K .......... .......... .......... .......... .......... 84% 89.7M 0s
    ##  67500K .......... .......... .......... .......... .......... 84% 81.2M 0s
    ##  67550K .......... .......... .......... .......... .......... 84% 91.0M 0s
    ##  67600K .......... .......... .......... .......... .......... 84% 70.2M 0s
    ##  67650K .......... .......... .......... .......... .......... 84% 70.2M 0s
    ##  67700K .......... .......... .......... .......... .......... 84% 62.2M 0s
    ##  67750K .......... .......... .......... .......... .......... 84% 88.0M 0s
    ##  67800K .......... .......... .......... .......... .......... 84% 24.2M 0s
    ##  67850K .......... .......... .......... .......... .......... 84% 59.9M 0s
    ##  67900K .......... .......... .......... .......... .......... 85% 90.1M 0s
    ##  67950K .......... .......... .......... .......... .......... 85% 65.6M 0s
    ##  68000K .......... .......... .......... .......... .......... 85% 66.2M 0s
    ##  68050K .......... .......... .......... .......... .......... 85% 74.2M 0s
    ##  68100K .......... .......... .......... .......... .......... 85% 62.6M 0s
    ##  68150K .......... .......... .......... .......... .......... 85% 96.4M 0s
    ##  68200K .......... .......... .......... .......... .......... 85% 64.7M 0s
    ##  68250K .......... .......... .......... .......... .......... 85% 55.9M 0s
    ##  68300K .......... .......... .......... .......... .......... 85% 78.0M 0s
    ##  68350K .......... .......... .......... .......... .......... 85% 71.2M 0s
    ##  68400K .......... .......... .......... .......... .......... 85% 63.6M 0s
    ##  68450K .......... .......... .......... .......... .......... 85% 79.2M 0s
    ##  68500K .......... .......... .......... .......... .......... 85% 80.6M 0s
    ##  68550K .......... .......... .......... .......... .......... 85% 85.2M 0s
    ##  68600K .......... .......... .......... .......... .......... 85% 76.3M 0s
    ##  68650K .......... .......... .......... .......... .......... 85%  103M 0s
    ##  68700K .......... .......... .......... .......... .......... 86% 78.1M 0s
    ##  68750K .......... .......... .......... .......... .......... 86% 58.3M 0s
    ##  68800K .......... .......... .......... .......... .......... 86% 60.5M 0s
    ##  68850K .......... .......... .......... .......... .......... 86% 73.4M 0s
    ##  68900K .......... .......... .......... .......... .......... 86% 84.2M 0s
    ##  68950K .......... .......... .......... .......... .......... 86%  126M 0s
    ##  69000K .......... .......... .......... .......... .......... 86% 76.2M 0s
    ##  69050K .......... .......... .......... .......... .......... 86% 60.7M 0s
    ##  69100K .......... .......... .......... .......... .......... 86% 64.8M 0s
    ##  69150K .......... .......... .......... .......... .......... 86%  110M 0s
    ##  69200K .......... .......... .......... .......... .......... 86% 93.5M 0s
    ##  69250K .......... .......... .......... .......... .......... 86%  103M 0s
    ##  69300K .......... .......... .......... .......... .......... 86% 37.4M 0s
    ##  69350K .......... .......... .......... .......... .......... 86% 91.4M 0s
    ##  69400K .......... .......... .......... .......... .......... 86% 27.4M 0s
    ##  69450K .......... .......... .......... .......... .......... 86%  103M 0s
    ##  69500K .......... .......... .......... .......... .......... 87% 48.4M 0s
    ##  69550K .......... .......... .......... .......... .......... 87% 89.9M 0s
    ##  69600K .......... .......... .......... .......... .......... 87% 88.4M 0s
    ##  69650K .......... .......... .......... .......... .......... 87% 80.8M 0s
    ##  69700K .......... .......... .......... .......... .......... 87%  107M 0s
    ##  69750K .......... .......... .......... .......... .......... 87% 15.6M 0s
    ##  69800K .......... .......... .......... .......... .......... 87% 55.8M 0s
    ##  69850K .......... .......... .......... .......... .......... 87%  120M 0s
    ##  69900K .......... .......... .......... .......... .......... 87% 93.4M 0s
    ##  69950K .......... .......... .......... .......... .......... 87% 84.4M 0s
    ##  70000K .......... .......... .......... .......... .......... 87% 98.1M 0s
    ##  70050K .......... .......... .......... .......... .......... 87%  109M 0s
    ##  70100K .......... .......... .......... .......... .......... 87% 48.8M 0s
    ##  70150K .......... .......... .......... .......... .......... 87% 37.8M 0s
    ##  70200K .......... .......... .......... .......... .......... 87% 63.5M 0s
    ##  70250K .......... .......... .......... .......... .......... 87% 58.0M 0s
    ##  70300K .......... .......... .......... .......... .......... 88% 79.8M 0s
    ##  70350K .......... .......... .......... .......... .......... 88%  123M 0s
    ##  70400K .......... .......... .......... .......... .......... 88% 84.8M 0s
    ##  70450K .......... .......... .......... .......... .......... 88%  120M 0s
    ##  70500K .......... .......... .......... .......... .......... 88% 36.1M 0s
    ##  70550K .......... .......... .......... .......... .......... 88% 99.9M 0s
    ##  70600K .......... .......... .......... .......... .......... 88% 61.5M 0s
    ##  70650K .......... .......... .......... .......... .......... 88% 62.3M 0s
    ##  70700K .......... .......... .......... .......... .......... 88% 25.2M 0s
    ##  70750K .......... .......... .......... .......... .......... 88%  114M 0s
    ##  70800K .......... .......... .......... .......... .......... 88% 98.3M 0s
    ##  70850K .......... .......... .......... .......... .......... 88% 62.6M 0s
    ##  70900K .......... .......... .......... .......... .......... 88% 85.6M 0s
    ##  70950K .......... .......... .......... .......... .......... 88%  105M 0s
    ##  71000K .......... .......... .......... .......... .......... 88% 94.4M 0s
    ##  71050K .......... .......... .......... .......... .......... 88%  114M 0s
    ##  71100K .......... .......... .......... .......... .......... 89% 80.7M 0s
    ##  71150K .......... .......... .......... .......... .......... 89% 55.9M 0s
    ##  71200K .......... .......... .......... .......... .......... 89% 75.3M 0s
    ##  71250K .......... .......... .......... .......... .......... 89% 79.6M 0s
    ##  71300K .......... .......... .......... .......... .......... 89% 29.3M 0s
    ##  71350K .......... .......... .......... .......... .......... 89% 83.1M 0s
    ##  71400K .......... .......... .......... .......... .......... 89% 38.4M 0s
    ##  71450K .......... .......... .......... .......... .......... 89% 83.9M 0s
    ##  71500K .......... .......... .......... .......... .......... 89% 30.2M 0s
    ##  71550K .......... .......... .......... .......... .......... 89% 88.5M 0s
    ##  71600K .......... .......... .......... .......... .......... 89%  106M 0s
    ##  71650K .......... .......... .......... .......... .......... 89% 97.4M 0s
    ##  71700K .......... .......... .......... .......... .......... 89% 89.2M 0s
    ##  71750K .......... .......... .......... .......... .......... 89% 60.3M 0s
    ##  71800K .......... .......... .......... .......... .......... 89% 91.1M 0s
    ##  71850K .......... .......... .......... .......... .......... 89%  109M 0s
    ##  71900K .......... .......... .......... .......... .......... 90% 33.3M 0s
    ##  71950K .......... .......... .......... .......... .......... 90% 36.2M 0s
    ##  72000K .......... .......... .......... .......... .......... 90% 74.0M 0s
    ##  72050K .......... .......... .......... .......... .......... 90% 73.9M 0s
    ##  72100K .......... .......... .......... .......... .......... 90% 77.9M 0s
    ##  72150K .......... .......... .......... .......... .......... 90% 53.8M 0s
    ##  72200K .......... .......... .......... .......... .......... 90%  116M 0s
    ##  72250K .......... .......... .......... .......... .......... 90%  119M 0s
    ##  72300K .......... .......... .......... .......... .......... 90% 47.2M 0s
    ##  72350K .......... .......... .......... .......... .......... 90% 87.1M 0s
    ##  72400K .......... .......... .......... .......... .......... 90%  106M 0s
    ##  72450K .......... .......... .......... .......... .......... 90%  113M 0s
    ##  72500K .......... .......... .......... .......... .......... 90% 34.2M 0s
    ##  72550K .......... .......... .......... .......... .......... 90% 68.9M 0s
    ##  72600K .......... .......... .......... .......... .......... 90% 83.8M 0s
    ##  72650K .......... .......... .......... .......... .......... 90% 65.6M 0s
    ##  72700K .......... .......... .......... .......... .......... 91% 65.9M 0s
    ##  72750K .......... .......... .......... .......... .......... 91% 77.0M 0s
    ##  72800K .......... .......... .......... .......... .......... 91% 83.9M 0s
    ##  72850K .......... .......... .......... .......... .......... 91% 98.3M 0s
    ##  72900K .......... .......... .......... .......... .......... 91%  110M 0s
    ##  72950K .......... .......... .......... .......... .......... 91% 56.3M 0s
    ##  73000K .......... .......... .......... .......... .......... 91% 67.9M 0s
    ##  73050K .......... .......... .......... .......... .......... 91% 32.2M 0s
    ##  73100K .......... .......... .......... .......... .......... 91% 72.2M 0s
    ##  73150K .......... .......... .......... .......... .......... 91% 64.6M 0s
    ##  73200K .......... .......... .......... .......... .......... 91% 80.0M 0s
    ##  73250K .......... .......... .......... .......... .......... 91% 41.0M 0s
    ##  73300K .......... .......... .......... .......... .......... 91% 77.9M 0s
    ##  73350K .......... .......... .......... .......... .......... 91%  109M 0s
    ##  73400K .......... .......... .......... .......... .......... 91% 70.6M 0s
    ##  73450K .......... .......... .......... .......... .......... 91% 84.2M 0s
    ##  73500K .......... .......... .......... .......... .......... 92% 62.5M 0s
    ##  73550K .......... .......... .......... .......... .......... 92% 98.7M 0s
    ##  73600K .......... .......... .......... .......... .......... 92% 30.5M 0s
    ##  73650K .......... .......... .......... .......... .......... 92% 58.0M 0s
    ##  73700K .......... .......... .......... .......... .......... 92% 32.8M 0s
    ##  73750K .......... .......... .......... .......... .......... 92% 90.2M 0s
    ##  73800K .......... .......... .......... .......... .......... 92%  103M 0s
    ##  73850K .......... .......... .......... .......... .......... 92%  100M 0s
    ##  73900K .......... .......... .......... .......... .......... 92%  104M 0s
    ##  73950K .......... .......... .......... .......... .......... 92% 24.1M 0s
    ##  74000K .......... .......... .......... .......... .......... 92% 95.0M 0s
    ##  74050K .......... .......... .......... .......... .......... 92% 94.7M 0s
    ##  74100K .......... .......... .......... .......... .......... 92% 35.4M 0s
    ##  74150K .......... .......... .......... .......... .......... 92%  122M 0s
    ##  74200K .......... .......... .......... .......... .......... 92% 64.8M 0s
    ##  74250K .......... .......... .......... .......... .......... 92% 20.3M 0s
    ##  74300K .......... .......... .......... .......... .......... 93% 61.8M 0s
    ##  74350K .......... .......... .......... .......... .......... 93% 96.4M 0s
    ##  74400K .......... .......... .......... .......... .......... 93% 99.5M 0s
    ##  74450K .......... .......... .......... .......... .......... 93%  116M 0s
    ##  74500K .......... .......... .......... .......... .......... 93% 95.9M 0s
    ##  74550K .......... .......... .......... .......... .......... 93% 88.7M 0s
    ##  74600K .......... .......... .......... .......... .......... 93% 99.4M 0s
    ##  74650K .......... .......... .......... .......... .......... 93% 24.4M 0s
    ##  74700K .......... .......... .......... .......... .......... 93% 88.1M 0s
    ##  74750K .......... .......... .......... .......... .......... 93%  112M 0s
    ##  74800K .......... .......... .......... .......... .......... 93% 87.3M 0s
    ##  74850K .......... .......... .......... .......... .......... 93% 48.5M 0s
    ##  74900K .......... .......... .......... .......... .......... 93% 67.4M 0s
    ##  74950K .......... .......... .......... .......... .......... 93% 88.8M 0s
    ##  75000K .......... .......... .......... .......... .......... 93%  106M 0s
    ##  75050K .......... .......... .......... .......... .......... 93% 49.0M 0s
    ##  75100K .......... .......... .......... .......... .......... 94% 54.2M 0s
    ##  75150K .......... .......... .......... .......... .......... 94%  116M 0s
    ##  75200K .......... .......... .......... .......... .......... 94% 57.8M 0s
    ##  75250K .......... .......... .......... .......... .......... 94% 79.5M 0s
    ##  75300K .......... .......... .......... .......... .......... 94%  124M 0s
    ##  75350K .......... .......... .......... .......... .......... 94% 96.9M 0s
    ##  75400K .......... .......... .......... .......... .......... 94% 92.6M 0s
    ##  75450K .......... .......... .......... .......... .......... 94% 66.7M 0s
    ##  75500K .......... .......... .......... .......... .......... 94% 84.6M 0s
    ##  75550K .......... .......... .......... .......... .......... 94% 93.9M 0s
    ##  75600K .......... .......... .......... .......... .......... 94% 21.2M 0s
    ##  75650K .......... .......... .......... .......... .......... 94% 76.4M 0s
    ##  75700K .......... .......... .......... .......... .......... 94% 75.6M 0s
    ##  75750K .......... .......... .......... .......... .......... 94%  119M 0s
    ##  75800K .......... .......... .......... .......... .......... 94% 59.9M 0s
    ##  75850K .......... .......... .......... .......... .......... 94% 99.3M 0s
    ##  75900K .......... .......... .......... .......... .......... 95%  110M 0s
    ##  75950K .......... .......... .......... .......... .......... 95% 34.5M 0s
    ##  76000K .......... .......... .......... .......... .......... 95% 82.8M 0s
    ##  76050K .......... .......... .......... .......... .......... 95% 83.4M 0s
    ##  76100K .......... .......... .......... .......... .......... 95% 32.3M 0s
    ##  76150K .......... .......... .......... .......... .......... 95% 54.9M 0s
    ##  76200K .......... .......... .......... .......... .......... 95% 63.3M 0s
    ##  76250K .......... .......... .......... .......... .......... 95% 77.6M 0s
    ##  76300K .......... .......... .......... .......... .......... 95% 42.5M 0s
    ##  76350K .......... .......... .......... .......... .......... 95% 28.5M 0s
    ##  76400K .......... .......... .......... .......... .......... 95%  103M 0s
    ##  76450K .......... .......... .......... .......... .......... 95%  108M 0s
    ##  76500K .......... .......... .......... .......... .......... 95% 94.0M 0s
    ##  76550K .......... .......... .......... .......... .......... 95% 35.9M 0s
    ##  76600K .......... .......... .......... .......... .......... 95% 97.6M 0s
    ##  76650K .......... .......... .......... .......... .......... 95%  125M 0s
    ##  76700K .......... .......... .......... .......... .......... 96% 86.2M 0s
    ##  76750K .......... .......... .......... .......... .......... 96% 86.3M 0s
    ##  76800K .......... .......... .......... .......... .......... 96% 80.8M 0s
    ##  76850K .......... .......... .......... .......... .......... 96%  104M 0s
    ##  76900K .......... .......... .......... .......... .......... 96% 97.0M 0s
    ##  76950K .......... .......... .......... .......... .......... 96% 22.9M 0s
    ##  77000K .......... .......... .......... .......... .......... 96% 83.0M 0s
    ##  77050K .......... .......... .......... .......... .......... 96% 95.7M 0s
    ##  77100K .......... .......... .......... .......... .......... 96% 49.6M 0s
    ##  77150K .......... .......... .......... .......... .......... 96% 96.7M 0s
    ##  77200K .......... .......... .......... .......... .......... 96% 66.4M 0s
    ##  77250K .......... .......... .......... .......... .......... 96% 83.0M 0s
    ##  77300K .......... .......... .......... .......... .......... 96% 90.2M 0s
    ##  77350K .......... .......... .......... .......... .......... 96% 82.3M 0s
    ##  77400K .......... .......... .......... .......... .......... 96%  106M 0s
    ##  77450K .......... .......... .......... .......... .......... 96%  109M 0s
    ##  77500K .......... .......... .......... .......... .......... 97% 83.3M 0s
    ##  77550K .......... .......... .......... .......... .......... 97% 86.4M 0s
    ##  77600K .......... .......... .......... .......... .......... 97% 80.2M 0s
    ##  77650K .......... .......... .......... .......... .......... 97% 27.7M 0s
    ##  77700K .......... .......... .......... .......... .......... 97% 75.2M 0s
    ##  77750K .......... .......... .......... .......... .......... 97% 79.1M 0s
    ##  77800K .......... .......... .......... .......... .......... 97% 95.2M 0s
    ##  77850K .......... .......... .......... .......... .......... 97% 44.2M 0s
    ##  77900K .......... .......... .......... .......... .......... 97%  102M 0s
    ##  77950K .......... .......... .......... .......... .......... 97% 76.7M 0s
    ##  78000K .......... .......... .......... .......... .......... 97% 48.9M 0s
    ##  78050K .......... .......... .......... .......... .......... 97% 48.3M 0s
    ##  78100K .......... .......... .......... .......... .......... 97% 75.5M 0s
    ##  78150K .......... .......... .......... .......... .......... 97% 83.8M 0s
    ##  78200K .......... .......... .......... .......... .......... 97% 32.6M 0s
    ##  78250K .......... .......... .......... .......... .......... 97% 69.3M 0s
    ##  78300K .......... .......... .......... .......... .......... 98% 59.6M 0s
    ##  78350K .......... .......... .......... .......... .......... 98%  122M 0s
    ##  78400K .......... .......... .......... .......... .......... 98% 55.2M 0s
    ##  78450K .......... .......... .......... .......... .......... 98% 98.6M 0s
    ##  78500K .......... .......... .......... .......... .......... 98% 75.3M 0s
    ##  78550K .......... .......... .......... .......... .......... 98%  102M 0s
    ##  78600K .......... .......... .......... .......... .......... 98% 44.1M 0s
    ##  78650K .......... .......... .......... .......... .......... 98% 47.7M 0s
    ##  78700K .......... .......... .......... .......... .......... 98% 59.0M 0s
    ##  78750K .......... .......... .......... .......... .......... 98%  104M 0s
    ##  78800K .......... .......... .......... .......... .......... 98% 34.2M 0s
    ##  78850K .......... .......... .......... .......... .......... 98% 94.1M 0s
    ##  78900K .......... .......... .......... .......... .......... 98% 59.5M 0s
    ##  78950K .......... .......... .......... .......... .......... 98% 97.4M 0s
    ##  79000K .......... .......... .......... .......... .......... 98% 88.8M 0s
    ##  79050K .......... .......... .......... .......... .......... 98% 79.3M 0s
    ##  79100K .......... .......... .......... .......... .......... 99% 77.6M 0s
    ##  79150K .......... .......... .......... .......... .......... 99% 44.9M 0s
    ##  79200K .......... .......... .......... .......... .......... 99% 79.4M 0s
    ##  79250K .......... .......... .......... .......... .......... 99%  118M 0s
    ##  79300K .......... .......... .......... .......... .......... 99% 50.0M 0s
    ##  79350K .......... .......... .......... .......... .......... 99%  102M 0s
    ##  79400K .......... .......... .......... .......... .......... 99% 39.9M 0s
    ##  79450K .......... .......... .......... .......... .......... 99% 82.0M 0s
    ##  79500K .......... .......... .......... .......... .......... 99% 97.4M 0s
    ##  79550K .......... .......... .......... .......... .......... 99% 83.0M 0s
    ##  79600K .......... .......... .......... .......... .......... 99% 80.4M 0s
    ##  79650K .......... .......... .......... .......... .......... 99% 57.9M 0s
    ##  79700K .......... .......... .......... .......... .......... 99% 65.5M 0s
    ##  79750K .......... .......... .......... .......... .......... 99% 99.0M 0s
    ##  79800K .......... .......... .......... .......... .......... 99% 73.2M 0s
    ##  79850K .......... .......... .......... .......... .......... 99%  120M 0s
    ##  79900K .......... .......... ..                              100%  120M=1.5s
    ## 
    ## 2020-11-27 16:31:48 (52.8 MB/s) - ‘silva_species_assignment_v138.fa.gz.5’ saved [81840166/81840166]

Ici nous créons une variable qui va recevoir les espèces obtenues grâce
à Silva

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
```

``` r
taxa <- addSpecies(taxa, "~/silva_species_assignment_v138.fa.gz")
```

On remarque donc après avoir affiché la table qu’on a créée on obtient
une majorité de Bacteroidetes ce qui est normal dans des échantillons
fécaux. D’autres espèces n’ont pas pu être assignées car on a peu de
données sur les bactéries des intestins des souris.

``` r
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum         Class         Order           Family          
    ## [1,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [2,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [3,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [4,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [5,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
    ## [6,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ##      Genus         Species
    ## [1,] NA            NA     
    ## [2,] NA            NA     
    ## [3,] NA            NA     
    ## [4,] NA            NA     
    ## [5,] "Bacteroides" NA     
    ## [6,] NA            NA

# Evaluate accuracy

Ici on cherche à comparer les variants donnés par la machine par rapport
à la composition de communauté attendue. Le premier code nous montre
qu’on s’attendait à avoir 20 souches différentes.

``` r
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
```

    ## DADA2 inferred 20 sample sequences present in the Mock community.

Dada2 a indentifié 20 ASVs. Cela correspond précisement à ce à quoi on
s’attendait, donc le taux d’erreur résiduel par dada2 est de 0%

``` r
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```

    ## Of those, 20 were exact matches to the expected reference sequences.

# Sauvegarde de l’environnement pour pouvoir jouer le code de phyloseq

``` r
save.image(file="03_DADA2_tutorial_FinalEnv")
```
