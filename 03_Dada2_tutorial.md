Dada2 tutorial
================

  - [Préparation de l’environnement](#préparation-de-lenvironnement)
      - [Chargement des librairies Phyloseq et Dada2 afin de pouvoir
        utiliser ses
        packages](#chargement-des-librairies-phyloseq-et-dada2-afin-de-pouvoir-utiliser-ses-packages)
      - [Création d’une variable ‘path’ dans laquelle on met l’objet
        “Miseq\_SOP”. Ensuite on les
        liste](#création-dune-variable-path-dans-laquelle-on-met-lobjet-miseq_sop.-ensuite-on-les-liste)
      - [Création de deux listes. Une avec les Forwards et une avec les
        Reverses qui correspondent respectivement au Read1 et Read2. La
        commande avec ‘sample.names’ permet de mettre toutes les
        séquences sous un même format pour harmoniser leurs
        noms.](#création-de-deux-listes.-une-avec-les-forwards-et-une-avec-les-reverses-qui-correspondent-respectivement-au-read1-et-read2.-la-commande-avec-sample.names-permet-de-mettre-toutes-les-séquences-sous-un-même-format-pour-harmoniser-leurs-noms.)
  - [Nous allons indiquer à la machine quels paramètres nous allons
    utiliser pour filtrer les séquences avant de les ranger. On indique
    ainsi dans la première fonction les deux fichiers d’où les séquences
    viendront (fnFs et fnRs), puis le fichier où elles seront rangées
    (filtFs et filtRs crées juste au dessus). Enuite nous allons
    indiquer où couper pour les deux sortes de
    séquences.](#nous-allons-indiquer-à-la-machine-quels-paramètres-nous-allons-utiliser-pour-filtrer-les-séquences-avant-de-les-ranger.-on-indique-ainsi-dans-la-première-fonction-les-deux-fichiers-doù-les-séquences-viendront-fnfs-et-fnrs-puis-le-fichier-où-elles-seront-rangées-filtfs-et-filtrs-crées-juste-au-dessus.-enuite-nous-allons-indiquer-où-couper-pour-les-deux-sortes-de-séquences.)
      - [Cette ligne nous permet de visualiser les erreurs qu’on a fait
        apprendre à la
        machine](#cette-ligne-nous-permet-de-visualiser-les-erreurs-quon-a-fait-apprendre-à-la-machine)
      - [Cette commande nous permet de visualiser le résultat global
        qu’on retrouve classé dans la liste dadaFs. Ils nous indiquent
        que sur les séquences on retrouve 128 séquences qui
        correspondent aux vrais variants, par rapport aux 1979
        séquences. Ils nous indiquent aussi les diagnostiques de
        qualité.](#cette-commande-nous-permet-de-visualiser-le-résultat-global-quon-retrouve-classé-dans-la-liste-dadafs.-ils-nous-indiquent-que-sur-les-séquences-on-retrouve-128-séquences-qui-correspondent-aux-vrais-variants-par-rapport-aux-1979-séquences.-ils-nous-indiquent-aussi-les-diagnostiques-de-qualité.)
      - [Ici on peut voir qu’on à 3% de séquences chimériques dans notre
        jeu de
        donnée.](#ici-on-peut-voir-quon-à-3-de-séquences-chimériques-dans-notre-jeu-de-donnée.)
      - [Ici nous créons une variable qui va recevoir les espèces
        obtenues grâce à
        Silva](#ici-nous-créons-une-variable-qui-va-recevoir-les-espèces-obtenues-grâce-à-silva)
      - [Ici on remarque donc après avoir affiché la table qu’on a créée
        on obtient une majorité de Bacteroidetes ce qui est normal dans
        des échantillons fécaux. D’autres espèces n’ont pas pu être
        assignées car on a peu de données sur les bactéries des
        intestins des
        souris.](#ici-on-remarque-donc-après-avoir-affiché-la-table-quon-a-créée-on-obtient-une-majorité-de-bacteroidetes-ce-qui-est-normal-dans-des-échantillons-fécaux.-dautres-espèces-nont-pas-pu-être-assignées-car-on-a-peu-de-données-sur-les-bactéries-des-intestins-des-souris.)
      - [Dada2 a indentifié 20 ASVs. Cela correspond précisement à ce à
        quoi on s’attendait, donc le taux d’erreur résiduel par dada2
        est de
        0%](#dada2-a-indentifié-20-asvs.-cela-correspond-précisement-à-ce-à-quoi-on-sattendait-donc-le-taux-derreur-résiduel-par-dada2-est-de-0)

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

## Création de deux listes. Une avec les Forwards et une avec les Reverses qui correspondent respectivement au Read1 et Read2. La commande avec ‘sample.names’ permet de mettre toutes les séquences sous un même format pour harmoniser leurs noms.

``` r
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

\#Inspect read quality profiles \#\# Permet d’inspecter la qualité des
différents reads et nous permet de voir à partir de quel nucléotide nous
allons avoir une baisse du score de qualité pour pouvoir par la suite
retirer les nucléotides de mauvaise qualité. Sur les forwards nous
allons retirer les 10 derniers nucléotides avec une commande retrouvée
plus bas dans le script.

``` r
plotQualityProfile(fnRs[1:4])
```

![](03_Dada2_tutorial_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
\#\# Ici nous faisons la même manipulation mais avec les Reverses qui
sont d’un peu moins bonne qualité que les Forwards, dû au sequençage
avec Illumina.Ici, nous retirerons tous les nucléotides à partir du 160e

``` r
plotQualityProfile(fnRs[1:4])
```

![](03_Dada2_tutorial_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
\#Filter and trim \#\# Ici nous créeons une nouvelle variable qui
recevra tous les fichiers filtrés, que ce soit pour les Forwards ou les
Reverses. Nous allons en même temps indiquer à la machine que les noms
utilisés pour ranger les séquences dans les dossiers seront les mêmes
que nous avons standardisé plus haut.

``` r
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

# Nous allons indiquer à la machine quels paramètres nous allons utiliser pour filtrer les séquences avant de les ranger. On indique ainsi dans la première fonction les deux fichiers d’où les séquences viendront (fnFs et fnRs), puis le fichier où elles seront rangées (filtFs et filtRs crées juste au dessus). Enuite nous allons indiquer où couper pour les deux sortes de séquences.

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

\#Learn the Error Rates \#\# Nous allons ici utiliser des lignes de
commandes qui vont permettre d’apprendre à la machine les différents
profils d’erreurs générées lors du séquençage. L’opération est faite sur
les deux types de séquence.

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

## Cette ligne nous permet de visualiser les erreurs qu’on a fait apprendre à la machine

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](03_Dada2_tutorial_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

\#Sample Inference \#\# Ici nous créons une autre variable “dadaFs” dans
laquelle nous mettons les fichiers obtenus après avoir filtré et
appliqué le profil d’erreur à nos séquences. Nous allons faire la même
chose avec dadaRS.

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

## Cette commande nous permet de visualiser le résultat global qu’on retrouve classé dans la liste dadaFs. Ils nous indiquent que sur les séquences on retrouve 128 séquences qui correspondent aux vrais variants, par rapport aux 1979 séquences. Ils nous indiquent aussi les diagnostiques de qualité.

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 128 sequence variants were inferred from 1979 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

\#Merge paired reads \#\# Ici nous voulons mettre en une seule séquence
les Forwards et les Reverses.Nous pouvons faire cette opération grâce
aux overlaps de 12 paires de base. Cela se fait grâce à un alignement
entre les forwards et les reverses qui vont permettre de contruire les
contigs.

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

\#Construct sequence table \#\# Ici nous allons construire une table des
variations de séquence dans les amplicons (ASV) qui permet une meilleure
résolution que les tables OTUs 97%

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

\#Remove chimeras \#\# Malgré qu’on ait pu appliquer les modèles
d’erreurs aux séquences, il reste des chimères. Ces chimères sont
facilement reconnaissables par la machine et peuvent etre réparées en y
rajoutant les parties droites et gauche des 2 séquences les plus
abondantes.

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 61 bimeras out of 293 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]  20 232

## Ici on peut voir qu’on à 3% de séquences chimériques dans notre jeu de donnée.

``` r
1-sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.03596257

\#Track reads through the pipeline \#\# Ce code nous permet de
visualiser le nombre de séquences obtenues à la suite de toutes nos
manipulations de filtrage. Ici nous pouvons voir qu’on a pu récupérer la
plupart de nos séquences brutes, ce qui est signe d’une bonne qualité de
séquençage.

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

\#Assign taxonomy \#\# Ici nous avons du récupérer silva afin d’analyser
et d’assigner les taxonomies.

``` bash
wget https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz
```

    ## --2020-11-26 19:48:39--  https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz
    ## Resolving zenodo.org (zenodo.org)... 137.138.76.77
    ## Connecting to zenodo.org (zenodo.org)|137.138.76.77|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 81840166 (78M) [application/octet-stream]
    ## Saving to: ‘silva_species_assignment_v138.fa.gz.1’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0% 7.12M 11s
    ##     50K .......... .......... .......... .......... ..........  0% 7.03M 11s
    ##    100K .......... .......... .......... .......... ..........  0% 13.0M 9s
    ##    150K .......... .......... .......... .......... ..........  0% 11.5M 9s
    ##    200K .......... .......... .......... .......... ..........  0% 11.6M 8s
    ##    250K .......... .......... .......... .......... ..........  0% 68.0M 7s
    ##    300K .......... .......... .......... .......... ..........  0% 86.3M 6s
    ##    350K .......... .......... .......... .......... ..........  0% 14.9M 6s
    ##    400K .......... .......... .......... .......... ..........  0% 86.3M 5s
    ##    450K .......... .......... .......... .......... ..........  0% 35.3M 5s
    ##    500K .......... .......... .......... .......... ..........  0% 91.6M 5s
    ##    550K .......... .......... .......... .......... ..........  0% 66.8M 4s
    ##    600K .......... .......... .......... .......... ..........  0%  102M 4s
    ##    650K .......... .......... .......... .......... ..........  0%  111M 4s
    ##    700K .......... .......... .......... .......... ..........  0% 89.3M 4s
    ##    750K .......... .......... .......... .......... ..........  1%  122M 4s
    ##    800K .......... .......... .......... .......... ..........  1% 28.6M 3s
    ##    850K .......... .......... .......... .......... ..........  1%  103M 3s
    ##    900K .......... .......... .......... .......... ..........  1% 43.6M 3s
    ##    950K .......... .......... .......... .......... ..........  1% 31.5M 3s
    ##   1000K .......... .......... .......... .......... ..........  1% 56.0M 3s
    ##   1050K .......... .......... .......... .......... ..........  1% 66.0M 3s
    ##   1100K .......... .......... .......... .......... ..........  1% 63.4M 3s
    ##   1150K .......... .......... .......... .......... ..........  1% 86.4M 3s
    ##   1200K .......... .......... .......... .......... ..........  1% 90.8M 3s
    ##   1250K .......... .......... .......... .......... ..........  1%  132M 3s
    ##   1300K .......... .......... .......... .......... ..........  1% 72.7M 3s
    ##   1350K .......... .......... .......... .......... ..........  1% 97.7M 3s
    ##   1400K .......... .......... .......... .......... ..........  1% 88.2M 2s
    ##   1450K .......... .......... .......... .......... ..........  1%  106M 2s
    ##   1500K .......... .......... .......... .......... ..........  1% 53.4M 2s
    ##   1550K .......... .......... .......... .......... ..........  2% 69.1M 2s
    ##   1600K .......... .......... .......... .......... ..........  2% 25.5M 2s
    ##   1650K .......... .......... .......... .......... ..........  2% 64.2M 2s
    ##   1700K .......... .......... .......... .......... ..........  2% 58.6M 2s
    ##   1750K .......... .......... .......... .......... ..........  2% 65.2M 2s
    ##   1800K .......... .......... .......... .......... ..........  2% 59.9M 2s
    ##   1850K .......... .......... .......... .......... ..........  2% 74.3M 2s
    ##   1900K .......... .......... .......... .......... ..........  2% 62.7M 2s
    ##   1950K .......... .......... .......... .......... ..........  2% 56.4M 2s
    ##   2000K .......... .......... .......... .......... ..........  2% 62.1M 2s
    ##   2050K .......... .......... .......... .......... ..........  2% 61.5M 2s
    ##   2100K .......... .......... .......... .......... ..........  2% 54.0M 2s
    ##   2150K .......... .......... .......... .......... ..........  2% 64.1M 2s
    ##   2200K .......... .......... .......... .......... ..........  2% 62.4M 2s
    ##   2250K .......... .......... .......... .......... ..........  2% 75.7M 2s
    ##   2300K .......... .......... .......... .......... ..........  2% 73.8M 2s
    ##   2350K .......... .......... .......... .......... ..........  3% 59.6M 2s
    ##   2400K .......... .......... .......... .......... ..........  3% 95.6M 2s
    ##   2450K .......... .......... .......... .......... ..........  3% 86.7M 2s
    ##   2500K .......... .......... .......... .......... ..........  3% 67.1M 2s
    ##   2550K .......... .......... .......... .......... ..........  3%  107M 2s
    ##   2600K .......... .......... .......... .......... ..........  3% 81.4M 2s
    ##   2650K .......... .......... .......... .......... ..........  3% 92.2M 2s
    ##   2700K .......... .......... .......... .......... ..........  3% 72.0M 2s
    ##   2750K .......... .......... .......... .......... ..........  3% 33.7M 2s
    ##   2800K .......... .......... .......... .......... ..........  3% 55.1M 2s
    ##   2850K .......... .......... .......... .......... ..........  3% 53.3M 2s
    ##   2900K .......... .......... .......... .......... ..........  3% 41.5M 2s
    ##   2950K .......... .......... .......... .......... ..........  3% 61.4M 2s
    ##   3000K .......... .......... .......... .......... ..........  3% 81.8M 2s
    ##   3050K .......... .......... .......... .......... ..........  3% 86.3M 2s
    ##   3100K .......... .......... .......... .......... ..........  3% 88.1M 2s
    ##   3150K .......... .......... .......... .......... ..........  4%  107M 2s
    ##   3200K .......... .......... .......... .......... ..........  4% 88.4M 2s
    ##   3250K .......... .......... .......... .......... ..........  4% 91.9M 2s
    ##   3300K .......... .......... .......... .......... ..........  4% 82.0M 2s
    ##   3350K .......... .......... .......... .......... ..........  4% 82.9M 2s
    ##   3400K .......... .......... .......... .......... ..........  4% 62.2M 2s
    ##   3450K .......... .......... .......... .......... ..........  4% 76.9M 2s
    ##   3500K .......... .......... .......... .......... ..........  4% 41.2M 2s
    ##   3550K .......... .......... .......... .......... ..........  4% 67.4M 2s
    ##   3600K .......... .......... .......... .......... ..........  4% 79.0M 2s
    ##   3650K .......... .......... .......... .......... ..........  4% 97.7M 2s
    ##   3700K .......... .......... .......... .......... ..........  4% 36.3M 2s
    ##   3750K .......... .......... .......... .......... ..........  4% 91.1M 2s
    ##   3800K .......... .......... .......... .......... ..........  4% 78.7M 2s
    ##   3850K .......... .......... .......... .......... ..........  4% 83.3M 2s
    ##   3900K .......... .......... .......... .......... ..........  4%  100M 2s
    ##   3950K .......... .......... .......... .......... ..........  5% 80.3M 2s
    ##   4000K .......... .......... .......... .......... ..........  5% 88.1M 2s
    ##   4050K .......... .......... .......... .......... ..........  5% 72.3M 2s
    ##   4100K .......... .......... .......... .......... ..........  5% 65.4M 2s
    ##   4150K .......... .......... .......... .......... ..........  5% 88.8M 2s
    ##   4200K .......... .......... .......... .......... ..........  5% 70.9M 2s
    ##   4250K .......... .......... .......... .......... ..........  5% 69.7M 2s
    ##   4300K .......... .......... .......... .......... ..........  5% 49.1M 2s
    ##   4350K .......... .......... .......... .......... ..........  5% 71.4M 2s
    ##   4400K .......... .......... .......... .......... ..........  5% 97.7M 2s
    ##   4450K .......... .......... .......... .......... ..........  5% 72.9M 2s
    ##   4500K .......... .......... .......... .......... ..........  5% 57.9M 2s
    ##   4550K .......... .......... .......... .......... ..........  5% 71.8M 2s
    ##   4600K .......... .......... .......... .......... ..........  5% 48.8M 2s
    ##   4650K .......... .......... .......... .......... ..........  5% 70.8M 2s
    ##   4700K .......... .......... .......... .......... ..........  5% 35.6M 2s
    ##   4750K .......... .......... .......... .......... ..........  6% 87.9M 2s
    ##   4800K .......... .......... .......... .......... ..........  6% 52.5M 1s
    ##   4850K .......... .......... .......... .......... ..........  6% 88.0M 1s
    ##   4900K .......... .......... .......... .......... ..........  6% 47.3M 1s
    ##   4950K .......... .......... .......... .......... ..........  6% 74.1M 1s
    ##   5000K .......... .......... .......... .......... ..........  6% 71.0M 1s
    ##   5050K .......... .......... .......... .......... ..........  6%  125M 1s
    ##   5100K .......... .......... .......... .......... ..........  6% 75.9M 1s
    ##   5150K .......... .......... .......... .......... ..........  6% 44.4M 1s
    ##   5200K .......... .......... .......... .......... ..........  6% 71.3M 1s
    ##   5250K .......... .......... .......... .......... ..........  6% 53.8M 1s
    ##   5300K .......... .......... .......... .......... ..........  6% 36.6M 1s
    ##   5350K .......... .......... .......... .......... ..........  6% 83.0M 1s
    ##   5400K .......... .......... .......... .......... ..........  6% 78.9M 1s
    ##   5450K .......... .......... .......... .......... ..........  6% 97.1M 1s
    ##   5500K .......... .......... .......... .......... ..........  6% 74.5M 1s
    ##   5550K .......... .......... .......... .......... ..........  7% 77.5M 1s
    ##   5600K .......... .......... .......... .......... ..........  7% 90.7M 1s
    ##   5650K .......... .......... .......... .......... ..........  7% 57.1M 1s
    ##   5700K .......... .......... .......... .......... ..........  7% 41.1M 1s
    ##   5750K .......... .......... .......... .......... ..........  7% 81.8M 1s
    ##   5800K .......... .......... .......... .......... ..........  7% 72.3M 1s
    ##   5850K .......... .......... .......... .......... ..........  7% 97.5M 1s
    ##   5900K .......... .......... .......... .......... ..........  7% 43.1M 1s
    ##   5950K .......... .......... .......... .......... ..........  7%  113M 1s
    ##   6000K .......... .......... .......... .......... ..........  7%  112M 1s
    ##   6050K .......... .......... .......... .......... ..........  7% 78.7M 1s
    ##   6100K .......... .......... .......... .......... ..........  7% 47.9M 1s
    ##   6150K .......... .......... .......... .......... ..........  7%  113M 1s
    ##   6200K .......... .......... .......... .......... ..........  7% 79.4M 1s
    ##   6250K .......... .......... .......... .......... ..........  7% 92.3M 1s
    ##   6300K .......... .......... .......... .......... ..........  7% 44.3M 1s
    ##   6350K .......... .......... .......... .......... ..........  8% 88.6M 1s
    ##   6400K .......... .......... .......... .......... ..........  8% 56.4M 1s
    ##   6450K .......... .......... .......... .......... ..........  8%  112M 1s
    ##   6500K .......... .......... .......... .......... ..........  8% 66.1M 1s
    ##   6550K .......... .......... .......... .......... ..........  8%  108M 1s
    ##   6600K .......... .......... .......... .......... ..........  8%  103M 1s
    ##   6650K .......... .......... .......... .......... ..........  8% 57.3M 1s
    ##   6700K .......... .......... .......... .......... ..........  8% 30.8M 1s
    ##   6750K .......... .......... .......... .......... ..........  8% 96.4M 1s
    ##   6800K .......... .......... .......... .......... ..........  8% 80.3M 1s
    ##   6850K .......... .......... .......... .......... ..........  8% 83.5M 1s
    ##   6900K .......... .......... .......... .......... ..........  8% 87.0M 1s
    ##   6950K .......... .......... .......... .......... ..........  8% 98.5M 1s
    ##   7000K .......... .......... .......... .......... ..........  8%  102M 1s
    ##   7050K .......... .......... .......... .......... ..........  8%  106M 1s
    ##   7100K .......... .......... .......... .......... ..........  8% 36.6M 1s
    ##   7150K .......... .......... .......... .......... ..........  9% 50.1M 1s
    ##   7200K .......... .......... .......... .......... ..........  9% 56.5M 1s
    ##   7250K .......... .......... .......... .......... ..........  9% 76.8M 1s
    ##   7300K .......... .......... .......... .......... ..........  9% 44.6M 1s
    ##   7350K .......... .......... .......... .......... ..........  9% 69.8M 1s
    ##   7400K .......... .......... .......... .......... ..........  9% 45.7M 1s
    ##   7450K .......... .......... .......... .......... ..........  9%  115M 1s
    ##   7500K .......... .......... .......... .......... ..........  9%  119M 1s
    ##   7550K .......... .......... .......... .......... ..........  9% 95.5M 1s
    ##   7600K .......... .......... .......... .......... ..........  9% 95.0M 1s
    ##   7650K .......... .......... .......... .......... ..........  9%  125M 1s
    ##   7700K .......... .......... .......... .......... ..........  9% 75.6M 1s
    ##   7750K .......... .......... .......... .......... ..........  9% 97.1M 1s
    ##   7800K .......... .......... .......... .......... ..........  9% 57.6M 1s
    ##   7850K .......... .......... .......... .......... ..........  9% 35.5M 1s
    ##   7900K .......... .......... .......... .......... ..........  9% 45.0M 1s
    ##   7950K .......... .......... .......... .......... .......... 10% 84.3M 1s
    ##   8000K .......... .......... .......... .......... .......... 10% 88.9M 1s
    ##   8050K .......... .......... .......... .......... .......... 10% 67.6M 1s
    ##   8100K .......... .......... .......... .......... .......... 10% 91.2M 1s
    ##   8150K .......... .......... .......... .......... .......... 10% 94.9M 1s
    ##   8200K .......... .......... .......... .......... .......... 10% 67.2M 1s
    ##   8250K .......... .......... .......... .......... .......... 10% 26.1M 1s
    ##   8300K .......... .......... .......... .......... .......... 10% 29.6M 1s
    ##   8350K .......... .......... .......... .......... .......... 10% 69.8M 1s
    ##   8400K .......... .......... .......... .......... .......... 10% 54.4M 1s
    ##   8450K .......... .......... .......... .......... .......... 10% 42.5M 1s
    ##   8500K .......... .......... .......... .......... .......... 10% 55.7M 1s
    ##   8550K .......... .......... .......... .......... .......... 10% 76.7M 1s
    ##   8600K .......... .......... .......... .......... .......... 10% 83.9M 1s
    ##   8650K .......... .......... .......... .......... .......... 10% 77.2M 1s
    ##   8700K .......... .......... .......... .......... .......... 10% 96.6M 1s
    ##   8750K .......... .......... .......... .......... .......... 11% 77.2M 1s
    ##   8800K .......... .......... .......... .......... .......... 11% 74.3M 1s
    ##   8850K .......... .......... .......... .......... .......... 11% 23.7M 1s
    ##   8900K .......... .......... .......... .......... .......... 11% 86.1M 1s
    ##   8950K .......... .......... .......... .......... .......... 11% 74.4M 1s
    ##   9000K .......... .......... .......... .......... .......... 11% 71.2M 1s
    ##   9050K .......... .......... .......... .......... .......... 11%  101M 1s
    ##   9100K .......... .......... .......... .......... .......... 11% 85.1M 1s
    ##   9150K .......... .......... .......... .......... .......... 11% 73.0M 1s
    ##   9200K .......... .......... .......... .......... .......... 11% 35.2M 1s
    ##   9250K .......... .......... .......... .......... .......... 11% 15.0M 1s
    ##   9300K .......... .......... .......... .......... .......... 11% 71.5M 1s
    ##   9350K .......... .......... .......... .......... .......... 11% 69.2M 1s
    ##   9400K .......... .......... .......... .......... .......... 11% 87.1M 1s
    ##   9450K .......... .......... .......... .......... .......... 11% 53.0M 1s
    ##   9500K .......... .......... .......... .......... .......... 11% 77.7M 1s
    ##   9550K .......... .......... .......... .......... .......... 12% 90.0M 1s
    ##   9600K .......... .......... .......... .......... .......... 12% 80.8M 1s
    ##   9650K .......... .......... .......... .......... .......... 12% 79.3M 1s
    ##   9700K .......... .......... .......... .......... .......... 12% 53.8M 1s
    ##   9750K .......... .......... .......... .......... .......... 12% 75.8M 1s
    ##   9800K .......... .......... .......... .......... .......... 12% 67.6M 1s
    ##   9850K .......... .......... .......... .......... .......... 12% 96.3M 1s
    ##   9900K .......... .......... .......... .......... .......... 12% 52.9M 1s
    ##   9950K .......... .......... .......... .......... .......... 12% 89.1M 1s
    ##  10000K .......... .......... .......... .......... .......... 12% 65.1M 1s
    ##  10050K .......... .......... .......... .......... .......... 12% 80.1M 1s
    ##  10100K .......... .......... .......... .......... .......... 12% 70.3M 1s
    ##  10150K .......... .......... .......... .......... .......... 12% 92.0M 1s
    ##  10200K .......... .......... .......... .......... .......... 12% 64.7M 1s
    ##  10250K .......... .......... .......... .......... .......... 12%  106M 1s
    ##  10300K .......... .......... .......... .......... .......... 12% 66.5M 1s
    ##  10350K .......... .......... .......... .......... .......... 13% 69.2M 1s
    ##  10400K .......... .......... .......... .......... .......... 13% 71.8M 1s
    ##  10450K .......... .......... .......... .......... .......... 13% 77.1M 1s
    ##  10500K .......... .......... .......... .......... .......... 13% 37.0M 1s
    ##  10550K .......... .......... .......... .......... .......... 13% 53.1M 1s
    ##  10600K .......... .......... .......... .......... .......... 13% 74.4M 1s
    ##  10650K .......... .......... .......... .......... .......... 13% 57.6M 1s
    ##  10700K .......... .......... .......... .......... .......... 13% 96.4M 1s
    ##  10750K .......... .......... .......... .......... .......... 13% 81.6M 1s
    ##  10800K .......... .......... .......... .......... .......... 13% 83.8M 1s
    ##  10850K .......... .......... .......... .......... .......... 13% 86.6M 1s
    ##  10900K .......... .......... .......... .......... .......... 13% 72.8M 1s
    ##  10950K .......... .......... .......... .......... .......... 13% 59.3M 1s
    ##  11000K .......... .......... .......... .......... .......... 13% 83.0M 1s
    ##  11050K .......... .......... .......... .......... .......... 13% 88.6M 1s
    ##  11100K .......... .......... .......... .......... .......... 13% 41.7M 1s
    ##  11150K .......... .......... .......... .......... .......... 14% 81.1M 1s
    ##  11200K .......... .......... .......... .......... .......... 14%  102M 1s
    ##  11250K .......... .......... .......... .......... .......... 14% 64.7M 1s
    ##  11300K .......... .......... .......... .......... .......... 14% 91.4M 1s
    ##  11350K .......... .......... .......... .......... .......... 14% 84.5M 1s
    ##  11400K .......... .......... .......... .......... .......... 14% 86.6M 1s
    ##  11450K .......... .......... .......... .......... .......... 14%  105M 1s
    ##  11500K .......... .......... .......... .......... .......... 14% 58.5M 1s
    ##  11550K .......... .......... .......... .......... .......... 14% 50.3M 1s
    ##  11600K .......... .......... .......... .......... .......... 14% 85.2M 1s
    ##  11650K .......... .......... .......... .......... .......... 14% 74.0M 1s
    ##  11700K .......... .......... .......... .......... .......... 14% 72.2M 1s
    ##  11750K .......... .......... .......... .......... .......... 14% 88.8M 1s
    ##  11800K .......... .......... .......... .......... .......... 14% 77.7M 1s
    ##  11850K .......... .......... .......... .......... .......... 14% 49.0M 1s
    ##  11900K .......... .......... .......... .......... .......... 14% 69.0M 1s
    ##  11950K .......... .......... .......... .......... .......... 15% 89.1M 1s
    ##  12000K .......... .......... .......... .......... .......... 15% 58.3M 1s
    ##  12050K .......... .......... .......... .......... .......... 15% 72.0M 1s
    ##  12100K .......... .......... .......... .......... .......... 15% 74.9M 1s
    ##  12150K .......... .......... .......... .......... .......... 15% 67.2M 1s
    ##  12200K .......... .......... .......... .......... .......... 15% 69.9M 1s
    ##  12250K .......... .......... .......... .......... .......... 15% 65.7M 1s
    ##  12300K .......... .......... .......... .......... .......... 15% 72.3M 1s
    ##  12350K .......... .......... .......... .......... .......... 15%  106M 1s
    ##  12400K .......... .......... .......... .......... .......... 15% 76.6M 1s
    ##  12450K .......... .......... .......... .......... .......... 15% 69.5M 1s
    ##  12500K .......... .......... .......... .......... .......... 15% 66.7M 1s
    ##  12550K .......... .......... .......... .......... .......... 15% 62.3M 1s
    ##  12600K .......... .......... .......... .......... .......... 15% 56.9M 1s
    ##  12650K .......... .......... .......... .......... .......... 15%  116M 1s
    ##  12700K .......... .......... .......... .......... .......... 15% 68.6M 1s
    ##  12750K .......... .......... .......... .......... .......... 16% 72.7M 1s
    ##  12800K .......... .......... .......... .......... .......... 16% 84.5M 1s
    ##  12850K .......... .......... .......... .......... .......... 16%  106M 1s
    ##  12900K .......... .......... .......... .......... .......... 16% 65.6M 1s
    ##  12950K .......... .......... .......... .......... .......... 16% 74.6M 1s
    ##  13000K .......... .......... .......... .......... .......... 16% 84.2M 1s
    ##  13050K .......... .......... .......... .......... .......... 16% 89.4M 1s
    ##  13100K .......... .......... .......... .......... .......... 16% 91.0M 1s
    ##  13150K .......... .......... .......... .......... .......... 16% 88.6M 1s
    ##  13200K .......... .......... .......... .......... .......... 16%  101M 1s
    ##  13250K .......... .......... .......... .......... .......... 16% 86.9M 1s
    ##  13300K .......... .......... .......... .......... .......... 16% 73.1M 1s
    ##  13350K .......... .......... .......... .......... .......... 16%  105M 1s
    ##  13400K .......... .......... .......... .......... .......... 16% 68.4M 1s
    ##  13450K .......... .......... .......... .......... .......... 16% 99.5M 1s
    ##  13500K .......... .......... .......... .......... .......... 16% 81.1M 1s
    ##  13550K .......... .......... .......... .......... .......... 17% 88.8M 1s
    ##  13600K .......... .......... .......... .......... .......... 17%  106M 1s
    ##  13650K .......... .......... .......... .......... .......... 17% 72.8M 1s
    ##  13700K .......... .......... .......... .......... .......... 17% 94.6M 1s
    ##  13750K .......... .......... .......... .......... .......... 17%  131M 1s
    ##  13800K .......... .......... .......... .......... .......... 17% 71.8M 1s
    ##  13850K .......... .......... .......... .......... .......... 17% 63.8M 1s
    ##  13900K .......... .......... .......... .......... .......... 17% 68.8M 1s
    ##  13950K .......... .......... .......... .......... .......... 17% 97.7M 1s
    ##  14000K .......... .......... .......... .......... .......... 17% 89.8M 1s
    ##  14050K .......... .......... .......... .......... .......... 17%  115M 1s
    ##  14100K .......... .......... .......... .......... .......... 17%  100M 1s
    ##  14150K .......... .......... .......... .......... .......... 17% 78.2M 1s
    ##  14200K .......... .......... .......... .......... .......... 17% 79.1M 1s
    ##  14250K .......... .......... .......... .......... .......... 17%  102M 1s
    ##  14300K .......... .......... .......... .......... .......... 17% 95.4M 1s
    ##  14350K .......... .......... .......... .......... .......... 18%  150M 1s
    ##  14400K .......... .......... .......... .......... .......... 18% 91.0M 1s
    ##  14450K .......... .......... .......... .......... .......... 18% 77.6M 1s
    ##  14500K .......... .......... .......... .......... .......... 18% 72.6M 1s
    ##  14550K .......... .......... .......... .......... .......... 18%  104M 1s
    ##  14600K .......... .......... .......... .......... .......... 18% 91.7M 1s
    ##  14650K .......... .......... .......... .......... .......... 18%  118M 1s
    ##  14700K .......... .......... .......... .......... .......... 18% 84.1M 1s
    ##  14750K .......... .......... .......... .......... .......... 18%  117M 1s
    ##  14800K .......... .......... .......... .......... .......... 18% 84.9M 1s
    ##  14850K .......... .......... .......... .......... .......... 18% 92.0M 1s
    ##  14900K .......... .......... .......... .......... .......... 18% 76.7M 1s
    ##  14950K .......... .......... .......... .......... .......... 18% 83.7M 1s
    ##  15000K .......... .......... .......... .......... .......... 18% 81.8M 1s
    ##  15050K .......... .......... .......... .......... .......... 18%  131M 1s
    ##  15100K .......... .......... .......... .......... .......... 18% 85.6M 1s
    ##  15150K .......... .......... .......... .......... .......... 19% 90.0M 1s
    ##  15200K .......... .......... .......... .......... .......... 19% 74.3M 1s
    ##  15250K .......... .......... .......... .......... .......... 19%  111M 1s
    ##  15300K .......... .......... .......... .......... .......... 19% 94.6M 1s
    ##  15350K .......... .......... .......... .......... .......... 19%  113M 1s
    ##  15400K .......... .......... .......... .......... .......... 19% 88.5M 1s
    ##  15450K .......... .......... .......... .......... .......... 19% 95.8M 1s
    ##  15500K .......... .......... .......... .......... .......... 19% 74.6M 1s
    ##  15550K .......... .......... .......... .......... .......... 19% 65.4M 1s
    ##  15600K .......... .......... .......... .......... .......... 19%  115M 1s
    ##  15650K .......... .......... .......... .......... .......... 19% 87.8M 1s
    ##  15700K .......... .......... .......... .......... .......... 19% 89.2M 1s
    ##  15750K .......... .......... .......... .......... .......... 19%  105M 1s
    ##  15800K .......... .......... .......... .......... .......... 19%  103M 1s
    ##  15850K .......... .......... .......... .......... .......... 19%  106M 1s
    ##  15900K .......... .......... .......... .......... .......... 19% 96.3M 1s
    ##  15950K .......... .......... .......... .......... .......... 20%  109M 1s
    ##  16000K .......... .......... .......... .......... .......... 20%  106M 1s
    ##  16050K .......... .......... .......... .......... .......... 20%  115M 1s
    ##  16100K .......... .......... .......... .......... .......... 20% 77.8M 1s
    ##  16150K .......... .......... .......... .......... .......... 20% 72.3M 1s
    ##  16200K .......... .......... .......... .......... .......... 20% 79.9M 1s
    ##  16250K .......... .......... .......... .......... .......... 20% 94.7M 1s
    ##  16300K .......... .......... .......... .......... .......... 20% 90.6M 1s
    ##  16350K .......... .......... .......... .......... .......... 20%  133M 1s
    ##  16400K .......... .......... .......... .......... .......... 20% 94.5M 1s
    ##  16450K .......... .......... .......... .......... .......... 20%  104M 1s
    ##  16500K .......... .......... .......... .......... .......... 20% 96.5M 1s
    ##  16550K .......... .......... .......... .......... .......... 20%  111M 1s
    ##  16600K .......... .......... .......... .......... .......... 20% 96.0M 1s
    ##  16650K .......... .......... .......... .......... .......... 20%  120M 1s
    ##  16700K .......... .......... .......... .......... .......... 20% 86.4M 1s
    ##  16750K .......... .......... .......... .......... .......... 21% 99.7M 1s
    ##  16800K .......... .......... .......... .......... .......... 21%  129M 1s
    ##  16850K .......... .......... .......... .......... .......... 21% 97.5M 1s
    ##  16900K .......... .......... .......... .......... .......... 21%  114M 1s
    ##  16950K .......... .......... .......... .......... .......... 21%  107M 1s
    ##  17000K .......... .......... .......... .......... .......... 21%  115M 1s
    ##  17050K .......... .......... .......... .......... .......... 21%  137M 1s
    ##  17100K .......... .......... .......... .......... .......... 21%  106M 1s
    ##  17150K .......... .......... .......... .......... .......... 21% 81.1M 1s
    ##  17200K .......... .......... .......... .......... .......... 21%  106M 1s
    ##  17250K .......... .......... .......... .......... .......... 21% 97.8M 1s
    ##  17300K .......... .......... .......... .......... .......... 21% 98.5M 1s
    ##  17350K .......... .......... .......... .......... .......... 21% 94.3M 1s
    ##  17400K .......... .......... .......... .......... .......... 21%  112M 1s
    ##  17450K .......... .......... .......... .......... .......... 21% 86.4M 1s
    ##  17500K .......... .......... .......... .......... .......... 21%  115M 1s
    ##  17550K .......... .......... .......... .......... .......... 22% 99.7M 1s
    ##  17600K .......... .......... .......... .......... .......... 22% 85.9M 1s
    ##  17650K .......... .......... .......... .......... .......... 22%  104M 1s
    ##  17700K .......... .......... .......... .......... .......... 22%  126M 1s
    ##  17750K .......... .......... .......... .......... .......... 22%  133M 1s
    ##  17800K .......... .......... .......... .......... .......... 22%  119M 1s
    ##  17850K .......... .......... .......... .......... .......... 22% 11.0M 1s
    ##  17900K .......... .......... .......... .......... .......... 22% 74.7M 1s
    ##  17950K .......... .......... .......... .......... .......... 22%  129M 1s
    ##  18000K .......... .......... .......... .......... .......... 22%  139M 1s
    ##  18050K .......... .......... .......... .......... .......... 22% 85.7M 1s
    ##  18100K .......... .......... .......... .......... .......... 22% 73.5M 1s
    ##  18150K .......... .......... .......... .......... .......... 22%  150M 1s
    ##  18200K .......... .......... .......... .......... .......... 22%  102M 1s
    ##  18250K .......... .......... .......... .......... .......... 22% 30.2M 1s
    ##  18300K .......... .......... .......... .......... .......... 22%  109M 1s
    ##  18350K .......... .......... .......... .......... .......... 23% 22.9M 1s
    ##  18400K .......... .......... .......... .......... .......... 23% 67.6M 1s
    ##  18450K .......... .......... .......... .......... .......... 23%  103M 1s
    ##  18500K .......... .......... .......... .......... .......... 23% 92.1M 1s
    ##  18550K .......... .......... .......... .......... .......... 23% 96.3M 1s
    ##  18600K .......... .......... .......... .......... .......... 23%  127M 1s
    ##  18650K .......... .......... .......... .......... .......... 23% 18.1M 1s
    ##  18700K .......... .......... .......... .......... .......... 23%  131M 1s
    ##  18750K .......... .......... .......... .......... .......... 23%  168M 1s
    ##  18800K .......... .......... .......... .......... .......... 23%  103M 1s
    ##  18850K .......... .......... .......... .......... .......... 23%  134M 1s
    ##  18900K .......... .......... .......... .......... .......... 23%  116M 1s
    ##  18950K .......... .......... .......... .......... .......... 23% 91.6M 1s
    ##  19000K .......... .......... .......... .......... .......... 23%  115M 1s
    ##  19050K .......... .......... .......... .......... .......... 23% 42.4M 1s
    ##  19100K .......... .......... .......... .......... .......... 23% 29.0M 1s
    ##  19150K .......... .......... .......... .......... .......... 24% 29.9M 1s
    ##  19200K .......... .......... .......... .......... .......... 24% 74.3M 1s
    ##  19250K .......... .......... .......... .......... .......... 24% 77.3M 1s
    ##  19300K .......... .......... .......... .......... .......... 24%  136M 1s
    ##  19350K .......... .......... .......... .......... .......... 24%  138M 1s
    ##  19400K .......... .......... .......... .......... .......... 24% 31.7M 1s
    ##  19450K .......... .......... .......... .......... .......... 24%  124M 1s
    ##  19500K .......... .......... .......... .......... .......... 24%  110M 1s
    ##  19550K .......... .......... .......... .......... .......... 24%  115M 1s
    ##  19600K .......... .......... .......... .......... .......... 24% 21.7M 1s
    ##  19650K .......... .......... .......... .......... .......... 24% 27.1M 1s
    ##  19700K .......... .......... .......... .......... .......... 24% 83.4M 1s
    ##  19750K .......... .......... .......... .......... .......... 24% 92.3M 1s
    ##  19800K .......... .......... .......... .......... .......... 24% 82.0M 1s
    ##  19850K .......... .......... .......... .......... .......... 24% 96.5M 1s
    ##  19900K .......... .......... .......... .......... .......... 24%  103M 1s
    ##  19950K .......... .......... .......... .......... .......... 25%  101M 1s
    ##  20000K .......... .......... .......... .......... .......... 25% 45.1M 1s
    ##  20050K .......... .......... .......... .......... .......... 25% 95.6M 1s
    ##  20100K .......... .......... .......... .......... .......... 25% 33.5M 1s
    ##  20150K .......... .......... .......... .......... .......... 25% 96.5M 1s
    ##  20200K .......... .......... .......... .......... .......... 25% 67.5M 1s
    ##  20250K .......... .......... .......... .......... .......... 25%  109M 1s
    ##  20300K .......... .......... .......... .......... .......... 25% 30.8M 1s
    ##  20350K .......... .......... .......... .......... .......... 25%  167M 1s
    ##  20400K .......... .......... .......... .......... .......... 25%  136M 1s
    ##  20450K .......... .......... .......... .......... .......... 25%  120M 1s
    ##  20500K .......... .......... .......... .......... .......... 25% 91.0M 1s
    ##  20550K .......... .......... .......... .......... .......... 25% 77.0M 1s
    ##  20600K .......... .......... .......... .......... .......... 25% 94.9M 1s
    ##  20650K .......... .......... .......... .......... .......... 25%  121M 1s
    ##  20700K .......... .......... .......... .......... .......... 25% 76.7M 1s
    ##  20750K .......... .......... .......... .......... .......... 26% 36.7M 1s
    ##  20800K .......... .......... .......... .......... .......... 26% 23.0M 1s
    ##  20850K .......... .......... .......... .......... .......... 26% 22.3M 1s
    ##  20900K .......... .......... .......... .......... .......... 26% 85.8M 1s
    ##  20950K .......... .......... .......... .......... .......... 26% 98.0M 1s
    ##  21000K .......... .......... .......... .......... .......... 26% 74.3M 1s
    ##  21050K .......... .......... .......... .......... .......... 26% 77.7M 1s
    ##  21100K .......... .......... .......... .......... .......... 26% 90.0M 1s
    ##  21150K .......... .......... .......... .......... .......... 26%  142M 1s
    ##  21200K .......... .......... .......... .......... .......... 26%  137M 1s
    ##  21250K .......... .......... .......... .......... .......... 26%  178M 1s
    ##  21300K .......... .......... .......... .......... .......... 26%  122M 1s
    ##  21350K .......... .......... .......... .......... .......... 26%  160M 1s
    ##  21400K .......... .......... .......... .......... .......... 26%  105M 1s
    ##  21450K .......... .......... .......... .......... .......... 26%  128M 1s
    ##  21500K .......... .......... .......... .......... .......... 26%  116M 1s
    ##  21550K .......... .......... .......... .......... .......... 27% 78.7M 1s
    ##  21600K .......... .......... .......... .......... .......... 27%  131M 1s
    ##  21650K .......... .......... .......... .......... .......... 27% 58.5M 1s
    ##  21700K .......... .......... .......... .......... .......... 27%  127M 1s
    ##  21750K .......... .......... .......... .......... .......... 27%  130M 1s
    ##  21800K .......... .......... .......... .......... .......... 27%  144M 1s
    ##  21850K .......... .......... .......... .......... .......... 27% 71.8M 1s
    ##  21900K .......... .......... .......... .......... .......... 27%  118M 1s
    ##  21950K .......... .......... .......... .......... .......... 27%  124M 1s
    ##  22000K .......... .......... .......... .......... .......... 27% 44.5M 1s
    ##  22050K .......... .......... .......... .......... .......... 27% 97.2M 1s
    ##  22100K .......... .......... .......... .......... .......... 27% 30.4M 1s
    ##  22150K .......... .......... .......... .......... .......... 27% 31.7M 1s
    ##  22200K .......... .......... .......... .......... .......... 27%  105M 1s
    ##  22250K .......... .......... .......... .......... .......... 27%  122M 1s
    ##  22300K .......... .......... .......... .......... .......... 27%  133M 1s
    ##  22350K .......... .......... .......... .......... .......... 28%  114M 1s
    ##  22400K .......... .......... .......... .......... .......... 28%  125M 1s
    ##  22450K .......... .......... .......... .......... .......... 28% 19.1M 1s
    ##  22500K .......... .......... .......... .......... .......... 28%  140M 1s
    ##  22550K .......... .......... .......... .......... .......... 28%  139M 1s
    ##  22600K .......... .......... .......... .......... .......... 28%  106M 1s
    ##  22650K .......... .......... .......... .......... .......... 28%  126M 1s
    ##  22700K .......... .......... .......... .......... .......... 28% 51.5M 1s
    ##  22750K .......... .......... .......... .......... .......... 28% 54.7M 1s
    ##  22800K .......... .......... .......... .......... .......... 28% 58.3M 1s
    ##  22850K .......... .......... .......... .......... .......... 28% 40.5M 1s
    ##  22900K .......... .......... .......... .......... .......... 28%  124M 1s
    ##  22950K .......... .......... .......... .......... .......... 28%  131M 1s
    ##  23000K .......... .......... .......... .......... .......... 28% 53.1M 1s
    ##  23050K .......... .......... .......... .......... .......... 28%  165M 1s
    ##  23100K .......... .......... .......... .......... .......... 28%  161M 1s
    ##  23150K .......... .......... .......... .......... .......... 29%  120M 1s
    ##  23200K .......... .......... .......... .......... .......... 29%  140M 1s
    ##  23250K .......... .......... .......... .......... .......... 29%  116M 1s
    ##  23300K .......... .......... .......... .......... .......... 29% 86.9M 1s
    ##  23350K .......... .......... .......... .......... .......... 29%  117M 1s
    ##  23400K .......... .......... .......... .......... .......... 29%  103M 1s
    ##  23450K .......... .......... .......... .......... .......... 29% 37.1M 1s
    ##  23500K .......... .......... .......... .......... .......... 29% 62.9M 1s
    ##  23550K .......... .......... .......... .......... .......... 29% 96.5M 1s
    ##  23600K .......... .......... .......... .......... .......... 29% 59.5M 1s
    ##  23650K .......... .......... .......... .......... .......... 29% 45.5M 1s
    ##  23700K .......... .......... .......... .......... .......... 29%  104M 1s
    ##  23750K .......... .......... .......... .......... .......... 29%  106M 1s
    ##  23800K .......... .......... .......... .......... .......... 29% 98.6M 1s
    ##  23850K .......... .......... .......... .......... .......... 29% 87.3M 1s
    ##  23900K .......... .......... .......... .......... .......... 29% 41.0M 1s
    ##  23950K .......... .......... .......... .......... .......... 30% 96.0M 1s
    ##  24000K .......... .......... .......... .......... .......... 30%  142M 1s
    ##  24050K .......... .......... .......... .......... .......... 30% 60.0M 1s
    ##  24100K .......... .......... .......... .......... .......... 30% 93.7M 1s
    ##  24150K .......... .......... .......... .......... .......... 30% 54.1M 1s
    ##  24200K .......... .......... .......... .......... .......... 30% 92.3M 1s
    ##  24250K .......... .......... .......... .......... .......... 30%  142M 1s
    ##  24300K .......... .......... .......... .......... .......... 30% 58.0M 1s
    ##  24350K .......... .......... .......... .......... .......... 30%  140M 1s
    ##  24400K .......... .......... .......... .......... .......... 30%  147M 1s
    ##  24450K .......... .......... .......... .......... .......... 30% 52.6M 1s
    ##  24500K .......... .......... .......... .......... .......... 30%  101M 1s
    ##  24550K .......... .......... .......... .......... .......... 30% 82.7M 1s
    ##  24600K .......... .......... .......... .......... .......... 30% 67.7M 1s
    ##  24650K .......... .......... .......... .......... .......... 30%  118M 1s
    ##  24700K .......... .......... .......... .......... .......... 30% 45.1M 1s
    ##  24750K .......... .......... .......... .......... .......... 31%  134M 1s
    ##  24800K .......... .......... .......... .......... .......... 31% 53.5M 1s
    ##  24850K .......... .......... .......... .......... .......... 31%  135M 1s
    ##  24900K .......... .......... .......... .......... .......... 31%  130M 1s
    ##  24950K .......... .......... .......... .......... .......... 31%  115M 1s
    ##  25000K .......... .......... .......... .......... .......... 31% 75.7M 1s
    ##  25050K .......... .......... .......... .......... .......... 31% 60.1M 1s
    ##  25100K .......... .......... .......... .......... .......... 31% 67.0M 1s
    ##  25150K .......... .......... .......... .......... .......... 31%  161M 1s
    ##  25200K .......... .......... .......... .......... .......... 31% 57.3M 1s
    ##  25250K .......... .......... .......... .......... .......... 31%  146M 1s
    ##  25300K .......... .......... .......... .......... .......... 31% 94.3M 1s
    ##  25350K .......... .......... .......... .......... .......... 31%  121M 1s
    ##  25400K .......... .......... .......... .......... .......... 31% 73.7M 1s
    ##  25450K .......... .......... .......... .......... .......... 31% 55.5M 1s
    ##  25500K .......... .......... .......... .......... .......... 31%  104M 1s
    ##  25550K .......... .......... .......... .......... .......... 32%  178M 1s
    ##  25600K .......... .......... .......... .......... .......... 32% 52.1M 1s
    ##  25650K .......... .......... .......... .......... .......... 32% 92.9M 1s
    ##  25700K .......... .......... .......... .......... .......... 32% 94.3M 1s
    ##  25750K .......... .......... .......... .......... .......... 32%  134M 1s
    ##  25800K .......... .......... .......... .......... .......... 32%  131M 1s
    ##  25850K .......... .......... .......... .......... .......... 32% 59.9M 1s
    ##  25900K .......... .......... .......... .......... .......... 32% 13.3M 1s
    ##  25950K .......... .......... .......... .......... .......... 32%  130M 1s
    ##  26000K .......... .......... .......... .......... .......... 32%  139M 1s
    ##  26050K .......... .......... .......... .......... .......... 32%  155M 1s
    ##  26100K .......... .......... .......... .......... .......... 32%  145M 1s
    ##  26150K .......... .......... .......... .......... .......... 32%  154M 1s
    ##  26200K .......... .......... .......... .......... .......... 32%  115M 1s
    ##  26250K .......... .......... .......... .......... .......... 32%  130M 1s
    ##  26300K .......... .......... .......... .......... .......... 32% 16.9M 1s
    ##  26350K .......... .......... .......... .......... .......... 33% 26.6M 1s
    ##  26400K .......... .......... .......... .......... .......... 33%  100M 1s
    ##  26450K .......... .......... .......... .......... .......... 33% 29.9M 1s
    ##  26500K .......... .......... .......... .......... .......... 33%  140M 1s
    ##  26550K .......... .......... .......... .......... .......... 33% 29.2M 1s
    ##  26600K .......... .......... .......... .......... .......... 33%  119M 1s
    ##  26650K .......... .......... .......... .......... .......... 33%  169M 1s
    ##  26700K .......... .......... .......... .......... .......... 33%  144M 1s
    ##  26750K .......... .......... .......... .......... .......... 33%  163M 1s
    ##  26800K .......... .......... .......... .......... .......... 33%  130M 1s
    ##  26850K .......... .......... .......... .......... .......... 33%  127M 1s
    ##  26900K .......... .......... .......... .......... .......... 33%  116M 1s
    ##  26950K .......... .......... .......... .......... .......... 33%  146M 1s
    ##  27000K .......... .......... .......... .......... .......... 33%  144M 1s
    ##  27050K .......... .......... .......... .......... .......... 33% 30.4M 1s
    ##  27100K .......... .......... .......... .......... .......... 33% 43.7M 1s
    ##  27150K .......... .......... .......... .......... .......... 34% 88.9M 1s
    ##  27200K .......... .......... .......... .......... .......... 34% 58.4M 1s
    ##  27250K .......... .......... .......... .......... .......... 34% 36.4M 1s
    ##  27300K .......... .......... .......... .......... .......... 34%  114M 1s
    ##  27350K .......... .......... .......... .......... .......... 34%  115M 1s
    ##  27400K .......... .......... .......... .......... .......... 34% 84.8M 1s
    ##  27450K .......... .......... .......... .......... .......... 34%  126M 1s
    ##  27500K .......... .......... .......... .......... .......... 34% 31.6M 1s
    ##  27550K .......... .......... .......... .......... .......... 34%  133M 1s
    ##  27600K .......... .......... .......... .......... .......... 34% 50.8M 1s
    ##  27650K .......... .......... .......... .......... .......... 34% 41.9M 1s
    ##  27700K .......... .......... .......... .......... .......... 34% 90.7M 1s
    ##  27750K .......... .......... .......... .......... .......... 34%  128M 1s
    ##  27800K .......... .......... .......... .......... .......... 34% 76.7M 1s
    ##  27850K .......... .......... .......... .......... .......... 34% 36.4M 1s
    ##  27900K .......... .......... .......... .......... .......... 34%  106M 1s
    ##  27950K .......... .......... .......... .......... .......... 35% 77.6M 1s
    ##  28000K .......... .......... .......... .......... .......... 35% 85.1M 1s
    ##  28050K .......... .......... .......... .......... .......... 35% 42.0M 1s
    ##  28100K .......... .......... .......... .......... .......... 35%  104M 1s
    ##  28150K .......... .......... .......... .......... .......... 35%  157M 1s
    ##  28200K .......... .......... .......... .......... .......... 35% 85.9M 1s
    ##  28250K .......... .......... .......... .......... .......... 35% 43.9M 1s
    ##  28300K .......... .......... .......... .......... .......... 35% 75.4M 1s
    ##  28350K .......... .......... .......... .......... .......... 35% 68.5M 1s
    ##  28400K .......... .......... .......... .......... .......... 35% 86.7M 1s
    ##  28450K .......... .......... .......... .......... .......... 35% 45.1M 1s
    ##  28500K .......... .......... .......... .......... .......... 35% 44.6M 1s
    ##  28550K .......... .......... .......... .......... .......... 35%  110M 1s
    ##  28600K .......... .......... .......... .......... .......... 35% 37.4M 1s
    ##  28650K .......... .......... .......... .......... .......... 35%  161M 1s
    ##  28700K .......... .......... .......... .......... .......... 35% 92.4M 1s
    ##  28750K .......... .......... .......... .......... .......... 36%  123M 1s
    ##  28800K .......... .......... .......... .......... .......... 36%  131M 1s
    ##  28850K .......... .......... .......... .......... .......... 36% 99.0M 1s
    ##  28900K .......... .......... .......... .......... .......... 36% 57.5M 1s
    ##  28950K .......... .......... .......... .......... .......... 36%  104M 1s
    ##  29000K .......... .......... .......... .......... .......... 36% 35.7M 1s
    ##  29050K .......... .......... .......... .......... .......... 36%  144M 1s
    ##  29100K .......... .......... .......... .......... .......... 36% 27.8M 1s
    ##  29150K .......... .......... .......... .......... .......... 36% 92.2M 1s
    ##  29200K .......... .......... .......... .......... .......... 36%  104M 1s
    ##  29250K .......... .......... .......... .......... .......... 36%  140M 1s
    ##  29300K .......... .......... .......... .......... .......... 36%  123M 1s
    ##  29350K .......... .......... .......... .......... .......... 36% 8.17M 1s
    ##  29400K .......... .......... .......... .......... .......... 36% 24.1M 1s
    ##  29450K .......... .......... .......... .......... .......... 36% 29.2M 1s
    ##  29500K .......... .......... .......... .......... .......... 36%  134M 1s
    ##  29550K .......... .......... .......... .......... .......... 37% 13.8M 1s
    ##  29600K .......... .......... .......... .......... .......... 37%  116M 1s
    ##  29650K .......... .......... .......... .......... .......... 37%  166M 1s
    ##  29700K .......... .......... .......... .......... .......... 37%  105M 1s
    ##  29750K .......... .......... .......... .......... .......... 37%  144M 1s
    ##  29800K .......... .......... .......... .......... .......... 37%  140M 1s
    ##  29850K .......... .......... .......... .......... .......... 37%  173M 1s
    ##  29900K .......... .......... .......... .......... .......... 37%  126M 1s
    ##  29950K .......... .......... .......... .......... .......... 37%  119M 1s
    ##  30000K .......... .......... .......... .......... .......... 37% 67.9M 1s
    ##  30050K .......... .......... .......... .......... .......... 37% 36.5M 1s
    ##  30100K .......... .......... .......... .......... .......... 37% 22.5M 1s
    ##  30150K .......... .......... .......... .......... .......... 37%  111M 1s
    ##  30200K .......... .......... .......... .......... .......... 37% 73.2M 1s
    ##  30250K .......... .......... .......... .......... .......... 37% 64.8M 1s
    ##  30300K .......... .......... .......... .......... .......... 37% 26.5M 1s
    ##  30350K .......... .......... .......... .......... .......... 38%  143M 1s
    ##  30400K .......... .......... .......... .......... .......... 38%  122M 1s
    ##  30450K .......... .......... .......... .......... .......... 38%  130M 1s
    ##  30500K .......... .......... .......... .......... .......... 38%  113M 1s
    ##  30550K .......... .......... .......... .......... .......... 38%  126M 1s
    ##  30600K .......... .......... .......... .......... .......... 38%  121M 1s
    ##  30650K .......... .......... .......... .......... .......... 38%  131M 1s
    ##  30700K .......... .......... .......... .......... .......... 38%  132M 1s
    ##  30750K .......... .......... .......... .......... .......... 38%  170M 1s
    ##  30800K .......... .......... .......... .......... .......... 38% 42.5M 1s
    ##  30850K .......... .......... .......... .......... .......... 38% 20.5M 1s
    ##  30900K .......... .......... .......... .......... .......... 38% 73.7M 1s
    ##  30950K .......... .......... .......... .......... .......... 38%  108M 1s
    ##  31000K .......... .......... .......... .......... .......... 38% 62.9M 1s
    ##  31050K .......... .......... .......... .......... .......... 38%  113M 1s
    ##  31100K .......... .......... .......... .......... .......... 38%  103M 1s
    ##  31150K .......... .......... .......... .......... .......... 39% 41.0M 1s
    ##  31200K .......... .......... .......... .......... .......... 39% 87.0M 1s
    ##  31250K .......... .......... .......... .......... .......... 39%  122M 1s
    ##  31300K .......... .......... .......... .......... .......... 39% 27.2M 1s
    ##  31350K .......... .......... .......... .......... .......... 39%  141M 1s
    ##  31400K .......... .......... .......... .......... .......... 39% 31.4M 1s
    ##  31450K .......... .......... .......... .......... .......... 39%  119M 1s
    ##  31500K .......... .......... .......... .......... .......... 39% 93.7M 1s
    ##  31550K .......... .......... .......... .......... .......... 39%  136M 1s
    ##  31600K .......... .......... .......... .......... .......... 39% 32.1M 1s
    ##  31650K .......... .......... .......... .......... .......... 39%  136M 1s
    ##  31700K .......... .......... .......... .......... .......... 39% 38.7M 1s
    ##  31750K .......... .......... .......... .......... .......... 39%  112M 1s
    ##  31800K .......... .......... .......... .......... .......... 39%  127M 1s
    ##  31850K .......... .......... .......... .......... .......... 39%  139M 1s
    ##  31900K .......... .......... .......... .......... .......... 39%  105M 1s
    ##  31950K .......... .......... .......... .......... .......... 40% 97.6M 1s
    ##  32000K .......... .......... .......... .......... .......... 40%  139M 1s
    ##  32050K .......... .......... .......... .......... .......... 40%  126M 1s
    ##  32100K .......... .......... .......... .......... .......... 40% 57.1M 1s
    ##  32150K .......... .......... .......... .......... .......... 40% 61.0M 1s
    ##  32200K .......... .......... .......... .......... .......... 40% 39.1M 1s
    ##  32250K .......... .......... .......... .......... .......... 40%  107M 1s
    ##  32300K .......... .......... .......... .......... .......... 40%  110M 1s
    ##  32350K .......... .......... .......... .......... .......... 40% 33.2M 1s
    ##  32400K .......... .......... .......... .......... .......... 40% 94.0M 1s
    ##  32450K .......... .......... .......... .......... .......... 40%  154M 1s
    ##  32500K .......... .......... .......... .......... .......... 40% 48.5M 1s
    ##  32550K .......... .......... .......... .......... .......... 40%  108M 1s
    ##  32600K .......... .......... .......... .......... .......... 40% 43.2M 1s
    ##  32650K .......... .......... .......... .......... .......... 40% 44.7M 1s
    ##  32700K .......... .......... .......... .......... .......... 40%  108M 1s
    ##  32750K .......... .......... .......... .......... .......... 41% 94.6M 1s
    ##  32800K .......... .......... .......... .......... .......... 41% 39.7M 1s
    ##  32850K .......... .......... .......... .......... .......... 41%  123M 1s
    ##  32900K .......... .......... .......... .......... .......... 41% 89.0M 1s
    ##  32950K .......... .......... .......... .......... .......... 41% 67.5M 1s
    ##  33000K .......... .......... .......... .......... .......... 41%  105M 1s
    ##  33050K .......... .......... .......... .......... .......... 41% 44.3M 1s
    ##  33100K .......... .......... .......... .......... .......... 41%  110M 1s
    ##  33150K .......... .......... .......... .......... .......... 41%  113M 1s
    ##  33200K .......... .......... .......... .......... .......... 41% 40.8M 1s
    ##  33250K .......... .......... .......... .......... .......... 41%  125M 1s
    ##  33300K .......... .......... .......... .......... .......... 41% 62.2M 1s
    ##  33350K .......... .......... .......... .......... .......... 41% 37.1M 1s
    ##  33400K .......... .......... .......... .......... .......... 41% 85.5M 1s
    ##  33450K .......... .......... .......... .......... .......... 41% 99.2M 1s
    ##  33500K .......... .......... .......... .......... .......... 41%  104M 1s
    ##  33550K .......... .......... .......... .......... .......... 42% 60.9M 1s
    ##  33600K .......... .......... .......... .......... .......... 42% 93.9M 1s
    ##  33650K .......... .......... .......... .......... .......... 42%  111M 1s
    ##  33700K .......... .......... .......... .......... .......... 42% 45.0M 1s
    ##  33750K .......... .......... .......... .......... .......... 42% 99.7M 1s
    ##  33800K .......... .......... .......... .......... .......... 42%  120M 1s
    ##  33850K .......... .......... .......... .......... .......... 42% 57.3M 1s
    ##  33900K .......... .......... .......... .......... .......... 42% 58.6M 1s
    ##  33950K .......... .......... .......... .......... .......... 42%  119M 1s
    ##  34000K .......... .......... .......... .......... .......... 42% 92.6M 1s
    ##  34050K .......... .......... .......... .......... .......... 42% 90.0M 1s
    ##  34100K .......... .......... .......... .......... .......... 42%  124M 1s
    ##  34150K .......... .......... .......... .......... .......... 42% 39.1M 1s
    ##  34200K .......... .......... .......... .......... .......... 42% 83.2M 1s
    ##  34250K .......... .......... .......... .......... .......... 42%  114M 1s
    ##  34300K .......... .......... .......... .......... .......... 42%  142M 1s
    ##  34350K .......... .......... .......... .......... .......... 43% 34.7M 1s
    ##  34400K .......... .......... .......... .......... .......... 43%  103M 1s
    ##  34450K .......... .......... .......... .......... .......... 43% 48.6M 1s
    ##  34500K .......... .......... .......... .......... .......... 43% 77.6M 1s
    ##  34550K .......... .......... .......... .......... .......... 43%  122M 1s
    ##  34600K .......... .......... .......... .......... .......... 43%  106M 1s
    ##  34650K .......... .......... .......... .......... .......... 43% 42.1M 1s
    ##  34700K .......... .......... .......... .......... .......... 43%  107M 1s
    ##  34750K .......... .......... .......... .......... .......... 43% 56.3M 1s
    ##  34800K .......... .......... .......... .......... .......... 43%  104M 1s
    ##  34850K .......... .......... .......... .......... .......... 43%  143M 1s
    ##  34900K .......... .......... .......... .......... .......... 43% 47.6M 1s
    ##  34950K .......... .......... .......... .......... .......... 43%  103M 1s
    ##  35000K .......... .......... .......... .......... .......... 43%  129M 1s
    ##  35050K .......... .......... .......... .......... .......... 43%  105M 1s
    ##  35100K .......... .......... .......... .......... .......... 43% 39.0M 1s
    ##  35150K .......... .......... .......... .......... .......... 44%  122M 1s
    ##  35200K .......... .......... .......... .......... .......... 44%  116M 1s
    ##  35250K .......... .......... .......... .......... .......... 44% 71.8M 1s
    ##  35300K .......... .......... .......... .......... .......... 44%  121M 1s
    ##  35350K .......... .......... .......... .......... .......... 44% 36.8M 1s
    ##  35400K .......... .......... .......... .......... .......... 44% 96.2M 1s
    ##  35450K .......... .......... .......... .......... .......... 44%  104M 1s
    ##  35500K .......... .......... .......... .......... .......... 44%  145M 1s
    ##  35550K .......... .......... .......... .......... .......... 44% 50.9M 1s
    ##  35600K .......... .......... .......... .......... .......... 44%  139M 1s
    ##  35650K .......... .......... .......... .......... .......... 44%  115M 1s
    ##  35700K .......... .......... .......... .......... .......... 44% 64.5M 1s
    ##  35750K .......... .......... .......... .......... .......... 44% 59.5M 1s
    ##  35800K .......... .......... .......... .......... .......... 44% 96.5M 1s
    ##  35850K .......... .......... .......... .......... .......... 44%  107M 1s
    ##  35900K .......... .......... .......... .......... .......... 44% 69.7M 1s
    ##  35950K .......... .......... .......... .......... .......... 45%  157M 1s
    ##  36000K .......... .......... .......... .......... .......... 45% 61.3M 1s
    ##  36050K .......... .......... .......... .......... .......... 45% 85.6M 1s
    ##  36100K .......... .......... .......... .......... .......... 45% 78.7M 1s
    ##  36150K .......... .......... .......... .......... .......... 45% 83.8M 1s
    ##  36200K .......... .......... .......... .......... .......... 45%  132M 1s
    ##  36250K .......... .......... .......... .......... .......... 45%  109M 1s
    ##  36300K .......... .......... .......... .......... .......... 45% 45.5M 1s
    ##  36350K .......... .......... .......... .......... .......... 45% 93.8M 1s
    ##  36400K .......... .......... .......... .......... .......... 45%  111M 1s
    ##  36450K .......... .......... .......... .......... .......... 45% 67.8M 1s
    ##  36500K .......... .......... .......... .......... .......... 45% 82.5M 1s
    ##  36550K .......... .......... .......... .......... .......... 45%  116M 1s
    ##  36600K .......... .......... .......... .......... .......... 45% 72.0M 1s
    ##  36650K .......... .......... .......... .......... .......... 45% 91.1M 1s
    ##  36700K .......... .......... .......... .......... .......... 45%  110M 1s
    ##  36750K .......... .......... .......... .......... .......... 46% 58.1M 1s
    ##  36800K .......... .......... .......... .......... .......... 46% 90.0M 1s
    ##  36850K .......... .......... .......... .......... .......... 46%  138M 1s
    ##  36900K .......... .......... .......... .......... .......... 46%  118M 1s
    ##  36950K .......... .......... .......... .......... .......... 46% 50.2M 1s
    ##  37000K .......... .......... .......... .......... .......... 46%  105M 1s
    ##  37050K .......... .......... .......... .......... .......... 46% 61.4M 1s
    ##  37100K .......... .......... .......... .......... .......... 46% 84.9M 1s
    ##  37150K .......... .......... .......... .......... .......... 46%  129M 1s
    ##  37200K .......... .......... .......... .......... .......... 46% 85.8M 1s
    ##  37250K .......... .......... .......... .......... .......... 46%  104M 1s
    ##  37300K .......... .......... .......... .......... .......... 46%  134M 1s
    ##  37350K .......... .......... .......... .......... .......... 46% 53.1M 1s
    ##  37400K .......... .......... .......... .......... .......... 46% 96.0M 1s
    ##  37450K .......... .......... .......... .......... .......... 46%  119M 1s
    ##  37500K .......... .......... .......... .......... .......... 46% 77.5M 1s
    ##  37550K .......... .......... .......... .......... .......... 47%  132M 1s
    ##  37600K .......... .......... .......... .......... .......... 47%  139M 1s
    ##  37650K .......... .......... .......... .......... .......... 47% 54.6M 1s
    ##  37700K .......... .......... .......... .......... .......... 47%  150M 1s
    ##  37750K .......... .......... .......... .......... .......... 47%  102M 1s
    ##  37800K .......... .......... .......... .......... .......... 47% 67.1M 1s
    ##  37850K .......... .......... .......... .......... .......... 47%  111M 1s
    ##  37900K .......... .......... .......... .......... .......... 47%  123M 1s
    ##  37950K .......... .......... .......... .......... .......... 47% 59.4M 1s
    ##  38000K .......... .......... .......... .......... .......... 47%  128M 1s
    ##  38050K .......... .......... .......... .......... .......... 47% 72.7M 1s
    ##  38100K .......... .......... .......... .......... .......... 47% 84.1M 1s
    ##  38150K .......... .......... .......... .......... .......... 47%  117M 1s
    ##  38200K .......... .......... .......... .......... .......... 47%  108M 1s
    ##  38250K .......... .......... .......... .......... .......... 47%  104M 1s
    ##  38300K .......... .......... .......... .......... .......... 47%  126M 1s
    ##  38350K .......... .......... .......... .......... .......... 48%  125M 1s
    ##  38400K .......... .......... .......... .......... .......... 48% 75.2M 1s
    ##  38450K .......... .......... .......... .......... .......... 48% 93.7M 1s
    ##  38500K .......... .......... .......... .......... .......... 48%  112M 1s
    ##  38550K .......... .......... .......... .......... .......... 48% 68.9M 1s
    ##  38600K .......... .......... .......... .......... .......... 48%  144M 1s
    ##  38650K .......... .......... .......... .......... .......... 48%  127M 1s
    ##  38700K .......... .......... .......... .......... .......... 48% 98.9M 1s
    ##  38750K .......... .......... .......... .......... .......... 48%  121M 1s
    ##  38800K .......... .......... .......... .......... .......... 48%  118M 1s
    ##  38850K .......... .......... .......... .......... .......... 48% 31.2M 1s
    ##  38900K .......... .......... .......... .......... .......... 48%  147M 1s
    ##  38950K .......... .......... .......... .......... .......... 48% 61.2M 1s
    ##  39000K .......... .......... .......... .......... .......... 48% 49.9M 1s
    ##  39050K .......... .......... .......... .......... .......... 48%  126M 1s
    ##  39100K .......... .......... .......... .......... .......... 48%  120M 1s
    ##  39150K .......... .......... .......... .......... .......... 49% 76.3M 1s
    ##  39200K .......... .......... .......... .......... .......... 49%  125M 1s
    ##  39250K .......... .......... .......... .......... .......... 49%  152M 1s
    ##  39300K .......... .......... .......... .......... .......... 49% 57.1M 1s
    ##  39350K .......... .......... .......... .......... .......... 49% 81.1M 1s
    ##  39400K .......... .......... .......... .......... .......... 49% 32.1M 1s
    ##  39450K .......... .......... .......... .......... .......... 49%  138M 1s
    ##  39500K .......... .......... .......... .......... .......... 49%  167M 1s
    ##  39550K .......... .......... .......... .......... .......... 49% 76.4M 1s
    ##  39600K .......... .......... .......... .......... .......... 49% 68.1M 1s
    ##  39650K .......... .......... .......... .......... .......... 49% 57.2M 1s
    ##  39700K .......... .......... .......... .......... .......... 49%  141M 1s
    ##  39750K .......... .......... .......... .......... .......... 49%  150M 1s
    ##  39800K .......... .......... .......... .......... .......... 49% 36.9M 1s
    ##  39850K .......... .......... .......... .......... .......... 49%  114M 1s
    ##  39900K .......... .......... .......... .......... .......... 49%  133M 1s
    ##  39950K .......... .......... .......... .......... .......... 50% 95.8M 1s
    ##  40000K .......... .......... .......... .......... .......... 50% 64.3M 1s
    ##  40050K .......... .......... .......... .......... .......... 50%  110M 1s
    ##  40100K .......... .......... .......... .......... .......... 50% 68.1M 1s
    ##  40150K .......... .......... .......... .......... .......... 50% 49.7M 1s
    ##  40200K .......... .......... .......... .......... .......... 50%  105M 1s
    ##  40250K .......... .......... .......... .......... .......... 50% 5.79M 1s
    ##  40300K .......... .......... .......... .......... .......... 50%  118M 1s
    ##  40350K .......... .......... .......... .......... .......... 50%  151M 1s
    ##  40400K .......... .......... .......... .......... .......... 50%  161M 1s
    ##  40450K .......... .......... .......... .......... .......... 50%  168M 1s
    ##  40500K .......... .......... .......... .......... .......... 50%  151M 1s
    ##  40550K .......... .......... .......... .......... .......... 50%  128M 1s
    ##  40600K .......... .......... .......... .......... .......... 50%  133M 1s
    ##  40650K .......... .......... .......... .......... .......... 50%  167M 1s
    ##  40700K .......... .......... .......... .......... .......... 50% 23.2M 1s
    ##  40750K .......... .......... .......... .......... .......... 51% 26.2M 1s
    ##  40800K .......... .......... .......... .......... .......... 51%  133M 1s
    ##  40850K .......... .......... .......... .......... .......... 51% 37.8M 1s
    ##  40900K .......... .......... .......... .......... .......... 51%  134M 1s
    ##  40950K .......... .......... .......... .......... .......... 51%  107M 1s
    ##  41000K .......... .......... .......... .......... .......... 51%  124M 1s
    ##  41050K .......... .......... .......... .......... .......... 51%  136M 1s
    ##  41100K .......... .......... .......... .......... .......... 51%  129M 1s
    ##  41150K .......... .......... .......... .......... .......... 51%  166M 1s
    ##  41200K .......... .......... .......... .......... .......... 51%  124M 1s
    ##  41250K .......... .......... .......... .......... .......... 51%  145M 1s
    ##  41300K .......... .......... .......... .......... .......... 51% 26.1M 1s
    ##  41350K .......... .......... .......... .......... .......... 51% 39.6M 1s
    ##  41400K .......... .......... .......... .......... .......... 51% 45.7M 1s
    ##  41450K .......... .......... .......... .......... .......... 51% 25.3M 1s
    ##  41500K .......... .......... .......... .......... .......... 51%  132M 1s
    ##  41550K .......... .......... .......... .......... .......... 52%  176M 1s
    ##  41600K .......... .......... .......... .......... .......... 52%  142M 1s
    ##  41650K .......... .......... .......... .......... .......... 52%  151M 1s
    ##  41700K .......... .......... .......... .......... .......... 52% 92.3M 1s
    ##  41750K .......... .......... .......... .......... .......... 52%  133M 1s
    ##  41800K .......... .......... .......... .......... .......... 52%  109M 1s
    ##  41850K .......... .......... .......... .......... .......... 52%  190M 1s
    ##  41900K .......... .......... .......... .......... .......... 52% 94.1M 1s
    ##  41950K .......... .......... .......... .......... .......... 52%  124M 1s
    ##  42000K .......... .......... .......... .......... .......... 52%  134M 1s
    ##  42050K .......... .......... .......... .......... .......... 52%  108M 1s
    ##  42100K .......... .......... .......... .......... .......... 52% 42.3M 1s
    ##  42150K .......... .......... .......... .......... .......... 52% 24.0M 1s
    ##  42200K .......... .......... .......... .......... .......... 52% 33.8M 1s
    ##  42250K .......... .......... .......... .......... .......... 52% 82.0M 1s
    ##  42300K .......... .......... .......... .......... .......... 52%  136M 1s
    ##  42350K .......... .......... .......... .......... .......... 53%  160M 1s
    ##  42400K .......... .......... .......... .......... .......... 53%  108M 1s
    ##  42450K .......... .......... .......... .......... .......... 53%  107M 1s
    ##  42500K .......... .......... .......... .......... .......... 53%  146M 1s
    ##  42550K .......... .......... .......... .......... .......... 53%  129M 1s
    ##  42600K .......... .......... .......... .......... .......... 53%  110M 1s
    ##  42650K .......... .......... .......... .......... .......... 53%  126M 1s
    ##  42700K .......... .......... .......... .......... .......... 53%  140M 1s
    ##  42750K .......... .......... .......... .......... .......... 53% 74.4M 1s
    ##  42800K .......... .......... .......... .......... .......... 53%  108M 1s
    ##  42850K .......... .......... .......... .......... .......... 53% 52.8M 1s
    ##  42900K .......... .......... .......... .......... .......... 53%  107M 1s
    ##  42950K .......... .......... .......... .......... .......... 53% 34.7M 1s
    ##  43000K .......... .......... .......... .......... .......... 53% 11.0M 1s
    ##  43050K .......... .......... .......... .......... .......... 53% 78.1M 1s
    ##  43100K .......... .......... .......... .......... .......... 53% 91.3M 1s
    ##  43150K .......... .......... .......... .......... .......... 54%  103M 1s
    ##  43200K .......... .......... .......... .......... .......... 54% 91.5M 1s
    ##  43250K .......... .......... .......... .......... .......... 54% 97.1M 1s
    ##  43300K .......... .......... .......... .......... .......... 54%  103M 1s
    ##  43350K .......... .......... .......... .......... .......... 54% 91.3M 1s
    ##  43400K .......... .......... .......... .......... .......... 54% 95.7M 1s
    ##  43450K .......... .......... .......... .......... .......... 54%  100M 1s
    ##  43500K .......... .......... .......... .......... .......... 54%  112M 1s
    ##  43550K .......... .......... .......... .......... .......... 54%  130M 1s
    ##  43600K .......... .......... .......... .......... .......... 54% 41.8M 1s
    ##  43650K .......... .......... .......... .......... .......... 54% 77.9M 1s
    ##  43700K .......... .......... .......... .......... .......... 54%  112M 1s
    ##  43750K .......... .......... .......... .......... .......... 54% 46.9M 1s
    ##  43800K .......... .......... .......... .......... .......... 54% 29.3M 1s
    ##  43850K .......... .......... .......... .......... .......... 54% 72.5M 1s
    ##  43900K .......... .......... .......... .......... .......... 54% 85.7M 1s
    ##  43950K .......... .......... .......... .......... .......... 55%  109M 1s
    ##  44000K .......... .......... .......... .......... .......... 55% 96.2M 1s
    ##  44050K .......... .......... .......... .......... .......... 55% 56.7M 1s
    ##  44100K .......... .......... .......... .......... .......... 55%  115M 1s
    ##  44150K .......... .......... .......... .......... .......... 55%  105M 1s
    ##  44200K .......... .......... .......... .......... .......... 55%  105M 1s
    ##  44250K .......... .......... .......... .......... .......... 55% 38.7M 1s
    ##  44300K .......... .......... .......... .......... .......... 55% 69.6M 1s
    ##  44350K .......... .......... .......... .......... .......... 55% 25.7M 1s
    ##  44400K .......... .......... .......... .......... .......... 55% 23.3M 1s
    ##  44450K .......... .......... .......... .......... .......... 55% 77.0M 1s
    ##  44500K .......... .......... .......... .......... .......... 55%  159M 1s
    ##  44550K .......... .......... .......... .......... .......... 55% 75.0M 1s
    ##  44600K .......... .......... .......... .......... .......... 55% 72.4M 1s
    ##  44650K .......... .......... .......... .......... .......... 55% 66.8M 1s
    ##  44700K .......... .......... .......... .......... .......... 55% 64.8M 1s
    ##  44750K .......... .......... .......... .......... .......... 56%  105M 1s
    ##  44800K .......... .......... .......... .......... .......... 56% 63.2M 1s
    ##  44850K .......... .......... .......... .......... .......... 56% 33.3M 1s
    ##  44900K .......... .......... .......... .......... .......... 56% 72.5M 1s
    ##  44950K .......... .......... .......... .......... .......... 56% 82.5M 1s
    ##  45000K .......... .......... .......... .......... .......... 56% 46.5M 1s
    ##  45050K .......... .......... .......... .......... .......... 56% 50.7M 1s
    ##  45100K .......... .......... .......... .......... .......... 56% 89.4M 1s
    ##  45150K .......... .......... .......... .......... .......... 56% 50.2M 1s
    ##  45200K .......... .......... .......... .......... .......... 56% 71.7M 1s
    ##  45250K .......... .......... .......... .......... .......... 56%  104M 1s
    ##  45300K .......... .......... .......... .......... .......... 56% 49.4M 1s
    ##  45350K .......... .......... .......... .......... .......... 56% 65.7M 1s
    ##  45400K .......... .......... .......... .......... .......... 56% 77.6M 1s
    ##  45450K .......... .......... .......... .......... .......... 56% 78.5M 1s
    ##  45500K .......... .......... .......... .......... .......... 56% 39.6M 1s
    ##  45550K .......... .......... .......... .......... .......... 57% 77.2M 1s
    ##  45600K .......... .......... .......... .......... .......... 57% 49.2M 1s
    ##  45650K .......... .......... .......... .......... .......... 57% 71.2M 1s
    ##  45700K .......... .......... .......... .......... .......... 57%  128M 0s
    ##  45750K .......... .......... .......... .......... .......... 57% 40.1M 0s
    ##  45800K .......... .......... .......... .......... .......... 57% 98.8M 0s
    ##  45850K .......... .......... .......... .......... .......... 57% 81.7M 0s
    ##  45900K .......... .......... .......... .......... .......... 57%  120M 0s
    ##  45950K .......... .......... .......... .......... .......... 57% 42.3M 0s
    ##  46000K .......... .......... .......... .......... .......... 57% 93.3M 0s
    ##  46050K .......... .......... .......... .......... .......... 57% 41.0M 0s
    ##  46100K .......... .......... .......... .......... .......... 57% 74.4M 0s
    ##  46150K .......... .......... .......... .......... .......... 57%  153M 0s
    ##  46200K .......... .......... .......... .......... .......... 57% 38.2M 0s
    ##  46250K .......... .......... .......... .......... .......... 57% 69.9M 0s
    ##  46300K .......... .......... .......... .......... .......... 57% 68.2M 0s
    ##  46350K .......... .......... .......... .......... .......... 58% 78.7M 0s
    ##  46400K .......... .......... .......... .......... .......... 58% 44.6M 0s
    ##  46450K .......... .......... .......... .......... .......... 58%  121M 0s
    ##  46500K .......... .......... .......... .......... .......... 58% 73.5M 0s
    ##  46550K .......... .......... .......... .......... .......... 58% 81.2M 0s
    ##  46600K .......... .......... .......... .......... .......... 58% 97.2M 0s
    ##  46650K .......... .......... .......... .......... .......... 58% 61.8M 0s
    ##  46700K .......... .......... .......... .......... .......... 58% 37.0M 0s
    ##  46750K .......... .......... .......... .......... .......... 58% 76.2M 0s
    ##  46800K .......... .......... .......... .......... .......... 58% 74.7M 0s
    ##  46850K .......... .......... .......... .......... .......... 58% 70.0M 0s
    ##  46900K .......... .......... .......... .......... .......... 58% 97.4M 0s
    ##  46950K .......... .......... .......... .......... .......... 58% 56.1M 0s
    ##  47000K .......... .......... .......... .......... .......... 58%  105M 0s
    ##  47050K .......... .......... .......... .......... .......... 58% 79.7M 0s
    ##  47100K .......... .......... .......... .......... .......... 58%  141M 0s
    ##  47150K .......... .......... .......... .......... .......... 59% 17.9M 0s
    ##  47200K .......... .......... .......... .......... .......... 59%  149M 0s
    ##  47250K .......... .......... .......... .......... .......... 59%  133M 0s
    ##  47300K .......... .......... .......... .......... .......... 59%  144M 0s
    ##  47350K .......... .......... .......... .......... .......... 59%  143M 0s
    ##  47400K .......... .......... .......... .......... .......... 59% 62.7M 0s
    ##  47450K .......... .......... .......... .......... .......... 59% 77.9M 0s
    ##  47500K .......... .......... .......... .......... .......... 59%  107M 0s
    ##  47550K .......... .......... .......... .......... .......... 59% 66.8M 0s
    ##  47600K .......... .......... .......... .......... .......... 59%  145M 0s
    ##  47650K .......... .......... .......... .......... .......... 59%  126M 0s
    ##  47700K .......... .......... .......... .......... .......... 59% 80.1M 0s
    ##  47750K .......... .......... .......... .......... .......... 59% 81.5M 0s
    ##  47800K .......... .......... .......... .......... .......... 59% 44.0M 0s
    ##  47850K .......... .......... .......... .......... .......... 59% 26.6M 0s
    ##  47900K .......... .......... .......... .......... .......... 59% 55.7M 0s
    ##  47950K .......... .......... .......... .......... .......... 60% 71.5M 0s
    ##  48000K .......... .......... .......... .......... .......... 60% 75.6M 0s
    ##  48050K .......... .......... .......... .......... .......... 60%  131M 0s
    ##  48100K .......... .......... .......... .......... .......... 60% 31.2M 0s
    ##  48150K .......... .......... .......... .......... .......... 60% 84.8M 0s
    ##  48200K .......... .......... .......... .......... .......... 60%  138M 0s
    ##  48250K .......... .......... .......... .......... .......... 60%  151M 0s
    ##  48300K .......... .......... .......... .......... .......... 60% 30.6M 0s
    ##  48350K .......... .......... .......... .......... .......... 60% 69.6M 0s
    ##  48400K .......... .......... .......... .......... .......... 60% 72.6M 0s
    ##  48450K .......... .......... .......... .......... .......... 60%  104M 0s
    ##  48500K .......... .......... .......... .......... .......... 60%  123M 0s
    ##  48550K .......... .......... .......... .......... .......... 60% 25.8M 0s
    ##  48600K .......... .......... .......... .......... .......... 60% 96.8M 0s
    ##  48650K .......... .......... .......... .......... .......... 60%  137M 0s
    ##  48700K .......... .......... .......... .......... .......... 60%  121M 0s
    ##  48750K .......... .......... .......... .......... .......... 61% 43.1M 0s
    ##  48800K .......... .......... .......... .......... .......... 61% 62.3M 0s
    ##  48850K .......... .......... .......... .......... .......... 61% 59.2M 0s
    ##  48900K .......... .......... .......... .......... .......... 61% 89.6M 0s
    ##  48950K .......... .......... .......... .......... .......... 61%  141M 0s
    ##  49000K .......... .......... .......... .......... .......... 61% 31.1M 0s
    ##  49050K .......... .......... .......... .......... .......... 61%  116M 0s
    ##  49100K .......... .......... .......... .......... .......... 61%  107M 0s
    ##  49150K .......... .......... .......... .......... .......... 61%  147M 0s
    ##  49200K .......... .......... .......... .......... .......... 61% 47.8M 0s
    ##  49250K .......... .......... .......... .......... .......... 61%  126M 0s
    ##  49300K .......... .......... .......... .......... .......... 61% 40.1M 0s
    ##  49350K .......... .......... .......... .......... .......... 61%  115M 0s
    ##  49400K .......... .......... .......... .......... .......... 61% 66.0M 0s
    ##  49450K .......... .......... .......... .......... .......... 61% 63.8M 0s
    ##  49500K .......... .......... .......... .......... .......... 61% 63.1M 0s
    ##  49550K .......... .......... .......... .......... .......... 62%  128M 0s
    ##  49600K .......... .......... .......... .......... .......... 62% 92.7M 0s
    ##  49650K .......... .......... .......... .......... .......... 62% 60.0M 0s
    ##  49700K .......... .......... .......... .......... .......... 62% 99.5M 0s
    ##  49750K .......... .......... .......... .......... .......... 62% 45.3M 0s
    ##  49800K .......... .......... .......... .......... .......... 62% 81.4M 0s
    ##  49850K .......... .......... .......... .......... .......... 62%  179M 0s
    ##  49900K .......... .......... .......... .......... .......... 62% 35.2M 0s
    ##  49950K .......... .......... .......... .......... .......... 62%  102M 0s
    ##  50000K .......... .......... .......... .......... .......... 62%  105M 0s
    ##  50050K .......... .......... .......... .......... .......... 62% 49.9M 0s
    ##  50100K .......... .......... .......... .......... .......... 62%  107M 0s
    ##  50150K .......... .......... .......... .......... .......... 62%  128M 0s
    ##  50200K .......... .......... .......... .......... .......... 62% 78.0M 0s
    ##  50250K .......... .......... .......... .......... .......... 62%  117M 0s
    ##  50300K .......... .......... .......... .......... .......... 62% 43.5M 0s
    ##  50350K .......... .......... .......... .......... .......... 63%  118M 0s
    ##  50400K .......... .......... .......... .......... .......... 63% 92.2M 0s
    ##  50450K .......... .......... .......... .......... .......... 63%  111M 0s
    ##  50500K .......... .......... .......... .......... .......... 63% 42.0M 0s
    ##  50550K .......... .......... .......... .......... .......... 63%  119M 0s
    ##  50600K .......... .......... .......... .......... .......... 63%  132M 0s
    ##  50650K .......... .......... .......... .......... .......... 63% 56.3M 0s
    ##  50700K .......... .......... .......... .......... .......... 63% 99.5M 0s
    ##  50750K .......... .......... .......... .......... .......... 63% 51.7M 0s
    ##  50800K .......... .......... .......... .......... .......... 63% 74.9M 0s
    ##  50850K .......... .......... .......... .......... .......... 63%  105M 0s
    ##  50900K .......... .......... .......... .......... .......... 63%  110M 0s
    ##  50950K .......... .......... .......... .......... .......... 63% 54.8M 0s
    ##  51000K .......... .......... .......... .......... .......... 63%  102M 0s
    ##  51050K .......... .......... .......... .......... .......... 63% 69.1M 0s
    ##  51100K .......... .......... .......... .......... .......... 63% 98.9M 0s
    ##  51150K .......... .......... .......... .......... .......... 64%  134M 0s
    ##  51200K .......... .......... .......... .......... .......... 64% 59.3M 0s
    ##  51250K .......... .......... .......... .......... .......... 64% 60.2M 0s
    ##  51300K .......... .......... .......... .......... .......... 64% 98.5M 0s
    ##  51350K .......... .......... .......... .......... .......... 64%  129M 0s
    ##  51400K .......... .......... .......... .......... .......... 64% 85.6M 0s
    ##  51450K .......... .......... .......... .......... .......... 64%  122M 0s
    ##  51500K .......... .......... .......... .......... .......... 64% 56.7M 0s
    ##  51550K .......... .......... .......... .......... .......... 64% 92.4M 0s
    ##  51600K .......... .......... .......... .......... .......... 64%  140M 0s
    ##  51650K .......... .......... .......... .......... .......... 64% 54.1M 0s
    ##  51700K .......... .......... .......... .......... .......... 64%  114M 0s
    ##  51750K .......... .......... .......... .......... .......... 64%  107M 0s
    ##  51800K .......... .......... .......... .......... .......... 64%  121M 0s
    ##  51850K .......... .......... .......... .......... .......... 64% 50.8M 0s
    ##  51900K .......... .......... .......... .......... .......... 65%  110M 0s
    ##  51950K .......... .......... .......... .......... .......... 65% 67.6M 0s
    ##  52000K .......... .......... .......... .......... .......... 65% 98.9M 0s
    ##  52050K .......... .......... .......... .......... .......... 65%  143M 0s
    ##  52100K .......... .......... .......... .......... .......... 65% 48.8M 0s
    ##  52150K .......... .......... .......... .......... .......... 65%  142M 0s
    ##  52200K .......... .......... .......... .......... .......... 65% 85.9M 0s
    ##  52250K .......... .......... .......... .......... .......... 65%  115M 0s
    ##  52300K .......... .......... .......... .......... .......... 65% 78.7M 0s
    ##  52350K .......... .......... .......... .......... .......... 65%  127M 0s
    ##  52400K .......... .......... .......... .......... .......... 65% 50.6M 0s
    ##  52450K .......... .......... .......... .......... .......... 65%  110M 0s
    ##  52500K .......... .......... .......... .......... .......... 65%  121M 0s
    ##  52550K .......... .......... .......... .......... .......... 65%  110M 0s
    ##  52600K .......... .......... .......... .......... .......... 65% 59.2M 0s
    ##  52650K .......... .......... .......... .......... .......... 65%  100M 0s
    ##  52700K .......... .......... .......... .......... .......... 66% 69.3M 0s
    ##  52750K .......... .......... .......... .......... .......... 66%  114M 0s
    ##  52800K .......... .......... .......... .......... .......... 66%  112M 0s
    ##  52850K .......... .......... .......... .......... .......... 66% 75.7M 0s
    ##  52900K .......... .......... .......... .......... .......... 66%  100M 0s
    ##  52950K .......... .......... .......... .......... .......... 66%  147M 0s
    ##  53000K .......... .......... .......... .......... .......... 66% 77.8M 0s
    ##  53050K .......... .......... .......... .......... .......... 66% 57.4M 0s
    ##  53100K .......... .......... .......... .......... .......... 66%  106M 0s
    ##  53150K .......... .......... .......... .......... .......... 66%  172M 0s
    ##  53200K .......... .......... .......... .......... .......... 66% 69.1M 0s
    ##  53250K .......... .......... .......... .......... .......... 66%  150M 0s
    ##  53300K .......... .......... .......... .......... .......... 66% 56.2M 0s
    ##  53350K .......... .......... .......... .......... .......... 66%  139M 0s
    ##  53400K .......... .......... .......... .......... .......... 66%  113M 0s
    ##  53450K .......... .......... .......... .......... .......... 66%  168M 0s
    ##  53500K .......... .......... .......... .......... .......... 67% 79.7M 0s
    ##  53550K .......... .......... .......... .......... .......... 67%  101M 0s
    ##  53600K .......... .......... .......... .......... .......... 67% 58.6M 0s
    ##  53650K .......... .......... .......... .......... .......... 67%  117M 0s
    ##  53700K .......... .......... .......... .......... .......... 67%  119M 0s
    ##  53750K .......... .......... .......... .......... .......... 67% 62.0M 0s
    ##  53800K .......... .......... .......... .......... .......... 67% 98.7M 0s
    ##  53850K .......... .......... .......... .......... .......... 67% 95.6M 0s
    ##  53900K .......... .......... .......... .......... .......... 67% 64.2M 0s
    ##  53950K .......... .......... .......... .......... .......... 67%  116M 0s
    ##  54000K .......... .......... .......... .......... .......... 67%  124M 0s
    ##  54050K .......... .......... .......... .......... .......... 67%  117M 0s
    ##  54100K .......... .......... .......... .......... .......... 67%  100M 0s
    ##  54150K .......... .......... .......... .......... .......... 67%  140M 0s
    ##  54200K .......... .......... .......... .......... .......... 67%  109M 0s
    ##  54250K .......... .......... .......... .......... .......... 67%  107M 0s
    ##  54300K .......... .......... .......... .......... .......... 68% 83.2M 0s
    ##  54350K .......... .......... .......... .......... .......... 68% 51.4M 0s
    ##  54400K .......... .......... .......... .......... .......... 68% 92.5M 0s
    ##  54450K .......... .......... .......... .......... .......... 68%  166M 0s
    ##  54500K .......... .......... .......... .......... .......... 68% 81.0M 0s
    ##  54550K .......... .......... .......... .......... .......... 68%  113M 0s
    ##  54600K .......... .......... .......... .......... .......... 68%  163M 0s
    ##  54650K .......... .......... .......... .......... .......... 68% 83.9M 0s
    ##  54700K .......... .......... .......... .......... .......... 68%  120M 0s
    ##  54750K .......... .......... .......... .......... .......... 68%  110M 0s
    ##  54800K .......... .......... .......... .......... .......... 68% 72.9M 0s
    ##  54850K .......... .......... .......... .......... .......... 68%  117M 0s
    ##  54900K .......... .......... .......... .......... .......... 68%  121M 0s
    ##  54950K .......... .......... .......... .......... .......... 68% 66.8M 0s
    ##  55000K .......... .......... .......... .......... .......... 68% 91.3M 0s
    ##  55050K .......... .......... .......... .......... .......... 68%  121M 0s
    ##  55100K .......... .......... .......... .......... .......... 69% 89.8M 0s
    ##  55150K .......... .......... .......... .......... .......... 69%  128M 0s
    ##  55200K .......... .......... .......... .......... .......... 69%  114M 0s
    ##  55250K .......... .......... .......... .......... .......... 69%  112M 0s
    ##  55300K .......... .......... .......... .......... .......... 69% 92.5M 0s
    ##  55350K .......... .......... .......... .......... .......... 69%  140M 0s
    ##  55400K .......... .......... .......... .......... .......... 69% 79.0M 0s
    ##  55450K .......... .......... .......... .......... .......... 69% 95.4M 0s
    ##  55500K .......... .......... .......... .......... .......... 69% 91.3M 0s
    ##  55550K .......... .......... .......... .......... .......... 69% 97.4M 0s
    ##  55600K .......... .......... .......... .......... .......... 69%  120M 0s
    ##  55650K .......... .......... .......... .......... .......... 69% 79.2M 0s
    ##  55700K .......... .......... .......... .......... .......... 69%  142M 0s
    ##  55750K .......... .......... .......... .......... .......... 69% 18.6M 0s
    ##  55800K .......... .......... .......... .......... .......... 69%  111M 0s
    ##  55850K .......... .......... .......... .......... .......... 69%  111M 0s
    ##  55900K .......... .......... .......... .......... .......... 70%  106M 0s
    ##  55950K .......... .......... .......... .......... .......... 70%  107M 0s
    ##  56000K .......... .......... .......... .......... .......... 70%  141M 0s
    ##  56050K .......... .......... .......... .......... .......... 70%  119M 0s
    ##  56100K .......... .......... .......... .......... .......... 70%  136M 0s
    ##  56150K .......... .......... .......... .......... .......... 70%  128M 0s
    ##  56200K .......... .......... .......... .......... .......... 70% 35.4M 0s
    ##  56250K .......... .......... .......... .......... .......... 70%  103M 0s
    ##  56300K .......... .......... .......... .......... .......... 70% 25.4M 0s
    ##  56350K .......... .......... .......... .......... .......... 70%  109M 0s
    ##  56400K .......... .......... .......... .......... .......... 70% 41.2M 0s
    ##  56450K .......... .......... .......... .......... .......... 70%  112M 0s
    ##  56500K .......... .......... .......... .......... .......... 70% 48.4M 0s
    ##  56550K .......... .......... .......... .......... .......... 70%  123M 0s
    ##  56600K .......... .......... .......... .......... .......... 70%  144M 0s
    ##  56650K .......... .......... .......... .......... .......... 70% 48.5M 0s
    ##  56700K .......... .......... .......... .......... .......... 71% 58.5M 0s
    ##  56750K .......... .......... .......... .......... .......... 71%  112M 0s
    ##  56800K .......... .......... .......... .......... .......... 71% 98.4M 0s
    ##  56850K .......... .......... .......... .......... .......... 71% 25.7M 0s
    ##  56900K .......... .......... .......... .......... .......... 71% 71.8M 0s
    ##  56950K .......... .......... .......... .......... .......... 71%  100M 0s
    ##  57000K .......... .......... .......... .......... .......... 71% 81.0M 0s
    ##  57050K .......... .......... .......... .......... .......... 71%  106M 0s
    ##  57100K .......... .......... .......... .......... .......... 71% 72.9M 0s
    ##  57150K .......... .......... .......... .......... .......... 71%  107M 0s
    ##  57200K .......... .......... .......... .......... .......... 71%  120M 0s
    ##  57250K .......... .......... .......... .......... .......... 71%  115M 0s
    ##  57300K .......... .......... .......... .......... .......... 71% 95.2M 0s
    ##  57350K .......... .......... .......... .......... .......... 71%  134M 0s
    ##  57400K .......... .......... .......... .......... .......... 71%  112M 0s
    ##  57450K .......... .......... .......... .......... .......... 71%  139M 0s
    ##  57500K .......... .......... .......... .......... .......... 72%  124M 0s
    ##  57550K .......... .......... .......... .......... .......... 72%  133M 0s
    ##  57600K .......... .......... .......... .......... .......... 72%  126M 0s
    ##  57650K .......... .......... .......... .......... .......... 72%  111M 0s
    ##  57700K .......... .......... .......... .......... .......... 72%  110M 0s
    ##  57750K .......... .......... .......... .......... .......... 72% 58.2M 0s
    ##  57800K .......... .......... .......... .......... .......... 72% 82.7M 0s
    ##  57850K .......... .......... .......... .......... .......... 72%  117M 0s
    ##  57900K .......... .......... .......... .......... .......... 72%  121M 0s
    ##  57950K .......... .......... .......... .......... .......... 72% 81.0M 0s
    ##  58000K .......... .......... .......... .......... .......... 72%  108M 0s
    ##  58050K .......... .......... .......... .......... .......... 72% 85.1M 0s
    ##  58100K .......... .......... .......... .......... .......... 72%  119M 0s
    ##  58150K .......... .......... .......... .......... .......... 72%  134M 0s
    ##  58200K .......... .......... .......... .......... .......... 72% 84.5M 0s
    ##  58250K .......... .......... .......... .......... .......... 72%  103M 0s
    ##  58300K .......... .......... .......... .......... .......... 73%  144M 0s
    ##  58350K .......... .......... .......... .......... .......... 73%  136M 0s
    ##  58400K .......... .......... .......... .......... .......... 73% 95.5M 0s
    ##  58450K .......... .......... .......... .......... .......... 73%  141M 0s
    ##  58500K .......... .......... .......... .......... .......... 73% 90.3M 0s
    ##  58550K .......... .......... .......... .......... .......... 73%  114M 0s
    ##  58600K .......... .......... .......... .......... .......... 73% 88.3M 0s
    ##  58650K .......... .......... .......... .......... .......... 73%  123M 0s
    ##  58700K .......... .......... .......... .......... .......... 73%  125M 0s
    ##  58750K .......... .......... .......... .......... .......... 73%  123M 0s
    ##  58800K .......... .......... .......... .......... .......... 73% 95.2M 0s
    ##  58850K .......... .......... .......... .......... .......... 73%  110M 0s
    ##  58900K .......... .......... .......... .......... .......... 73%  145M 0s
    ##  58950K .......... .......... .......... .......... .......... 73%  149M 0s
    ##  59000K .......... .......... .......... .......... .......... 73%  105M 0s
    ##  59050K .......... .......... .......... .......... .......... 73%  171M 0s
    ##  59100K .......... .......... .......... .......... .......... 74% 70.3M 0s
    ##  59150K .......... .......... .......... .......... .......... 74%  119M 0s
    ##  59200K .......... .......... .......... .......... .......... 74% 67.4M 0s
    ##  59250K .......... .......... .......... .......... .......... 74%  115M 0s
    ##  59300K .......... .......... .......... .......... .......... 74%  109M 0s
    ##  59350K .......... .......... .......... .......... .......... 74%  142M 0s
    ##  59400K .......... .......... .......... .......... .......... 74%  103M 0s
    ##  59450K .......... .......... .......... .......... .......... 74%  126M 0s
    ##  59500K .......... .......... .......... .......... .......... 74%  125M 0s
    ##  59550K .......... .......... .......... .......... .......... 74%  162M 0s
    ##  59600K .......... .......... .......... .......... .......... 74%  137M 0s
    ##  59650K .......... .......... .......... .......... .......... 74% 89.7M 0s
    ##  59700K .......... .......... .......... .......... .......... 74% 86.3M 0s
    ##  59750K .......... .......... .......... .......... .......... 74%  122M 0s
    ##  59800K .......... .......... .......... .......... .......... 74%  101M 0s
    ##  59850K .......... .......... .......... .......... .......... 74%  133M 0s
    ##  59900K .......... .......... .......... .......... .......... 75%  115M 0s
    ##  59950K .......... .......... .......... .......... .......... 75% 78.4M 0s
    ##  60000K .......... .......... .......... .......... .......... 75% 78.7M 0s
    ##  60050K .......... .......... .......... .......... .......... 75% 84.1M 0s
    ##  60100K .......... .......... .......... .......... .......... 75%  109M 0s
    ##  60150K .......... .......... .......... .......... .......... 75%  141M 0s
    ##  60200K .......... .......... .......... .......... .......... 75%  127M 0s
    ##  60250K .......... .......... .......... .......... .......... 75%  156M 0s
    ##  60300K .......... .......... .......... .......... .......... 75%  113M 0s
    ##  60350K .......... .......... .......... .......... .......... 75%  129M 0s
    ##  60400K .......... .......... .......... .......... .......... 75%  116M 0s
    ##  60450K .......... .......... .......... .......... .......... 75%  131M 0s
    ##  60500K .......... .......... .......... .......... .......... 75%  113M 0s
    ##  60550K .......... .......... .......... .......... .......... 75%  101M 0s
    ##  60600K .......... .......... .......... .......... .......... 75%  109M 0s
    ##  60650K .......... .......... .......... .......... .......... 75%  116M 0s
    ##  60700K .......... .......... .......... .......... .......... 76%  126M 0s
    ##  60750K .......... .......... .......... .......... .......... 76% 92.2M 0s
    ##  60800K .......... .......... .......... .......... .......... 76%  113M 0s
    ##  60850K .......... .......... .......... .......... .......... 76%  102M 0s
    ##  60900K .......... .......... .......... .......... .......... 76% 91.8M 0s
    ##  60950K .......... .......... .......... .......... .......... 76% 98.7M 0s
    ##  61000K .......... .......... .......... .......... .......... 76%  144M 0s
    ##  61050K .......... .......... .......... .......... .......... 76%  101M 0s
    ##  61100K .......... .......... .......... .......... .......... 76%  101M 0s
    ##  61150K .......... .......... .......... .......... .......... 76%  122M 0s
    ##  61200K .......... .......... .......... .......... .......... 76%  113M 0s
    ##  61250K .......... .......... .......... .......... .......... 76%  140M 0s
    ##  61300K .......... .......... .......... .......... .......... 76%  147M 0s
    ##  61350K .......... .......... .......... .......... .......... 76%  102M 0s
    ##  61400K .......... .......... .......... .......... .......... 76%  105M 0s
    ##  61450K .......... .......... .......... .......... .......... 76%  126M 0s
    ##  61500K .......... .......... .......... .......... .......... 77%  107M 0s
    ##  61550K .......... .......... .......... .......... .......... 77%  121M 0s
    ##  61600K .......... .......... .......... .......... .......... 77%  143M 0s
    ##  61650K .......... .......... .......... .......... .......... 77%  117M 0s
    ##  61700K .......... .......... .......... .......... .......... 77%  117M 0s
    ##  61750K .......... .......... .......... .......... .......... 77%  119M 0s
    ##  61800K .......... .......... .......... .......... .......... 77%  127M 0s
    ##  61850K .......... .......... .......... .......... .......... 77%  141M 0s
    ##  61900K .......... .......... .......... .......... .......... 77%  136M 0s
    ##  61950K .......... .......... .......... .......... .......... 77%  104M 0s
    ##  62000K .......... .......... .......... .......... .......... 77%  120M 0s
    ##  62050K .......... .......... .......... .......... .......... 77%  134M 0s
    ##  62100K .......... .......... .......... .......... .......... 77%  110M 0s
    ##  62150K .......... .......... .......... .......... .......... 77%  119M 0s
    ##  62200K .......... .......... .......... .......... .......... 77%  129M 0s
    ##  62250K .......... .......... .......... .......... .......... 77%  106M 0s
    ##  62300K .......... .......... .......... .......... .......... 78% 94.5M 0s
    ##  62350K .......... .......... .......... .......... .......... 78%  119M 0s
    ##  62400K .......... .......... .......... .......... .......... 78%  101M 0s
    ##  62450K .......... .......... .......... .......... .......... 78%  113M 0s
    ##  62500K .......... .......... .......... .......... .......... 78%  151M 0s
    ##  62550K .......... .......... .......... .......... .......... 78%  107M 0s
    ##  62600K .......... .......... .......... .......... .......... 78%  123M 0s
    ##  62650K .......... .......... .......... .......... .......... 78%  120M 0s
    ##  62700K .......... .......... .......... .......... .......... 78%  125M 0s
    ##  62750K .......... .......... .......... .......... .......... 78%  147M 0s
    ##  62800K .......... .......... .......... .......... .......... 78%  100M 0s
    ##  62850K .......... .......... .......... .......... .......... 78%  127M 0s
    ##  62900K .......... .......... .......... .......... .......... 78%  103M 0s
    ##  62950K .......... .......... .......... .......... .......... 78%  117M 0s
    ##  63000K .......... .......... .......... .......... .......... 78%  126M 0s
    ##  63050K .......... .......... .......... .......... .......... 78%  145M 0s
    ##  63100K .......... .......... .......... .......... .......... 79%  144M 0s
    ##  63150K .......... .......... .......... .......... .......... 79%  115M 0s
    ##  63200K .......... .......... .......... .......... .......... 79%  110M 0s
    ##  63250K .......... .......... .......... .......... .......... 79%  101M 0s
    ##  63300K .......... .......... .......... .......... .......... 79% 89.9M 0s
    ##  63350K .......... .......... .......... .......... .......... 79% 98.9M 0s
    ##  63400K .......... .......... .......... .......... .......... 79%  115M 0s
    ##  63450K .......... .......... .......... .......... .......... 79%  124M 0s
    ##  63500K .......... .......... .......... .......... .......... 79% 89.8M 0s
    ##  63550K .......... .......... .......... .......... .......... 79%  116M 0s
    ##  63600K .......... .......... .......... .......... .......... 79%  122M 0s
    ##  63650K .......... .......... .......... .......... .......... 79%  155M 0s
    ##  63700K .......... .......... .......... .......... .......... 79% 97.0M 0s
    ##  63750K .......... .......... .......... .......... .......... 79%  114M 0s
    ##  63800K .......... .......... .......... .......... .......... 79% 99.5M 0s
    ##  63850K .......... .......... .......... .......... .......... 79%  111M 0s
    ##  63900K .......... .......... .......... .......... .......... 80%  123M 0s
    ##  63950K .......... .......... .......... .......... .......... 80%  141M 0s
    ##  64000K .......... .......... .......... .......... .......... 80%  101M 0s
    ##  64050K .......... .......... .......... .......... .......... 80%  141M 0s
    ##  64100K .......... .......... .......... .......... .......... 80%  122M 0s
    ##  64150K .......... .......... .......... .......... .......... 80%  116M 0s
    ##  64200K .......... .......... .......... .......... .......... 80%  112M 0s
    ##  64250K .......... .......... .......... .......... .......... 80%  122M 0s
    ##  64300K .......... .......... .......... .......... .......... 80%  103M 0s
    ##  64350K .......... .......... .......... .......... .......... 80%  117M 0s
    ##  64400K .......... .......... .......... .......... .......... 80%  100M 0s
    ##  64450K .......... .......... .......... .......... .......... 80%  141M 0s
    ##  64500K .......... .......... .......... .......... .......... 80%  131M 0s
    ##  64550K .......... .......... .......... .......... .......... 80%  146M 0s
    ##  64600K .......... .......... .......... .......... .......... 80%  137M 0s
    ##  64650K .......... .......... .......... .......... .......... 80%  125M 0s
    ##  64700K .......... .......... .......... .......... .......... 81%  100M 0s
    ##  64750K .......... .......... .......... .......... .......... 81%  124M 0s
    ##  64800K .......... .......... .......... .......... .......... 81%  108M 0s
    ##  64850K .......... .......... .......... .......... .......... 81%  142M 0s
    ##  64900K .......... .......... .......... .......... .......... 81% 88.2M 0s
    ##  64950K .......... .......... .......... .......... .......... 81%  171M 0s
    ##  65000K .......... .......... .......... .......... .......... 81%  116M 0s
    ##  65050K .......... .......... .......... .......... .......... 81% 93.4M 0s
    ##  65100K .......... .......... .......... .......... .......... 81%  102M 0s
    ##  65150K .......... .......... .......... .......... .......... 81%  120M 0s
    ##  65200K .......... .......... .......... .......... .......... 81% 26.1M 0s
    ##  65250K .......... .......... .......... .......... .......... 81%  102M 0s
    ##  65300K .......... .......... .......... .......... .......... 81%  120M 0s
    ##  65350K .......... .......... .......... .......... .......... 81%  142M 0s
    ##  65400K .......... .......... .......... .......... .......... 81% 65.5M 0s
    ##  65450K .......... .......... .......... .......... .......... 81% 76.9M 0s
    ##  65500K .......... .......... .......... .......... .......... 82%  112M 0s
    ##  65550K .......... .......... .......... .......... .......... 82%  102M 0s
    ##  65600K .......... .......... .......... .......... .......... 82%  120M 0s
    ##  65650K .......... .......... .......... .......... .......... 82%  170M 0s
    ##  65700K .......... .......... .......... .......... .......... 82%  126M 0s
    ##  65750K .......... .......... .......... .......... .......... 82% 61.5M 0s
    ##  65800K .......... .......... .......... .......... .......... 82%  122M 0s
    ##  65850K .......... .......... .......... .......... .......... 82%  108M 0s
    ##  65900K .......... .......... .......... .......... .......... 82% 62.1M 0s
    ##  65950K .......... .......... .......... .......... .......... 82% 95.1M 0s
    ##  66000K .......... .......... .......... .......... .......... 82% 92.4M 0s
    ##  66050K .......... .......... .......... .......... .......... 82%  103M 0s
    ##  66100K .......... .......... .......... .......... .......... 82%  143M 0s
    ##  66150K .......... .......... .......... .......... .......... 82%  158M 0s
    ##  66200K .......... .......... .......... .......... .......... 82%  119M 0s
    ##  66250K .......... .......... .......... .......... .......... 82% 11.6M 0s
    ##  66300K .......... .......... .......... .......... .......... 83%  127M 0s
    ##  66350K .......... .......... .......... .......... .......... 83%  184M 0s
    ##  66400K .......... .......... .......... .......... .......... 83%  148M 0s
    ##  66450K .......... .......... .......... .......... .......... 83%  128M 0s
    ##  66500K .......... .......... .......... .......... .......... 83%  119M 0s
    ##  66550K .......... .......... .......... .......... .......... 83%  114M 0s
    ##  66600K .......... .......... .......... .......... .......... 83%  117M 0s
    ##  66650K .......... .......... .......... .......... .......... 83%  115M 0s
    ##  66700K .......... .......... .......... .......... .......... 83% 30.3M 0s
    ##  66750K .......... .......... .......... .......... .......... 83%  140M 0s
    ##  66800K .......... .......... .......... .......... .......... 83% 23.5M 0s
    ##  66850K .......... .......... .......... .......... .......... 83%  115M 0s
    ##  66900K .......... .......... .......... .......... .......... 83% 30.9M 0s
    ##  66950K .......... .......... .......... .......... .......... 83%  130M 0s
    ##  67000K .......... .......... .......... .......... .......... 83% 26.3M 0s
    ##  67050K .......... .......... .......... .......... .......... 83%  176M 0s
    ##  67100K .......... .......... .......... .......... .......... 84%  123M 0s
    ##  67150K .......... .......... .......... .......... .......... 84%  127M 0s
    ##  67200K .......... .......... .......... .......... .......... 84%  184M 0s
    ##  67250K .......... .......... .......... .......... .......... 84%  129M 0s
    ##  67300K .......... .......... .......... .......... .......... 84%  145M 0s
    ##  67350K .......... .......... .......... .......... .......... 84%  149M 0s
    ##  67400K .......... .......... .......... .......... .......... 84%  131M 0s
    ##  67450K .......... .......... .......... .......... .......... 84%  141M 0s
    ##  67500K .......... .......... .......... .......... .......... 84%  126M 0s
    ##  67550K .......... .......... .......... .......... .......... 84% 60.6M 0s
    ##  67600K .......... .......... .......... .......... .......... 84%  115M 0s
    ##  67650K .......... .......... .......... .......... .......... 84% 20.1M 0s
    ##  67700K .......... .......... .......... .......... .......... 84% 81.7M 0s
    ##  67750K .......... .......... .......... .......... .......... 84%  135M 0s
    ##  67800K .......... .......... .......... .......... .......... 84% 64.5M 0s
    ##  67850K .......... .......... .......... .......... .......... 84% 38.9M 0s
    ##  67900K .......... .......... .......... .......... .......... 85%  112M 0s
    ##  67950K .......... .......... .......... .......... .......... 85%  113M 0s
    ##  68000K .......... .......... .......... .......... .......... 85% 38.8M 0s
    ##  68050K .......... .......... .......... .......... .......... 85%  108M 0s
    ##  68100K .......... .......... .......... .......... .......... 85% 61.7M 0s
    ##  68150K .......... .......... .......... .......... .......... 85% 83.0M 0s
    ##  68200K .......... .......... .......... .......... .......... 85%  103M 0s
    ##  68250K .......... .......... .......... .......... .......... 85%  117M 0s
    ##  68300K .......... .......... .......... .......... .......... 85% 30.9M 0s
    ##  68350K .......... .......... .......... .......... .......... 85%  154M 0s
    ##  68400K .......... .......... .......... .......... .......... 85% 82.4M 0s
    ##  68450K .......... .......... .......... .......... .......... 85% 24.8M 0s
    ##  68500K .......... .......... .......... .......... .......... 85%  114M 0s
    ##  68550K .......... .......... .......... .......... .......... 85%  135M 0s
    ##  68600K .......... .......... .......... .......... .......... 85% 61.9M 0s
    ##  68650K .......... .......... .......... .......... .......... 85%  113M 0s
    ##  68700K .......... .......... .......... .......... .......... 86% 31.5M 0s
    ##  68750K .......... .......... .......... .......... .......... 86% 83.6M 0s
    ##  68800K .......... .......... .......... .......... .......... 86%  102M 0s
    ##  68850K .......... .......... .......... .......... .......... 86%  167M 0s
    ##  68900K .......... .......... .......... .......... .......... 86% 47.6M 0s
    ##  68950K .......... .......... .......... .......... .......... 86%  126M 0s
    ##  69000K .......... .......... .......... .......... .......... 86% 82.3M 0s
    ##  69050K .......... .......... .......... .......... .......... 86% 59.1M 0s
    ##  69100K .......... .......... .......... .......... .......... 86%  101M 0s
    ##  69150K .......... .......... .......... .......... .......... 86% 34.8M 0s
    ##  69200K .......... .......... .......... .......... .......... 86% 85.6M 0s
    ##  69250K .......... .......... .......... .......... .......... 86%  106M 0s
    ##  69300K .......... .......... .......... .......... .......... 86%  114M 0s
    ##  69350K .......... .......... .......... .......... .......... 86% 64.3M 0s
    ##  69400K .......... .......... .......... .......... .......... 86%  102M 0s
    ##  69450K .......... .......... .......... .......... .......... 86%  133M 0s
    ##  69500K .......... .......... .......... .......... .......... 87% 60.8M 0s
    ##  69550K .......... .......... .......... .......... .......... 87%  109M 0s
    ##  69600K .......... .......... .......... .......... .......... 87% 44.5M 0s
    ##  69650K .......... .......... .......... .......... .......... 87% 99.7M 0s
    ##  69700K .......... .......... .......... .......... .......... 87% 47.5M 0s
    ##  69750K .......... .......... .......... .......... .......... 87%  101M 0s
    ##  69800K .......... .......... .......... .......... .......... 87%  108M 0s
    ##  69850K .......... .......... .......... .......... .......... 87% 38.5M 0s
    ##  69900K .......... .......... .......... .......... .......... 87%  112M 0s
    ##  69950K .......... .......... .......... .......... .......... 87% 92.6M 0s
    ##  70000K .......... .......... .......... .......... .......... 87%  121M 0s
    ##  70050K .......... .......... .......... .......... .......... 87% 43.5M 0s
    ##  70100K .......... .......... .......... .......... .......... 87% 96.1M 0s
    ##  70150K .......... .......... .......... .......... .......... 87%  133M 0s
    ##  70200K .......... .......... .......... .......... .......... 87% 63.9M 0s
    ##  70250K .......... .......... .......... .......... .......... 87% 72.5M 0s
    ##  70300K .......... .......... .......... .......... .......... 88% 97.9M 0s
    ##  70350K .......... .......... .......... .......... .......... 88% 50.7M 0s
    ##  70400K .......... .......... .......... .......... .......... 88% 93.5M 0s
    ##  70450K .......... .......... .......... .......... .......... 88%  139M 0s
    ##  70500K .......... .......... .......... .......... .......... 88% 59.1M 0s
    ##  70550K .......... .......... .......... .......... .......... 88% 91.1M 0s
    ##  70600K .......... .......... .......... .......... .......... 88% 63.8M 0s
    ##  70650K .......... .......... .......... .......... .......... 88%  111M 0s
    ##  70700K .......... .......... .......... .......... .......... 88% 89.5M 0s
    ##  70750K .......... .......... .......... .......... .......... 88% 92.5M 0s
    ##  70800K .......... .......... .......... .......... .......... 88% 47.3M 0s
    ##  70850K .......... .......... .......... .......... .......... 88%  109M 0s
    ##  70900K .......... .......... .......... .......... .......... 88%  115M 0s
    ##  70950K .......... .......... .......... .......... .......... 88% 64.9M 0s
    ##  71000K .......... .......... .......... .......... .......... 88%  114M 0s
    ##  71050K .......... .......... .......... .......... .......... 88% 61.9M 0s
    ##  71100K .......... .......... .......... .......... .......... 89%  125M 0s
    ##  71150K .......... .......... .......... .......... .......... 89% 74.1M 0s
    ##  71200K .......... .......... .......... .......... .......... 89%  113M 0s
    ##  71250K .......... .......... .......... .......... .......... 89% 60.2M 0s
    ##  71300K .......... .......... .......... .......... .......... 89% 91.8M 0s
    ##  71350K .......... .......... .......... .......... .......... 89%  133M 0s
    ##  71400K .......... .......... .......... .......... .......... 89% 60.2M 0s
    ##  71450K .......... .......... .......... .......... .......... 89%  124M 0s
    ##  71500K .......... .......... .......... .......... .......... 89% 90.1M 0s
    ##  71550K .......... .......... .......... .......... .......... 89%  108M 0s
    ##  71600K .......... .......... .......... .......... .......... 89% 50.7M 0s
    ##  71650K .......... .......... .......... .......... .......... 89%  114M 0s
    ##  71700K .......... .......... .......... .......... .......... 89% 90.6M 0s
    ##  71750K .......... .......... .......... .......... .......... 89%  136M 0s
    ##  71800K .......... .......... .......... .......... .......... 89%  100M 0s
    ##  71850K .......... .......... .......... .......... .......... 89%  101M 0s
    ##  71900K .......... .......... .......... .......... .......... 90% 85.7M 0s
    ##  71950K .......... .......... .......... .......... .......... 90% 78.0M 0s
    ##  72000K .......... .......... .......... .......... .......... 90% 94.6M 0s
    ##  72050K .......... .......... .......... .......... .......... 90%  136M 0s
    ##  72100K .......... .......... .......... .......... .......... 90% 65.5M 0s
    ##  72150K .......... .......... .......... .......... .......... 90% 87.3M 0s
    ##  72200K .......... .......... .......... .......... .......... 90% 77.1M 0s
    ##  72250K .......... .......... .......... .......... .......... 90% 51.7M 0s
    ##  72300K .......... .......... .......... .......... .......... 90%  131M 0s
    ##  72350K .......... .......... .......... .......... .......... 90% 41.1M 0s
    ##  72400K .......... .......... .......... .......... .......... 90%  109M 0s
    ##  72450K .......... .......... .......... .......... .......... 90%  135M 0s
    ##  72500K .......... .......... .......... .......... .......... 90%  126M 0s
    ##  72550K .......... .......... .......... .......... .......... 90%  120M 0s
    ##  72600K .......... .......... .......... .......... .......... 90%  119M 0s
    ##  72650K .......... .......... .......... .......... .......... 90% 55.5M 0s
    ##  72700K .......... .......... .......... .......... .......... 91% 36.0M 0s
    ##  72750K .......... .......... .......... .......... .......... 91% 58.1M 0s
    ##  72800K .......... .......... .......... .......... .......... 91% 58.3M 0s
    ##  72850K .......... .......... .......... .......... .......... 91% 30.0M 0s
    ##  72900K .......... .......... .......... .......... .......... 91% 21.1M 0s
    ##  72950K .......... .......... .......... .......... .......... 91%  119M 0s
    ##  73000K .......... .......... .......... .......... .......... 91% 41.9M 0s
    ##  73050K .......... .......... .......... .......... .......... 91%  109M 0s
    ##  73100K .......... .......... .......... .......... .......... 91%  108M 0s
    ##  73150K .......... .......... .......... .......... .......... 91% 55.4M 0s
    ##  73200K .......... .......... .......... .......... .......... 91%  112M 0s
    ##  73250K .......... .......... .......... .......... .......... 91% 61.9M 0s
    ##  73300K .......... .......... .......... .......... .......... 91% 47.4M 0s
    ##  73350K .......... .......... .......... .......... .......... 91% 89.3M 0s
    ##  73400K .......... .......... .......... .......... .......... 91% 69.5M 0s
    ##  73450K .......... .......... .......... .......... .......... 91% 33.3M 0s
    ##  73500K .......... .......... .......... .......... .......... 92% 83.9M 0s
    ##  73550K .......... .......... .......... .......... .......... 92% 81.6M 0s
    ##  73600K .......... .......... .......... .......... .......... 92% 51.0M 0s
    ##  73650K .......... .......... .......... .......... .......... 92% 95.5M 0s
    ##  73700K .......... .......... .......... .......... .......... 92%  112M 0s
    ##  73750K .......... .......... .......... .......... .......... 92% 52.9M 0s
    ##  73800K .......... .......... .......... .......... .......... 92% 59.9M 0s
    ##  73850K .......... .......... .......... .......... .......... 92%  127M 0s
    ##  73900K .......... .......... .......... .......... .......... 92% 35.2M 0s
    ##  73950K .......... .......... .......... .......... .......... 92% 90.9M 0s
    ##  74000K .......... .......... .......... .......... .......... 92% 64.7M 0s
    ##  74050K .......... .......... .......... .......... .......... 92% 61.6M 0s
    ##  74100K .......... .......... .......... .......... .......... 92% 62.0M 0s
    ##  74150K .......... .......... .......... .......... .......... 92% 65.3M 0s
    ##  74200K .......... .......... .......... .......... .......... 92% 91.7M 0s
    ##  74250K .......... .......... .......... .......... .......... 92% 97.1M 0s
    ##  74300K .......... .......... .......... .......... .......... 93% 50.5M 0s
    ##  74350K .......... .......... .......... .......... .......... 93% 98.4M 0s
    ##  74400K .......... .......... .......... .......... .......... 93% 80.3M 0s
    ##  74450K .......... .......... .......... .......... .......... 93% 38.3M 0s
    ##  74500K .......... .......... .......... .......... .......... 93% 80.5M 0s
    ##  74550K .......... .......... .......... .......... .......... 93% 98.4M 0s
    ##  74600K .......... .......... .......... .......... .......... 93% 54.0M 0s
    ##  74650K .......... .......... .......... .......... .......... 93% 44.7M 0s
    ##  74700K .......... .......... .......... .......... .......... 93% 95.1M 0s
    ##  74750K .......... .......... .......... .......... .......... 93% 87.9M 0s
    ##  74800K .......... .......... .......... .......... .......... 93% 83.3M 0s
    ##  74850K .......... .......... .......... .......... .......... 93%  106M 0s
    ##  74900K .......... .......... .......... .......... .......... 93% 75.9M 0s
    ##  74950K .......... .......... .......... .......... .......... 93% 37.8M 0s
    ##  75000K .......... .......... .......... .......... .......... 93% 60.8M 0s
    ##  75050K .......... .......... .......... .......... .......... 93%  123M 0s
    ##  75100K .......... .......... .......... .......... .......... 94% 88.0M 0s
    ##  75150K .......... .......... .......... .......... .......... 94% 58.6M 0s
    ##  75200K .......... .......... .......... .......... .......... 94% 98.6M 0s
    ##  75250K .......... .......... .......... .......... .......... 94%  116M 0s
    ##  75300K .......... .......... .......... .......... .......... 94% 55.3M 0s
    ##  75350K .......... .......... .......... .......... .......... 94% 66.4M 0s
    ##  75400K .......... .......... .......... .......... .......... 94%  102M 0s
    ##  75450K .......... .......... .......... .......... .......... 94% 47.8M 0s
    ##  75500K .......... .......... .......... .......... .......... 94% 63.9M 0s
    ##  75550K .......... .......... .......... .......... .......... 94% 96.4M 0s
    ##  75600K .......... .......... .......... .......... .......... 94%  105M 0s
    ##  75650K .......... .......... .......... .......... .......... 94% 85.3M 0s
    ##  75700K .......... .......... .......... .......... .......... 94% 95.3M 0s
    ##  75750K .......... .......... .......... .......... .......... 94%  154M 0s
    ##  75800K .......... .......... .......... .......... .......... 94% 44.1M 0s
    ##  75850K .......... .......... .......... .......... .......... 94% 73.1M 0s
    ##  75900K .......... .......... .......... .......... .......... 95%  117M 0s
    ##  75950K .......... .......... .......... .......... .......... 95% 57.4M 0s
    ##  76000K .......... .......... .......... .......... .......... 95% 58.6M 0s
    ##  76050K .......... .......... .......... .......... .......... 95%  145M 0s
    ##  76100K .......... .......... .......... .......... .......... 95%  132M 0s
    ##  76150K .......... .......... .......... .......... .......... 95% 45.2M 0s
    ##  76200K .......... .......... .......... .......... .......... 95% 87.1M 0s
    ##  76250K .......... .......... .......... .......... .......... 95%  111M 0s
    ##  76300K .......... .......... .......... .......... .......... 95%  118M 0s
    ##  76350K .......... .......... .......... .......... .......... 95%  126M 0s
    ##  76400K .......... .......... .......... .......... .......... 95% 90.5M 0s
    ##  76450K .......... .......... .......... .......... .......... 95%  103M 0s
    ##  76500K .......... .......... .......... .......... .......... 95% 58.8M 0s
    ##  76550K .......... .......... .......... .......... .......... 95%  104M 0s
    ##  76600K .......... .......... .......... .......... .......... 95% 96.1M 0s
    ##  76650K .......... .......... .......... .......... .......... 95% 43.5M 0s
    ##  76700K .......... .......... .......... .......... .......... 96% 86.6M 0s
    ##  76750K .......... .......... .......... .......... .......... 96%  136M 0s
    ##  76800K .......... .......... .......... .......... .......... 96%  125M 0s
    ##  76850K .......... .......... .......... .......... .......... 96% 76.4M 0s
    ##  76900K .......... .......... .......... .......... .......... 96% 63.1M 0s
    ##  76950K .......... .......... .......... .......... .......... 96%  110M 0s
    ##  77000K .......... .......... .......... .......... .......... 96%  127M 0s
    ##  77050K .......... .......... .......... .......... .......... 96% 78.3M 0s
    ##  77100K .......... .......... .......... .......... .......... 96%  111M 0s
    ##  77150K .......... .......... .......... .......... .......... 96% 51.8M 0s
    ##  77200K .......... .......... .......... .......... .......... 96% 88.0M 0s
    ##  77250K .......... .......... .......... .......... .......... 96% 89.1M 0s
    ##  77300K .......... .......... .......... .......... .......... 96%  153M 0s
    ##  77350K .......... .......... .......... .......... .......... 96% 57.8M 0s
    ##  77400K .......... .......... .......... .......... .......... 96%  102M 0s
    ##  77450K .......... .......... .......... .......... .......... 96%  116M 0s
    ##  77500K .......... .......... .......... .......... .......... 97% 84.2M 0s
    ##  77550K .......... .......... .......... .......... .......... 97% 97.7M 0s
    ##  77600K .......... .......... .......... .......... .......... 97%  124M 0s
    ##  77650K .......... .......... .......... .......... .......... 97% 48.2M 0s
    ##  77700K .......... .......... .......... .......... .......... 97% 88.8M 0s
    ##  77750K .......... .......... .......... .......... .......... 97%  102M 0s
    ##  77800K .......... .......... .......... .......... .......... 97%  115M 0s
    ##  77850K .......... .......... .......... .......... .......... 97% 67.9M 0s
    ##  77900K .......... .......... .......... .......... .......... 97%  102M 0s
    ##  77950K .......... .......... .......... .......... .......... 97%  100M 0s
    ##  78000K .......... .......... .......... .......... .......... 97%  120M 0s
    ##  78050K .......... .......... .......... .......... .......... 97% 70.4M 0s
    ##  78100K .......... .......... .......... .......... .......... 97%  101M 0s
    ##  78150K .......... .......... .......... .......... .......... 97%  160M 0s
    ##  78200K .......... .......... .......... .......... .......... 97% 56.8M 0s
    ##  78250K .......... .......... .......... .......... .......... 97% 80.3M 0s
    ##  78300K .......... .......... .......... .......... .......... 98%  124M 0s
    ##  78350K .......... .......... .......... .......... .......... 98%  152M 0s
    ##  78400K .......... .......... .......... .......... .......... 98% 68.3M 0s
    ##  78450K .......... .......... .......... .......... .......... 98% 76.4M 0s
    ##  78500K .......... .......... .......... .......... .......... 98%  132M 0s
    ##  78550K .......... .......... .......... .......... .......... 98%  111M 0s
    ##  78600K .......... .......... .......... .......... .......... 98% 88.3M 0s
    ##  78650K .......... .......... .......... .......... .......... 98% 90.7M 0s
    ##  78700K .......... .......... .......... .......... .......... 98% 55.2M 0s
    ##  78750K .......... .......... .......... .......... .......... 98%  107M 0s
    ##  78800K .......... .......... .......... .......... .......... 98%  159M 0s
    ##  78850K .......... .......... .......... .......... .......... 98%  102M 0s
    ##  78900K .......... .......... .......... .......... .......... 98% 63.4M 0s
    ##  78950K .......... .......... .......... .......... .......... 98%  149M 0s
    ##  79000K .......... .......... .......... .......... .......... 98% 63.6M 0s
    ##  79050K .......... .......... .......... .......... .......... 98% 72.2M 0s
    ##  79100K .......... .......... .......... .......... .......... 99%  144M 0s
    ##  79150K .......... .......... .......... .......... .......... 99%  138M 0s
    ##  79200K .......... .......... .......... .......... .......... 99% 80.8M 0s
    ##  79250K .......... .......... .......... .......... .......... 99%  160M 0s
    ##  79300K .......... .......... .......... .......... .......... 99%  105M 0s
    ##  79350K .......... .......... .......... .......... .......... 99% 63.0M 0s
    ##  79400K .......... .......... .......... .......... .......... 99% 85.6M 0s
    ##  79450K .......... .......... .......... .......... .......... 99%  122M 0s
    ##  79500K .......... .......... .......... .......... .......... 99% 60.9M 0s
    ##  79550K .......... .......... .......... .......... .......... 99% 72.6M 0s
    ##  79600K .......... .......... .......... .......... .......... 99% 75.9M 0s
    ##  79650K .......... .......... .......... .......... .......... 99%  128M 0s
    ##  79700K .......... .......... .......... .......... .......... 99%  117M 0s
    ##  79750K .......... .......... .......... .......... .......... 99%  104M 0s
    ##  79800K .......... .......... .......... .......... .......... 99%  126M 0s
    ##  79850K .......... .......... .......... .......... .......... 99%  158M 0s
    ##  79900K .......... .......... ..                              100%  132M=1.1s
    ## 
    ## 2020-11-26 19:48:41 (72.6 MB/s) - ‘silva_species_assignment_v138.fa.gz.1’ saved [81840166/81840166]

## Ici nous créons une variable qui va recevoir les espèces obtenues grâce à Silva

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
```

``` r
taxa <- addSpecies(taxa, "~/silva_species_assignment_v138.fa.gz")
```

## Ici on remarque donc après avoir affiché la table qu’on a créée on obtient une majorité de Bacteroidetes ce qui est normal dans des échantillons fécaux. D’autres espèces n’ont pas pu être assignées car on a peu de données sur les bactéries des intestins des souris.

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

\#Evaluate accuracy \#\# Ici on cherche à comparer les variants donnés
par la machine par rapport à la composition de communauté attendue. Le
premier code nous montre qu’on s’attendait à avoir 20 souches
différentes.

``` r
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
```

    ## DADA2 inferred 20 sample sequences present in the Mock community.

## Dada2 a indentifié 20 ASVs. Cela correspond précisement à ce à quoi on s’attendait, donc le taux d’erreur résiduel par dada2 est de 0%

``` r
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```

    ## Of those, 20 were exact matches to the expected reference sequences.

\#Sauvegarde de l’environnement pour pouvoir jouer le code de phyloseq

``` r
save.image(file="03_DADA2_tutorial_FinalEnv")
```
