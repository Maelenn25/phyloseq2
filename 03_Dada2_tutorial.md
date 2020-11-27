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
  - [Inspect read quality profiles](#inspect-read-quality-profiles)
      - [Permet d’inspecter la qualité des différents reads et nous
        permet de voir à partir de quel nucléotide nous allons avoir une
        baisse du score de qualité pour pouvoir par la suite retirer les
        nucléotides de mauvaise qualité. Sur les forwards nous allons
        retirer les 10 derniers nucléotides avec une commande retrouvée
        plus bas dans le
        script.](#permet-dinspecter-la-qualité-des-différents-reads-et-nous-permet-de-voir-à-partir-de-quel-nucléotide-nous-allons-avoir-une-baisse-du-score-de-qualité-pour-pouvoir-par-la-suite-retirer-les-nucléotides-de-mauvaise-qualité.-sur-les-forwards-nous-allons-retirer-les-10-derniers-nucléotides-avec-une-commande-retrouvée-plus-bas-dans-le-script.)
  - [Nous allons indiquer à la machine quels paramètres nous allons
    utiliser pour filtrer les séquences avant de les ranger. On indique
    ainsi dans la première fonction les deux fichiers d’où les séquences
    viendront (fnFs et fnRs), puis le fichier où elles seront rangées
    (filtFs et filtRs crées juste au dessus). Enuite nous allons
    indiquer où couper pour les deux sortes de
    séquences.](#nous-allons-indiquer-à-la-machine-quels-paramètres-nous-allons-utiliser-pour-filtrer-les-séquences-avant-de-les-ranger.-on-indique-ainsi-dans-la-première-fonction-les-deux-fichiers-doù-les-séquences-viendront-fnfs-et-fnrs-puis-le-fichier-où-elles-seront-rangées-filtfs-et-filtrs-crées-juste-au-dessus.-enuite-nous-allons-indiquer-où-couper-pour-les-deux-sortes-de-séquences.)
  - [Learn the Error Rates](#learn-the-error-rates)
      - [Nous allons ici utiliser des lignes de commandes qui vont
        permettre d’apprendre à la machine les différents profils
        d’erreurs générées lors du séquençage. L’opération est faite
        sur les deux types de
        séquence.](#nous-allons-ici-utiliser-des-lignes-de-commandes-qui-vont-permettre-dapprendre-à-la-machine-les-différents-profils-derreurs-générées-lors-du-séquençage.-lopération-est-faite-sur-les-deux-types-de-séquence.)
      - [Cette ligne nous permet de visualiser les erreurs qu’on a fait
        apprendre à la
        machine](#cette-ligne-nous-permet-de-visualiser-les-erreurs-quon-a-fait-apprendre-à-la-machine)
  - [Sample Inference](#sample-inference)
      - [Ici nous créons une autre variable “dadaFs” dans laquelle nous
        mettons les fichiers obtenus après avoir filtré et appliqué le
        profil d’erreur à nos séquences. Nous allons faire la même chose
        avec
        dadaRS.](#ici-nous-créons-une-autre-variable-dadafs-dans-laquelle-nous-mettons-les-fichiers-obtenus-après-avoir-filtré-et-appliqué-le-profil-derreur-à-nos-séquences.-nous-allons-faire-la-même-chose-avec-dadars.)
      - [Cette commande nous permet de visualiser le résultat global
        qu’on retrouve classé dans la liste dadaFs. Ils nous indiquent
        que sur les séquences on retrouve 128 séquences qui
        correspondent aux vrais variants, par rapport aux 1979
        séquences. Ils nous indiquent aussi les diagnostiques de
        qualité.](#cette-commande-nous-permet-de-visualiser-le-résultat-global-quon-retrouve-classé-dans-la-liste-dadafs.-ils-nous-indiquent-que-sur-les-séquences-on-retrouve-128-séquences-qui-correspondent-aux-vrais-variants-par-rapport-aux-1979-séquences.-ils-nous-indiquent-aussi-les-diagnostiques-de-qualité.)
  - [Merge paired reads](#merge-paired-reads)
      - [Ici nous voulons mettre en une seule séquence les Forwards et
        les Reverses.Nous pouvons faire cette opération grâce aux
        overlaps de 12 paires de base. Cela se fait grâce à un
        alignement entre les forwards et les reverses qui vont permettre
        de contruire les
        contigs.](#ici-nous-voulons-mettre-en-une-seule-séquence-les-forwards-et-les-reverses.nous-pouvons-faire-cette-opération-grâce-aux-overlaps-de-12-paires-de-base.-cela-se-fait-grâce-à-un-alignement-entre-les-forwards-et-les-reverses-qui-vont-permettre-de-contruire-les-contigs.)
  - [Construct sequence table](#construct-sequence-table)
      - [Ici nous allons construire une table des variations de séquence
        dans les amplicons (ASV) qui permet une meilleure résolution que
        les tables OTUs
        97%](#ici-nous-allons-construire-une-table-des-variations-de-séquence-dans-les-amplicons-asv-qui-permet-une-meilleure-résolution-que-les-tables-otus-97)
  - [Remove chimeras](#remove-chimeras)
      - [Malgré qu’on ait pu appliquer les modèles d’erreurs aux
        séquences, il reste des chimères. Ces chimères sont facilement
        reconnaissables par la machine et peuvent etre réparées en y
        rajoutant les parties droites et gauche des 2 séquences les plus
        abondantes.](#malgré-quon-ait-pu-appliquer-les-modèles-derreurs-aux-séquences-il-reste-des-chimères.-ces-chimères-sont-facilement-reconnaissables-par-la-machine-et-peuvent-etre-réparées-en-y-rajoutant-les-parties-droites-et-gauche-des-2-séquences-les-plus-abondantes.)
      - [Ici on peut voir qu’on à 3% de séquences chimériques dans notre
        jeu de
        donnée.](#ici-on-peut-voir-quon-à-3-de-séquences-chimériques-dans-notre-jeu-de-donnée.)
  - [Track reads through the
    pipeline](#track-reads-through-the-pipeline)
      - [Ce code nous permet de visualiser le nombre de séquences
        obtenues à la suite de toutes nos manipulations de filtrage. Ici
        nous pouvons voir qu’on a pu récupérer la plupart de nos
        séquences brutes, ce qui est signe d’une bonne qualité de
        séquençage.](#ce-code-nous-permet-de-visualiser-le-nombre-de-séquences-obtenues-à-la-suite-de-toutes-nos-manipulations-de-filtrage.-ici-nous-pouvons-voir-quon-a-pu-récupérer-la-plupart-de-nos-séquences-brutes-ce-qui-est-signe-dune-bonne-qualité-de-séquençage.)
  - [Assign taxonomy](#assign-taxonomy)
      - [Ici nous avons du récupérer silva afin d’analyser et d’assigner
        les
        taxonomies.](#ici-nous-avons-du-récupérer-silva-afin-danalyser-et-dassigner-les-taxonomies.)
      - [Ici nous créons une variable qui va recevoir les espèces
        obtenues grâce à
        Silva](#ici-nous-créons-une-variable-qui-va-recevoir-les-espèces-obtenues-grâce-à-silva)
      - [Ici on remarque donc après avoir affiché la table qu’on a créée
        on obtient une majorité de Bacteroidetes ce qui est normal dans
        des échantillons fécaux. D’autres espèces n’ont pas pu être
        assignées car on a peu de données sur les bactéries des
        intestins des
        souris.](#ici-on-remarque-donc-après-avoir-affiché-la-table-quon-a-créée-on-obtient-une-majorité-de-bacteroidetes-ce-qui-est-normal-dans-des-échantillons-fécaux.-dautres-espèces-nont-pas-pu-être-assignées-car-on-a-peu-de-données-sur-les-bactéries-des-intestins-des-souris.)
  - [Evaluate accuracy](#evaluate-accuracy)
      - [Ici on cherche à comparer les variants donnés par la machine
        par rapport à la composition de communauté attendue. Le premier
        code nous montre qu’on s’attendait à avoir 20 souches
        différentes.](#ici-on-cherche-à-comparer-les-variants-donnés-par-la-machine-par-rapport-à-la-composition-de-communauté-attendue.-le-premier-code-nous-montre-quon-sattendait-à-avoir-20-souches-différentes.)
      - [Dada2 a indentifié 20 ASVs. Cela correspond précisement à ce à
        quoi on s’attendait, donc le taux d’erreur résiduel par dada2
        est de
        0%](#dada2-a-indentifié-20-asvs.-cela-correspond-précisement-à-ce-à-quoi-on-sattendait-donc-le-taux-derreur-résiduel-par-dada2-est-de-0)
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

## Création de deux listes. Une avec les Forwards et une avec les Reverses qui correspondent respectivement au Read1 et Read2. La commande avec ‘sample.names’ permet de mettre toutes les séquences sous un même format pour harmoniser leurs noms.

``` r
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

# Inspect read quality profiles

## Permet d’inspecter la qualité des différents reads et nous permet de voir à partir de quel nucléotide nous allons avoir une baisse du score de qualité pour pouvoir par la suite retirer les nucléotides de mauvaise qualité. Sur les forwards nous allons retirer les 10 derniers nucléotides avec une commande retrouvée plus bas dans le script.

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

![](03_Dada2_tutorial_files/figure-gfm/unnamed-chunk-5-1.png)<!-- --> \#
Filter and trim \#\# Ici nous créeons une nouvelle variable qui recevra
tous les fichiers filtrés, que ce soit pour les Forwards ou les
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

# Learn the Error Rates

## Nous allons ici utiliser des lignes de commandes qui vont permettre d’apprendre à la machine les différents profils d’erreurs générées lors du séquençage. L’opération est faite sur les deux types de séquence.

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

# Sample Inference

## Ici nous créons une autre variable “dadaFs” dans laquelle nous mettons les fichiers obtenus après avoir filtré et appliqué le profil d’erreur à nos séquences. Nous allons faire la même chose avec dadaRS.

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

# Merge paired reads

## Ici nous voulons mettre en une seule séquence les Forwards et les Reverses.Nous pouvons faire cette opération grâce aux overlaps de 12 paires de base. Cela se fait grâce à un alignement entre les forwards et les reverses qui vont permettre de contruire les contigs.

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

## Ici nous allons construire une table des variations de séquence dans les amplicons (ASV) qui permet une meilleure résolution que les tables OTUs 97%

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

## Malgré qu’on ait pu appliquer les modèles d’erreurs aux séquences, il reste des chimères. Ces chimères sont facilement reconnaissables par la machine et peuvent etre réparées en y rajoutant les parties droites et gauche des 2 séquences les plus abondantes.

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

# Track reads through the pipeline

## Ce code nous permet de visualiser le nombre de séquences obtenues à la suite de toutes nos manipulations de filtrage. Ici nous pouvons voir qu’on a pu récupérer la plupart de nos séquences brutes, ce qui est signe d’une bonne qualité de séquençage.

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

## Ici nous avons du récupérer silva afin d’analyser et d’assigner les taxonomies.

``` bash
wget https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz
```

    ## --2020-11-27 13:48:16--  https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz
    ## Resolving zenodo.org (zenodo.org)... 137.138.76.77
    ## Connecting to zenodo.org (zenodo.org)|137.138.76.77|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 81840166 (78M) [application/octet-stream]
    ## Saving to: ‘silva_species_assignment_v138.fa.gz.3’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0% 7.20M 11s
    ##     50K .......... .......... .......... .......... ..........  0% 12.7M 8s
    ##    100K .......... .......... .......... .......... ..........  0% 8.46M 9s
    ##    150K .......... .......... .......... .......... ..........  0% 13.5M 8s
    ##    200K .......... .......... .......... .......... ..........  0% 74.8M 7s
    ##    250K .......... .......... .......... .......... ..........  0% 39.2M 6s
    ##    300K .......... .......... .......... .......... ..........  0% 17.6M 6s
    ##    350K .......... .......... .......... .......... ..........  0% 67.2M 5s
    ##    400K .......... .......... .......... .......... ..........  0% 18.4M 5s
    ##    450K .......... .......... .......... .......... ..........  0% 36.0M 5s
    ##    500K .......... .......... .......... .......... ..........  0% 26.7M 5s
    ##    550K .......... .......... .......... .......... ..........  0% 36.8M 4s
    ##    600K .......... .......... .......... .......... ..........  0% 60.4M 4s
    ##    650K .......... .......... .......... .......... ..........  0% 29.3M 4s
    ##    700K .......... .......... .......... .......... ..........  0% 22.2M 4s
    ##    750K .......... .......... .......... .......... ..........  1% 27.4M 4s
    ##    800K .......... .......... .......... .......... ..........  1% 51.7M 4s
    ##    850K .......... .......... .......... .......... ..........  1% 59.1M 4s
    ##    900K .......... .......... .......... .......... ..........  1% 17.0M 4s
    ##    950K .......... .......... .......... .......... ..........  1% 46.7M 4s
    ##   1000K .......... .......... .......... .......... ..........  1% 16.4M 4s
    ##   1050K .......... .......... .......... .......... ..........  1% 53.3M 3s
    ##   1100K .......... .......... .......... .......... ..........  1% 79.1M 3s
    ##   1150K .......... .......... .......... .......... ..........  1% 18.5M 3s
    ##   1200K .......... .......... .......... .......... ..........  1% 69.9M 3s
    ##   1250K .......... .......... .......... .......... ..........  1% 78.1M 3s
    ##   1300K .......... .......... .......... .......... ..........  1% 16.8M 3s
    ##   1350K .......... .......... .......... .......... ..........  1% 64.9M 3s
    ##   1400K .......... .......... .......... .......... ..........  1% 12.2M 3s
    ##   1450K .......... .......... .......... .......... ..........  1% 65.7M 3s
    ##   1500K .......... .......... .......... .......... ..........  1% 15.9M 3s
    ##   1550K .......... .......... .......... .......... ..........  2% 51.3M 3s
    ##   1600K .......... .......... .......... .......... ..........  2% 89.6M 3s
    ##   1650K .......... .......... .......... .......... ..........  2% 96.4M 3s
    ##   1700K .......... .......... .......... .......... ..........  2% 24.5M 3s
    ##   1750K .......... .......... .......... .......... ..........  2% 18.7M 3s
    ##   1800K .......... .......... .......... .......... ..........  2% 58.3M 3s
    ##   1850K .......... .......... .......... .......... ..........  2% 58.3M 3s
    ##   1900K .......... .......... .......... .......... ..........  2% 21.7M 3s
    ##   1950K .......... .......... .......... .......... ..........  2% 55.3M 3s
    ##   2000K .......... .......... .......... .......... ..........  2% 53.3M 3s
    ##   2050K .......... .......... .......... .......... ..........  2% 20.0M 3s
    ##   2100K .......... .......... .......... .......... ..........  2% 33.8M 3s
    ##   2150K .......... .......... .......... .......... ..........  2% 28.0M 3s
    ##   2200K .......... .......... .......... .......... ..........  2% 27.1M 3s
    ##   2250K .......... .......... .......... .......... ..........  2% 68.3M 3s
    ##   2300K .......... .......... .......... .......... ..........  2% 28.8M 3s
    ##   2350K .......... .......... .......... .......... ..........  3% 23.7M 3s
    ##   2400K .......... .......... .......... .......... ..........  3% 31.9M 3s
    ##   2450K .......... .......... .......... .......... ..........  3% 19.5M 3s
    ##   2500K .......... .......... .......... .......... ..........  3% 89.7M 3s
    ##   2550K .......... .......... .......... .......... ..........  3% 36.3M 3s
    ##   2600K .......... .......... .......... .......... ..........  3% 23.2M 3s
    ##   2650K .......... .......... .......... .......... ..........  3% 27.5M 3s
    ##   2700K .......... .......... .......... .......... ..........  3% 77.2M 3s
    ##   2750K .......... .......... .......... .......... ..........  3% 30.1M 3s
    ##   2800K .......... .......... .......... .......... ..........  3% 25.2M 3s
    ##   2850K .......... .......... .......... .......... ..........  3% 24.4M 3s
    ##   2900K .......... .......... .......... .......... ..........  3% 30.0M 3s
    ##   2950K .......... .......... .......... .......... ..........  3% 68.7M 3s
    ##   3000K .......... .......... .......... .......... ..........  3% 23.6M 3s
    ##   3050K .......... .......... .......... .......... ..........  3% 32.6M 3s
    ##   3100K .......... .......... .......... .......... ..........  3% 22.2M 3s
    ##   3150K .......... .......... .......... .......... ..........  4% 33.0M 3s
    ##   3200K .......... .......... .......... .......... ..........  4% 73.5M 3s
    ##   3250K .......... .......... .......... .......... ..........  4% 29.0M 3s
    ##   3300K .......... .......... .......... .......... ..........  4% 25.9M 3s
    ##   3350K .......... .......... .......... .......... ..........  4% 29.6M 3s
    ##   3400K .......... .......... .......... .......... ..........  4% 24.5M 3s
    ##   3450K .......... .......... .......... .......... ..........  4% 76.8M 3s
    ##   3500K .......... .......... .......... .......... ..........  4% 37.9M 3s
    ##   3550K .......... .......... .......... .......... ..........  4% 20.1M 3s
    ##   3600K .......... .......... .......... .......... ..........  4% 44.4M 3s
    ##   3650K .......... .......... .......... .......... ..........  4%  101M 3s
    ##   3700K .......... .......... .......... .......... ..........  4% 16.9M 3s
    ##   3750K .......... .......... .......... .......... ..........  4% 74.4M 3s
    ##   3800K .......... .......... .......... .......... ..........  4% 13.3M 3s
    ##   3850K .......... .......... .......... .......... ..........  4% 68.5M 3s
    ##   3900K .......... .......... .......... .......... ..........  4% 38.2M 3s
    ##   3950K .......... .......... .......... .......... ..........  5%  105M 3s
    ##   4000K .......... .......... .......... .......... ..........  5% 49.7M 3s
    ##   4050K .......... .......... .......... .......... ..........  5% 29.2M 3s
    ##   4100K .......... .......... .......... .......... ..........  5% 73.5M 3s
    ##   4150K .......... .......... .......... .......... ..........  5%  104M 3s
    ##   4200K .......... .......... .......... .......... ..........  5% 50.8M 3s
    ##   4250K .......... .......... .......... .......... ..........  5% 31.3M 3s
    ##   4300K .......... .......... .......... .......... ..........  5% 49.5M 3s
    ##   4350K .......... .......... .......... .......... ..........  5% 75.3M 3s
    ##   4400K .......... .......... .......... .......... ..........  5% 70.3M 2s
    ##   4450K .......... .......... .......... .......... ..........  5% 28.0M 2s
    ##   4500K .......... .......... .......... .......... ..........  5% 50.4M 2s
    ##   4550K .......... .......... .......... .......... ..........  5% 65.3M 2s
    ##   4600K .......... .......... .......... .......... ..........  5% 56.6M 2s
    ##   4650K .......... .......... .......... .......... ..........  5% 73.5M 2s
    ##   4700K .......... .......... .......... .......... ..........  5% 26.4M 2s
    ##   4750K .......... .......... .......... .......... ..........  6% 58.0M 2s
    ##   4800K .......... .......... .......... .......... ..........  6% 54.2M 2s
    ##   4850K .......... .......... .......... .......... ..........  6% 79.5M 2s
    ##   4900K .......... .......... .......... .......... ..........  6% 28.0M 2s
    ##   4950K .......... .......... .......... .......... ..........  6% 36.3M 2s
    ##   5000K .......... .......... .......... .......... ..........  6% 74.1M 2s
    ##   5050K .......... .......... .......... .......... ..........  6%  100M 2s
    ##   5100K .......... .......... .......... .......... ..........  6% 38.4M 2s
    ##   5150K .......... .......... .......... .......... ..........  6% 45.0M 2s
    ##   5200K .......... .......... .......... .......... ..........  6% 44.8M 2s
    ##   5250K .......... .......... .......... .......... ..........  6% 26.7M 2s
    ##   5300K .......... .......... .......... .......... ..........  6% 63.6M 2s
    ##   5350K .......... .......... .......... .......... ..........  6% 62.7M 2s
    ##   5400K .......... .......... .......... .......... ..........  6% 53.1M 2s
    ##   5450K .......... .......... .......... .......... ..........  6% 42.9M 2s
    ##   5500K .......... .......... .......... .......... ..........  6% 37.9M 2s
    ##   5550K .......... .......... .......... .......... ..........  7% 77.7M 2s
    ##   5600K .......... .......... .......... .......... ..........  7% 60.3M 2s
    ##   5650K .......... .......... .......... .......... ..........  7% 34.7M 2s
    ##   5700K .......... .......... .......... .......... ..........  7% 59.5M 2s
    ##   5750K .......... .......... .......... .......... ..........  7% 67.7M 2s
    ##   5800K .......... .......... .......... .......... ..........  7% 76.2M 2s
    ##   5850K .......... .......... .......... .......... ..........  7% 25.7M 2s
    ##   5900K .......... .......... .......... .......... ..........  7% 23.9M 2s
    ##   5950K .......... .......... .......... .......... ..........  7%  110M 2s
    ##   6000K .......... .......... .......... .......... ..........  7% 85.5M 2s
    ##   6050K .......... .......... .......... .......... ..........  7% 70.4M 2s
    ##   6100K .......... .......... .......... .......... ..........  7% 20.8M 2s
    ##   6150K .......... .......... .......... .......... ..........  7% 49.9M 2s
    ##   6200K .......... .......... .......... .......... ..........  7% 90.5M 2s
    ##   6250K .......... .......... .......... .......... ..........  7%  129M 2s
    ##   6300K .......... .......... .......... .......... ..........  7% 40.3M 2s
    ##   6350K .......... .......... .......... .......... ..........  8% 25.6M 2s
    ##   6400K .......... .......... .......... .......... ..........  8% 85.9M 2s
    ##   6450K .......... .......... .......... .......... ..........  8% 79.0M 2s
    ##   6500K .......... .......... .......... .......... ..........  8%  107M 2s
    ##   6550K .......... .......... .......... .......... ..........  8% 79.1M 2s
    ##   6600K .......... .......... .......... .......... ..........  8% 53.4M 2s
    ##   6650K .......... .......... .......... .......... ..........  8% 50.0M 2s
    ##   6700K .......... .......... .......... .......... ..........  8% 31.2M 2s
    ##   6750K .......... .......... .......... .......... ..........  8% 85.3M 2s
    ##   6800K .......... .......... .......... .......... ..........  8% 70.0M 2s
    ##   6850K .......... .......... .......... .......... ..........  8% 72.7M 2s
    ##   6900K .......... .......... .......... .......... ..........  8% 31.7M 2s
    ##   6950K .......... .......... .......... .......... ..........  8% 85.8M 2s
    ##   7000K .......... .......... .......... .......... ..........  8% 58.3M 2s
    ##   7050K .......... .......... .......... .......... ..........  8% 91.0M 2s
    ##   7100K .......... .......... .......... .......... ..........  8% 89.3M 2s
    ##   7150K .......... .......... .......... .......... ..........  9% 23.7M 2s
    ##   7200K .......... .......... .......... .......... ..........  9% 70.5M 2s
    ##   7250K .......... .......... .......... .......... ..........  9% 77.6M 2s
    ##   7300K .......... .......... .......... .......... ..........  9% 64.5M 2s
    ##   7350K .......... .......... .......... .......... ..........  9% 37.2M 2s
    ##   7400K .......... .......... .......... .......... ..........  9% 56.7M 2s
    ##   7450K .......... .......... .......... .......... ..........  9%  106M 2s
    ##   7500K .......... .......... .......... .......... ..........  9% 74.6M 2s
    ##   7550K .......... .......... .......... .......... ..........  9% 27.3M 2s
    ##   7600K .......... .......... .......... .......... ..........  9% 45.5M 2s
    ##   7650K .......... .......... .......... .......... ..........  9% 94.7M 2s
    ##   7700K .......... .......... .......... .......... ..........  9% 89.9M 2s
    ##   7750K .......... .......... .......... .......... ..........  9% 56.9M 2s
    ##   7800K .......... .......... .......... .......... ..........  9% 35.5M 2s
    ##   7850K .......... .......... .......... .......... ..........  9% 43.3M 2s
    ##   7900K .......... .......... .......... .......... ..........  9% 73.9M 2s
    ##   7950K .......... .......... .......... .......... .......... 10% 98.2M 2s
    ##   8000K .......... .......... .......... .......... .......... 10% 30.0M 2s
    ##   8050K .......... .......... .......... .......... .......... 10%  104M 2s
    ##   8100K .......... .......... .......... .......... .......... 10% 45.3M 2s
    ##   8150K .......... .......... .......... .......... .......... 10% 72.8M 2s
    ##   8200K .......... .......... .......... .......... .......... 10% 40.2M 2s
    ##   8250K .......... .......... .......... .......... .......... 10% 82.9M 2s
    ##   8300K .......... .......... .......... .......... .......... 10% 22.8M 2s
    ##   8350K .......... .......... .......... .......... .......... 10% 74.7M 2s
    ##   8400K .......... .......... .......... .......... .......... 10% 39.3M 2s
    ##   8450K .......... .......... .......... .......... .......... 10% 56.5M 2s
    ##   8500K .......... .......... .......... .......... .......... 10% 2.54M 2s
    ##   8550K .......... .......... .......... .......... .......... 10% 45.7M 2s
    ##   8600K .......... .......... .......... .......... .......... 10% 62.4M 2s
    ##   8650K .......... .......... .......... .......... .......... 10% 76.8M 2s
    ##   8700K .......... .......... .......... .......... .......... 10% 52.8M 2s
    ##   8750K .......... .......... .......... .......... .......... 11% 19.7M 2s
    ##   8800K .......... .......... .......... .......... .......... 11% 47.9M 2s
    ##   8850K .......... .......... .......... .......... .......... 11% 51.9M 2s
    ##   8900K .......... .......... .......... .......... .......... 11% 83.3M 2s
    ##   8950K .......... .......... .......... .......... .......... 11% 54.4M 2s
    ##   9000K .......... .......... .......... .......... .......... 11% 40.6M 2s
    ##   9050K .......... .......... .......... .......... .......... 11% 56.8M 2s
    ##   9100K .......... .......... .......... .......... .......... 11% 42.5M 2s
    ##   9150K .......... .......... .......... .......... .......... 11% 64.8M 2s
    ##   9200K .......... .......... .......... .......... .......... 11% 67.6M 2s
    ##   9250K .......... .......... .......... .......... .......... 11% 87.1M 2s
    ##   9300K .......... .......... .......... .......... .......... 11% 49.9M 2s
    ##   9350K .......... .......... .......... .......... .......... 11% 59.3M 2s
    ##   9400K .......... .......... .......... .......... .......... 11% 53.3M 2s
    ##   9450K .......... .......... .......... .......... .......... 11% 41.9M 2s
    ##   9500K .......... .......... .......... .......... .......... 11% 50.4M 2s
    ##   9550K .......... .......... .......... .......... .......... 12% 86.9M 2s
    ##   9600K .......... .......... .......... .......... .......... 12% 65.9M 2s
    ##   9650K .......... .......... .......... .......... .......... 12% 20.4M 2s
    ##   9700K .......... .......... .......... .......... .......... 12% 64.1M 2s
    ##   9750K .......... .......... .......... .......... .......... 12% 72.0M 2s
    ##   9800K .......... .......... .......... .......... .......... 12% 66.5M 2s
    ##   9850K .......... .......... .......... .......... .......... 12% 85.3M 2s
    ##   9900K .......... .......... .......... .......... .......... 12% 39.3M 2s
    ##   9950K .......... .......... .......... .......... .......... 12% 47.4M 2s
    ##  10000K .......... .......... .......... .......... .......... 12% 27.2M 2s
    ##  10050K .......... .......... .......... .......... .......... 12% 74.8M 2s
    ##  10100K .......... .......... .......... .......... .......... 12% 88.3M 2s
    ##  10150K .......... .......... .......... .......... .......... 12% 84.4M 2s
    ##  10200K .......... .......... .......... .......... .......... 12% 56.9M 2s
    ##  10250K .......... .......... .......... .......... .......... 12% 60.9M 2s
    ##  10300K .......... .......... .......... .......... .......... 12% 36.7M 2s
    ##  10350K .......... .......... .......... .......... .......... 13% 68.8M 2s
    ##  10400K .......... .......... .......... .......... .......... 13% 71.0M 2s
    ##  10450K .......... .......... .......... .......... .......... 13% 63.5M 2s
    ##  10500K .......... .......... .......... .......... .......... 13% 55.9M 2s
    ##  10550K .......... .......... .......... .......... .......... 13% 34.4M 2s
    ##  10600K .......... .......... .......... .......... .......... 13% 66.8M 2s
    ##  10650K .......... .......... .......... .......... .......... 13% 69.6M 2s
    ##  10700K .......... .......... .......... .......... .......... 13% 79.5M 2s
    ##  10750K .......... .......... .......... .......... .......... 13% 67.9M 2s
    ##  10800K .......... .......... .......... .......... .......... 13% 57.4M 2s
    ##  10850K .......... .......... .......... .......... .......... 13% 40.9M 2s
    ##  10900K .......... .......... .......... .......... .......... 13% 63.7M 2s
    ##  10950K .......... .......... .......... .......... .......... 13% 72.3M 2s
    ##  11000K .......... .......... .......... .......... .......... 13% 50.7M 2s
    ##  11050K .......... .......... .......... .......... .......... 13% 95.4M 2s
    ##  11100K .......... .......... .......... .......... .......... 13% 78.8M 2s
    ##  11150K .......... .......... .......... .......... .......... 14% 55.6M 2s
    ##  11200K .......... .......... .......... .......... .......... 14% 58.5M 2s
    ##  11250K .......... .......... .......... .......... .......... 14% 72.8M 2s
    ##  11300K .......... .......... .......... .......... .......... 14% 75.6M 2s
    ##  11350K .......... .......... .......... .......... .......... 14% 68.0M 2s
    ##  11400K .......... .......... .......... .......... .......... 14% 79.5M 2s
    ##  11450K .......... .......... .......... .......... .......... 14% 33.1M 2s
    ##  11500K .......... .......... .......... .......... .......... 14% 67.2M 2s
    ##  11550K .......... .......... .......... .......... .......... 14% 60.0M 2s
    ##  11600K .......... .......... .......... .......... .......... 14% 62.7M 2s
    ##  11650K .......... .......... .......... .......... .......... 14% 19.8M 2s
    ##  11700K .......... .......... .......... .......... .......... 14% 72.7M 2s
    ##  11750K .......... .......... .......... .......... .......... 14% 58.1M 2s
    ##  11800K .......... .......... .......... .......... .......... 14% 58.9M 2s
    ##  11850K .......... .......... .......... .......... .......... 14% 70.4M 2s
    ##  11900K .......... .......... .......... .......... .......... 14% 66.7M 2s
    ##  11950K .......... .......... .......... .......... .......... 15% 51.6M 2s
    ##  12000K .......... .......... .......... .......... .......... 15% 62.8M 2s
    ##  12050K .......... .......... .......... .......... .......... 15% 55.2M 2s
    ##  12100K .......... .......... .......... .......... .......... 15% 64.9M 2s
    ##  12150K .......... .......... .......... .......... .......... 15% 53.3M 2s
    ##  12200K .......... .......... .......... .......... .......... 15% 58.8M 2s
    ##  12250K .......... .......... .......... .......... .......... 15% 62.5M 2s
    ##  12300K .......... .......... .......... .......... .......... 15% 61.5M 2s
    ##  12350K .......... .......... .......... .......... .......... 15% 76.6M 2s
    ##  12400K .......... .......... .......... .......... .......... 15% 61.7M 2s
    ##  12450K .......... .......... .......... .......... .......... 15% 70.3M 2s
    ##  12500K .......... .......... .......... .......... .......... 15% 85.5M 2s
    ##  12550K .......... .......... .......... .......... .......... 15% 68.7M 2s
    ##  12600K .......... .......... .......... .......... .......... 15% 80.4M 2s
    ##  12650K .......... .......... .......... .......... .......... 15% 47.1M 2s
    ##  12700K .......... .......... .......... .......... .......... 15% 62.1M 2s
    ##  12750K .......... .......... .......... .......... .......... 16% 54.3M 2s
    ##  12800K .......... .......... .......... .......... .......... 16% 73.0M 2s
    ##  12850K .......... .......... .......... .......... .......... 16%  110M 2s
    ##  12900K .......... .......... .......... .......... .......... 16% 50.7M 2s
    ##  12950K .......... .......... .......... .......... .......... 16% 61.7M 2s
    ##  13000K .......... .......... .......... .......... .......... 16% 82.1M 2s
    ##  13050K .......... .......... .......... .......... .......... 16% 76.5M 2s
    ##  13100K .......... .......... .......... .......... .......... 16% 76.1M 2s
    ##  13150K .......... .......... .......... .......... .......... 16% 57.8M 2s
    ##  13200K .......... .......... .......... .......... .......... 16% 48.9M 2s
    ##  13250K .......... .......... .......... .......... .......... 16% 90.2M 2s
    ##  13300K .......... .......... .......... .......... .......... 16% 73.5M 2s
    ##  13350K .......... .......... .......... .......... .......... 16% 58.2M 2s
    ##  13400K .......... .......... .......... .......... .......... 16% 88.5M 2s
    ##  13450K .......... .......... .......... .......... .......... 16% 57.0M 2s
    ##  13500K .......... .......... .......... .......... .......... 16% 50.6M 2s
    ##  13550K .......... .......... .......... .......... .......... 17% 32.8M 2s
    ##  13600K .......... .......... .......... .......... .......... 17% 66.6M 2s
    ##  13650K .......... .......... .......... .......... .......... 17% 87.6M 2s
    ##  13700K .......... .......... .......... .......... .......... 17% 74.7M 2s
    ##  13750K .......... .......... .......... .......... .......... 17% 60.3M 2s
    ##  13800K .......... .......... .......... .......... .......... 17% 57.9M 2s
    ##  13850K .......... .......... .......... .......... .......... 17% 71.4M 2s
    ##  13900K .......... .......... .......... .......... .......... 17% 68.1M 2s
    ##  13950K .......... .......... .......... .......... .......... 17% 82.9M 2s
    ##  14000K .......... .......... .......... .......... .......... 17% 38.2M 2s
    ##  14050K .......... .......... .......... .......... .......... 17% 79.4M 2s
    ##  14100K .......... .......... .......... .......... .......... 17% 57.2M 2s
    ##  14150K .......... .......... .......... .......... .......... 17% 91.7M 2s
    ##  14200K .......... .......... .......... .......... .......... 17% 28.2M 2s
    ##  14250K .......... .......... .......... .......... .......... 17% 65.8M 2s
    ##  14300K .......... .......... .......... .......... .......... 17% 53.2M 2s
    ##  14350K .......... .......... .......... .......... .......... 18% 72.3M 2s
    ##  14400K .......... .......... .......... .......... .......... 18% 16.5M 2s
    ##  14450K .......... .......... .......... .......... .......... 18% 66.4M 2s
    ##  14500K .......... .......... .......... .......... .......... 18% 85.1M 2s
    ##  14550K .......... .......... .......... .......... .......... 18% 80.9M 2s
    ##  14600K .......... .......... .......... .......... .......... 18% 84.3M 2s
    ##  14650K .......... .......... .......... .......... .......... 18% 71.9M 2s
    ##  14700K .......... .......... .......... .......... .......... 18% 61.0M 2s
    ##  14750K .......... .......... .......... .......... .......... 18% 44.6M 2s
    ##  14800K .......... .......... .......... .......... .......... 18% 54.8M 2s
    ##  14850K .......... .......... .......... .......... .......... 18% 88.2M 2s
    ##  14900K .......... .......... .......... .......... .......... 18% 93.1M 2s
    ##  14950K .......... .......... .......... .......... .......... 18% 62.3M 2s
    ##  15000K .......... .......... .......... .......... .......... 18% 55.0M 2s
    ##  15050K .......... .......... .......... .......... .......... 18% 46.1M 2s
    ##  15100K .......... .......... .......... .......... .......... 18% 79.8M 2s
    ##  15150K .......... .......... .......... .......... .......... 19% 25.2M 2s
    ##  15200K .......... .......... .......... .......... .......... 19% 81.6M 2s
    ##  15250K .......... .......... .......... .......... .......... 19% 89.4M 2s
    ##  15300K .......... .......... .......... .......... .......... 19% 91.9M 2s
    ##  15350K .......... .......... .......... .......... .......... 19% 90.7M 2s
    ##  15400K .......... .......... .......... .......... .......... 19% 17.0M 2s
    ##  15450K .......... .......... .......... .......... .......... 19% 76.8M 2s
    ##  15500K .......... .......... .......... .......... .......... 19% 56.8M 2s
    ##  15550K .......... .......... .......... .......... .......... 19% 90.8M 2s
    ##  15600K .......... .......... .......... .......... .......... 19% 86.0M 2s
    ##  15650K .......... .......... .......... .......... .......... 19% 72.1M 2s
    ##  15700K .......... .......... .......... .......... .......... 19% 81.2M 2s
    ##  15750K .......... .......... .......... .......... .......... 19% 35.5M 2s
    ##  15800K .......... .......... .......... .......... .......... 19%  104M 2s
    ##  15850K .......... .......... .......... .......... .......... 19% 71.0M 2s
    ##  15900K .......... .......... .......... .......... .......... 19% 66.1M 1s
    ##  15950K .......... .......... .......... .......... .......... 20% 84.9M 1s
    ##  16000K .......... .......... .......... .......... .......... 20% 49.7M 1s
    ##  16050K .......... .......... .......... .......... .......... 20% 34.0M 1s
    ##  16100K .......... .......... .......... .......... .......... 20% 79.5M 1s
    ##  16150K .......... .......... .......... .......... .......... 20% 72.9M 1s
    ##  16200K .......... .......... .......... .......... .......... 20% 73.3M 1s
    ##  16250K .......... .......... .......... .......... .......... 20% 72.5M 1s
    ##  16300K .......... .......... .......... .......... .......... 20% 74.6M 1s
    ##  16350K .......... .......... .......... .......... .......... 20% 64.8M 1s
    ##  16400K .......... .......... .......... .......... .......... 20% 92.5M 1s
    ##  16450K .......... .......... .......... .......... .......... 20% 36.8M 1s
    ##  16500K .......... .......... .......... .......... .......... 20% 89.3M 1s
    ##  16550K .......... .......... .......... .......... .......... 20% 59.8M 1s
    ##  16600K .......... .......... .......... .......... .......... 20%  108M 1s
    ##  16650K .......... .......... .......... .......... .......... 20% 87.4M 1s
    ##  16700K .......... .......... .......... .......... .......... 20% 28.0M 1s
    ##  16750K .......... .......... .......... .......... .......... 21% 50.3M 1s
    ##  16800K .......... .......... .......... .......... .......... 21% 75.9M 1s
    ##  16850K .......... .......... .......... .......... .......... 21% 62.4M 1s
    ##  16900K .......... .......... .......... .......... .......... 21% 98.9M 1s
    ##  16950K .......... .......... .......... .......... .......... 21% 60.1M 1s
    ##  17000K .......... .......... .......... .......... .......... 21% 33.8M 1s
    ##  17050K .......... .......... .......... .......... .......... 21% 98.7M 1s
    ##  17100K .......... .......... .......... .......... .......... 21% 39.8M 1s
    ##  17150K .......... .......... .......... .......... .......... 21% 89.5M 1s
    ##  17200K .......... .......... .......... .......... .......... 21% 66.4M 1s
    ##  17250K .......... .......... .......... .......... .......... 21% 43.3M 1s
    ##  17300K .......... .......... .......... .......... .......... 21% 71.6M 1s
    ##  17350K .......... .......... .......... .......... .......... 21%  131M 1s
    ##  17400K .......... .......... .......... .......... .......... 21% 47.5M 1s
    ##  17450K .......... .......... .......... .......... .......... 21% 65.4M 1s
    ##  17500K .......... .......... .......... .......... .......... 21% 49.2M 1s
    ##  17550K .......... .......... .......... .......... .......... 22%  112M 1s
    ##  17600K .......... .......... .......... .......... .......... 22% 58.3M 1s
    ##  17650K .......... .......... .......... .......... .......... 22% 62.9M 1s
    ##  17700K .......... .......... .......... .......... .......... 22% 66.4M 1s
    ##  17750K .......... .......... .......... .......... .......... 22% 48.9M 1s
    ##  17800K .......... .......... .......... .......... .......... 22% 48.8M 1s
    ##  17850K .......... .......... .......... .......... .......... 22% 81.5M 1s
    ##  17900K .......... .......... .......... .......... .......... 22%  100M 1s
    ##  17950K .......... .......... .......... .......... .......... 22% 53.1M 1s
    ##  18000K .......... .......... .......... .......... .......... 22% 46.5M 1s
    ##  18050K .......... .......... .......... .......... .......... 22% 97.1M 1s
    ##  18100K .......... .......... .......... .......... .......... 22% 46.8M 1s
    ##  18150K .......... .......... .......... .......... .......... 22% 93.9M 1s
    ##  18200K .......... .......... .......... .......... .......... 22% 46.9M 1s
    ##  18250K .......... .......... .......... .......... .......... 22%  105M 1s
    ##  18300K .......... .......... .......... .......... .......... 22% 41.6M 1s
    ##  18350K .......... .......... .......... .......... .......... 23% 87.6M 1s
    ##  18400K .......... .......... .......... .......... .......... 23% 58.9M 1s
    ##  18450K .......... .......... .......... .......... .......... 23% 58.0M 1s
    ##  18500K .......... .......... .......... .......... .......... 23% 92.8M 1s
    ##  18550K .......... .......... .......... .......... .......... 23% 58.6M 1s
    ##  18600K .......... .......... .......... .......... .......... 23% 49.0M 1s
    ##  18650K .......... .......... .......... .......... .......... 23% 79.3M 1s
    ##  18700K .......... .......... .......... .......... .......... 23% 45.7M 1s
    ##  18750K .......... .......... .......... .......... .......... 23% 18.2M 1s
    ##  18800K .......... .......... .......... .......... .......... 23% 98.7M 1s
    ##  18850K .......... .......... .......... .......... .......... 23%  102M 1s
    ##  18900K .......... .......... .......... .......... .......... 23% 57.0M 1s
    ##  18950K .......... .......... .......... .......... .......... 23% 39.0M 1s
    ##  19000K .......... .......... .......... .......... .......... 23% 23.3M 1s
    ##  19050K .......... .......... .......... .......... .......... 23% 58.6M 1s
    ##  19100K .......... .......... .......... .......... .......... 23%  102M 1s
    ##  19150K .......... .......... .......... .......... .......... 24%  118M 1s
    ##  19200K .......... .......... .......... .......... .......... 24% 88.2M 1s
    ##  19250K .......... .......... .......... .......... .......... 24% 97.9M 1s
    ##  19300K .......... .......... .......... .......... .......... 24% 95.9M 1s
    ##  19350K .......... .......... .......... .......... .......... 24% 34.9M 1s
    ##  19400K .......... .......... .......... .......... .......... 24%  103M 1s
    ##  19450K .......... .......... .......... .......... .......... 24%  120M 1s
    ##  19500K .......... .......... .......... .......... .......... 24% 90.0M 1s
    ##  19550K .......... .......... .......... .......... .......... 24% 68.3M 1s
    ##  19600K .......... .......... .......... .......... .......... 24% 78.3M 1s
    ##  19650K .......... .......... .......... .......... .......... 24% 91.7M 1s
    ##  19700K .......... .......... .......... .......... .......... 24% 44.5M 1s
    ##  19750K .......... .......... .......... .......... .......... 24%  104M 1s
    ##  19800K .......... .......... .......... .......... .......... 24% 20.6M 1s
    ##  19850K .......... .......... .......... .......... .......... 24% 83.2M 1s
    ##  19900K .......... .......... .......... .......... .......... 24%  107M 1s
    ##  19950K .......... .......... .......... .......... .......... 25%  125M 1s
    ##  20000K .......... .......... .......... .......... .......... 25%  112M 1s
    ##  20050K .......... .......... .......... .......... .......... 25% 37.8M 1s
    ##  20100K .......... .......... .......... .......... .......... 25% 56.4M 1s
    ##  20150K .......... .......... .......... .......... .......... 25% 61.0M 1s
    ##  20200K .......... .......... .......... .......... .......... 25% 78.1M 1s
    ##  20250K .......... .......... .......... .......... .......... 25% 73.4M 1s
    ##  20300K .......... .......... .......... .......... .......... 25% 35.9M 1s
    ##  20350K .......... .......... .......... .......... .......... 25% 98.6M 1s
    ##  20400K .......... .......... .......... .......... .......... 25%  100M 1s
    ##  20450K .......... .......... .......... .......... .......... 25% 81.8M 1s
    ##  20500K .......... .......... .......... .......... .......... 25% 66.3M 1s
    ##  20550K .......... .......... .......... .......... .......... 25% 36.7M 1s
    ##  20600K .......... .......... .......... .......... .......... 25% 86.7M 1s
    ##  20650K .......... .......... .......... .......... .......... 25% 74.0M 1s
    ##  20700K .......... .......... .......... .......... .......... 25% 59.4M 1s
    ##  20750K .......... .......... .......... .......... .......... 26% 81.2M 1s
    ##  20800K .......... .......... .......... .......... .......... 26% 41.8M 1s
    ##  20850K .......... .......... .......... .......... .......... 26% 66.1M 1s
    ##  20900K .......... .......... .......... .......... .......... 26% 86.0M 1s
    ##  20950K .......... .......... .......... .......... .......... 26% 99.2M 1s
    ##  21000K .......... .......... .......... .......... .......... 26%  112M 1s
    ##  21050K .......... .......... .......... .......... .......... 26% 25.6M 1s
    ##  21100K .......... .......... .......... .......... .......... 26% 73.7M 1s
    ##  21150K .......... .......... .......... .......... .......... 26%  112M 1s
    ##  21200K .......... .......... .......... .......... .......... 26% 82.7M 1s
    ##  21250K .......... .......... .......... .......... .......... 26% 96.9M 1s
    ##  21300K .......... .......... .......... .......... .......... 26% 32.9M 1s
    ##  21350K .......... .......... .......... .......... .......... 26% 53.0M 1s
    ##  21400K .......... .......... .......... .......... .......... 26% 29.9M 1s
    ##  21450K .......... .......... .......... .......... .......... 26% 77.2M 1s
    ##  21500K .......... .......... .......... .......... .......... 26% 74.5M 1s
    ##  21550K .......... .......... .......... .......... .......... 27% 67.9M 1s
    ##  21600K .......... .......... .......... .......... .......... 27% 76.4M 1s
    ##  21650K .......... .......... .......... .......... .......... 27% 73.9M 1s
    ##  21700K .......... .......... .......... .......... .......... 27% 6.75M 1s
    ##  21750K .......... .......... .......... .......... .......... 27% 37.3M 1s
    ##  21800K .......... .......... .......... .......... .......... 27% 72.4M 1s
    ##  21850K .......... .......... .......... .......... .......... 27% 82.3M 1s
    ##  21900K .......... .......... .......... .......... .......... 27% 89.7M 1s
    ##  21950K .......... .......... .......... .......... .......... 27%  112M 1s
    ##  22000K .......... .......... .......... .......... .......... 27% 63.8M 1s
    ##  22050K .......... .......... .......... .......... .......... 27% 61.6M 1s
    ##  22100K .......... .......... .......... .......... .......... 27% 63.0M 1s
    ##  22150K .......... .......... .......... .......... .......... 27% 74.2M 1s
    ##  22200K .......... .......... .......... .......... .......... 27%  102M 1s
    ##  22250K .......... .......... .......... .......... .......... 27% 85.3M 1s
    ##  22300K .......... .......... .......... .......... .......... 27% 80.9M 1s
    ##  22350K .......... .......... .......... .......... .......... 28% 91.2M 1s
    ##  22400K .......... .......... .......... .......... .......... 28% 74.5M 1s
    ##  22450K .......... .......... .......... .......... .......... 28% 67.9M 1s
    ##  22500K .......... .......... .......... .......... .......... 28% 39.7M 1s
    ##  22550K .......... .......... .......... .......... .......... 28% 76.0M 1s
    ##  22600K .......... .......... .......... .......... .......... 28% 38.2M 1s
    ##  22650K .......... .......... .......... .......... .......... 28% 65.6M 1s
    ##  22700K .......... .......... .......... .......... .......... 28% 79.8M 1s
    ##  22750K .......... .......... .......... .......... .......... 28% 75.8M 1s
    ##  22800K .......... .......... .......... .......... .......... 28% 89.0M 1s
    ##  22850K .......... .......... .......... .......... .......... 28% 73.5M 1s
    ##  22900K .......... .......... .......... .......... .......... 28% 61.6M 1s
    ##  22950K .......... .......... .......... .......... .......... 28% 66.1M 1s
    ##  23000K .......... .......... .......... .......... .......... 28% 65.3M 1s
    ##  23050K .......... .......... .......... .......... .......... 28% 50.5M 1s
    ##  23100K .......... .......... .......... .......... .......... 28% 78.5M 1s
    ##  23150K .......... .......... .......... .......... .......... 29% 80.1M 1s
    ##  23200K .......... .......... .......... .......... .......... 29% 99.4M 1s
    ##  23250K .......... .......... .......... .......... .......... 29% 4.09M 1s
    ##  23300K .......... .......... .......... .......... .......... 29% 90.8M 1s
    ##  23350K .......... .......... .......... .......... .......... 29%  109M 1s
    ##  23400K .......... .......... .......... .......... .......... 29%  107M 1s
    ##  23450K .......... .......... .......... .......... .......... 29%  111M 1s
    ##  23500K .......... .......... .......... .......... .......... 29% 36.7M 1s
    ##  23550K .......... .......... .......... .......... .......... 29% 56.4M 1s
    ##  23600K .......... .......... .......... .......... .......... 29% 28.5M 1s
    ##  23650K .......... .......... .......... .......... .......... 29% 64.7M 1s
    ##  23700K .......... .......... .......... .......... .......... 29% 71.8M 1s
    ##  23750K .......... .......... .......... .......... .......... 29% 77.9M 1s
    ##  23800K .......... .......... .......... .......... .......... 29% 74.1M 1s
    ##  23850K .......... .......... .......... .......... .......... 29% 25.7M 1s
    ##  23900K .......... .......... .......... .......... .......... 29% 71.3M 1s
    ##  23950K .......... .......... .......... .......... .......... 30% 48.6M 1s
    ##  24000K .......... .......... .......... .......... .......... 30%  106M 1s
    ##  24050K .......... .......... .......... .......... .......... 30% 98.4M 1s
    ##  24100K .......... .......... .......... .......... .......... 30% 85.6M 1s
    ##  24150K .......... .......... .......... .......... .......... 30% 64.4M 1s
    ##  24200K .......... .......... .......... .......... .......... 30% 58.2M 1s
    ##  24250K .......... .......... .......... .......... .......... 30% 54.3M 1s
    ##  24300K .......... .......... .......... .......... .......... 30% 62.2M 1s
    ##  24350K .......... .......... .......... .......... .......... 30% 91.3M 1s
    ##  24400K .......... .......... .......... .......... .......... 30% 65.3M 1s
    ##  24450K .......... .......... .......... .......... .......... 30% 58.4M 1s
    ##  24500K .......... .......... .......... .......... .......... 30% 52.0M 1s
    ##  24550K .......... .......... .......... .......... .......... 30% 64.9M 1s
    ##  24600K .......... .......... .......... .......... .......... 30% 75.7M 1s
    ##  24650K .......... .......... .......... .......... .......... 30%  102M 1s
    ##  24700K .......... .......... .......... .......... .......... 30% 22.7M 1s
    ##  24750K .......... .......... .......... .......... .......... 31% 62.4M 1s
    ##  24800K .......... .......... .......... .......... .......... 31% 52.0M 1s
    ##  24850K .......... .......... .......... .......... .......... 31% 91.1M 1s
    ##  24900K .......... .......... .......... .......... .......... 31% 80.5M 1s
    ##  24950K .......... .......... .......... .......... .......... 31% 55.8M 1s
    ##  25000K .......... .......... .......... .......... .......... 31% 53.4M 1s
    ##  25050K .......... .......... .......... .......... .......... 31% 69.2M 1s
    ##  25100K .......... .......... .......... .......... .......... 31% 41.0M 1s
    ##  25150K .......... .......... .......... .......... .......... 31%  111M 1s
    ##  25200K .......... .......... .......... .......... .......... 31% 86.4M 1s
    ##  25250K .......... .......... .......... .......... .......... 31% 71.0M 1s
    ##  25300K .......... .......... .......... .......... .......... 31% 47.6M 1s
    ##  25350K .......... .......... .......... .......... .......... 31% 77.2M 1s
    ##  25400K .......... .......... .......... .......... .......... 31% 38.0M 1s
    ##  25450K .......... .......... .......... .......... .......... 31% 86.3M 1s
    ##  25500K .......... .......... .......... .......... .......... 31% 79.4M 1s
    ##  25550K .......... .......... .......... .......... .......... 32% 71.5M 1s
    ##  25600K .......... .......... .......... .......... .......... 32% 52.5M 1s
    ##  25650K .......... .......... .......... .......... .......... 32% 82.1M 1s
    ##  25700K .......... .......... .......... .......... .......... 32% 41.5M 1s
    ##  25750K .......... .......... .......... .......... .......... 32% 72.8M 1s
    ##  25800K .......... .......... .......... .......... .......... 32% 49.2M 1s
    ##  25850K .......... .......... .......... .......... .......... 32%  104M 1s
    ##  25900K .......... .......... .......... .......... .......... 32% 58.6M 1s
    ##  25950K .......... .......... .......... .......... .......... 32% 57.8M 1s
    ##  26000K .......... .......... .......... .......... .......... 32% 57.6M 1s
    ##  26050K .......... .......... .......... .......... .......... 32%  100M 1s
    ##  26100K .......... .......... .......... .......... .......... 32% 77.2M 1s
    ##  26150K .......... .......... .......... .......... .......... 32% 31.4M 1s
    ##  26200K .......... .......... .......... .......... .......... 32% 89.8M 1s
    ##  26250K .......... .......... .......... .......... .......... 32%  116M 1s
    ##  26300K .......... .......... .......... .......... .......... 32%  105M 1s
    ##  26350K .......... .......... .......... .......... .......... 33% 77.2M 1s
    ##  26400K .......... .......... .......... .......... .......... 33% 53.0M 1s
    ##  26450K .......... .......... .......... .......... .......... 33% 59.7M 1s
    ##  26500K .......... .......... .......... .......... .......... 33% 48.9M 1s
    ##  26550K .......... .......... .......... .......... .......... 33% 96.1M 1s
    ##  26600K .......... .......... .......... .......... .......... 33% 54.5M 1s
    ##  26650K .......... .......... .......... .......... .......... 33% 61.6M 1s
    ##  26700K .......... .......... .......... .......... .......... 33% 46.5M 1s
    ##  26750K .......... .......... .......... .......... .......... 33% 95.6M 1s
    ##  26800K .......... .......... .......... .......... .......... 33% 72.8M 1s
    ##  26850K .......... .......... .......... .......... .......... 33% 58.5M 1s
    ##  26900K .......... .......... .......... .......... .......... 33% 17.5M 1s
    ##  26950K .......... .......... .......... .......... .......... 33%  115M 1s
    ##  27000K .......... .......... .......... .......... .......... 33%  118M 1s
    ##  27050K .......... .......... .......... .......... .......... 33%  123M 1s
    ##  27100K .......... .......... .......... .......... .......... 33%  125M 1s
    ##  27150K .......... .......... .......... .......... .......... 34%  108M 1s
    ##  27200K .......... .......... .......... .......... .......... 34% 86.9M 1s
    ##  27250K .......... .......... .......... .......... .......... 34% 32.9M 1s
    ##  27300K .......... .......... .......... .......... .......... 34% 42.8M 1s
    ##  27350K .......... .......... .......... .......... .......... 34% 48.1M 1s
    ##  27400K .......... .......... .......... .......... .......... 34% 55.0M 1s
    ##  27450K .......... .......... .......... .......... .......... 34% 44.9M 1s
    ##  27500K .......... .......... .......... .......... .......... 34% 98.2M 1s
    ##  27550K .......... .......... .......... .......... .......... 34% 40.7M 1s
    ##  27600K .......... .......... .......... .......... .......... 34% 85.0M 1s
    ##  27650K .......... .......... .......... .......... .......... 34% 98.3M 1s
    ##  27700K .......... .......... .......... .......... .......... 34%  112M 1s
    ##  27750K .......... .......... .......... .......... .......... 34% 44.8M 1s
    ##  27800K .......... .......... .......... .......... .......... 34% 31.5M 1s
    ##  27850K .......... .......... .......... .......... .......... 34% 96.6M 1s
    ##  27900K .......... .......... .......... .......... .......... 34% 61.2M 1s
    ##  27950K .......... .......... .......... .......... .......... 35% 83.6M 1s
    ##  28000K .......... .......... .......... .......... .......... 35% 74.1M 1s
    ##  28050K .......... .......... .......... .......... .......... 35% 97.8M 1s
    ##  28100K .......... .......... .......... .......... .......... 35% 23.5M 1s
    ##  28150K .......... .......... .......... .......... .......... 35% 20.1M 1s
    ##  28200K .......... .......... .......... .......... .......... 35%  103M 1s
    ##  28250K .......... .......... .......... .......... .......... 35% 87.0M 1s
    ##  28300K .......... .......... .......... .......... .......... 35%  120M 1s
    ##  28350K .......... .......... .......... .......... .......... 35% 94.7M 1s
    ##  28400K .......... .......... .......... .......... .......... 35% 58.8M 1s
    ##  28450K .......... .......... .......... .......... .......... 35% 57.6M 1s
    ##  28500K .......... .......... .......... .......... .......... 35% 53.7M 1s
    ##  28550K .......... .......... .......... .......... .......... 35% 55.2M 1s
    ##  28600K .......... .......... .......... .......... .......... 35% 62.7M 1s
    ##  28650K .......... .......... .......... .......... .......... 35%  100M 1s
    ##  28700K .......... .......... .......... .......... .......... 35%  124M 1s
    ##  28750K .......... .......... .......... .......... .......... 36% 41.7M 1s
    ##  28800K .......... .......... .......... .......... .......... 36% 48.0M 1s
    ##  28850K .......... .......... .......... .......... .......... 36% 57.3M 1s
    ##  28900K .......... .......... .......... .......... .......... 36% 50.2M 1s
    ##  28950K .......... .......... .......... .......... .......... 36%  114M 1s
    ##  29000K .......... .......... .......... .......... .......... 36%  103M 1s
    ##  29050K .......... .......... .......... .......... .......... 36% 43.7M 1s
    ##  29100K .......... .......... .......... .......... .......... 36%  110M 1s
    ##  29150K .......... .......... .......... .......... .......... 36% 63.7M 1s
    ##  29200K .......... .......... .......... .......... .......... 36% 43.3M 1s
    ##  29250K .......... .......... .......... .......... .......... 36%  103M 1s
    ##  29300K .......... .......... .......... .......... .......... 36% 45.7M 1s
    ##  29350K .......... .......... .......... .......... .......... 36%  118M 1s
    ##  29400K .......... .......... .......... .......... .......... 36%  101M 1s
    ##  29450K .......... .......... .......... .......... .......... 36% 77.8M 1s
    ##  29500K .......... .......... .......... .......... .......... 36% 52.8M 1s
    ##  29550K .......... .......... .......... .......... .......... 37% 75.9M 1s
    ##  29600K .......... .......... .......... .......... .......... 37% 67.7M 1s
    ##  29650K .......... .......... .......... .......... .......... 37% 69.6M 1s
    ##  29700K .......... .......... .......... .......... .......... 37% 47.9M 1s
    ##  29750K .......... .......... .......... .......... .......... 37% 99.3M 1s
    ##  29800K .......... .......... .......... .......... .......... 37% 49.0M 1s
    ##  29850K .......... .......... .......... .......... .......... 37%  120M 1s
    ##  29900K .......... .......... .......... .......... .......... 37% 45.8M 1s
    ##  29950K .......... .......... .......... .......... .......... 37% 49.6M 1s
    ##  30000K .......... .......... .......... .......... .......... 37% 81.8M 1s
    ##  30050K .......... .......... .......... .......... .......... 37% 51.7M 1s
    ##  30100K .......... .......... .......... .......... .......... 37% 79.0M 1s
    ##  30150K .......... .......... .......... .......... .......... 37% 62.3M 1s
    ##  30200K .......... .......... .......... .......... .......... 37% 40.4M 1s
    ##  30250K .......... .......... .......... .......... .......... 37%  114M 1s
    ##  30300K .......... .......... .......... .......... .......... 37% 91.6M 1s
    ##  30350K .......... .......... .......... .......... .......... 38% 46.1M 1s
    ##  30400K .......... .......... .......... .......... .......... 38% 57.3M 1s
    ##  30450K .......... .......... .......... .......... .......... 38% 62.6M 1s
    ##  30500K .......... .......... .......... .......... .......... 38% 90.9M 1s
    ##  30550K .......... .......... .......... .......... .......... 38% 97.7M 1s
    ##  30600K .......... .......... .......... .......... .......... 38% 72.5M 1s
    ##  30650K .......... .......... .......... .......... .......... 38% 33.2M 1s
    ##  30700K .......... .......... .......... .......... .......... 38% 96.3M 1s
    ##  30750K .......... .......... .......... .......... .......... 38%  125M 1s
    ##  30800K .......... .......... .......... .......... .......... 38% 30.7M 1s
    ##  30850K .......... .......... .......... .......... .......... 38% 77.0M 1s
    ##  30900K .......... .......... .......... .......... .......... 38%  107M 1s
    ##  30950K .......... .......... .......... .......... .......... 38% 83.7M 1s
    ##  31000K .......... .......... .......... .......... .......... 38%  106M 1s
    ##  31050K .......... .......... .......... .......... .......... 38% 65.8M 1s
    ##  31100K .......... .......... .......... .......... .......... 38% 33.2M 1s
    ##  31150K .......... .......... .......... .......... .......... 39% 91.2M 1s
    ##  31200K .......... .......... .......... .......... .......... 39% 87.1M 1s
    ##  31250K .......... .......... .......... .......... .......... 39%  113M 1s
    ##  31300K .......... .......... .......... .......... .......... 39% 85.0M 1s
    ##  31350K .......... .......... .......... .......... .......... 39% 64.2M 1s
    ##  31400K .......... .......... .......... .......... .......... 39% 79.7M 1s
    ##  31450K .......... .......... .......... .......... .......... 39% 48.0M 1s
    ##  31500K .......... .......... .......... .......... .......... 39% 83.9M 1s
    ##  31550K .......... .......... .......... .......... .......... 39% 74.6M 1s
    ##  31600K .......... .......... .......... .......... .......... 39% 56.6M 1s
    ##  31650K .......... .......... .......... .......... .......... 39% 74.8M 1s
    ##  31700K .......... .......... .......... .......... .......... 39%  112M 1s
    ##  31750K .......... .......... .......... .......... .......... 39% 69.0M 1s
    ##  31800K .......... .......... .......... .......... .......... 39% 60.5M 1s
    ##  31850K .......... .......... .......... .......... .......... 39% 54.0M 1s
    ##  31900K .......... .......... .......... .......... .......... 39% 68.1M 1s
    ##  31950K .......... .......... .......... .......... .......... 40% 89.9M 1s
    ##  32000K .......... .......... .......... .......... .......... 40% 87.0M 1s
    ##  32050K .......... .......... .......... .......... .......... 40% 73.1M 1s
    ##  32100K .......... .......... .......... .......... .......... 40% 28.8M 1s
    ##  32150K .......... .......... .......... .......... .......... 40%  116M 1s
    ##  32200K .......... .......... .......... .......... .......... 40% 99.4M 1s
    ##  32250K .......... .......... .......... .......... .......... 40%  120M 1s
    ##  32300K .......... .......... .......... .......... .......... 40% 67.4M 1s
    ##  32350K .......... .......... .......... .......... .......... 40%  103M 1s
    ##  32400K .......... .......... .......... .......... .......... 40%  110M 1s
    ##  32450K .......... .......... .......... .......... .......... 40% 32.1M 1s
    ##  32500K .......... .......... .......... .......... .......... 40% 83.1M 1s
    ##  32550K .......... .......... .......... .......... .......... 40% 61.1M 1s
    ##  32600K .......... .......... .......... .......... .......... 40% 36.6M 1s
    ##  32650K .......... .......... .......... .......... .......... 40% 86.4M 1s
    ##  32700K .......... .......... .......... .......... .......... 40% 91.1M 1s
    ##  32750K .......... .......... .......... .......... .......... 41%  101M 1s
    ##  32800K .......... .......... .......... .......... .......... 41% 53.2M 1s
    ##  32850K .......... .......... .......... .......... .......... 41%  106M 1s
    ##  32900K .......... .......... .......... .......... .......... 41% 37.5M 1s
    ##  32950K .......... .......... .......... .......... .......... 41% 64.6M 1s
    ##  33000K .......... .......... .......... .......... .......... 41% 41.1M 1s
    ##  33050K .......... .......... .......... .......... .......... 41% 76.9M 1s
    ##  33100K .......... .......... .......... .......... .......... 41% 59.2M 1s
    ##  33150K .......... .......... .......... .......... .......... 41% 89.9M 1s
    ##  33200K .......... .......... .......... .......... .......... 41% 78.8M 1s
    ##  33250K .......... .......... .......... .......... .......... 41% 56.0M 1s
    ##  33300K .......... .......... .......... .......... .......... 41% 88.9M 1s
    ##  33350K .......... .......... .......... .......... .......... 41% 44.9M 1s
    ##  33400K .......... .......... .......... .......... .......... 41% 73.9M 1s
    ##  33450K .......... .......... .......... .......... .......... 41% 99.4M 1s
    ##  33500K .......... .......... .......... .......... .......... 41% 53.0M 1s
    ##  33550K .......... .......... .......... .......... .......... 42% 88.4M 1s
    ##  33600K .......... .......... .......... .......... .......... 42% 86.9M 1s
    ##  33650K .......... .......... .......... .......... .......... 42% 50.8M 1s
    ##  33700K .......... .......... .......... .......... .......... 42% 78.0M 1s
    ##  33750K .......... .......... .......... .......... .......... 42% 62.9M 1s
    ##  33800K .......... .......... .......... .......... .......... 42% 53.4M 1s
    ##  33850K .......... .......... .......... .......... .......... 42% 78.9M 1s
    ##  33900K .......... .......... .......... .......... .......... 42% 70.1M 1s
    ##  33950K .......... .......... .......... .......... .......... 42% 72.9M 1s
    ##  34000K .......... .......... .......... .......... .......... 42% 63.0M 1s
    ##  34050K .......... .......... .......... .......... .......... 42% 70.9M 1s
    ##  34100K .......... .......... .......... .......... .......... 42% 89.0M 1s
    ##  34150K .......... .......... .......... .......... .......... 42% 51.0M 1s
    ##  34200K .......... .......... .......... .......... .......... 42% 63.0M 1s
    ##  34250K .......... .......... .......... .......... .......... 42% 74.3M 1s
    ##  34300K .......... .......... .......... .......... .......... 42% 51.3M 1s
    ##  34350K .......... .......... .......... .......... .......... 43%  109M 1s
    ##  34400K .......... .......... .......... .......... .......... 43% 55.7M 1s
    ##  34450K .......... .......... .......... .......... .......... 43% 68.5M 1s
    ##  34500K .......... .......... .......... .......... .......... 43% 69.4M 1s
    ##  34550K .......... .......... .......... .......... .......... 43% 99.2M 1s
    ##  34600K .......... .......... .......... .......... .......... 43% 66.3M 1s
    ##  34650K .......... .......... .......... .......... .......... 43%  103M 1s
    ##  34700K .......... .......... .......... .......... .......... 43% 51.2M 1s
    ##  34750K .......... .......... .......... .......... .......... 43%  139M 1s
    ##  34800K .......... .......... .......... .......... .......... 43% 59.5M 1s
    ##  34850K .......... .......... .......... .......... .......... 43% 75.2M 1s
    ##  34900K .......... .......... .......... .......... .......... 43% 67.8M 1s
    ##  34950K .......... .......... .......... .......... .......... 43% 56.0M 1s
    ##  35000K .......... .......... .......... .......... .......... 43% 68.2M 1s
    ##  35050K .......... .......... .......... .......... .......... 43%  108M 1s
    ##  35100K .......... .......... .......... .......... .......... 43% 60.3M 1s
    ##  35150K .......... .......... .......... .......... .......... 44% 73.3M 1s
    ##  35200K .......... .......... .......... .......... .......... 44% 68.1M 1s
    ##  35250K .......... .......... .......... .......... .......... 44% 98.2M 1s
    ##  35300K .......... .......... .......... .......... .......... 44% 71.7M 1s
    ##  35350K .......... .......... .......... .......... .......... 44% 72.1M 1s
    ##  35400K .......... .......... .......... .......... .......... 44% 75.4M 1s
    ##  35450K .......... .......... .......... .......... .......... 44% 74.7M 1s
    ##  35500K .......... .......... .......... .......... .......... 44% 91.8M 1s
    ##  35550K .......... .......... .......... .......... .......... 44% 68.2M 1s
    ##  35600K .......... .......... .......... .......... .......... 44% 64.2M 1s
    ##  35650K .......... .......... .......... .......... .......... 44% 69.2M 1s
    ##  35700K .......... .......... .......... .......... .......... 44%  115M 1s
    ##  35750K .......... .......... .......... .......... .......... 44%  110M 1s
    ##  35800K .......... .......... .......... .......... .......... 44% 72.6M 1s
    ##  35850K .......... .......... .......... .......... .......... 44% 78.1M 1s
    ##  35900K .......... .......... .......... .......... .......... 44% 59.6M 1s
    ##  35950K .......... .......... .......... .......... .......... 45%  126M 1s
    ##  36000K .......... .......... .......... .......... .......... 45% 87.5M 1s
    ##  36050K .......... .......... .......... .......... .......... 45% 62.5M 1s
    ##  36100K .......... .......... .......... .......... .......... 45% 73.9M 1s
    ##  36150K .......... .......... .......... .......... .......... 45%  132M 1s
    ##  36200K .......... .......... .......... .......... .......... 45% 81.9M 1s
    ##  36250K .......... .......... .......... .......... .......... 45% 79.8M 1s
    ##  36300K .......... .......... .......... .......... .......... 45% 82.0M 1s
    ##  36350K .......... .......... .......... .......... .......... 45% 61.1M 1s
    ##  36400K .......... .......... .......... .......... .......... 45% 75.0M 1s
    ##  36450K .......... .......... .......... .......... .......... 45%  104M 1s
    ##  36500K .......... .......... .......... .......... .......... 45% 70.6M 1s
    ##  36550K .......... .......... .......... .......... .......... 45%  101M 1s
    ##  36600K .......... .......... .......... .......... .......... 45% 82.4M 1s
    ##  36650K .......... .......... .......... .......... .......... 45% 67.1M 1s
    ##  36700K .......... .......... .......... .......... .......... 45% 97.9M 1s
    ##  36750K .......... .......... .......... .......... .......... 46% 80.3M 1s
    ##  36800K .......... .......... .......... .......... .......... 46% 85.4M 1s
    ##  36850K .......... .......... .......... .......... .......... 46% 82.0M 1s
    ##  36900K .......... .......... .......... .......... .......... 46%  106M 1s
    ##  36950K .......... .......... .......... .......... .......... 46% 62.6M 1s
    ##  37000K .......... .......... .......... .......... .......... 46% 96.7M 1s
    ##  37050K .......... .......... .......... .......... .......... 46% 60.1M 1s
    ##  37100K .......... .......... .......... .......... .......... 46%  100M 1s
    ##  37150K .......... .......... .......... .......... .......... 46% 82.7M 1s
    ##  37200K .......... .......... .......... .......... .......... 46% 90.8M 1s
    ##  37250K .......... .......... .......... .......... .......... 46% 76.6M 1s
    ##  37300K .......... .......... .......... .......... .......... 46%  112M 1s
    ##  37350K .......... .......... .......... .......... .......... 46% 22.3M 1s
    ##  37400K .......... .......... .......... .......... .......... 46% 9.94M 1s
    ##  37450K .......... .......... .......... .......... .......... 46%  118M 1s
    ##  37500K .......... .......... .......... .......... .......... 46%  142M 1s
    ##  37550K .......... .......... .......... .......... .......... 47%  145M 1s
    ##  37600K .......... .......... .......... .......... .......... 47% 65.6M 1s
    ##  37650K .......... .......... .......... .......... .......... 47%  124M 1s
    ##  37700K .......... .......... .......... .......... .......... 47% 45.4M 1s
    ##  37750K .......... .......... .......... .......... .......... 47% 53.3M 1s
    ##  37800K .......... .......... .......... .......... .......... 47% 29.7M 1s
    ##  37850K .......... .......... .......... .......... .......... 47% 87.4M 1s
    ##  37900K .......... .......... .......... .......... .......... 47% 96.1M 1s
    ##  37950K .......... .......... .......... .......... .......... 47% 76.4M 1s
    ##  38000K .......... .......... .......... .......... .......... 47% 36.4M 1s
    ##  38050K .......... .......... .......... .......... .......... 47%  148M 1s
    ##  38100K .......... .......... .......... .......... .......... 47%  109M 1s
    ##  38150K .......... .......... .......... .......... .......... 47% 98.7M 1s
    ##  38200K .......... .......... .......... .......... .......... 47% 33.7M 1s
    ##  38250K .......... .......... .......... .......... .......... 47% 86.1M 1s
    ##  38300K .......... .......... .......... .......... .......... 47% 97.4M 1s
    ##  38350K .......... .......... .......... .......... .......... 48% 67.3M 1s
    ##  38400K .......... .......... .......... .......... .......... 48% 49.8M 1s
    ##  38450K .......... .......... .......... .......... .......... 48%  100M 1s
    ##  38500K .......... .......... .......... .......... .......... 48%  118M 1s
    ##  38550K .......... .......... .......... .......... .......... 48% 69.4M 1s
    ##  38600K .......... .......... .......... .......... .......... 48%  109M 1s
    ##  38650K .......... .......... .......... .......... .......... 48% 53.4M 1s
    ##  38700K .......... .......... .......... .......... .......... 48% 59.0M 1s
    ##  38750K .......... .......... .......... .......... .......... 48% 90.1M 1s
    ##  38800K .......... .......... .......... .......... .......... 48% 98.0M 1s
    ##  38850K .......... .......... .......... .......... .......... 48% 59.6M 1s
    ##  38900K .......... .......... .......... .......... .......... 48% 54.3M 1s
    ##  38950K .......... .......... .......... .......... .......... 48% 19.5M 1s
    ##  39000K .......... .......... .......... .......... .......... 48%  116M 1s
    ##  39050K .......... .......... .......... .......... .......... 48%  168M 1s
    ##  39100K .......... .......... .......... .......... .......... 48% 49.2M 1s
    ##  39150K .......... .......... .......... .......... .......... 49%  125M 1s
    ##  39200K .......... .......... .......... .......... .......... 49%  138M 1s
    ##  39250K .......... .......... .......... .......... .......... 49% 17.1M 1s
    ##  39300K .......... .......... .......... .......... .......... 49% 17.8M 1s
    ##  39350K .......... .......... .......... .......... .......... 49% 59.2M 1s
    ##  39400K .......... .......... .......... .......... .......... 49% 77.9M 1s
    ##  39450K .......... .......... .......... .......... .......... 49% 76.8M 1s
    ##  39500K .......... .......... .......... .......... .......... 49% 65.5M 1s
    ##  39550K .......... .......... .......... .......... .......... 49% 80.3M 1s
    ##  39600K .......... .......... .......... .......... .......... 49% 94.0M 1s
    ##  39650K .......... .......... .......... .......... .......... 49%  103M 1s
    ##  39700K .......... .......... .......... .......... .......... 49% 83.6M 1s
    ##  39750K .......... .......... .......... .......... .......... 49% 63.3M 1s
    ##  39800K .......... .......... .......... .......... .......... 49% 79.4M 1s
    ##  39850K .......... .......... .......... .......... .......... 49% 73.5M 1s
    ##  39900K .......... .......... .......... .......... .......... 49% 49.3M 1s
    ##  39950K .......... .......... .......... .......... .......... 50% 60.4M 1s
    ##  40000K .......... .......... .......... .......... .......... 50% 54.8M 1s
    ##  40050K .......... .......... .......... .......... .......... 50% 71.1M 1s
    ##  40100K .......... .......... .......... .......... .......... 50% 83.0M 1s
    ##  40150K .......... .......... .......... .......... .......... 50% 84.0M 1s
    ##  40200K .......... .......... .......... .......... .......... 50% 71.1M 1s
    ##  40250K .......... .......... .......... .......... .......... 50% 60.0M 1s
    ##  40300K .......... .......... .......... .......... .......... 50% 50.6M 1s
    ##  40350K .......... .......... .......... .......... .......... 50% 69.9M 1s
    ##  40400K .......... .......... .......... .......... .......... 50%  103M 1s
    ##  40450K .......... .......... .......... .......... .......... 50% 82.0M 1s
    ##  40500K .......... .......... .......... .......... .......... 50% 61.4M 1s
    ##  40550K .......... .......... .......... .......... .......... 50% 70.4M 1s
    ##  40600K .......... .......... .......... .......... .......... 50% 68.6M 1s
    ##  40650K .......... .......... .......... .......... .......... 50% 57.2M 1s
    ##  40700K .......... .......... .......... .......... .......... 50% 82.8M 1s
    ##  40750K .......... .......... .......... .......... .......... 51% 92.7M 1s
    ##  40800K .......... .......... .......... .......... .......... 51% 52.7M 1s
    ##  40850K .......... .......... .......... .......... .......... 51% 63.5M 1s
    ##  40900K .......... .......... .......... .......... .......... 51% 53.8M 1s
    ##  40950K .......... .......... .......... .......... .......... 51% 72.4M 1s
    ##  41000K .......... .......... .......... .......... .......... 51% 40.7M 1s
    ##  41050K .......... .......... .......... .......... .......... 51% 62.0M 1s
    ##  41100K .......... .......... .......... .......... .......... 51% 79.9M 1s
    ##  41150K .......... .......... .......... .......... .......... 51% 83.2M 1s
    ##  41200K .......... .......... .......... .......... .......... 51% 65.4M 1s
    ##  41250K .......... .......... .......... .......... .......... 51% 75.3M 1s
    ##  41300K .......... .......... .......... .......... .......... 51% 45.1M 1s
    ##  41350K .......... .......... .......... .......... .......... 51% 78.4M 1s
    ##  41400K .......... .......... .......... .......... .......... 51% 72.6M 1s
    ##  41450K .......... .......... .......... .......... .......... 51% 70.6M 1s
    ##  41500K .......... .......... .......... .......... .......... 51% 63.4M 1s
    ##  41550K .......... .......... .......... .......... .......... 52% 78.9M 1s
    ##  41600K .......... .......... .......... .......... .......... 52% 48.0M 1s
    ##  41650K .......... .......... .......... .......... .......... 52% 81.7M 1s
    ##  41700K .......... .......... .......... .......... .......... 52% 62.4M 1s
    ##  41750K .......... .......... .......... .......... .......... 52% 93.0M 1s
    ##  41800K .......... .......... .......... .......... .......... 52% 48.3M 1s
    ##  41850K .......... .......... .......... .......... .......... 52% 43.1M 1s
    ##  41900K .......... .......... .......... .......... .......... 52% 67.9M 1s
    ##  41950K .......... .......... .......... .......... .......... 52% 71.6M 1s
    ##  42000K .......... .......... .......... .......... .......... 52% 66.8M 1s
    ##  42050K .......... .......... .......... .......... .......... 52% 94.7M 1s
    ##  42100K .......... .......... .......... .......... .......... 52% 83.0M 1s
    ##  42150K .......... .......... .......... .......... .......... 52% 54.9M 1s
    ##  42200K .......... .......... .......... .......... .......... 52% 63.1M 1s
    ##  42250K .......... .......... .......... .......... .......... 52% 60.7M 1s
    ##  42300K .......... .......... .......... .......... .......... 52% 59.3M 1s
    ##  42350K .......... .......... .......... .......... .......... 53% 87.6M 1s
    ##  42400K .......... .......... .......... .......... .......... 53% 49.1M 1s
    ##  42450K .......... .......... .......... .......... .......... 53% 62.5M 1s
    ##  42500K .......... .......... .......... .......... .......... 53% 66.7M 1s
    ##  42550K .......... .......... .......... .......... .......... 53% 78.5M 1s
    ##  42600K .......... .......... .......... .......... .......... 53% 83.1M 1s
    ##  42650K .......... .......... .......... .......... .......... 53% 97.0M 1s
    ##  42700K .......... .......... .......... .......... .......... 53% 55.9M 1s
    ##  42750K .......... .......... .......... .......... .......... 53% 63.7M 1s
    ##  42800K .......... .......... .......... .......... .......... 53% 74.0M 1s
    ##  42850K .......... .......... .......... .......... .......... 53% 90.5M 1s
    ##  42900K .......... .......... .......... .......... .......... 53% 82.5M 1s
    ##  42950K .......... .......... .......... .......... .......... 53% 64.7M 1s
    ##  43000K .......... .......... .......... .......... .......... 53% 74.2M 1s
    ##  43050K .......... .......... .......... .......... .......... 53% 83.2M 1s
    ##  43100K .......... .......... .......... .......... .......... 53% 56.0M 1s
    ##  43150K .......... .......... .......... .......... .......... 54% 85.6M 1s
    ##  43200K .......... .......... .......... .......... .......... 54% 74.1M 1s
    ##  43250K .......... .......... .......... .......... .......... 54% 94.5M 1s
    ##  43300K .......... .......... .......... .......... .......... 54% 59.2M 1s
    ##  43350K .......... .......... .......... .......... .......... 54% 81.3M 1s
    ##  43400K .......... .......... .......... .......... .......... 54% 55.5M 1s
    ##  43450K .......... .......... .......... .......... .......... 54% 91.6M 1s
    ##  43500K .......... .......... .......... .......... .......... 54% 89.3M 1s
    ##  43550K .......... .......... .......... .......... .......... 54% 77.5M 1s
    ##  43600K .......... .......... .......... .......... .......... 54% 55.4M 1s
    ##  43650K .......... .......... .......... .......... .......... 54% 86.2M 1s
    ##  43700K .......... .......... .......... .......... .......... 54% 62.5M 1s
    ##  43750K .......... .......... .......... .......... .......... 54% 93.1M 1s
    ##  43800K .......... .......... .......... .......... .......... 54% 75.6M 1s
    ##  43850K .......... .......... .......... .......... .......... 54% 91.1M 1s
    ##  43900K .......... .......... .......... .......... .......... 54% 66.4M 1s
    ##  43950K .......... .......... .......... .......... .......... 55% 58.3M 1s
    ##  44000K .......... .......... .......... .......... .......... 55% 98.2M 1s
    ##  44050K .......... .......... .......... .......... .......... 55% 90.0M 1s
    ##  44100K .......... .......... .......... .......... .......... 55% 87.5M 1s
    ##  44150K .......... .......... .......... .......... .......... 55% 82.5M 1s
    ##  44200K .......... .......... .......... .......... .......... 55% 85.1M 1s
    ##  44250K .......... .......... .......... .......... .......... 55% 74.7M 1s
    ##  44300K .......... .......... .......... .......... .......... 55% 69.5M 1s
    ##  44350K .......... .......... .......... .......... .......... 55% 75.1M 1s
    ##  44400K .......... .......... .......... .......... .......... 55% 90.5M 1s
    ##  44450K .......... .......... .......... .......... .......... 55% 84.4M 1s
    ##  44500K .......... .......... .......... .......... .......... 55% 69.6M 1s
    ##  44550K .......... .......... .......... .......... .......... 55% 83.3M 1s
    ##  44600K .......... .......... .......... .......... .......... 55% 81.5M 1s
    ##  44650K .......... .......... .......... .......... .......... 55% 74.9M 1s
    ##  44700K .......... .......... .......... .......... .......... 55% 88.8M 1s
    ##  44750K .......... .......... .......... .......... .......... 56% 74.4M 1s
    ##  44800K .......... .......... .......... .......... .......... 56% 87.4M 1s
    ##  44850K .......... .......... .......... .......... .......... 56%  111M 1s
    ##  44900K .......... .......... .......... .......... .......... 56% 72.0M 1s
    ##  44950K .......... .......... .......... .......... .......... 56%  103M 1s
    ##  45000K .......... .......... .......... .......... .......... 56% 79.6M 1s
    ##  45050K .......... .......... .......... .......... .......... 56% 83.9M 1s
    ##  45100K .......... .......... .......... .......... .......... 56%  119M 1s
    ##  45150K .......... .......... .......... .......... .......... 56% 22.8M 1s
    ##  45200K .......... .......... .......... .......... .......... 56% 46.6M 1s
    ##  45250K .......... .......... .......... .......... .......... 56%  103M 1s
    ##  45300K .......... .......... .......... .......... .......... 56% 65.6M 1s
    ##  45350K .......... .......... .......... .......... .......... 56% 91.4M 1s
    ##  45400K .......... .......... .......... .......... .......... 56% 88.4M 1s
    ##  45450K .......... .......... .......... .......... .......... 56% 72.4M 1s
    ##  45500K .......... .......... .......... .......... .......... 56%  101M 1s
    ##  45550K .......... .......... .......... .......... .......... 57% 22.6M 1s
    ##  45600K .......... .......... .......... .......... .......... 57% 30.3M 1s
    ##  45650K .......... .......... .......... .......... .......... 57%  120M 1s
    ##  45700K .......... .......... .......... .......... .......... 57%  118M 1s
    ##  45750K .......... .......... .......... .......... .......... 57%  116M 1s
    ##  45800K .......... .......... .......... .......... .......... 57%  121M 1s
    ##  45850K .......... .......... .......... .......... .......... 57%  112M 1s
    ##  45900K .......... .......... .......... .......... .......... 57% 23.0M 1s
    ##  45950K .......... .......... .......... .......... .......... 57% 54.0M 1s
    ##  46000K .......... .......... .......... .......... .......... 57% 82.1M 1s
    ##  46050K .......... .......... .......... .......... .......... 57% 86.3M 1s
    ##  46100K .......... .......... .......... .......... .......... 57%  119M 1s
    ##  46150K .......... .......... .......... .......... .......... 57%  134M 1s
    ##  46200K .......... .......... .......... .......... .......... 57% 84.5M 1s
    ##  46250K .......... .......... .......... .......... .......... 57% 38.0M 1s
    ##  46300K .......... .......... .......... .......... .......... 57% 22.3M 1s
    ##  46350K .......... .......... .......... .......... .......... 58%  105M 1s
    ##  46400K .......... .......... .......... .......... .......... 58%  119M 1s
    ##  46450K .......... .......... .......... .......... .......... 58%  149M 1s
    ##  46500K .......... .......... .......... .......... .......... 58%  111M 1s
    ##  46550K .......... .......... .......... .......... .......... 58% 99.3M 1s
    ##  46600K .......... .......... .......... .......... .......... 58%  103M 1s
    ##  46650K .......... .......... .......... .......... .......... 58% 49.0M 1s
    ##  46700K .......... .......... .......... .......... .......... 58% 67.5M 1s
    ##  46750K .......... .......... .......... .......... .......... 58%  109M 1s
    ##  46800K .......... .......... .......... .......... .......... 58% 33.7M 1s
    ##  46850K .......... .......... .......... .......... .......... 58%  107M 1s
    ##  46900K .......... .......... .......... .......... .......... 58%  106M 1s
    ##  46950K .......... .......... .......... .......... .......... 58% 19.2M 1s
    ##  47000K .......... .......... .......... .......... .......... 58%  109M 1s
    ##  47050K .......... .......... .......... .......... .......... 58%  108M 1s
    ##  47100K .......... .......... .......... .......... .......... 58%  111M 1s
    ##  47150K .......... .......... .......... .......... .......... 59% 97.3M 1s
    ##  47200K .......... .......... .......... .......... .......... 59% 64.2M 1s
    ##  47250K .......... .......... .......... .......... .......... 59%  116M 1s
    ##  47300K .......... .......... .......... .......... .......... 59% 41.1M 1s
    ##  47350K .......... .......... .......... .......... .......... 59%  113M 1s
    ##  47400K .......... .......... .......... .......... .......... 59% 27.4M 1s
    ##  47450K .......... .......... .......... .......... .......... 59% 97.7M 1s
    ##  47500K .......... .......... .......... .......... .......... 59% 31.3M 1s
    ##  47550K .......... .......... .......... .......... .......... 59%  116M 1s
    ##  47600K .......... .......... .......... .......... .......... 59%  107M 1s
    ##  47650K .......... .......... .......... .......... .......... 59% 78.8M 1s
    ##  47700K .......... .......... .......... .......... .......... 59%  104M 1s
    ##  47750K .......... .......... .......... .......... .......... 59%  145M 1s
    ##  47800K .......... .......... .......... .......... .......... 59% 30.8M 1s
    ##  47850K .......... .......... .......... .......... .......... 59% 89.7M 1s
    ##  47900K .......... .......... .......... .......... .......... 59% 83.7M 1s
    ##  47950K .......... .......... .......... .......... .......... 60%  148M 1s
    ##  48000K .......... .......... .......... .......... .......... 60% 92.3M 1s
    ##  48050K .......... .......... .......... .......... .......... 60%  104M 1s
    ##  48100K .......... .......... .......... .......... .......... 60% 93.4M 1s
    ##  48150K .......... .......... .......... .......... .......... 60% 12.9M 1s
    ##  48200K .......... .......... .......... .......... .......... 60% 48.2M 1s
    ##  48250K .......... .......... .......... .......... .......... 60% 25.0M 1s
    ##  48300K .......... .......... .......... .......... .......... 60% 60.3M 1s
    ##  48350K .......... .......... .......... .......... .......... 60% 64.0M 1s
    ##  48400K .......... .......... .......... .......... .......... 60% 56.5M 1s
    ##  48450K .......... .......... .......... .......... .......... 60% 70.6M 1s
    ##  48500K .......... .......... .......... .......... .......... 60% 58.0M 1s
    ##  48550K .......... .......... .......... .......... .......... 60% 83.3M 1s
    ##  48600K .......... .......... .......... .......... .......... 60% 69.6M 1s
    ##  48650K .......... .......... .......... .......... .......... 60% 34.6M 1s
    ##  48700K .......... .......... .......... .......... .......... 60% 36.7M 1s
    ##  48750K .......... .......... .......... .......... .......... 61% 33.4M 1s
    ##  48800K .......... .......... .......... .......... .......... 61% 36.1M 1s
    ##  48850K .......... .......... .......... .......... .......... 61% 83.3M 1s
    ##  48900K .......... .......... .......... .......... .......... 61% 67.9M 1s
    ##  48950K .......... .......... .......... .......... .......... 61% 48.5M 1s
    ##  49000K .......... .......... .......... .......... .......... 61% 44.4M 1s
    ##  49050K .......... .......... .......... .......... .......... 61% 97.5M 1s
    ##  49100K .......... .......... .......... .......... .......... 61% 75.6M 1s
    ##  49150K .......... .......... .......... .......... .......... 61% 68.9M 1s
    ##  49200K .......... .......... .......... .......... .......... 61% 75.4M 1s
    ##  49250K .......... .......... .......... .......... .......... 61% 88.8M 1s
    ##  49300K .......... .......... .......... .......... .......... 61% 76.2M 1s
    ##  49350K .......... .......... .......... .......... .......... 61% 89.6M 1s
    ##  49400K .......... .......... .......... .......... .......... 61% 77.6M 1s
    ##  49450K .......... .......... .......... .......... .......... 61%  105M 1s
    ##  49500K .......... .......... .......... .......... .......... 61% 49.1M 1s
    ##  49550K .......... .......... .......... .......... .......... 62% 89.2M 1s
    ##  49600K .......... .......... .......... .......... .......... 62% 30.2M 1s
    ##  49650K .......... .......... .......... .......... .......... 62% 68.3M 1s
    ##  49700K .......... .......... .......... .......... .......... 62% 13.9M 1s
    ##  49750K .......... .......... .......... .......... .......... 62% 47.7M 1s
    ##  49800K .......... .......... .......... .......... .......... 62% 57.4M 1s
    ##  49850K .......... .......... .......... .......... .......... 62% 96.8M 1s
    ##  49900K .......... .......... .......... .......... .......... 62% 66.2M 1s
    ##  49950K .......... .......... .......... .......... .......... 62% 86.4M 1s
    ##  50000K .......... .......... .......... .......... .......... 62% 93.1M 1s
    ##  50050K .......... .......... .......... .......... .......... 62% 69.6M 1s
    ##  50100K .......... .......... .......... .......... .......... 62% 82.5M 1s
    ##  50150K .......... .......... .......... .......... .......... 62%  111M 1s
    ##  50200K .......... .......... .......... .......... .......... 62% 56.5M 1s
    ##  50250K .......... .......... .......... .......... .......... 62% 80.9M 1s
    ##  50300K .......... .......... .......... .......... .......... 62% 96.5M 1s
    ##  50350K .......... .......... .......... .......... .......... 63% 56.4M 1s
    ##  50400K .......... .......... .......... .......... .......... 63% 96.6M 1s
    ##  50450K .......... .......... .......... .......... .......... 63% 76.8M 1s
    ##  50500K .......... .......... .......... .......... .......... 63% 84.1M 1s
    ##  50550K .......... .......... .......... .......... .......... 63% 92.8M 1s
    ##  50600K .......... .......... .......... .......... .......... 63% 77.0M 1s
    ##  50650K .......... .......... .......... .......... .......... 63% 35.0M 1s
    ##  50700K .......... .......... .......... .......... .......... 63% 95.3M 1s
    ##  50750K .......... .......... .......... .......... .......... 63%  101M 1s
    ##  50800K .......... .......... .......... .......... .......... 63% 52.8M 1s
    ##  50850K .......... .......... .......... .......... .......... 63% 84.2M 1s
    ##  50900K .......... .......... .......... .......... .......... 63% 88.4M 1s
    ##  50950K .......... .......... .......... .......... .......... 63% 36.4M 1s
    ##  51000K .......... .......... .......... .......... .......... 63% 26.2M 1s
    ##  51050K .......... .......... .......... .......... .......... 63% 74.8M 1s
    ##  51100K .......... .......... .......... .......... .......... 63% 24.9M 1s
    ##  51150K .......... .......... .......... .......... .......... 64% 78.0M 1s
    ##  51200K .......... .......... .......... .......... .......... 64% 33.9M 1s
    ##  51250K .......... .......... .......... .......... .......... 64% 42.1M 1s
    ##  51300K .......... .......... .......... .......... .......... 64% 75.7M 1s
    ##  51350K .......... .......... .......... .......... .......... 64% 64.0M 1s
    ##  51400K .......... .......... .......... .......... .......... 64% 75.5M 1s
    ##  51450K .......... .......... .......... .......... .......... 64% 67.4M 1s
    ##  51500K .......... .......... .......... .......... .......... 64% 92.7M 1s
    ##  51550K .......... .......... .......... .......... .......... 64% 86.4M 1s
    ##  51600K .......... .......... .......... .......... .......... 64% 66.6M 1s
    ##  51650K .......... .......... .......... .......... .......... 64% 78.8M 1s
    ##  51700K .......... .......... .......... .......... .......... 64% 94.6M 1s
    ##  51750K .......... .......... .......... .......... .......... 64% 79.9M 1s
    ##  51800K .......... .......... .......... .......... .......... 64% 81.5M 1s
    ##  51850K .......... .......... .......... .......... .......... 64% 89.2M 1s
    ##  51900K .......... .......... .......... .......... .......... 65% 72.7M 1s
    ##  51950K .......... .......... .......... .......... .......... 65% 79.7M 1s
    ##  52000K .......... .......... .......... .......... .......... 65% 98.3M 1s
    ##  52050K .......... .......... .......... .......... .......... 65% 80.8M 1s
    ##  52100K .......... .......... .......... .......... .......... 65% 56.7M 1s
    ##  52150K .......... .......... .......... .......... .......... 65% 4.40M 1s
    ##  52200K .......... .......... .......... .......... .......... 65% 91.3M 1s
    ##  52250K .......... .......... .......... .......... .......... 65%  132M 1s
    ##  52300K .......... .......... .......... .......... .......... 65%  107M 1s
    ##  52350K .......... .......... .......... .......... .......... 65%  108M 1s
    ##  52400K .......... .......... .......... .......... .......... 65%  118M 1s
    ##  52450K .......... .......... .......... .......... .......... 65% 97.6M 1s
    ##  52500K .......... .......... .......... .......... .......... 65% 16.6M 1s
    ##  52550K .......... .......... .......... .......... .......... 65% 97.2M 1s
    ##  52600K .......... .......... .......... .......... .......... 65% 39.3M 1s
    ##  52650K .......... .......... .......... .......... .......... 65% 75.1M 1s
    ##  52700K .......... .......... .......... .......... .......... 66% 95.4M 1s
    ##  52750K .......... .......... .......... .......... .......... 66%  115M 0s
    ##  52800K .......... .......... .......... .......... .......... 66%  101M 0s
    ##  52850K .......... .......... .......... .......... .......... 66%  110M 0s
    ##  52900K .......... .......... .......... .......... .......... 66% 37.4M 0s
    ##  52950K .......... .......... .......... .......... .......... 66% 98.5M 0s
    ##  53000K .......... .......... .......... .......... .......... 66% 22.4M 0s
    ##  53050K .......... .......... .......... .......... .......... 66% 73.1M 0s
    ##  53100K .......... .......... .......... .......... .......... 66% 82.4M 0s
    ##  53150K .......... .......... .......... .......... .......... 66%  101M 0s
    ##  53200K .......... .......... .......... .......... .......... 66% 87.1M 0s
    ##  53250K .......... .......... .......... .......... .......... 66% 84.7M 0s
    ##  53300K .......... .......... .......... .......... .......... 66%  103M 0s
    ##  53350K .......... .......... .......... .......... .......... 66% 35.8M 0s
    ##  53400K .......... .......... .......... .......... .......... 66% 83.6M 0s
    ##  53450K .......... .......... .......... .......... .......... 66% 69.6M 0s
    ##  53500K .......... .......... .......... .......... .......... 67% 43.7M 0s
    ##  53550K .......... .......... .......... .......... .......... 67% 67.6M 0s
    ##  53600K .......... .......... .......... .......... .......... 67% 93.7M 0s
    ##  53650K .......... .......... .......... .......... .......... 67%  112M 0s
    ##  53700K .......... .......... .......... .......... .......... 67% 48.4M 0s
    ##  53750K .......... .......... .......... .......... .......... 67% 99.9M 0s
    ##  53800K .......... .......... .......... .......... .......... 67% 55.1M 0s
    ##  53850K .......... .......... .......... .......... .......... 67% 74.0M 0s
    ##  53900K .......... .......... .......... .......... .......... 67%  122M 0s
    ##  53950K .......... .......... .......... .......... .......... 67% 48.4M 0s
    ##  54000K .......... .......... .......... .......... .......... 67% 24.3M 0s
    ##  54050K .......... .......... .......... .......... .......... 67%  114M 0s
    ##  54100K .......... .......... .......... .......... .......... 67%  111M 0s
    ##  54150K .......... .......... .......... .......... .......... 67% 96.1M 0s
    ##  54200K .......... .......... .......... .......... .......... 67% 78.4M 0s
    ##  54250K .......... .......... .......... .......... .......... 67% 28.1M 0s
    ##  54300K .......... .......... .......... .......... .......... 68% 57.0M 0s
    ##  54350K .......... .......... .......... .......... .......... 68%  102M 0s
    ##  54400K .......... .......... .......... .......... .......... 68%  108M 0s
    ##  54450K .......... .......... .......... .......... .......... 68%  134M 0s
    ##  54500K .......... .......... .......... .......... .......... 68%  105M 0s
    ##  54550K .......... .......... .......... .......... .......... 68%  131M 0s
    ##  54600K .......... .......... .......... .......... .......... 68% 78.2M 0s
    ##  54650K .......... .......... .......... .......... .......... 68% 96.5M 0s
    ##  54700K .......... .......... .......... .......... .......... 68% 48.3M 0s
    ##  54750K .......... .......... .......... .......... .......... 68% 31.7M 0s
    ##  54800K .......... .......... .......... .......... .......... 68% 90.8M 0s
    ##  54850K .......... .......... .......... .......... .......... 68% 92.4M 0s
    ##  54900K .......... .......... .......... .......... .......... 68% 77.9M 0s
    ##  54950K .......... .......... .......... .......... .......... 68% 58.4M 0s
    ##  55000K .......... .......... .......... .......... .......... 68% 91.6M 0s
    ##  55050K .......... .......... .......... .......... .......... 68%  109M 0s
    ##  55100K .......... .......... .......... .......... .......... 69%  115M 0s
    ##  55150K .......... .......... .......... .......... .......... 69% 68.3M 0s
    ##  55200K .......... .......... .......... .......... .......... 69% 31.6M 0s
    ##  55250K .......... .......... .......... .......... .......... 69% 89.5M 0s
    ##  55300K .......... .......... .......... .......... .......... 69% 96.7M 0s
    ##  55350K .......... .......... .......... .......... .......... 69% 96.6M 0s
    ##  55400K .......... .......... .......... .......... .......... 69% 52.8M 0s
    ##  55450K .......... .......... .......... .......... .......... 69% 83.0M 0s
    ##  55500K .......... .......... .......... .......... .......... 69% 93.0M 0s
    ##  55550K .......... .......... .......... .......... .......... 69% 37.7M 0s
    ##  55600K .......... .......... .......... .......... .......... 69% 84.5M 0s
    ##  55650K .......... .......... .......... .......... .......... 69% 73.0M 0s
    ##  55700K .......... .......... .......... .......... .......... 69% 59.4M 0s
    ##  55750K .......... .......... .......... .......... .......... 69%  104M 0s
    ##  55800K .......... .......... .......... .......... .......... 69% 47.0M 0s
    ##  55850K .......... .......... .......... .......... .......... 69% 65.6M 0s
    ##  55900K .......... .......... .......... .......... .......... 70% 51.2M 0s
    ##  55950K .......... .......... .......... .......... .......... 70%  116M 0s
    ##  56000K .......... .......... .......... .......... .......... 70% 57.7M 0s
    ##  56050K .......... .......... .......... .......... .......... 70% 42.8M 0s
    ##  56100K .......... .......... .......... .......... .......... 70% 96.6M 0s
    ##  56150K .......... .......... .......... .......... .......... 70%  119M 0s
    ##  56200K .......... .......... .......... .......... .......... 70%  114M 0s
    ##  56250K .......... .......... .......... .......... .......... 70% 96.4M 0s
    ##  56300K .......... .......... .......... .......... .......... 70%  113M 0s
    ##  56350K .......... .......... .......... .......... .......... 70%  151M 0s
    ##  56400K .......... .......... .......... .......... .......... 70% 13.3M 0s
    ##  56450K .......... .......... .......... .......... .......... 70% 79.4M 0s
    ##  56500K .......... .......... .......... .......... .......... 70% 65.6M 0s
    ##  56550K .......... .......... .......... .......... .......... 70% 80.3M 0s
    ##  56600K .......... .......... .......... .......... .......... 70%  119M 0s
    ##  56650K .......... .......... .......... .......... .......... 70%  118M 0s
    ##  56700K .......... .......... .......... .......... .......... 71% 83.2M 0s
    ##  56750K .......... .......... .......... .......... .......... 71% 94.2M 0s
    ##  56800K .......... .......... .......... .......... .......... 71% 73.6M 0s
    ##  56850K .......... .......... .......... .......... .......... 71% 40.4M 0s
    ##  56900K .......... .......... .......... .......... .......... 71% 45.1M 0s
    ##  56950K .......... .......... .......... .......... .......... 71% 92.1M 0s
    ##  57000K .......... .......... .......... .......... .......... 71%  106M 0s
    ##  57050K .......... .......... .......... .......... .......... 71% 99.5M 0s
    ##  57100K .......... .......... .......... .......... .......... 71% 80.0M 0s
    ##  57150K .......... .......... .......... .......... .......... 71% 79.4M 0s
    ##  57200K .......... .......... .......... .......... .......... 71% 79.9M 0s
    ##  57250K .......... .......... .......... .......... .......... 71% 72.5M 0s
    ##  57300K .......... .......... .......... .......... .......... 71% 60.9M 0s
    ##  57350K .......... .......... .......... .......... .......... 71% 39.3M 0s
    ##  57400K .......... .......... .......... .......... .......... 71% 37.8M 0s
    ##  57450K .......... .......... .......... .......... .......... 71% 79.2M 0s
    ##  57500K .......... .......... .......... .......... .......... 72% 95.2M 0s
    ##  57550K .......... .......... .......... .......... .......... 72%  106M 0s
    ##  57600K .......... .......... .......... .......... .......... 72% 89.0M 0s
    ##  57650K .......... .......... .......... .......... .......... 72% 91.7M 0s
    ##  57700K .......... .......... .......... .......... .......... 72%  106M 0s
    ##  57750K .......... .......... .......... .......... .......... 72% 66.4M 0s
    ##  57800K .......... .......... .......... .......... .......... 72% 69.7M 0s
    ##  57850K .......... .......... .......... .......... .......... 72% 35.8M 0s
    ##  57900K .......... .......... .......... .......... .......... 72% 44.0M 0s
    ##  57950K .......... .......... .......... .......... .......... 72% 82.0M 0s
    ##  58000K .......... .......... .......... .......... .......... 72% 82.1M 0s
    ##  58050K .......... .......... .......... .......... .......... 72% 89.6M 0s
    ##  58100K .......... .......... .......... .......... .......... 72% 53.9M 0s
    ##  58150K .......... .......... .......... .......... .......... 72% 38.6M 0s
    ##  58200K .......... .......... .......... .......... .......... 72% 52.5M 0s
    ##  58250K .......... .......... .......... .......... .......... 72% 62.0M 0s
    ##  58300K .......... .......... .......... .......... .......... 73% 98.0M 0s
    ##  58350K .......... .......... .......... .......... .......... 73%  109M 0s
    ##  58400K .......... .......... .......... .......... .......... 73% 84.6M 0s
    ##  58450K .......... .......... .......... .......... .......... 73%  143M 0s
    ##  58500K .......... .......... .......... .......... .......... 73%  110M 0s
    ##  58550K .......... .......... .......... .......... .......... 73%  113M 0s
    ##  58600K .......... .......... .......... .......... .......... 73% 75.9M 0s
    ##  58650K .......... .......... .......... .......... .......... 73% 74.2M 0s
    ##  58700K .......... .......... .......... .......... .......... 73% 72.4M 0s
    ##  58750K .......... .......... .......... .......... .......... 73%  131M 0s
    ##  58800K .......... .......... .......... .......... .......... 73% 84.6M 0s
    ##  58850K .......... .......... .......... .......... .......... 73% 90.0M 0s
    ##  58900K .......... .......... .......... .......... .......... 73% 83.8M 0s
    ##  58950K .......... .......... .......... .......... .......... 73% 83.6M 0s
    ##  59000K .......... .......... .......... .......... .......... 73% 45.5M 0s
    ##  59050K .......... .......... .......... .......... .......... 73% 58.2M 0s
    ##  59100K .......... .......... .......... .......... .......... 74% 67.8M 0s
    ##  59150K .......... .......... .......... .......... .......... 74% 62.6M 0s
    ##  59200K .......... .......... .......... .......... .......... 74% 68.7M 0s
    ##  59250K .......... .......... .......... .......... .......... 74% 87.0M 0s
    ##  59300K .......... .......... .......... .......... .......... 74% 71.1M 0s
    ##  59350K .......... .......... .......... .......... .......... 74% 99.7M 0s
    ##  59400K .......... .......... .......... .......... .......... 74% 66.5M 0s
    ##  59450K .......... .......... .......... .......... .......... 74% 49.9M 0s
    ##  59500K .......... .......... .......... .......... .......... 74% 58.5M 0s
    ##  59550K .......... .......... .......... .......... .......... 74%  114M 0s
    ##  59600K .......... .......... .......... .......... .......... 74% 44.6M 0s
    ##  59650K .......... .......... .......... .......... .......... 74% 85.5M 0s
    ##  59700K .......... .......... .......... .......... .......... 74% 94.0M 0s
    ##  59750K .......... .......... .......... .......... .......... 74% 77.8M 0s
    ##  59800K .......... .......... .......... .......... .......... 74%  106M 0s
    ##  59850K .......... .......... .......... .......... .......... 74%  131M 0s
    ##  59900K .......... .......... .......... .......... .......... 75% 55.5M 0s
    ##  59950K .......... .......... .......... .......... .......... 75% 98.5M 0s
    ##  60000K .......... .......... .......... .......... .......... 75% 94.4M 0s
    ##  60050K .......... .......... .......... .......... .......... 75% 52.6M 0s
    ##  60100K .......... .......... .......... .......... .......... 75% 22.1M 0s
    ##  60150K .......... .......... .......... .......... .......... 75%  141M 0s
    ##  60200K .......... .......... .......... .......... .......... 75% 46.4M 0s
    ##  60250K .......... .......... .......... .......... .......... 75% 56.9M 0s
    ##  60300K .......... .......... .......... .......... .......... 75%  106M 0s
    ##  60350K .......... .......... .......... .......... .......... 75%  132M 0s
    ##  60400K .......... .......... .......... .......... .......... 75% 68.1M 0s
    ##  60450K .......... .......... .......... .......... .......... 75%  144M 0s
    ##  60500K .......... .......... .......... .......... .......... 75%  144M 0s
    ##  60550K .......... .......... .......... .......... .......... 75% 69.1M 0s
    ##  60600K .......... .......... .......... .......... .......... 75% 40.1M 0s
    ##  60650K .......... .......... .......... .......... .......... 75%  100M 0s
    ##  60700K .......... .......... .......... .......... .......... 76%  122M 0s
    ##  60750K .......... .......... .......... .......... .......... 76%  145M 0s
    ##  60800K .......... .......... .......... .......... .......... 76% 73.1M 0s
    ##  60850K .......... .......... .......... .......... .......... 76%  101M 0s
    ##  60900K .......... .......... .......... .......... .......... 76% 35.9M 0s
    ##  60950K .......... .......... .......... .......... .......... 76%  110M 0s
    ##  61000K .......... .......... .......... .......... .......... 76% 26.7M 0s
    ##  61050K .......... .......... .......... .......... .......... 76%  130M 0s
    ##  61100K .......... .......... .......... .......... .......... 76%  110M 0s
    ##  61150K .......... .......... .......... .......... .......... 76% 69.8M 0s
    ##  61200K .......... .......... .......... .......... .......... 76% 90.0M 0s
    ##  61250K .......... .......... .......... .......... .......... 76%  131M 0s
    ##  61300K .......... .......... .......... .......... .......... 76% 98.3M 0s
    ##  61350K .......... .......... .......... .......... .......... 76% 74.4M 0s
    ##  61400K .......... .......... .......... .......... .......... 76%  116M 0s
    ##  61450K .......... .......... .......... .......... .......... 76% 20.3M 0s
    ##  61500K .......... .......... .......... .......... .......... 77% 69.0M 0s
    ##  61550K .......... .......... .......... .......... .......... 77%  120M 0s
    ##  61600K .......... .......... .......... .......... .......... 77%  126M 0s
    ##  61650K .......... .......... .......... .......... .......... 77%  132M 0s
    ##  61700K .......... .......... .......... .......... .......... 77%  105M 0s
    ##  61750K .......... .......... .......... .......... .......... 77%  119M 0s
    ##  61800K .......... .......... .......... .......... .......... 77% 90.2M 0s
    ##  61850K .......... .......... .......... .......... .......... 77% 84.8M 0s
    ##  61900K .......... .......... .......... .......... .......... 77% 36.9M 0s
    ##  61950K .......... .......... .......... .......... .......... 77%  112M 0s
    ##  62000K .......... .......... .......... .......... .......... 77% 47.8M 0s
    ##  62050K .......... .......... .......... .......... .......... 77%  111M 0s
    ##  62100K .......... .......... .......... .......... .......... 77%  101M 0s
    ##  62150K .......... .......... .......... .......... .......... 77%  161M 0s
    ##  62200K .......... .......... .......... .......... .......... 77%  116M 0s
    ##  62250K .......... .......... .......... .......... .......... 77% 85.1M 0s
    ##  62300K .......... .......... .......... .......... .......... 78%  136M 0s
    ##  62350K .......... .......... .......... .......... .......... 78%  157M 0s
    ##  62400K .......... .......... .......... .......... .......... 78% 35.4M 0s
    ##  62450K .......... .......... .......... .......... .......... 78% 60.3M 0s
    ##  62500K .......... .......... .......... .......... .......... 78% 84.8M 0s
    ##  62550K .......... .......... .......... .......... .......... 78% 40.2M 0s
    ##  62600K .......... .......... .......... .......... .......... 78% 58.0M 0s
    ##  62650K .......... .......... .......... .......... .......... 78% 98.5M 0s
    ##  62700K .......... .......... .......... .......... .......... 78%  102M 0s
    ##  62750K .......... .......... .......... .......... .......... 78% 84.3M 0s
    ##  62800K .......... .......... .......... .......... .......... 78% 62.1M 0s
    ##  62850K .......... .......... .......... .......... .......... 78% 82.1M 0s
    ##  62900K .......... .......... .......... .......... .......... 78%  111M 0s
    ##  62950K .......... .......... .......... .......... .......... 78% 45.1M 0s
    ##  63000K .......... .......... .......... .......... .......... 78% 91.2M 0s
    ##  63050K .......... .......... .......... .......... .......... 78% 53.0M 0s
    ##  63100K .......... .......... .......... .......... .......... 79% 67.8M 0s
    ##  63150K .......... .......... .......... .......... .......... 79%  151M 0s
    ##  63200K .......... .......... .......... .......... .......... 79%  104M 0s
    ##  63250K .......... .......... .......... .......... .......... 79% 79.3M 0s
    ##  63300K .......... .......... .......... .......... .......... 79% 62.2M 0s
    ##  63350K .......... .......... .......... .......... .......... 79% 48.0M 0s
    ##  63400K .......... .......... .......... .......... .......... 79%  110M 0s
    ##  63450K .......... .......... .......... .......... .......... 79% 66.6M 0s
    ##  63500K .......... .......... .......... .......... .......... 79% 69.8M 0s
    ##  63550K .......... .......... .......... .......... .......... 79%  170M 0s
    ##  63600K .......... .......... .......... .......... .......... 79% 40.6M 0s
    ##  63650K .......... .......... .......... .......... .......... 79% 71.3M 0s
    ##  63700K .......... .......... .......... .......... .......... 79% 53.8M 0s
    ##  63750K .......... .......... .......... .......... .......... 79% 4.84M 0s
    ##  63800K .......... .......... .......... .......... .......... 79% 30.0M 0s
    ##  63850K .......... .......... .......... .......... .......... 79% 65.7M 0s
    ##  63900K .......... .......... .......... .......... .......... 80% 64.2M 0s
    ##  63950K .......... .......... .......... .......... .......... 80% 55.2M 0s
    ##  64000K .......... .......... .......... .......... .......... 80% 50.8M 0s
    ##  64050K .......... .......... .......... .......... .......... 80% 53.5M 0s
    ##  64100K .......... .......... .......... .......... .......... 80% 26.8M 0s
    ##  64150K .......... .......... .......... .......... .......... 80% 71.6M 0s
    ##  64200K .......... .......... .......... .......... .......... 80% 72.0M 0s
    ##  64250K .......... .......... .......... .......... .......... 80% 70.9M 0s
    ##  64300K .......... .......... .......... .......... .......... 80% 66.1M 0s
    ##  64350K .......... .......... .......... .......... .......... 80% 67.0M 0s
    ##  64400K .......... .......... .......... .......... .......... 80% 94.6M 0s
    ##  64450K .......... .......... .......... .......... .......... 80%  108M 0s
    ##  64500K .......... .......... .......... .......... .......... 80% 90.3M 0s
    ##  64550K .......... .......... .......... .......... .......... 80% 36.7M 0s
    ##  64600K .......... .......... .......... .......... .......... 80% 67.8M 0s
    ##  64650K .......... .......... .......... .......... .......... 80% 91.0M 0s
    ##  64700K .......... .......... .......... .......... .......... 81% 71.3M 0s
    ##  64750K .......... .......... .......... .......... .......... 81% 70.4M 0s
    ##  64800K .......... .......... .......... .......... .......... 81% 74.8M 0s
    ##  64850K .......... .......... .......... .......... .......... 81% 81.2M 0s
    ##  64900K .......... .......... .......... .......... .......... 81% 69.0M 0s
    ##  64950K .......... .......... .......... .......... .......... 81% 73.4M 0s
    ##  65000K .......... .......... .......... .......... .......... 81% 57.2M 0s
    ##  65050K .......... .......... .......... .......... .......... 81% 50.2M 0s
    ##  65100K .......... .......... .......... .......... .......... 81% 65.5M 0s
    ##  65150K .......... .......... .......... .......... .......... 81% 39.1M 0s
    ##  65200K .......... .......... .......... .......... .......... 81% 99.8M 0s
    ##  65250K .......... .......... .......... .......... .......... 81% 86.2M 0s
    ##  65300K .......... .......... .......... .......... .......... 81% 67.8M 0s
    ##  65350K .......... .......... .......... .......... .......... 81% 51.1M 0s
    ##  65400K .......... .......... .......... .......... .......... 81% 96.2M 0s
    ##  65450K .......... .......... .......... .......... .......... 81% 92.7M 0s
    ##  65500K .......... .......... .......... .......... .......... 82% 71.0M 0s
    ##  65550K .......... .......... .......... .......... .......... 82% 61.5M 0s
    ##  65600K .......... .......... .......... .......... .......... 82% 29.1M 0s
    ##  65650K .......... .......... .......... .......... .......... 82%  125M 0s
    ##  65700K .......... .......... .......... .......... .......... 82% 42.7M 0s
    ##  65750K .......... .......... .......... .......... .......... 82% 69.3M 0s
    ##  65800K .......... .......... .......... .......... .......... 82% 73.8M 0s
    ##  65850K .......... .......... .......... .......... .......... 82% 9.98M 0s
    ##  65900K .......... .......... .......... .......... .......... 82%  118M 0s
    ##  65950K .......... .......... .......... .......... .......... 82%  114M 0s
    ##  66000K .......... .......... .......... .......... .......... 82%  123M 0s
    ##  66050K .......... .......... .......... .......... .......... 82%  146M 0s
    ##  66100K .......... .......... .......... .......... .......... 82%  151M 0s
    ##  66150K .......... .......... .......... .......... .......... 82%  165M 0s
    ##  66200K .......... .......... .......... .......... .......... 82%  144M 0s
    ##  66250K .......... .......... .......... .......... .......... 82%  188M 0s
    ##  66300K .......... .......... .......... .......... .......... 83%  143M 0s
    ##  66350K .......... .......... .......... .......... .......... 83%  157M 0s
    ##  66400K .......... .......... .......... .......... .......... 83%  108M 0s
    ##  66450K .......... .......... .......... .......... .......... 83% 42.7M 0s
    ##  66500K .......... .......... .......... .......... .......... 83% 46.5M 0s
    ##  66550K .......... .......... .......... .......... .......... 83% 34.1M 0s
    ##  66600K .......... .......... .......... .......... .......... 83% 81.7M 0s
    ##  66650K .......... .......... .......... .......... .......... 83% 81.0M 0s
    ##  66700K .......... .......... .......... .......... .......... 83%  110M 0s
    ##  66750K .......... .......... .......... .......... .......... 83% 94.6M 0s
    ##  66800K .......... .......... .......... .......... .......... 83%  103M 0s
    ##  66850K .......... .......... .......... .......... .......... 83%  107M 0s
    ##  66900K .......... .......... .......... .......... .......... 83% 14.9M 0s
    ##  66950K .......... .......... .......... .......... .......... 83% 87.6M 0s
    ##  67000K .......... .......... .......... .......... .......... 83%  124M 0s
    ##  67050K .......... .......... .......... .......... .......... 83%  138M 0s
    ##  67100K .......... .......... .......... .......... .......... 84%  116M 0s
    ##  67150K .......... .......... .......... .......... .......... 84%  127M 0s
    ##  67200K .......... .......... .......... .......... .......... 84%  118M 0s
    ##  67250K .......... .......... .......... .......... .......... 84% 32.7M 0s
    ##  67300K .......... .......... .......... .......... .......... 84% 32.4M 0s
    ##  67350K .......... .......... .......... .......... .......... 84% 27.2M 0s
    ##  67400K .......... .......... .......... .......... .......... 84%  106M 0s
    ##  67450K .......... .......... .......... .......... .......... 84% 62.9M 0s
    ##  67500K .......... .......... .......... .......... .......... 84% 88.4M 0s
    ##  67550K .......... .......... .......... .......... .......... 84%  148M 0s
    ##  67600K .......... .......... .......... .......... .......... 84% 38.5M 0s
    ##  67650K .......... .......... .......... .......... .......... 84%  149M 0s
    ##  67700K .......... .......... .......... .......... .......... 84%  143M 0s
    ##  67750K .......... .......... .......... .......... .......... 84%  102M 0s
    ##  67800K .......... .......... .......... .......... .......... 84%  114M 0s
    ##  67850K .......... .......... .......... .......... .......... 84%  122M 0s
    ##  67900K .......... .......... .......... .......... .......... 85%  103M 0s
    ##  67950K .......... .......... .......... .......... .......... 85% 49.4M 0s
    ##  68000K .......... .......... .......... .......... .......... 85%  159M 0s
    ##  68050K .......... .......... .......... .......... .......... 85% 47.9M 0s
    ##  68100K .......... .......... .......... .......... .......... 85% 80.1M 0s
    ##  68150K .......... .......... .......... .......... .......... 85% 54.1M 0s
    ##  68200K .......... .......... .......... .......... .......... 85% 60.7M 0s
    ##  68250K .......... .......... .......... .......... .......... 85%  104M 0s
    ##  68300K .......... .......... .......... .......... .......... 85% 33.7M 0s
    ##  68350K .......... .......... .......... .......... .......... 85% 25.5M 0s
    ##  68400K .......... .......... .......... .......... .......... 85%  128M 0s
    ##  68450K .......... .......... .......... .......... .......... 85%  147M 0s
    ##  68500K .......... .......... .......... .......... .......... 85%  134M 0s
    ##  68550K .......... .......... .......... .......... .......... 85%  151M 0s
    ##  68600K .......... .......... .......... .......... .......... 85%  112M 0s
    ##  68650K .......... .......... .......... .......... .......... 85% 12.1M 0s
    ##  68700K .......... .......... .......... .......... .......... 86%  115M 0s
    ##  68750K .......... .......... .......... .......... .......... 86%  142M 0s
    ##  68800K .......... .......... .......... .......... .......... 86%  111M 0s
    ##  68850K .......... .......... .......... .......... .......... 86%  140M 0s
    ##  68900K .......... .......... .......... .......... .......... 86%  146M 0s
    ##  68950K .......... .......... .......... .......... .......... 86%  120M 0s
    ##  69000K .......... .......... .......... .......... .......... 86%  111M 0s
    ##  69050K .......... .......... .......... .......... .......... 86% 29.3M 0s
    ##  69100K .......... .......... .......... .......... .......... 86% 17.4M 0s
    ##  69150K .......... .......... .......... .......... .......... 86% 43.6M 0s
    ##  69200K .......... .......... .......... .......... .......... 86% 63.3M 0s
    ##  69250K .......... .......... .......... .......... .......... 86% 45.1M 0s
    ##  69300K .......... .......... .......... .......... .......... 86% 76.6M 0s
    ##  69350K .......... .......... .......... .......... .......... 86% 88.0M 0s
    ##  69400K .......... .......... .......... .......... .......... 86% 74.1M 0s
    ##  69450K .......... .......... .......... .......... .......... 86% 83.5M 0s
    ##  69500K .......... .......... .......... .......... .......... 87% 74.8M 0s
    ##  69550K .......... .......... .......... .......... .......... 87% 93.7M 0s
    ##  69600K .......... .......... .......... .......... .......... 87% 97.5M 0s
    ##  69650K .......... .......... .......... .......... .......... 87% 90.8M 0s
    ##  69700K .......... .......... .......... .......... .......... 87% 79.7M 0s
    ##  69750K .......... .......... .......... .......... .......... 87% 92.5M 0s
    ##  69800K .......... .......... .......... .......... .......... 87% 45.5M 0s
    ##  69850K .......... .......... .......... .......... .......... 87% 75.4M 0s
    ##  69900K .......... .......... .......... .......... .......... 87% 30.5M 0s
    ##  69950K .......... .......... .......... .......... .......... 87% 62.5M 0s
    ##  70000K .......... .......... .......... .......... .......... 87% 77.6M 0s
    ##  70050K .......... .......... .......... .......... .......... 87% 73.5M 0s
    ##  70100K .......... .......... .......... .......... .......... 87% 23.0M 0s
    ##  70150K .......... .......... .......... .......... .......... 87% 77.0M 0s
    ##  70200K .......... .......... .......... .......... .......... 87% 83.7M 0s
    ##  70250K .......... .......... .......... .......... .......... 87% 87.9M 0s
    ##  70300K .......... .......... .......... .......... .......... 88% 93.0M 0s
    ##  70350K .......... .......... .......... .......... .......... 88%  102M 0s
    ##  70400K .......... .......... .......... .......... .......... 88% 93.6M 0s
    ##  70450K .......... .......... .......... .......... .......... 88% 95.4M 0s
    ##  70500K .......... .......... .......... .......... .......... 88% 40.3M 0s
    ##  70550K .......... .......... .......... .......... .......... 88%  100M 0s
    ##  70600K .......... .......... .......... .......... .......... 88% 76.0M 0s
    ##  70650K .......... .......... .......... .......... .......... 88% 69.6M 0s
    ##  70700K .......... .......... .......... .......... .......... 88% 88.3M 0s
    ##  70750K .......... .......... .......... .......... .......... 88% 78.3M 0s
    ##  70800K .......... .......... .......... .......... .......... 88% 37.8M 0s
    ##  70850K .......... .......... .......... .......... .......... 88% 92.0M 0s
    ##  70900K .......... .......... .......... .......... .......... 88% 85.0M 0s
    ##  70950K .......... .......... .......... .......... .......... 88%  105M 0s
    ##  71000K .......... .......... .......... .......... .......... 88% 95.5M 0s
    ##  71050K .......... .......... .......... .......... .......... 88% 92.0M 0s
    ##  71100K .......... .......... .......... .......... .......... 89%  103M 0s
    ##  71150K .......... .......... .......... .......... .......... 89%  104M 0s
    ##  71200K .......... .......... .......... .......... .......... 89%  102M 0s
    ##  71250K .......... .......... .......... .......... .......... 89% 74.1M 0s
    ##  71300K .......... .......... .......... .......... .......... 89% 48.5M 0s
    ##  71350K .......... .......... .......... .......... .......... 89% 60.0M 0s
    ##  71400K .......... .......... .......... .......... .......... 89% 74.5M 0s
    ##  71450K .......... .......... .......... .......... .......... 89% 73.7M 0s
    ##  71500K .......... .......... .......... .......... .......... 89% 31.5M 0s
    ##  71550K .......... .......... .......... .......... .......... 89% 79.2M 0s
    ##  71600K .......... .......... .......... .......... .......... 89% 97.2M 0s
    ##  71650K .......... .......... .......... .......... .......... 89% 67.3M 0s
    ##  71700K .......... .......... .......... .......... .......... 89% 44.4M 0s
    ##  71750K .......... .......... .......... .......... .......... 89% 63.4M 0s
    ##  71800K .......... .......... .......... .......... .......... 89% 96.7M 0s
    ##  71850K .......... .......... .......... .......... .......... 89% 75.5M 0s
    ##  71900K .......... .......... .......... .......... .......... 90% 35.2M 0s
    ##  71950K .......... .......... .......... .......... .......... 90% 76.8M 0s
    ##  72000K .......... .......... .......... .......... .......... 90% 70.4M 0s
    ##  72050K .......... .......... .......... .......... .......... 90% 89.3M 0s
    ##  72100K .......... .......... .......... .......... .......... 90% 47.0M 0s
    ##  72150K .......... .......... .......... .......... .......... 90% 99.8M 0s
    ##  72200K .......... .......... .......... .......... .......... 90% 90.0M 0s
    ##  72250K .......... .......... .......... .......... .......... 90% 37.7M 0s
    ##  72300K .......... .......... .......... .......... .......... 90% 59.9M 0s
    ##  72350K .......... .......... .......... .......... .......... 90% 81.7M 0s
    ##  72400K .......... .......... .......... .......... .......... 90% 78.7M 0s
    ##  72450K .......... .......... .......... .......... .......... 90% 85.7M 0s
    ##  72500K .......... .......... .......... .......... .......... 90% 69.1M 0s
    ##  72550K .......... .......... .......... .......... .......... 90% 74.7M 0s
    ##  72600K .......... .......... .......... .......... .......... 90% 57.9M 0s
    ##  72650K .......... .......... .......... .......... .......... 90% 69.8M 0s
    ##  72700K .......... .......... .......... .......... .......... 91% 55.6M 0s
    ##  72750K .......... .......... .......... .......... .......... 91%  104M 0s
    ##  72800K .......... .......... .......... .......... .......... 91% 81.7M 0s
    ##  72850K .......... .......... .......... .......... .......... 91% 33.1M 0s
    ##  72900K .......... .......... .......... .......... .......... 91%  107M 0s
    ##  72950K .......... .......... .......... .......... .......... 91% 69.8M 0s
    ##  73000K .......... .......... .......... .......... .......... 91% 97.3M 0s
    ##  73050K .......... .......... .......... .......... .......... 91% 82.0M 0s
    ##  73100K .......... .......... .......... .......... .......... 91% 85.1M 0s
    ##  73150K .......... .......... .......... .......... .......... 91% 89.1M 0s
    ##  73200K .......... .......... .......... .......... .......... 91% 76.0M 0s
    ##  73250K .......... .......... .......... .......... .......... 91% 35.3M 0s
    ##  73300K .......... .......... .......... .......... .......... 91%  107M 0s
    ##  73350K .......... .......... .......... .......... .......... 91% 92.1M 0s
    ##  73400K .......... .......... .......... .......... .......... 91% 57.0M 0s
    ##  73450K .......... .......... .......... .......... .......... 91% 84.0M 0s
    ##  73500K .......... .......... .......... .......... .......... 92% 97.8M 0s
    ##  73550K .......... .......... .......... .......... .......... 92% 94.4M 0s
    ##  73600K .......... .......... .......... .......... .......... 92% 45.2M 0s
    ##  73650K .......... .......... .......... .......... .......... 92% 87.7M 0s
    ##  73700K .......... .......... .......... .......... .......... 92% 91.7M 0s
    ##  73750K .......... .......... .......... .......... .......... 92%  113M 0s
    ##  73800K .......... .......... .......... .......... .......... 92% 62.9M 0s
    ##  73850K .......... .......... .......... .......... .......... 92% 79.1M 0s
    ##  73900K .......... .......... .......... .......... .......... 92% 85.7M 0s
    ##  73950K .......... .......... .......... .......... .......... 92%  104M 0s
    ##  74000K .......... .......... .......... .......... .......... 92% 48.8M 0s
    ##  74050K .......... .......... .......... .......... .......... 92% 79.6M 0s
    ##  74100K .......... .......... .......... .......... .......... 92% 85.0M 0s
    ##  74150K .......... .......... .......... .......... .......... 92%  128M 0s
    ##  74200K .......... .......... .......... .......... .......... 92% 52.8M 0s
    ##  74250K .......... .......... .......... .......... .......... 92%  102M 0s
    ##  74300K .......... .......... .......... .......... .......... 93% 70.5M 0s
    ##  74350K .......... .......... .......... .......... .......... 93% 81.1M 0s
    ##  74400K .......... .......... .......... .......... .......... 93% 63.8M 0s
    ##  74450K .......... .......... .......... .......... .......... 93% 95.3M 0s
    ##  74500K .......... .......... .......... .......... .......... 93% 96.9M 0s
    ##  74550K .......... .......... .......... .......... .......... 93%  113M 0s
    ##  74600K .......... .......... .......... .......... .......... 93% 52.5M 0s
    ##  74650K .......... .......... .......... .......... .......... 93% 96.2M 0s
    ##  74700K .......... .......... .......... .......... .......... 93% 89.1M 0s
    ##  74750K .......... .......... .......... .......... .......... 93% 70.5M 0s
    ##  74800K .......... .......... .......... .......... .......... 93% 75.5M 0s
    ##  74850K .......... .......... .......... .......... .......... 93% 95.1M 0s
    ##  74900K .......... .......... .......... .......... .......... 93% 93.2M 0s
    ##  74950K .......... .......... .......... .......... .......... 93% 46.5M 0s
    ##  75000K .......... .......... .......... .......... .......... 93% 91.5M 0s
    ##  75050K .......... .......... .......... .......... .......... 93% 76.1M 0s
    ##  75100K .......... .......... .......... .......... .......... 94% 84.5M 0s
    ##  75150K .......... .......... .......... .......... .......... 94% 98.0M 0s
    ##  75200K .......... .......... .......... .......... .......... 94%  103M 0s
    ##  75250K .......... .......... .......... .......... .......... 94% 91.1M 0s
    ##  75300K .......... .......... .......... .......... .......... 94% 87.9M 0s
    ##  75350K .......... .......... .......... .......... .......... 94%  119M 0s
    ##  75400K .......... .......... .......... .......... .......... 94% 92.8M 0s
    ##  75450K .......... .......... .......... .......... .......... 94% 88.4M 0s
    ##  75500K .......... .......... .......... .......... .......... 94% 76.2M 0s
    ##  75550K .......... .......... .......... .......... .......... 94% 83.5M 0s
    ##  75600K .......... .......... .......... .......... .......... 94% 84.9M 0s
    ##  75650K .......... .......... .......... .......... .......... 94%  105M 0s
    ##  75700K .......... .......... .......... .......... .......... 94% 83.7M 0s
    ##  75750K .......... .......... .......... .......... .......... 94% 96.0M 0s
    ##  75800K .......... .......... .......... .......... .......... 94% 98.1M 0s
    ##  75850K .......... .......... .......... .......... .......... 94%  109M 0s
    ##  75900K .......... .......... .......... .......... .......... 95% 60.0M 0s
    ##  75950K .......... .......... .......... .......... .......... 95%  118M 0s
    ##  76000K .......... .......... .......... .......... .......... 95% 53.5M 0s
    ##  76050K .......... .......... .......... .......... .......... 95%  142M 0s
    ##  76100K .......... .......... .......... .......... .......... 95% 95.7M 0s
    ##  76150K .......... .......... .......... .......... .......... 95% 82.0M 0s
    ##  76200K .......... .......... .......... .......... .......... 95% 87.2M 0s
    ##  76250K .......... .......... .......... .......... .......... 95%  129M 0s
    ##  76300K .......... .......... .......... .......... .......... 95%  118M 0s
    ##  76350K .......... .......... .......... .......... .......... 95% 71.9M 0s
    ##  76400K .......... .......... .......... .......... .......... 95% 58.9M 0s
    ##  76450K .......... .......... .......... .......... .......... 95% 44.6M 0s
    ##  76500K .......... .......... .......... .......... .......... 95% 85.0M 0s
    ##  76550K .......... .......... .......... .......... .......... 95%  130M 0s
    ##  76600K .......... .......... .......... .......... .......... 95%  110M 0s
    ##  76650K .......... .......... .......... .......... .......... 95% 95.6M 0s
    ##  76700K .......... .......... .......... .......... .......... 96% 44.2M 0s
    ##  76750K .......... .......... .......... .......... .......... 96% 43.9M 0s
    ##  76800K .......... .......... .......... .......... .......... 96% 99.8M 0s
    ##  76850K .......... .......... .......... .......... .......... 96% 58.7M 0s
    ##  76900K .......... .......... .......... .......... .......... 96% 29.0M 0s
    ##  76950K .......... .......... .......... .......... .......... 96% 87.1M 0s
    ##  77000K .......... .......... .......... .......... .......... 96%  100M 0s
    ##  77050K .......... .......... .......... .......... .......... 96%  118M 0s
    ##  77100K .......... .......... .......... .......... .......... 96%  131M 0s
    ##  77150K .......... .......... .......... .......... .......... 96% 94.6M 0s
    ##  77200K .......... .......... .......... .......... .......... 96%  107M 0s
    ##  77250K .......... .......... .......... .......... .......... 96% 82.1M 0s
    ##  77300K .......... .......... .......... .......... .......... 96%  103M 0s
    ##  77350K .......... .......... .......... .......... .......... 96% 97.5M 0s
    ##  77400K .......... .......... .......... .......... .......... 96% 84.6M 0s
    ##  77450K .......... .......... .......... .......... .......... 96% 30.2M 0s
    ##  77500K .......... .......... .......... .......... .......... 97% 69.8M 0s
    ##  77550K .......... .......... .......... .......... .......... 97% 90.8M 0s
    ##  77600K .......... .......... .......... .......... .......... 97% 62.1M 0s
    ##  77650K .......... .......... .......... .......... .......... 97% 87.4M 0s
    ##  77700K .......... .......... .......... .......... .......... 97%  119M 0s
    ##  77750K .......... .......... .......... .......... .......... 97% 93.5M 0s
    ##  77800K .......... .......... .......... .......... .......... 97% 36.4M 0s
    ##  77850K .......... .......... .......... .......... .......... 97% 58.0M 0s
    ##  77900K .......... .......... .......... .......... .......... 97% 82.5M 0s
    ##  77950K .......... .......... .......... .......... .......... 97% 95.0M 0s
    ##  78000K .......... .......... .......... .......... .......... 97% 53.9M 0s
    ##  78050K .......... .......... .......... .......... .......... 97% 69.5M 0s
    ##  78100K .......... .......... .......... .......... .......... 97% 90.8M 0s
    ##  78150K .......... .......... .......... .......... .......... 97%  107M 0s
    ##  78200K .......... .......... .......... .......... .......... 97% 62.5M 0s
    ##  78250K .......... .......... .......... .......... .......... 97% 61.5M 0s
    ##  78300K .......... .......... .......... .......... .......... 98% 99.9M 0s
    ##  78350K .......... .......... .......... .......... .......... 98% 46.2M 0s
    ##  78400K .......... .......... .......... .......... .......... 98% 77.0M 0s
    ##  78450K .......... .......... .......... .......... .......... 98% 90.2M 0s
    ##  78500K .......... .......... .......... .......... .......... 98%  108M 0s
    ##  78550K .......... .......... .......... .......... .......... 98% 62.4M 0s
    ##  78600K .......... .......... .......... .......... .......... 98% 68.4M 0s
    ##  78650K .......... .......... .......... .......... .......... 98% 84.1M 0s
    ##  78700K .......... .......... .......... .......... .......... 98% 94.8M 0s
    ##  78750K .......... .......... .......... .......... .......... 98% 69.7M 0s
    ##  78800K .......... .......... .......... .......... .......... 98% 52.6M 0s
    ##  78850K .......... .......... .......... .......... .......... 98%  108M 0s
    ##  78900K .......... .......... .......... .......... .......... 98%  130M 0s
    ##  78950K .......... .......... .......... .......... .......... 98% 37.5M 0s
    ##  79000K .......... .......... .......... .......... .......... 98%  101M 0s
    ##  79050K .......... .......... .......... .......... .......... 98% 90.8M 0s
    ##  79100K .......... .......... .......... .......... .......... 99%  135M 0s
    ##  79150K .......... .......... .......... .......... .......... 99% 93.9M 0s
    ##  79200K .......... .......... .......... .......... .......... 99%  142M 0s
    ##  79250K .......... .......... .......... .......... .......... 99% 28.0M 0s
    ##  79300K .......... .......... .......... .......... .......... 99%  112M 0s
    ##  79350K .......... .......... .......... .......... .......... 99%  146M 0s
    ##  79400K .......... .......... .......... .......... .......... 99%  111M 0s
    ##  79450K .......... .......... .......... .......... .......... 99% 72.7M 0s
    ##  79500K .......... .......... .......... .......... .......... 99%  115M 0s
    ##  79550K .......... .......... .......... .......... .......... 99% 47.1M 0s
    ##  79600K .......... .......... .......... .......... .......... 99%  102M 0s
    ##  79650K .......... .......... .......... .......... .......... 99%  144M 0s
    ##  79700K .......... .......... .......... .......... .......... 99% 25.1M 0s
    ##  79750K .......... .......... .......... .......... .......... 99%  111M 0s
    ##  79800K .......... .......... .......... .......... .......... 99%  131M 0s
    ##  79850K .......... .......... .......... .......... .......... 99%  116M 0s
    ##  79900K .......... .......... ..                              100% 99.3M=1.4s
    ## 
    ## 2020-11-27 13:48:20 (57.0 MB/s) - ‘silva_species_assignment_v138.fa.gz.3’ saved [81840166/81840166]

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

# Evaluate accuracy

## Ici on cherche à comparer les variants donnés par la machine par rapport à la composition de communauté attendue. Le premier code nous montre qu’on s’attendait à avoir 20 souches différentes.

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

# Sauvegarde de l’environnement pour pouvoir jouer le code de phyloseq

``` r
save.image(file="03_DADA2_tutorial_FinalEnv")
```
