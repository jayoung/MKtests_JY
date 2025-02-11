test_alns_with_gaps.Rmd
================
Janet Young

2025-02-10

I was getting an error when all members of one population have a gap for
a whole codon. This showed up when Sage used my code to analyze the
mod(mdg4) gene.

For testing, I made a small example alignment, that has a gap in every
pop1 sequences at some positions (file =
`test_data/testData_for_troubleshooting/test_otherMiscProblemAlns/test_onePopAllGaps.fa`)

``` r
## example alignment WITHOUT gaps in all members of one population
alnFile <- here("test_data/MKTwebsiteExample/MKTwebsite_testAln.fa")
aln <- readBStringSet(alnFile)
names(aln) <- sapply(strsplit(names(aln), " "), "[[", 1)
alnSpecies <- gsub("\\d","",names(aln))
alnPopulations <- split(names(aln), alnSpecies)

## example alignment WITH gaps in all members of one population
alnWithGapsFile <- here("test_data/testData_for_troubleshooting/test_otherMiscProblemAlns/test_onePopAllGaps.fa")
alnWithGaps <- readBStringSet(alnWithGapsFile)
names(alnWithGaps) <- sapply(strsplit(names(alnWithGaps), " "), "[[", 1)
alnWithGapsSpecies <- gsub("\\d","",names(alnWithGaps))
alnWithGapsPopulations <- split(names(alnWithGaps), alnWithGapsSpecies)
```

I tracked down the error, to a call within `doMKtest()`, to the
`getCodonChangeCounts()` function (I’ve fixed it now). This code was the
minimal example to show the error:

``` r
getCodonChangeCounts("---","CTA")
```

    ##   codon before after ns s tot
    ## 1     1    ---   CTA  0 0   0

Run `doMKtest` on alignment without gaps:

``` r
MKresults_noGaps <- doMKtest(alnFile, 
                      pop1seqs=alnWithGapsPopulations[["pongo"]],
                      pop2seqs=alnWithGapsPopulations[["trachy"]],
                      quiet=TRUE)
```

    ## Warning in doMKtest(alnFile, pop1seqs = alnWithGapsPopulations[["pongo"]], : this alignment has a codon at the end with only stop codons - we will strip it out and not count any changes in this codon

Run `doMKtest` on alignment with gaps:

``` r
MKresults <- doMKtest(alnWithGapsFile, 
                      pop1seqs=alnWithGapsPopulations[["pongo"]],
                      pop2seqs=alnWithGapsPopulations[["trachy"]],
                      quiet=TRUE)
```

    ## Warning in doMKtest(alnWithGapsFile, pop1seqs = alnWithGapsPopulations[["pongo"]], : this alignment has a codon at the end with only stop codons - we will strip it out and not count any changes in this codon

    ## Warning: there are nucleotide positions where (in one population) every sequence has a gap or N. Positions:  16,17,18,28,29,30

    ## Warning: there are nucleotide positions where (in one population) every sequence has a gap or N. Positions:  40,41,42

Other than the `input` column (tells us the input alignment file name),
the results are identical. That’s good, because the codons where I
introduced artificial gaps did not contain any evolutionary changes.

``` r
identical(MKresults_noGaps[["summary"]] %>% dplyr::select(-input), 
          MKresults[["summary"]] %>% dplyr::select(-input))
```

    ## [1] TRUE

Show MK results (alignment with gaps):

``` r
MKresults[["summary"]] %>% 
    dplyr::select(-input) %>% 
    kable() %>% 
    kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
num_seqs
</th>
<th style="text-align:right;">
num_seqs_pop1
</th>
<th style="text-align:right;">
num_seqs_pop2
</th>
<th style="text-align:right;">
num_seqs_outgroup
</th>
<th style="text-align:left;">
seqs_not_used
</th>
<th style="text-align:right;">
length_NT
</th>
<th style="text-align:right;">
length_AA
</th>
<th style="text-align:right;">
region_start_coord_codons
</th>
<th style="text-align:right;">
region_end_coord_codons
</th>
<th style="text-align:left;">
filter_rare_alleles
</th>
<th style="text-align:right;">
rare_allele_freq_threshold
</th>
<th style="text-align:left;">
parameter_combining_approach
</th>
<th style="text-align:left;">
first_pop1_seq
</th>
<th style="text-align:left;">
first_pop2_seq
</th>
<th style="text-align:right;">
pop1_vs_pop2_chiSq_pVal
</th>
<th style="text-align:right;">
pop1_vs_pop2_FET_pVal
</th>
<th style="text-align:right;">
pop1_vs_pop2_Dn
</th>
<th style="text-align:right;">
pop1_vs_pop2_Ds
</th>
<th style="text-align:right;">
pop1_vs_pop2_Pn
</th>
<th style="text-align:right;">
pop1_vs_pop2_Ps
</th>
<th style="text-align:right;">
pop1_vs_pop2_NI
</th>
<th style="text-align:right;">
pop1_vs_pop2_alpha
</th>
<th style="text-align:left;">
pop1_vs_pop2_result
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
pop1_vs_pop2_Dn
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
</td>
<td style="text-align:right;">
396
</td>
<td style="text-align:right;">
132
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
132
</td>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
conservative
</td>
<td style="text-align:left;">
pongo1
</td>
<td style="text-align:left;">
trachy1
</td>
<td style="text-align:right;">
0.1184649
</td>
<td style="text-align:right;">
0.2307692
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
not signif
</td>
</tr>
</tbody>
</table>

Show positions table for one of the gap codons:

``` r
MKresults[["positions"]] %>% 
    # mutate is because kable turns hyphens into bullets
    # https://github.com/haozhu233/kableExtra/issues/223
    mutate(pop1_anc=gsub("-","--",pop1_anc)) %>% 
    mutate(pop2_anc=gsub("-","--",pop2_anc)) %>% 
    as_tibble() %>% 
    filter(pos %in% 16:18) %>% 
    kable() %>% 
    kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
pos
</th>
<th style="text-align:right;">
codon
</th>
<th style="text-align:right;">
codon_pos
</th>
<th style="text-align:left;">
pop1_anc
</th>
<th style="text-align:left;">
pop1_poly
</th>
<th style="text-align:right;">
pop1_num_alleles
</th>
<th style="text-align:left;">
pop2_anc
</th>
<th style="text-align:left;">
pop2_poly
</th>
<th style="text-align:right;">
pop2_num_alleles
</th>
<th style="text-align:left;">
fixed_difference
</th>
<th style="text-align:right;">
pop1_A
</th>
<th style="text-align:right;">
pop1_C
</th>
<th style="text-align:right;">
pop1_G
</th>
<th style="text-align:right;">
pop1_T
</th>
<th style="text-align:right;">
pop1\_-
</th>
<th style="text-align:right;">
pop1_N
</th>
<th style="text-align:left;">
pop1_major_allele
</th>
<th style="text-align:left;">
pop1_minor_alleles
</th>
<th style="text-align:left;">
pop1_major_allele_freqs
</th>
<th style="text-align:left;">
pop1_minor_allele_freqs
</th>
<th style="text-align:right;">
pop2_A
</th>
<th style="text-align:right;">
pop2_C
</th>
<th style="text-align:right;">
pop2_G
</th>
<th style="text-align:right;">
pop2_T
</th>
<th style="text-align:right;">
pop2\_-
</th>
<th style="text-align:right;">
pop2_N
</th>
<th style="text-align:left;">
pop2_major_allele
</th>
<th style="text-align:left;">
pop2_minor_alleles
</th>
<th style="text-align:left;">
pop2_major_allele_freqs
</th>
<th style="text-align:left;">
pop2_minor_allele_freqs
</th>
<th style="text-align:right;">
pop1_vs_pop2_Dn
</th>
<th style="text-align:right;">
pop1_vs_pop2_Ds
</th>
<th style="text-align:right;">
pop1_codon
</th>
<th style="text-align:left;">
pop1_codonPos
</th>
<th style="text-align:right;">
pop1_Pn
</th>
<th style="text-align:right;">
pop1_Ps
</th>
<th style="text-align:right;">
pop2_codon
</th>
<th style="text-align:left;">
pop2_codonPos
</th>
<th style="text-align:right;">
pop2_Pn
</th>
<th style="text-align:right;">
pop2_Ps
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
–
</td>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
C
</td>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">

-   </td>
    <td style="text-align:left;">
    NA
    </td>
    <td style="text-align:left;">
    1
    </td>
    <td style="text-align:left;">
    NA
    </td>
    <td style="text-align:right;">
    0
    </td>
    <td style="text-align:right;">
    2
    </td>
    <td style="text-align:right;">
    0
    </td>
    <td style="text-align:right;">
    0
    </td>
    <td style="text-align:right;">
    0
    </td>
    <td style="text-align:right;">
    0
    </td>
    <td style="text-align:left;">
    C
    </td>
    <td style="text-align:left;">
    NA
    </td>
    <td style="text-align:left;">
    1
    </td>
    <td style="text-align:left;">
    NA
    </td>
    <td style="text-align:right;">
    0
    </td>
    <td style="text-align:right;">
    0
    </td>
    <td style="text-align:right;">
    6
    </td>
    <td style="text-align:left;">
    1
    </td>
    <td style="text-align:right;">
    0
    </td>
    <td style="text-align:right;">
    0
    </td>
    <td style="text-align:right;">
    6
    </td>
    <td style="text-align:left;">
    1
    </td>
    <td style="text-align:right;">
    0
    </td>
    <td style="text-align:right;">
    0
    </td>
    </tr>
    <tr>
    <td style="text-align:right;">
    17
    </td>
    <td style="text-align:right;">
    6
    </td>
    <td style="text-align:right;">
    2
    </td>
    <td style="text-align:left;">
    –
    </td>
    <td style="text-align:left;">
    FALSE
    </td>
    <td style="text-align:right;">
    1
    </td>
    <td style="text-align:left;">
    T
    </td>
    <td style="text-align:left;">
    FALSE
    </td>
    <td style="text-align:right;">
    1
    </td>
    <td style="text-align:left;">
    TRUE
    </td>
    <td style="text-align:right;">
    0
    </td>
    <td style="text-align:right;">
    0
    </td>
    <td style="text-align:right;">
    0
    </td>
    <td style="text-align:right;">
    0
    </td>
    <td style="text-align:right;">
    2
    </td>
    <td style="text-align:right;">
    0
    </td>
    <td style="text-align:left;">

    -   </td>
        <td style="text-align:left;">
        NA
        </td>
        <td style="text-align:left;">
        1
        </td>
        <td style="text-align:left;">
        NA
        </td>
        <td style="text-align:right;">
        0
        </td>
        <td style="text-align:right;">
        0
        </td>
        <td style="text-align:right;">
        0
        </td>
        <td style="text-align:right;">
        2
        </td>
        <td style="text-align:right;">
        0
        </td>
        <td style="text-align:right;">
        0
        </td>
        <td style="text-align:left;">
        T
        </td>
        <td style="text-align:left;">
        NA
        </td>
        <td style="text-align:left;">
        1
        </td>
        <td style="text-align:left;">
        NA
        </td>
        <td style="text-align:right;">
        NA
        </td>
        <td style="text-align:right;">
        NA
        </td>
        <td style="text-align:right;">
        6
        </td>
        <td style="text-align:left;">
        2
        </td>
        <td style="text-align:right;">
        0
        </td>
        <td style="text-align:right;">
        0
        </td>
        <td style="text-align:right;">
        6
        </td>
        <td style="text-align:left;">
        2
        </td>
        <td style="text-align:right;">
        0
        </td>
        <td style="text-align:right;">
        0
        </td>
        </tr>
        <tr>
        <td style="text-align:right;">
        18
        </td>
        <td style="text-align:right;">
        6
        </td>
        <td style="text-align:right;">
        3
        </td>
        <td style="text-align:left;">
        –
        </td>
        <td style="text-align:left;">
        FALSE
        </td>
        <td style="text-align:right;">
        1
        </td>
        <td style="text-align:left;">
        A
        </td>
        <td style="text-align:left;">
        FALSE
        </td>
        <td style="text-align:right;">
        1
        </td>
        <td style="text-align:left;">
        TRUE
        </td>
        <td style="text-align:right;">
        0
        </td>
        <td style="text-align:right;">
        0
        </td>
        <td style="text-align:right;">
        0
        </td>
        <td style="text-align:right;">
        0
        </td>
        <td style="text-align:right;">
        2
        </td>
        <td style="text-align:right;">
        0
        </td>
        <td style="text-align:left;">

        -   </td>
            <td style="text-align:left;">
            NA
            </td>
            <td style="text-align:left;">
            1
            </td>
            <td style="text-align:left;">
            NA
            </td>
            <td style="text-align:right;">
            2
            </td>
            <td style="text-align:right;">
            0
            </td>
            <td style="text-align:right;">
            0
            </td>
            <td style="text-align:right;">
            0
            </td>
            <td style="text-align:right;">
            0
            </td>
            <td style="text-align:right;">
            0
            </td>
            <td style="text-align:left;">
            A
            </td>
            <td style="text-align:left;">
            NA
            </td>
            <td style="text-align:left;">
            1
            </td>
            <td style="text-align:left;">
            NA
            </td>
            <td style="text-align:right;">
            NA
            </td>
            <td style="text-align:right;">
            NA
            </td>
            <td style="text-align:right;">
            6
            </td>
            <td style="text-align:left;">
            3
            </td>
            <td style="text-align:right;">
            0
            </td>
            <td style="text-align:right;">
            0
            </td>
            <td style="text-align:right;">
            6
            </td>
            <td style="text-align:left;">
            3
            </td>
            <td style="text-align:right;">
            0
            </td>
            <td style="text-align:right;">
            0
            </td>
            </tr>
            </tbody>
            </table>

# Finished

``` r
sessionInfo()
```

    ## R version 4.4.0 (2024-04-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 18.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: FlexiBLAS OPENBLAS;  LAPACK version 3.11.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] Biostrings_2.72.0   GenomeInfoDb_1.40.0 XVector_0.44.0     
    ##  [4] IRanges_2.38.0      S4Vectors_0.42.0    BiocGenerics_0.50.0
    ##  [7] kableExtra_1.4.0    lubridate_1.9.3     forcats_1.0.0      
    ## [10] stringr_1.5.1       dplyr_1.1.4         purrr_1.0.2        
    ## [13] readr_2.1.5         tidyr_1.3.1         tibble_3.2.1       
    ## [16] ggplot2_3.5.1       tidyverse_2.0.0     here_1.0.1         
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] utf8_1.2.4              generics_0.1.3          xml2_1.3.6             
    ##  [4] stringi_1.8.4           hms_1.1.3               digest_0.6.35          
    ##  [7] magrittr_2.0.3          evaluate_0.23           grid_4.4.0             
    ## [10] timechange_0.3.0        fastmap_1.2.0           jsonlite_1.8.8         
    ## [13] rprojroot_2.0.4         httr_1.4.7              fansi_1.0.6            
    ## [16] UCSC.utils_1.0.0        viridisLite_0.4.2       scales_1.3.0           
    ## [19] cli_3.6.2               crayon_1.5.2            rlang_1.1.4            
    ## [22] munsell_0.5.1           withr_3.0.0             yaml_2.3.8             
    ## [25] tools_4.4.0             tzdb_0.4.0              colorspace_2.1-0       
    ## [28] GenomeInfoDbData_1.2.12 vctrs_0.6.5             R6_2.5.1               
    ## [31] lifecycle_1.0.4         zlibbioc_1.50.0         pkgconfig_2.0.3        
    ## [34] pillar_1.9.0            gtable_0.3.5            glue_1.7.0             
    ## [37] systemfonts_1.1.0       highr_0.10              xfun_0.44              
    ## [40] tidyselect_1.2.1        rstudioapi_0.16.0       knitr_1.46             
    ## [43] htmltools_0.5.8.1       rmarkdown_2.26          svglite_2.1.3          
    ## [46] compiler_4.4.0
