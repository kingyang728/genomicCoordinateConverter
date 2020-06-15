
<!-- README.md is generated from README.Rmd. Please edit that file -->

# genomicCoordinateConverter

<!-- badges: start -->

<!-- badges: end -->

The goal of genomicCoordinateConverter is to provide functionalities
supporting convert coordinate to amino acid information.The package
contains 6 main converter functionalities: <br/> ‘*’ genomic location
(default genome: HG38), nucleotide exchange -\> HUGO Symbol, transcript
ID, position in primary transcript, exchanged amino acid <br/> ’*’
genomic location (default genome: HG38), nucleotide exchange -\> HUGO
Symbol, for all transcripts :transcript ID, position exchanged amino
acid <br/> ‘*’ genomic location (default genome: HG38) -\> HUGO Symbol,
list of all transcript ID <br/> ’*’ genomic location (default genome:
HG38), nucleotide exchange -\> HUGO symbol, exon number <br/> ‘*’
transcript ID, position , target genome -\> genomic location, triplet
and number in triplet <br/> ’*’ Genomeversion, genomic location, taget
genomeversion -\> genomic location <br/>

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kingyang728/genomicCoordinateConverter")
```

## Example

This is a basic example which shows you how to convert genomic
coordinate to target dataframe:

``` r
library(genomicCoordinateConverter)
#> Warning: replacing previous import 'R.utils::header' by
#> 'VariantAnnotation::header' when loading 'genomicCoordinateConverter'
#> 
#> 
## basic example code
#genomic location,  nucleotide exchange -> HUGO Symbol, transcript ID, position in primary transcript, exchanged amino acid
Coordinate_Covnerter1("chrX",48823056,48823056,"G","C")
#> Import genomic features from the file as a GRanges object ...
#> OK
#> Prepare the 'metadata' data frame ... OK
#> Make the TxDb object ...
#> Warning in .extract_exons_from_GRanges(exon_IDX, gr, mcols0, tx_IDX, feature = "exon", : 1922 exons couldn't be linked to a transcript so were dropped (showing
#>   only the first 6):
#>   seqid   start     end strand   ID            Name                     Parent
#> 1     1 2566410 2566605      + <NA> ENSE00001889931 transcript:ENST00000456687
#> 2     1 2567866 2568054      + <NA> ENSE00001659318 transcript:ENST00000456687
#> 3     1 2568977 2569888      + <NA> ENSE00001624410 transcript:ENST00000456687
#> 4     1 3205988 3208664      + <NA> ENSE00003757578 transcript:ENST00000624175
#> 5     1 9826289 9828271      - <NA> ENSE00003807109 transcript:ENST00000639753
#> 6     1 9827610 9828258      - <NA> ENSE00003759368 transcript:ENST00000623863
#>   Parent_type
#> 1        <NA>
#> 2        <NA>
#> 3        <NA>
#> 4        <NA>
#> 5        <NA>
#> 6        <NA>
#> Warning in .extract_exons_from_GRanges(cds_IDX, gr, mcols0, tx_IDX, feature = "cds", : 752 CDS couldn't be linked to a transcript so were dropped (showing
#>   only the first 6):
#>   seqid    start      end strand                  ID            Name
#> 1    14 21621904 21621946      + CDS:ENSP00000446309 ENSP00000446309
#> 2    14 21622285 21622567      + CDS:ENSP00000446309 ENSP00000446309
#> 3    14 21642973 21643015      + CDS:ENSP00000439668 ENSP00000439668
#> 4    14 21643304 21643578      + CDS:ENSP00000439668 ENSP00000439668
#> 5    14 21712331 21712394      + CDS:ENSP00000438195 ENSP00000438195
#> 6    14 21712570 21712843      + CDS:ENSP00000438195 ENSP00000438195
#>                       Parent Parent_type
#> 1 transcript:ENST00000542354        <NA>
#> 2 transcript:ENST00000542354        <NA>
#> 3 transcript:ENST00000390423        <NA>
#> 4 transcript:ENST00000390423        <NA>
#> 5 transcript:ENST00000390424        <NA>
#> 6 transcript:ENST00000390424        <NA>
#> OK
#> Warning in valid.GenomicRanges.seqinfo(x, suggest.trim = TRUE): GRanges object contains 25 out-of-bound ranges located on sequences
#>   219562, 219565, 219570, 219571, and 219573. Note that ranges located on
#>   a sequence whose length is unknown (NA) or on a circular sequence are
#>   not considered out-of-bound (use seqlengths() and isCircular() to get
#>   the lengths and circularity flags of the underlying sequences). You can
#>   use trim() to trim these ranges. See ?`trim,GenomicRanges-method` for
#>   more information.
#> [1] "ENST00000376619" "ENST00000334136" "ENST00000643374" "ENST00000644068"
#> [5] "ENST00000643934"
#>   seqnames    start      end Hugo_Symbol   transcript_id  TXstart    TXend
#> 1     chrX 48823056 48823056       HDAC6 ENST00000334136 48802067 48824976
#> 2     chrX 48823056 48823056       HDAC6 ENST00000376619 48802034 48824976
#> 3     chrX 48823056 48823056       HDAC6 ENST00000643374 48802159 48824956
#> 4     chrX 48823056 48823056       HDAC6 ENST00000643934 48802187 48824958
#> 5     chrX 48823056 48823056       HDAC6 ENST00000644068 48802168 48824976
#>   CDS_position RefAminoAcid VarAminoAcid
#> 1         2657            S            S
#> 2         2657            S            S
#> 3         2657            S            S
#> 4         2492            S            S
#> 5         2657            S            S
#genomic location, nucleotide exchange -> HUGO Symbol, for all transcripts :transcript ID, position exchanged amino acid
Coordinate_Covnerter2("chrX",48823056,48823056,"G","C")
#> Warning in valid.GenomicRanges.seqinfo(x, suggest.trim = TRUE): GRanges object contains 25 out-of-bound ranges located on sequences
#>   219562, 219565, 219570, 219571, and 219573. Note that ranges located on
#>   a sequence whose length is unknown (NA) or on a circular sequence are
#>   not considered out-of-bound (use seqlengths() and isCircular() to get
#>   the lengths and circularity flags of the underlying sequences). You can
#>   use trim() to trim these ranges. See ?`trim,GenomicRanges-method` for
#>   more information.
#> [1] "ENST00000376619" "ENST00000334136" "ENST00000643374" "ENST00000644068"
#> [5] "ENST00000643934"
#>   seqnames    start      end Hugo_Symbol PROTEINLOC   transcript_id
#> 1     chrX 48823056 48823056       HDAC6        886 ENST00000376619
#> 2     chrX 48823056 48823056       HDAC6        886 ENST00000334136
#> 3     chrX 48823056 48823056       HDAC6        886 ENST00000643374
#> 4     chrX 48823056 48823056       HDAC6        886 ENST00000644068
#> 5     chrX 48823056 48823056       HDAC6        831 ENST00000643934
#genomic location -> HUGO Symbol, list of all transcript ID
Coordinate_Covnerter3("chrX",48823056,48823056)
#>         seqnames    start      end   transcript_id Hugo_Symbol
#> 2851349     chrX 48802034 48824976 ENST00000376619       HDAC6
#> 2851410     chrX 48802055 48824982 ENST00000477528       HDAC6
#> 2851454     chrX 48802067 48824976 ENST00000334136       HDAC6
#> 2851555     chrX 48802133 48824860 ENST00000643659       HDAC6
#> 2851585     chrX 48802159 48824956 ENST00000643374       HDAC6
#> 2851646     chrX 48802168 48824976 ENST00000644068       HDAC6
#> 2851722     chrX 48802187 48824958 ENST00000643934       HDAC6
#> 2851831     chrX 48802888 48824967 ENST00000647398       HDAC6
#> 2851907     chrX 48813353 48824976 ENST00000486227       HDAC6
#> 2851953     chrX 48816470 48824947 ENST00000645643       HDAC6
#> 2851998     chrX 48818071 48823109 ENST00000480525       HDAC6
#> 2852006     chrX 48818313 48823579 ENST00000498808       HDAC6
#> 2852012     chrX 48819447 48824976 ENST00000488543       HDAC6
#genomic location (default genome: HG38), nucleotide exchange ->  HUGO symbol, exon number
Coordinate_Covnerter4("chrX",48823056,48823056)
#>         seqnames    start      end         exon_id rank
#> 2851399     chrX 48822912 48823588 ENSE00003465573   25
#> 2851434     chrX 48822912 48823588 ENSE00003634146   24
#> 2851504     chrX 48822912 48823588 ENSE00003465573   25
#> 2851580     chrX 48822912 48823588 ENSE00003634146   25
#> 2851635     chrX 48822912 48823588 ENSE00003465573   25
#> 2851696     chrX 48822912 48823588 ENSE00003465573   25
#> 2851770     chrX 48822912 48823588 ENSE00003465573   24
#> 2851920     chrX 48822912 48823588 ENSE00003634146   13
#> 2851965     chrX 48822912 48823588 ENSE00003634146    6
#> 2852004     chrX 48822912 48823109 ENSE00001839278    3
#> 2852015     chrX 48822912 48823588 ENSE00003634146    3
#transcript ID, position , target genome -> genomic location, triplet and number in triplet
##transcript id passed in should be ensembl transcript id like "ENST00000334136"
Coordinate_Covnerter5_2("ENST00000334136",2657)
#> Warning: `select_()` is deprecated as of dplyr 0.7.0.
#> Please use `select()` instead.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_warnings()` to see where this warning was generated.
#> Warning: `filter_()` is deprecated as of dplyr 0.7.0.
#> Please use `filter()` instead.
#> See vignette('programming') for more help
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_warnings()` to see where this warning was generated.
#> Cache found
#>   seqnames    start      end condon condonNumber
#> 1     chrX 48823056 48823056    TCA            2
#Genomeversion, genomic location, taget genomeversion -> genomic location
Coordinate_Covnerter6("chrX",45060024,45060024,"hg16","hg18")
#>   seqnames    start      end width strand
#> 1     chrX 45898986 45898986     1      *
```
