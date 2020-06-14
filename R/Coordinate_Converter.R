#' Coordinate Converter 1
#'
#' This function take genomic location,  nucleotide exchange as input and convert these to a dataframe which contains
#' HUGO Symbol, transcript ID, position in primary transcript, exchanged amino acid.
#'
#' @param chromosome sequence name like "chrX"
#' @param genomic_start genomic start position
#' @param genomic_end genomic end position
#' @param refAllele  reference allele
#' @param varAllele  variant allele
#' @param strand    sequence strand "+" or "-"
#'
#' @return dataframe which contians corresponding HUGO Symbol, transcript ID, position in primary transcript, exchanged amino acid.
#' @export
#'
#' @import GenomicFeatures
#' @import VariantAnnotation
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import biomaRt
#' @import rtracklayer
#' @import AnnotationDbi
#' @import Biostrings
#' @import gwascat
#' @import liftOver
#'
#'
#' @examples
Coordinate_Covnerter1<-function(chromosome,genomic_start,genomic_end,refAllele,varAllele,strand){

  ################
  ######Load Data
  if(!exists("gff_File")){
    GFFURI ="ftp://ftp.ensembl.org/pub/release-100/gff3/homo_sapiens/Homo_sapiens.GRCh38.100.gff3.gz"
    gff_File<<-download_GFF_fromURI(GFFURI)
    GFFDataload(gff_File)
  }
  else {
    GFFDataload(gff_File)

  }


  #Coordinate_input_Table <- read.table(Coordinate_input,sep="\t", header=TRUE)
  if(missing(strand)){
    Coordinate_Granges <- GRanges(seqnames = Rle(chromosome),
                                  ranges = IRanges(start = genomic_start, end = genomic_end),
                                  REF = refAllele)

  } else {
    Coordinate_Granges <- GRanges(seqnames = Rle(chromosome),
                                  ranges = IRanges(start = genomic_start, end = genomic_end),
                                  strand = strand,
                                  REF = refAllele)

  }

  seqlevelsStyle(txdb) <- "UCSC"
  seqlevelsStyle(Coordinate_Granges) <- "UCSC"

  VarAllele<-DNAStringSet(varAllele)

  coding <- predictCoding(Coordinate_Granges, txdb, Hsapiens,VarAllele)

  codingTB <- data.frame(coding)

  Transcript_Coding_TB<-Coordinate_to_codingTranscriptTB(codingTB,gff_DataTable)
  if(nrow(Transcript_Coding_TB)==0){

    DF<-dplyr::select(Transcript_Coding_TB, seqnames, start, end,GENEID,transcript_id,TXstart,TXend,CDSLOC.start)
    DF$REFAA <- character()
    DF$VARAA <- character()
    names(DF)[names(DF) == "GENEID"] <- "Hugo_Symbol"
    names(DF)[names(DF) == "REFAA"] <- "RefAminoAcid"
    names(DF)[names(DF) == "VARAA"] <- "VarAminoAcid"
    names(DF)[names(DF) == "CDSLOC.start"] <- "CDS_position"
    return(DF)
  } else{

    DF<-dplyr::select(Transcript_Coding_TB, seqnames, start, end,GENEID,transcript_id,TXstart,TXend,CDSLOC.start,REFAA,VARAA)
    names(DF)[names(DF) == "GENEID"] <- "Hugo_Symbol"
    names(DF)[names(DF) == "REFAA"] <- "RefAminoAcid"
    names(DF)[names(DF) == "VARAA"] <- "VarAminoAcid"
    names(DF)[names(DF) == "CDSLOC.start"] <- "CDS_position"
    return(DF)
  }



}



#' Coordinate Converter 2
#'
#'This function take genomic location,  nucleotide exchange as input and convert these to a dataframe which contains
#'HUGO Symbol, transcript ID, position exchanged amino acid.
#'
#'
#' @param chromosome    sequence name like "chrX"
#' @param genomic_start   genomic start position
#' @param genomic_end   genomic end position
#' @param refAllele   reference allele
#' @param varAllele   variant allele
#' @param strand     sequence strand "+" or "-"
#'
#' @return  dataframe which contians corresponding HUGO Symbol, transcript ID, position exchanged amino acid.
#' @export
#'
#' @examples
Coordinate_Covnerter2<-function(chromosome,genomic_start,genomic_end,refAllele,varAllele,strand){
  ################
  ######Load Data
  if(!exists("gff_File")){
    GFFURI ="ftp://ftp.ensembl.org/pub/release-100/gff3/homo_sapiens/Homo_sapiens.GRCh38.100.gff3.gz"
    gff_File<<-download_GFF_fromURI(GFFURI)
    GFFDataload(gff_File)
  }
  else {
    GFFDataload(gff_File)

  }

  #Coordinate_input_Table <- read.table(Coordinate_input,sep="\t", header=TRUE)
  if(missing(strand)){
    Coordinate_Granges <- GRanges(seqnames = Rle(chromosome),
                                  ranges = IRanges(start = genomic_start, end = genomic_end),
                                  REF = refAllele)

  } else {
    Coordinate_Granges <- GRanges(seqnames = Rle(chromosome),
                                  ranges = IRanges(start = genomic_start, end = genomic_end),
                                  strand = strand,
                                  REF = refAllele)

  }

  seqlevelsStyle(txdb) <- "UCSC"
  seqlevelsStyle(Coordinate_Granges) <- "UCSC"

  VarAllele<-DNAStringSet(varAllele)

  coding <- predictCoding(Coordinate_Granges, txdb, Hsapiens,VarAllele)

  codingTB <- data.frame(coding)

  Transcript_Coding_TB<-Coordinate_to_transcriptTB(codingTB,GFFCDSTable)

  DF<-dplyr::select(Transcript_Coding_TB, seqnames, start, end,GENEID,PROTEINLOC,transcript_id)
  names(DF)[names(DF) == "GENEID"] <- "Hugo_Symbol"

  return(DF)
}


#' Coordinate Converter 3
#'
#' This function take genomic location as input and convert these to a dataframe which contains
#' corresponding HUGO Symbol, list of all transcript ID.
#'
#'
#' @param chromosome
#' @param genomic_start
#' @param genomic_end
#'
#' @return  dataframe which contians corresponding HUGO Symbol, list of all transcript ID.
#' @export
#'
#' @examples
Coordinate_Covnerter3<-function(chromosome,genomic_start,genomic_end){
  ################
  ######Load Data
  if(!exists("gff_File")){
    GFFURI ="ftp://ftp.ensembl.org/pub/release-100/gff3/homo_sapiens/Homo_sapiens.GRCh38.100.gff3.gz"
    gff_File<<-download_GFF_fromURI(GFFURI)
    GFFDataload(gff_File)
  }
  else {
    GFFDataload(gff_File)

  }

  #Coordinate_Mapping_gffTB()
  tempTB<-Coordinate_Mapping_gffTB(chromosome,genomic_start,genomic_end,gff_DataTable)
  Hugo_Symbol <- tempTB[tempTB$type=="gene",]$Name

  tempTB<-tempTB[!is.na(tempTB$transcript_id),]
  DF<-dplyr::select(tempTB, seqnames, start, end,transcript_id,Name)
  names(DF)[names(DF) == "Name"] <- "Hugo_Symbol"
  DF$Hugo_Symbol<-Hugo_Symbol


  return(DF)
}

#' Coordinate Converter 4
#'
#' This function take genomic location as input and convert these to a dataframe which contains
#' corresponding HUGO symbol and exon number.
#'
#'
#' @param chromosome  sequence name like "chrX"
#' @param genomic_start   genomic start position
#' @param genomic_end   genomic end position
#'
#' @return dataframe which contians corresponding HUGO symbol, exon number.
#' @export
#'
#' @examples
Coordinate_Covnerter4<-function(chromosome,genomic_start,genomic_end){
  ################
  ######Load Data
  if(!exists("gff_File")){
    GFFURI ="ftp://ftp.ensembl.org/pub/release-100/gff3/homo_sapiens/Homo_sapiens.GRCh38.100.gff3.gz"
    gff_File<<-download_GFF_fromURI(GFFURI)
    GFFDataload(gff_File)
  }
  else {
    GFFDataload(gff_File)

  }

  #Coordinate_Mapping_gffTB()
  tempTB<-Coordinate_Mapping_gffTB(chromosome,genomic_start,genomic_end,gff_DataTable)
  Hugo_Symbol <- tempTB[tempTB$type=="gene",]$Name

  tempTB<-tempTB[!is.na(tempTB$rank),]
  DF<-dplyr::select(tempTB, seqnames, start, end,exon_id,rank)

  #DF$Hugo_Symbol<-Hugo_Symbol


  return(DF)
}

#' Coordinate Converter 5_1
#'
#' This function take transcript id and transcript location as input and convert these to a dataframe which contains
#' corresponding genomic position, codon and codon's number.
#'
#' @param transcript_id   transcript passed in should be external transcript name like "HDAC6-201"
#' @param TX_start  transcript start position like 2657
#' @param TX_end    transcript end position like 2657
#'
#' @return  dataframe which contians corresponding genomic position, codon and codon's number.
#' @export
#'
#' @examples
Coordinate_Covnerter5_1<-function(transcript_id,TX_start,TX_end){
  ################
  ######Load Data
  if(!exists("gff_File")){
    GFFURI ="ftp://ftp.ensembl.org/pub/release-100/gff3/homo_sapiens/Homo_sapiens.GRCh38.100.gff3.gz"
    gff_File<<-download_GFF_fromURI(GFFURI)
    GFFDataload(gff_File)
  }
  else {
    GFFDataload(gff_File)

  }

  #transcript_id <- "HDAC6-203"
  #TX_start=2657
  #TX_end=2657



  #genomeTX <-transcripts(txdb)
  #names (genomeTX) <- id2name (txdb, "tx")
  #levelsActivate <- seqlevels(txdb)
  #levelsActivate<-head(levelsActivate,25)
  #seqlevels(txdb) <- levelsActivate
  seqlevelsStyle(txdb) <- "UCSC"
    cdsTX<-cdsBy(txdb, "tx", use.names=TRUE)[transcript_id]

    tx_coord <- GRanges(transcript_id, IRanges(TX_start, TX_end))
    genomeCoord <- mapFromTranscripts(tx_coord, cdsTX)
    genomeCoord <- granges(genomeCoord)
    genomicCoordTB<-data.frame(genomeCoord)



    ####convert genomic coordinate to codon and codon number
    DF<-Coordinate_to_codons(as.character(genomicCoordTB$seqnames),genomicCoordTB$start,genomicCoordTB$end,GFFCDSTable)


    return(DF)
}

#' Coordinate Converter 5_2
#'
#' This function take transcript id and transcript location as input and convert these to a dataframe which contains
#' corresponding genomic position, codon and codon's number.
#'
#'
#' @param transcript_id   transcript id passed in should be ensembl transcript id like "ENST00000334136"
#' @param TX_position   transcript position like 2657
#'
#' @return  dataframe which contians corresponding genomic position, codon and codon's number.
#' @export
#'
#' @examples
Coordinate_Covnerter5_2<-function(transcript_id,TX_position){
  ################
  ######Load Data
  if(!exists("gff_File")){
    GFFURI ="ftp://ftp.ensembl.org/pub/release-100/gff3/homo_sapiens/Homo_sapiens.GRCh38.100.gff3.gz"
    gff_File<<-download_GFF_fromURI(GFFURI)
    GFFDataload(gff_File)
  }
  else {
    GFFDataload(gff_File)

  }

  #transcript_id="ENST00000486554"
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  #grch37
  ###mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice",dataset="hsapiens_gene_ensembl")


  results <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id","hgnc_symbol","transcript_start","external_transcript_name"),
                   filters = c("ensembl_transcript_id"),
                   values = transcript_id,
                   mart = mart)
  transcript_id_name <- results$external_transcript_name

  DF<-Coordinate_Covnerter5_1(transcript_id_name,TX_position,TX_position)

  #transcript_id <- "HDAC6-203"
  #TX_start=2657
  #TX_end=2657
  #TX_position=2657
  #transcript_id<- "ENST00000334136"
  #proteinLoc<-ceiling(TX_position/3)

  #GFFCDSTable$transcript_id<-gsub(pattern = "^.*:",replacement = "",GFFCDSTable$Parent)
  #targetTB<-GFFCDSTable[GFFCDSTable$transcript_id==transcript_id,]
  #targetTB$index <-targetTB$start+2657

  #***genomicCoordTB<-gff_DataTable[!is.na(gff_DataTable$transcript_id),]
  #***genomicCoordTB<-genomicCoordTB[genomicCoordTB$transcript_id == transcript_id,]
  #***TXtoGenomicindex <- genomicCoordTB$start + TX_position
  #TXtoGenomicindex<-48823056


  ####convert genomic coordinate to codon and codon number
  #***DF<-Coordinate_to_codons(as.character(genomicCoordTB$seqnames),TXtoGenomicindex,TXtoGenomicindex,GFFCDSTable)



  return(DF)

}

#' Coordinate Converter 6
#'
#' This function take genomic location, orginal and target genome version as input and convert these to a dataframe which contains
#' corresponding genomic position.
#'
#' @param chromosome    sequence name like "chrX"
#' @param genomic_start    genomic start position
#' @param genomic_end    genomic end position
#' @param genomeVersion   original genome version
#' @param targetGenomeVersion  target genome version which want to transfer
#'
#' @return  dataframe which contians corresponding converted target genome coordinate.
#' @export
#'
#'
#' @examples
Coordinate_Covnerter6<-function(chromosome,genomic_start,genomic_end,genomeVersion,targetGenomeVersion){

  #genomeVersion="hg38"
  #targetGenomeVersion="hg19"
  Coordinate_Granges <- GRanges(seqnames = Rle(chromosome),
                                ranges = IRanges(start = genomic_start, end = genomic_end))
  genome(Coordinate_Granges) = genomeVersion

  tempTarget<-paste(toupper(substr(targetGenomeVersion, 1, 1)), substr(targetGenomeVersion, 2, nchar(targetGenomeVersion)), sep="")
  ChainFileName <- paste(genomeVersion,"To",tempTarget,".over.chain",sep = "")

  path = system.file(package="liftOver", "extdata")
  ChainFilepath = paste(path,"/",ChainFileName,sep = "")
  if(file.exists(ChainFilepath)){ch = import.chain(ChainFilepath)}
  else if(genomeVersion=="hg16"){
    if(targetGenomeVersion=="hg17"){
      ChainURI = "https://hgdownload.soe.ucsc.edu/goldenPath/hg16/liftOver/hg16ToHg17.over.chain.gz"

      ChainGZfile = paste(path,"/",basename(ChainURI),sep = "")
      download.file(ChainURI,ChainGZfile)
      gunzip(ChainGZfile)
      ChainFilepath <- tools::file_path_sans_ext(ChainGZfile)

    } else if(targetGenomeVersion=="hg18"){
      ChainURI = "https://hgdownload.soe.ucsc.edu/goldenPath/hg16/liftOver/hg16ToHg18.over.chain.gz"
      ChainGZfile =  paste(path,"/",basename(ChainURI),sep = "")
      download.file(ChainURI,ChainGZfile)
      gunzip(ChainGZfile)
      ChainFilepath <- tools::file_path_sans_ext(ChainGZfile)

    }else if(targetGenomeVersion=="hg19"){
      ChainURI = "https://hgdownload.soe.ucsc.edu/goldenPath/hg16/liftOver/hg16ToHg19.over.chain.gz"
      ChainGZfile =  paste(path,"/",basename(ChainURI),sep = "")
      download.file(ChainURI,ChainGZfile)
      gunzip(ChainGZfile)
      ChainFilepath <- tools::file_path_sans_ext(ChainGZfile)

    }else if(targetGenomeVersion=="hg38"){
      ChainURI = "https://hgdownload.soe.ucsc.edu/goldenPath/hg16/liftOver/hg16ToHg38.over.chain.gz"
      ChainGZfile =  paste(path,"/",basename(ChainURI),sep = "")
      download.file(ChainURI,ChainGZfile)
      gunzip(ChainGZfile)
      ChainFilepath <- tools::file_path_sans_ext(ChainGZfile)

    } else {print("chain file to target version missing,please download target file first ")
      return(NULL)}
  }
  else if(genomeVersion=="hg17"){
    if(targetGenomeVersion=="hg16"){
      ChainURI = "https://hgdownload.soe.ucsc.edu/goldenPath/hg17/liftOver/hg17ToHg16.over.chain.gz"

      ChainGZfile = paste(path,"/",basename(ChainURI),sep = "")
      download.file(ChainURI,ChainGZfile)
      gunzip(ChainGZfile)
      ChainFilepath <- tools::file_path_sans_ext(ChainGZfile)

    } else if(targetGenomeVersion=="hg18"){
      ChainURI = "https://hgdownload.soe.ucsc.edu/goldenPath/hg17/liftOver/hg17ToHg18.over.chain.gz"
      ChainGZfile =  paste(path,"/",basename(ChainURI),sep = "")
      download.file(ChainURI,ChainGZfile)
      gunzip(ChainGZfile)
      ChainFilepath <- tools::file_path_sans_ext(ChainGZfile)

    }else if(targetGenomeVersion=="hg19"){
      ChainURI = "https://hgdownload.soe.ucsc.edu/goldenPath/hg17/liftOver/hg17ToHg19.over.chain.gz"
      ChainGZfile =  paste(path,"/",basename(ChainURI),sep = "")
      download.file(ChainURI,ChainGZfile)
      gunzip(ChainGZfile)
      ChainFilepath <- tools::file_path_sans_ext(ChainGZfile)

    }else if(targetGenomeVersion=="hg38"){
      ChainURI = "https://hgdownload.soe.ucsc.edu/goldenPath/hg17/liftOver/hg17ToHg38.over.chain.gz"
      ChainGZfile =  paste(path,"/",basename(ChainURI),sep = "")
      download.file(ChainURI,ChainGZfile)
      gunzip(ChainGZfile)
      ChainFilepath <- tools::file_path_sans_ext(ChainGZfile)

    } else {print("chain file to target version missing,please download target file first ")
      return(NULL)}
  }
  else if(genomeVersion=="hg18"){
    if(targetGenomeVersion=="hg19"){
      ChainURI = "http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg19.over.chain.gz"

      ChainGZfile = paste(path,"/",basename(ChainURI),sep = "")
      download.file(ChainURI,ChainGZfile)
      gunzip(ChainGZfile)
      ChainFilepath <- tools::file_path_sans_ext(ChainGZfile)

    } else if(targetGenomeVersion=="hg17"){
      ChainURI = "http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg17.over.chain.gz"
      ChainGZfile =  paste(path,"/",basename(ChainURI),sep = "")
      download.file(ChainURI,ChainGZfile)
      gunzip(ChainGZfile)
      ChainFilepath <- tools::file_path_sans_ext(ChainGZfile)

    }else if(targetGenomeVersion=="hg38"){
      ChainURI = "http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg38.over.chain.gz"
      ChainGZfile = paste(path,"/",basename(ChainURI),sep = "")
      download.file(ChainURI,ChainGZfile)
      gunzip(ChainGZfile)
      ChainFilepath <- tools::file_path_sans_ext(ChainGZfile)
    } else {print("chain file to target version missing,please download target file first ")
      return(NULL)}

  }
  else if(genomeVersion=="hg19"){
    if(targetGenomeVersion=="hg18"){
      ChainURI = "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg18.over.chain.gz"
      ChainGZfile =  paste(path,"/",basename(ChainURI),sep = "")
      download.file(ChainURI,ChainGZfile)
      gunzip(ChainGZfile)
      ChainFilepath <- tools::file_path_sans_ext(ChainGZfile)

    } else if(targetGenomeVersion=="hg17"){
      ChainURI = "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg17.over.chain.gz"
      ChainGZfile = paste(path,"/", basename(ChainURI),sep = "")
      download.file(ChainURI,ChainGZfile)
      gunzip(ChainGZfile)
      ChainFilepath <- tools::file_path_sans_ext(ChainGZfile)

    }else if(targetGenomeVersion=="hg38"){
      ChainURI = "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz"
      ChainGZfile = paste(path,"/", basename(ChainURI),sep = "")
      download.file(ChainURI,ChainGZfile)
      gunzip(ChainGZfile)
      ChainFilepath <- tools::file_path_sans_ext(ChainGZfile)
    } else {print("chain file to target version missing,please download target file first ")
      return(NULL)}
  }
  else if(genomeVersion=="hg38"){
    if(targetGenomeVersion=="hg19"){
      ChainURI = "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz"
      ChainGZfile = paste(path,"/",  basename(ChainURI),sep = "")
      download.file(ChainURI,ChainGZfile)
      gunzip(ChainGZfile)
      ChainFilepath <- tools::file_path_sans_ext(ChainGZfile)
    } else {print("chain file to target version missing,please download target file first ")
      return(NULL)}
  } else {
    print("Genome conversion chain file missing,please download target file first ")
    return(NULL)
  }

  ch = import.chain(ChainFilepath)

  seqlevelsStyle(Coordinate_Granges) = "UCSC"  # necessary
  liftOver_Coordinate_GR = liftOver(Coordinate_Granges, ch)
  liftOver_Coordinate_GR = unlist(liftOver_Coordinate_GR)

  genome(liftOver_Coordinate_GR) = targetGenomeVersion
  DF <-data.frame(liftOver_Coordinate_GR)


  #Coordinate_Mapping_gffTB()
  #DF$Hugo_Symbol<-Hugo_Symbol


  return(DF)
}
