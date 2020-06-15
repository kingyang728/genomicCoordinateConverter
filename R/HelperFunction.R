#' Download gff file from URI
#'
#' @param GFFURI the gff file URI link
#'
#' @return gff file downloaded path.
#' @export
#'
#' @examples
#' gff_File<-download_GFF_fromURI("ftp://ftp.ensembl.org/pub/release-100/gff3/homo_sapiens/Homo_sapiens.GRCh38.100.gff3.gz")
download_GFF_fromURI<-function(GFFURI){

  path = system.file(package="genomicCoordinateConverter", "extdata")
  if(path==""){
    path="."
  }
  gff_File <<- paste(path,"/",basename(GFFURI),sep = "")
  #gfffilename <- basename(GFFURI)
  if(!file.exists(gff_File)){
    download.file(GFFURI,gff_File)
  }

  return(gff_File)
}

#' Generate TxDB, gff dataframe, cds_gff dataframe based on gff file
#'
#' @param gff_File gff file located path
#'
#' @return
#' @export
#'
#'
GFFDataload <- function(gff_File){

  if(exists("txDBPath") & file.exists("txDBPath")  & !exists("txdb")){
    txdb <<- AnnotationDbi::loadDb(txDBPath)
    seqlevelsStyle(txdb) <- "UCSC"
  } else if(exists("txDBPath") & !file.exists("txDBPath")  & !exists("txdb")){
    txdb <<- makeTxDbFromGFF(gff_File,organism="Homo sapiens")
    seqlevelsStyle(txdb) <- "UCSC"
    txdbName <- tools::file_path_sans_ext(gff_File)
    txDBPath <<-paste(txdbName,"txdb.sqlite",sep = "_")
    saveDb(txdb, file=txDBPath)

  } else if(! exists("txDBPath") & ! exists("txdb")){
    txdb <<- makeTxDbFromGFF(gff_File,organism="Homo sapiens")
    seqlevelsStyle(txdb) <- "UCSC"
    txdbName <- tools::file_path_sans_ext(gff_File)
    txDBPath <<-paste(txdbName,"txdb.sqlite",sep = "_")
    saveDb(txdb, file=txDBPath)
  }

  if(! exists("gff_DataTable")){
    gff_Granges <- import(gff_File)
    seqlevelsStyle(gff_Granges) <- "UCSC"
    gff_DataTable<<-data.frame(gff_Granges)
  }
  if(exists("gff_DataTable") & ! exists("GFFCDSTable")){
    GFFCDSTable <<- gff_DataTable[gff_DataTable$type=="CDS",]
  }
}

#' Coordinate select function for gff dataframe Table.
#'
#' @param chr sequence name like "chrX"
#' @param start_pos genomic start position
#' @param end_pos genomic end position
#' @param gff_Table gff data frame generated from gff file
#'
#' @return selected dataframe according to chr, genomic start and end coordinate.
#' @export
#'
#'
Coordinate_Mapping_gffTB <- function(chr,start_pos,end_pos,gff_Table){
  #start_pos<-78755835
  #end_pos<-78755835
  #gff_Table<-XCDSTable
  #chr="X"
  subsetgffTable<-gff_Table[which(gff_Table$seqnames==chr & gff_Table$start<=start_pos & gff_Table$end>=end_pos),]
  #subsetgffTable<-subset(gff_Table, seqnames==chr&start <= start_pos & end >= end_pos)
  #subsetgffTable$start<-start_pos
  #subsetgffTable$end<-end_pos
  #temp<-geneMapping("X",78755835,78755835,XCDSTable)

  return(subsetgffTable)
}


#' Coordinate coding table to transcript coding table
#'
#' @param CodingTable AA coding table generated from predictCoding()
#' @param GFFCDSTable cds_gff table extracted from gff table.
#'
#' @return coding table with corresponding transcript id
#' @export
#'
#'
Coordinate_to_transcriptTB<-function(CodingTable,GFFCDSTable){

  GFFCDSTable$transcript_id<-gsub(pattern = "^.*:",replacement = "",GFFCDSTable$Parent)
  #CodingTable<-codingTB
  seqnames<-CodingTable$seqnames
  start<-CodingTable$start
  end <- CodingTable$end
  seqRegionTable<-data.frame(seqnames,start,end)
  seqRegionTable <-unique(seqRegionTable)

  #GFFCDSTable <- GFFCDSTable

  TranscriptID<-c()
  for(row in 1:nrow(seqRegionTable)){

    tempTable <-Coordinate_Mapping_gffTB(as.character(seqRegionTable[row,1]),seqRegionTable[row,2],seqRegionTable[row,3],GFFCDSTable)
    #tempTable<-select(tempTable, seqnames, start, end,Transcript_id)
    print(as.character(tempTable$transcript_id))
    TranscriptID<-append(TranscriptID,as.character(tempTable$transcript_id))
  }
  CodingTable$transcript_id<-TranscriptID
  return(CodingTable)
}


#' Coordinate coding table to full transcript coding table
#'
#' @param CodingTable AA coding table generated from predictCoding()
#' @param gff_DataTable gff table generated from gff file
#'
#' @return full merged transcript-coding table
#' @export
#'
#'
Coordinate_to_codingTranscriptTB<-function(CodingTable,gff_DataTable){
  #CodingTable<-codingTB
  seqnames<-CodingTable$seqnames
  start<-CodingTable$start
  end <- CodingTable$end
  seqRegionTable<-data.frame(seqnames,start,end)
  seqRegionTable <-unique(seqRegionTable)
  #GFFCDSTable <- GFFCDSTable

  TranscriptID<-c()
  PremRNATB<-gff_DataTable[FALSE,]   # create empty data frame .
  for(row in 1:nrow(seqRegionTable)){

    tempTable <-Coordinate_Mapping_gffTB(as.character(seqRegionTable[row,1]),seqRegionTable[row,2],seqRegionTable[row,3],gff_DataTable)
    CDCtempTB <- tempTable[tempTable$type=="CDS",]
    transcript_id<-gsub(pattern = "^.*:",replacement = "",CDCtempTB$Parent)
    mRNAtempTB<-tempTable[tempTable$type=="mRNA",]
    preMRNAtempTB<-mRNAtempTB[mRNAtempTB$transcript_id  %in% transcript_id,]

    #tempTable<-select(tempTable, seqnames, start, end,Transcript_id)
    print(transcript_id)
    TranscriptID<-append(TranscriptID,transcript_id)
    PremRNATB<-rbind(PremRNATB,preMRNAtempTB)
  }
  CodingTable$transcript_id<-TranscriptID

  PremRNATB$seqnames=NULL
  names(PremRNATB)[names(PremRNATB) == "start"] <- "TXstart"
  names(PremRNATB)[names(PremRNATB) == "end"] <- "TXend"
  names(PremRNATB)[names(PremRNATB) == "width"] <- "TXwidth"
  DF<-merge(CodingTable,PremRNATB,by="transcript_id")

  #GFFCDSTable <- GFFCDSTable
  return(DF)
}



#' Coordinate to codons table
#'
#' @param chr  sequence name like "chrX"
#' @param start_pos  genomic start position
#' @param end_pos  genomic end position
#' @param GFFCDSTable  cds_gff table extracted from gff table.
#'
#' @return table with coorespoding codons and codons number
#' @export
#'
#'
Coordinate_to_codons<-function(chr,start_pos,end_pos,GFFCDSTable){
  #chr = genomicCoordTB$seqnames
  #start_pos=genomicCoordTB$start
  #end_pos=genomicCoordTB$end

  GFFCDSTable$transcript_id<-gsub(pattern = "^.*:",replacement = "",GFFCDSTable$Parent)

  tempTB<-Coordinate_Mapping_gffTB(as.character(chr),start_pos,end_pos,GFFCDSTable)

  #tempTB<-Coordinate_Mapping_gffTB("chrX",48823056,48823056,GFFCDSTable )
  tempTB<-dplyr::select(tempTB,seqnames,start,end,width,strand,phase)
  tempTB<-unique(tempTB)
  tempTB$width <- tempTB$width-tempTB$phase
  if(length(tempTB$strand)==0){
    index = NULL

  } else if(tempTB$strand=="-"){
    tempTB$end <- tempTB$end - tempTB$phase

    index = tempTB$end -start_pos

  }else if(tempTB$strand=="+" | tempTB$strand=="*"){
    tempTB$start <- tempTB$start+tempTB$phase

    index = start_pos-tempTB$start

  }


  remaider = index %% 3
  division = index %/% 3

  cdsTargetRegion<-makeGRangesListFromDataFrame(tempTB,split.field = "width")
  cds_seqs<-extractTranscriptSeqs(Hsapiens, cdsTargetRegion)

  i=1
  while(i <= length(cds_seqs)){
    targetCodon <- as.character(codons(cds_seqs[[i]])[division+1])
    i=i+1
  }
  #targetCodon <- as.character(codons(cds_seqs[[1]])[division+1])
  if(!exists("targetCodon")){
    targetCodon<- ""
  }
  if(length(remaider)==0){
    remaider<-""
  }

  DF<-data.frame(chr,start_pos,end_pos)
  if(nrow(DF)==0){
    DF$condon= character()
    DF$condonNumber=integer()

    colnames(DF) <- c("seqnames", "start","end","condon","condonNumber")
    return(DF)
  } else{
    DF$condon <- targetCodon
    DF$condonNumber <- remaider+1
    colnames(DF) <- c("seqnames", "start","end","condon","condonNumber")
    return(DF)
  }

}

