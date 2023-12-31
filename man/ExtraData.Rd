% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExtraData.R
\docType{data}
\name{ExtraData}
\alias{ExtraData}
\alias{ARIscores}
\alias{BioCartagenesetsID}
\alias{BiocartaGenesetsSymbol}
\alias{GOgenesetsID}
\alias{GOgenesetsSymbol}
\alias{GSVAResults}
\alias{HallmarkgenesetsID}
\alias{HallmarkgenesetsSymbol}
\alias{KEGGgenesetsID}
\alias{KEGGgenesetsSymbol}
\alias{KICH}
\alias{KICHcli}
\alias{KICHfilter}
\alias{ReactomegenesetsID}
\alias{ReactomegenesetsSymbol}
\alias{WikiPathwaysgenesetsID}
\alias{WikipathwaysSymbol}
\title{Data about PathActMarker}
\format{
A 2-column dataframe, the first column being the method names and the second column being the ARI scores.

BioCartagenesetsID:pathway names; list:292 lists of genes (gene id) in each pathway.

BiocartaGenesetsSymbol:pathway names; list:292 lists of genes (gene symbol) in each pathway.

GOgenesetsID:pathway names; list:10561 lists of genes (gene id) in each pathway.

GOgenesetsID:pathway names; list:10561 lists of genes (gene symbol) in each pathway.

The n x m pathway activity matrix, where n is the pathways and m the samples.

HallmarkgenesetsID:pathway names; list:50 lists of genes (gene id) in each pathway.

HallmarkgenesetsSymbol:pathway names; list:50 lists of genes (gene symbol) in each pathway.

KEGGgenesetsID:pathway names; list:352 lists of genes (gene id) in each pathway.

KEGGgenesetsSymbol:pathway names; list:352 lists of genes (gene symbol) in each pathway.

The tbl_df of kich gene expression matrix.

A data frame with 3 rows(OS:survival or not, OS.time:survival times , label:the prognosis ) and n columns(samples).

The n x m gene expression matrix, where n is the genes and m the samples.

ReactomegenesetsID:pathway names; list:1635 lists of genes (gene id) in each pathway.

ReactomegenesetsSymbol:pathway names; list:1635 lists of genes (gene symbol) in each pathway.

ReactomegenesetsID:pathway names; list:1635 lists of genes (gene id) in each pathway.

WikipathwaysSymbol:pathway names; list:1635 lists of genes (gene symbol) in each pathway.
}
\usage{
ARIscores

BioCartagenesetsID

BiocartaGenesetsSymbol

GOgenesetsID

GOgenesetsSymbol

GSVAResults

HallmarkgenesetsID

HallmarkgenesetsSymbol

KEGGgenesetsID

KEGGgenesetsSymbol

KICH

KICHcli

KICHfilter

ReactomegenesetsID

ReactomegenesetsSymbol

WikiPathwaysgenesetsID

WikipathwaysSymbol
}
\description{
This dataset includes method names and method names.

This dataset includes pathway names and gene id in each pathway.

This dataset includes pathway names and gene symbol in each pathway.

This dataset includes pathway names and gene id in each pathway.

This dataset includes pathway names and gene symbol in each pathway.

The pathway activity matrix was obtained by GSVA algorithm.

This dataset includes pathway names and gene id in each pathway.

This dataset includes pathway names and gene symbol in each pathway.

This dataset includes pathway names and gene id in each pathway.

This dataset includes pathway names and gene symbol in each pathway.

This dataset is the Kidney Chromophobe(KICH) gene expression matrix.

This dataset includes Survival information.

This dataset includes standardized kich gene expression matrix.

This dataset includes pathway names and gene id in each pathway.

This dataset includes pathway names and gene symbol in each pathway.

This dataset includes pathway names and gene id in each pathway.

This dataset includes pathway names and gene symbol in each pathway.
}
\seealso{
Other related datasets or functions

Other related datasets or functions

Other related datasets or functions

Other related datasets or functions

Other related datasets or functions

Other related datasets or functions

Other related datasets or functions

Other related datasets or functions

Other related datasets or functions

Other related datasets or functions

Other related datasets or functions

Other related datasets or functions

Other related datasets or functions

Other related datasets or functions

Other related datasets or functions

Other related datasets or functions

Other related datasets or functions
}
\keyword{dataset}
