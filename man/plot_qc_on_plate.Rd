% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{plot_qc_on_plate}
\alias{plot_qc_on_plate}
\title{Plot QC on Plate}
\usage{
plot_qc_on_plate(
  adata,
  features,
  output_dir = "plots",
  width = 8,
  height = 6,
  dpi = 300,
  facet_var = "Plate"
)
}
\arguments{
\item{adata}{An object containing metadata with `WellID` and `Plate` information.}

\item{features}{A character vector representing the feature names, used for coloring the heatmap and naming the files.}

\item{output_dir}{A character string specifying the directory where the files will be saved (default: "plots").}

\item{width}{Numeric value for the width of the saved plot in inches (default: 8).}

\item{height}{Numeric value for the height of the saved plot in inches (default: 6).}

\item{dpi}{Numeric value specifying the resolution in dots per inch (default: 300).}

\item{facet_var}{A character string specifying the column name to be used for faceting.}
}
\value{
The function does not return a value; it saves the plots as PNG files.
}
\description{
This function generates heatmap plots from metadata for multiple features and saves them as PNG files.
The filenames are determined based on the feature names provided.
}
\examples{
# plot_qc_on_plate(adata, c("percentMT", "nFeature_RNA"), facet_var = "Plate")

}
