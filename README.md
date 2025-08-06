# MSI Clustering Suite

This is a **Shiny** application for the analysis and clustering of **Mass Spectrometry Imaging (MSI)** data.

The app allows you to:
- Load data in `.rds`, `.csv`, `.txt`, `.xlsx`, or `.imzML` format  
- Choose between different distance metrics and clustering algorithms  
- Visualize and classify new tissue regions  
- Save results as `.csv` reports

## 📁 Folder contents

- `lancio_app.R` – Main Shiny app script  
- `FinalMatrix_Peptides.RDS` – Example intensity matrix  
- `FinalMatrix_Coor.RDS` – Example XY coordinates

## ▶️ Run the app

```r
shiny::runApp("lancio_app.R")
library(shiny)
library(ggplot2)
library(viridis)
library(dplyr)
library(proxy)
library(readxl)
library(vegan)
library(plotly)
library(data.table)
library(tools)
library(bslib)
library(gridExtra)
library(Cardinal)
library(SummarizedExperiment)
library(MALDIquant)
library(MALDIquantForeign)

