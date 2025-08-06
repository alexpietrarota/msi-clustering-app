# MSI Clustering Suite

Questa è un'applicazione **Shiny** per l'analisi e il clustering di dati di **spettrometria di massa da immagini (MSI)**.

L'app permette di:
- Caricare dati in formato `.rds` o `.imzML`
- Scegliere tra diverse metriche di distanza e algoritmi di clustering
- Visualizzare e classificare nuovi tessuti
- Salvare i risultati come report `.csv`

## 📁 Contenuto della cartella

- `lancio_app.R` – Codice principale dell'app Shiny
- `FinalMatrix_Peptides.RDS` – Esempio di matrice intensità
- `FinalMatrix_Coor.RDS` – Esempio di coordinate XY

## ▶️ Esecuzione

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

