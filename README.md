# Integrative single-cell and spatial mapping of oxidative stress response uncovers GCLC+ mesenchymal tumor cell state linked with favorable outcomes in triple negative breast cancer

---

## TABLE OF CONTENTS <br>

1. GENERAL INFORMATION
2. DATA PROCESSING
3. FIGURE-TO-CODE MAPPING <br>
---
## :book: General Information

### Integrative single-cell and spatial mapping of oxidative stress response uncovers GCLC+ mesenchymal tumor cell state linked with favorable outcomes in triple negative breast cancer <br><br>
**Authors:**

Nomeda Girnius1,2,3*, Tuulia Vallius2,3,4*, Wenqing Chen5, Inga-Maria Launonen5, Sara Palomino5,6, Jia-Ren Lin2,3, Caitlin E. Mills3, Silja Kauppila5, Pauliina Kronqvist7, Antti Ellonen8, Merja Perala9, Eloise Withnell, Yu-An Chen2, Maria Secrier, Sandro Santagata2,3,4,10, Peter K. Sorger2,3,4, Anniina Farkkila#

**Affiliations**
1.	Department of Cell Biology, Ludwig Center at Harvard, Harvard Medical School, Boston, MA, USA;
2.	Ludwig Center at Harvard;
3.	Laboratory of Systems Pharmacology, Harvard Program in Therapeutic Science, Harvard Medical School, 200 Longwood Avenue, Boston, MA 02115;
4.	Department of Systems Biology, Harvard Medical School, 200 Longwood Avenue, Boston, MA 02115
5.	Research Program in Systems Oncology, University of Helsinki, Helsinki, Finland
6.	Institute for Molecular Medicine Finland, FIMM, University of Helsinki, Helsinki, Finland
7.	Institute of Biomedicine, University of Turku, and Department of Pathology, Turku University Hospital, Turku, Finland
8.	Department of Oncology, Turku University Central Hospital, Turku, Finland
9.	Auria Biobank, University of Turku and Turku University
10.	Department of Pathology, Brigham and Women’s Hospital, Harvard Medical School, Boston, Massachusetts
    
*These authors contributed equally

#Corresponding author
<br><br><br>

## 🧹 Data Processing

### 🔹 t-CycIF
| Script                                       | Description           |
|----------------------------------------------|-----------------------|
| `Codes/t-CycIF/Data_rename.ipynb`            | Re-name the celltypes |
| `Codes/t-CycIF/Scimap_Spatial-Count.ipynb`   | Get RCNs using Scimap |
         
### 🔹 scRNA-seq
| Script                                       | Description           |
|----------------------------------------------|-----------------------|
| `Codes/scRNA-seq/01_scRNAseq_preparation.R`            | scRNAseq data exploration and preparation |
| `Codes/scRNA-seq/02_Antiox_signature.R`   | Compute Antiox signature at single cell level |
| `Codes/scRNA-seq/03_GCLCVIM_signature.R`   | Compute GCLCVIM signature at single cell level |


### 🔹 spatial_transcriptomics
| Script                                       | Description           |
|----------------------------------------------|-----------------------|
| `Codes/spatial_transcriptomics/Spatial_signatures.R` | Explore Visium data and compute signature at spatial level |
| `Codes/spatial_transcriptomics/TNBCs_spatial_preparation.ipynb`   | Deconvolution exploration and preparation for hotpsot analysis |
| `Codes/spatial_transcriptomics/Hotspot_calculation_and_slides_plots.ipynb`   | Hotspot calculation |

<br><br>
## 📊 Figure-to-Code Mapping

### Main Figures

| Figure        | Script                                                |
|---------------|-------------------------------------------------------|
| Figure 1b     | `Codes/scRNA-seq/Oxstress_plots.R` & `Codes/spatial_transcriptomics/Hotspot_calculation_and_slides_plots.ipynb` |
| Figure 1c     | `Codes/scRNA-seq/Oxstress_plots.R`                    |
| Figure 1d     | `Codes/scRNA-seq/Oxstress_plots.R`                    |
| Figure 1e     | `Codes/scRNA-seq/UMAP.R`                              |
| Figure 1f     | `Codes/scRNA-seq/GSEA_Volcano_plots.R`                |
| Figure 1g     | `Codes/spatial_transcriptomics/Hotspot_calculation_and_slides_plots.ipynb`   |
| Figure 1h     | `Codes/spatial_transcriptomics/Spatial_figures.R`  |
| Figure 1i     | `Codes/spatial_transcriptomics/Spatial_figures.R`  |
| Figure 1j     | `Codes/spatial_transcriptomics/Spatial_figures.R`  |
| Figure 2e     | `Codes/t-CycIF/Correlation_plot.R`                    |
| Figure 2h     | `Codes/t-CycIF/GCLC_VIM_plots.ipynb`                  |
| Figure 2i     | `Codes/t-CycIF/GCLC_VIM_plots.ipynb`                  | 
| Figure 2j     | `Codes/scRNA-seq/GCLC_VIM_plots.R`                    | 
| Figure 2k     | `Codes/scRNA-seq/GCLC_VIM_plots.R`                    |
| Figure 3a     | `Codes/t-CycIF/Spatial_plots.ipynb`                   | 
| Figure 3b     | `Codes/t-CycIF/Spatial_plots.ipynb`                   |
| Figure 3d     | `Codes/t-CycIF/Spatial_plots.ipynb`                   |
| Figure 3e     | `Codes/t-CycIF/GCLC_VIM_plots.ipynb`                  | 
| Figure 3g     | `Codes/t-CycIF/Scimap_pScore.ipynb`                   |
| Figure 4a     | `Codes/scRNA-seq/GSEA_Volcano_plots.R`                | 
| Figure 4b     | `Codes/scRNA-seq/GSEA_Volcano_plots.R`                | 
| Figure 4c     | `Codes/scRNA-seq/Proliferation_anastasis_score.ipynb` | 
| Figure 4d     | `Codes/scRNA-seq/Proliferation_anastasis_score.ipynb` | 
| Figure 4e     | `Codes/scRNA-seq/Proliferation_anastasis_score.ipynb` | 
| Figure 4f     | `Codes/scRNA-seq/Proliferation_anastasis_score.ipynb` | 
| Figure 4g     | `Codes/scRNA-seq/UMAP.R`                              | 
| Figure 4h     | `Codes/scRNA-seq/Oxstress_plots.R`                    |
| Figure 4i     | `Codes/scRNA-seq/UMAP.R`                              |
| Figure 4j     | `Codes/spatial_transcriptomics/Spatial_figures.R`     |
| Figure 5d     | `Codes/t-CycIF/Survival.R`                            | 
| Figure 5e     | `Codes/t-CycIF/Survival.R`                            | 



### Supplementary Figures

| Figure        | Script                                |
|---------------|---------------------------------------|
| SFigure 1a    | `Codes/scRNA-seq/Oxstress_plots.R`    |
| SFigure 1g    | `Codes/spatial_transcriptomics/Spatial_figures.R`|
| SFigure 1e    | `Codes/spatial_transcriptomics/Spatial_figures.R`|
| SFigure 2c    | `Codes/scRNA-seq/Antiox_expression.R` |
| SFigure 2d    | `Codes/t-CycIF/Spatial_plots.ipynb`   |
| SFigure 3b    | `Codes/t-CycIF/Moran_violin.ipynb`    | 
| SFigure 3d    | `Codes/t-CycIF/Spatial_plots.ipynb`   | 
| SFigure 4a    | `Codes/scRNA-seq/Oxstress_plots.R`    |
| SFigure 4b    | `Codes/spatial_transcriptomics/Hotspot_calculation_and_slides_plots.ipynb`   |  
| SFigure 4c    | `Codes/spatial_transcriptomics/Spatial_figures.R`  |
| SFigure 4d    | `Codes/spatial_transcriptomics/Spatial_figures.R`  |

---
