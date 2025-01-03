# Reconstructing-East-Asian-Data-Assimilation 

Code to reproduce figures from the paper "Reconstructing East Asian Temperatures from 1368 to 1911 Using Historical Documents, Climate Models, and Data Assimilation" by Eric Sun, Kuan-hui Elaine Lin, Wan-Ling Tseng, Pao K. Wang and Hsin-Cheng Huang. 

## Organization

### Data

REACHES data can be downloaded on https://www.ncei.noaa.gov/access/paleo-search/study/37720.

LME data can be download on https://www.cesm.ucar.edu/community-projects/lme, before download the data, please register the account on Climate Data Gateway.

GHCN data can be download on https://www.ncei.noaa.gov/pub/data/ghcn/v4/.

The directory `Getdata` contains instructions and code to download and aggregate the data used for the analysis. See detailed instructions in `Getdata/readme.md`. For convenience, the final aggregated datasets are already provided in the `Data` folder.


### Needed data to plot Figures  

The code `Code/Figure1&2.R` uses the REACHES data to produce Figures 1-2. 

The code `Code/Figure3.R` uses the LME data to produce Figure 3.  (LME)

Figure 4. (REACHES)

Figure 5-6. (REACHES)

Figures 7-9. (REACHES, LME)

Figure 10-12. (REACHES, LME)

Figure 13-15. (REACHES, LME, GHCN)

The R scripts used in main manuscript may be found in the directory `scripts`.

### simulation-results  

Contains the results from running the code in the simulation-code folder as described above. 

### figures  

Produces Figures 1-6 in the paper and Figures S1-S3 in the supplement. All figures except Figures 3 and 6 depend on the simulation-results folder. 

### figures-code  

Contains the results from calling the code in the figures-code folder.

### real-data-code  

The code to run the real data analysis in Section 6.1 is in penguins.R. Instructions on how to download the data can be found [here](https://allisonhorst.github.io/palmerpenguins/articles/download.html). 

The code to run the real data analysis in Section 6.2 is in zheng.R. The T-cell data can be downloaded [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/memory_t), the B-cell data can be downloaded [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/b_cells) and the monocyte data can be downloaded [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/cd14_monocytes).  You'll need to put the three resulting folders in the raw directory inside real-data-code and rename the three folders to 'filtered_matrices_mex_memory', 'filtered_matrices_mex_bcell', and 'filtered_matrices_mex_mono'. 
