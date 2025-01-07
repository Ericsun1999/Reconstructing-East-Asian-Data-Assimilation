# Reconstructing-East-Asian-Data-Assimilation 

Code to reproduce figures from the paper "Reconstructing East Asian Temperatures from 1368 to 1911 Using Historical Documents, Climate Models, and Data Assimilation" by Eric Sun, Kuan-hui Elaine Lin, Wan-Ling Tseng, Pao K. Wang and Hsin-Cheng Huang ([arxiv link](http://arxiv.org/abs/2410.21790)). 

## Organization

### Data

We have three different datasets:

REACHES Data: This dataset is reconstructed based on historical document records. It includes information on location (longitude, latitude), levels, and times (years).

LME Data: This is a simulated dataset containing location (longitude, latitude), temperatures, and times (months).

GHCN Data: This dataset consists of the earliest instrumental records, including location (longitude, lati- tude), temperatures, and times (months).

REACHES data can be downloaded on https://www.ncei.noaa.gov/access/paleo-search/study/37720.

LME data can be download on https://www.cesm.ucar.edu/community-projects/lme, before download the data, please register the account on Climate Data Gateway.

GHCN data can be download on https://www.ncei.noaa.gov/pub/data/ghcn/v4/, please click ghcnm.tavg.latest.qcf.tar.gz for downlad.

The directory `Data/Get_data` contains instructions and code to download and aggregate the data used for the analysis. See detailed instructions in `Data/Get_data/README.md`. For convenience, the final aggregated datasets are already provided in the `Data` folder.


### Code 

Detail are described in `Code/README.md`.

The code `Code/Figure1&2.R` uses the REACHES data to produce Figures 1-2. 

The code `Code/Figure3.R` uses the LME data to produce Figure 3.

The code `Code/Figure4.R` uses the REACHES data to produce Figure 4.

The code `Code/Figure5&6.R` uses the REACHES data to produce Figure 5-6.

The code `Code/Figure7-9.R` uses the REACHES data and the LME data to produce Figures 7-9.

The code `Code/Figure10&11.R` uses the LME data to produce Figure 10-11.

The code `Code/Figure12.R` uses the REACHES data and the LME data to produce Figure 12.

The code `Code/Figure13.R` uses the REACHES data, the LME data and the GHCN data to produce Figure 13.

The code `Code/Figure14.R` uses the REACHES data, the LME data and the GHCN data to produce Figure 14.

The code `Code/Figure15.R` uses the REACHES data, the LME data and the GHCN data to produce Figure 15.


### Figure  

Produces Figures 1-15 in the paper.


