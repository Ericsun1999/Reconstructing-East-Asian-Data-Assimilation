# Reconstructing-East-Asian-Data-Assimilation 

Code to reproduce figures from the paper "Reconstructing East Asian Temperatures from 1368 to 1911 Using Historical Documents, Climate Models, and Data Assimilation" by Eric Sun, Kuan-hui Elaine Lin, Wan-Ling Tseng, Pao K. Wang and Hsin-Cheng Huang ([arxiv link](http://arxiv.org/abs/2410.21790)). 

## Organization

### Code 

Detailed descriptions of the code and usage instructions are available in `Code/README.md`.

### Data

Detailed descriptions of the data and usage instructions are available in `Data/README.md`.

#### REACHES Data:

Reconstructed from historical document records, this dataset includes information on location (longitude, latitude), levels, and time (years).

Download Link: [NOAA Paleo Search - Study 37720](https://www.ncei.noaa.gov/access/paleo-search/study/37720.)

#### LME Data:

A simulated dataset containing information on location (longitude, latitude), temperatures, and time (months).

Download Link: [CESM Community Projects - LME](https://www.cesm.ucar.edu/community-projects/lme)
(Note: You need to register for an account on the Climate Data Gateway before downloading.)

#### GHCN Data: 

Composed of the earliest instrumental records, this dataset includes location (longitude, latitude), temperatures, and time (months).

Download Link:[ NOAA GHCN v4](https://www.ncei.noaa.gov/pub/data/ghcn/v4/)
(Click on ghcnm.tavg.latest.qcf.tar.gz to download.)

#### Data Aggregation

The directory `Data/Get_data` contains instructions and code to download and aggregate the data used for the analysis. See detailed instructions in `Data/Get_data/README.md`. 

For convenience, the final aggregated datasets are already provided in the `Data` folder.


### Figure  

The provided code generates Figures 1–15 as presented in the paper.


