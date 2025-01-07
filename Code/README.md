The code `Figure1&2.R` uses the REACHES data to produce Figures 1-2. 

The code `Figure3.R` uses the LME data to produce Figure 3.

The code `Figure4.R` uses the REACHES data to produce Figure 4.

The code `Figure5&6.R` uses the REACHES data to produce Figure 5-6.

The code `Get_tempe_all_data.R` uses the REACHES data to produce the kriged REACHES data. (It can be skip for reproducing the Figures due to time issues)

The code `Figure7-9.R` uses the kriged REACHES data and the LME data to produce Figures 7-9.

The code `Choosetuning.R` uses the LME data to estimate the tuning parameters in section 4.1. (It can be skip for reproducing the Figures due to time issues, this code should be executed after running `Figure7-9.R`)

The code `Fused.R` uses the LME data to estimate the parameters in section 4.1, the result of parameters have been saved in `Data\par`. (It can be skip for reproducing the Figures due to time issues, this code should be executed after running `Figure7-9.R`)

The code `Figure10&11.R` uses the estimated Parameters for Beijing (LME data) to produce Figure 10-11.

The code `Figure12.R` uses the Mean-Square Prediction Error of Kriged REACHES Data, the kriged REACHES data and the estimated Parameters for Beijing (LME data) to produce Figure 12. (This code should be executed after running `Figure7-9.R`)

The code `Figure13.R` uses the kriged REACHES data, the LME data and the GHCN data to produce Figure 13. (This code should be executed after running `Figure12.R`)

The code `Figure14.R` uses the kriged REACHES data, the LME data and the GHCN data to produce Figure 14. (This code should be executed after running `Figure12.R`)

The code `Figure15.R` uses the kriged REACHES data, the LME data and the GHCN data to produce Figure 15. (This code should be executed after running `Figure12.R`)
