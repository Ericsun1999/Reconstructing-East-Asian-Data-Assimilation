## Code Descriptions

Before running the code, download the required data from the Data folder, and refer to the descriptions in `Data/README.md`.

Follow the code below to produce figures for this paper. 

1. `Figure1&2.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) to produce Figures 1 and 2.

2. `Figure3.R` uses the LME data (`Data/LME data`) to produce Figure 3.

3. `Figure4.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) to produce Figure 4.

4. `Figure5&6.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) to produce Figures 5 and 6.

5. `Figure7-8.R` uses the kriged REACHES data (`Data/tempe_all.csv`) and the LME data (`Data/LME data`) to produce Figures 7 and 8.

6. `Figure9.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) and the LME data (`Data/LME data/Figure9`) to produce Figure 9.

7. `Figure10.R` uses the kriged REACHES data (`Data/tempe_all.csv`) and the LME data (`Data/LME data`) to produce Figure 10. (It should be executed after running `Figure7-8.R`)

8. `Figure11&12.R` uses the estimated parameters for Beijing (`Data/par`) to produce Figures 11 and 12.

9. `Figure13.R` uses the Mean-Square Prediction Error of the kriged REACHES data (`Data/std`), the kriged REACHES data itself (`Data/tempe_all.csv`), and the estimated parameters for Beijing (`Data/par`) to produce Figure 13. (It should be executed after running `Figure7-8.R`)
   
10. `Figure14.R` uses the kriged REACHES data (`Data/tempe_all.csv`), the LME data (`Data/LME data`), and the GHCN data (`Data/GHCNv4.xlsx`) to produce Figure 14. (It should be executed after running `Figure13.R`)

11. `Figure15.R` uses the kriged REACHES data (`Data/tempe_all.csv`), the LME data (`Data/LME data`), and the GHCN data (`Data/GHCNv4.xlsx`) to produce Figure 15. (It should be executed after running `Figure13.R`)

12. `Figure16.R` uses the kriged REACHES data (`Data/tempe_all.csv`), the LME data (`Data/LME data`), and the GHCN data (`Data/GHCNv4.xlsx`) to produce Figure 16. (It should be executed after running `Figure13.R`)

The following code is not essential for reproducing the figures presented in the paper. Skipping these will reduce the total runtime to less than one hour.

`Get_tempe_all_data.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) to produce the kriged REACHES data.

`Choosetuning.R` uses the LME data (`Data/LME data`) to estimate the tuning parameters described in Section 4.1. (It should be executed after running `Figure7-8.R`)

`Fused.R` uses the LME data (`Data/LME data`) to estimate the parameters described in Section 4.1. The resulting parameters are saved in `Data/par`. (It should be executed after running `Figure7-8.R`)
