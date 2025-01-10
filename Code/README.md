## Code Descriptions

Before running the code, download the required data from the Data folder, and refer to the descriptions in `Data/README.md`.

Follow below code to produce figures for this paper. 

1. `Figure1&2.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) to produce Figures 1 and 2.

2. `Figure3.R` uses the LME data (`Data/LME data`) to produce Figure 3.

3. `Figure4.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) to produce Figure 4.

4. `Figure5&6.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) to produce Figure 5 and 6.

5. `Figure7-9.R` uses the kriged REACHES data (`Data/tempe_all.csv`) and the LME data (`Data/LME data`) to produce Figures 7, 8, and 9.

6. `Figure10&11.R` uses the estimated parameters for Beijing (`Data/par`) to produce Figures 10 and 11.

7. `Figure12.R` uses the Mean-Square Prediction Error of the kriged REACHES data (`Data/std`), the kriged REACHES data itself (`Data/tempe_all.csv`), and the estimated parameters for Beijing (`Data/par`) to produce Figure 12. (It should be executed after running `Figure7-9.R`)
   
8. `Figure13.R` uses the kriged REACHES data (`Data/tempe_all.csv`), the LME data (`Data/LME data`), and the GHCN data (`Data/GHCNv4.xlsx`) to produce Figure 13. (It should be executed after running `Figure12.R`)

9. `Figure14.R` uses the kriged REACHES data (`Data/tempe_all.csv`), the LME data (`Data/LME data`), and the GHCN data (`Data/GHCNv4.xlsx`) to produce Figure 14. (It should be executed after running `Figure12.R`)

10. `Figure15.R` uses the kriged REACHES data (`Data/tempe_all.csv`), the LME data (`Data/LME data`), and the GHCN data (`Data/GHCNv4.xlsx`) to produce Figure 15. (It should be executed after running `Figure12.R`)

The following code are not essential for reproducing the figures presented in the paper. Skipping these will reduce the total runtime to less than one hour.

`Get_tempe_all_data.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) to produce the kriged REACHES data.
(This step can be skipped when reproducing the figures due to time constraints.)

`Choosetuning.R` uses the LME data (`Data/LME data`) to estimate the tuning parameters described in Section 4.1. (This step can be skipped when reproducing the figures due to time constraints. It should be executed after running `Figure7-9.R`.)

`Fused.R` uses the LME data (`Data/LME data`) to estimate the parameters described in Section 4.1. The resulting parameters are saved in `Data/par`. (This step can be skipped when reproducing the figures due to time constraints. It should be executed after running `Figure7-9.R`.)
