## Code Descriptions

Before running the code, you can download data on `Data`, and the description in `Data/README.md`.

You can follow the order below to run the code for this paper. However, If you want to just reproduce the figure in the paper, you should skip the code `Get_tempe_all_data.R`, `Choosetuning.R`, `Fused.R`, the expected running times would be cut to less than 1 hour.

### 1. `Figure1&2.R`
Uses the REACHES data (`Data/temperature index value.va.xlsx`) to produce Figures 1 and 2.

### 2. `Figure3.R`
Uses the LME data (`Data/LME data`) to produce Figure 3.

### 3. `Figure4.R`
Uses the REACHES data (`Data/temperature index value.va.xlsx`) to produce Figure 4.

### 4. `Figure5&6.R`
Uses the REACHES data (`Data/temperature index value.va.xlsx`) to produce Figure 5 and 6.

### 5. `Get_tempe_all_data.R`
Uses the REACHES data (`Data/temperature index value.va.xlsx`) to produce the kriged REACHES data.
(This step can be skipped when reproducing the figures due to time constraints.)

### 6. `Figure7-9.R`
Uses the kriged REACHES data (`Data/tempe_all.csv`) and the LME data (`Data/LME data`) to produce Figures 7, 8, and 9.

### 7. `Choosetuning.R`
Uses the LME data (`Data/LME data`) to estimate the tuning parameters described in Section 4.1. (This step can be skipped when reproducing the figures due to time constraints. It should be executed after running `Figure7-9.R`.)

### 8. `Fused.R`
Uses the LME data (`Data/LME data`) to estimate the parameters described in Section 4.1. The resulting parameters are saved in `Data\par`. (This step can be skipped when reproducing the figures due to time constraints. It should be executed after running `Figure7-9.R`.)

### 9. `Figure10&11.R`
Uses the estimated parameters for Beijing (`Data\par`) to produce Figures 10 and 11.

### 10. `Figure12.R`
Uses the Mean-Square Prediction Error of the kriged REACHES data (`Data\std`), the kriged REACHES data itself (`Data/tempe_all.csv`), and the estimated parameters for Beijing (`Data\par`) to produce Figure 12. (It should be executed after running `Figure7-9.R`)

### 11. `Figure13.R`
Uses the kriged REACHES data (`Data/tempe_all.csv`), the LME data (`Data/LME data`), and the GHCN data (`Data/GHCNv4.xlsx`) to produce Figure 13. (It should be executed after running `Figure12.R`)

### 12. `Figure14.R`
Uses the kriged REACHES data (`Data/tempe_all.csv`), the LME data (`Data/LME data`), and the GHCN data (`Data/GHCNv4.xlsx`) to produce Figure 14. (It should be executed after running `Figure12.R`)

### 13. `Figure15.R`
Uses the kriged REACHES data (`Data/tempe_all.csv`), the LME data (`Data/LME data`), and the GHCN data (`Data/GHCNv4.xlsx`) to produce Figure 15. (It should be executed after running `Figure12.R`)

