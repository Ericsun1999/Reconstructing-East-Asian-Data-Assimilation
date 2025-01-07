## Code Descriptions

### 1. `Figure1&2.R`
Uses the REACHES data to produce Figures 1 and 2.

### 2. `Figure3.R`
Uses the LME data to produce Figure 3.

### 3. `Figure4.R`
Uses the REACHES data to produce Figure 4.

### 4. `Figure5&6.R`
Uses the REACHES data to produce Figure 5 and 6.

### 5. `Get_tempe_all_data.R`
Uses the REACHES data to produce the kriged REACHES data.
(This step can be skipped when reproducing the figures due to time constraints.)

### 6. `Figure7-9.R`
Uses the kriged REACHES data and the LME data to produce Figures 7, 8, and 9.

### 7. `Choosetuning.R`
Uses the LME data to estimate the tuning parameters described in Section 4.1. (This step can be skipped when reproducing the figures due to time constraints. It should be executed after running `Figure7-9.R`.)

### 8. `Fused.R`
Uses the LME data to estimate the parameters described in Section 4.1. The resulting parameters are saved in `Data\par`. (This step can be skipped when reproducing the figures due to time constraints. It should be executed after running `Figure7-9.R`.)

### 9. `Figure10&11.R`
Uses the estimated parameters for Beijing (from the LME data) to produce Figures 10 and 11.

### 10. `Figure12.R`
Uses the Mean-Square Prediction Error of the kriged REACHES data, the kriged REACHES data itself, and the estimated parameters for Beijing (from the LME data) to produce Figure 12. (It should be executed after running `Figure7-9.R`)

### 11. `Figure13.R`
Uses the kriged REACHES data, the LME data, and the GHCN data to produce Figure 13. (It should be executed after running `Figure12.R`)

### 12. `Figure14.R`
Uses the kriged REACHES data, the LME data, and the GHCN data to produce Figure 14. (It should be executed after running `Figure12.R`)

### 13. `Figure15.R`
Uses the kriged REACHES data, the LME data, and the GHCN data to produce Figure 15. (It should be executed after running `Figure12.R`)

