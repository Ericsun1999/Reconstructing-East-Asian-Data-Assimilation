## Code Descriptions

Before running the code, download the required data from the Data folder, and refer to the descriptions in `Data/README.md`.

Follow the code below to produce figures for this paper. 

1. `Figure2.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) to produce Figures 2(a-c).

2. `Figure3.R` uses the LME data (`Data/LME data/`) to produce Figure 3(a-c).

3. `Figure4.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) to produce Figure 4(a-c).

4. `Figure5.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) to produce Figures 5(a1-b3). (It should be executed after running `Figure4.R`)

5. `Figure6-7.R` uses the Mean-Square Prediction Error of the kriged REACHES data (`Data/tempe_all_std.csv`), the kriged REACHES data (`Data/tempe_all_v3.csv`) and the LME data (`Data/LME data/`) to produce Figures 6(a-d) and 7(a-d).

6. `Figure6e.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) and the LME data (`Data/LME data/Figure6e/b0.csv`) to produce Figure 6(e). (It should be executed after running Figure4.R)

7. `Figure8abc.R` uses the estimated parameters for Beijing (`Data/par/`) to produce Figures 8(a-c).

8. `Figure8d.R` uses the Mean-Square Prediction Error of the kriged REACHES data (`Data/tempe_all_std.csv`), the kriged REACHES data (`Data/tempe_all_v3.csv`) and the estimated parameters for Beijing (`Data/par/`) to produce Figures 8d.

9.  `Figure8d.R` also produce the dataset used to plot Figure9.
   
10.  `Figure9(a1)-(c1).R` uses the Beijing's Valid data (`Data/Valid/tempBv5.csv`) and the GHCN data (`Data/GHCNv4.xlsx`) to produce Figure9(a1)-(c1).

11. `Figure9(a2)-(c2).R` uses the Shanghai's Valid data (`Data/Valid/tempSv5.csv`) and the GHCN data (`Data/GHCNv4.xlsx`) to produce Figure9(a2)-(c2).

12. `Figure9(a3)-(c3).R` uses the Hong Kong's Valid data (`Data/Valid/tempHv5.csv`) and the GHCN data (`Data/GHCNv4.xlsx`) to produce Figure9(a3)-(c3).

The following code is not essential for reproducing the figures presented in the paper. Skipping these will reduce the total runtime to less than one hour.

`Get_tempe_all_data.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) to produce the kriged REACHES data and the Mean-Square Prediction Error of the kriged REACHES data. (It should be executed after running `Figure4.R`)

`Prior.R` uses the LME data (`Data/LME data/`) to estimate the tuning parameters described in Section 4.1.

