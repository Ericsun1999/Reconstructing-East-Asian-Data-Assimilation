## Code Descriptions

Before running the code, download the required data from the Data folder, and refer to the descriptions in `Data/README.md`.

Follow the code below to produce figures for this paper. 

1. `Figure2.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) to produce Figures 2(a)-(c).

2. `Figure3.R` uses the LME data (`Data/LME data/`) to produce Figures 3(a)-(c).

3. `Figure4.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) to produce Figures 4(a)-(c).

4. `Figure5.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) to produce Figures 5(a1)-(b3). (It should be executed after running `Figure4.R`)

5. `Figure6.R` uses the kriged REACHES data (`Data/tempe_figure6.csv`) to produce Figures 6.

6. `Figure7-8.R` uses the Mean-Square Prediction Error of the kriged REACHES data (`Data/tempe_all_std.csv`), the kriged REACHES data (`Data/tempe_all_v3.csv`) and the LME data (`Data/LME data/`) to produce Figures 7(a)-(d) and 8(a)-(d).

7. `Figure7e.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) and the LME data (`Data/LME data/Figure6e/b0.csv`) to produce Figure 6(e). (It should be executed after running Figure4.R)

8. `Figure9abc.R` uses the estimated parameters for Beijing (`Data/par/`) to produce Figures 8(a)-(c).

9. `Figure9d.R` uses the Mean-Square Prediction Error of the kriged REACHES data (`Data/tempe_all_std.csv`), the kriged REACHES data (`Data/tempe_all_v3.csv`) and the estimated parameters for Beijing (`Data/par/`) to produce Figure 8(d).

10.  `Figure9d.R` also produce the summary temperature dataset used to plot Figures 9.
   
11.  `Figure10a.R` uses the Beijing's summary data (`Data/Valid/tempBv5.csv`) and the GHCN data (`Data/GHCNv4.xlsx`) to produce Figures 10(a).

12. `Figure10b.R` uses the Shanghai's summary data (`Data/Valid/tempSv5.csv`) and the GHCN data (`Data/GHCNv4.xlsx`) to produce Figures 10(b).

13. `Figure10c.R` uses the Hong Kong's summary data (`Data/Valid/tempHv5.csv`) and the GHCN data (`Data/GHCNv4.xlsx`) to produce Figures 10(c).

The following code is not essential for reproducing the figures presented in the paper. Skipping these will reduce the total runtime to less than one hour.

`Get_tempe_all_data.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) to produce the kriged REACHES data and the Mean-Square Prediction Error of the kriged REACHES data. (It should be executed after running `Figure4.R`)

`Prior.R` uses the LME data (`Data/LME data/`) to estimate the tuning parameters described in Section 4.1.

Follow the code below to produce figures for the supplementary.

1. `Supplementary/FigureS1.R` uses the REACHES data (`Data/temperature index value.v1.xlsx`) to produce Figure S1(a)-(e).

2. `Supplementary/FigureS2.R` uses the summary data (`Data/Valid/`) and the GHCN data (`Data/GHCNv4.xlsx`) to produce Figure S2(a)-(c).

3.  `Figure9abc.R` uses the estimated parameters for Shanghai and Hong Kong (`Data/par/`) to produce Figures S3(a)-(c) and S4(a)-(c).

4.  `Figure9d.R` uses the Mean-Square Prediction Error of the kriged REACHES data (`Data/tempe_all_std.csv`), the kriged REACHES data (`Data/tempe_all_v3.csv`) and the estimated parameters for Shanghai and Hong Kong (`Data/par/`) to produce Figure S3(d) and S4(d).
