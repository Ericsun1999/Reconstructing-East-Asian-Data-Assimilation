### Data Aggregation

The directory `Data/Get_data` contains instructions and code to download and aggregate the data used for the analysis. For detailed instructions, refer to `Data/Get_data/README.md`. 

For convenience, the following final aggregated datasets datasets are available for analysis:

### 1. REACHES Data:

The temperature index values are provided in `temperature index value.v1.xlsx`.

### 2. LME Data:

Hong Kong: `LME data/d1.csv`.

Shanghai: `LME data/d2.csv`.

Beijing: `LME data/d3.csv`.

Figure7e: All 13 LME time series for produce Figure7(e).

### 3. GHCN Data:

Data is provided in `GHCNv4.xlsx`.

### 4. Kriged REACHES Data (Section 3.1):

Data is provided in `tempe_all_v3.csv` & `tempe_figure6.csv`.

### 5. Mean-Square Prediction Error of Kriged REACHES Data (Section 3.1):

Data is provided in `tempe_all_std.csv`.

### 6. Estimated Parameters for Beijing (Section 4.1):

Parameter mt: `par/mtB.csv`. 

Parameter mu: `par/muB.csv`.

Parameter rt: `par/rtB.csv`.

### 7. Estimated Parameters for Shanghai & Hong Kong (Section 4.1):

These follow the same naming convention as the Beijing parameters, with the last character indicating the location:

For Shanghai, replace the last character with "S" (e.g., `par/mtS.csv`).

For Hong Kong, replace the last character with "H" (e.g., `par/mtH.csv`).

### 8. Summary temperature (Section 4.3):

It includes assimilated, LME, and kriged REACHES temperatures for Beijing, Shanghai, and Hong Kong.

Beijing:   `Valid/tempB5.csv`.

Shanghai:  `Valid/tempS5.csv`.

Hong Kong: `Valid/tempH5.csv`.


