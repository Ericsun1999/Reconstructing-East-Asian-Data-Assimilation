### Data Aggregation

The directory `Data/Get_data` contains instructions and code to download and aggregate the data used for the analysis. For detailed instructions, refer to `Data/Get_data/README.md`. 

For convenience, the following final aggregated datasets datasets are available for analysis:

### 1. REACHES Data:

The temperature index values are provided in `temperature index value.v1.xlsx`.

### 2. LME Data:

Hong Kong: `LME data/a1.csv`.

Shanghai: `LME data/a2.csv`.

Beijing: `LME data/a3.csv`.

Figure9: All 13 LME time series for produce Figure9.

### 3. GHCN Data:

Data is provided in `GHCNv4.xlsx`.

### 4. Kriged REACHES Data (Section 3.1):

Data is provided in `tempe_all.csv`.

### 5. Mean-Square Prediction Error of Kriged REACHES Data (Section 3.1):

Beijing: `std/Bstd.csv`.

Shanghai: `std/BstdS.csv`.

Hong Kong: `std/BstdH.csv`.

### 6. Estimated Parameters for Beijing (Section 4.1):

Parameter mt: `par/mtB.csv`. 

Parameter mu: `par/muB.csv`.

Parameter rt: `par/rtB.csv`.

### 7. Estimated Parameters for Shanghai & Hong Kong (Section 4.1):

These follow the same naming convention as the Beijing parameters, with the last character indicating the location:

For Shanghai, replace the last character with "S" (e.g., `par/mtS.csv`).

For Hong Kong, replace the last character with "H" (e.g., `par/mtH.csv`).


