# BRIDGE

**BRIDGE** is a framework for deconvolving bulk expression data into molecular subtype fractions and using these profiles to predict therapy response.
All data and code to reproduce the paper can be found at https://hpc.nih.gov/~Lab_ruppin/BRIDGE_main_scripts.zip.

---

## Installation

You can install BRIDGE directly from GitHub:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install BRIDGE
devtools::install_github("cantorethomas/BRIDGE")
```

## How to run 

### 1. Derive subtype abundances `BRIDGEdeconv()`
Please Note: `BRIDGEdeconv()` can be run using 
- `reference='PAM50'`
- `reference='TNBC'`
- `reference='INTCLUST'`

```r
library(BRIDGE)

# Synthetic expression matrix (genes x samples)

set.seed(3)
gsel <- sample(rownames(reference_BRIDGE_PAM50), 600)
expr <- matrix(round(runif(600 * 3), 3) * 1000, nrow = 600,
               dimnames = list(gsel, paste0("Sample", 1:3)))

# Run BRIDGE deconvolution
res <- BRIDGEdeconv(expr_matrix = expr, reference = "PAM50")
```

The function returns a **numeric matrix** with dimensions *(samples × subtypes)*.

- **Rows**: individual samples from the input expression matrix  
- **Columns**: molecular subtypes defined by the chosen reference (e.g., PAM50)  
- **Values**: estimated relative abundance of each subtype within each sample  

Example:

| Sample | BASA | LUMA | LUMB | HER2 | 
|--------|-------|------|------|------|
| P1     | 0.70  | 0.10 | 0.05 | 0.10 |
| P2     | 0.15  | 0.40 | 0.30 | 0.10 |
| P3     | 0.05  | 0.50 | 0.35 | 0.05 |

---

### 2. Deconvolution and response prediction `BRIDGEpredict()`

> ⚠️ **Research use only.** BRIDGE is intended strictly for research purposes and has not been validated for clinical decision-making.

`BRIDGEpredict()` requires two arguments: `subtype` and `therapy`. The following combinations are currently supported:

| `subtype`       | `therapy`    | Model type    |
|-----------------|--------------|---------------|
| `ERpos_HER2neg` | `CHEMO`      | Main          |
| `ERpos_HER2neg` | `IMMUNO`     | Exploratory   |
| `ERpos`         | `ENDO`       | Exploratory   |
| `HER2pos`       | `ANTI_HER2`  | Main          |
| `TNBC`          | `CHEMO`      | Main          |
| `TNBC`          | `IMMUNO`     | Exploratory   |

> **Main** models have been trained and validated on multiple cohorts. **Exploratory** models are based on limited data and should be interpreted with caution.

```r
# Run BRIDGE predictive model for a given subtype and therapy
res_pred <- BRIDGEpredict(expr_matrix = expr,
                          subtype     = "ERpos_HER2neg",
                          therapy     = "CHEMO")

# Access results
fractions <- res_pred$fractions
scores    <- res_pred$BRIDGE_SCORE

head(fractions)
head(scores)
```

The function returns a **list** with two elements:

1. **`fractions`**  
   Same output as `BRIDGEdeconv()` (samples × subtypes matrix). Useful for exploring intra-tumoral heterogeneity independently of the predictive model.

2. **`BRIDGE_SCORE`**  
   A data frame with three columns:

   | Column      | Description                                              |
   |-------------|----------------------------------------------------------|
   | `SCORE`     | Continuous response score (0–1)                         |
   | `CLASS`     | Predicted class: `"high"` or `"low"` response           |
   | `SUBCLASS`  | Quantile-based sub-class (`q1`–`q5`); `NA` if unavailable |

   Example output:

   |        | SCORE | CLASS | SUBCLASS |
   |--------|-------|-------|----------|
   | P1     | 0.82  | high  | q5       |
   | P2     | 0.47  | low   | q3       |
   | P3     | 0.65  | high  | q4       |

