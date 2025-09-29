# BRIDGE

**BRIDGE** is a framework for deconvolving bulk expression data into molecular subtype fractions and using these profiles to predict therapy response.

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

```r
# Run BRIDGE predictive model for a given clinical subtype
res_pred <- BRIDGEpredict(expr_matrix = expr,
                          bcor = FALSE,
                          subtype = "ERpos_HER2neg")

# Access results
fractions <- res_pred$fractions
scores    <- res_pred$BRIDGE_SCORE

head(fractions)
head(scores)
```
The function returns a **list** with two elements:

1. **fractions**  
   - Same output as `BRIDGEdeconv()` (samples × subtypes matrix).  
   - Useful for exploring intra-tumoral heterogeneity.  

2. **BRIDGE_SCORE**  
   - A two-column data frame:  
     - **Column 1**: continuous response score (numeric)  
     - **Column 2**: predicted class label (e.g., `"Responder"`, `"Non-responder"`)  

Example:

| Sample | Score  | Prediction   |
|--------|--------|--------------|
| P1     | 0.82   | Responder    |
| P2     | 0.47   | Non-responder|
| P3     | 0.65   | Responder    |

