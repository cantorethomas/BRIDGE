# BRIDGE

**BRIDGE** is a framework for deconvolving bulk RNA-seq data into molecular subtype fractions and using these profiles to predict therapy response.

---

## Installation

Clone the repository:

```bash
git clone git@github.com:<USERNAME>/BRIDGE.git
cd BRIDGE
```

## Quick Start

Run a minimal example:

```r
# Install required dependencies
install.packages(c("glmnet", "e1071", "Seurat", "GSVA"))   # or see DESCRIPTION

# Load BRIDGE
devtools::load_all(".")

# Example input: expression matrix (genes x samples)
expr <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(expr) <- paste0("Gene", 1:100)
colnames(expr) <- paste0("Sample", 1:10)

# Run deconvolution
bridge_res <- BRIDGE::deconvolve(expr) # Per-sample subtypes abundances 

# Predict therapy response
pred <- BRIDGE::predict_response(bridge_res) # Per-sample prediction of score 

print(pred)
```
