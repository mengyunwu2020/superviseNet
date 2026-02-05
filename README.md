# Supplementary Material for Manuscript

**Title:** Supervised Bayesian Joint Graphical Model for Simultaneous Network Estimation and Subgroup Identification

**Package Name:** superviseNet 

---

## 1. Directory Structure
The unzipped folder is organized as follows:

* **`README.md`**: This instruction file.
* **`superviseNet_0.1.0.tar.gz`**: The source code of the R package.
* **`code_reproduction/`**: Folder containing R scripts to reproduce results.
    * `00_install_package.R`: Script to install dependencies and the package itself.
    * `01_simulation_main.R`: Main script to run simulations.
    * `02_analyze_results.R`: Script to analyze results. 
    
## 2. Instructions for Reproduction

1. Open R or RStudio.
2. Set the working directory to the root of the unzipped directory.
3. Install the package and dependencies:
   ```r
   source("code_reproduction/00_install_package.R")
   ```
4. Run the simulation study:
   ```r
   source("code_reproduction/01_simulation_main.R")
   ```
5. Analyze the results:
   ```r
   source("code_reproduction/02_analyze_results.R")
   ```

