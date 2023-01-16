# PCA
The data files and R scripts for the Nature Reviews Methods Primer: Principal Component analysis, by Michael Greenacre, Trevor Hastie, Patrick Groenen, Alfonso Iodice d'Enza, Angelos Markos and Elena Tuzhilina (2022)

### Data files

**world-happiness-report-2021.csv**: CSV file of the World Happiness data

(The Khan cancer data is provided in package ISLR2)

**BarentsFish.csv**: CSV file of the Barents Sea fish data


### R scripts

**Nature_PCA_happy.R**: script for analysing (PCA) the World Happiness data

**Nature_PCA_cancer.R**: script for analysing (PCA and sparse PCA) the Khan cancer data

**Nature_PCA_fish.R**: script for analysing (CA) the Barents Sea fish data

**Nature_happy_missing.R**: script for missing value imputation exercise


### Videos

**Video1_PCA_Centroid_3D.gif**: A three-dimensional animation of the centroid analysis of the four tumour groups.

**Video2_PCA_Regular_To_Centroid.gif**: A dynamic transition from the regular PCA to the PCA of the four tumour group centroids, as weight is transferred from the individual tumours to the tumour group centroids. This shows how the centroid analysis separates the groups better in the two-dimensional PCA solution, as well as how the highly contributing genes change.

**Video3_PCA_Centroid_Regular_To_Sparse.gif**: A dynamic transition from the PCA of the group centroids to the corresponding sparse PCA solution. This shows how most genes are shrunk to the origin, and are thus eliminated, while the others are generally shrunk to the axes, which means they are contributing to only one PC. A few genes still contribute to both PCs.
