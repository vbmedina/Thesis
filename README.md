# Deep-learning for the Discovery of Small-Molecule Antimalarial Drug Leads

This repo contains all the codes used in the mioengineering masters thesis: Deep-learning for the Discovery of Small-Molecule Antimalarial Drug Leads. This repo is meant to help researchers reproduce what has been done in my thesis.

<img width="4200" height="3000" alt="missingno" src="https://github.com/user-attachments/assets/fe43d4ef-16f5-43f7-9872-1aba8f51d12c" />

## ChEMBL Dataset
This study recommends the following pipeline to streamline preprocessing time of the CHMBL364 dataset, especially in scenarios where P. falciparum strains are needed.

<img width="794" height="111" alt="Screenshot 2025-08-31 at 10 24 47 pm" src="https://github.com/user-attachments/assets/6bff7c54-1978-4c48-a033-36d67b006016" />

Out of the 39,624 remaining molecules from the preprocessing pipeline, 18,875 unique molecules and 263 unique strains. We used ChEMBL_35 which was accessed in February of 2025.

## Introduction to the Splits
To simulate real world virtual screening scenarios, this study used the following spliting methods on different D-MPNNs (Chemprop) to evalute cross-stage generalization.

- Random Split: Given the emphasis our study puts on strain preservation, we recommend a random stratified split strategy when keeping molecules with many IC50 values. This is done to benchmark random splitting methods to other studies using unique molecules.
- Scaffold Split: Scaffold Splitting is based on the chemical scaffolds of the molecules (Bemis-Murko).
- Butina Split: The Butina Split uses a distance-based method to groups compounds into clusters based on chemical similarity.
- UMAP-based Clustering Split: UMAP-based Clustering uses Uniform Manifold Approximation and Projection (UMAP) for dimensionality reduction. This study used UMAP twice, once with k-means clustering and another with ward clustering.

<img width="3600" height="1800" alt="max_tanimoto_violins_box_2048" src="https://github.com/user-attachments/assets/68175a68-c66c-4f92-9608-17c7fd3f093f" />

## Models were trained on an 80/20 split, using 5-fold per model.

<img width="4150" height="911" alt="splits_data_vis" src="https://github.com/user-attachments/assets/3bacfb41-4e26-4396-a79a-deb06a108058" />

## Installation
These instructions will guide you through setting up the Conda environment for the project.

### Prerequisites
Make sure you have Conda installed on your system. If not, you can download and install it from [here](https://www.anaconda.com/download).

### Clone the Repository
Clone this repository to your local machine using the following command:

```bash
git clone https://github.com/vbmedina/Thesis.git
cd Thesis
```

### Set Up Conda Environment
Create a Conda environment using the provided `requirements.txt` file. Run the following commands in the project root:

```bash
conda create --name my_environment
conda activate dmpnn
conda install --file requirements.txt
```

### Download and Extract the data
All the data are available on the repo under p0_all_csvs

## Authors

- [@vbmedina](https://www.github.com/vbmedina)

## Contributing

Let me know what you think! 

## Acknowledgments

[1] Guo, Q., Hernandez-Hernandez, S., & Ballester, P. J. “UMAP-based clustering split for rigorous evaluation of AI models for virtual screening on cancer cell lines.” *Journal of Cheminformatics* 17, 94 (2025). https://doi.org/10.1186/s13321-025-01039-8

[2] Stokes, J. M., Yang, K., Swanson, K., Jin, W., *et al.* “A Deep Learning Approach to Antibiotic Discovery.” *Cell* 180(4), 688–702.e13 (2020). https://doi.org/10.1016/j.cell.2020.01.021  (ScienceDirect: https://www.sciencedirect.com/science/article/pii/S0092867420301021)

[3] Carracedo-Reboredo, P., Liñares-Blanco, J., Rodríguez-Fernández, N., Cedrón, F., *et al.* “A review on machine learning approaches and trends in drug discovery.” *Computational and Structural Biotechnology Journal* 19, 4538–4558 (2021). https://doi.org/10.1016/j.csbj.2021.08.011  (ScienceDirect: https://www.sciencedirect.com/science/article/pii/S2001037021003421)
