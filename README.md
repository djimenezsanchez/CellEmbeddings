# Unsupervised Learning of Contextual Information in Multiplex Immunofluorescence Tissue Cytometry

This is the official implementation (Tensorflow and MATLAB) of: 

**Unsupervised Learning of Contextual Information in Multiplex Immunofluorescence Tissue Cytometry**  
Daniel Jiménez-Sánchez, Mikel Ariz, Carlos Ortiz-de-Solórzano.  <a href="https://ieeexplore.ieee.org/document/9098352">Paper</a>

<div align="center">
  <img width="90%" alt="Contextual cell embeddings" src="https://github.com/djimenezsanchez/Contextual_cell_embeddings/blob/main/Method_Description.gif">
</div>

### Abstract:

New machine learning models designed to capture the histopathology of tissues should account not only for the phenotype and morphology of the cells but also learn complex spatial relationships between them. To achieve this, we represent the tissue as an interconnected graph, where previously segmented cells become nodes of the graph. Then the relationships between cells are learned and embedded into a low-dimensional vector, using a Graph Neural Network. This is a fully unsupervised method that learns how to optimally encode cell phenotypes, morphologies, and cell-to-cell interactions from histological tissues labeled using multiplex immunohistochemistry. We provide real multispectral images of human lung adenocarcinoma tissue samples and the code that generates contextual embeddings for each cell tissue. This framework solves limitations of previous image understanding methods whose data representation strategies are either inappropriate for Machine Learning or do not properly capture the interaction between cells. This may render next-generation image understanding deep learning frameworks predicting clinically relevant patient-level labels.
 
 
### Data download

To replicate the paper's experiments on Tissue Cytometry data, first download the images following the <a href="http://doi.org/10.5281/zenodo.4965746">link</a>
We suggest downloading the complete folder and maintaining the original folder structure to prevent image loading errors.

### Installation
This package is compatible with Windows. 
1. Install conda (https://www.anaconda.com/products/individual)
2. Install tensorflow and create a conda environment using the following command: `conda create --name tf_gpu tensorflow-gpu==1.15 networkx==1.11 joblib`

### Usage
Execute main.m 

### Cite
If you make use of this code in your own work, please cite our paper:
```
@inproceedings{jimenez2020unsupervised,
  title={Unsupervised Learning of Contextual Information in Multiplex Immunofluorescence Tissue Cytometry},
  author={Jiménez-Sánchez, Daniel and Ariz, Mikel and Ortiz-de-Solórzano, Carlos},
  booktitle={2020 IEEE 17th International Symposium on Biomedical Imaging (ISBI)},
  pages={1275--1279},
  year={2020},
  organization={IEEE}
}
```
