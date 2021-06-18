# Unsupervised Learning of Contextual Information in Multiplex Immunofluorescence Tissue Cytometry

This is the official implementation (Tensorflow and MATLAB) of CellEmbeddings: 

**Unsupervised Learning of Contextual Information in Multiplex Immunofluorescence Tissue Cytometry**  
Daniel Jiménez-Sánchez, Mikel Ariz, Carlos Ortiz-de-Solórzano.  <a href="https://ieeexplore.ieee.org/document/9098352">Paper</a>.

<div align="center">
  <img width="90%" alt="Contextual cell embeddings" src="https://github.com/djimenezsanchez/Contextual_cell_embeddings/blob/main/Method_Description.gif">
</div>
<div align="center">
  An illustration of CellEmbeddings. 
</div>

### Abstract:

New machine learning models designed to capture the histopathology of tissues should account not only for the phenotype and morphology of the cells but also learn complex spatial relationships between them. To achieve this, we represent the tissue as an interconnected graph, where previously segmented cells become nodes of the graph. Then the relationships between cells are learned and embedded into a low-dimensional vector, using a Graph Neural Network. CellEmbeddings is a fully unsupervised method that learns how to optimally encode cell phenotypes, morphologies, and cell-to-cell interactions from histological tissues labeled using multiplex immunohistochemistry. This framework solves limitations of previous image understanding methods whose data representation strategies are either inappropriate for Machine Learning or do not properly capture the interaction between cells. This may render next-generation image understanding deep learning frameworks predicting clinically relevant patient-level labels.
 
 
### Data download

To replicate the paper's experiments on lung cancer multiplexed data, first download the images following the <a href="http://doi.org/10.5281/zenodo.4965746">link</a>.
Then, add the images to the folder 'Images/Original/'. As a working example, there are two images available in this directory. 

### Installation
This package is compatible with Windows. You'll need a cuda-enabled GPU in your computer system to run CellEmbeddings, visit this <a href="https://medium.com/analytics-vidhya/install-tensorflow-gpu-cuda-in-windows-10-with-easy-to-follow-instructions-614d79782d26">tutorial</a>. You can check that your cuda is installed using the `nvidia-smi` command. Then, follow the instructions:
1. Install conda (https://www.anaconda.com/products/individual)
2. Create a conda environment using the following command: `conda create --name tf_gpu tensorflow-gpu==1.15 networkx==1.11 joblib`

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
