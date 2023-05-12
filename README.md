# Using Transformer to Extract Structure Information for Multi-Sequences to Boost Function Learning

**Complete Capstone Project by: Aihan Liu, Hsuehyi Lu**  
**Master in Data Science, May 2023**

## Abstract
This repository contains the complete capstone project that focuses on utilizing transformers to extract structure information from multi-sequences in order to enhance function learning. The project explores the application of CoT-Transfer Learning, deepBreaks, and GCN (Graph Convolutional Networks) for RNA structure prediction and function learning. The research investigates the impact of contact matrices, RNA cuts, and feature importance in understanding RNA structure and function.

## Methodology

### Concepts & Background
- **RNA**: RNA (Ribonucleic acid) utilizes four bases (A, U, C, G), and RNA sequences are combinations of these nucleotides.
- **RNA contact**: Tertiary nucleotide-nucleotide interactions play a critical role in determining RNA structure and function. Contacts are defined as physical distances less than 10Å.
- **CoT-Transfer Learning**: Transfer learning from protein-based models can improve RNA contact prediction, addressing the data scarcity issue in RNA structural prediction.
- **deepBreaks**: A machine learning tool that identifies and prioritizes important positions in genotype-phenotype associations. It fits multiple models, selects the best one based on cross-validation score, and predicts the phenotype based on provided sequences.
- **GCN**: Graph Convolutional Networks designed to work with graph-structured data. They leverage convolutional operations to aggregate information from a node's local neighborhood, enabling effective learning and prediction on graphs.

### Future Direction
- **Optimize CoT-Transfer learning**: Improve the mapping method and the transfer learning network to overcome the limitations on sequence length.
- **Further study of cuts**: Research the appropriate definition of cuts in RNA structure and investigate the relationship between the number of cuts and the complexity of RNA structure.
- **GCN structure**: Explore the use of another network to extract a "feature map" for the sequence to enhance prediction results.
- **GCN explainability**: Develop methods, such as gradient-based contrast or class activation mapping, to identify important nodes in GCN for improved interpretability.

## Results
The project achieved high accuracy in contact prediction and provided insights into important positions in the RNA sequence. The GCN model captured complex interactions within RNA structures, enhancing the accuracy of predictions. Results from different datasets and models are shown below:

| Data Description | Model Name                 | Training Size | Accuracy |
|------------------|----------------------------|---------------|----------|
| HIV-1 based on V3 | Logistic Regression        | 35,424, 105  | 0.9925   |
|                  | GCN (real contact)         | 35,424, 105  | 0.9924   |
|                  | GCN (random contact)       | 35,424, 105  | 0.9883   |
|                  | GCN (real contact)         | 100, 105     | 1.0000   |
|                  | GCN (random contact)       | 100, 105     | 0.8000   |
| SARS-CoV-2       | Extra Trees Classifier      | 900, 3822    | 0.9745   |
|                  | GCN (real contact)         | 900, 500     | 0.9495   |
|                  | GCN (random contact)       | 900, 500     | 0.9376   |
|                  | GCN (real contact)        
