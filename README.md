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
|                  | GCN (real contact)         | 900, 500     | 0.9376   |
|                  | GCN (random contact)      | 300, 500      | 0.3469   |
|------------------|--------------------------|---------------|----------|

**Note:** The best ML method provided by deepBreaks was used for the logistic regression and extra trees classifier models.

## Structure and Function Relationship
The research indicates that with a large number of samples for training, the contact information has a minimal impact on prediction accuracy. However, the type and order of nucleotides become more critical. On the other hand, with smaller sample sizes, the contact information significantly influences the classification accuracy. The number and distribution of RNA cuts do not exhibit a significant pattern, and further investigation is required.

## Feature Importance
Certain positions in the sequence appear to be more relevant than others. Some of these positions align with the contact matrix, as shown by the red points in the figure. This provides an explainable result for RNA structure prediction. The research suggests that studying the feature importance can offer valuable insights into RNA structure and function.

## Future Work
The project opens up new avenues for further research and exploration. Future directions include:
- Optimizing CoT-Transfer learning by improving the mapping method and transfer learning network to overcome sequence length limitations.
- Conducting additional studies on RNA cuts to define appropriate criteria and examine the relationship between cuts and RNA structure complexity.
- Exploring alternative network architectures to extract a "feature map" for the sequence, potentially improving prediction results.
- Enhancing GCN explainability by employing gradient-based contrast or class activation mapping methods to identify important nodes in the network.

## Acknowledgements
We would like to express our sincere gratitude to our supervisor, Professor Chen Zeng, for his invaluable guidance and support throughout this project. We would also like to thank Professor Edwin Lo for providing valuable insights into our research. Special thanks to our classmates for their contributions and suggestions.

## References
1. Jian et al (Forthcoming), "Knowledge from Large-Scale Protein Contact Prediction Models can be Transferred to the Data-Scarce RNA Contact Prediction Task," Nature Machine Intelligence (submitted).
2. Rahnavard, A., Baghbanzadeh, M., Dawson, T., Sayoldin, B., Oakley, T., & Crandall, K. (2023). "deepBreaks: A Machine Learning Tool for Identifying and Prioritizing Genotype-Phenotype Associations."
3. Kipf, T. N., & Welling, M. (2016). "Semi-Supervised Classification with Graph Convolutional Networks." arXiv preprint arXiv:1609.02907.

---

## Installation for using deepBreaks ##

## Installation ##
* First install *conda*  
Go to the [Anaconda website](https://www.anaconda.com/) and download the latest version for your operating system.  
* For Windows users: do not forget to add `conda` to your system `path`
* Second is to check for conda availability  
open a terminal (or command line for Windows users) and run:
```
conda --version
```
it should out put something like:
```
conda 4.9.2
```
if not, you must make *conda* available to your system for further steps.
if you have problems adding conda to PATH, you can find instructions
[here](https://docs.anaconda.com/anaconda/user-guide/faq/).  

### Windows Linux Mac ###
If you are using an **Apple M1/M2 MAC** please go to the [Apple M1/M2 MAC](#apple-m1m2-mac) for installation
instructions.  
If you have a working conda on your system, you can safely skip to step three.  
If you are using windows, please make sure you have both git and Microsoft Visual C++ 14.0 or greater installed.
install [git](https://gitforwindows.org/)
[Microsoft C++ build tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/)
In case you face issues with this step, [this link](https://github.com/pycaret/pycaret/issues/1254) may help you.
1) Create a new conda environment (let's call it deepBreaks_env) with the following command:
```
conda create --name deepBreaks_env python=3.9
```
2) Activate your conda environment:
```commandline
conda activate deepBreaks_env 
```
3) Install *deepBreaks*:
install with pip:
```commandline
pip install deepBreaks
```
or you can directly install if from GitHub:
```commandline
python -m pip install git+https://github.com/omicsEye/deepbreaks
```
### Apple M1/M2 MAC ###
1) Update/install Xcode Command Line Tools
  ```commandline
  xcode-select --install
  ```
2) Install [Brew](https://brew.sh/index_fr)
  ```commandline
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
  ```
3) Install libraries for brew
  ```commandline
  brew install cmake libomp
  ```
4) Install miniforge
  ```commandline
  brew install miniforge
  ```
5) Close the current terminal and open a new terminal
6) Create a new conda environment (let's call it deepBreaks_env) with the following command:
  ```commandline
  conda create --name deepBreaks_env python=3.9
  ```
7) Activate the conda environment
  ```commandline
  conda activate deepBreaks_env
  ```
8) Install packages from Conda
  ```commandline
  conda install lightgbm
  pip install xgboost
  ```
9) Finally, install *deepBreaks*:
install with pip:
```commandline
pip install deepBreaks
```
or you can directly install if from GitHub:
```commandline
python -m pip install git+https://github.com/omicsEye/deepbreaks
```


Please feel free to customize the template as needed and add any additional sections or information that may be relevant to your specific project.

