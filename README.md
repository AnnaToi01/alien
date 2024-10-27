# Alien Experiments

Group project for Applied Bioinformatics: Microbiome

**Group member** (lexigographical order): Anna Toidze, Marc Kesselring, Minghang Li, Willem Fütterer

## Project Description

In summer 2024, a number of cancer patients mysteriously disappeared from the University Hospital Basel. They inexplicably turned up a single day later, with few memories of what happened. The only common thread in their story: Everyone believes they had been abducted by aliens. These friendly aliens supposedly administered a special, gut microbiota-modifying treatment to these patients in efforts to cure their disease. Lo and behold, some of these patients do appear to miraculously recover.

As a member of the clinical diagnostics team at the University of Zürich Institute of Medical Microbiology, you are tasked with investigating the potential ramifications of these “treatments”. To this end, you are given a dataset containing **fecal microbiome sequences** (“Sequences”) and a metadata table with additional information on the samples (“Sample metadata”) from these patients. To identify the effects of alien abduction, you should examine the relationship between **the fecal microbiome**, **patient covariates**, and **recovery times**.

The patients and their families await your (positive and negative) results with bated breath. The scientific community is already very interested in learning about your approach.

## File Structure

The file structure follows the recommended structure as described in `/ref/project_structure_handout.pdf`.

```text
.
.
├── README.md
├── data
│   ├── processed [29 processed *.qza and *.qzv files]
│   └── raw
│       └── metadata.tsv
├── docs
├── environment.yml
├── ref
├── results [8 visualizations]
├── scripts
│   ├── 00_EDA.ipynb
│   ├── 01_data_download_and_preprocessing.ipynb
│   └── 02_taxonomy.ipynb
└── src
```

- `README.md`: This file.
- `data`: contains all data in the project.
  - `raw`: contains the original data obtained via downloading from polybox. _This folder **should not** be modified._
  - `processed`: contains intermediate processed data, like `.qzv` or `.qza` files. It's not ignored to avoid re-running analysis pipelines. But ideally the synced files should not be larger than `100MB`.
- `docs`: contains all documentation and reports.
- `environment.yml`: a `conda` environment file to recreate the environment used in the project.
- `ref`: contains all reference materials, including the project description, the project structure handout, and any other relevant information related to the group project.
  - Contents in the reference folder might be copied into README.md for better accessibility.
- `results`: contains the final visualizations (`.png`) and tables (`.tsv` or `.csv` files) to put into the report.
- **`scripts`: contains all scripts used in the project.**
  - **`00_EDA.ipynb`**: Jupyter notebook for exploratory data analysis on metadata.
  - **`01_data_download_and_preprocessing.ipynb`**: Jupyter notebook for downloading and preprocessing (denoising, clustering etc) data.
  - **`02_taxonomy.ipynb`**: Jupyter notebook for taxonomy classification via SILVA database.
- `src`: a folder with custom-made modules that will be reused in scripts and notebooks.

## Tasks

### Common Tasks

Each of these tasks should be accompanied by one or more visualizations or tables (as appropriate)

> **Note**: Task 01 should be finished **before midterm**.

- [ ] **01. Sequence quality control.**
  - Appropriate quality control and filtering procedures should be applied.
  - Appropriate denoising and/or clustering techniques should be applied, except if a group is given explicitly pre-processed data.

---

- [ ] **02. Taxonomy classification.**
  - Appropriate taxonomy classification techniques should be applied and explained.
  - Try different methods and databases — you do not need to include this in your report, but you could if the differences are meaningful (e.g., one database clearly gives better results).
- [ ] **03. Alpha diversity analysis.**
  - Should be estimated for all samples, and compared across the primary sample categories(e.g., sample types or other main groups) and/or gradients* (e.g., age, space, time, pH, or other continuous sample metadata).
  - Test out a few different metrics (including a phylogenetic metric) to see what they reveal. Apply appropriate statistical tests.
- [ ] **04. Beta diversity analysis.**
  - Should be measured and compared between the primary sample categories or gradients.
  - Test out a few different metrics (including a phylogenetic metric) to see what they reveal. Apply appropriate statistical tests.
- [ ] **05. Feature differential abundance analysis.**
  - Which features are more/less abundant in different groups?
  - Apply appropriate statistical tests and/or supervised learning methods to answer this question. Consider carefully what types of features you want to compare between groups — maybe use different feature definitions (e.g., collapsed on taxonomy or functional information).
- [ ] 06. [_Optional_] Functional prediction.
  - Use complementary methods to assess microbial functions in a real research/clinical/industrial context (e.g., targeted methods, culture-based methods, metagenomics, etc).
  - Use `q2-picrust2` to predict the metagenome composition of your samples and use this to compare your primary categories or gradients.
  - Repeat beta diversity analysis based on predicted metagenome composition.
- [ ] **07. Do something new.** (_Extra_ work for our group)
  - Pick ~~a~~ **three methods** that we did not learn about in class, and which will add value to your group project (especially if it relates to one of the specific objectives of your group project).
  - Could be a QIIME 2 plugin or action that we did not use, or different statistical or visualization packages for Python or even R.

### Group Project Specific Tasks

- [ ] Do you find any **associations within your metadata** (e.g. between stool consistency and sample day)?
  - Hint: pandas summary statistics on metadata
- [ ]  Explore the **microbial communities of the alleged alien abductees**. What characteristics do you observe? Are you able to identify **predictive biomarkers** for the speed of recovery?
- [ ]  **Bonus**: Some patients gave **more than one stool sample**. How do **microbial communities change** between these timepoints?

## Development Tips

### GitHub in JupyterHub

Please follow the instructions in the [course's GitHub guide](ref/W1_Git_and_GitHub.pdf). Some tips:

- Recommend to generate **classic** OAuth token instead of fine-grained ones, as the latter option requires significant familiarity with GitHub API.
- The OAuth token should at least have `repo`, `workflow` and `read:org` in `admin:org` subsection! -- that's a combination of `gh auth login` requirement and the course's GitHub guide.

### Local environment setup

It is _highly_ recommended to use `libmamba` and the `conda` solver to speed up the environment creation.

#### Install `libmamba`

```bash
# make sure it's at least conda 22.11
conda update -n base conda

# install libmamba
conda install -n base conda-libmamba-solver
conda config --set solver libmamba
```

#### Create `conda` environment

You can create a local `conda` environment using

```bash
conda env create -f environment.yml
```

By default it will create an environment named `qiime2`, and has the exact copy of the packages installed in the course's provided environemnt.

> The configuration file `environment.yml` was almost the same as [QIIME2 Amplicon Distribution](https://docs.qiime2.org/2024.5/install/native/#qiime-2-amplicon-distribution), with `picrust2` installed.
> This is also almost the `conda` environment used in the cluster, except that:
>
> 1. Metagenomics tools (`q2-assembly`, `q2-fondue` etc) are not included
> 2. `mamba` and `jupyterhub` related dependencies are not included
