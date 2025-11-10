# 🧬 TSProm: Tissue-Specific Promoter DNA LLM

![TSProm Pipeline](files/graphical.png)

---

## Overview
**TSProm** is a framework that fine-tunes DNA foundation models (e.g., **DNABERT2**) to decipher the **tissue-specific regulatory grammar** encoded in promoter DNA sequences.  
By comparing models specialized for general vs. tissue-specific promoter logic, TSProm isolates sequence features and motifs that uniquely define **tissue identity**.

---

## Abstract
Characterizing tissue-specific (TSp) gene expression is crucial for understanding development and disease. Traditional expression-based methods, however, often overlook the latent *regulatory grammar* embedded in non-coding DNA, particularly in distal promoter regions.

Here, we introduce **TSProm**, a framework that specializes a DNA foundation model (DNABERT2) to decipher the long-range regulatory logic of TSp promoters at the gene isoform level.

The contributions of this work are twofold:

1. **Comparative Model Design** — Two distinct models are trained:  
   - **Model A:** for general promoter biology  
   - **Model B:** for tissue-specific promoter regulation  
   This comparative strategy isolates motifs that uniquely define tissue identity.

2. **Integrated Explainable AI (xAI)** — Combines attention-based motif discovery with SHAP-based interpretability for robust, cross-validated feature attribution.

Applying TSProm to human **brain**, **liver**, and **testis** promoters, we identify clinically relevant transcription factors (e.g., *SP1, MYC, HES6*) and validate their known roles in diseases such as gliomas and neuroblastomas.  
Our analysis further reveals that **C2H2 Zinc Finger proteins** dominate the global landscape of tissue-specific gene regulation.

**TSProm** provides an interpretable computational framework for identifying regulatory elements that drive tissue-specific gene expression — offering new insights into normal and disease states.

**Keywords:** Explainable AI, Tissue Specificity, DNA Language Models, Interpretable Genomics

---

## Repository Structure
```
TSProm/
├── files/
│   ├── Figure1.pdf
│   ├── DNABERT2_TSP_vs_TransTExRest.csv
│   ├── GENA-LM_TSP_vs_TransTExRest.csv
│   ├── ModelA_best.csv
│   ├── ModelB_best.csv
│   ├── comparison_Models_A_B.svg
│   ├── comparison_methods.svg
│   ├── jobs.json
│   └── jobs1.json
│
├── notebooks/
│   ├── Fig2.ipynb
│   ├── Fig3.ipynb
│   └── Fig4.ipynb
│   ├── runs/
│   │   ├── 1_attention.....ipynb
│   │   └── 2....ipynb
│   │   └── 3_TransSHap.ipynb
│   │   └── 2C_Biclustering.ipynb
│
├── run_scripts/
│   ├── 0_generate_data_run.sh
│   ├── 1_dnabert2_finetune.sh
│   ├── 1_GENA-LM_finetune.sh
│   ├── 2_run_predict.sh
│   └── 3_run_attention.sh
│
├── src/
│   ├── 0_generate_data/
│   │   ├── make_data.R
│   │   ├── make_modelA_negSet.sh
│   │   └── make_modelB_nullSeqs.R
│   │
│   ├── 1_finetune/
│   │   ├── DNABERT2_AttentionExtracted.py
│   │   └── FineTune_GENALM.py
│   │
│   ├── 2_predict/
│   │   ├── 1_predict.py
│   │   └── 1_runAll_predict.py
│   │
│   └── 3_attention/
│       ├── 1_raw_attention_extract-gpu.py
│       ├── 2A_save_meme.py
│       └── 3_SHAP.py
│
└── README.md
```

---

## Quick Setup

Once you clone this repository, run the following commands to install all dependencies required for both **Python** and **R** environments.

```bash
# 1️⃣ Clone the repository
git clone https://github.com/pallavisurana1/TSProm.git
cd TSProm

# 2️⃣ Set up Python environment
python3 -m venv venv
source venv/bin/activate
pip install -r requirements_LLM.txt

# 3️⃣ Set up R environment
# (Installs CRAN + Bioconductor packages automatically)
Rscript install_R_requirements.R

# 4️⃣ Verify installations
python -c "import torch, transformers, pandas; print('✅ Python OK')"
Rscript -e "library(GenomicRanges); library(BSgenome.Hsapiens.UCSC.hg38); cat('✅ R OK\n')"
```

---

## 🚀 Typical Pipeline Flow

| Step  | Script                                             | Description                                                                                                                                                           |
| ----- | -------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **0** | `0_generate_data_run.sh`                           | Generates benchmark datasets and Model A/B inputs. Positive sets can be obtained from [Human TransTEx groupings](https://bmi.cewit.stonybrook.edu/transtexdb/search). |
| **1** | `1_dnabert2_finetune.sh` / `1_GENA-LM_finetune.sh` | Fine-tunes DNABERT2 and GENA-LM models. This step runs Model A and B training for the selected tissue type.                                                           |
| **2** | `2_run_predict.sh`                                 | Runs prediction across train/dev/test splits. Subsets only true predictions for downstream motif analysis.                                                            |
| **3** | `3_run_attention.sh`                               | Executes attention-based motif discovery and SHAP-based interpretability module using enriched JASPAR motifs and TF analyses.                                         |

---

## Example Outputs
- Example CSV outputs and visualizations are provided under the `files/` and `notebooks/`  directories.  `notebooks/runs` directory has information about the attention extraction module runs from the fine-tuned models.
- You can regenerate figures from the paper using `Fig2.ipynb`, `Fig3.ipynb`, and `Fig4.ipynb`.

---

## Model Checkpoints
Fine-tuned **DNABERT2** weights for **Model A** and **Model B** (Brain, Liver, Testis) will be publicly released on **Zenodo** after journal acceptance.  

For pre-release access, please email:  
📧 *pallavi.surana@stonybrook.edu*  
📧 *ramana.davuluri@stonybrookmedicine.edu*  
Include a brief statement of intended use.

---

## Citation
If you use **TSProm** in your research, please cite:


```
Surana, Pallavi, Pratik Dutta, Nimisha Papineni, Rekha Sathian, Zhihan Zhou, Han Liu, and Ramana V. Davuluri. 2025.
“TSProm: Deciphering the Genomic Context of Tissue Specificity.” bioRxiv. https://doi.org/10.1101/2025.10.30.685161
```

---

## License
This project is released under the **MIT License**.  
See [LICENSE](LICENSE) for details.
