# 🧬 Spatial Transcriptomics with 10x Genomics

**Submitted by:** Fakhruddin Ahmad
**Course:** Special Topics in Bioinformatics
**Deadline:** 26th April 2026

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue)](https://www.python.org/)
[![Scanpy](https://img.shields.io/badge/Scanpy-1.9%2B-green)](https://scanpy.readthedocs.io/)
[![Squidpy](https://img.shields.io/badge/Squidpy-1.2%2B-orange)](https://squidpy.readthedocs.io/)
[![Platform](https://img.shields.io/badge/Platform-Google%20Colab-yellow)](https://colab.research.google.com/)
[![License](https://img.shields.io/badge/License-MIT-red)](LICENSE)

A complete spatial transcriptomics analysis project following the 10x Genomics pipeline tutorials covering Scanpy basics, Visium fluorescence, Visium H&E, and Xenium single-cell resolution analysis.

---

## 📖 What Is Spatial Transcriptomics?

Spatial transcriptomics answers a fundamental question in biology: **not just which genes are expressed, but *where* they are expressed inside a tissue.** Traditional scRNA-seq destroys spatial information during dissociation. Spatial transcriptomics preserves the physical location of each measurement, allowing us to map gene expression onto tissue architecture.

---

## 🗂️ Project Structure
spatial-transcriptomics-10x/
├── notebooks/
│   ├── 01_scanpy_basics.ipynb
│   ├── 02_visium_fluorescence.ipynb
│   ├── 03_visium_hne.ipynb
│   └── 04_xenium_analysis.ipynb
├── figures/
│   ├── 01_scanpy/
│   ├── 02_visium_fluorescence/
│   ├── 03_visium_hne/
│   └── 04_xenium/
├── data/
│   └── README.md
├── docs/
│   └── methodology.md
├── requirements.txt
└── README.md

---

## 📓 Notebooks Summary

| Notebook | Dataset | Key Analysis |
|----------|---------|-------------|
| 01 — Scanpy Basics | Human Lymph Node (Visium) | QC, clustering, UMAP, spatial plots |
| 02 — Visium Fluorescence | Mouse Brain (fluorescence) | Image segmentation, DAPI/NEUN/GFAP features |
| 03 — Visium H&E | Mouse Brain (H&E) | Spatial graphs, ligand-receptor, Moran's I |
| 04 — Xenium Analysis | Synthetic Lung Cancer | Single-cell resolution, Delaunay graph |

---

## 🔬 Key Concepts Covered

- Visium spot-level analysis
- Image feature extraction
- Spatial graph construction
- Neighborhood enrichment
- Co-occurrence analysis
- Ligand-receptor interactions
- Moran's I spatial autocorrelation
- Delaunay triangulation
- Centrality scores

---

## 📦 Dependencies

```bash
pip install -r requirements.txt
```

| Package | Purpose |
|---------|---------|
| scanpy | Core single-cell analysis |
| squidpy | Spatial transcriptomics |
| anndata | Data format |
| pandas | Data manipulation |
| numpy | Numerical computing |
| matplotlib | Plotting |

---

## 📊 Datasets

| Dataset | Source |
|---------|--------|
| Human Lymph Node | 10x Genomics |
| Mouse Brain Fluorescence | Squidpy built-in |
| Mouse Brain H&E | Squidpy built-in |
| Xenium Lung Cancer | Synthetic |

---

## 🙏 Acknowledgements

- [Scanpy](https://scanpy.readthedocs.io/)
- [Squidpy](https://squidpy.readthedocs.io/)
- [10x Genomics](https://www.10xgenomics.com/)
