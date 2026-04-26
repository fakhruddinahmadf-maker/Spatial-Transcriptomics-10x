# Methodology

## Spatial Transcriptomics Analysis Pipeline

**Project:** 10x Genomics Spatial Methods  
**Tools:** Scanpy, Squidpy, AnnData  
**Platforms analyzed:** Visium (spot-based), Xenium (single-cell)

---

## 1. Overview

This project follows a four-stage analysis pipeline that progressively increases in analytical complexity and spatial resolution. Each stage builds on the tools and concepts established in the previous one.

```
Stage 1: Core toolkit (Scanpy)
    ↓
Stage 2: Image feature analysis (Visium Fluorescence)
    ↓
Stage 3: Spatial graph statistics (Visium H&E)
    ↓
Stage 4: Single-cell resolution (Xenium)
```

---

## 2. Data Formats

### AnnData (`.h5ad`)
The central data structure used throughout. An AnnData object contains:

| Slot | Contents |
|---|---|
| `adata.X` | Expression matrix (cells/spots × genes) |
| `adata.obs` | Cell/spot metadata (QC metrics, cluster labels) |
| `adata.var` | Gene metadata (highly variable flags, Moran's I results) |
| `adata.obsm["spatial"]` | 2D spatial coordinates |
| `adata.obsm["X_pca"]` | PCA embedding |
| `adata.obsm["X_umap"]` | UMAP embedding |
| `adata.obsp["spatial_connectivities"]` | Spatial connectivity graph |
| `adata.uns["spatial"]` | Tissue image data |
| `adata.uns["moranI"]` | Moran's I results per gene |

### SpatialData (`.zarr`) — Xenium only
A multi-modal container format for complex spatial datasets:

| Component | Contents |
|---|---|
| Images | Morphology focus image (multiscale) |
| Labels | Cell and nucleus segmentation masks |
| Points | Individual transcript locations (x, y, z) |
| Shapes | Cell and nucleus boundary polygons |
| Tables | AnnData with count matrix and cell metadata |

---

## 3. Quality Control (QC)

### Visium QC (Notebooks 01, 02, 03)

| Metric | Description | Filter threshold |
|---|---|---|
| `total_counts` | Total RNA counts per spot | min=5,000; max=35,000 |
| `n_genes_by_counts` | Unique genes detected per spot | Visualized only |
| `pct_counts_mt` | % mitochondrial gene counts | < 20% |
| Gene minimum cells | Genes appearing in fewer spots | min=10 spots |

**Mitochondrial genes** (prefixed `MT-` in human) are flagged because high MT% indicates cellular stress or apoptosis — the cell membrane has lysed but mitochondria (enclosed in a double membrane) remain intact, artificially inflating their representation.

### Xenium QC (Notebook 04)

| Metric | Description | Filter threshold |
|---|---|---|
| `total_counts` | Total transcripts per cell | min=10 |
| `n_genes_by_counts` | Unique genes per cell | Visualized |
| `cell_area` | Segmented cell area (µm²) | Visualized |
| `nucleus_area / cell_area` | Nucleus-to-cell ratio | Visualized (healthy: 0.2–0.8) |
| `control_probe_counts %` | Negative DNA probe rate | Should be < 0.1% |
| `control_codeword_counts %` | Negative decoding error rate | Should be < 0.1% |

---

## 4. Preprocessing Pipeline

Applied in all four notebooks:

### Step 1: Normalization
```python
sc.pp.normalize_total(adata, inplace=True)
```
Scales each cell/spot so that total counts equal 10,000 (library-size normalization). This removes technical variation in sequencing depth, making cells comparable.

### Step 2: Log-transformation
```python
sc.pp.log1p(adata)
```
Applies log(x+1) transformation. This:
- Reduces skew in count distributions
- Stabilizes variance across the expression range
- Makes the data approximately normally distributed for PCA

### Step 3: Highly Variable Gene Selection (Visium only)
```python
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
```
Selects the 2,000 genes with the most variable expression across spots using the Seurat dispersion method. These genes carry the most biological signal.

### Step 4: PCA
```python
sc.pp.pca(adata)
```
Linear dimensionality reduction. Compresses the gene expression matrix from thousands of dimensions into 50 principal components that capture most of the variance.

### Step 5: KNN Graph
```python
sc.pp.neighbors(adata)
```
Constructs a k-nearest neighbor graph in PCA space. Cells/spots that are transcriptionally similar are connected as neighbors. This graph is the basis for both UMAP and Leiden clustering.

### Step 6: UMAP
```python
sc.tl.umap(adata)
```
Non-linear dimensionality reduction for visualization. Projects the high-dimensional KNN graph into 2D while preserving local neighborhood structure.

### Step 7: Leiden Clustering
```python
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)
```
Graph-based community detection. Identifies groups of cells/spots that are more densely connected to each other than to the rest of the graph. Unlike k-means, does not require pre-specifying the number of clusters.

---

## 5. Image Feature Extraction (Notebook 02)

### Image Segmentation
1. **Smoothing** (`sq.im.process`): Gaussian smoothing to reduce shot noise before segmentation
2. **Watershed segmentation** (`sq.im.segment`): Applied to the DAPI channel. Treats pixel intensity as a topographic surface, "floods" from local maxima (seed points) to identify individual nuclei. Each nucleus is assigned a unique integer label.

### Feature Types

| Feature | Description | Biological use |
|---|---|---|
| **Summary** | Mean, standard deviation, percentiles of pixel intensity per spot | General image brightness and contrast |
| **Histogram** | Binned distribution of pixel intensities | Intensity profile of the spot |
| **Texture** | Gray-Level Co-occurrence Matrix (GLCM) features: contrast, homogeneity, energy, correlation | Tissue microstructure and texture patterns |
| **Segmentation** | Count of segmented objects (cells) per spot; mean intensity per channel per object | Cell density; cell-type proportion proxy |

### Multi-Scale Feature Extraction
Features are extracted at multiple scales (0.25×, 1.0×) and spatial contexts (with/without circle mask) to capture both local and regional tissue properties.

---

## 6. Spatial Graph Analysis (Notebooks 03, 04)

### Spatial Neighbors Graph

**Visium** (`sq.gr.spatial_neighbors`): Builds a connectivity graph based on physical proximity of Visium spots on the array grid.

**Xenium** (`sq.gr.spatial_neighbors(..., coord_type="generic", delaunay=True)`): Uses Delaunay triangulation — an algorithm that connects each point to its geometrically nearest neighbors without pre-specifying a radius, resulting in a triangulation that maximizes the minimum angle of all triangles.

### Neighborhood Enrichment (`sq.gr.nhood_enrichment`)
Tests whether clusters are spatially adjacent more often than expected by chance.

**Method:**
1. For each cluster pair (A, B), count how many spots/cells of cluster B are neighbors of spots/cells of cluster A
2. Permute cluster labels 1,000 times to build a null distribution
3. Compute z-score of observed count vs permuted counts

**Interpretation:** Positive score = enriched spatial proximity. Negative score = spatial avoidance.

### Co-occurrence Analysis (`sq.gr.co_occurrence`)
Measures how the conditional probability of observing cluster B changes with increasing distance from cluster A.

$$\text{Score}(r) = \frac{p(\text{cluster B} \mid \text{within radius } r \text{ of cluster A})}{p(\text{cluster B})}$$

Computed across a range of radii. Score > 1 at radius r means cluster B is enriched within distance r of cluster A.

### Ligand-Receptor Analysis (`sq.gr.ligrec`)
Re-implementation of CellPhoneDB using the OmniPath database of curated ligand-receptor pairs.

**Method:**
1. For each annotated ligand-receptor pair in the database, compute the mean expression of the ligand in the source cluster and the receptor in the target cluster
2. Permute cluster labels 100 times to build a null distribution
3. Compute p-value and adjusted p-value for each pair
4. Report enriched pairs (high mean expression, low adjusted p-value)

### Spatially Variable Genes — Moran's I (`sq.gr.spatial_autocorr`)
Global spatial autocorrelation statistic.

$$I = \frac{N}{\sum_{i}\sum_{j} w_{ij}} \cdot \frac{\sum_{i}\sum_{j} w_{ij}(x_i - \bar{x})(x_j - \bar{x})}{\sum_{i}(x_i - \bar{x})^2}$$

Where:
- N = number of spots/cells
- w_ij = spatial weight between spots i and j (from connectivity graph)
- x_i = expression of gene in spot i

**Interpretation:**
- I ≈ +1: expression is spatially clustered (spatially variable gene)
- I ≈ 0: random spatial pattern
- I ≈ −1: expression is spatially dispersed (checkerboard pattern)

Significance is assessed by comparing the observed I to a permutation-based null distribution (100 permutations).

---

## 7. Centrality Scores (Notebook 04)

Computed for each Leiden cluster using the spatial connectivity graph:

| Score | Definition | Biological meaning |
|---|---|---|
| **Closeness centrality** | Mean shortest path length from a cluster to all other nodes | How spatially accessible/central is this cluster? |
| **Clustering coefficient** | Fraction of a node's neighbors that are also neighbors of each other | How tightly do cells of this type group together? |
| **Degree centrality** | Fraction of all other cells connected to this cluster | What proportion of the tissue does this cluster interface with? |

---

## 8. Software Versions

| Package | Version tested |
|---|---|
| Python | 3.9 |
| scanpy | 1.9.x |
| squidpy | 1.2.x |
| anndata | 0.8.x |
| numpy | 1.23.x |
| pandas | 1.5.x |
| matplotlib | 3.6.x |
| seaborn | 0.12.x |
| igraph | 0.10.x |
| leidenalg | 0.9.x |

---

## 9. References

1. Wolf, F.A., Angerer, P. & Theis, F.J. **SCANPY: large-scale single-cell gene expression data analysis.** *Genome Biology* 19, 15 (2018). https://doi.org/10.1186/s13059-017-1382-0

2. Palla, G., Spitzer, H., Klein, M. et al. **Squidpy: a scalable framework for spatial omics analysis.** *Nature Methods* 19, 171–178 (2022). https://doi.org/10.1038/s41592-021-01358-2

3. Efremova, M., Vento-Tormo, M., Teichmann, S.A. & Vento-Tormo, R. **CellPhoneDB: inferring cell–cell communication from combined expression of multi-subunit ligand–receptor complexes.** *Nature Protocols* 15, 1484–1506 (2020).

4. Türei, D., Korcsmáros, T. & Saez-Rodriguez, J. **OmniPath: guidelines and gateway for literature-curated signaling pathway resources.** *Nature Methods* 13, 966–967 (2016).

5. 10x Genomics. **Visium Spatial Gene Expression.** https://www.10xgenomics.com/products/spatial-gene-expression

6. 10x Genomics. **Xenium In Situ.** https://www.10xgenomics.com/products/xenium-in-situ
