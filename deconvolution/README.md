# Cell type deconvolution

To identify cell-type-specific spatially variable genes (SVGs) while considering the cell type mixtures in the spatial transcriptomic (ST) data, we first need to determine the cell type proportions by using deconvolution methods.
Through our numerical experiments, we found that the estimated cell type proportions play an important role in the cell-type-specific SVG identification, and thus keeping the deconvolution and SVG models consistent in terms of the definition of cell type proportion is essential to controlling the false positive rates and yielding reliable and interpretable results.

Within the **MMM** model, cell type proportion represents the ratio of the transcript count from a certain cell type to all transcripts at the spot.
Here, we highlight its difference with another definition based on the number of cells rather than transcripts.
Therefore, the two most suitable deconvolution methods for **MMM** are [`RCTD`](https://github.com/dmcable/spacexr)
and [`STitch3D`](https://github.com/YangLabHKUST/STitch3D).
Specifically, `STitch3D` is a deep learning-based method that can infer 3D spatial distributions of fine-grained cell types in tissues.
In this repository, we provide a modified version of `STitch3D` with new features for sparsifying the estimated cell type proportions and writing the full deconvolution results to disk.

## References

- Dylan M. Cable, Evan Murray, Luli S. Zou, Aleksandrina Goeva, Evan Z. Macosko, Fei Chen, and Rafael A. Irizarry. Robust decomposition of cell type mixtures in spatial transcriptomics. *Nature Biotechnology* 40, 517–526 (2022). <https://doi.org/10.1038/s41587-021-00830-w>

- Gefei Wang, Jia Zhao, Yan Yan, Yang Wang, Angela Ruohao Wu, and Can Yang. Construction of a 3D whole organism spatial atlas by joint modelling of multiple slices with deep neural networks. *Nature Machine Intelligence* 5, 1200–1213 (2023). <https://doi.org/10.1038/s42256-023-00734-1>
