---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.1'
      jupytext_version: 1.2.4
  kernelspec:
    display_name: (Fertig) Python 3.7
    language: python
    name: fertig_python_3_7
---

```python
from IPython.core.display import display, HTML
display(HTML("<style>.container { width:90% !important; }</style>"))
%matplotlib inline

import scanpy as sc
import scvelo as scv
scv.logging.print_version()

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.set_figure_params('scvelo')  # for beautified visualization
```

```python
adata = scv.read('processed_adata.h5ad')
scv.utils.show_proportions(adata)
adata
```

https://scvelo-notebooks.readthedocs.io/DentateGyrus.html

```python
scv.pp.filter_and_normalize(adata, min_shared_counts=30, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
```

```python
adata
```

```python
scv.tl.velocity(adata)
```

```python
scv.tl.velocity_graph(adata)
adata
```

```python
sc.tl.umap(adata)
```

```python
adata.write_h5ad('adataWithVelocity.h5ad')
```

```python
scv.pl.velocity_embedding_stream(adata, basis='umap')
```
