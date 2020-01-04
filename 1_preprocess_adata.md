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
import scanpy as sc
import os
import re
import scipy as scp
import pandas as pd
import scvelo as scv
```

```python
#Read the downloaded output files from the folders
adataNames = os.listdir('./adataObjs/')
adataList = []
for aN in adataNames:
    adataList.append(sc.read_h5ad('./adataObjs/'+aN))
adataList
```

```python

```
