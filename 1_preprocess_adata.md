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
samples = []
for an in adataNames:
    samples.append(re.search('S.*[0-9]_',an).group()[:-1])
#for z in zip(adataNames,samples):
    #print(z)
```

```python
for i,adata in enumerate(adataList):
    adata.obs["folder"] = samples[i]
adataList
```

```python
# check if the barcodes are unique
uniqueCheck = []
for adata in adataList:
    uniqueCheck.append(len(set(adata.obs.index)) == len(adata.obs.index))
all(uniqueCheck)
```

```python
#filter cells and genes
for adata in adataList:
    sc.pp.filter_cells(adata,min_genes=200)
    sc.pp.filter_genes(adata,min_cells=3)
adataList
```

```python
adata = adataList[0].concatenate(adataList[1:],join='outer')
adata
```

```python
obs = []
for barcode in adata.obs.index:
    obs.append(barcode)
len(obs) == len(set(obs))
```

```python
#store spliced and unspliced for concatenation separately because it will get concat
s = []
u = []
for ad in adataList:
    s.append(scp.sparse.csr_matrix(ad.layers['spliced']))
    u.append(scp.sparse.csr_matrix(ad.layers['unspliced']))
```

```python
#Convert the spliced and unspliced matrices to pandas dataframe
pdS = []
pdU = []
for i,sm in enumerate(s):
    pdS.append(pd.DataFrame(sm.toarray(),index=adataList[i].obs.index,columns=adataList[i].var.index))
for i,um in enumerate(u): 
    pdU.append(pd.DataFrame(um.toarray(),index=adataList[i].obs.index,columns=adataList[i].var.index))
```

```python
#concatenate the data
pdSCon = pd.concat(pdS,sort=False)
pdUCon = pd.concat(pdU,sort=False)
```

```python
# Rearrange the gene names (cols) to match with the adata concatenation
pdSCon = pdSCon.loc[:,adata.var.index]
pdUCon = pdUCon.loc[:,adata.var.index]
#fill na with 0
pdSCon.fillna(0,inplace=True)
pdUCon.fillna(0,inplace=True)
```

```python
equal = 0
for i,val in enumerate(pdSCon.columns):
    if val == adata.var.index[i]:
        equal +=1
equal == len(pdSCon.columns)
```

```python
scv.utils.show_proportions(adata)
adata
```

```python
adata.obs.head()
```

```python
adata.var.head()
```

```python
adata.obs.drop('batch',axis=1,inplace=True)
adata.obs.head()
```

```python
#convert na to 0
adata.var.fillna(0,inplace=True)
```

```python
df = pd.DataFrame(data=adata.var.apply(sum,1,result_type='reduce'),columns=['n_cells'])
df
```

```python
adata.var = df
adata.var.head()
```

```python
print(adata)
print(adata.obs.head())
print(adata.var.head())
adata.layers['spliced']
```

```python
adata.write('processed_adata.h5ad')
```
