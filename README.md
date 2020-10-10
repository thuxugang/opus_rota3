# OPUS-Rota3

Side-chain modeling is critical for protein structure prediction since the uniqueness of the protein structure is largely determined by its side-chain packing conformation. In this paper, differing from most approaches which reply on rotamer library sampling, we first propose a novel side-chain rotamer prediction method based on deep neural networks, named OPUS-RotaNN. Then, on the basis of our previous work OPUS-Rota2, we propose an open-source side-chain modeling framework, OPUS-Rota3, which integrates the results of different methods into its rotamer library as the sampling candidates. By including OPUS-RotaNN into OPUS-Rota3, we conduct our experiments on three native backbone test sets and one non-native backbone test set. On the native backbone test set CAMEO-Hard61 for example, OPUS-Rota3 successfully predicts 51.14% of all side-chain dihedral angles with a tolerance criterion of 20°, outperforms OSCAR-star (50.87%), SCWRL4 (50.40%) and FASPR (49.85%). On the non-native backbone test set DB379-ITASSER, the accuracy of OPUS-Rota3 is 52.49%, better than OSCAR-star (48.95%), FASPR (48.69%) and SCWRL4 (48.29%).

## Datasets

### Train & Val:

The training and validation sets for OPUS-RotaNN are the same as that for OPUS-TASS, which are hosted on [Baidu Drive](https://pan.baidu.com/s/1L6w_qBIKvleO2uFr1Ekevw) with password `zmc1`. Also, they can be downloaded directly from [Here](http://ma-lab.rice.edu/MaLab/dist/opus_tass_datasets.zip). More details can be found [here](https://github.com/thuxugang/opus_tass).

### Test2016 & CASP-FM & CAMEO-Hard61

The side-chain dihedrals of these three sets can be found [here](https://github.com/thuxugang/opus_rota3/tree/master/testsets).

 .dihedrals (resid resname x1 x2 x3 x4): 181: Missing atoms for calculation 182: X doesn't exist

### DB379 & DB379-ITASSER:

These two sets can be downloaded from FASPR paper.
