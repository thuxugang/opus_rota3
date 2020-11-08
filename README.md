# OPUS-Rota3

Side-chain modeling is critical for protein structure prediction since the uniqueness of the protein structure is largely determined by its side-chain packing conformation. In this paper, differing from most approaches which reply on rotamer library sampling, we first propose a novel side-chain rotamer prediction method based on deep neural networks, named OPUS-RotaNN. Then, on the basis of our previous work OPUS-Rota2, we propose an open-source side-chain modeling framework, OPUS-Rota3, which integrates the results of different methods into its rotamer library as the sampling candidates. By including OPUS-RotaNN into OPUS-Rota3, we conduct our experiments on three native backbone test sets and one non-native backbone test set. On the native backbone test set CAMEO-Hard61 for example, OPUS-Rota3 successfully predicts 51.14% of all side-chain dihedral angles with a tolerance criterion of 20°, outperforms OSCAR-star (50.87%), SCWRL4 (50.40%) and FASPR (49.85%). On the non-native backbone test set DB379-ITASSER, the accuracy of OPUS-Rota3 is 52.49%, better than OSCAR-star (48.95%), FASPR (48.69%) and SCWRL4 (48.29%).

## Datasets

### Train & Val:

The training and validation sets for OPUS-RotaNN are the same as that for OPUS-TASS, which are hosted on [Baidu Drive](https://pan.baidu.com/s/1L6w_qBIKvleO2uFr1Ekevw) with password `zmc1`. Also, they can be downloaded directly from [Here](http://ma-lab.rice.edu/MaLab/dist/opus_tass_datasets.zip). More details can be found in [OPUS-TASS](https://github.com/thuxugang/opus_tass) repo.

The CSF120 features for the training and validation sets are hosted on [Baidu Drive](https://pan.baidu.com/s/15pfFzL4kaYdaV8A37BwqPw) with password `7bmi`. Also, they can be downloaded directly from [Here](http://ma-lab.rice.edu/MaLab/dist/rota3_csf120.zip).

### Test2016 & CASP-FM & CAMEO-Hard61

The pdb files of these three sets are hosted on [Baidu Drive](https://pan.baidu.com/s/1mwsG6OeuOwzmHsWkN1EAJQ) with password `em4s`. Also, they can be downloaded directly from [Here](http://ma-lab.rice.edu/MaLab/dist/rota3_testsets.zip).

### DB379 & DB379-ITASSER:

These two sets can be downloaded from FASPR paper.

## OPUS-RotaNN

### Dependency

```
Python 3.7
TensorFlow v2.0
xssp-3.0.10
```

The training & inference codes of OPUS-RotaNN can be found [Here](https://github.com/thuxugang/opus_rota3/tree/master/OPUS-RotaNN).

The standalone version of OPUS-RotaNN (including pre-trained models, the tmp_files and the results for CASP-FM) is hosted on [Baidu Drive](https://pan.baidu.com/s/11UO508bMR9rOfUYLBA_2fA) with password `3l8k`. Also, it can be downloaded directly from [Here](http://ma-lab.rice.edu/MaLab/dist/OPUS-RotaNN.zip).

Note that for higher speed, we simplify the calculation in xssp. The pipeline for recompiling the xssp is as following:

```
cd mkdssp
tar zxvf xssp-3.0.10.tar.gz
cp dssp.cpp xssp-3.0.10/src/
cd xssp-3.0.10/
./autogen.sh
./configure
make mkdssp
```

## OPUS-Rota3v and OPUS-Rota3

### Dependency

For your convenience, the packages we used can be found in *OPUS-Rota3/dependency* folder.

```
redis
hiredis
```

### Third-Party Programs

#### OSCAR-star

The executable oscar-star file we downloaded from [OSCAR-star website](https://sysimm.ifrec.osaka-u.ac.jp/OSCAR/) can be found in *OPUS-Rota3/oscar* folder.

#### FASPR

The source code we downloaded from [FASPR](https://github.com/tommyhuangthu/FASPR) can be found in *OPUS-Rota3/FASPR* folder. For higher speed, we modified the source code. You can recompile it following their instruction.

#### OPUS-Rota2

The information of OPUS-Rota2 can be found [here](https://github.com/thuxugang/opus_rota2).

### Usage

1. Database preparation.

   1. Install Redis.

      ```
      cd redis-3.2.12
      make
      cd src
      sudo make install
      ```

   2. Download OPUS-DASF lookup tables and rotamer library data, which are hosted on [Baidu Drive](https://pan.baidu.com/s/1XxYz3HpYdpv_DPcZBsMgQA) with password `9dqu`. Also, they can be downloaded directly from [Here](http://ma-lab.rice.edu/dist/rota2_db.zip). Unzip *dump.zip* and put it into *redis-3.2.12*.


2. Install eigen3

   ```
   sudo apt install libeigen3-dev
   cd /usr/include/eigen3
   sudo cp Eigen/ .. -R
   ```
3. Install hiredis

   ```
   cd hiredis
   make
   sudo make install
   sudo ldconfig /usr/local/lib
   ```

4. Build OPUS-Rota3 ***(gcc/g++ 4.8 is necessary to reproduce our results)***

   ```
   cd src
   mkdir build
   cd build
   cmake ..
   make
   cp rota3 ../../
   ```
 
5. Run Redis. 
   ```
   cd /home/redis-3.2.12
   redis-server redis.conf
   ```   

   This will take approximately 1min until *"The server is now ready to accept connections on port 6379"* appeared in the console. 
 
4. Run OPUS-Rota3.

   Set the parameters in *rota3.ini*. The examples for OPUS-Rota3v and OPUS-Rota3 are *rota3v_example.ini* and *rota3_example.ini*, you can replace *rota3.ini* with what you need. Note that ***absolute path*** is required.
   
   1. `optimized`.
   
   `optimized=0` means the sampling process will not be performed. If you want to use OPUS-Rota3v, please set `optimized=0`. For OPUS-Rota3, please set `optimized=1`.
   
   2. `use_faspr`
   
   If you want to use the results from FASPR, please set `use_faspr=1`. For both OPUS-Rota3v and OPUS-Rota3, `use_faspr` is 1.
   
   3. `use_oscar`

   If you want to use the results from OCAR-star, please set `use_oscar=1`. 

   4. `use_others`
   
   Note that OCAR-star will provide stochastic results. Therefore, we provide another option that will not run OCAR-star on-the-fly, but introduce their existing predicted results. In this case, please set `use_oscar=0` and `use_others=1`. The `others_list.txt` in `tmp_dir` should be maintained for establishing the relations between backbone files and OCAR-star predicted results. Instead of introducing the results from OCAR-star, the results from other side-chain modeling programs can also be included in this way. 
  `.dihedrals`: real resid in backbone file, resname, chi1, chi2, chi3, chi4 (181: Missing atoms for calculation 182: chi doesn't exist).
  
   5. `use_rotann`
   
   If you want to use the results from OPUS-RotaNN, please set `use_rotann=1` and maintain the `rotann_list.txt` in `tmp_dir` as `use_others`. Note that the different between `use_others` and `use_rotann` is `.dihedrals` uses real resid in backbone file and `use_rotann` uses the resid starting from 1.
   `.rota`: resid starting from 1, resname, chi1, chi2, chi3, chi4 (181: Missing atoms for calculation 182: chi doesn't exist).

## Performance

### Performance of different predictors on TEST2016

|Predictors |MAE (χ1)	|MAE (χ2)	|MAE (χ3)	|MAE (χ4)	|ACC|
|:----:|:----:|:----:|:----:|:----:|:----:|
|OPUS-RotaNN	|24.75 	|36.64 	|51.62 	|51.65 	|46.26%|
|RotamerLib |34.72 	|46.34 	|53.60 	|54.82 	|46.73%|
|SCWRL4|22.12 |	37.37 |	49.76 |	53.97 	|57.72%|
|FASPR|21.79 |	36.97 |	49.45 |	54.86 |	57.80%|
|OSCAR-star |19.74 	|35.85 |	47.90 	|53.30 	|58.99%|

### Performance of different predictors on DB379

|Predictors |MAE (χ1)	|MAE (χ2)	|MAE (χ3)	|MAE (χ4)	|ACC|
|:----:|:----:|:----:|:----:|:----:|:----:|
|RotamerLib |34.47 |	45.38 |	56.50 |	53.28 	|46.82%|
|FASPR |21.52 |	36.41 |	53.09 	|53.45 |	57.98%|
|SCWRL4 |21.99 |	36.89 	|52.84 |	52.55 |	58.03%|
|OSCAR-star |19.48 |	35.71 |	51.67 	|52.05 |	59.02%|
|OPUS-Rota3v	|18.98 	|33.37 |	49.80 |	52.06 |	60.50%|
|OPUS-Rota3	|17.52 |	32.85 |	48.96 |	50.42 |	61.52%|


### Performance of different predictors on CASP-FM

|Predictors |MAE (χ1)	|MAE (χ2)	|MAE (χ3)	|MAE (χ4)	|ACC|
|:----:|:----:|:----:|:----:|:----:|:----:|
|OPUS-RotaNN|	29.41 |	38.95| 	53.26 	|49.19 |	42.86%|
|RotamerLib |38.07 |	47.53 |	55.65 |	52.86 |	44.95%|
|FASPR |26.63 |	39.75 |	53.40 	|54.81 |	53.11%|
|SCWRL4|27.09 |	40.44 |	52.67 |	54.61 |	53.17%|
|OSCAR-star|24.63 |	37.68 |	50.61 |	53.35 |	54.84%|
|OPUS-Rota3v|	23.05 |	36.05 |	50.82 |	53.63 |	56.65%|
|OPUS-Rota3|	21.38 |	34.50 |	49.07 |	51.51 |	58.05%|
|OPUS-Rota3 (w/ RotaNN)|	21.30 	|33.61 	|49.07 	|51.48 	|58.56%|


### Performance of different predictors on CAMEO-Hard61

|Predictors |MAE (χ1)	|MAE (χ2)	|MAE (χ3)	|MAE (χ4)	|ACC|
|:----:|:----:|:----:|:----:|:----:|:----:|
|OPUS-RotaNN	|32.73 	|42.18| 	56.37 	|51.01 |	38.51%|
|RotamerLib |41.39 |	50.25 |	58.52 	|54.44| 	40.44%|
|FASPR |28.47 	|42.03 	|55.52 |	56.91| 	49.85%|
|OPUS-Rota3v|	28.15 |	41.97 	|55.76 |	57.07 |	50.24%|
|SCWRL4 |28.26 |	42.35| 	55.96 	|55.97 |	50.40%|
|OSCAR-star |26.33 	|41.47| 	55.07 |	56.44 |	50.87%|
|OPUS-Rota3 (w/ RotaNN)|	26.28 |	41.43 |	55.21 |	56.46 |	51.14%|
|OPUS-Rota3	|26.29 |	41.48 |	55.26 |	56.49 |	51.16%|

### Performance of different predictors on DB379-ITASSER (***Non-Native Backbone Test Set***)

|Predictors |MAE (χ1)	|MAE (χ2)	|MAE (χ3)	|MAE (χ4)	|ACC|
|:----:|:----:|:----:|:----:|:----:|:----:|
|RotamerLib |43.10 |	47.64 |	56.94 |	53.27 	|40.61%|
|SCWRL4 |33.29 	|42.74 |	56.17 |	55.52 	|48.29%|
|FASPR |32.59 	|42.05 	|57.05 |	55.69 	|48.69%|
|OSCAR-star |30.78 	|42.45 	|56.07 	|55.74 	|48.95%|
|OPUS-Rota3v	|28.33 	|38.84 |	53.64 	|54.49 	|51.88%|
|OPUS-Rota3	|27.05 	|39.08 |	53.41 	|54.30 |	52.28%|


## Reference 
```bibtex
@article{xu2020opus4,
  title={OPUS-Rota3: Improving Protein Side-Chain Modeling by Deep Neural Networks and Ensemble Methods},
  author={Xu, Gang and Wang, Qinghua and Ma, Jianpeng},
  journal={Journal of Chemical Information and Modeling},
  year={2020},
  publisher={ACS Publications}
}
```
