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

## Performance (To be done)

