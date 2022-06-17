# Hybrid Biometric Template Protection: Resolving the Agony of Choice between Bloom Filters and Homomorphic Encryption

## Description

This repository contains the code for experimentally comparing Bloom filter-based against homomorphic encryption-based biometric template protection schemes (BTPs) for the case study of iris recognition (in the `BiometricPerformance` folder) and proposing the proof-of-concept implementation for the hybrid BTP scheme that combines both approaches (in the `HybridBTP` folder) .
There are two categories of BF-based BTPs: the 1st BF-based BTP [RBB13] applies a XOR with a key to the blocks before generating the BF, while the 2nd BF-based BTP [GRG+16] applies a row permutation to the blocks.

Our proposed hybrid BTP scheme combines both approaches: the HE and the 1st BF in the user-specific key setting that requires the system to use a different key for each subject to generate its template.
Hence, the hybrid BTP requires two keys: a key for the subject to transform its iris-code into a BF representation and another key for the system to learn the recognition outcome.
In this hybrid BTP, BF is used as a representation of the iris-code that leads to an accurate biometric comparison while the HE brings confidentiality, unlinkablility and irreversibility to the system.

## Dependencies per folder

### BiometricPerformance folder

This folder contains a Python 3.9 implementation that requires the following packages:

- [`NumPy`](https://numpy.org/)  
- [`SciPy`](https://scipy.org/)

### HybridBTP folder

This folder contains a C++ implementation that requires the following libraries:

- [`PALISADE version v1.11.5`](https://gitlab.com/palisade/palisade-release)
- [`OpenMP`](https://www.openmp.org/)

## Dataset

For the study case of iris recognition, we evaluate the biometric performance of those BTP approaches on the [IITD iris database (Version 1.0)](http://www4.comp.polyu.edu.hk/~csajaykr/IITD/Database_Iris.htm/) that comprises $224$ different subjects and $5$ samples of each eye per subject; we consider the left eye only.
The iris-codes are of dimension $20 \times 512 $ and were extracted using the Log-Gabor (LG) [M03] feature extraction algorithm from the [Iris Toolkit](https://www.wavelab.at/sources/) [RUW+16].
The iris-codes should be placed the `./data` folder.

## Experiments per folder

### BiometricPerformance folder

The following experiments measure the biometric performance of the baseline iris recognition (that corresponds to the HE BTP) and the performance of the 1st BF BTP and the 2nd BF BTP in both the application-specific and user-specific settings.
Note that the biometric performance of the Hybrid BTP corresponds to the performance of the 1st BF BTP in the user-specific setting.
All the performances are evaluated over the [IITD iris database (Version 1.0)](http://www4.comp.polyu.edu.hk/~csajaykr/IITD/Database_Iris.htm/).

- The baseline iris recognition is in file `baselineHD.py`
- The 1st BF in the application-specific setting is in file `appSpec-1stBF.py`
- The 1st BF in the user-specific setting is in file `userSpec-1stBF.py`
- The 2nd BF in the application-specific setting is in file `appSpec-2ndBF.py`
- The 2nd BF in the user-specific setting is in file `userSpec-2ndBF.py`

### HybridBTP folder

The following experiments measure the runtime of the 1st BF BTP, the 2nd BF BTP, the HE BTP, and the Hybrid BTP.

- The 1st BF BTP in the user-specific setting is in file `exp1stBF.cpp`
- The 2nd BF BTP in the user-specific setting is in file `exp2ndBF.cpp`
- The HE BTP for three security levels is in file `expHEPacked.cpp`
- The Hybrid BTP for three security levels is in file `expHybridPacked.cpp`

The following experiments measure the template size of the 1st BF BTP, the 2nd BF BTP, the HE BTP, and the Hybrid BTP, by storing them in binary files.
- The 1st BF template size is in file `serial1stBFTemplate.cpp`
- The 2nd BF template size is in file `serial2ndBFTemplate.cpp`
- The HE template size for three security levels is in file `serialHETemplate.cpp`
- The Hybrid template size for three security levels are in files `serialHybridTemplate.cpp`

Before launching the above experiments, execute the following commands

```
sudo mkdir build
cd build
sudo cmake ..
```

To run any of the above experiments, execute the following commands

```

cd build
sudo make <experiment>
./experiment

```

Commands for HE and Hybrid experiments, add nBits = {128, 192, 256}
```
cd build
sudo make <experiment>
./experiment nBits
```

## References

[[ RBB13 ]](https://ieeexplore.ieee.org/abstract/document/6612976)
[[ GRG+16 ]](https://www.sciencedirect.com/science/article/pii/S0020025516304753)
[[ RUW+16 ]](https://link.springer.com/chapter/10.1007/978-1-4471-6784-6_16)
[[ M03 ]](https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.90.5112&rep=rep1&type=pdf)

## Bibtex Citation

```
Bassit,A., et al.: Hybrid biometric template protection: Resolving the agony of choice between Bloom filters and homomorphic encryption. IETBiome.1–15(2022). https://doi.org/10.1049/bme2.12075 
```
