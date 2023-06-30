<p align="center">
	<img src="logo.jpeg" />
</p>

This repository is built upon the TFHE library (https://github.com/tfhe/tfhe).

#HOW TO BUILD

1. mkdir build
2. cd build
3. cmake ../src -DENABLE_TESTS=on -DENABLE_FFTW=off -DCMAKE_BUILD_TYPE=optim
4. make
