# test_gamma_pdf

Create gamma PDF using a numerical method (gamma_func.f90) used in [mizuRoute](https://github.com/nmizukami/mizuRoute "mizuRoute main")
and output in text file.

All the f90 files are from mizuRoute with minor modifications.


To compile,

1. cd ./build

2. cp Makefile Makefile.local (so git does not track your personal makefile)

3. edit Makefile.local for FC, EXE, F_MASTER (see inline comments in Makefile for details)

4. make -f Makefile.local

5. cd ../

6. ./<executable_name>


TODO
 - provide some script(s) to compare with analytical one 
