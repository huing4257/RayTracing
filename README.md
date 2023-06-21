usage:
```shell
#mac
clang++ -std=c++14 -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp src/main.cpp -O3 -o test
./test
```
```shell
#linu
g++ -O2 -fopenmp src/main.cpp src/mesh.cpp  -o test
./test
```