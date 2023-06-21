usage:
```shell
#mac
clang++ -std=c++14 -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp src/main.cpp -O3 -o test
./test
```
```shell
#windows
g++ -O3 -fopenmp smallpt.cpp -o smallpt
```