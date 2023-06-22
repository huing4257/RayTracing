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

### 参数曲面求交

$$
x(t)cos\theta = \frac{(y(t)-o_y)}{d_y}d_x+o_x \\
x(t)sin\theta = \frac{(y(t)-o_y)}{d_y}d_z+o_z \\
\Rightarrow x^2(t) = ({ \frac{(y(t)-o_y)}{d_y}})^2(d_x^2+d_z^2)+o_x^2+o_z^2+2(o_x+o_z)*\frac{(y(t)-o_y)}{d_y} \\
\Rightarrow F(t)=({ \frac{(y(t)-o_y)}{d_y}})^2(d_x^2+d_z^2)+o_x^2+o_z^2+2(o_x+o_z)*\frac{(y(t)-o_y)}{d_y}-x^2(t)\\
F'(t)=\frac{2(y(t)-o_y)}{d_y^2}(d_x^2+d_z^2)y'(t) + \frac{2(o_x+o_z)}{d_y}y'(t)-2x(t)x'(t)
$$

