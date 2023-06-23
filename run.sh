cmake -B build
cmake --build build
./build/RT 2000 && convert news.ppm news.jpg 
cd ..