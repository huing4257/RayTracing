cmake -B build
cmake --build build
./build/RT 1000 && convert news.ppm news.jpg 
cd ..