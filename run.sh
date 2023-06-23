cmake -B build
cmake --build build
./build/RT 10 && convert news.ppm news.jpg 
cd ..