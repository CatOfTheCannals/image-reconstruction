g++ tests/*.cpp codigo/*.cpp codigo/ppmloader/ppmloader.cpp --std=c++11 -O3 -o build/tests -Wall -Wno-sign-compare && cd build && ./tests && cd ..
