all:
	mkdir -p build
	g++ codigo/*.cpp codigo/ppmloader/ppmloader.cpp --std=c++11 -O3 -o build/tp3 -Wall -mavx -DNDEBUG

verificador:
	mkdir -p build
	g++ codigo/verificador/*.cpp codigo/ppmloader/ppmloader.cpp --std=c++11 -O3 -o build/verificar -Wall #-mavx -DNDEBUG

tests:
	mkdir -p build
	g++ tests/*.cpp codigo/*.cpp codigo/ppmloader/ppmloader.cpp --std=c++11 -O3 -o build/tests -Wall -Wno-sign-compare #-mavx -DNDEBUG

clean:
	rm -rf build
