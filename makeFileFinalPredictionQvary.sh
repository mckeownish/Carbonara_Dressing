#/bin/sh
rm src/*.o

g++ -c -O3 -std=gnu++14 -o src/point.o  src/point.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/polyHelix.o src/polyHelix.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/randomMolGen.o src/randomMolGen.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/ktlMoleculeRandom.o src/ktlMoleculeRandom.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/experimentalData.o src/experimentalData.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/hydrationShellRandom.o src/hydrationShellRandom.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/writheFP.o src/writheFP.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/moleculeFitAndState.o src/moleculeFitAndState.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/mainPredictionFinalQvar.o src/mainPredictionFinalQvar.cpp -fopenmp

g++ -c -O3 -std=gnu++14 -o src/Flexible_generator.o src/Flexible_generator.cpp -fopenmp

g++ -O3 -std=gnu++14 -o predictStructureQvary src/point.o src/polyHelix.o src/randomMolGen.o src/ktlMoleculeRandom.o src/experimentalData.o src/hydrationShellRandom.o src/writheFP.o src/moleculeFitAndState.o src/mainPredictionFinalQvar.o -fopenmp

g++ -O3 -std=gnu++14 -o generate_structure src/point.o src/polyHelix.o src/randomMolGen.o src/ktlMoleculeRandom.o src/experimentalData.o src/hydrationShellRandom.o src/writheFP.o src/moleculeFitAndState.o src/Flexible_generator.o -fopenmp
