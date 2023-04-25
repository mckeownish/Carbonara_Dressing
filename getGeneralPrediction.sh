#!/bin/bash
# $1 is the moelcule name $2 the output file $3 the restart file if sied -> molXSubstep_Y_Z.dat  X = run (1..30), Y = substep, z = particular mol (this case usually 1)
ScatterFile=newFitData/$1/$1Saxs.dat
fileLocs=newFitData/$1/
initialCoordsFile=newFitData/$1/coordinates.dat
initialCoordsFile=False
noStructures=1
pairedPredictions=none
fixedsections=newFitData/$1/varyingSectionSecondary.dat 
crystalSymmetry=none
withinMonomerHydroCover=none
betweenMonomerHydroCover=none
kmin=0.022;
kmax=0.25;
maxNoFitSteps=20000;

mkdir newFitData/$1/$2


for i in {1..30}
do
    ./predictStructureQvary $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $betweenMonomerHydroCover $kmin $kmax $maxNoFitSteps newFitData/$1/$2/mol$i newFitData/$1/$2/scatter$i.dat newFitData/$1/mixtureFile.dat newFitData/$1/$3 newFitData/$1/$2/fitLog$i.dat
done
 
 
