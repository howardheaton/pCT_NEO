#!/bin/bash
if [ -e ProtonCTSimulation ]
then
    echo "Removing current 'Main' executable file." 
    rm ProtonCTSimulation
fi

make
#Thsi is ac
g++ Main.cpp -o ProtonCTSimulation   #Compile the main file into an executable
./ProtonCTSimulation
