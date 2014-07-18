#!/bin/bash

if [ -e SimulationData/My2dSim ] # Check to see if the directory exists
then
    echo "The directory /SimulationData/My2dSim/ already exists."
  
else
    echo "Creating directory  /SimulationData/My2dSim/."
    mkdir SimulationData/My2dSim            # Create each directory
    mkdir SimulationData/My2dSim/DetectData
    mkdir SimulationData/My2dSim/MatrixData
fi

echo " " # This space is included to create an extra line gap
         # before outputting that the directory's been created.
echo "Directory creation successful."
