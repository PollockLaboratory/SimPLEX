#!/bin/env bash

# A basic script to complile and install SimPLEX.

INSTALL_LOCATION=$HOME/.local/bin
echo "Installing SimPLEX."

if test -f $INSTALL_LOCATION/SimPLEX; then
    echo "Removing old version of SimPLEX from" $INSTALL_LOCATION/SimPLEX
    rm -v $INSTALL_LOCATION/SimPLEX
fi

echo "Creating build directory."
mkdir -vp build/
cd build/

echo "Making build scripts."
cmake .. -DCMAKE_BUILD_TYPE=RELEASE

echo -e "Compiling simPLEX."
make 

if test -f ./bin/SimPLEX; then
    echo -e "Successful compile, copying binary to" $INSTALL_LOCATION/SimPLEX
    install ./bin/SimPLEX $INSTALL_LOCATION/SimPLEX
fi




