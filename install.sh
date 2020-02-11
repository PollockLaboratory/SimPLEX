#!/bin/bash

# A basic install script for simPLEX.

echo "Installing simPLEX."

if test -f /usr/bin/simPLEX; then
    echo "Removing old version of simPLEX from /usr/bin"
    rm -v /usr/bin/simPLEX
fi

echo "Creating build directory."
mkdir -vp build/
cd build/

echo "Making build scripts."
cmake .. -DCMAKE_BUILD_TYPE=RELEASE

echo -e "Compiling simPLEX."
make 

if test -f ./bin/simPLEX; then
    echo -e "Successful compile, copying binary to /usr/bin."
    cp -v ./bin/simPLEX /usr/bin/simPLEX
fi




