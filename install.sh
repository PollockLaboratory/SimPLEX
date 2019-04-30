#!/bin/bash

# A basic install script for simPLEX.

echo "$(tput setaf 4)Installing simPLEX.$(tput setaf 0)"

if test -f /usr/bin/simPLEX; then
    echo "$(tput setaf 4)Removing old version of simPLEX from /usr/bin$(tput setaf 0)"
    rm -v /usr/bin/simPLEX
fi

echo "$(tput setaf 4)Creating build directory.$(tput setaf 0)"
mkdir -vp build/
cd build/

echo "$(tput setaf 4)Making build scripts.$(tput setaf 0)"
cmake .. -DCMAKE_BUILD_TYPE=RELEASE

echo -e "$(tput setaf 4)Compiling simPLEX.$(tput setaf 0)"
make 

if test -f ./bin/simPLEX; then
    echo -e "$(tput setaf 4)Successful compile, copying binary to /usr/bin.$(tput setaf 0)"
    cp -v ./bin/simPLEX /usr/bin/simPLEX
fi




