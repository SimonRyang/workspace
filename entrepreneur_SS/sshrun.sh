#!/bin/bash

ulimit -s unlimited
cd ~/workspace/
git pull origin master
cd ~/workspace/include/
cp * ~/workspace/entrepreneur_SS/
cd ~/workspace/entrepreneur_SS/
rm *.o
rm *.mod
make
