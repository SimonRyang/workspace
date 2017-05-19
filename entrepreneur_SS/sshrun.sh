#!/bin/bash
ulimit -s unlimited
source ~/intel/compilers_and_libraries_2017/linux/bin/compilervars.sh intel64
cd ~/workspace/entrepreneur_SS/
$1 make
