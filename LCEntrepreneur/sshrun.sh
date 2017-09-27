#!/bin/bash
ulimit -s unlimited
source ~/intel/compilers_and_libraries_2017/linux/bin/compilervars.sh intel64
cd ~/workspace/LCEntrepreneur/
$1 make -B
