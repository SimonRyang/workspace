#!/bin/bash

clear

git commit -a -m "edit"
git push origin master

ssh -t -t flh71wr@$1 << EOF
ulimit -s unlimited
cd workspace/include/
git pull origin master
cd workspace/worker/
git pull origin master
cp *.f90 ~/workspace/include/*.f90 ~/workspace/worker
rm -r Build
rm *.mod
$2 make
exit
EOF
