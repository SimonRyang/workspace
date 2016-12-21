#!/bin/bash

clear

git commit -a -m "edit"
git push origin master

ssh -t -t flh71wr@$1 << EOF
ulimit -s unlimited
cd ~/workspace/
git pull origin master
cd ~/workspace/include/
cp * ~/workspace/entrepreneur/
cd ~/workspace/entrepreneur/
rm *.mod
rm *.o
$2 make
exit
EOF
