#!/bin/bash

clear

git commit -a -m "edit"
git push origin master

ssh -t -t flh71wr@$1 << EOF
ulimit -s unlimited
cd workspace/
git pull origin master
cd ~/workspace/include/
cp * ~/workspace/worker/
cd ~/workspace/worker/
rm -f *.mod
rm -f *.o
$2 make
exit
EOF
