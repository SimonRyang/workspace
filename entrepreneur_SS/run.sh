#!/bin/bash

clear

git commit -a -m "edit"
git push origin master

ssh -t -t -X flh71wr@$1 << EOF
ulimit -s unlimited
cd ~/workspace/
git pull origin master
cd ~/workspace/include/
cp * ~/workspace/entrepreneur_SS/
cd ~/workspace/entrepreneur_SS/
$2 make

EOF
