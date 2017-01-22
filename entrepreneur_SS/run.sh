#!/bin/bash

clear

git commit -a -m "edit"
git push origin master

ssh -t -t -X flh71wr@$1 << EOF
ulimit -s unlimited
cd ~/workspace/
git pull origin master
$2 make
exit
EOF
