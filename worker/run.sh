#!/bin/bash

clear

git commit -a -m "edit"
git push origin master

ssh -t -t flh71wr@$1 << EOF
ulimit -s unlimited
cd workspace/worker/
git pull origin master
rm -r Build
rm *.mod
$2 make
exit
EOF
