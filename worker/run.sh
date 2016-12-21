#!/bin/bash

clear

git commit -a -m "edit"
git push origin master

ssh -t -t flh71wr@$1 << EOF
ulimit -s unlimited
cd ~/workspace/include/

EOF
