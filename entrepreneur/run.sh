#!/bin/bash

clear

git commit -a -m "edit"
git push origin master

ssh -t -t -X flh71wr@$1 << EOF
ulimit -s unlimited
cd ~/workspace/
git pull origin master
cd ~/workspace/include/
cp * ~/workspace/entrepreneur/
cd ~/workspace/entrepreneur/
rm *.o
rm *.mod
#make
exit
EOF

ssh -t -t -X flh71wr@$1 "/home/flh71wr/workspace/entrepreneur/sshrun.sh"
