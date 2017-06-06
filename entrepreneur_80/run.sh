#!/bin/bash

clear

git commit -a -m "edit"
git push origin master

ssh -t -t -X flh71wr@$1 << EOF
cd ~/workspace/
git pull origin master
cd ~/workspace/include/
cp * ~/workspace/entrepreneur_80/
cd ~/workspace/entrepreneur_80/
rm *.o
rm *.mod
exit
EOF

ssh -t -t -X flh71wr@$1 "/home/flh71wr/workspace/entrepreneur_80/sshrun.sh $2"
