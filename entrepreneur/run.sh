#!/bin/bash

clear

git commit -a -m "edit"
git push origin master

ssh -t -t -X flh71wr@$1 << EOF
cd ~/workspace/
git pull origin master
cd ~/workspace/entrepreneur/
#rm globals.o
#rm globals.mod
exit
EOF

ssh -t -t -X flh71wr@$1 "/home/flh71wr/workspace/entrepreneur/sshrun.sh $2"
