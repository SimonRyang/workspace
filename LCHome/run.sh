#!/bin/bash

clear

git commit -a -m "edit"
git push origin master

ssh -t -t -X flh71wr@$1 << EOF
cd ~/workspace/
git pull origin master
exit
EOF

ssh -t -t -X flh71wr@$1 "/home/flh71wr/workspace/LCHome/sshrun.sh $2"
