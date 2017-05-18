#!/bin/bash

clear

git commit -a -m "edit"
git push origin master

ssh -t -t -X flh71wr@$1 << EOF
cd ~/workspace/
git pull origin master
cd ~/workspace/include/
cp * ~/workspace/test/
cd ~/workspace/test/
rm *.o
rm *.mod
exit
EOF

ssh -t -t -X flh71wr@$1 "/home/flh71wr/workspace/test/sshrun.sh"
