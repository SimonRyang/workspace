#!/bin/bash

clear

git commit -a -m "edit"
git push origin master

ssh -t -t -X flh71wr@$1 "/home/flh71wr/workspace/entrepreneur_SS/sshrun.sh"
