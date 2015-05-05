#!/bin/bash
 
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
pushd . >/dev/null
cd $DIR
echo $DIR
/groups/flyprojects/home/leea30/bin/git-graph.sh | head -n 40
#git rev-list HEAD | head -n 1
git status --porcelain
popd >/dev/null