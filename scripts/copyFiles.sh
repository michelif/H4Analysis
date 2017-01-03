#!/bin/bash

user=micheli
folder=/tmp/micheli/

echo the following items will be copied "in" folder $folder by user $user
for i in $( ls ntuples/); do
    echo  $i
done

read -p "Are you sure? " -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
  echo daje
else
  echo exiting  
fi

scp -r ntuples $user@lxplus046.cern.ch:$folder