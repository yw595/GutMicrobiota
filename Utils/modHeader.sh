#!/bin/bash

doTemp=${3:-"false"}

if [ $doTemp == "false" ]; then
    cat $1 | sed -r "s/\s+/|/" | sed -r "s/\s+/_/g" > $2
else
    cat $1 | cut -d"[" -f1 | sed -r "s/\s+/_/g" | grep -P -v "(^>$|^$)" > $2
fi
