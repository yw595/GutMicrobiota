#!/bin/bash

# note had to change fields to 1,1, add grep to remove empty lines
sed -e '/^>/s/$/@/' -e 's/^>/#/' $1 | tr -d '\n' | tr "#" "\n" | tr "@" "\t" | sort -u -t ' ' -f -k 1,1 | sed -e 's/^/>/' -e 's/\t/\n/' | grep -P -v "(>$|^$)" > $2
#rm $1
#mv $2 $1
