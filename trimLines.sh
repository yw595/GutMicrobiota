#!/bin/bash

awk 'length>255 && $0 ~ />/{len=length;$0=substr($0,0,255)};1' $1 > $2
mv $2 $1
