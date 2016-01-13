#!/bin/bash
tr "\000" "@" < $1 | grep -n "@"
