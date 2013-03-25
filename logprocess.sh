#!/bin/sh
awk 'BEGIN{count=1} {printf("%d %s %s\n", count, $4, $8); count++} END{}' $1 > $1_1
