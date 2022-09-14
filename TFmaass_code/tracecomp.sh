#!/bin/sh
time ./trace 14
cat ./tfdata/* > ./tfdata/merged.out
rm -f ./tfdata/*.tmp
awk -F\, '{print>"./tfdata/"$1".txt"}' ./tfdata/merged.out
rm ./tfdata/merged.out
