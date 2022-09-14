#!/bin/sh
q=67            # number of processes; must be odd
x=1000000000    # compute up to this limit
prec=100        # output precision (in decimal)
for a in `seq 1 $q`; do
	echo "read(\"L.gp\");go($a,$q,$x,\"L.tmp.$a\",$prec);" |gp -q &
done
wait
time sort -n -k1,1 -S 100G -T /home/fo19175/ L.tmp.* > 0_to_1000000000.out
