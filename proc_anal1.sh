#!/bin/sh
#dir=/projects/jpma/library/pdb83
dir=/shared.scratch/lzw/tls
dirwork=/Users/alncat/work/tls
listfile=list
for i in `cat $listfile`;do
    echo $i
    cd $i
    $dirwork/twin.py -p tlsout.pdb -t nmout.pdb  > comparelog1
    echo tls
    sleep 2
    $dirwork/twin.py -p aniso.pdb -t tlsout.pdb  > tlsaniso1
    echo nm
    sleep 2
    $dirwork/twin.py -p aniso.pdb -t nmout.pdb  > nmaniso1
    cd $dirwork
done

