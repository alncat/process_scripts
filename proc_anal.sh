#!/bin/sh
#dir=/projects/jpma/library/pdb83
dir=/shared.scratch/lzw/tls
dirwork=/Users/alncat/work/tls
listfile=list
for i in `cat $listfile`;do
    echo $i
    cd $i
    $dirwork/self.py -p tlsout.pdb -t phenix -o tlsout > tlslog
    sleep 5
    $dirwork/self.py -p nmout.pdb -t phenix -o nmout > nmlog
    sleep 5
    $dirwork/self.py -p aniso.pdb -t phenix -o origout > origlog
    sleep 5
    cd $dirwork
done

