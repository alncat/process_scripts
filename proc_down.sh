#!/bin/sh
#dir=/projects/jpma/library/pdb83
dir=/shared.scratch/lzw/tls
dirwork=/Users/alncat/work/tls
listfile=list
ssh -fN -L 1234:sugar.rice.edu:22 tianwu@tigger.rice.edu
for i in `cat $listfile`;do
    echo $i
    cd $i
    scp -P 1234 tz4@127.0.0.1:$dir/$i/stat_final .
    awk '{print $NF, $4, $1}' stat_final > stat_final_p
    awk 'BEGIN{count=1} {if(count==1 && $2 < 10) printf("%s 0%s %s\n", $1, $2, $3);else printf("%s %s %s\n", $1, $2, $3);count++} {}' stat_final_p > stat_final_f
    awk '{printf("%s\/%s\/out.pdb%s\n", $1, $2, $3)}' stat_final_f > filetod
    test=0
    for j in `cat filetod`;do
        if [ $test -eq 0 ]; then
            scp -P 1234 tz4@127.0.0.1:$dir/$i/$j tlsout.pdb
            test=1
        else
            scp -P 1234 tz4@127.0.0.1:$dir/$i/$j nmout.pdb
        fi
    done
    test=0
    for j in `cat stat_final_f`;do
        if [ $test -eq 1 ];then
            scp -P 1234 tz4@127.0.0.1:$dir/tlsin/$i"_"$j.phenix phenix
            scp -P 1234 tz4@127.0.0.1:$dir/tlsin/$i"_"tls.pdb pdb
            break
        fi
        ((test=test + 1))
    done
#    $dirwork/self.py -p tlsout.pdb -t phenix -o tlsout > tlslog
#    $dirwork/self.py -p nmout.pdb -t phenix -o nmout > nmlog
#    $dirwork/self.py -p pdb -t phenix -o origout > origlog
    cd $dirwork
done

