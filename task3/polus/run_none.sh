module load SpectrumMPI
module load OpenMPI

for L in 1.0 3.14
do
    for N in 128 256 512
    do
        for P in 1 2 3 4
        do
        echo $L"_"$N"_"$P
        mpisubmit.pl --stdout polus/none/out/$L"_"$N"_"$P.out --stderr polus/none/err/$L"_"$N"_"$P.err -p $P -w 00:15 main -- -N $N -L $L -T 0.025
        done
    done
done