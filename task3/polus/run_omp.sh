module load SpectrumMPI
module load OpenMPI

for L in 1.0 3.14
do
    for N in 128 256 512
    do
        for P in 01 10 20 40
        do
        echo $L"_"$N"_"$P
        mpisubmit.pl --stdout polus/omp/out/$L"_"$N"_"$P.out --stderr polus/omp/err/$L"_"$N"_"$P.err -p $P  -t 4 -w 00:15 omp_main -- -N $N -L $L -T 0.025
        done
    done
done