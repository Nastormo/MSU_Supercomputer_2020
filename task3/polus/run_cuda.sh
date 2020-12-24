module load SpectrumMPI

for L in 1.0 3.14
do
    for N in 128 256 512
    do
        for P in 01 10 20 40
        do
        echo $L"_"$N"_"$P
        bsub -oo polus/cuda/out/$L"_"$N"_"$P.out -eo polus/cuda/err/$L"_"$N"_"$P.err -n $P -gpu "num=1" mpiexec ./cuda_main -N $N -L $L -T 0.025
        done
    done
done
