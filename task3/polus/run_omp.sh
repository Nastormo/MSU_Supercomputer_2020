module load SpectrumMPI
module load OpenMPI

for L in 1.0 3.14
do
    for N in 128 256 512
    do
        for P in 1 2 3 4
        do
        echo $L"_"$N"_"$P
        bsub < "./polus/omp/lsf/"$L"_"$N"_"$P".lsf"
        done
    done
done