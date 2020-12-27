module load SpectrumMPI
module load OpenMPI

for L in 1.0 3.14
do
    for N in 128 256 512
    do
        for P in 1 2 3 4
        do
        echo "source /polusfs/setenv/setup.SMPI" > "./polus/omp/lsf/"$L"_"$N"_"$P".lsf"
        echo "#BSUB -n "$P >> "./polus/omp/lsf/"$L"_"$N"_"$P".lsf"
        echo "#BSUB -W 00:15" >> "./polus/omp/lsf/"$L"_"$N"_"$P".lsf"
        echo "#BSUB -o polus/omp/out/"$L"_"$N"_"$P".out" >> "./polus/omp/lsf/"$L"_"$N"_"$P".lsf"
        echo "#BSUB -e polus/omp/err/"$L"_"$N"_"$P".err" >> "./polus/omp/lsf/"$L"_"$N"_"$P".lsf"
        echo "#BSUB -R \"span[ptile=1]\"" >> "./polus/omp/lsf/"$L"_"$N"_"$P".lsf"
        echo "OMP_NUM_THREADS=8 mpiexec omp_main -N "$N" -L "$L" -T 0.025" >> "./polus/omp/lsf/"$L"_"$N"_"$P".lsf"
        done
    done
done