module load SpectrumMPI
# mpisubmit.pl --gpu -p 16 -w 00:15 cuda_main
# bsub -oo cuda_profile_2.txt -eo cuda_profile_2.err -n 1 -R "span[ptile=1]" -gpu "num=1" mpiexec nvprof ./cuda_main
bsub -oo polus/result/test.out -eo polus/error/test.err -n 5 -gpu "num=1" mpiexec ./cuda_main
bsub -oo polus/result/test1.out -eo polus/error/test1.err -n 5 mpiexec ./omp_main