mpicc -o dgemm_benchmark dgemm_benchmark.c \
      -I/path/to/scalapack/include \
      -I${MKLROOT}/include \
      -L/path/to/scalapack/lib \
      -L${MKLROOT}/lib/intel64 \
      -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a \
                        ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a \
                        ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
                        ${MKLROOT}/lib/intel64/libmkl_sequential.a \
                        ${MKLROOT}/lib/intel64/libmkl_core.a \
      -Wl,--end-group -L${MKLROOT}/lib/intel64/ -lpthread -lm -ldl 

#mpirun -np <number_of_processes> ./dgemm_benchmark <matrix_size> <block_size>
