#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scalapack.h"

#define IDX2C(i,j,ld) (((j)*(ld))+(i))

void initialize_matrix(double *A, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            A[IDX2C(i, j, rows)] = drand48();
        }
    }
}

int main(int argc, char **argv) {
    int ictxt, nprow, npcol, myrow, mycol;
    int myrank, nprocs;
    int n, nb;
    int info, descA[9], descB[9], descC[9];
    double alpha = 1.0, beta = 0.0;
    double *A, *B, *C;
    int mpA, nqA, mpB, nqB, mpC, nqC;
    double start, end, local_time, max_time;
    int IZERO = 0, IONE=1, ITWO=2;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc != 3) {
        if (myrank == 0) {
            fprintf(stderr, "Usage: %s <matrix_size> <block_size>\n", argv[0]);
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    n = atoi(argv[1]);
    nb = atoi(argv[2]);

    nprow = (int)sqrt((double)nprocs);
    npcol = nprocs / nprow;

    if (myrank == 0) {
        printf("Running DGEMM benchmark with matrix size %d x %d and block size %d\n", n, n, nb);
        printf("Process grid: %d x %d\n", nprow, npcol);
    }

    Cblacs_pinfo(&myrank, &nprocs);
    Cblacs_get(-1, 0, &ictxt);
    Cblacs_gridinit(&ictxt, "Row", nprow, npcol);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    mpA = numroc_(&n, &nb, &myrow, &IZERO, &nprow);
    nqA = numroc_(&n, &nb, &mycol, &IZERO, &npcol);
    mpB = numroc_(&n, &nb, &myrow, &IZERO, &nprow);
    nqB = numroc_(&n, &nb, &mycol, &IZERO, &npcol);
    mpC = numroc_(&n, &nb, &myrow, &IZERO, &nprow);
    nqC = numroc_(&n, &nb, &mycol, &IZERO, &npcol);

    A = (double *)malloc(mpA * nqA * sizeof(double));
    B = (double *)malloc(mpB * nqB * sizeof(double));
    C = (double *)malloc(mpC * nqC * sizeof(double));

    initialize_matrix(A, mpA, nqA);
    initialize_matrix(B, mpB, nqB);
    initialize_matrix(C, mpC, nqC);

    descinit_(descA, &n, &n, &nb, &nb, &IZERO, &IZERO, &ictxt, &mpA, &info);
    descinit_(descB, &n, &n, &nb, &nb, &IZERO, &IZERO, &ictxt, &mpB, &info);
    descinit_(descC, &n, &n, &nb, &nb, &IZERO, &IZERO, &ictxt, &mpC, &info);

    start = MPI_Wtime();
    pdgemm_("N", "N", &n, &n, &n, &alpha, A, &IONE, &IONE, descA, B, &IONE, &IONE, descB, &beta, C, &IONE, &IONE, descC);
    end = MPI_Wtime();

    local_time = end - start;
    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (myrank == 0) {
        printf("DGEMM completed in %f seconds (max time across all processes)\n", max_time);
    }

    free(A);
    free(B);
    free(C);
    Cblacs_gridexit(0);
    MPI_Finalize();
    return 0;
}
