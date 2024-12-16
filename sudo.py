#!/usr/bin/env python


# generic code for distributing workload
# does not showcase actual setup of MPI
num_ranks = 0
my_rank = 0

total_alignments = num_seqs - 1
# init MPI
init_mpi(&num_ranks, &my_rank)


if my_rank == 0:
    #num_ranks must be % 2 == 0
        alignments_per_rank = total_alignments / 2

mpi_distribute(sequences, alignments_per_rank)
mpi_barrier()
mpi_gather()
