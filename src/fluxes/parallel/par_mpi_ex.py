#!/usr/bin/env/ python
import numpy as np
import math
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nproc = comm.Get_size()

def main():
    if rank == 0:
        comm.send(rank+1,dest=rank+1,tag=1)
    elif rank == 4:
        comm.send(rank-1,dest=rank-1,tag=2)
    else:
        comm.send(rank+1,dest=rank+1,tag=1)
        comm.send(rank-1,dest=rank-1,tag=2) 
    if rank == 0:
        F = comm.recv(source=rank+1,tag=2)
        B = 0
    elif rank == 4:
        B = comm.recv(source=rank-1,tag=1)
        F = 0
    else:
        F = comm.recv(source=rank+1,tag=2)
        B = comm.recv(source=rank-1,tag=1)
    file = open("Check_%s"%(rank),'a')
    file.write("F:"+str(F)+"B:"+str(B))
    file.close()

main()
