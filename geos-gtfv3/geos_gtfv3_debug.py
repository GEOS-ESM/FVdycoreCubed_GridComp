import numpy as np

def write_sum_of_vars(comm, data):

    rank = comm.Get_rank()
    nranks = comm.Get_size()

    for i in range(nranks):
        if i == rank:
            print()
            print('P: rank:', rank, flush=True)
            for varname in data:
                print('P:', varname, ':', np.sum(data[varname]), flush=True)
        comm.Barrier()
