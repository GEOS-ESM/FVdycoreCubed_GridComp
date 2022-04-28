import numpy as np

def write_sum_of_vars(comm, data, var_list=None):

    rank = comm.Get_rank()
    nranks = comm.Get_size()

    for i in range(nranks):
        if i == rank:
            print()
            print('P: rank:', rank, flush=True)
            if var_list is None:
                var_list = data.keys()
            for varname in var_list:
                print('P:', varname, ':', data[varname].shape, np.sum(data[varname]), flush=True)
        comm.Barrier()
