import numpy as np

def write_sum_of_untranslated_arrays(comm, data):

    rank = comm.Get_rank()
    nranks = comm.Get_size()

    for i in range(nranks):
        if i == rank:
            print('P: rank:', rank, flush=True)
            print('P: u:',
                  np.sum(data['u']),
                  np.sum(data['v']),
                  np.sum(data['w']),
                  np.sum(data['delz']),
                  flush=True)
            print('P: pt:',
                  np.sum(data['pt']),
                  np.sum(data['delp']),
                  np.sum(data['qvapor']),
                  np.sum(data['qliquid']),
                  np.sum(data['qice']),
                  np.sum(data['qrain']),
                  np.sum(data['qsnow']),
                  np.sum(data['qgraupel']),
                  np.sum(data['qcld']),
                  flush=True)
            print('P: ps:',
                  np.sum(data['ps']),
                  np.sum(data['pe']),
                  np.sum(data['pk']),
                  np.sum(data['peln']),
                  np.sum(data['pkz']),
                  flush=True)
            print('P: phis:',
                  np.sum(data['phis']),
                  np.sum(data['q_con']),
                  np.sum(data['omga']),
                  flush=True)
            print('P: ua:',
                  np.sum(data['ua']),
                  np.sum(data['va']),
                  np.sum(data['uc']),
                  np.sum(data['vc']),
                  flush=True)
            print('P: mfx:',
                  np.sum(data['mfxd']),
                  np.sum(data['mfyd']),
                  np.sum(data['cxd']),
                  np.sum(data['cyd']),
                  flush=True)
            print('P: diss_est:', np.sum(data['diss_estd']), flush=True)
        comm.Barrier()
