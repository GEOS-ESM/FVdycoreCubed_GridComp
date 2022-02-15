import numpy as np

def write_sum_of_untranslated_arrays(comm, output_data):

    rank = comm.Get_rank()
    nranks = comm.Get_size()

    for i in range(nranks):
        if i == rank:
            print('P: rank:', rank, flush=True)
            print('P: u:',
                  np.sum(output_data['u']),
                  np.sum(output_data['v']),
                  np.sum(output_data['w']),
                  np.sum(output_data['delz']),
                  flush=True)
            print('P: pt:',
                  np.sum(output_data['pt']),
                  np.sum(output_data['delp']),
                  np.sum(output_data['qvapor']),
                  np.sum(output_data['qliquid']),
                  np.sum(output_data['qice']),
                  np.sum(output_data['qrain']),
                  np.sum(output_data['qsnow']),
                  np.sum(output_data['qgraupel']),
                  np.sum(output_data['qcld']),
                  flush=True)
            print('P: ps:',
                  np.sum(output_data['ps']),
                  np.sum(output_data['pe']),
                  np.sum(output_data['pk']),
                  np.sum(output_data['peln']),
                  np.sum(output_data['pkz']),
                  flush=True)
            print('P: phis:',
                  np.sum(output_data['phis']),
                  np.sum(output_data['q_con']),
                  np.sum(output_data['omga']),
                  flush=True)
            print('P: ua:',
                  np.sum(output_data['ua']),
                  np.sum(output_data['va']),
                  np.sum(output_data['uc']),
                  np.sum(output_data['vc']),
                  flush=True)
            print('P: mfx:',
                  np.sum(output_data['mfxd']),
                  np.sum(output_data['mfyd']),
                  np.sum(output_data['cxd']),
                  np.sum(output_data['cyd']),
                  flush=True)
            print('P: diss_est:', np.sum(output_data['diss_estd']), flush=True)
        comm.Barrier()
