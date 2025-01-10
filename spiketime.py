import FSIN as sim
import numpy as np
import matplotlib.pyplot as plt
import sys

method = "euler";
dt = 5.e-4;
sim_time=200.;

gms=1.65

kneu = 43

with open('observation.txt','w') as file:
    for x in range(60):
        dv=0.1*x
        rates, spMon = sim.gewnet(kneu=kneu,sim_time=sim_time,return_state=False,Nrec=100,dt=dt,method=method,act_syns=True,mod_gL=False, gjs=False,gext=7.,sin_cond=False,const_cond=True,tau_rise=.3,homo_neurons=True,gms=gms,dist_syns=False,read_delays=False, read_gms=False,dt_rec=.1,dmin=dv,dmax=dv,connectivity="FID",return_syns=False,read_syns=False,USEm=1.,std=False, Es=-55.,randics=False,one_perturbed=True); 
        spike_trains = spMon.spike_trains();
#with open('output.txt', 'w') as file:
#   for neuron_index, spikes in spike_trains.items():
#        np.savetxt(file, spikes, header=f"Neuron {neuron_index}", fmt='%.18e')

#   with open('spike_trains.txt','w') as file:
#        for neuron_index, spikes in spike_trains.items():
#            spike_times_str = ' '.join(map(str, spikes/sim.ms))
#            file.write(f'{spike_times_str}\n')


        ISIs = {index: np.diff(spikes)
        for index, spikes in spMon.spike_trains().items()}
        synchrony=0.0
        for index, spikes in spMon.spike_trains().items():
            synchrony= synchrony+abs((spike_trains[0][-1]-spike_trains[index][-1])/sim.ms)
        isi=(ISIs[0][-1]/sim.ms)
        freq= 1000.0/(ISIs[0][-1]/sim.ms)
        freq2= 1000.0/(ISIs[0][-2]/sim.ms)
        print (synchrony)
        if (synchrony<0.7) and (abs(freq-freq2)<1.0):
            file.write(f"{dv:.1f} {freq:.3f} \n")
        print (dv,isi)
         


