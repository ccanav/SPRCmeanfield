import FSIN as sim
import numpy as np
import matplotlib.pyplot as plt
import sys

method = "euler";
dt = 5.e-4;
sim_time=1000.;

gms=1.65

kneu = 43

dv=4.9  
rates, spMon = sim.gewnet(kneu=kneu,sim_time=sim_time,return_state=False,Nrec=100,dt=dt,method=method,act_syns=True,mod_gL=False, gjs=False,gext=7.,sin_cond=False,const_cond=True,tau_rise=.3,homo_neurons=True,gms=gms,dist_syns=False,read_delays=False, read_gms=False,dt_rec=.1,dmin=dv,dmax=dv,connectivity="FID",return_syns=False,read_syns=False,USEm=1.,std=False, Es=-75.,randics=True,one_perturbed=False); 
spike_trains = spMon.spike_trains();
#with open('output.txt', 'w') as file:
#   for neuron_index, spikes in spike_trains.items():
#        np.savetxt(file, spikes, header=f"Neuron {neuron_index}", fmt='%.18e')

with open('spike_trains.txt','w') as file:
    for neuron_index, spikes in spMon.spike_trains().items():
         spike_times_str = ' '.join(map(str, spikes/sim.ms))
         file.write(f'{spike_times_str}\n')


ISIs = {index: np.diff(spikes)
for index, spikes in spMon.spike_trains().items()}
synchrony=0.0
for index, spikes in spMon.spike_trains().items():
   synchrony= synchrony+abs((spike_trains[0][-1]-spike_trains[index][-1])/sim.ms)

freq= 1.0/(ISIs[0][-1])
freq2= 1.0/(ISIs[0][-2])
print (synchrony)
print (freq)
print (freq2)


width = 5.2 # in inches
height = 7.# in inches
tmp = plt.figure(figsize=[width,height]);

tmp1 = plt.subplot(3,1,1); tmp = [ plt.plot(spMon.spike_trains()[kneu]/sim.ms, kneu*np.ones(len(spMon.spike_trains()[kneu])), "ok", ms=1.) for kneu in range(20) ]; tmp = plt.xlim(900.,1000.);

ISIses_heteact = sim.b2.diff(spMon.spike_trains()[0]/sim.ms);

print(ISIses_heteact[-1]-ISIses_heteact[-2])

#plt.title("B. $E_{syn}$=-75, $\delta$=0.8 ms, random")
tmp1.spines['top'].set_visible(False)
tmp1.spines['right'].set_visible(False)
tmp1.spines['bottom'].set_visible(False)
tmp1.spines['left'].set_visible(False)
tmp1.set_yticks([])
tmp1.set_yticklabels([])
plt.savefig("Fig5D.png");
plt.savefig("Fig5d.eps"); plt.clf();

         


