import ROOT
from ROOT import TH1D, TH2D
from ROOT import TFile, TTree
import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#root file 
f  = ROOT.TFile("mcv6.4.gsg_anumu-CCHEDIS_1e2-1e8GeV.sirene.jterbr.jchain.aashower.dst.merged_9635_10005_inline.root")


#-------Comparing run ids with events and mc_run_id----------------------
test_list1=[]
test_list2=[]
test_list3=[]
with uproot.open("mcv6.4.gsg_anumu-CCHEDIS_1e2-1e8GeV.sirene.jterbr.jchain.aashower.dst.merged_9635_10005_inline.root")["E"] as ifile:
   x = ifile['run_id']
   test_list1.append(x.arrays(library='pd'))
   y = ifile['id']
   test_list2.append(y.arrays(library = 'pd'))
   z = ifile['mc_run_id']
   test_list3.append(z.arrays(library = 'pd'))

#RUN_ID
run_id = []
for i in range(len(test_list1[0]['run_id'])):
   run_id.append(test_list1[0]['run_id'][i])
#print(run_id)
print(len(run_id))

#EVENTS
events = []
for i in range(len(test_list2[0]['id'])):
   events.append(test_list2[0]['id'][i])
#print(events)
print(len(events))

#MC_RUN_ID
mc_run_id = []
for i in range(len(test_list3[0]['mc_run_id'])):
   mc_run_id.append(test_list3[0]['mc_run_id'][i])
#print(mc_run_id)
print(len(mc_run_id))

#We have created the lists run_id, events and mc_run_id containing the run ids, the events and the mc run ids. After checking that the length of the lists run_id, events and mc_run_id is the same, we can proceed:

#-------------------------------------

#Read the branches that  we will need and fill the lists:

test_list_track_dir_z = [] #for zenith angle
test_list_likelihood = [] #for likelihood

with uproot.open("mcv6.4.gsg_anumu-CCHEDIS_1e2-1e8GeV.sirene.jterbr.jchain.aashower.dst.merged_9635_10005_inline.root")["E"] as ifile:
   dirz  = ifile['mc_trks.dir.z']
   lkld = ifile['trks.lik']
   test_list_track_dir_z.append(dirz.arrays(library='pd'))
   test_list_likelihood.append(lkld.arrays(library = 'pd'))

#-------track_dir_z-------
track_dir_z = []
for i in range(len(test_list_track_dir_z[0]['mc_trks.dir.z'])):
   track_dir_z.append(test_list_track_dir_z[0]['mc_trks.dir.z'][i][0])
#print(track_dir_z)

#-------likelihood------
likelihood = []
for i in range(len(test_list_likelihood[0]['trks.lik'])):
   likelihood.append(test_list_likelihood[0]['trks.lik'][i][0])
#print(likelihood)


#transform lists into dataframes
#df_run_id = pd.DataFrame(run_id , columns = ['Run Id']) #for run ids
#print(df_run_id)

#df_event = pd.DataFrame(events, columns = ['Events']) #for events
#print(df_event)

#df_mc_run_id = pd.DataFrame(mc_run_id, columns = ['Mc Run Id'])
#print(df_mc_run_ids)


# The following loop will give the likelihoods, track.dir.z and whatever else we need:

for i in range(len(run_id)):
   if run_id[i] == 9986 and events[i] == 1 and mc_run_id[i] == 254:
      print(likelihood[i])
      print(track_dir_z)
   else:
      quit 
   

 
#------------------------------------------------
#-------------Histograms-------------------
#run id
plt.hist(run_id, bins = 100)
plt.show()

#events
plt.hist(events, bins = 100)
plt.show()

#mc run id
plt.hist(mc_run_id, bins = 100)
plt.show()

#track_dir_z
plt.hist(track_dir_z, bins = 100)
plt.show()

#likelihood
plt.hist(likelihood, bins = 100)
plt.show()


                           








  









