import sys
sys.path.append('/home/jra21/work/source/repos/opep')

from openep import io as openep_io
from openep import case as openep_case
from openep import mesh_routines as openep_mesh
from openep import case_routines as case_routines
from openep import draw_routines as draw
import numpy as np
import matplotlib.pyplot as plt



filename = '/home/jra21/work/source/repos/opep/examples/data/new_dataset_2.mat'

ep_case = openep_io.load_case(filename)

# DrawVoltage Map
# hsurf = draw.DrawMap(ep_case,
#                     volt='bip',
#                     freeboundary_color='black',
#                     freeboundary_width=5,
#                     cmap='jet_r',
#                     minval=0,
#                     maxval=2,
#                     volt_below_color='brown', 
#                     volt_above_color='magenta', 
#                     nan_color='gray',
#                     plot=True)

# GetAnatomicalStructure and DrawFreeBoundaries on mesh
# hsurf1 = draw.getAnatomicalStructures(ep_case,plot=True)

def get_egms_at_points(mesh_case,*args):
    '''
    Access electrogram from stored in the openep data format
    
    ''' 
    n_standard_args = 1
    i_egm = np.array(':')
    n_arg_in = len(args)+1

    electric_egm_bip = []
    electric_egm_uni = []
    egm_traces = []

    if n_arg_in>n_standard_args:
        for i in range(0,(n_arg_in-n_standard_args),2):
            if np.char.lower(args[i]) == 'iegm':
                i_egm = args[i+1]

    if (type(i_egm)==str) and (np.char.compare_chararrays(i_egm,':','==',True)):
        electric_egm_bip = mesh_case.electric["egm"].T
        electric_egm_uni = mesh_case.electric["egmUni"].T

    else:
        electric_egm_bip_raw = mesh_case.electric["egm"].T
        electric_egm_bip_raw = electric_egm_bip_raw.astype(float)
        electric_egm_uni_raw = mesh_case.electric["egmUni"].T
        electric_egm_uni_raw = electric_egm_uni_raw.astype(float)
        
        if len(i_egm)==2:
            start_indx = i_egm[0]
            end_indx = i_egm[1]+1

            if (end_indx > len(electric_egm_bip_raw)):
                raise IndexError('egm index out of range')

            else:
                for i in range(start_indx,end_indx):
                    electric_egm_bip.append(electric_egm_bip_raw[i,:])
                    electric_egm_uni.append(electric_egm_uni_raw[i,:])            

        else:
            if (i_egm[0] > (len(electric_egm_bip_raw)-1)):
                raise IndexError('egm index out of range')

            else:
                electric_egm_bip.append(electric_egm_bip_raw[i_egm])
                electric_egm_uni.append(electric_egm_uni_raw[i_egm])
                    
    electric_egm_bip = np.array(electric_egm_bip).astype(float)
    electric_egm_uni = np.array(electric_egm_uni).astype(float)

    
    egm_traces.append(electric_egm_bip)
    egm_traces.append(electric_egm_uni[0][0][:,0])
    egm_traces.append(electric_egm_uni[0][0][:,1])

    return {'egm_traces':egm_traces,'electric_egm_bip':electric_egm_bip}


egm = get_egms_at_points(ep_case,"iegm",[0])

egm_bip = egm["electric_egm_bip"][0]
size_of_egm = egm_bip.shape[1]

t = np.arange(0,size_of_egm,1)
# t_lin = np.linspace(0,size_of_egm,num=size_of_egm)
print("t-shape\n",t.shape)
print("egm_bip\n",egm_bip[0])

plt.plot(t,egm_bip[0])
plt.ylabel("egm")
plt.show()