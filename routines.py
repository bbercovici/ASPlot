# MIT License

# Copyright (c) 2018 Benjamin Bercovici and Andrew French

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.




import matplotlib
import sys
import os
import multiprocessing

# Switching backend on fortuna
print ("Platform : " + sys.platform)
if sys.platform == "linux":
    print ("Switching backend to 'agg'")

    matplotlib.use('Agg')

import numpy as np



import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.cm as cmaps
import os 

# files can be retrieved from 
# scp -r bebe0705@fortuna.colorado.edu:/home/anfr8485/FALCON/Filter/output/ outputs/


'''
Returns direction spanning input vector
Inputs:
------
x : np.array 
Outputs:
-------
u : unit-vector spanning x : x = |x| * u
'''
def normalize(x):
    if np.linalg.norm(x) < 1e-8:
        Raise(ValueError("Norm of provided vector is suspiciously small"))

    return x / np.linalg.norm(x)


'''
Loads the data from the specified directory
Inputs:
------
- inputpath: directory from where to pull data
- only_covs : if true, will only load epoch, cov and ref_traj
- kept : determines the number of measurements to kept between two consecutive kept measurements (kept == -1: no discarding)
Ouputs:
------
- epoch : vector of times (days)
- stm : vector of flattened STM s
- cov : vector of flattened covariances
- deviations : vector of state deviations
- ref_traj : vector of reference trajectory states
- mnvr : vector of maneuver epochs

'''
def load_data(inputpath,only_covs = False,kept = -1):

    epoch = np.loadtxt(inputpath + "epoch.txt")

    if kept > 0:
        kept_indices = np.linspace(0,len(epoch) -1 ,kept,dtype = int)
    else:
        kept_indices = np.linspace(0,len(epoch) -1 ,len(epoch),dtype = int)


    print("- Kept " + str(len(kept_indices)) + " observations from " + str(len(epoch)))

    epoch = np.copy(epoch[kept_indices])


    cov = np.loadtxt(inputpath + "cov.txt")[kept_indices,:]
    ref_traj = np.loadtxt(inputpath + "state.txt")[kept_indices,:]

    if not only_covs:
        stm = np.loadtxt(inputpath + "stm.txt")[kept_indices,:]
        deviations = np.loadtxt(inputpath + "dev.txt")[kept_indices,:]
    else:
        stm = None
        deviations = None

    mnvr = np.loadtxt(inputpath + "mnvr.txt")

    return epoch,stm,cov,deviations,ref_traj,mnvr



'''
Generates results plots from files in provided directory
Inputs:
-------
- args : tuple of inputs

'''
def plot_results_mlp(args):

    inputpath,convert_to_RTN,savepath,kept,title,log_scale = args 

    try:
        epoch,stm,cov,deviations,ref_traj,mnvr = load_data(inputpath,only_covs = False,kept = kept)

        covEME = cov
        # If need be, state deviations and covariances are converted to the RTN frame
        if convert_to_RTN :
            cov,deviations = convert2RTN(cov,deviations,ref_traj)

        # Plot labels are created
        labels = create_labels(convert_to_RTN)

        # The state deviations are plotted along with the covariances
        plot_everything(epoch,labels,None,cov,mnvr,savepath,title,log_scale)

        map_to_mnvrs(mnvr, epoch, covEME, stm, ref_traj, savepath)
    
    except:
        print("Filter did not converge")


def map_to_mnvrs(mnvr_times, epoch, cov, stm, ref, savepath):

    N = ref[0].shape[0]
    cov = [np.reshape(c, (N, N)) for c in cov]
    stm = [np.reshape(s, (N, N)) for s in stm]

    # DCO in days prior to burn
    # Given by Adv. Space
    mnvr_labels = ['TCM01', 'TCM02', 'TCM03', 'GSI01', 'GSI02', 'GSI03', 'GSI04']
    mnvr_dco_offset = [1., 3., 3., 1., 1., 1., 1.]

       # Mapping Covariance

    # Finding index of dco (finds latest epoch before dco)
    dco_time_index = []
    for i in range(len(mnvr_times)):
        dco_time = mnvr_times[i] - mnvr_dco_offset[i]
        epoch_minus = [x - dco_time for x in epoch]
        epoch_minus_temp_index = []
        
        for j in range(len(epoch_minus)):
            if epoch_minus[j] <= 0.:
                epoch_minus_temp_index.append(j)

        dco_time_index.append(epoch_minus_temp_index[-1])

    # Finding index of manuevers (finds latest epoch before manuever)
    mnvr_time_index = []
    for i in range(len(mnvr_times)):
        mnvr_time = mnvr_times[i]
        epoch_minus = [x - mnvr_time for x in epoch]
        epoch_minus_temp_index = []
        
        for j in range(len(epoch_minus)):
            if epoch_minus[j] <= 0.:
                epoch_minus_temp_index.append(j)

        mnvr_time_index.append(epoch_minus_temp_index[-1])

    mapFile = os.path.join(savepath, 'sigmas.map')

    # Mapping convariances
    with open(mapFile, 'w') as f:
        cov_mapped = []
        for i in range(len(dco_time_index)):
            dco_index = dco_time_index[i]
            mnvr_index = mnvr_time_index[i]

            # Creating STM from dco to mnvr
            stm_dco = stm[dco_index]
            stm_mnvr = stm[mnvr_index]
            stm_map = stm_mnvr @ np.linalg.inv(stm_dco)

            # Mapping covariance
            cov_map = stm_map @ cov[dco_index] @ stm_map.transpose()
            cov_map = np.reshape(cov_map, (1, N*N))

            # Rotating mapped covariance
            cov_map_rot, _ = convert2RTN(cov_map, None, np.reshape(ref[dco_index], (1, N)))
            cov_mapped.append(cov_map_rot)

            sigmas = np.sqrt(np.diag(np.reshape(cov_map_rot, (N, N))))
            sigmas = ['{:0.3e}'.format(s) for s in sigmas]

            line = [mnvr_labels[i]] + sigmas[:6]
            f.write(' & '.join(line) + ' \\\\\n')



'''
Generates results plots from files in provided directory
Inputs:
-------
- inputpath : path to directory where epoch.txt, cov.txt, state.txt , ref_traj.txt and stm.txt are
- convert_to_RTN : True if the provided state/covariances must be expressed in the RTN frame, False otherwise
- savepath : path to folder where to save generated plots, must be terminated by "/". If None provided then will just 
show plots
- kept : determines the number of measurements to kept between two consecutive kept measurements (kept == 1: no discarding)

'''
def plot_results(inputpath,convert_to_RTN = False,savepath = None,kept = -1,title = None,log_scale = False):

    epoch,stm,cov,deviations,ref_traj,mnvr = load_data(inputpath,only_covs = False,kept = kept)

    covEME = cov
    # If need be, state deviations and covariances are converted to the RTN frame
    if convert_to_RTN :
        cov,deviations = convert2RTN(cov,deviations,ref_traj)

    # Plot labels are created
    labels = create_labels(convert_to_RTN)

    # The state deviations are plotted along with the covariances
    plot_everything(epoch,labels,None,cov,mnvr,savepath,title,log_scale)

    map_to_mnvrs(mnvr, epoch, covEME, stm, ref_traj, savepath)


'''
Converts the provided deviations from the EME2000 frame to the RTN frame
Inputs:
------
- deviations : np.arrays of state errors (pxN) (pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,mass,Cr). Position/velocities in EME2000 (km,km/s)
- cov : np.arrays of covariances on state errors (pxN) (pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,mass,Cr). Position/velocities in EME2000 (km,km/s)
Outputs:
-------

- deviations_RTN : np.arrays of state errors (pxN) (pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,mass,Cr). Position/velocities in RTN frame (km,km/s)
- cov_RTN : np.arrays of covariances on state errors (pxN) (pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,mass,Cr). Position/velocities in RTN frame (km,km/s)

'''
def convert2RTN(cov,deviations,ref_traj):

    if deviations is not None:
        deviations_RTN = np.zeros(deviations.shape)
        N = deviations.shape[1] 

    else:
        deviations_RTN = None
        N = (int)(np.sqrt(cov.shape[1]))

    cov_RTN = np.zeros(cov.shape)


    for i in range (ref_traj.shape[0]):

        # The dcm between the inertial frame N and RTN frame R is assembled (RN):
        RN = np.zeros([3,3])
        e_r = normalize(ref_traj[i,0:3])
        e_n = normalize(np.cross(ref_traj[i,0:3] ,ref_traj[i,3:6] ))
        e_t = np.cross(e_n,e_r)

        RN[0,:] = e_r
        RN[1,:] = e_t
        RN[2,:] = e_n

        # The transformation matrix rotating the frame but leaving mass and Cr invariant is formed
        transform_mat = np.eye(N,N)
        transform_mat[0:3,0:3] = RN
        transform_mat[3:6,3:6] = RN

        # The deviations are transformed
        if deviations is not None:
            deviations_RTN[i,:] = transform_mat.dot(deviations[i,:])

        # The covariances are transformed
        cov_RTN[i,:] = (transform_mat.dot(np.reshape(cov[i,:],[N,N])).dot(transform_mat.T)).flatten()

    return cov_RTN,deviations_RTN


'''    
Plots the state errors at each of the provided epoch along with their covariance envelopes
Inputs:
-------
- epoch : np.arrays of epochs (px1)
- labels : dictionary of labels
- deviations : np.arrays of state errors (pxN) (pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,mass,Cr). Position/velocities in EME2000 (km,km/s)
- cov : np.arrays of covariances on state errors (pxN) (pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,mass,Cr). Position/velocities in EME2000 (km,km/s)
- savepath : path to folder where to save generated plots, must be terminated by "/". If None provided then will just 
show plots
- title : title to be added to each plot
- log_scale : true if results must be shown in logarithm space
'''
def plot_everything(epoch,labels,deviations = None,cov = None,mnvr = None, savepath = None,title = None,log_scale = False):

    # Shifting the epoch
    if mnvr is not None:
        mnvr = np.copy(mnvr - epoch[0])
    epoch = np.copy(epoch - epoch[0])

    # Extracting state dimension
    if deviations is not None:
        N = deviations.shape[1]
    elif cov is not None:
        N = (int)(np.sqrt(cov.shape[1]))
    else:
        Raise(ValueError("Need at minimum a vector of state deviations or flattened covariances"))

    # Extracting covariances bounds if provided 
    if cov is not None:
        sd_states = [np.sqrt(np.diag(np.reshape(cov[i,:],[N,N]))) for i in range(epoch.shape[0])]
        sd_states = np.vstack(sd_states)

    # Position
    if log_scale:
        if deviations is not None and cov is None:
            plt.semilogy(epoch,deviations[:,0],label = labels["e1"])
            plt.semilogy(epoch,deviations[:,1],label = labels["e2"])
            plt.semilogy(epoch,deviations[:,2],label = labels["e3"])
        elif cov is not None and deviations is None:
            plt.semilogy(epoch,3 * sd_states[:,0] ,".",label = labels["e1"])
            plt.semilogy(epoch,3 * sd_states[:,1] ,".",label = labels["e2"])
            plt.semilogy(epoch,3 * sd_states[:,2] ,".",label = labels["e3"])
            
        else:
            plt.semilogy(epoch,deviations[:,0],label = labels["e1"])
            plt.semilogy(epoch,deviations[:,1],label = labels["e2"])
            plt.semilogy(epoch,deviations[:,2],label = labels["e3"])
            plt.gca().set_color_cycle(None)

            plt.semilogy(epoch,3 * sd_states[:,0] ,".")
            plt.semilogy(epoch,3 * sd_states[:,1] ,".")
            plt.semilogy(epoch,3 * sd_states[:,2] ,".")
    else:
        if deviations is not None and cov is None:
            plt.plot(epoch,deviations[:,0],label = labels["e1"])
            plt.plot(epoch,deviations[:,1],label = labels["e2"])
            plt.plot(epoch,deviations[:,2],label = labels["e3"])
        elif cov is not None and deviations is None:
            plt.plot(epoch,3 * sd_states[:,0] ,".",label = labels["e1"])
            plt.plot(epoch,3 * sd_states[:,1] ,".",label = labels["e2"])
            plt.plot(epoch,3 * sd_states[:,2] ,".",label = labels["e3"])
            
        else:
            plt.plot(epoch,deviations[:,0],label = labels["e1"])
            plt.plot(epoch,deviations[:,1],label = labels["e2"])
            plt.plot(epoch,deviations[:,2],label = labels["e3"])
            plt.gca().set_color_cycle(None)

            plt.plot(epoch,3 * sd_states[:,0] ,".")
            plt.plot(epoch,3 * sd_states[:,1] ,".")
            plt.plot(epoch,3 * sd_states[:,2] ,".")

    if mnvr is not None:
        plot_maneuvers(mnvr)

    plt.xlabel("Days since Epoch")
    plt.ylabel("Position (km)")
    plt.title(title)
    plt.legend(loc = "center right",bbox_to_anchor = (1.25,0.5))
    # plt.tight_layout()
    if savepath is None:
        plt.show()
    else:
        plt.savefig(savepath + "positions.pdf",bbox_inches='tight')

    plt.clf()

    # Velocity
    if log_scale:
        if deviations is not None and cov is None:
            plt.semilogy(epoch,deviations[:,3],label = labels["e1"])
            plt.semilogy(epoch,deviations[:,4],label = labels["e2"])
            plt.semilogy(epoch,deviations[:,5],label = labels["e3"])
        elif cov is not None and deviations is None:
            plt.semilogy(epoch,3 * sd_states[:,3] ,".",label = labels["e1"])
            plt.semilogy(epoch,3 * sd_states[:,4] ,".",label = labels["e2"])
            plt.semilogy(epoch,3 * sd_states[:,5] ,".",label = labels["e3"])
            
        else:
            plt.semilogy(epoch,deviations[:,3],label = labels["e1"])
            plt.semilogy(epoch,deviations[:,4],label = labels["e2"])
            plt.semilogy(epoch,deviations[:,5],label = labels["e3"])
            plt.gca().set_color_cycle(None)

            plt.semilogy(epoch,3 * sd_states[:,3] ,".")
            plt.semilogy(epoch,3 * sd_states[:,4] ,".")
            plt.semilogy(epoch,3 * sd_states[:,5] ,".")
    else:
        if deviations is not None and cov is None:
            plt.plot(epoch,deviations[:,3],label = labels["e1"])
            plt.plot(epoch,deviations[:,4],label = labels["e2"])
            plt.plot(epoch,deviations[:,5],label = labels["e3"])
        elif cov is not None and deviations is None:
            plt.plot(epoch,3 * sd_states[:,3] ,".",label = labels["e1"])
            plt.plot(epoch,3 * sd_states[:,4] ,".",label = labels["e2"])
            plt.plot(epoch,3 * sd_states[:,5] ,".",label = labels["e3"])
        
        else:
            plt.plot(epoch,deviations[:,3],label = labels["e1"])
            plt.plot(epoch,deviations[:,4],label = labels["e2"])
            plt.plot(epoch,deviations[:,5],label = labels["e3"])
            plt.gca().set_color_cycle(None)

            plt.plot(epoch,3 * sd_states[:,3] ,".")
            plt.plot(epoch,3 * sd_states[:,4] ,".")
            plt.plot(epoch,3 * sd_states[:,5] ,".")
        

    
    if mnvr is not None:
        plot_maneuvers(mnvr)

    plt.xlabel("Days since Epoch")
    plt.ylabel("Velocity (km/s)")
    plt.title(title)

    plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))

    plt.legend(loc = "center right",bbox_to_anchor = (1.25,0.5))
    
    # plt.tight_layout()

    if savepath is None:
        plt.show()
    else:
        plt.savefig(savepath + "velocities.pdf",bbox_inches='tight')

    plt.clf()


    # # Mass
    # if deviations is not None :
    #     plt.plot(epoch,deviations[:,6])
    # if cov is not None:
    #     plt.gca().set_color_cycle(None)

    #     plt.plot(epoch,3 * sd_states[:,6] ,".")
    #     plt.gca().set_color_cycle(None)
    #     plt.plot(epoch,- 3 * sd_states[:,6] ,".")
    
    # plt.xlabel("Days since Epoch")
    # plt.ylabel("Mass (kg)")
    # plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
    # plt.legend(loc = "best")
    # plt.tight_layout()

    # if savepath is None:
    #     plt.show()
    # else:
    #     plt.savefig(savepath + "mass.pdf")

    # plt.clf()

    # # Cr
    # if deviations is not None :
    #     plt.plot(epoch,deviations[:,7])
    # if cov is not None:
    #     plt.gca().set_color_cycle(None)

    #     plt.plot(epoch,3 * sd_states[:,7] ,".")
    #     plt.gca().set_color_cycle(None)
    #     plt.plot(epoch,- 3 * sd_states[:,7] ,".")
   

    # plt.xlabel("Days since Epoch")
    # plt.ylabel("Cr (-)")
    # plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))

    # plt.legend(loc = "best")
    # plt.tight_layout()

    # if savepath is None:
    #     plt.show()
    # else:
    #     plt.savefig(savepath + "Cr.pdf")

    # plt.clf()


'''
Fills up a dictionary with labels for each plotted quantity
Inputs:
------
convert_to_RTN : True if position/velocity are expressed in the RTN frame, false otherwise
'''

def create_labels(convert_to_RTN):

    if convert_to_RTN: 
        labels = {"e1" : r"R","e2" : r"T","e3" : r"N","frame" : "RTN"}

    else:
        labels = {"e1" : r"X","e2" : r"Y","e3" : r"Z","frame" : "EME2000"}

        
    return labels

'''
Marks maneuvers epochs and labels them
Inputs:
-------
- mnvr : array of maneuver epochs
'''
def plot_maneuvers(mnvr):

    labels = ["TCM_1","TCM_2","TCM_3","GSI_1","GSI_2","GSI_3","GSI_4"]

    cmap =cmaps.get_cmap('viridis')

    for i in range(mnvr.shape[0]):
        plt.gca().axvline(mnvr[i],label = labels[i],linestyle = '--',color = cmap(float(i)/(mnvr.shape[0] -1 )))


'''
Plots covariance envelope 
for each state component given the provided
observation schedules (WORK IN PROGRESS)
Inputs:
-------
- inputfolder : path to enclosing folder for all cases
- convert_to_RTN : True if computed states must be converted to RTN
- log_scale : True if covariances must be plotted in semilogy scale
- kept : determines the number of measurements to kept between two consecutive kept measurements (kept == -1: no discarding)
'''
def plot_covariance_overlay_from_enclosing_folder(inputfolder,convert_to_RTN,log_scale = True,outputname = None,kept = - 1):
    covs = []
    cases = []
    ref_trajs = []
    epochs = []
    sd_states_all_cases = []

    labels =  create_labels(convert_to_RTN)

    # Creating the subplots
    ax_x_pos = plt.subplot(321)
    plt.ylabel(labels["e1"] + " position (km)")

    # Plot e2 pos component
    ax_y_pos = plt.subplot(323, sharex= ax_x_pos)
    plt.ylabel(labels["e2"] + " position (km)")

    # Plot e3 pos component
    ax_z_pos = plt.subplot(325, sharex= ax_x_pos)
    
    plt.xlabel("Days since Epoch")
    plt.ylabel(labels["e3"] + " position (km)")

    # Plot e1 vel component
    ax_x_vel = plt.subplot(322, sharex= ax_x_pos)
    plt.ylabel(labels["e1"] + " velocity (km/s)")
    
    # Plot e2 vel component
    ax_y_vel = plt.subplot(324, sharex= ax_x_pos)
    plt.ylabel(labels["e2"] + " velocity (km/s)")
    
    # Plot e3 vel component
    ax_z_vel = plt.subplot(326, sharex= ax_x_pos)
    plt.ylabel(labels["e3"] + " velocity (km/s)")

    plt.xlabel("Days since Epoch")

    # Loading and plotting
    for folder in os.walk(inputfolder) :
        for subfolder in folder[1]:
            foldername = folder[0] + subfolder + "/"
            print("Loading case " + subfolder)

            # Loading
            epoch,stm,cov,deviations,ref_traj,mnvr = load_data(foldername,only_covs = True,kept = kept)
            epoch = epoch - epoch[0]


            # Converting
            if convert_to_RTN:
                dummy = None
                cov,dummy = convert2RTN(cov,dummy,ref_traj)

            # Extracting SDs
            N = (int)(np.sqrt(cov.shape[1]))
            sd_states = [np.sqrt(np.diag(np.reshape(cov[i,:],[N,N]))) for i in range(epoch.shape[0])]
            sd_states = np.vstack(sd_states)

            # Plot e1 pos component
            if log_scale:
                ax_x_pos.semilogy(epoch,3 * sd_states[:,0],'.')
            else:
                ax_x_pos.plot(epoch,3 * sd_states[:,0],'.')


            # Plot e2 pos component
            if log_scale:
                ax_y_pos.semilogy(epoch,3 * sd_states[:,1],'.')
            else:
                ax_y_pos.plot(epoch,3 * sd_states[:,1],'.')

            # Plot e3 pos component
            if log_scale:
                ax_z_pos.semilogy(epoch,3 * sd_states[:,2],'.')
            else:
                ax_z_pos.plot(epoch,3 * sd_states[:,2],'.')

           
            # Plot e1 vel component
            if log_scale:
                ax_x_vel.semilogy(epoch,3 * sd_states[:,3],'.')
            else:
                ax_x_vel.plot(epoch,3 * sd_states[:,3],'.')
            
            # Plot e2 vel component
            if log_scale:
                ax_y_vel.semilogy(epoch,3 * sd_states[:,4],'.')
            else:
                ax_y_vel.plot(epoch,3 * sd_states[:,4],'.')
            
            # Plot e3 vel component
            if log_scale:
                ax_z_vel.semilogy(epoch,3 * sd_states[:,5],'.')
            else:
                ax_z_vel.plot(epoch,3 * sd_states[:,5],'.')


    if outputname is not None:
        plt.savefig(outputname)
    else:
        plt.show()





'''
Generates the covariance overlays for the required combination of cases
'''
def plot_covariance_overlays():

    #########################
    # Hours-per-day tradoff #
    #########################

    source_folder = "/home/anfr8485/FALCON/Filter/test/advspc_results/"

    folder_list_AFSCN_ONLY  = [
    source_folder + "AFSCN_ONLY_24_HRS_R_3_RR_6_PN_12/",
    source_folder + "AFSCN_ONLY_8_HRS_R_3_RR_6_PN_12/",
    source_folder + "AFSCN_ONLY_4_HRS_R_3_RR_6_PN_12/"] 

    plot_covariance_overlay_from_list_of_folder(
    folder_list_AFSCN_ONLY,
    source_folder,
    True,
    kept = 4000,
    outputname = "AFSCN_ONLY",
    outputfolder = source_folder) 

    folder_list_AFSCN_DSN  = [
    source_folder + "AFSCN_DSN_24_HRS_R_3_RR_6_PN_12/",
    source_folder + "AFSCN_DSN_8_HRS_R_3_RR_6_PN_12/",
    source_folder + "AFSCN_DSN_4_HRS_R_3_RR_6_PN_12/"]

    plot_covariance_overlay_from_list_of_folder(
    folder_list_AFSCN_DSN,
    source_folder,
    True,
    kept = 4000,
    outputname = "AFSCN_DSN",
    outputfolder = source_folder)


    ###################
    # R Noise tradoff #
    ###################

    folder_list_AFSCN_DSN_8_HRS_R_NOISE_LVL  = [
    source_folder + "AFSCN_DSN_8_HRS_R_5_RR_6_PN_12/",
    source_folder + "AFSCN_DSN_8_HRS_R_4_RR_6_PN_12/",
    source_folder + "AFSCN_DSN_8_HRS_R_3_RR_6_PN_12/",
    source_folder + "AFSCN_DSN_8_HRS_R_2_RR_6_PN_12/"]  

    plot_covariance_overlay_from_list_of_folder(
    folder_list_AFSCN_DSN_8_HRS_R_NOISE_LVL,
    source_folder,
    True,
    kept = 4000,
    outputname = "AFSCN_DSN_8_HRS_R_NOISE_LVL",
    outputfolder = source_folder)


    folder_list_AFSCN_ONLY_8_HRS_R_NOISE_LVL  = [
    source_folder + "AFSCN_ONLY_8_HRS_R_5_RR_6_PN_12/",
    source_folder + "AFSCN_ONLY_8_HRS_R_4_RR_6_PN_12/",
    source_folder + "AFSCN_ONLY_8_HRS_R_3_RR_6_PN_12/",
    source_folder + "AFSCN_ONLY_8_HRS_R_2_RR_6_PN_12/"]  

    plot_covariance_overlay_from_list_of_folder(
    folder_list_AFSCN_ONLY_8_HRS_R_NOISE_LVL,
    source_folder,
    True,
    kept = 4000,
    outputname = "AFSCN_ONLY_8_HRS_R_NOISE_LVL",
    outputfolder = source_folder)


    ####################
    # RR Noise tradoff #
    ####################

    folder_list_AFSCN_DSN_8_HRS_RR_NOISE_LVL  = [
    source_folder + "AFSCN_DSN_8_HRS_R_3_RR_8_PN_12/",
    source_folder + "AFSCN_DSN_8_HRS_R_3_RR_7_PN_12/",
    source_folder + "AFSCN_DSN_8_HRS_R_3_RR_6_PN_12/",
    source_folder + "AFSCN_DSN_8_HRS_R_3_RR_5_PN_12/"]  

    plot_covariance_overlay_from_list_of_folder(
    folder_list_AFSCN_DSN_8_HRS_RR_NOISE_LVL,
    source_folder,
    True,
    kept = 4000,
    outputname = "AFSCN_DSN_8_HRS_RR_NOISE_LVL",
    outputfolder = source_folder)

    folder_list_AFSCN_8_HRS_ONLY_RR_NOISE_LVL  = [
    source_folder + "AFSCN_ONLY_8_HRS_R_3_RR_8_PN_12/",
    source_folder + "AFSCN_ONLY_8_HRS_R_3_RR_7_PN_12/",
    source_folder + "AFSCN_ONLY_8_HRS_R_3_RR_6_PN_12/",
    source_folder + "AFSCN_ONLY_8_HRS_R_3_RR_5_PN_12/"]  

    plot_covariance_overlay_from_list_of_folder(
    folder_list_AFSCN_8_HRS_ONLY_RR_NOISE_LVL,
    source_folder,
    True,
    kept = 4000,
    outputname = "AFSCN_8_HRS_ONLY_RR_NOISE_LVL",
    outputfolder = source_folder)


    ####################
    # Process Noise tradoff #
    ####################

    folder_list_AFSCN_DSN_8_HRS_PN_LVL  = [
    source_folder + "AFSCN_DSN_8_HRS_R_3_RR_6_PN_12/",
    source_folder + "AFSCN_DSN_8_HRS_R_3_RR_6_PN_10/",
    source_folder + "AFSCN_DSN_8_HRS_R_3_RR_6_PN_8/",
    source_folder + "AFSCN_DSN_8_HRS_R_3_RR_6_PN_6/"]

    plot_covariance_overlay_from_list_of_folder(
    folder_list_AFSCN_DSN_8_HRS_PN_LVL,
    source_folder,
    True,
    kept = 4000,
    outputname = "AFSCN_DSN_8_HRS_PN_LVL",
    outputfolder = source_folder)

    folder_list_AFSCN_8_HRS_ONLY_PN_LVL  = [
    source_folder + "AFSCN_ONLY_8_HRS_R_3_RR_6_PN_12/",
    source_folder + "AFSCN_ONLY_8_HRS_R_3_RR_6_PN_10/",
    source_folder + "AFSCN_ONLY_8_HRS_R_3_RR_6_PN_8/",
    source_folder + "AFSCN_ONLY_8_HRS_R_3_RR_6_PN_6/"]  

    plot_covariance_overlay_from_list_of_folder(
    folder_list_AFSCN_8_HRS_ONLY_PN_LVL,
    source_folder,
    True,
    kept = 4000,
    outputname = "AFSCN_8_HRS_ONLY_PN_LVL",
    outputfolder = source_folder)


'''
Plots covariance envelope 
for each state component given the provided
observation schedules (WORK IN PROGRESS)
Inputs:
-------
- folder_list : list of enclosing folders we want to have overlaid
- inputfolder : path to enclosing folder for all cases
- convert_to_RTN : True if computed states must be converted to RTN
- log_scale : True if covariances must be plotted in semilogy scale
- kept : determines the number of measurements to kept between two consecutive kept measurements (kept == -1: no discarding)
- plot_RSS : if True, will only plot position and velocity RSS
- outputfolder : folder where to save figure
'''
def plot_covariance_overlay_from_list_of_folder(folder_list,
    inputfolder,
    convert_to_RTN,
    log_scale = True,
    outputname = None,
    kept = - 1,
    plot_RSS = True,
    outputfolder = None):
    covs = []
    cases = []
    ref_trajs = []
    epochs = []
    sd_states_all_cases = []

    labels = create_labels(convert_to_RTN)
    plt.figure()
    plt.suptitle(r"3$\sigma$ RSS for " + outputname)

    if plot_RSS:


        # Creating the subplots
        ax_pos_RSS = plt.subplot(211)
        plt.ylabel("Position RSS (km)")

        # Plot e3 vel component
        ax_vel_RSS = plt.subplot(212, sharex= ax_pos_RSS)
        plt.ylabel("Velocity RSS (km/s)")
        plt.xlabel("Days since Epoch")

    else:
        # Creating the subplots
        ax_x_pos = plt.subplot(321)
        plt.ylabel(labels["e1"] + " position (km)")

        # Plot e2 pos component
        ax_y_pos = plt.subplot(323, sharex= ax_x_pos)
        plt.ylabel(labels["e2"] + " position (km)")

        # Plot e3 pos component
        ax_z_pos = plt.subplot(325, sharex= ax_x_pos)
        
        plt.xlabel("Days since Epoch")
        plt.ylabel(labels["e3"] + " position (km)")

        # Plot e1 vel component
        ax_x_vel = plt.subplot(322, sharex= ax_x_pos)
        plt.ylabel(labels["e1"] + " velocity (km/s)")
        
        # Plot e2 vel component
        ax_y_vel = plt.subplot(324, sharex= ax_x_pos)
        plt.ylabel(labels["e2"] + " velocity (km/s)")
        
        # Plot e3 vel component
        ax_z_vel = plt.subplot(326, sharex= ax_x_pos)
        plt.ylabel(labels["e3"] + " velocity (km/s)")

        plt.xlabel("Days since Epoch")

    # Loading and plotting
    ncol = 0


    for foldername_from_list in folder_list:

        for folder in os.walk(inputfolder) :
            for subfolder in folder[1]:

                foldername = folder[0] + subfolder + "/"

                if foldername != foldername_from_list:
                    continue

                print("Loading case " + subfolder)

                # Loading
                epoch,stm,cov,deviations,ref_traj,mnvr = load_data(foldername,only_covs = True,kept = kept)
                epoch = epoch - epoch[0]

                # Converting
                if convert_to_RTN:
                    dummy = None
                    cov,dummy = convert2RTN(cov,dummy,ref_traj)

                # Extracting SDs
                N = (int)(np.sqrt(cov.shape[1]))

                if plot_RSS:

                    sd_states_pos = np.array([np.sqrt(np.trace(np.reshape(cov[i,:],[N,N])[0:3,0:3])) for i in range(epoch.shape[0])])
                    sd_states_vel = np.array([np.sqrt(np.trace(np.reshape(cov[i,:],[N,N])[3:6,3:6])) for i in range(epoch.shape[0])])

                    # Extracting legend
                    if outputname == "AFSCN_ONLY" or outputname == "AFSCN_DSN":
                        if "_4_HRS" in foldername:
                            label = "4 hours"
                        elif "8_HRS" in foldername:
                            label = "8 hours"
                        elif "24_HRS" in foldername:
                            label = "24 hours"
                        ncol = ncol + 1

                    elif outputname == "AFSCN_DSN_8_HRS_R_NOISE_LVL" or outputname == "AFSCN_ONLY_8_HRS_R_NOISE_LVL":
                        if "_R_5" in foldername:
                            label = r"$10^{-5}\ \mathrm{km}$"
                        elif "_R_4" in foldername:
                            label = r"$10^{-4}\ \mathrm{km}$"
                        elif "_R_3" in foldername:
                            label = r"$10^{-3}\ \mathrm{km}$"
                        elif "_R_2" in foldername:
                            label = r"$10^{-2}\ \mathrm{km}$"
                        ncol = ncol + 1


                    elif outputname == "AFSCN_DSN_8_HRS_RR_NOISE_LVL" or outputname == "AFSCN_8_HRS_ONLY_RR_NOISE_LVL":
                        if "RR_8" in foldername:
                            label = r"$10^{-8}\ \mathrm{km/s}$"
                        elif "RR_7" in foldername:
                            label = r"$10^{-7}\ \mathrm{km/s}$"
                        elif "RR_6" in foldername:
                            label = r"$10^{-6}\ \mathrm{km/s}$"
                        elif "RR_5" in foldername:
                            label = r"$10^{-5}\ \mathrm{km/s}$"
                        ncol = ncol + 1


                    elif outputname == "AFSCN_DSN_8_HRS_PN_LVL" or outputname == "AFSCN_8_HRS_ONLY_PN_LVL":
                        if "PN_12" in foldername:
                            label = r"$10^{-12}\ \mathrm{km/s^2}$"
                        elif "PN_10" in foldername:
                            label = r"$10^{-10}\ \mathrm{km/s^2}$"
                        elif "PN_8" in foldername:
                            label = r"$10^{-8}\ \mathrm{km/s^2}$"
                        elif "PN_6" in foldername:
                            label = r"$10^{-6}\ \mathrm{km/s^2}$"
                        ncol = ncol + 1



                    # Plot e1 pos component
                    if log_scale:
                        ax_pos_RSS.semilogy(epoch,3 * sd_states_pos,'.',label = label)
                    else:
                        ax_pos_RSS.plot(epoch,3 * sd_states_pos,'.',label = label)




                    if log_scale:
                        ax_vel_RSS.semilogy(epoch,3 * sd_states_vel,'.')
                    else:
                        ax_vel_RSS.plot(epoch,3 * sd_states_vel,'.')

            
                else:

                    sd_states = [np.sqrt(np.diag(np.reshape(cov[i,:],[N,N]))) for i in range(epoch.shape[0])]
                    sd_states = np.vstack(sd_states)

                    # Plot e1 pos component
                    if log_scale:
                        ax_x_pos.semilogy(epoch,3 * sd_states[:,0],'.')
                    else:
                        ax_x_pos.plot(epoch,3 * sd_states[:,0],'.')

                    # Plot e2 pos component
                    if log_scale:
                        ax_y_pos.semilogy(epoch,3 * sd_states[:,1],'.')
                    else:
                        ax_y_pos.plot(epoch,3 * sd_states[:,1],'.')

                    # Plot e3 pos component
                    if log_scale:
                        ax_z_pos.semilogy(epoch,3 * sd_states[:,2],'.')
                    else:
                        ax_z_pos.plot(epoch,3 * sd_states[:,2],'.')

                    # Plot e1 vel component
                    if log_scale:
                        ax_x_vel.semilogy(epoch,3 * sd_states[:,3],'.')
                    else:
                        ax_x_vel.plot(epoch,3 * sd_states[:,3],'.')
                    
                    # Plot e2 vel component
                    if log_scale:
                        ax_y_vel.semilogy(epoch,3 * sd_states[:,4],'.')
                    else:
                        ax_y_vel.plot(epoch,3 * sd_states[:,4],'.')
                    
                    # Plot e3 vel component
                    if log_scale:
                        ax_z_vel.semilogy(epoch,3 * sd_states[:,5],'.')
                    else:
                        ax_z_vel.plot(epoch,3 * sd_states[:,5],'.')



    ax_pos_RSS.legend(loc = "upper center",bbox_to_anchor = (0.5,1.2),ncol = ncol)  
             
    if outputname is not None:
        plt.savefig(outputfolder + outputname + ".pdf",bbox_inches='tight')
    else:
        plt.show()



'''
Generates results plots for all simulation cases in the enclosing folder
Inputs:
-------
- inputpath : path to directory where epoch.txt, cov.txt, state.txt , ref_traj.txt and stm.txt are
- convert_to_RTN : True if the provided state/covariances must be expressed in the RTN frame, False otherwise
- savepath : path to folder where to save generated plots, must be terminated by "/". If None provided then will just 
show plots
- kept : determines the number of measurements to kept between two consecutive kept measurements (kept == -1: no discarding)
- log_scale : true if results must be shown in logarithm space
'''
def plot_results_from_enclosing_folder(inputfolder,convert_to_RTN = False,savepath = None,kept = -1,log_scale = True,mltpro = False):


    inputfolders = []
    savepaths = []
    converts = []
    kepts = []
    titles = []
    log_scales = []

    for folder in os.walk(inputfolder) :
        for subfolder in folder[1]:
            foldername = folder[0] + subfolder + "/"
            if mltpro is False:
                print("Loading case " + subfolder)
                plot_results(foldername,convert_to_RTN = convert_to_RTN,savepath = foldername,kept = kept,title = subfolder,log_scale = log_scale)
                
            else:
                inputfolders += [foldername]
                converts += [convert_to_RTN]
                savepaths += [foldername]
                kepts += [kept]
                titles += [subfolder]
                log_scales += [log_scale]





    if mltpro:
        pool = multiprocessing.Pool()
        input = zip(inputfolders, converts,savepaths,kepts,titles,log_scales)
        pool.map(plot_results_mlp, input)



if sys.platform != "linux":
    plot_results("/Users/bbercovici/Desktop/AFSCN_DSN_8_HRS_R_4_RR_7/",convert_to_RTN = False,
        savepath = None,kept = -1,title  = "AFSCN_DSN_8_HRS_R_4_RR_7")






