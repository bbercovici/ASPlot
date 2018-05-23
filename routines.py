# MIT License

# Copyright (c) 2018 Benjamin Bercovici and Jay McMahon

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
        input = zip(foldernames, converts,savepaths,kepts,titles,log_scales)
        pool.map(plot, input)

'''
Generates results plots from files in provided directory
Inputs:
-------
- args : tuple of inputs

'''
def plot_results(args):

    inputpath,convert_to_RTN,savepath,kept,title,log_scale = args 

    epoch,stm,cov,deviations,ref_traj,mnvr = load_data(inputpath,only_covs = True,kept = kept)

    # If need be, state deviations and covariances are converted to the RTN frame
    if convert_to_RTN :
        cov,deviations = convert2RTN(cov,deviations,ref_traj)

    # Plot labels are created
    labels = create_labels(convert_to_RTN)

    # The state deviations are plotted along with the covariances
    plot_everything(epoch,labels,None,cov,mnvr,savepath,title,log_scale)






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

    epoch,stm,cov,deviations,ref_traj,mnvr = load_data(inputpath,only_covs = True,kept = kept)

    # If need be, state deviations and covariances are converted to the RTN frame
    if convert_to_RTN :
        cov,deviations = convert2RTN(cov,deviations,ref_traj)

    # Plot labels are created
    labels = create_labels(convert_to_RTN)

    # The state deviations are plotted along with the covariances
    plot_everything(epoch,labels,None,cov,mnvr,savepath,title,log_scale)


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

def plot_covariance_schedule_from_enclosing_folder(inputfolder,convert_to_RTN,log_scale = True,outputname = None,kept = - 1):
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


if sys.platform != "linux":
    plot_results("/Users/bbercovici/Desktop/AFSCN_DSN_8_HRS_R_4_RR_7/",convert_to_RTN = False,
        savepath = None,kept = -1,title  = "AFSCN_DSN_8_HRS_R_4_RR_7")






