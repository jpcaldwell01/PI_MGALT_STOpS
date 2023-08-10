# FORM: var = clusterPruning(dates,opts,var,alg,num_mig,I)
#
# |-----------------------------------------------------------------------
# |
# | NOTES:
# |     -This function is the Cluster Pruning function. It uses the DBSCAN
# |      Algorithm to search for clusters in the Departure Date and ToF 
# |      Windows. It is N-Dimensional and works with any number of transfers.
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -dates      (npop*ngen,len(bod['bodies'])+1)       [array]        [unitless]
# |         An array containing the Departure Date, Times of Flight, and Flyby Date of each population member
# |     -opts               (17)       [dict]        [unitless]
# |         A dict containing the mission options
# |     -var                (10)       [dict]        [unitless]
# |         A dict containing the variable limits
# |     -alg                (1,1)      [string]      [unitless]
# |         Island name to perform cluster pruning ('DE','GA','PSO','MBH')
# |     -num_mig            (1,1)       [int]           [unitless]
# |         The current migration number
# |     -I           (len(bod['bodies'])+1)       [list]           [unitless]
# |         Cluster Pruning Indexes
# |     -plot_me            (1,1)       [int]        [unitless]
# |         Variable determining whether or not to plot the data. A value of 1 produces plots.
# |         
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |     -var          (11)       [dict]        [unitless]
# |         A dict containing the variable limits (with cluster pruning limits changed)
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |     Plotting capabilities do exist. However, they only exist for 2 transfer trajectories.
# |
# |-----------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import DBSCAN
import pylab as p
import time
from IPython import get_ipython
from julian import from_jd


def clusterPruning(dates,opts,var,alg,num_mig,I,plot_me):
    
    arr = dates[:,-1]+dates[:,-2]
    dates = np.hstack((dates,arr[:,None]))
    clustMig = 1
     
    if num_mig == clustMig and alg == 'DE' and opts['cluster_Pruning'] == 'y':
           
        if plot_me == 1:
            # Plot Dep and ToF
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            boundingBoxRect(var['low'][I],var['high'][I],'k--',opts,ax)
            boundingBoxRect(var['lowC'][I],var['highC'][I],'b--',opts,ax)
            Legend = {}



        # Cluster Pruning
        get_ipython().magic('clear')
        print("Cluster Pruning...")
        xL = np.zeros((np.shape(dates)[1]-var['transfers'],))
        xU = np.zeros((np.shape(dates)[1]-var['transfers'],))
        
        
        clustering = DBSCAN(eps=7.5, min_samples=10,n_jobs=-1).fit(dates[:,:-var['transfers']])
        DBSCAN_dataset = dates.copy()
        DBSCAN_dataset = np.hstack((DBSCAN_dataset,clustering.labels_[:,None]))
        count = 0
        sigData = np.empty((0,np.shape(dates)[1]))
        for it in np.unique(DBSCAN_dataset[:,-1]):
            if np.shape(DBSCAN_dataset[DBSCAN_dataset[:,-1]==it])[0] >= .1*np.shape(dates)[0] and it != -1:
                if plot_me == 1:
                    ax.scatter(DBSCAN_dataset[DBSCAN_dataset[:,-1]==it][:,0],DBSCAN_dataset[DBSCAN_dataset[:,-1]==it][:,1],DBSCAN_dataset[DBSCAN_dataset[:,-1]==it][:,2],s=5)
                    Legend[count] = str(it)
                count = count+1
                sigData = np.append(sigData,DBSCAN_dataset[DBSCAN_dataset[:,-1]==it][:,:-1],axis=0)
    
            else:
                if plot_me == 1:
                    ax.scatter(DBSCAN_dataset[DBSCAN_dataset[:,-1]==it][:,0],DBSCAN_dataset[DBSCAN_dataset[:,-1]==it][:,1],DBSCAN_dataset[DBSCAN_dataset[:,-1]==it][:,2],c='k',s=5,label='_nolegend_')
        
        if len(sigData) > 0:
            xL = np.zeros((np.shape(dates)[1]-var['transfers'],))
            xU = np.zeros((np.shape(dates)[1]-var['transfers'],))
            for it in range(np.shape(dates)[1]-var['transfers']):
                xL[it] = np.min(sigData[:,it])
                xU[it] = np.max(sigData[:,it])
            
                I2 = I[:var['transfers']+1]
                var['lowC'][I2] = xL
                var['highC'][I2] = xU
            if plot_me == 1:
                boundingBoxRect(xL,xU,'g--',opts,ax)
                plt.legend(Legend)
                
            print('Cluster pruning was a success!')
            if count == 1:
                print('1 cluster was found containing more than 10 percent of the data. Using this cluster, the new departure and time of flight bounds are:')
            else:
                print('%i clusters were found each containing more than 10 percent of the data. Using these clusters, the new departure and time of flight bounds are:' % count)
            
            var['clusterResults'] = np.zeros((var['transfers']+2,2))
            
            date1 = from_jd(xL[0])
            date2 = from_jd(xU[0])
            var['clusterResults'][0,:] = [xL[0],xU[0]]
            print('     Dep Date: ',date1.strftime("%d/%m/%Y"),' --> ',date2.strftime("%d/%m/%Y"))
            
            for it2 in range(len(xL)-1):
                print('     ToF %i: %.2f --> %.2f days' % (it2+1, xL[it2+1],xU[it2+1]))
                var['clusterResults'][it2+1,:] = [xL[it2+1],xU[it2+1]]

            print('This information will also be shown at the end of optimization.')
            print('Pausing for 15 seconds...')
            
            
            
            
        else:
            print("No clusters found. Returning to optimization in 15 seconds.")
            var['clusterResults'] = []
        time.sleep(15)
        
        if plot_me == 1:

            ax.scatter(dates[:,0],dates[:,1],dates[:,2],c='k',s=5,label='_nolegend_')
        
            ax.set_xlabel('Departure Date')
            ax.set_zlabel('ToF 2')
            ax.set_ylabel('ToF 1')
            plt.title(alg + ' Migration %i, Dep & ToF Scatter Plot' % num_mig)
            plt.tight_layout()
            plt.show()
            
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            boundingBoxDiag(var['low'],var['high'],'k--',opts,ax)
    
        for it in np.unique(DBSCAN_dataset[:,-1]):
            if np.shape(DBSCAN_dataset[DBSCAN_dataset[:,-1]==it])[0] >= .1*np.shape(dates)[0] and it != -1:
                if plot_me == 1:
                    ax.scatter(DBSCAN_dataset[DBSCAN_dataset[:,-1]==it][:,0],DBSCAN_dataset[DBSCAN_dataset[:,-1]==it][:,-3],DBSCAN_dataset[DBSCAN_dataset[:,-1]==it][:,-2],s=5)
                    
                count = count+1
                sigData = np.append(sigData,DBSCAN_dataset[DBSCAN_dataset[:,-1]==it][:,:-1],axis=0)
    
            else:
                if plot_me == 1:
                    ax.scatter(DBSCAN_dataset[DBSCAN_dataset[:,-1]==it][:,0],DBSCAN_dataset[DBSCAN_dataset[:,-1]==it][:,-3],DBSCAN_dataset[DBSCAN_dataset[:,-1]==it][:,-2],c='k',s=5,label='_nolegend_')
        if plot_me == 1:

            plt.legend(Legend)
            ax.scatter(dates[:,0],dates[:,-2],dates[:,-1],c='k',s=5,label='_nolegend_')
            
            ax.set_xlabel('Departure Date')
            ax.set_ylabel('Flyby Date')
            ax.set_zlabel('Arrival Date')
            plt.title(alg + ' Migration %i, Dep, Fly, & Arr Scatter Plot' % num_mig)
            plt.tight_layout()
            plt.show()
            
            p.waitforbuttonpress(1)


    return var

def boundingBoxDiag(low,high,color,opts,ax):
    # Bounding Box
    sortVarLow = np.sort(low)
    sortVarHigh = np.sort(high)
    
    x = [sortVarLow[-2],sortVarLow[-2],sortVarLow[-2],sortVarLow[-2],sortVarHigh[-2],sortVarHigh[-2],sortVarHigh[-2],sortVarHigh[-2]]
    y = [sortVarLow[-1],sortVarLow[-1],sortVarLow[-1]+np.sum(opts['tof_margin'][0]),sortVarLow[-1]+np.sum(opts['tof_margin'][0]),sortVarHigh[-1]-np.sum(opts['tof_margin'][0]),sortVarHigh[-1]-np.sum(opts['tof_margin'][0]),sortVarHigh[-1],sortVarHigh[-1]]
    z = [sortVarLow[-1]+sortVarLow[-3],sortVarLow[-1]+sortVarHigh[-3],sortVarLow[-1]+sortVarLow[-3]+np.sum(opts['tof_margin'][0]),sortVarLow[-1]+sortVarHigh[-3]+np.sum(opts['tof_margin'][0]),sortVarHigh[-1]+sortVarHigh[-3]-np.sum(opts['tof_margin'][0])-np.sum(opts['tof_margin'][1]),sortVarHigh[-1]+sortVarHigh[-3]-np.sum(opts['tof_margin'][0]),sortVarHigh[-1] + sortVarHigh[-3] - np.sum(opts['tof_margin'][1]),sortVarHigh[-1]+sortVarHigh[-3]]
    
    ax.plot([x[0],x[4]],[y[0],y[4]],[z[0],z[4]],color,label='_nolegend_')
    ax.plot([x[1],x[5]],[y[1],y[5]],[z[1],z[5]],color,label='_nolegend_')
    ax.plot([x[2],x[6]],[y[2],y[6]],[z[2],z[6]],color,label='_nolegend_')
    ax.plot([x[3],x[7]],[y[3],y[7]],[z[3],z[7]],color,label='_nolegend_')
    
    ax.plot([x[0],x[1]],[y[0],y[1]],[z[0],z[1]],color,label='_nolegend_')
    ax.plot([x[1],x[3]],[y[1],y[3]],[z[1],z[3]],color,label='_nolegend_')
    ax.plot([x[2],x[3]],[y[2],y[3]],[z[2],z[3]],color,label='_nolegend_')
    ax.plot([x[0],x[2]],[y[0],y[2]],[z[0],z[2]],color,label='_nolegend_')

    ax.plot([x[4],x[5]],[y[4],y[5]],[z[4],z[5]],color,label='_nolegend_')
    ax.plot([x[5],x[7]],[y[5],y[7]],[z[5],z[7]],color,label='_nolegend_')
    ax.plot([x[6],x[7]],[y[6],y[7]],[z[6],z[7]],color,label='_nolegend_')
    ax.plot([x[4],x[6]],[y[4],y[6]],[z[4],z[6]],color,label='_nolegend_')
    return

def boundingBoxRect(low,high,color,opts,ax):
    # Bounding Box
    xL = 1*low
    xU = 1*high
    
    ax.plot([xL[0],xL[0]],[xL[1],xL[1]],[xL[2],xU[2]],color,label='_nolegend_') 
    ax.plot([xL[0],xL[0]],[xL[1],xU[1]],[xU[2],xU[2]],color,label='_nolegend_') 
    ax.plot([xL[0],xL[0]],[xU[1],xU[1]],[xL[2],xU[2]],color,label='_nolegend_') 
    ax.plot([xL[0],xL[0]],[xL[1],xU[1]],[xL[2],xL[2]],color,label='_nolegend_') 
    
    ax.plot([xU[0],xU[0]],[xL[1],xL[1]],[xL[2],xU[2]],color,label='_nolegend_') 
    ax.plot([xU[0],xU[0]],[xL[1],xU[1]],[xU[2],xU[2]],color,label='_nolegend_') 
    ax.plot([xU[0],xU[0]],[xU[1],xU[1]],[xL[2],xU[2]],color,label='_nolegend_') 
    ax.plot([xU[0],xU[0]],[xL[1],xU[1]],[xL[2],xL[2]],color,label='_nolegend_') 
    
    ax.plot([xL[0],xU[0]],[xL[1],xL[1]],[xL[2],xL[2]],color,label='_nolegend_') 
    ax.plot([xL[0],xU[0]],[xL[1],xL[1]],[xU[2],xU[2]],color,label='_nolegend_') 
    ax.plot([xL[0],xU[0]],[xU[1],xU[1]],[xU[2],xU[2]],color,label='_nolegend_') 
    ax.plot([xL[0],xU[0]],[xU[1],xU[1]],[xL[2],xL[2]],color,label='_nolegend_')
    return