# FORM: [] = plotOrbits(bod,const,opts,var,plot_vars)
#
# |-----------------------------------------------------------------------
# | NOTES:
# |     -Function designed to plot MGALT transfers
# |
# |     -This funciton plots the departure planet, flyby planet(s), and 
# |     the target planet along with the spacecraft's trajectories between 
# |     all of them. The funciton also plots the thruster pointing angle 
# |     for all of the transfers
# |
# |-----------------------------------------------------------------------
# |
# | INPUTS:
# |     -bod                (1,1)       [struct]        [unitless]
# |         A struct containing information pertaining to the planetary
# |         bodies. Contains list of bodies, launch windows and ToF, and 
# |         planetary R/V/JD vectors. This struct has dynamic fields and 
# |         will adapt to contain only the necesary information
# |     -const              (1,1)       [struct]        [unitless]
# |         A struct containing constants used in the calcs. Contains
# |         values for AU, TU, Sun (rad/mu/rp) and (rad/mu/rp/SOI/per) 
# |         for any bodies used in the optsimization scheme. This is a 
# |         dynamic struct and will adapt to contain only the necesary 
# |         information
# |     -opts                (1,1)       [struct]        [unitless]
# |         A struct containing constants user optsions. Contains the save 
# |         folder, ToF values, and more structs containing informaiton 
# |         for the island model, cost parameters, weighting parameters, 
# |         and all of the islands used in the optsimization process
# |     -var                (1,1)       [struct]        [unitless]
# |         A struct containing the variable limits
# |     -plot_vars          (1,1)       [struct]     	[unitless]
# |         An object containing a lot of information about the 
# |         optsimization parameters including: transfers(t and y ode 
# |         outputs), thrust values, thruster pointing angles, transfer 
# |         starting position, planet start/end locations for each 
# |         transfer, JD of each transfer, and tspans of each transfer
# |
# |-----------------------------------------------------------------------
# |
# | OUTPUTS:
# |
# |-----------------------------------------------------------------------
# |
# | MISC:
# |
# |-----------------------------------------------------------------------

import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.integrate import odeint
from Other_Functions.orbit3D import orbit3D
import pylab as p

 
def plotOrbits(bod,const,opts,var,plot_vars):
    #%%
    ## User Choices
    
    # Size of Sun and Planets
    planet_size = 20
    sun_size = 40
    
    # General line sizes and shapes
    line_width_orbit = 2
    line_width_thruster_angle = 2
    line_width_terminal_point = 5
    line_shape_terminal_point = 'kd'
    
    # Direct method points and arrows
    point_size = 10
    arrow_length = 0.3
    arrow_head_size = 3
    arrow_width = .002
    scale_arrow_length = .1
    
    # Font
    font_size = 24
    font_style = 'Times'
    font_terminal = 'Terminal Points'
    
    # Colors
    sun_color = (0.95,0.95,0)
    arrow_color = (0,0.5,0)
    
    ## constants
    
    # Get constants
    mew_sun = const['Sun_mu'][0]   # km^3/s^2
    AU = const['AU'][0]              # km/AU
    SOI_divisions = np.arange(0,2*math.pi,math.pi/180,dtype=float)
    orbit_margin = 1.01
    
    # Plot colors
    plot_col = plt.rcParams['axes.prop_cycle'].by_key()['color'][:7]
    
    ## Plot Orbits
    #%%
    if opts['solver'] == 'LT_DIR_FSM_2D':
        '''
        #***FIGURE 1: ORBITAL TRANSFERS***
        figure(1)
        
        # Grid Basics
        plot(0,0,'color',sun_color,'marker','.','markersize',sun_size)  # Plot the Sun
        hold on
        grid on
        Legend1{1} = 'Sun'
            
        # Plot the Departure body Position
        body_data = plot_vars.planetary_conditions(:,1)
        body_X = body_data(1)/AU	# DU
        body_Y = body_data(2)/AU 	# DU
        plot(body_X,body_Y,'Color',plot_col(1,:),'marker','.','markersize',planet_size)    # body Position

        # Plot the Departure body's SOI
        body_SOI = const.(strcat(bod.bodies{1},"_SOI"))/AU
        body_SOI_X = (body_SOI/2)*cos(SOI_divisions) + body_X
        body_SOI_Y = (body_SOI/2)*sin(SOI_divisions) + body_Y
        body_SOI_line(1,1) = plot(body_SOI_X,body_SOI_Y,'k--')
        set(get(get(body_SOI_line(1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

        # Plot the Departure body's Orbit
        body_period = const.(strcat(bod.bodies{1},"_per"))*orbit_margin
        body_prop_period = linspace(0,body_period*86400,500)
        [~,body_orbit] = ode45(@orbit3D,body_prop_period,plot_vars.planetary_conditions(:,1),opts.ode,mew_sun)
        body_orbit_line(1,1) = plot(body_orbit(:,1)/AU,body_orbit(:,2)/AU,'k--')
        set(get(get(body_orbit_line(1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')      # Exclude line from legend
        thrust_plot_dist = max([norm(body_orbit(1,:)),norm(body_orbit(2,:))])                          	# Used for plotting, the length of the thrust pointing line
        clear body_period body_prop_period body_orbit

        # Plot the Spacecraft Transfer Arc
        sc_data = plot_vars.transfers{1,1}
        sc_pos_rad = sc_data(:,2)        	# DU
        sc_pos_ang = sc_data(:,3)*pi/180	# rad
        [sc_X,sc_Y] = pol2cart(sc_pos_ang,sc_pos_rad)      # Convert to Cartesian
        plot(sc_X,sc_Y,'Color',plot_col(2,:),'linewidth',line_width_orbit)      # Spacecraft Trajectory

        # Calculate the Thrust Vectors
        thrust_switch = plot_vars.thrust_switch{1,1}                   # binary
        thrust_angle = plot_vars.thrust_phi{1,1}                       # deg
        thrust_R = plot_vars.seg_start{1,1}(:,1)                       # DU
        thrust_arc = plot_vars.seg_start{1,1}(:,2)                     # deg
        [thrust_X,thrust_Y] = pol2cart(thrust_arc*pi/180,thrust_R)     # points for the thrust vectors
        thrust_mag = arrow_length*ones(length(thrust_switch),1)       	# magnitude of thrust vectors for visual
        thrust_ang_global = thrust_arc + 90 - thrust_angle'         	# global angle for thrust vectors
        [thrust_X_global,thrust_Y_global] = pol2cart(thrust_ang_global*pi/180,thrust_mag)	# points for the thrust vectors 

        # Plot the Thrust Vectors
        quiver_count = 1
        for i2 = 1:length(thrust_X)

            if thrust_switch(i2)
                thrust_point_line(1,i2) = plot(thrust_X(i2),thrust_Y(i2),...
                    'Color',arrow_color,'marker','.','markersize',point_size)
                thrust_quiver_line(1,quiver_count) = quiver(thrust_X(i2),thrust_Y(i2),thrust_X_global(i2),thrust_Y_global(i2),...
                    'color',[0,0.5,0],'MaxHeadSize',(arrow_head_size*AU),'LineWidth',arrow_width)
                set(get(get(thrust_quiver_line(1,quiver_count),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')	# Exclude line from legend
                quiver_count = quiver_count+1
            else
                thrust_point_line(1,i2) = plot(thrust_X(i2),thrust_Y(i2),...
                    'Color',arrow_color,'marker','o','markersize',point_size/4)
            end

            set(get(get(thrust_point_line(1,i2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')              # Exclude line from legend

        end

        # Labels
        Legend1{2} = bod.bodies{1}
        transfer_string1 = strcat("Transfer ",bod.bodies{1}," to ",bod.bodies{end})
        Legend1{3} = (transfer_string1)
        
        # Plot the Target body
        body_data = plot_vars.planetary_conditions(:,end)
        body_X = body_data(1)/AU	# DU
        body_Y = body_data(2)/AU 	# DU
        plot(body_X,body_Y,'Color',plot_col(2,:),'marker','.','markersize',planet_size)    # body Position
        Legend1{4} = bod.bodies{end}
        
        # Plot the Target body's SOI
        body_SOI = const.(strcat(bod.bodies{end},"_SOI"))/AU
        body_SOI_X = (body_SOI/2)*cos(SOI_divisions) + body_X
        body_SOI_Y = (body_SOI/2)*sin(SOI_divisions) + body_Y
        body_SOI_line(end,1) = plot(body_SOI_X,body_SOI_Y,'k--')
        set(get(get(body_SOI_line(end,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
        
        # Plot the Target body's Orbit
        body_period = const.(strcat(bod.bodies{2},"_per"))*orbit_margin
        body_prop_period = linspace(0,body_period*86400,500)
        [~,body_orbit] = ode45(@orbit3D,body_prop_period,plot_vars.planetary_conditions(:,end),opts.ode,mew_sun)
        body_orbit_line(end,1) = plot(body_orbit(:,1)/AU,body_orbit(:,2)/AU,'k--')
        set(get(get(body_orbit_line(end,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
        
        # Final Graph Elements
        axis square
        hold off
        xlabel('X Position (AU)')
        ylabel('Y Positon (AU)')
        legend(string(Legend1),'location','eastoutside')
        title('Low-Thrust Orbit Transfer with Direct FSM')
        set(gca,'FontSize',font_size)
        set(gca,'FontName',font_style)
        
        
        
        #***FIGURE 2: THRUSTER ANGLE***
        figure(2)
        
        # Grid Basics
        hold on
        grid on
            
        # Plot Thruster Angle
        for i2 = 1:size(plot_vars.tspan{1,1},1)

            # If full or dashed
            if plot_vars.thrust_switch{1,1}(i2)
                thrust_segment_line(1,i2) = plot(plot_vars.tspan{1,1}(i2,:),...
                    [plot_vars.thrust_phi{1,1}(i2),plot_vars.thrust_phi{1,1}(i2)],...
                    'Color',plot_col(2,:),'linewidth',line_width_thruster_angle)
                set(get(get(thrust_segment_line(1,i2),'Annotation'),...
                    'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
            else
                thrust_segment_line(1,i2) = plot(plot_vars.tspan{1,1}(i2,:),...
                    [plot_vars.thrust_phi{1,1}(i2),plot_vars.thrust_phi{1,1}(i2)],...
                    '--','Color',plot_col(2,:),'linewidth',line_width_thruster_angle)
                set(get(get(thrust_segment_line(1,i2),'Annotation'),...
                    'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
            end

        end

        # Plot connecting lines
        for i3 = 2:size(plot_vars.tspan{1,1},1)

            # If full or dashed
            if plot_vars.thrust_switch{1,1}(i3)
                thrust_vertical_line(1,i3) = plot([plot_vars.tspan{1,1}(i3,1),plot_vars.tspan{1,1}(i3,1)],...
                    [plot_vars.thrust_phi{1,1}(i3-1),plot_vars.thrust_phi{1,1}(i3)],...
                    'Color',plot_col(2,:),'linewidth',line_width_thruster_angle)
                set(get(get(thrust_vertical_line(1,i3),'Annotation'),...
                    'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
            else
                thrust_vertical_line(1,i3) = plot([plot_vars.tspan{1,1}(i3,1),plot_vars.tspan{1,1}(i3,1)],...
                    [plot_vars.thrust_phi{1,1}(i3-1),plot_vars.thrust_phi{1,1}(i3)],...
                    '--','Color',plot_col(2,:),'linewidth',line_width_thruster_angle)
                set(get(get(thrust_vertical_line(1,i3),'Annotation'),...
                    'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
            end

        end

        # Plot Thruster Angle
        plot([0,0.00001],[0,00001],'Color',plot_col(2,:),'linewidth',line_width_thruster_angle)
        transfer_string2 = strcat("Transfer ",bod.bodies{1}," to ",bod.bodies{end})
        Legend2{1} = (transfer_string2)

        # Final Graph Elements
        hold off
        xlabel('Time (TU)')
        ylabel('\phi (deg)')
        legend(string(Legend2),'location','eastoutside')
        title('Thruster Pointing Angle \phi for Direct FSM')
        set(gca,'FontSize',font_size)
        set(gca,'FontName',font_style)
        ylim([0 360])
        yticks([0,45,90,135,180,225,270,315,360])
        '''
    elif opts['solver'] == 'LT_IN_FSM_2D':
        '''
        #***FIGURE 1: ORBITAL TRANSFERS***
        figure(1)
        
        # Grid Basics
        plot(0,0,'color',sun_color,'marker','.','markersize',sun_size)  # Plot the Sun
        hold on
        grid on
        Legend1{1} = 'Sun'
            
        # Plot the Departure body Position
        body_data = plot_vars.planetary_conditions(:,1)
        body_X = body_data(1)/AU	# DU
        body_Y = body_data(2)/AU 	# DU
        plot(body_X,body_Y,'Color',plot_col(1,:),'marker','.','markersize',planet_size)    # body Position

        # Plot the Departure body's SOI
        body_SOI = const.(strcat(bod.bodies{1},"_SOI"))/AU
        body_SOI_X = (body_SOI/2)*cos(SOI_divisions) + body_X
        body_SOI_Y = (body_SOI/2)*sin(SOI_divisions) + body_Y
        body_SOI_line(1,1) = plot(body_SOI_X,body_SOI_Y,'k--')
        set(get(get(body_SOI_line(1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

        # Plot the Departure body's Orbit
        body_period = const.(strcat(bod.bodies{1},"_per"))*orbit_margin
        body_prop_period = linspace(0,body_period*86400,500)
        [~,body_orbit] = ode45(@orbit3D,body_prop_period,plot_vars.planetary_conditions(:,1),opts.ode,mew_sun)
        body_orbit_line(1,1) = plot(body_orbit(:,1)/AU,body_orbit(:,2)/AU,'k--')
        set(get(get(body_orbit_line(1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
        clear body_period body_prop_period body_orbit

        # Plot the Spacecraft Transfer Arc
        sc_data = plot_vars.transfers{1,1}
        sc_pos_rad = sc_data(:,5)            # DU
        sc_pos_ang = sc_data(:,6)*pi/180   # rad
        [sc_X,sc_Y] = pol2cart(sc_pos_ang,sc_pos_rad)    # Convert to Cartesian
        plot(sc_X,sc_Y,'Color',plot_col(2,:),'linewidth',line_width_orbit)	# Spacecraft Trajectory

        # Labels
        Legend1{2} = bod.bodies{1}
        transfer_string1 = strcat("Transfer ",bod.bodies{1}," to ",bod.bodies{end})
        Legend1{3} = (transfer_string1)
        
        # Plot the Target body
        body_data = plot_vars.planetary_conditions(:,end)
        body_X = body_data(1)/AU	# DU
        body_Y = body_data(2)/AU 	# DU
        plot(body_X,body_Y,'Color',plot_col(2,:),'marker','.','markersize',planet_size)    # body Position
        Legend1{4} = bod.bodies{end}
        
        # Plot the Target body's SOI
        body_SOI = const.(strcat(bod.bodies{end},"_SOI"))/AU
        body_SOI_X = (body_SOI/2)*cos(SOI_divisions) + body_X
        body_SOI_Y = (body_SOI/2)*sin(SOI_divisions) + body_Y
        body_SOI_line(end,1) = plot(body_SOI_X,body_SOI_Y,'k--')
        set(get(get(body_SOI_line(end,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
        
        # Plot the Target body's Orbit
        body_period = const.(strcat(bod.bodies{2},"_per"))*orbit_margin
        body_prop_period = linspace(0,body_period*86400,500)
        [~,body_orbit] = ode45(@orbit3D,body_prop_period,plot_vars.planetary_conditions(:,end),opts.ode,mew_sun)
        body_orbit_line(end,1) = plot(body_orbit(:,1)/AU,body_orbit(:,2)/AU,'k--')
        set(get(get(body_orbit_line(end,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
        
        # Final Graph Elements
        axis square
        hold off
        xlabel('X Position (AU)')
        ylabel('Y Positon (AU)')
        legend(string(Legend1),'location','eastoutside')
        title('Low-Thrust Orbit Transfer with Indirect FSM')
        set(gca,'FontSize',font_size)
        set(gca,'FontName',font_style)
        
        
        
        #***FIGURE 2: THRUSTER ANGLE***
        figure(2)
        
        # Grid Basics
        hold on
        grid on
        time = 0
        phi = 0

        # Thruster Data
        thruster_data = plot_vars.transfers{1,1}
        time = thruster_data(:,1)+time(end)
        lambda_1 = thruster_data(:,2)
        lambda_2 = thruster_data(:,3)

        # Run through all thruster angles
        for i2 = 1:length(lambda_1)
            phi(i2) = atan2d(-lambda_1(i2)/sqrt(lambda_1(i2)^2 + lambda_2(i2)^2),...
                            -lambda_2(i2)/sqrt(lambda_1(i2)^2 + lambda_2(i2)^2))
            if phi(i2) < 0
                phi(i2) = 360 + phi(i2)
            end
        end

        # Plot Thruster Angle
        plot(time,phi,'Color',plot_col(2,:),'linewidth',line_width_thruster_angle)
        transfer_string2 = strcat("Transfer ",bod.bodies{1}," to ",bod.bodies{end})
        Legend2{1} = (transfer_string2)
        
        # Final Graph Elements
        hold off
        xlabel('Time (TU)')
        ylabel('\phi (deg)')
        legend(string(Legend2),'location','eastoutside')
        title('Thruster Pointing Angle \phi for Indirect FSM')
        set(gca,'FontSize',font_size)
        set(gca,'FontName',font_style)
        ylim([0 360])
        yticks([0,45,90,135,180,225,270,315,360])
        '''
    elif opts['solver'] == 'MGALT_DIR_FBSM_2D':
        
        #***FIGURE 1: ORBITAL TRANSFERS***
        fig = plt.figure()
        ax = fig.add_subplot() 
        
        # Grid Basics
        plt.grid()
        ax.plot(0,0,marker='.',color=sun_color,markersize=sun_size)  # Plot the Sun
        count = 0
        
        Legend1 = [[] for i in range(2*var['transfers']+3)]
        Legend1[count] = 'Sun'        
        
        # Run through all the transfers
        body_SOI_line = [[] for i in range(var['transfers'])]
        body_orbit_line = [[] for i in range(var['transfers'])]
        sc_fs_terminal_point = [[] for i in range(var['transfers'])]
        sc_bs_terminal_point = [[] for i in range(var['transfers'])]
        sc_transfer_line = [[] for i in range(var['transfers'])]
        thrust_point_line_fs = [[[] for i in range(np.shape(plot_vars['seg_start_fs'])[0])] for ii in range(var['transfers'])]
        thrust_quiver_line_fs = [[[] for i in range(np.shape(plot_vars['seg_start_fs'])[0])] for ii in range(var['transfers'])]
        thrust_point_line_bs = [[[] for i in range(np.shape(plot_vars['seg_start_fs'])[0])] for ii in range(var['transfers'])]
        thrust_quiver_line_bs = [[[] for i in range(np.shape(plot_vars['seg_start_fs'])[0])] for ii in range(var['transfers'])]

        for i1 in range(var['transfers']):
            
            count = count+1
            
            # Plot the Departure body Position
            body_data = plot_vars['planetary_conditions'][:,i1]
            body_X = body_data[0]/AU	# DU
            body_Y = body_data[1]/AU 	# DU
            ax.plot(body_X,body_Y,color=plot_col[i1],marker='.',markersize=planet_size)   # body Position
            
            # Plot the Departure body's SOI
            body_SOI = const[bod['bodies'][i1]+"_SOI"][0]/AU
            body_SOI_X = [(body_SOI/2)*math.cos(SOI_divisions[i]) + body_X for i in range(len(SOI_divisions))]
            body_SOI_Y = [(body_SOI/2)*math.sin(SOI_divisions[i]) + body_Y for i in range(len(SOI_divisions))]
            body_SOI_line[i1] = ax.plot(body_SOI_X,body_SOI_Y,'k--',label='_nolegend_')
            # set(get(get(body_SOI_line(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
            
            # Plot the Departure body's Orbit            
            body_period = const[bod['bodies'][i1]+"_per"][0]*orbit_margin
            body_prop_period = np.linspace(0,body_period*86400,500)
            body_orbit = odeint(orbit3D,plot_vars['planetary_conditions'][:,i1],body_prop_period,args=(mew_sun,))[:,:2]
            body_orbit_line[i1] = ax.plot(body_orbit[:,0]/AU,body_orbit[:,1]/AU,'k--', label='_nolegend_')
            thrust_plot_dist = max([np.linalg.norm(body_orbit[0,:]),np.linalg.norm(body_orbit[1,:])])         # Used for plotting, the length of the thrust pointing line
            # set(get(get(body_orbit_line(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
            # clear body_period body_prop_period body_orbit
            
            # Plot the Spacecraft Forward Transfer Arc
            sc_data_fs = plot_vars['transfers_fs'][:,i1*int(np.shape(plot_vars['transfers'])[1]/var['transfers']):(i1+1)*int(np.shape(plot_vars['transfers'])[1]/var['transfers'])]
            sc_pos_rad_fs = sc_data_fs[:,1]            # DU
            sc_pos_ang_fs = sc_data_fs[:,2]*math.pi/180   # rad
            sc_X_fs = []
            sc_Y_fs = []
            for iter1 in range(len(sc_pos_rad_fs)):           
                sc_X_fs = np.append(sc_X_fs,sc_pos_rad_fs[iter1]*math.cos(sc_pos_ang_fs[iter1]))    # Convert to Cartesian
                sc_Y_fs = np.append(sc_Y_fs,sc_pos_rad_fs[iter1]*math.sin(sc_pos_ang_fs[iter1]))    # Convert to Cartesian
            ax.plot(sc_X_fs,sc_Y_fs,color=plot_col[i1+1],linewidth=line_width_orbit)	# Spacecraft Trajectory
            sc_fs_terminal_point[i1] = ax.plot(sc_X_fs[-1],sc_Y_fs[-1],line_shape_terminal_point,linewidth=line_width_terminal_point, label='_nolegend_')
            # set(get(get(sc_fs_terminal_point(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
            
            # Plot the Spacecraft Backward Transfer Arc
            sc_data_bs = plot_vars['transfers_bs'][:,i1*int(np.shape(plot_vars['transfers'])[1]/var['transfers']):(i1+1)*int(np.shape(plot_vars['transfers'])[1]/var['transfers'])]
            sc_pos_rad_bs = sc_data_bs[:,1]            # DU
            sc_pos_ang_bs = sc_data_bs[:,2]*math.pi/180   # rad
            sc_X_bs = []
            sc_Y_bs = []
            for iter1 in range(len(sc_pos_rad_bs)):           
                sc_X_bs = np.append(sc_X_bs,sc_pos_rad_bs[iter1]*math.cos(sc_pos_ang_bs[iter1]))    # Convert to Cartesian
                sc_Y_bs = np.append(sc_Y_bs,sc_pos_rad_bs[iter1]*math.sin(sc_pos_ang_bs[iter1]))    # Convert to Cartesian
            sc_transfer_line[i1] = ax.plot(sc_X_bs,sc_Y_bs,color=plot_col[i1+1],linewidth=line_width_orbit, label='_nolegend_')	# Spacecraft Trajectory
            sc_bs_terminal_point[i1] = ax.plot(sc_X_bs[-1],sc_Y_bs[-1],line_shape_terminal_point,linewidth=line_width_terminal_point, label='_nolegend_')
            # set(get(get(sc_transfer_line(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
            # set(get(get(sc_bs_terminal_point(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
            
            # Calculate the Thrust Vectors
            thrust_switch_fs = plot_vars['thrust_switch_fs'][0,int(i1*np.shape(plot_vars['thrust_switch_fs'])[1]/var['transfers']):int((i1+1)*np.shape(plot_vars['thrust_switch_fs'])[1]/var['transfers'])]    # binary
            thrust_switch_bs = plot_vars['thrust_switch_bs'][0,int(i1*np.shape(plot_vars['thrust_switch_fs'])[1]/var['transfers']):int((i1+1)*np.shape(plot_vars['thrust_switch_fs'])[1]/var['transfers'])]    # binary
            thrust_angle_fs = plot_vars['thrust_phi_fs'][0,int(i1*np.shape(plot_vars['thrust_switch_fs'])[1]/var['transfers']):int((i1+1)*np.shape(plot_vars['thrust_switch_fs'])[1]/var['transfers'])]       # deg
            thrust_angle_bs = plot_vars['thrust_phi_bs'][0,int(i1*np.shape(plot_vars['thrust_switch_fs'])[1]/var['transfers']):int((i1+1)*np.shape(plot_vars['thrust_switch_fs'])[1]/var['transfers'])]    # binary
            thrust_R_fs = plot_vars['seg_start_fs'][:,int(i1*5)]       	# DU
            thrust_R_bs = plot_vars['seg_start_bs'][:,int(i1*5)]
            thrust_arc_fs = plot_vars['seg_start_fs'][:,int(i1*5)+1]     # deg
            thrust_arc_bs = plot_vars['seg_start_bs'][:,int(i1*5)+1]
            thrust_X_fs = []
            thrust_Y_fs = []
            for iter1 in range(len(thrust_R_fs)):
                thrust_X_fs = np.append(thrust_X_fs,math.cos(math.radians(thrust_arc_fs[iter1]))*thrust_R_fs[iter1])
                thrust_Y_fs = np.append(thrust_Y_fs,math.sin(math.radians(thrust_arc_fs[iter1]))*thrust_R_fs[iter1])
            thrust_X_bs = []
            thrust_Y_bs = []
            for iter1 in range(len(thrust_R_bs)):
                thrust_X_bs = np.append(thrust_X_bs,math.cos(math.radians(thrust_arc_bs[iter1]))*thrust_R_bs[iter1])
                thrust_Y_bs = np.append(thrust_Y_bs,math.sin(math.radians(thrust_arc_bs[iter1]))*thrust_R_bs[iter1])
            thrust_mag_fs = arrow_length*np.ones(len(thrust_switch_fs))
            thrust_mag_bs = arrow_length*np.ones(len(thrust_switch_bs))
            thrust_ang_global_fs = thrust_arc_fs + 90 - thrust_angle_fs       	# global angle for thrust vectors
            thrust_ang_global_bs = thrust_arc_bs + 90 - thrust_angle_bs
            thrust_X_global_fs = []
            thrust_Y_global_fs = []
            thrust_X_global_bs = []
            thrust_Y_global_bs = []
            for iter1 in range(len(thrust_X_fs)):           
                thrust_X_global_fs = np.append(thrust_X_global_fs,thrust_mag_fs[iter1]*math.cos(math.radians(thrust_ang_global_fs[iter1])))    # Convert to Cartesian
                thrust_Y_global_fs = np.append(thrust_Y_global_fs,thrust_mag_fs[iter1]*math.sin(math.radians(thrust_ang_global_fs[iter1])))    # Convert to Cartesian
                thrust_X_global_bs = np.append(thrust_X_global_bs,thrust_mag_bs[iter1]*math.cos(math.radians(thrust_ang_global_bs[iter1])))    # Convert to Cartesian
                thrust_Y_global_bs = np.append(thrust_Y_global_bs,thrust_mag_bs[iter1]*math.sin(math.radians(thrust_ang_global_bs[iter1])))    # Convert to Cartesian
            
            # Plot the Forward Shooting Thrust Vectors
            quiver_count_fs = 0
            for i2 in range(len(thrust_X_fs)):
                
                if thrust_switch_fs[i2]:
                    thrust_point_line_fs[i1][i2] = ax.plot(thrust_X_fs[i2],thrust_Y_fs[i2],color=arrow_color,marker='.',markersize=point_size, label='_nolegend_')
                    thrust_quiver_line_fs[i1][quiver_count_fs] = ax.quiver(thrust_X_fs[i2],thrust_Y_fs[i2],thrust_X_global_fs[i2]*scale_arrow_length,thrust_Y_global_fs[i2]*scale_arrow_length,color=(0,0.5,0),scale_units='inches',scale=.1,headwidth=(arrow_head_size),width=arrow_width, label='_nolegend_')
                    # set(get(get(thrust_quiver_line_fs(i1,quiver_count_fs),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')                              # Exclude line from legend
                    quiver_count_fs = quiver_count_fs+1
                else:
                    thrust_point_line_fs[i1][i2] = ax.plot(thrust_X_fs[i2],thrust_Y_fs[i2],color=arrow_color,marker='o',markersize=point_size/2)
                
                # set(get(get(thrust_point_line_fs(i1,i2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')             # Exclude line from legend
            
            # Plot the Backward Shooting Thrust Vectors
            quiver_count_bs = 0
            for i3 in range(len(thrust_X_bs)):
                
                if thrust_switch_bs[i3]:
                    thrust_point_line_bs[i1][i3] = ax.plot(thrust_X_bs[i3],thrust_Y_bs[i3],color=arrow_color,marker='.',markersize=point_size, label='_nolegend_')
                    thrust_quiver_line_bs[i1][quiver_count_bs] = ax.quiver(thrust_X_bs[i3],thrust_Y_bs[i3],thrust_X_global_bs[i3]*scale_arrow_length,thrust_Y_global_bs[i3]*scale_arrow_length,color=(0,0.5,0),scale_units='inches',scale=.1,headwidth=(arrow_head_size),width=arrow_width, label='_nolegend_')
                    # set(get(get(thrust_quiver_line_fs(i1,quiver_count_fs),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')                              # Exclude line from legend
                    quiver_count_bs = quiver_count_bs+1
                else:
                    thrust_point_line_bs[i1][i3] = ax.plot(thrust_X_bs[i3],thrust_Y_bs[i3],color=arrow_color,marker='o',markersize=point_size/2)

            
            # Labels
            Legend1[count] = bod['bodies'][i1]
            transfer_string1 = "Transfer "+bod['bodies'][i1]+" to "+bod['bodies'][i1+1]
            count = count+1
            Legend1[count] = transfer_string1
        
        # Plot the Target body
        body_data = plot_vars['planetary_conditions'][:,-1]
        body_X = body_data[0]/AU	# DU
        body_Y = body_data[1]/AU 	# DU
        ax.plot(body_X,body_Y,color=plot_col[i1+1],marker='.',markersize=planet_size)    # body Position
        Legend1[count+1] = bod['bodies'][-1]
        
        # Plot the Target body's SOI
        body_SOI = const[bod['bodies'][-1]+"_SOI"][0]/AU
        body_SOI_X = [(body_SOI/2)*math.cos(SOI_divisions[i]) + body_X for i in range(len(SOI_divisions))]
        body_SOI_Y = [(body_SOI/2)*math.sin(SOI_divisions[i]) + body_Y for i in range(len(SOI_divisions))]
        body_SOI_line[i1] = ax.plot(body_SOI_X,body_SOI_Y,'k--',label='_nolegend_')
        
        # Plot the Target body's Orbit
        body_period = const[bod['bodies'][-1]+"_per"][0]*orbit_margin
        body_prop_period = np.linspace(0,body_period*86400,500)
        body_orbit = odeint(orbit3D,plot_vars['planetary_conditions'][:,-1],body_prop_period,args=(mew_sun,))[:,:2]
        body_orbit_line[i1] = ax.plot(body_orbit[:,0]/AU,body_orbit[:,1]/AU,'k--', label='_nolegend_')
        thrust_plot_dist = max([np.linalg.norm(body_orbit[0,:]),np.linalg.norm(body_orbit[1,:])])         # Used for plotting, the length of the thrust pointing line
        
        # Have 1 match point displayed
        sc_bs_terminal_point[i1][0] = ax.plot(sc_X_bs[-1],sc_Y_bs[-1],line_shape_terminal_point,linewidth=line_width_terminal_point)
        Legend1[count+2] = font_terminal
        
        # Final Graph Elements
        ax.set_aspect(1)
        plt.rc('axes', titlesize=font_size)
        plt.rc('axes', labelsize=font_size)
        plt.xlabel('X Position (AU)')
        plt.ylabel('Y Positon (AU)')
        plt.legend(Legend1,loc='center left', bbox_to_anchor=(1, 0.5),prop={'size': font_size})
        plt.title('Low-Thrust Orbit Transfer with Direct FBSM')
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized() 
        plt.show()
        
        #***FIGURE 2: THRUSTER ANGLE***
        fig = plt.figure()
        ax = fig.add_subplot()
        
        # Grid Basics
        plt.grid()
        
        time_fs = np.zeros(np.shape(plot_vars['tspan_fs']))
        time_bs = np.zeros(np.shape(plot_vars['tspan_bs']))
        # # time_fs = 0
        # # time_bs = 0
        count = 0
        
        # Run through all the transfers
        thrust_segment_line_fs = [[[] for i in range(np.shape(plot_vars['tspan_fs'][i1])[0])] for ii in range(var['transfers'])]
        thrust_vertical_line_fs = [[[] for i in range(np.shape(plot_vars['tspan_fs'][i1])[0])] for ii in range(var['transfers'])]
        thrust_segment_line_bs = [[[] for i in range(np.shape(plot_vars['tspan_bs'][i1])[0])] for ii in range(var['transfers'])]
        thrust_vertical_line_bs = [[[] for i in range(np.shape(plot_vars['tspan_bs'][i1])[0])] for ii in range(var['transfers'])]
        Legend2 = [[] for i in range(2*var['transfers']+3)]

        for i1 in range(var['transfers']):
            
            # Thruster Data
            if var['transfers'] == 1:
                time_fs = plot_vars['tspan_fs'] + time_bs[-1,-1]
                time_bs = plot_vars['tspan_bs']
                time_bs = np.flipud(time_bs) + time_fs[-1,-1]
            else:
                time_fs = plot_vars['tspan_fs'][i1] + time_bs[-1,0]
                time_bs = plot_vars['tspan_bs'][i1]
                time_bs = np.flipud(time_bs) + time_fs[-1,-1]
            
            # Plot Forward Thruster Angle
            for i2 in range(int(opts['thrust']['Nseg']/2)):
                
                # If full or dashed
                if plot_vars['thrust_switch_fs'][0,i1*int(opts['thrust']['Nseg']/2)+i2]:
                    thrust_segment_line_fs[i1][i2] = ax.plot([time_fs[i2,0],time_fs[i2,-1]],[plot_vars['thrust_phi_fs'][0,i1*int(opts['thrust']['Nseg']/2)+i2],plot_vars['thrust_phi_fs'][0,i1*int(opts['thrust']['Nseg']/2)+i2]],color=plot_col[i1+1],linewidth=line_width_thruster_angle, label='_nolegend_')
                    # set(get(get(thrust_segment_line_fs(i1,i2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
                else:
                    thrust_segment_line_fs[i1][i2] = ax.plot([time_fs[i2,0],time_fs[i2,-1]],[plot_vars['thrust_phi_fs'][0,i1*int(opts['thrust']['Nseg']/2)+i2],plot_vars['thrust_phi_fs'][0,i1*int(opts['thrust']['Nseg']/2)+i2]],'--',color=plot_col[i1+1],linewidth=line_width_thruster_angle, label='_nolegend_')
                    # set(get(get(thrust_segment_line_fs(i1,i2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
            
            # Plot Forward Connecting Lines
            for i3 in range(1,int(opts['thrust']['Nseg']/2)):
                
                # If full or dashed
                if plot_vars['thrust_switch_fs'][0,i1*int(opts['thrust']['Nseg']/2)+i3]:
                    thrust_vertical_line_fs[i1][i3] = ax.plot([time_fs[i3,0],time_fs[i3,0]],[plot_vars['thrust_phi_fs'][0,i1*int(opts['thrust']['Nseg']/2)+i3-1],plot_vars['thrust_phi_fs'][0,i1*int(opts['thrust']['Nseg']/2)+i3]],color=plot_col[i1+1],linewidth=line_width_thruster_angle, label='_nolegend_')

                else:
                    thrust_vertical_line_fs[i1][i3] = ax.plot([time_fs[i3,0],time_fs[i3,0]],[plot_vars['thrust_phi_fs'][0,i1*int(opts['thrust']['Nseg']/2)+i3-1],plot_vars['thrust_phi_fs'][0,i1*int(opts['thrust']['Nseg']/2)+i3]],'--',color=plot_col[i1+1],linewidth=line_width_thruster_angle, label='_nolegend_')
            
            # Plot Backward Thruster Angle
            for i4 in range(1,int(opts['thrust']['Nseg']/2)):
                
                # If full or dashed
                if plot_vars['thrust_switch_bs'][0,i1*int(opts['thrust']['Nseg']/2)+i4]:
                    thrust_segment_line_bs[i1][i4] = ax.plot([time_bs[-i4-1,0],time_bs[-i4-1,-1]],[plot_vars['thrust_phi_bs'][0,i1*int(opts['thrust']['Nseg']/2)+i4],plot_vars['thrust_phi_bs'][0,i1*int(opts['thrust']['Nseg']/2)+i4]],color=plot_col[i1+1],linewidth=line_width_thruster_angle, label='_nolegend_')

                else:
                    thrust_segment_line_bs[i1][i4] = ax.plot([time_bs[-i4-1,0],time_bs[-i4-1,-1]],[plot_vars['thrust_phi_bs'][0,i1*int(opts['thrust']['Nseg']/2)+i4],plot_vars['thrust_phi_bs'][0,i1*int(opts['thrust']['Nseg']/2)+i4]],'--',color=plot_col[i1+1],linewidth=line_width_thruster_angle, label='_nolegend_')

            # Plot Backward Connecting Lines
            for i5 in range(1,int(opts['thrust']['Nseg']/2)):
                
                # If full or dashed
                if plot_vars['thrust_switch_bs'][0,i1*int(opts['thrust']['Nseg']/2)+i5]:
                    thrust_vertical_line_bs[i1][i5] = ax.plot([time_bs[-i5-1,0],time_bs[-i5-1,0]],[plot_vars['thrust_phi_bs'][0,i1*int(opts['thrust']['Nseg']/2)+i5-1],plot_vars['thrust_phi_bs'][0,i1*int(opts['thrust']['Nseg']/2)+i5]],color=plot_col[i1+1],linewidth=line_width_thruster_angle, label='_nolegend_')

                else:
                    thrust_vertical_line_fs[i1][i5] = ax.plot([time_bs[-i5-1,0],time_bs[-i5-1,0]],[plot_vars['thrust_phi_bs'][0,i1*int(opts['thrust']['Nseg']/2)+i5-1],plot_vars['thrust_phi_bs'][0,i1*int(opts['thrust']['Nseg']/2)+i5]],'--',color=plot_col[i1+1],linewidth=line_width_thruster_angle, label='_nolegend_')
            
            # Plot Match Points
            sc_fs_terminal_point[i1] = ax.plot(time_fs[i2,-1],plot_vars['thrust_phi_fs'][0,i1*int(opts['thrust']['Nseg']/2)+i2],line_shape_terminal_point,linewidth=line_width_terminal_point, label='_nolegend_')
            sc_bs_terminal_point[i1] = ax.plot(time_bs[0,-1],plot_vars['thrust_phi_bs'][0,i1*int(opts['thrust']['Nseg']/2) + int(opts['thrust']['Nseg']/2)-1],line_shape_terminal_point,linewidth=line_width_terminal_point, label='_nolegend_')
            
            # Plot Thruster Angle
            ax.plot([0,0.00001],[0,0.00001],color=plot_col[i1+1],linewidth=line_width_thruster_angle)
            transfer_string2 = "Transfer " + bod['bodies'][i1]+" to "+bod['bodies'][i1+1]
            Legend2[count] = transfer_string2
            count = count+1

        # Final Graph Elements
        ax.plot(time_bs[0,-1],plot_vars['thrust_phi_bs'][0,i1*int(opts['thrust']['Nseg']/2) + int(opts['thrust']['Nseg']/2)-1],line_shape_terminal_point,linewidth=line_width_terminal_point, label='_nolegend_')
        Legend2[count] = font_terminal

        plt.xlabel('Time (TU)')
        plt.ylabel('$\phi$ (deg)')
        plt.legend(Legend2,loc='center left', bbox_to_anchor=(1, 0.5),prop={'size': font_size/2})
        plt.title('Thruster Pointing Angle $\phi$ for Direct FBSM')
        plt.ylim((-45,360))
        plt.yticks((-45,0,45,90,135,180,225,275,315,360))
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized() 
        plt.tight_layout()
        plt.show()
        #%%
    elif opts['solver'] == 'MGALT_IN_FBSM_2D':
        
        #***FIGURE 1: ORBITAL TRANSFERS***
        fig = plt.figure()
        ax = fig.add_subplot()
        
        # Grid Basics
        plt.grid()
        ax.plot(0,0,marker='.',color=sun_color,markersize=sun_size)  # Plot the Sun
        count = 0
        
        Legend1 = [[] for i in range(2*var['transfers']+3)]
        Legend1[count] = 'Sun'
        
        # Run through all the transfers
        body_SOI_line = [[] for i in range(var['transfers'])]
        body_orbit_line = [[] for i in range(var['transfers'])]
        sc_fs_terminal_point = [[] for i in range(var['transfers'])]
        sc_bs_terminal_point = [[] for i in range(var['transfers'])]
        sc_transfer_line = [[] for i in range(var['transfers'])]
        
        for i1 in range(var['transfers']):
            
            count = count+1
            
            # Plot the Departure body Position
            body_data = plot_vars['planetary_conditions'][:,i1]
            body_X = body_data[0]/AU	# DU
            body_Y = body_data[1]/AU 	# DU
            ax.plot(body_X,body_Y,color=plot_col[i1],marker='.',markersize=planet_size)   # body Position
            
            # Plot the Departure body's SOI
            body_SOI = const[bod['bodies'][i1]+"_SOI"][0]/AU
            body_SOI_X = [(body_SOI/2)*math.cos(SOI_divisions[i]) + body_X for i in range(len(SOI_divisions))]
            body_SOI_Y = [(body_SOI/2)*math.sin(SOI_divisions[i]) + body_Y for i in range(len(SOI_divisions))]
            body_SOI_line[i1] = ax.plot(body_SOI_X,body_SOI_Y,'k--',label='_nolegend_')
            # set(get(get(body_SOI_line(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
            
            # Plot the Departure body's Orbit
            body_period = const[bod['bodies'][i1]+"_per"][0]*orbit_margin
            body_prop_period = np.linspace(0,body_period*86400,500)
            body_orbit = odeint(orbit3D,plot_vars['planetary_conditions'][:,i1],body_prop_period,args=(mew_sun,))[:,:2]
            body_orbit_line[i1] = ax.plot(body_orbit[:,0]/AU,body_orbit[:,1]/AU,'k--', label='_nolegend_')
            # set(get(get(body_orbit_line(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
            # clear body_period body_prop_period body_orbit
            
            # Plot the Spacecraft Forward Transfer Arc
            sc_data_fs = plot_vars['transfers_fs'][:,i1*int(np.shape(plot_vars['transfers'])[1]/var['transfers']):(i1+1)*int(np.shape(plot_vars['transfers'])[1]/var['transfers'])]
            sc_pos_rad_fs = sc_data_fs[:,4]            # DU
            sc_pos_ang_fs = sc_data_fs[:,5]*math.pi/180   # rad
            sc_X_fs = []
            sc_Y_fs = []
            for iter1 in range(len(sc_pos_rad_fs)):           
                sc_X_fs = np.append(sc_X_fs,sc_pos_rad_fs[iter1]*math.cos(sc_pos_ang_fs[iter1]))    # Convert to Cartesian
                sc_Y_fs = np.append(sc_Y_fs,sc_pos_rad_fs[iter1]*math.sin(sc_pos_ang_fs[iter1]))    # Convert to Cartesian
            ax.plot(sc_X_fs,sc_Y_fs,color=plot_col[i1+1],linewidth=line_width_orbit)	# Spacecraft Trajectory
            sc_fs_terminal_point[i1] = ax.plot(sc_X_fs[-1],sc_Y_fs[-1],line_shape_terminal_point,linewidth=line_width_terminal_point, label='_nolegend_')
            # set(get(get(sc_fs_terminal_point(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
            
            # Plot the Spacecraft Backwards Transfer Arc
            sc_data_bs = plot_vars['transfers_bs'][:,i1*int(np.shape(plot_vars['transfers'])[1]/var['transfers']):(i1+1)*int(np.shape(plot_vars['transfers'])[1]/var['transfers'])]
            sc_pos_rad_bs = sc_data_bs[:,4]            # DU
            sc_pos_ang_bs = sc_data_bs[:,5]*math.pi/180   # rad
            sc_X_bs = []
            sc_Y_bs = []
            for iter1 in range(len(sc_pos_rad_bs)):           
                sc_X_bs = np.append(sc_X_bs,sc_pos_rad_bs[iter1]*math.cos(sc_pos_ang_bs[iter1]))    # Convert to Cartesian
                sc_Y_bs = np.append(sc_Y_bs,sc_pos_rad_bs[iter1]*math.sin(sc_pos_ang_bs[iter1]))    # Convert to Cartesian
            sc_transfer_line[i1] = ax.plot(sc_X_bs,sc_Y_bs,color=plot_col[i1+1],linewidth=line_width_orbit, label='_nolegend_')	# Spacecraft Trajectory
            sc_bs_terminal_point[i1] = ax.plot(sc_X_bs[-1],sc_Y_bs[-1],line_shape_terminal_point,linewidth=line_width_terminal_point, label='_nolegend_')
            # set(get(get(sc_transfer_line(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
            # set(get(get(sc_bs_terminal_point(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
            
            # Labels
            Legend1[count] = bod['bodies'][i1]
            transfer_string1 = "Transfer "+bod['bodies'][i1]+" to "+bod['bodies'][i1+1]
            count = count+1
            Legend1[count] = transfer_string1
            
        
        # Plot the Target body
        body_data = plot_vars['planetary_conditions'][:,-1]
        body_X = body_data[0]/AU	# DU
        body_Y = body_data[1]/AU 	# DU
        ax.plot(body_X,body_Y,color=plot_col[i1+1],marker='.',markersize=planet_size)    # body Position
        Legend1[count+1] = bod['bodies'][-1]
        
        # Plot the Target body's SOI
        body_SOI = const[bod['bodies'][-1]+"_SOI"][0]/AU
        body_SOI_X = [(body_SOI/2)*math.cos(SOI_divisions[i]) + body_X for i in range(len(SOI_divisions))]
        body_SOI_Y = [(body_SOI/2)*math.sin(SOI_divisions[i]) + body_Y for i in range(len(SOI_divisions))]
        body_SOI_line[-1] = ax.plot(body_SOI_X,body_SOI_Y,'k--', label='_nolegend_')
        # set(get(get(body_SOI_line(end,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
        
        # Plot the Target body's Orbit
        body_period = const[bod['bodies'][-1]+"_per"][0]*orbit_margin
        body_prop_period = np.linspace(0,body_period*86400,500)
        body_orbit = odeint(orbit3D,plot_vars['planetary_conditions'][:,-1],body_prop_period,args=(mew_sun,))[:,:2]
        body_orbit_line[-1] = ax.plot(body_orbit[:,0]/AU,body_orbit[:,1]/AU,'k--', label='_nolegend_')
        # set(get(get(body_orbit_line(end,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
        
        # Have 1 match point displayed
        sc_bs_terminal_point[i1] = ax.plot(sc_X_bs[-1],sc_Y_bs[-1],line_shape_terminal_point,linewidth=line_width_terminal_point)
        Legend1[count+2] = font_terminal
        
        # Final Graph Elements
        ax.set_aspect(1)
        plt.rc('axes', titlesize=font_size)
        plt.rc('axes', labelsize=font_size)
        plt.xlabel('X Position (AU)')
        plt.ylabel('Y Positon (AU)')
        plt.legend(Legend1,loc='center left', bbox_to_anchor=(1, 0.5),prop={'size': font_size})
        plt.title('Low-Thrust Orbit Transfer with Indirect FBSM')
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized() 
        plt.show()
        
        #***FIGURE 2: THRUSTER ANGLE***
        fig = plt.figure()
        ax = fig.add_subplot()
        
        # Grid Basics
        plt.grid()
        
        time_fs = np.zeros((len(plot_vars['transfers_fs']),))
        time_bs = np.zeros((len(plot_vars['transfers_bs']),))
        count = 0
        
        # Run through all the transfers
        sc_bs_line = [[] for i in range(var['transfers'])]
        Legend2 = [[] for i in range(2*var['transfers']+3)]
        for i1 in range(var['transfers']):
            
            # Thruster Data
            thruster_fs = plot_vars['transfers_fs'][:,i1*int(np.shape(plot_vars['transfers'])[1]/var['transfers']):(i1+1)*int(np.shape(plot_vars['transfers'])[1]/var['transfers'])]
            time_fs = thruster_fs[:,0]+time_bs[-1]
            lambda_1_fs = thruster_fs[:,1]
            lambda_2_fs = thruster_fs[:,2]
            
            thruster_bs = plot_vars['transfers_bs'][:,i1*int(np.shape(plot_vars['transfers'])[1]/var['transfers']):(i1+1)*int(np.shape(plot_vars['transfers'])[1]/var['transfers'])]
            thruster_bs = np.flipud(thruster_bs)
            time_bs = thruster_bs[:,0]+time_fs[-1]
            lambda_1_bs = thruster_bs[:,1]
            lambda_2_bs = thruster_bs[:,2]
            
            # Preallocate
            phi_fs = np.zeros((len(time_fs),))
            phi_bs = np.zeros((len(time_bs),))
            
            # Run through all thruster angles
            for i2 in range(len(lambda_1_fs)):
                phi_fs[i2] = math.degrees(math.atan2(-lambda_1_fs[i2]/math.sqrt(lambda_1_fs[i2]**2 + lambda_2_fs[i2]**2),-lambda_2_fs[i2]/math.sqrt(lambda_1_fs[i2]**2 + lambda_2_fs[i2]**2)))
                if phi_fs[i2] < 0:
                    phi_fs[i2] = 360 + phi_fs[i2]
            
            for i3 in range(len(lambda_1_bs)):
                phi_bs[i3] = math.degrees(math.atan2(-lambda_1_bs[i3]/math.sqrt(lambda_1_bs[i3]**2 + lambda_2_bs[i3]**2),-lambda_2_bs[i3]/math.sqrt(lambda_1_bs[i3]**2 + lambda_2_bs[i3]**2)))
                if phi_bs[i3] < 0:
                    phi_bs[i3] = 360 + phi_bs[i3]
            
            # Plot Thruster Angle
            ax.plot(time_fs,phi_fs,color=plot_col[i1+1],linewidth=line_width_thruster_angle)
            sc_bs_line[i1] = ax.plot(time_bs,phi_bs,color=plot_col[i1+1],linewidth=line_width_thruster_angle,label='_nolegend_')
            # set(get(get(sc_bs_line(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
            
            sc_fs_terminal_point[i1] = ax.plot(time_fs[-1],phi_fs[-1],line_shape_terminal_point,linewidth=line_width_terminal_point,label='_nolegend_')
            # set(get(get(sc_fs_terminal_point(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
            sc_bs_terminal_point[i1] = ax.plot(time_bs[0],phi_bs[0],line_shape_terminal_point,linewidth=line_width_terminal_point,label='_nolegend_')
            # set(get(get(sc_bs_terminal_point(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')     #Exclude line from legend
            
            transfer_string2 = "Transfer "+bod['bodies'][i1]+" to "+bod['bodies'][i1+1]
            Legend2[count] = transfer_string2
            count = count+1
        
        # Final Graph Elements
        sc_bs_terminal_point[i1] = ax.plot(time_bs[0],phi_bs[0],line_shape_terminal_point,linewidth=line_width_terminal_point)
        Legend2[count] = font_terminal

        plt.xlabel('Time (TU)')
        plt.ylabel('$\phi$ (deg)')
        plt.legend(Legend2,loc='center left', bbox_to_anchor=(1, 0.5),prop={'size': font_size/2})
        plt.title('Thruster Pointing Angle $\phi$ for Indirect FBSM')
        plt.ylim((0,360))
        plt.yticks((0,45,90,135,180,225,270,315,360))
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized() 
        plt.show()
        
        p.waitforbuttonpress(1)
        
    else:
        
        raise Exception("Incorrect Solver Selection")
    
    return
