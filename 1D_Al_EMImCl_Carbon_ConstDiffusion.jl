"""
Created on Mon Dec 31 9pm 2018

Porous electrode theory for Al:EMImCl carbon powder electrode intercalating AlCl4-

This file executes simulations of the solid phase intercalation of a typical pcl
at each locaiton in the electrode, and also a simulation of the diffusion of ions
in the electrolyte.



Simulation of the solid-phase intercalation eqns --- see Damon's notes in Zhang2000
the vector's row/matrix format is setup as shwon below.
the row direction is into the particles
the column direction is along the electrode
                      ppppppppppppppppppppppppppppppppppppppppppppppppp
                      ccccccccccccccccccccccccccccccccccccccccccccccccc
                      lllllllllllllllllllllllllllllllllllllllllllllllll
                      1234567891011121314151617181920212223242526272829
  current collector | ---electrolyte------ y=0 -----electrolyte----- y=1|
                                       pcl center
f
For 1-D linear diffusion equation (dm/dt = -D d2m/dx_pcl2) the stability criteria is dt < dx_pcl^2/D

@author: damon
"""

function main()
    ####### EQUATIONS FOR FREE ENERGY OF THE ELECTRON IN THE ELECTROCHEMICAL REACTION ############
    #Li-ion intercalation plateaus from Zhang2000
    #The equation from Dolye1996 and Zhang2000 is typo-fixed error-fixed:
    #Ueq = 4.19829 + 0.0565661*tanh(-14.5546*y[-1] + 8.60942) - 0.0275479 * ((0.998432 - y[-1])**-0.492465 - 1.90111) - 0.157123*exp(-0.04738*y[-1]**8) + 0.810239 * exp(-40*y[-1] + 5.355)
    #One intercalation plateau (see Excel file)
    #Ueq=25*exp(-((1-yy)+0.02)/0.04)+2-35*exp(((1-yy)-1.03)/0.04)-(1-yy)*0.05
    #Three intercalation plateaus (see Excel file)
    #Ueq=25*exp(-((1-yy)+0.04)/0.02)-0.025*(1+erf(((1-yy)-0.66)/(0.03*sqrt(2))))-0.025*(1+erf(((1-yy)-0.34)/(0.03*sqrt(2))))+2-25*exp(((1-yy)-1.03)/0.02)-(1-yy)*0.03
    #Sloping Three Peaks
    #Ueq=25*exp(-((1-yy)+0.04)/0.02)-0.015*(1+erf(((1-yy)-0.66)/(0.03*sqrt(2))))-0.015*(1+erf(((1-yy)-0.34)/(0.03*sqrt(2))))+2.3-25*exp(((1-yy)-1.03)/0.02)-(1-yy)*0.5
    #Sloping
    #Ueq=25*exp(-((1-yy)+0.04)/0.02)+2.3-25*exp(((1-yy)-1.03)/0.02)-(1-yy)*0.5
    #Experimental curves of Ueq vs y[-1]
    #Ueq=interp(y[-1],imported_y,imported_Ueq)
    ######## INSERT ONE OF THESE EQNS INTO THE NEXT LINE ############
    Ueq(yy) = 25*exp(-((1-yy)+0.04)/0.02)+2.3-25*exp(((1-yy)-1.03)/0.02)-(1-yy)*0.5
    T =    298.0    # Kelvin
    F =    96500.0  # Faradays Constant
    R =    8.3      #(J/mole/Celcius)
    beta = 0.5      #Butler-Volmer transfer coefficient
    cl =   1.0      #(mol / L)  Li+ concentration in liquid phase  (constant!)
    r0 =   15E-4    #(cm)  radius of spherical particle
    c0 =   2.28     #(mol / L) initial concentration of intercalated-Li throughout pcl
    ct =   24       #(mol / L) total concentration of intercalation sites (occupied or unoccpied)
    Utop=  2.3      #(mV /s)  CV scan rate
    Ubottom=1.5     #(mV /s)  CV scan rate
    z=     -1.0     #charge of the intercalating ion
    kb =   0.05     ###(cm^5/2 s^-1 mol^1/2)
    D =    1.0E-8   ###(cm^2 /s)
    dtau=  0.0001   #dtau=dt*D/r0/r0 ##
    tau_nodes=Array(0:4000000)*dtau####
    saved_time_spacing = 30000 ####
    vreal= 0.00001 ###(V /s) scan rate
    v=     vreal*(r0*r0/D)

    #Create arrays to hold data
    #dx_real = 1E-4      # in units of cm.  dx is 100 nm
    dx_pcl=       0.05 #  non-dimensional x for the pcl simulation
    pcl_x_nodes=  Array(0:20)*dx_pcl  # 0 to 0.1 mm
    dx_electrode=       0.05 #  non-dimensional x for the pcl simulation
    electrode_x_nodes = Array(0:20)*dx_electrode  # 0 to 0.1 mm
    #dt=       0.001
    dtau=     dtau      #It was defined above
    tau_nodes=tau_nodes #It was defined above
    print("dt is ", dtau*r0*r0/D, " seconds\n")
    print("D dt/dx_pcl^2 is ",dtau/dx_pcl/dx_pcl, " and must be less than 0.5\n")
    y=        ones(size(pcl_x_nodes)[1],size(electrode_x_nodes)[1])*c0/ct  # array of y , y=c/ct
    Ustart=Ueq(y[end,1])# the initial applied potential
    dydx_pcl=     zeros(size(pcl_x_nodes)[1],size(electrode_x_nodes)[1])  # array of y , y=c/ct
    d2ydx_pcl2=   zeros(size(pcl_x_nodes)[1],size(electrode_x_nodes)[1])  # array of y , y=c/ct
    dydtau=     zeros(size(pcl_x_nodes)[1],size(electrode_x_nodes)[1])

    #STABILITY dt < dx_pcl^2/D  , in my case here dt = 0.001  and  dx_pcl^2=0.00001 and D=1
    j=1
    cs=y[end,:]*ct
    eta=zeros(size(cs)[1])
    Ucollector=Ustart
    eta = Ucollector.-Ueq.(y[end,:]) #eta is overpotential
    dydx_pcl[end,:]=r0*z*kb*(ct.-cs)/cl/ct/D.*(exp.(-beta*R*T/F*eta) - exp.((1-beta)*R*T/F*eta) )
    dydx_pcl[1,:].=0
    d2ydx_pcl2[1,:]=(2*y[2,:] - 2*y[1,:]) / dx_pcl / dx_pcl
    d2ydx_pcl2[end,:]=(dydx_pcl[end,:]-dydx_pcl[end-1,:])/dx_pcl
    dydtau = d2ydx_pcl2 #+ 2/pcl_x_nodes*dydx_pcl
    #print(j, Ucollector,eta,cs,dydx_pcl[-1],d2ydx_pcl2[-1])

    #Create the arrays that will save the chronological results to disk, and viewed later for insights
    saved_time_spacing = saved_time_spacing  #It was defined above
    tau_array_saved = cat([1.,2.,4.],6:saved_time_spacing:size(tau_nodes)[1],dims=1)
    y_saved=zeros(size(tau_array_saved)[1],size(pcl_x_nodes)[1],size(electrode_x_nodes)[1])
    i_saved=zeros((size(tau_array_saved)[1]),size(electrode_x_nodes)[1])
    Ucollector_saved=zeros((size(tau_array_saved)[1]))
    y_saved[1,:,:]=y
    i_saved[1,:]=-F*z*D*dydx_pcl[end,:]*ct/r0
    Ucollector_saved[1]=Ucollector
    k=1

    for j in 2:size(tau_nodes)[1]   #loop over time
        Ucollector=Ucollector+v*dtau
        if Ucollector >= Utop && sign(v)==1
            v=-v
            Ucollector=Ucollector+v*dtau*2
        end
        if Ucollector <= Ubottom && sign(v)==-1
            v=-v
            Ucollector=Ucollector+v*dtau*2
        end
        y=y.+dydtau*dtau                #update the concentration field inside the pcl
        cs=y[end,:]*ct                 #calculate the surface concentraction of intercalated ions
        eta = Ucollector.-Ueq.(y[end,:])
        for i in 2:size(pcl_x_nodes)[1]-1
            dydx_pcl[i,:] = (y[i+1,:] - y[i-1,:]) / 2 / dx_pcl
        end
        #I used this following equation for deriving the dydx_pcl[-1] equation seen below: Ueq=4.0 + R*T/F*log((ct/cs-1)/cl)
        dydx_pcl[end,:]=r0*z*kb*(ct.-cs)/cl/ct/D.*(exp.(-beta*R*T/F*eta) - exp.((1-beta)*R*T/F*eta) )
        #print(Ucollector,dydx_pcl[-1])
        dydx_pcl[1,:].=0
        for i in 2:size(pcl_x_nodes)[1]-1
            d2ydx_pcl2[i,:] = (y[i-1,:] - 2*y[i,:] + y[i+1,:]) / dx_pcl / dx_pcl
        end
        d2ydx_pcl2[end,:]=(dydx_pcl[end,:]-dydx_pcl[end-1,:])/dx_pcl
        d2ydx_pcl2[1,:]=(2*y[2,:] - 2*y[1,:]) / dx_pcl / dx_pcl
        dydtau = d2ydx_pcl2 #+ 2/pcl_x_nodes*dydx_pcl

        if any(j.==tau_array_saved)
            k=k+1
            println(k)
            #print(Ucollector,eta,y[1],y[end-2],y[end])
            #println(j, Ucollector,eta,cs,dydx_pcl[end],d2ydx_pcl2[end])
            y_saved[k,:,:]=y
            i_saved[k,:]=-F*z*D*dydx_pcl[end,:]*ct/r0
            Ucollector_saved[k]=Ucollector
        end
    end
    print("elapsed time is ",r0*r0/D*tau_nodes[end], " seconds\n")
    plot(Ucollector_saved,i_saved[:,2])
    #println(Ucollector_saved)
    #println(i_saved[:,1])
end

using Plots
main()
