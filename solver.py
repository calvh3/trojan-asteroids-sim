"""
This Programme solves the restricted three-body problem for a particle placed 
in the Jupiter/Mars systems. The chosen reference frame has origin at the Sun-Planet 
barycentre and rotates at the planets angular velocity. 
"""

import numpy as np
import scipy.integrate
import cmath
import matplotlib.pyplot as plt
#The Planet to be simulated is defined as a class 
#Define a Planet by Planet Mass (Mp) and Sun-Planet Radius R.
class Planet:
    def __init__(self, Mp, R):
        self.Mp = Mp                                #Planet Mass (in Solar Masses)
        self.R = R                                  #Sun-Planet Radius (in Au)
        
        self.setup(Mp,R)
        
    #setup defines/calculates all of the required constants associated with Planet.
    def setup(self, Mp, R):
        
        #Define Uiversal System Constants
        self.G = 4*np.pi**2                         #Gravitational Constant
        self.Ms = 1.0                               #Relatvie mass of Sun
        
        #Calculate Planet Specific Constants
        self.Rs = self.R*self.Mp/(self.Ms+self.Mp)  #Sun - Barycentre (Origin) Distance
        self.Rp = self.R*self.Ms/(self.Ms+self.Mp)  #Planet - Barycentre (Origin) Distance
        self.T = self.Ms*self.R**1.5                #Planet orbital time period (From Kepler III)
        self.Omega = 2*np.pi/self.T                 #Planet Angular Velocity
        
        #L4 Lagrangian Point Coordiantes
        self.Lx = self.Rp-self.R/2                  #Lx (x coordinate)
        self.Ly = np.sqrt(3)*self.R/2               #Ly (y coordiante)

class Trojan:
    def __init__(self,Planet,X0,width=0,res=1,Polar=False):
        #Define a Trojan from its start position at t=0, and planet class in trojan system.
        self.Planet = Planet
        self.res = res
        self.Polar = Polar
        if self.res == 0:
            self.res = 1
            print("cannot have zero resolution (set: res=1)")
        self.res2 = self.res*self.res
        if self.res==1:
            self.width = 0
        else: 
            self.width = width
        if Polar == False:
            self.X0 = X0
            delta = self.width/self.res              #Step size 
            (self.x0, self.y0)=np.meshgrid(self.X0[0]-self.width/2+delta*np.arange(0,self.res,1),self.X0[1]-self.width/2+delta*np.arange(0,self.res,1))
            
            (Vx,Vy)=(np.zeros(self.res2),np.zeros(self.res2))
        
            self.X00 = np.concatenate((self.x0.reshape(self.res2),self.y0.reshape(self.res2),Vx.reshape(self.res2),Vy.reshape(self.res2)))

        if Polar == True:
            Rin = np.sqrt(self.Planet.Lx**2+self.Planet.Ly**2)-width/2
            Rout = Rin+width
            (self.theta,self.r)=np.meshgrid(np.linspace(0.1,np.pi,self.res),np.linspace(Rin,Rout,self.res))
            (self.x0,self.y0)=(self.r*np.cos(self.theta),self.r*np.sin(self.theta))
            self.X00 = np.concatenate((self.x0.reshape(self.res2),self.y0.reshape(self.res2),np.zeros(self.res2),np.zeros(self.res2)))


    def OribtSolve(self,orbits,n=1,Precision=100,solve_ivp=False):
        #Solve Trojan trajetory over given number of oribts.
        self.orbits = orbits
        self.t = np.linspace(0,self.orbits*self.Planet.T,self.orbits*Precision)
        def OrbitSolve (X0, t, n):
        #Solve the ODE to find the trajectory of the orbit,
        #Starting position X0, over a set of time points, t
            def Derivative (t,X):
            #Equations of Motion
                Rx, Ry, Vx, Vy = (X[0:n],X[n:2*n],X[2*n:3*n],X[3*n:])
                
                Rst = np.clip(np.sqrt((self.Planet.Rs+Rx)**2+Ry**2),0.05,None) #Sun-Trojan Distance
                Rpt = np.clip(np.sqrt((Rx-self.Planet.Rp)**2+Ry**2),0.05,None) #self.Planet-Trojan Distance
                #Clipping is to stop errors in the integrator for orbits that 'hit' each body.
                
                dXdt= (Vx, Vy, 
                       -self.Planet.G*((Rx+self.Planet.Rs)*self.Planet.Ms/(Rst**3)+
                                       (Rx-self.Planet.Rp)*self.Planet.Mp/(Rpt**3))+
                                        Rx*self.Planet.Omega**2+2*self.Planet.Omega*Vy,
                       -self.Planet.G*(Ry*self.Planet.Ms/(Rst**3)+Ry*self.Planet.Mp/(Rpt**3))+
                                        Ry*self.Planet.Omega**2-2*self.Planet.Omega*Vx)
                return np.concatenate(dXdt)
            
            '''
            ######################################################################################################
            Choose Integrator:
            odeint, uses  LSODA from the FORTRAN library odepack, switches stiff/nonstiff solver.
            solve_ivp, is set to use RK45, use with 2D heatmaps.
            '''
            if(solve_ivp==True):
                X = scipy.integrate.solve_ivp(Derivative,(self.t[0],self.t[-1]),X0,
                                              t_eval=self.t,method='RK45',vectorized=True,max_step=0.5).y
            else:
                X = scipy.integrate.odeint(Derivative,X0,t,tfirst=True).T
            
            
            '''
            ######################################################################################################
            '''                
            return X
        self.X= OrbitSolve(self.X00, self.t,n)
        return self.X
        
    def OrbitQuickPlot(self):
        '''
        Plot a single Orbit
        '''
        try:
            if isinstance(self.X,np.ndarray):
                #Plot the orbit in 2D showing the Sun and Planet
                plt.plot(self.X[0],self.X[1],'k',linewidth=0.75)
                #Mark L4 point
                plt.text(self.Planet.Lx+0.1, self.Planet.Ly+0.1, 'L$_4$',color='blue')
                plt.plot(self.Planet.Lx, self.Planet.Ly, 'b+')
                #Mark the Sun and self.Planet
                plt.text(-self.Planet.Rs*1.08, 0.3, 'Sun')
                plt.plot(-self.Planet.Rs, 0,'o' ,color='orange', markersize=15)
                plt.text(self.Planet.Rp*0.88, 0.15, 'Planet')
                plt.plot(self.Planet.Rp, 0, 'ro', markersize=5)
                #Label Axis
                plt.xlabel('x / Au')
                plt.ylabel('y / Au')
                plt.tight_layout()
                plt.axis("equal")
                return plt.show()
        except AttributeError:
            print("Compute an Orbit first using .OribtSolve")
            
    def Wander(self,x,y):
        #Computes the wander distance of orbit, ie max absolute distance in space from the L4 Lagrangain Point (Lx,Ly). 
        try:
            if isinstance(self.X,np.ndarray):
                self.wander = np.max(np.sqrt((x-self.Planet.Lx)**2+(y-self.Planet.Ly)**2),axis=0)
                return self.wander
        except AttributeError:
            print("Compute an Orbit first using .OribtSolve")

    def Libration(self,x,y):
        #Computes the libration angle of orbit, ie the max-min angle the asteroid makes about the origin.
        try:
            if isinstance(self.X,np.ndarray):
                theta = np.angle(x+y*1j); theta[theta<0] += 2*np.pi
                self.libration = np.max(theta,axis=0)-np.min(theta,axis=0)
                return self.libration
        except AttributeError:
            print("Compute an Orbit first using .OribtSolve")

    def WanderGrid(self,orbits):
        X = self.OribtSolve(orbits,n=self.res2,solve_ivp=True).T
        self.wandergrid = self.Wander(X[:,0:self.res2],X[:,self.res2:2*self.res2]).reshape((self.res,self.res))
        return self.wandergrid
    
    def LibrationGrid(self,orbits):
        X = self.OribtSolve(orbits,n=self.res2,solve_ivp=True).T
        self.librationgrid = self.Libration(X[:,0:self.res2],X[:,self.res2:2*self.res2]).reshape((self.res,self.res))                
        return self.librationgrid
    

