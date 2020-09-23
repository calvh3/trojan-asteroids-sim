from solver import *

'''
Libration Grid Example for JUPITER:
    -Plot a 2D polar grid showing how far test particles librate in position space.
    -Resolution is set at 50x50 and each point is computed over 30 orbits.
    -The plot is scaled radially around the computed region to show more detail.
'''

#Define the planet Jupiter
Jupiter = Planet(Mp=0.001,R=5.2)

x0 = [Jupiter.Lx,Jupiter.Ly,0,0]

trojan3 = Trojan(Jupiter,x0,width=0.5,res=50,Polar=True)
librationgrid=trojan3.LibrationGrid(30)


#Lambda function radially scales data by a factor of 'scale' about the sun-planet barycentre.
scale = 10
rp=lambda r:((r-Jupiter.Rp)*scale+Jupiter.Rp)
(Tscale,Rscale)=(rp(trojan3.r)*np.cos(trojan3.theta),rp(trojan3.r)*np.sin(trojan3.theta))

#Render the heat map from the data using MatPlotLib, Wander (Libration) in units of pi
plt.axis("equal")
plt.contourf(Tscale,Rscale,librationgrid/np.pi,levels=np.linspace(0,1.9,80),cmap='inferno')
plt.colorbar(ticks=np.linspace(0,2,5)).set_label('Angular Libration /units of $\\pi$ Radians')
plt.contour(Tscale,Rscale,librationgrid, levels=[0.2*np.pi,0.4*np.pi,0.6*np.pi], 
            colors="white", linewidths=0.15) #Add contours 

#Mark the L4 point on the plot
plt.text(Jupiter.Lx*1.19, Jupiter.Ly*1.11, r'L$_4$')
plt.plot(Jupiter.Lx, Jupiter.Ly, 'w+')

#Mark the Sun and planet on the plot
plt.plot(-Jupiter.Rs, 0,'o' ,color='orange', markersize=15)     
plt.text(-Jupiter.Rs*0.92, 0.058*Jupiter.R, 'Sun')
plt.plot(Jupiter.Rp, 0, 'ro', markersize=5)
plt.text(Jupiter.Rp*0.88, 0.029*Jupiter.R, 'Planet')

#Label Axis
plt.xlabel('x / Au')
plt.ylabel('y / Au')

#Plot line of constant Radius at Planet Orbital Radius
Jx = np.linspace(-Jupiter.R-Jupiter.Rs, Jupiter.Rp, 100)
Jy = np.sqrt(abs(Jupiter.R**2-(Jx+Jupiter.Rs)**2))
plt.plot(Jx, Jy, 'r--', linewidth=0.35)

#Plot inner and outer bounding radii to make clear the radial scaling used in the plot
Rinbound = Jupiter.R-0.1
Routbound = Jupiter.R+0.1

Bx = np.linspace(-rp(Rinbound), rp(Rinbound), 100)
By = np.sqrt(rp(Rinbound)**2-(Bx)**2)
plt.plot(Bx, By, 'k--',linewidth=0.75)
Bx = np.linspace(-rp(Routbound), rp(Routbound), 100)
By = np.sqrt(rp(Routbound)**2-(Bx)**2)
plt.plot(Bx, By, 'k--',linewidth=0.75)

#Plot (double headed) arrow to show the distance between the two bounding lines.
inner=-rp(Rinbound)
outer=-rp(Routbound)
mid=0.5*(inner+outer)
plt.arrow(mid,0,inner-mid,0,head_length=0.2,head_width=0.1,
          length_includes_head=True,color='black',linewidth=0.75)
plt.arrow(mid,0,outer-mid,0,head_length=0.2,head_width=0.1,
          length_includes_head=True,color='black',linewidth=0.75)
plt.text(-Jupiter.R, -0.096*Jupiter.R, '$\\Delta R=%.3f\\mathrm{Au}$'%(Routbound-Rinbound))
plt.text(-Jupiter.R, -0.2*Jupiter.R, '(x%.0f Scaling)'%(scale))

#Plot an arrow to show the distance between the Sun and planet.
plt.arrow(0,0,Jupiter.R,0,head_length=0.2,head_width=0.1,
          length_includes_head=True,color='black',linewidth=0.75)
plt.text(Jupiter.R/2, -0.096*Jupiter.R, 'R={}Au'.format(Jupiter.R))
plt.show()
