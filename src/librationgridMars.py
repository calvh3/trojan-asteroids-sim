from solver import *

'''
Libration Grid Example for MARS:
    -Plot a 2D polar grid showing how far test particles librate in position space.
    -Resolution is set at 50x50 and each point is computed over 120 orbits.
    -The plot is scaled radially around the computed region to show more detail.
'''

#Define the planet Mars
Mars = Planet(Mp=3.2e-7,R=1.5)

x0 = [Mars.Lx,Mars.Ly,0,0]

trojan3 = Trojan(Mars,x0,width=0.006,res=50,Polar=True)
librationgrid=trojan3.LibrationGrid(120)


#Lambda function radially scales data by a factor of 'scale' about the sun-planet barycentre.
scale = 150
rp=lambda r:((r-Mars.Rp)*scale+Mars.Rp)
(Tscale,Rscale)=(rp(trojan3.r)*np.cos(trojan3.theta),rp(trojan3.r)*np.sin(trojan3.theta))

#Render the heat map from the data using MatPlotLib, Wander (Libration) in units of pi
plt.axis("equal")
plt.contourf(Tscale,Rscale,librationgrid/np.pi,levels=np.linspace(0,1.9,80),cmap='inferno')
plt.colorbar(ticks=np.linspace(0,2,5)).set_label('Angular Libration /units of $\\pi$ Radians')
plt.contour(Tscale,Rscale,librationgrid, levels=[0.2*np.pi,0.4*np.pi,0.6*np.pi], 
            colors="white", linewidths=0.15) #Add contours 

#Mark the L4 point on the plot
plt.text(Mars.Lx*1.19, Mars.Ly*1.11, r'L$_4$')
plt.plot(Mars.Lx, Mars.Ly, 'w+')

#Mark the Sun and planet on the plot
plt.plot(-Mars.Rs, 0,'o' ,color='orange', markersize=15)     
plt.text(-Mars.Rs*0.92, 0.058*Mars.R, 'Sun')
plt.plot(Mars.Rp, 0, 'ro', markersize=5)
plt.text(Mars.Rp*0.88, 0.029*Mars.R, 'Planet')

#Label Axis
plt.xlabel('x / Au')
plt.ylabel('y / Au')

#Plot line of constant Radius at Planet Orbital Radius
Jx = np.linspace(-Mars.R-Mars.Rs, Mars.Rp, 100)
Jy = np.sqrt(abs(Mars.R**2-(Jx+Mars.Rs)**2))
plt.plot(Jx, Jy, 'r--', linewidth=0.35)

#Plot inner and outer bounding radii to make clear the radial scaling used in the plot
Rinbound = Mars.R-0.003
Routbound = Mars.R+0.003

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
plt.text(-Mars.R, -0.096*Mars.R, '$\\Delta R=%.3f\\mathrm{Au}$'%(Routbound-Rinbound))
plt.text(-Mars.R, -0.2*Mars.R, '(x%.0f Scaling)'%(scale))

#Plot an arrow to show the distance between the Sun and planet.
plt.arrow(0,0,Mars.R,0,head_length=0.2,head_width=0.1,
          length_includes_head=True,color='black',linewidth=0.75)
plt.text(Mars.R/2, -0.096*Mars.R, 'R={}Au'.format(Mars.R))

plt.show()
