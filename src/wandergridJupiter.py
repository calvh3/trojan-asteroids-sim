from solver import *

'''
Distance Wander Grid Example for JUPITER:
    -Plot a 2D grid showing how far test particles wander in position space away from the Lagrange Point.
    -Example uses a small 50x50 grid, each point is computed over 30 orbits
'''

#Define the planet Jupiter
Jupiter = Planet(Mp=0.001,R=5.2)

x0 = [Jupiter.Lx,Jupiter.Ly,0,0]
trojan2 = Trojan(Jupiter,x0,width=0.2,res=50)
wandergrid = trojan2.WanderGrid(orbits=30)

R=Jupiter.R 
plt.axis("equal")
plt.contourf(trojan2.x0,trojan2.y0,wandergrid/R,levels=np.linspace(0,2,75),cmap='inferno')
plt.colorbar(ticks=np.linspace(0,2,5)).set_label('Wander from L4 /units of R (%.1fAu)'%(R))
plt.xlabel('x/Au')
plt.ylabel('y/Au')
plt.contour(trojan2.x0,trojan2.y0,wandergrid,[0.2*R,0.4*R,0.6*R,0.8*R,R], 
            colors="white", linewidths=0.25) #Add contours to plot

#Mark the L4 point on the plot
plt.text(Jupiter.Lx*1.19, Jupiter.Ly*1.11, 'L$_4$')
plt.plot(Jupiter.Lx, Jupiter.Ly, 'w+')
plt.show()

'''
Note: the method is vectorised which leads to fast runtimes for small grids but is extremely 
memory intensive for larger grids. To produce larger images, it is best to 'stitch together'
lots of the smaller images.
'''
