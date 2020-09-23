from solver import *

'''
Basic Orbit examples:
-Test particle moved 0.5% along x/y from Lagrange 4 point.
'''

#Define the Planets Jupiter/Mars 
Jupiter = Planet(Mp=0.001,R=5.2)
Mars = Planet(Mp=3.2e-7,R=1.5)

#Jupiter plot
x0=[Jupiter.Lx*1.005,Jupiter.Ly*1.005,0,0]
trojan1 = Trojan(Jupiter,x0)
trojan1.OribtSolve(30)
trojan1.OrbitQuickPlot()

#Mars plot
x0=[Mars.Lx*1.005,Mars.Ly*1.005,0,0]
trojan2 = Trojan(Mars,x0)
trojan2.OribtSolve(30)
trojan2.OrbitQuickPlot()
