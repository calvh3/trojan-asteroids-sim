from solver import *

'''
Libration vs Planet mass plot
    -Simple plot showing how libration changes with planet-sun mass ratio
    -Shows the analytical result of lagrange point stability breaking down at around mp/ms = 0.02
    -(Once libration approaches 2pi (marked in plot), particle is no longer bound).
'''



W = []
Mps = []
Mps = np.logspace(-10, np.log10(0.025), 100)

for Mp in Mps:
     #print(Mp,'Computed')
     p = Planet(Mp,5)
     t = Trojan(p,[p.Lx,p.Ly,0,0])
     t.OribtSolve(100)
     wander = t.Wander(t.X[0],t.X[1])
     W.append(wander)
     print("Mass = ", Mp,"wander = ",wander)
     
     
plt.plot(Mps/(1+Mps),W,'k',linewidth=1)
plt.plot([1e-10,0.025],[2*np.pi,2*np.pi],'k--',linewidth=0.75)
plt.text(1e-10,2.2*np.pi,'libration = 2π')
plt.xlabel("Mp/(Ms+Mp)")
plt.ylabel("Libration / rad")
plt.xscale('log')
plt.yscale('log')
plt.show()
