import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator , FormatStrFormatter
import pylab
#% matplotlib inline




ic= complex (0 ,1)

Lx =4.0
Nx =10000
Nt =100
tmax= - ic * 10.0
dt=tmax/Nt
gap =10
numplots =Nt/gap




x = [i * 2.0 * np.pi * ( Lx/Nx) for i in range( - Nx //2 ,1+ Nx //2)]

k_x = (1.0/ Lx )*  np.array ([ic * n for n in list(range (0,Nx //2)) + [0] + list(range( - Nx //2+1 ,0))])

Xmin=min(x)
Xmax=max(x)
dx =( Xmax - Xmin )/Nx


k2xm=np.zeros ((Nx), dtype=float )
xx=np.zeros ((Nx), dtype=float )


for i in range (Nx ):
    k2xm[i] = abs(np.real( k_x [i ] ** 2))
    xx[i]=x[i]

pot=np.zeros ((Nx), dtype= float)
u=np.zeros ((Nx), dtype= complex )

una=np.zeros ((Nx), dtype= complex )
unb=np.zeros ((Nx), dtype= complex )
v=np.zeros ((Nx), dtype= complex )
vna=np.zeros ((Nx), dtype= complex )
vnb=np.zeros ((Nx), dtype= complex )


alpha =1.0
ampl =1.0
uinit =ampl * np.exp(( - alpha * xx ** 2)/2)
usinit =abs(uinit ) #** 2

norm2i =0
for i in range (Nx):
    norm2i += usinit [i ] * dx

u= uinit /( norm2i) #** 0.5)


fig4=plt.figure()
ax= fig4.add_subplot(111)
ax.plot(xx ,abs(u ) ** 2 , 'g-', linewidth =3)
plt.title('probability density')
plt.xlabel ('x[a.u.]')
plt.ylabel ('| Psi | ^ 2')
plt.grid()
plt.legend ()
plt.show ()



tt= dt * np.array(range (Nt))
#iXX ,iTT = pylab.meshgrid (xx ,tt)
#phi=np.zeros_like (iXX , dtype= complex )


usquared =abs(u) ** 2
plotnum =0
Etott =[]; norm2t =[]; temps =[];

K=0.8
t=ic*0.0

vi=2
if vi==1:
    g=0
elif vi == 2:
   
    g=1.5
else:
    g=10
   
fr=10
for nt in range ( 10 ):
    itim= fr
    for n in range (gap):
        uphi=u; itim= gap * nt + n
       
        pot1 =g*(abs(u) ** 2)
        pot2=K*np.cos(xx)
       
        pot= pot1 + pot2

        una=np.exp( complex (0, - 0.5) * dt * pot ) * u
     
        vna=np.fft.fftn(una)
       
        vnb=np.exp( complex (0, -1) * dt * k2xm ) * vna
       
        unb=np.fft.ifftn(vnb)
       
        v=np.exp( complex (0, - 0.5) * dt * pot ) * unb
       
        u=v
       
        up2=abs(np.fft.fftn(u )) ** 2
       
        ux2=abs(u)**2
       
        Ekin =0; Epot =0; norm2 =0;
       
        T1=np.fft.fftn(u); T2=np.fft.ifftn(k2xm * T1); T3=np.conjugate (u)*T2;
       
       
        for i in range (Nx):
            Ekin += T3[i ] * dx
            Epot += pot[i ] * ux2[i ] * dx
            norm2 += ux2[i ] * dx
        #    phi[itim ,i] = uphi[i]
       
        u=u/np.sqrt(norm2)
        Etot=np.real(Epot+Ekin)
        Etott.append (Etot)
     
        norm2t.append (norm2 ** 0.5)
        temps.append (abs(t))
        t+=dt
        plotnum +=1

ktime=np.linspace(0,100,100)



fig12 =plt.figure ()
ax= fig12.add_subplot (111)
plt.cla ()
#ax.plot(xx ,np.imag(u ) ** 2+ np.real(u ) ** 2 , 'g', linewidth =3)
ax.plot(xx ,abs(u)**2, 'r',linewidth =3)
plt.title('Probability density')
plt.xlabel ('x [a.u.]')
plt.ylabel ('$ | \ psi | ^ 2 $')

plt.grid()
plt.show ()

fig12 =plt.figure ()
ax= fig12.add_subplot (111)
plt.cla ()
#ax.plot(xx ,np.imag(u ) ** 2+ np.real(u ) ** 2 , 'g', linewidth =3)
ax.plot(ktime,Etott,'r',linewidth =3)
plt.title('Average Energy plot')
plt.xlabel ('kick time')
plt.ylabel ('Average Energy')
plt.grid()
plt.show ()




