import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator , FormatStrFormatter
import pylab
#% matplotlib inline

Nx = 512*2
#Defining arrays to store values
k2xm=np.zeros ((Nx), dtype=float )
xx=np.zeros ((Nx), dtype=float )
pot=np.zeros ((Nx), dtype= float)
u=np.zeros ((Nx), dtype= complex )
una=np.zeros ((Nx), dtype= complex )
unb=np.zeros ((Nx), dtype= complex )
v=np.zeros ((Nx), dtype= complex )
vna=np.zeros ((Nx), dtype= complex )
vnb=np.zeros ((Nx), dtype= complex )
psi =np.zeros ((Nx), dtype= complex )


#defining Parameter
ic= complex (0,1)
Lx =16.0
Nt =1000
tmax= - ic * 10.0
dt=tmax/Nt
gap =100
numplots =Nt/gap

#Defining x and k
x = [i * 2.0 * np.pi*( Lx/Nx) for i in range( - Nx //2 ,1+ Nx //2)]
k_x = (1.0/ Lx )*  np.array ([ic * n for n in list(range (0,Nx //2)) + [0] + list(range( - Nx //2+1 ,0))])
Xmin=min(x)
Xmax=max(x)
dx =( Xmax - Xmin )/Nx
for i in range (Nx ):
    k2xm[i] = abs(np.real( k_x [i ] ** 2))
    xx[i]=x[i]


#Defining initial wave function
alpha =1.0/20
ampl =1.0
uinit =ampl * np.exp(( - alpha * xx ** 2)/2.0)
usinit =abs(uinit ) ** 2
norm2i =0
for i in range (Nx):
    norm2i += usinit [i ] * dx
u= uinit /( norm2i) **0.5

#plot the initial wave function
fig4=plt.figure()
ax= fig4.add_subplot(111)
ax.plot(xx ,abs(u ) ** 2 , 'g-', linewidth =3)
plt.title('probability density')
plt.xlabel ('x[a.u.]')
plt.ylabel ('| Psi | ^ 2')
plt.grid()
plt.show ()

#wave function
itt= dt*np.array(range(Nt))
iXX,iTT = pylab.meshgrid(xx,itt)
phi=np.zeros_like(iXX, dtype=complex)
#parameter to store values
Etott =[]; norm2t =[]; temps =[];

t=0.0
fr = 10
p=1
plotgap = 10
for nt in range ( 10 ):
    for n in range (gap):
        uphi=u; itim=plotgap*nt+n
        pot =(0.5*xx**2)
        una=np.exp( complex (0, - 0.5) * dt * pot )*u
        vna=np.fft.fftn(una)
        vnb=np.exp( complex (0, -0.5) * dt * k2xm ) * vna
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
            phi[itim,i] = uphi[i]
        u = u/np.sqrt(norm2)
        Etot=np.real(Epot+Ekin)
        Etott.append (Etot)
        norm2t.append (norm2 ** 0.5)
        t+=dt
        temps.append (abs(t))
       
fig12 =plt.figure ()
ax= fig12.add_subplot (111)
ax.plot(xx ,abs(u)**2, 'r',linewidth =3)
ax.plot(xx ,abs(pot)**2, 'r',linewidth =3)
plt.ylim(0,1)
plt.title('Probability density')
plt.xlabel ('x [a.u.]')
plt.ylabel ('$ | \ psi | ^ 2 $')
plt.grid()
plt.show ()

fig11 =plt.figure ()
ax= fig11.add_subplot (111)
ax.plot(temps,Etott,'r',linewidth =3)
plt.title('Average Energy plot')
plt.xlabel ('kick time')
plt.ylabel ('Average Energy')
plt.grid()
plt.show ()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator , FormatStrFormatter
import pylab
#% matplotlib inline

Nx = 512*2
#Defining arrays to store values
k2xm=np.zeros ((Nx), dtype=float )
xx=np.zeros ((Nx), dtype=float )
pot=np.zeros ((Nx), dtype= float)
u=np.zeros ((Nx), dtype= complex )
una=np.zeros ((Nx), dtype= complex )
unb=np.zeros ((Nx), dtype= complex )
v=np.zeros ((Nx), dtype= complex )
vna=np.zeros ((Nx), dtype= complex )
vnb=np.zeros ((Nx), dtype= complex )
psi =np.zeros ((Nx), dtype= complex )
psiii=np.zeros ((Nx), dtype=float )

#defining Parameter
ic= complex (0,1)
Lx =16.0
Nt =1000
tmax= - ic * 10.0
dt=tmax/Nt
gap =100
numplots =Nt/gap

#Defining x and k
x = [i * 2.0 * np.pi * ( Lx/Nx) for i in range( - Nx //2 ,1+ Nx //2)]
k_x = (1.0/ Lx )*  np.array ([ic * n for n in list(range (0,Nx //2)) + [0] + list(range( - Nx //2+1 ,0))])
Xmin=min(x)
Xmax=max(x)
dx =( Xmax - Xmin )/Nx
for i in range (Nx ):
    k2xm[i] = abs(np.real( k_x [i ] ** 2))
    xx[i]=x[i]


#Defining initial wave function
alpha =1.0/20
ampl =1.0
uinit =ampl * np.exp(( - alpha * xx ** 2)/2)
usinit =abs(uinit ) ** 2
norm2i =0
for i in range (Nx):
    norm2i += usinit [i ] * dx
u= uinit /( norm2i) **0.5
psiii=uinit /( norm2i) **0.5
#plot the initial wave function
fig4=plt.figure()
ax= fig4.add_subplot(111)
ax.plot(xx ,abs(u ) ** 2 , 'g-', linewidth =3)
plt.title('probability density')
plt.xlabel ('x[a.u.]')
plt.ylabel ('| Psi | ^ 2')
plt.grid()
plt.show ()

#wave function
itt= dt*np.array(range(Nt))
iXX,iTT = pylab.meshgrid(xx,itt)
phi=np.zeros_like(iXX, dtype=complex)
#parameter to store values
Etott =[]; norm2t =[]; temps =[];

t=0.0
fr = 10
p=1
plotgap = 10
for nt in range ( 10 ):
    for n in range (gap):
        uphi=u; itim=plotgap*nt+n
        pot = 0.5*xx**2
        una=np.exp( complex (0, - 0.5) * dt * pot )*u
        vna=np.fft.fftn(una)
        vnb=np.exp( complex (0, -0.5) * dt * k2xm ) * vna
        unb=np.fft.ifftn(vnb)
        v=np.exp( complex (0, -0.5) * dt * pot ) * unb  
        u=v
        up2=abs(np.fft.fftn(u )) ** 2
        ux2=abs(u)**2
   
        Ekin =0; Epot =0; norm2 =0;
        T1=np.fft.fftn(psi); T2=np.fft.ifftn(k2xm * T1); T3=np.conjugate (psi)*T2;
        for i in range (Nx):
            Ekin += T3[i ] * dx
            Epot += pot[i ] * ux2[i ] * dx
            norm2 += ux2[i ] * dx  
            phi[itim,i] = uphi[i]
        psi = u/np.sqrt(norm2)
        Etot=np.real(Epot+Ekin)
        Etott.append (Etot)
        norm2t.append (norm2 ** 0.5)
        t+=dt
        temps.append (abs(t))

fig12 =plt.figure ()
ax= fig12.add_subplot (111)
ax.plot(xx ,abs(psi)**2, 'r',linewidth =3)
ax.plot(xx ,abs(pot)**2, 'r',linewidth =3)
plt.ylim(0,1)
plt.xlim(-20,20)
plt.title('Probability density')
plt.xlabel ('x [a.u.]')
plt.ylabel ('$ | \ psi | ^ 2 $')
plt.grid()
plt.show ()

fig11 =plt.figure ()
ax= fig11.add_subplot (111)
ax.plot(temps,Etott,'r',linewidth =3)
plt.title('Average Energy plot')
plt.xlabel ('kick time')
plt.ylabel ('Average Energy')
plt.grid()
