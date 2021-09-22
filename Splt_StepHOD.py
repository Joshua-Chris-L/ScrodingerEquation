% python code 2 Split Cell

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib qt
from scipy.fftpack import fft, ifft
import pylab
import numpy as np
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import cm

Nx = 512
kk=np.zeros((Nx), dtype=complex)
xx=np.zeros((Nx), dtype = float)
#Definining Parameters
I = complex (0,1)
Lx = 8
tmax = - I*5
Nt = 1000
dt = -I*1.0/Nt

#Defining x and k
x = [i * 2.0 * np.pi * ( Lx/Nx) for i in range(-Nx//2 ,1+ Nx//2)]
k = (1.0/ Lx )*  np.array ([I * n for n in list(range(0,Nx//2)) + [0] + list(range(- Nx//2+1 ,0))])
for i in range (Nx ):
    kk[i] = abs(np.real( k[i] ** 2))/2.0
    xx[i]=x[i]
Xmin=min(x)
Xmax=max(x)
dx =( Xmax - Xmin )/(Nx-1)

#Defining parameters for psi_0
A = 1
alpha =1/2
f = 1
if f==1:
    f = 1
elif f==2:
    f=2*np.sqrt(alpha)*xx
else:
    f=2-4*alpha*xx**2
   
#Initial wave function
psi_0 = A*f*np.exp(-0.5*alpha*xx**2)
v=np.zeros((Nx), dtype=complex)
#normalising the initial wave function
psi_bar= abs(psi_0)**2
norm = 0
for i in range(Nx):
    norm += psi_bar[i]*dx
Npsi_0 = psi_0/(norm)**0.5

#plot for initial wave function
fig1 =plt.figure()
ax= fig1.add_subplot(111)
ax.plot(xx ,abs(Npsi_0) ** 2 , 'g-', linewidth =3)
plt.title('probability density')
plt.xlabel (r'$x[a.u.]$')
plt.ylabel (r'$| Psi | ^ 2$')
plt.grid()
plt.show ()

#definining the potential

def pot(xx,p):
    return p*(0.5*xx**2)

V= pot(xx,0.03)
I= complex (0,1)

def splitstep(Npsi_0):
    for i in range(Nx):
        Npsi_0[i] = np.exp(-I*dt*V[i]*0.5)*Npsi_0[i]
        c = np.fft.fftshift(np.fft.fft(Npsi_0))
    for i in range(Nx):
        c[i] = c[i]*np.exp(-I*dt*kk[i])
    # Convert back to physical space
        N_psi = np.fft.ifft(np.fft.fftshift(c));
        psi = np.exp(-I*dt*V*0.5)*N_psi
    return psi

#function for 2d plot at every time interval
def plot2D(Y, L, N):
    plt.plot(np.linspace(0, L, num=N), abs(Y), 'r-', linewidth = 3)
plt.show()

#set up plot at everytime interval
psiEv = np.zeros(shape=(Nt, Nx))
for i in range(Nt):
    psi = splitstep(Npsi_0)
    psiEv[i] = abs(psi)
    plot2D(psi, Lx, Nx) # Plot 2D graph at every time step
    plt.show()

T=1
dtt = 0.01

def plot3D(psiEv, Lx, Nx, T, dtt):
    fig = plt.figure(figsize=(20,10))
    ax = fig.gca(projection='3d')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#plot3D(psiEv, Lx, Nx, T, dtt)
fig = plt.figure(figsize=(20,10))
ax = fig.gca(projection='3d')

# Make data.
x = np.arange(0, Lx, Nx)
y = np.arange(0, T, Nx)
X, Y = np.meshgrid(x, y)

# Plot the surface.
surf = ax.plot_surface(X, Y, psiEv, cmap=cm.jet, linewidth=10, antialiased=True, rstride=1,cstride=1)

# Customize the z axis.
#ax.set_zlim(0, 1)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.view_init(30, 190)

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=10)
plt.show()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HOD WOrk

##Propagating the non-linear schrodinger equation using the split operator method
import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator , FormatStrFormatter
import pylab
#% matplotlib inline

Nx = 512*2
#Defining arrays to store values
kk=np.zeros ((Nx), dtype=float )
xx=np.zeros ((Nx), dtype=float )
psi =np.zeros ((Nx), dtype= complex )
pot=np.zeros ((Nx), dtype= float)
U1=np.zeros ((Nx), dtype= complex )
V1=np.zeros ((Nx), dtype= complex )
V2=np.zeros ((Nx), dtype= complex )
U2=np.zeros ((Nx), dtype= complex )
psi_t=np.zeros ((Nx), dtype= complex )
psiii=np.zeros ((Nx), dtype= complex )



#defining Parameter
I= complex (0,1)
Lx =16.0
Nt =1000
tmax= - I * 10.0
dt=tmax/Nt
gap =100
numplots =Nt/gap

#Defining x and k
x = [i * 2.0 * np.pi * ( Lx/Nx) for i in range( - Nx //2 ,1+ Nx //2)]
k = (1.0/ Lx )*  np.array ([I*n for n in list(range (0,Nx //2)) + [0] + list(range( - Nx //2+1 ,0))])
Xmin=min(x)
Xmax=max(x)
dx =( Xmax - Xmin )/Nx
for i in range (Nx ):
    kk[i] = abs(np.real( k[i ] ** 2)/2.0)
    xx[i]=x[i]


#Defining initial wave function
alpha =1.0/20
A=1
gWf= A*np.exp((-alpha * xx ** 2)/2)
agWf =abs(gWf)**2
norm =0


for i in range (Nx):
    norm += agWf[i ]*dx
psi = gWf/(norm) **0.5

#plot the initial wave function
fig4=plt.figure()
ax= fig4.add_subplot(111)
ax.plot(xx ,abs(psi) ** 2 , 'g-', linewidth =3)
plt.title('probability density')
plt.xlabel ('x[a.u.]')
plt.ylabel (r'$| Psi | ^ 2$')
plt.grid()
plt.show ()

#array to store values at every time interval
itt= dt*np.array(range(Nt))
iXX,iTT = pylab.meshgrid(xx,itt)for nt in range ( 10 ):
    for n in range (gap):
        uphi=psi; itim=plotgap*nt+n
        pot =(0.5*xx**2) #+ abs(psi)**2
        U1=np.exp( complex (0, - 0.5) * dt * pot )*psi
        V1=np.fft.fftshift(np.fft.fft(U1))
        V2=np.exp( complex (0, -1) * dt * kk ) * V1
        U2=np.fft.ifft(np.fft.ifftshift(V2))
        psi_t=np.exp(complex (0, -0.5) * dt * pot ) * U2
        psi=psi_t
        ux2=abs(psi)**2
       
        Ekin =0; Epot =0; norm2 =0;
        T1=np.fft.fftshift(np.fft.fft(psi)); T2=np.fft.ifft(np.fft.ifftshift(kk * T1)); T3=np.conjugate(psi)*T2;
        for i in range (Nx):
            Ekin += T3[i ] * dx
            Epot += pot[i ] * ux2[i ] * dx
            norm2 += ux2[i ] * dx  
           # phi[itim,i] = uphi[i]
        psi = psi/np.sqrt(norm2)
        Etot=np.real(Epot+Ekin)
        Etott.append (Etot)
        norm2t.append (norm2 ** 0.5)
        t+=dt
        temps.append (abs(t))
psiii = psi/np.sqrt(norm2)

phi=np.zeros_like(iXX, dtype=complex)
#parameter to store values
Etott =[]; norm2t =[]; temps =[];

plotgap = 10;t=0*I

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
plt.ylim(0.4,1)
plt.title('Average Energy plot')
plt.xlabel ('kick time')
plt.ylabel ('Average Energy')
plt.grid()
plt.show ()

%#Real time propagation with a time dependent potential
gapp=10
plotgap = 10; t=0*I
dtt=0.0001
#array to store values at every time interval
itt= dt*np.array(range(Nt))
iXX,iTT = pylab.meshgrid(xx,itt)
phii=np.zeros_like(iXX, dtype=complex)

#Defining arrays to store values
pot_t=np.zeros ((Nx), dtype= float)
U11=np.zeros ((Nx), dtype= complex )
V11=np.zeros ((Nx), dtype= complex )
V12=np.zeros ((Nx), dtype= complex )
U12=np.zeros ((Nx), dtype= complex )
psi_tt=np.zeros ((Nx), dtype= complex )
psii=np.zeros ((Nx), dtype= complex )
psim=np.zeros ((Nx), dtype= complex )
psimm=np.zeros ((Nx), dtype= complex )

VV=np.zeros((Nx), dtype=float)

def pott(xx,t,phy):
    F_0=7
    w=800
    V = 0.5*xx**2
    F_t=V-xx*F_0*np.sin(w*t + phy)
    return F_t

t=10.0
Etottt =[]; norm2tt =[]; tempss =[];
Np=10
for nt in range ( Np ):
    for n in range (gapp):
        uphi=psiii; itim=plotgap*nt+n
        pot_t = pott(xx,10,10)
        U11=np.exp( complex (0, - 0.5) * dtt * pot_t )*psiii
        V11=np.fft.fftshift(np.fft.fft(U11))
        V12=np.exp( complex (0, -1) * dtt * kk ) * V11
        U12=np.fft.ifft(np.fft.ifftshift(V12))
        psi_tt=np.exp(complex (0, -0.5) * dtt * pot_t ) * U12
        psiii=psi_tt
        psim=np.fft.ifft(np.fft.ifftshift(psi_tt))
        psimm=psim
        ux2=abs(psiii)**2
       
        Ekin =0; Epot =0; norm2 =0;
        T1=np.fft.fftshift(np.fft.fft(psii)); T2=np.fft.ifft(np.fft.ifftshift(kk * T1)); T3=np.conjugate(psii)*T2;
        for i in range (Nx):
            Ekin += T3[i ] * dx
            Epot += pot[i ] * ux2[i ] * dx
            norm2 += ux2[i ] * dx  
            phii[itim,i] = uphi[i]
        psii = psiii/np.sqrt(norm2)
        psim = psim/np.sqrt(norm2)
        Etot=np.real(Epot+Ekin)
        Etottt.append (Etot)
        norm2tt.append (norm2 ** 0.5)
        t+=dt
        tempss.append (abs(t))

fig12 =plt.figure ()
ax= fig12.add_subplot (111)
ax.plot(xx ,abs(psii)**2, 'r',linewidth =3)
#ax.plot(xx ,pot_t, 'r',linewidth =3)
#plt.ylim(0,1)
plt.xlim(-20,20)
plt.title('Probability density in position space')
plt.xlabel ('x [a.u.]')
plt.ylabel ('$ | \ psi | ^ 2 $')
plt.grid()
plt.show ()

fig13 =plt.figure ()
ax= fig13.add_subplot (111)
ax.plot(xx ,abs(psim)**2, 'r',linewidth =3)
#ax.plot(xx ,abs(pot_t)**2, 'r',linewidth =3)
#plt.ylim(0,1)
plt.xlim(-20,20)
plt.title('Probability density in momentum space')
plt.xlabel ('x [a.u.]')
plt.ylabel ('$ | \ psi | ^ 2 $')
plt.grid()
plt.show ()


fig11 =plt.figure ()
ax= fig11.add_subplot (111)
ax.plot(np.linspace(0,1,100),Etottt,'r',linewidth =3)
#plt.ylim(0.4,1)
plt.title('Average Energy plot')
plt.xlabel ('time')
plt.ylabel ('Average Energy')
plt.grid()

plt.show ()




