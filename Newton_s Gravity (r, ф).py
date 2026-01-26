import math
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.integrate import quad

# known values and constants (revolving around TON 618 - supermassive black hole {hereinafter BH})

G = 6.67430 * 10 ** (-11) #Gravitational constant in MS
M = (6.6 * 10 ** 10) * 1.98847 * 10 ** 30 #BH mass in kg
m0 = 555                 #body's mass in kg
a_u = 149.597870700 * 10**9 #astronomical unit in m
r_btw = 10 * 2600 * a_u #distance to BH
T = 15000 * 24 * 3600 #body's revolving period in sec

# data input and calculations

print('input original speed (vector`s module) and eccentricity of the orbit')
v0 = float(input())
e = float(input())

E = -m0*(math.pi*G*M/2/T)**(2/3)
L = ((e**2-1)*(G*M)**2*m0**3/(2*E))**0.5

# finding r_min and r_max

def find_root(r):
    return 2/m0*(E+G*m0*M/r)-L**2/(m0*r)**2
  
# finding r_min

r_initial=G*m0*M/(2*abs(E))*(1-e)
sol=optimize.root(find_root,r_initial,method='lm',options={'xtol':1e-20,'ftol':1e-40,'maxiter': 200000})
r_min=sol.x[0]
print('r_min =', r_min)

# finding r_max

r_initial=G*m0*M/(2*abs(E))*(1+e)
sol=optimize.root(find_root,r_initial,method='lm',options={'xtol':1e-20,'ftol':1e-40,'maxiter': 200000})
r_max=sol.x[0]
print('r_max = ', r_max)


# nondimensionalization when integrating r of r_min and constant's input for optimisation

alpha = 2*E/m0
beta = 2*G*M/r_min
gamma = L/m0/r_min
print('nondimensionalization paramters')
print(alpha,beta,gamma, r_max/r_min)

#when integrating  r is nondimensionalised into r_min

def Phi_R(r1):
    return gamma/(r1**2)/(alpha+beta/r1-gamma**2/r1**2)**0.5
  
def T_R(r1):
    return 1/(alpha+beta/r1-gamma**2/r1**2)**0.5


r=[]
phi=[]
t=[]
N=1000 # number of intervals
dr=(r_max/r_min-1)/N  # nondimensionalised step 
K=10 #number of half-turns

for i in range(0,K*N+1):
    r.append(r_max/r_min)
    phi.append(0)
    t.append(0)

q=1
for j in range(1,K+1):
    if j % 2 !=0:
        for i in range(1,N+1):
            r[q]=r_max/r_min-dr*i
            s=quad(Phi_R,r_max/r_min, r_max/r_min-dr*i, epsabs=1.49e-10, epsrel=1.49e-10)
            phi[q]=-s[0]+phi[(j-1)*N]
            p=quad(T_R,r_max/r_min, r_max/r_min-dr*i, epsabs=1.49e-10, epsrel=1.49e-10)
            t[q]=(-r_min*p[0])+t[(j-1)*N]
            q=q+1
    if j % 2 ==0:
        for i in range(1,N+1):
            r[q]=(r_max/r_min-dr*(N-i))
            s=quad(Phi_R,r_max/r_min-dr*(N-i), r_max/r_min-dr*N, epsabs=1.49e-10, epsrel=1.49e-10)
            phi[q]=-s[0]+phi[(j-1)*N]
            p=quad(T_R,r_max/r_min-dr*(N-i), r_max/r_min-dr*N, epsabs=1.49e-10, epsrel=1.49e-10)
            t[q]=(-r_min*p[0])+t[(j-1)*N]
            q=q+1

print('comparing true number of K * pie, and the one from integration:')
print((phi[len(phi)-1]-K*math.pi)/(K*math.pi))
print('comparing true K * {half a period}, and the one from integration:')
print((t[len(phi)-1]-K*0.5*math.pi*G*m0*M*(m0/2/(abs(E))**3)**0.5)/(K*0.5*math.pi*G*m0*M*(m0/2/(abs(E))**3)**0.5))

#data output: graph of dependence of nondimentionalised r of polar angle

plt.plot(phi,r, 'r',linewidth=1)
plt.minorticks_on()
plt.xlabel(u'$\phi $ ',fontsize=19)
plt.ylabel(u'$r/r_{min} $',fontsize=19) 
plt.tick_params(direction='in',length=6, width=1,bottom=1, top=1, left=1, right=1,labelsize=16)
plt.tick_params(which='minor',direction='in',length=3, width=1,bottom=1, top=1, left=1, right=1)
plt.grid()
plt.show()

#data output: orbit graph in two-dimensional coordinate system (x, y)

print(len(r),len(phi))
x=[]
y=[]
for i in range (0, len(r)):
    x.append(r[i] * math.cos(phi[i]) * r_min)
    y.append(r[i] * math.sin(phi[i]) * r_min)

plt.scatter(0,0)
plt.plot(x,y, 'r',linewidth=1)
plt.minorticks_on()
plt.xlabel(u'$x [m] $ ',fontsize=19)
plt.ylabel(u'$y [m] $',fontsize=19) 
plt.tick_params(direction='in',length=6, width=1,bottom=1, top=1, left=1, right=1,labelsize=16)
plt.tick_params(which='minor',direction='in',length=3, width=1,bottom=1, top=1, left=1, right=1)
plt.grid()
plt.show()

# clearing data

r.clear()
phi.clear()
t.clear()
x.clear()
y.clear()