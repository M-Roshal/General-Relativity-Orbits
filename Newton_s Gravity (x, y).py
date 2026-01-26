import matplotlib.pyplot as plt
import numpy as np

# known constants and data input

G = 6.67430 * 10**(-11) #m3·с−2·kg−1
print('input the original speed of the object (m / sec) for the two axes respectively')
v0_x = int(input())
v0_y = int(input()) # original speed in m/sec
a_u = 149.597870700 * 10**9 #astronomical unit in m
M = 1.98847 * 10**30 #star mass in kg

print('The particle orbit in Newtonian Gravity is illustrated in the output. Single segment - 10^11 m. Graph is two-dimensional')

# variables

v2 = np.empty(3) # objects' speed
r2 = np.empty(3) # objects' coordinates

# variables for integration

dt = 3600 #sec, step of integration (measures every hour)
t0 = 0 #sec, current time
T = 5*365*24*3600 #sed, time when measurements stop

#Original coordinates for Runge-Kutta integration method (4th order)
r2[0] = a_u #x
r2[1] = 0 #y

#Original speeds for Runge-Kutta integration method (4th order)
v2[0] = v0_x #x
v2[1] = v0_y #y

#integration via Runge-Kutta method (4th order)

def GravityForceX(x,y,vx,vy):
    return -G*M*x/(x**2+y**2)**1.5
def GravityForceY(x,y,vx,vy):
    return -G*M*y/(x**2+y**2)**1.5

def a_rX(x,y,vx,vy):
    return vx*dt
def a_rY(x,y,vx,vy):
    return vy*dt

def a_vX(x,y,vx,vy):
    return dt*GravityForceX(x,y,vx,vy)
def a_vY(x,y,vx,vy):
    return dt*GravityForceY(x,y,vx,vy)

def b_rX(x,y,vx,vy):
    return (vx+0.5*a_vX(x,y,vx,vy))*dt
def b_rY(x,y,vx,vy):
    return (vy+0.5*a_vY(x,y,vx,vy))*dt

def b_vX(x,y,vx,vy):
    return dt*GravityForceX(x+0.5*a_rX(x,y,vx,vy), y+0.5*a_rY(x,y,vx,vy),vx+0.5*a_vX(x,y,vx,vy),vy+0.5*a_vY(x,y,vx,vy))
def b_vY(x,y,vx,vy):
    return dt*GravityForceY(x+0.5*a_rX(x,y,vx,vy), y+0.5*a_rY(x,y,vx,vy),vx+0.5*a_vX(x,y,vx,vy),vy+0.5*a_vY(x,y,vx,vy))

def c_rX(x,y,vx,vy):
    return (vx+0.5*b_vX(x,y,vx,vy))*dt
def c_rY(x,y,vx,vy):
    return (vy+0.5*b_vY(x,y,vx,vy))*dt
def c_vX(x,y,vx,vy):
    return dt*GravityForceX(x+0.5*b_rX(x,y,vx,vy), y+0.5*b_rY(x,y,vx,vy),vx+0.5*b_vX(x,y,vx,vy),vy+0.5*b_vY(x,y,vx,vy))
def c_vY(x,y,vx,vy):
    return dt*GravityForceY(x+0.5*b_rX(x,y,vx,vy), y+0.5*b_rY(x,y,vx,vy),vx+0.5*b_vX(x,y,vx,vy),vy+0.5*b_vY(x,y,vx,vy))

def d_rX(x,y,vx,vy):
    return (vx+c_vX(x,y,vx,vy))*dt
def d_rY(x,y,vx,vy):
    return (vy+c_vY(x,y,vx,vy))*dt
def d_vX(x,y,vx,vy):
    return dt*GravityForceX(x+c_rX(x,y,vx,vy), y+c_rY(x,y,vx,vy),vx+c_vX(x,y,vx,vy),vy+c_vY(x,y,vx,vy))
def d_vY(x,y,vx,vy):
    return dt*GravityForceY(x+c_rX(x,y,vx,vy), y+c_rY(x,y,vx,vy),vx+c_vX(x,y,vx,vy),vy+c_vY(x,y,vx,vy))


filename2 = 'Object_orbit_Classic_Theory_of_Gravitaty.dat' 
outfile2 = open(filename2, 'w')

#saving data and output

while t0<T:
    outfile2.write ('%4.9f %4.9f %4.3f \n' % (r2[0], r2[1],  t0))

    l2=r2[0]
    m2=r2[1]
    v_2=v2[0]
    vy2=v2[1]
    r2[0]=r2[0]+(a_rX(l2,m2,v_2,vy2)+2*b_rX(l2,m2,v_2,vy2)+2*c_rX(l2,m2,v_2,vy2)+d_rX(l2,m2,v_2,vy2))/6
    r2[1]=r2[1]+(a_rY(l2,m2,v_2,vy2)+2*b_rY(l2,m2,v_2,vy2)+2*c_rY(l2,m2,v_2,vy2)+d_rY(l2,m2,v_2,vy2))/6

    v2[0]=v2[0]+(a_vX(l2,m2,v_2,vy2)+2*b_vX(l2,m2,v_2,vy2)+2*c_vX(l2,m2,v_2,vy2)+d_vX(l2,m2,v_2,vy2))/6
    v2[1]=v2[1]+(a_vY(l2,m2,v_2,vy2)+2*b_vY(l2,m2,v_2,vy2)+2*c_vY(l2,m2,v_2,vy2)+d_vY(l2,m2,v_2,vy2))/6
    
    t0=t0+dt


outfile2.close()

x2=[]
y2=[]
f=open('Object_orbit_Classic_Theory_of_Gravitaty.dat', 'r')
for line in f: 
    words = line.split() # разбиение строки на слова  
    x2.append(float(words[0]))
    y2.append(float(words[1]))
f.close()

plt.scatter(0,0)
plt.title('Bodies Trajectory in Newtonian Gravity')
plt.plot(x2,y2,linewidth=1)

plt.minorticks_on()
plt.xlabel(u'$x [m]$ ',fontsize=19)
plt.ylabel(u'$y [m]$',fontsize=19) 
plt.tick_params(direction='in',length=6, width=1,bottom=1, top=1, left=1, right=1,labelsize=16)
plt.tick_params(which='minor',direction='in',length=3, width=1,bottom=1, top=1, left=1, right=1)
plt.grid()
plt.show()