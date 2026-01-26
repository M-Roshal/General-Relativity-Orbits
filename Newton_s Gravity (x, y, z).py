import matplotlib.pyplot as plt
import numpy as np

ax = plt.axes(projection = '3d')

# known constants and input


G = 6.67430 * 10**(-11) #m3·s−2·kg−1
print('Input original speed of the body by 3 axes respectively')
v0_x = int(input())
v0_y = int(input()) 
v0_z = int(input()) # original speed in m/sec in 3 axes
a_u = 149.597870700 * 10**9 # astronomical unit in meters
M = 1.98847 * 10**30 # star mass in kg

print('The graph of the trajectory of the body in NEwtonian Gravity is shown. Single segment - 10^11 м. Graph is three-dimensional')

# variables

v2=np.empty(3) # body's mass
r2=np.empty(3) # body's coordinates

# values for integration

dt = 3600 #sec, integration step (measures every hour)
t0 = 0 #sec, current time
T = 5*365*24*3600 #sec, time when measurements stop

# Original coordinates for Runge-Kutta integration method (4th order)
r2[0] = a_u #x
r2[1] = 0 #y
r2[2] = 0 #z

# Original speeds for Runge-Kutta integration method (4th order)
v2[0] = v0_x #x
v2[1] = v0_y #y
v2[2] = v0_z #z

# integration via Runge-Kutta method (4th order)

def GravityForceX(x,y,z,vx,vy,vz):
    return -G*M*x/(x**2+y**2+z**2)**1.5
def GravityForceY(x,y,z,vx,vy,vz):
    return -G*M*y/(x**2+y**2+z**2)**1.5
def GravityForceZ(x,y,z,vx,vy,vz):
    return -G*M*z/(x**2+y**2+z**2)**1.5

def a_rX(x,y,z,vx,vy,vz):
    return vx*dt
def a_rY(x,y,z,vx,vy,vz):
    return vy*dt
def a_rZ(x,y,z,vx,vy,vz):
    return vz*dt

def a_vX(x,y,z,vx,vy,vz):
    return dt*GravityForceX(x,y,z,vx,vy,vz)
def a_vY(x,y,z,vx,vy,vz):
    return dt*GravityForceY(x,y,z,vx,vy,vz)
def a_vZ(x,y,z,vx,vy,vz):
    return dt*GravityForceZ(x,y,z,vx,vy,vz)

def b_rX(x,y,z,vx,vy,vz):
    return (vx+0.5*a_vX(x,y,z,vx,vy,vz))*dt
def b_rY(x,y,z,vx,vy,vz):
    return (vy+0.5*a_vY(x,y,z,vx,vy,vz))*dt
def b_rZ(x,y,z,vx,vy,vz):
    return (vz+0.5*a_vZ(x,y,z,vx,vy,vz))*dt

def b_vX(x,y,z,vx,vy,vz):
    return dt*GravityForceX(x+0.5*a_rX(x,y,z,vx,vy,vz), y+0.5*a_rY(x,y,z,vx,vy,vz), z+0.5*a_rZ(x,y,z,vx,vy,vz), vx+0.5*a_vX(x,y,z,vx,vy,vz), vy+0.5*a_vY(x,y,z,vx,vy,vz), vz+0.5*a_vZ(x,y,z,vx,vy,vz))
def b_vY(x,y,z,vx,vy,vz):
    return dt*GravityForceY(x+0.5*a_rX(x,y,z,vx,vy,vz), y+0.5*a_rY(x,y,z,vx,vy,vz), z+0.5*a_rZ(x,y,z,vx,vy,vz), vx+0.5*a_vX(x,y,z,vx,vy,vz), vy+0.5*a_vY(x,y,z,vx,vy,vz), vz+0.5*a_vZ(x,y,z,vx,vy,vz))
def b_vZ(x,y,z,vx,vy,vz):
    return dt*GravityForceZ(x+0.5*a_rX(x,y,z,vx,vy,vz), y+0.5*a_rY(x,y,z,vx,vy,vz), z+0.5*a_rZ(x,y,z,vx,vy,vz), vx+0.5*a_vX(x,y,z,vx,vy,vz), vy+0.5*a_vY(x,y,z,vx,vy,vz), vz+0.5*a_vZ(x,y,z,vx,vy,vz))

def c_rX(x,y,z,vx,vy,vz):
    return (vx+0.5*b_vX(x,y,z,vx,vy,vz))*dt
def c_rY(x,y,z,vx,vy,vz):
    return (vy+0.5*b_vY(x,y,z,vx,vy,vz))*dt
def c_rZ(x,y,z,vx,vy,vz):
    return (vz+0.5*b_vZ(x,y,z,vx,vy,vz))*dt

def c_vX(x,y,z,vx,vy,vz):
    return dt*GravityForceX(x+0.5*b_rX(x,y,z,vx,vy,vz), y+0.5*b_rY(x,y,z,vx,vy,vz), z+0.5*b_rZ(x,y,z,vx,vy,vz), vx+0.5*b_vX(x,y,z,vx,vy,vz), vy+0.5*b_vY(x,y,z,vx,vy,vz), vz+0.5*b_vZ(x,y,z,vx,vy,vz))
def c_vY(x,y,z,vx,vy,vz):
    return dt*GravityForceY(x+0.5*b_rX(x,y,z,vx,vy,vz), y+0.5*b_rY(x,y,z,vx,vy,vz), z+0.5*b_rZ(x,y,z,vx,vy,vz), vx+0.5*b_vX(x,y,z,vx,vy,vz), vy+0.5*b_vY(x,y,z,vx,vy,vz), vz+0.5*b_vZ(x,y,z,vx,vy,vz))
def c_vZ(x,y,z,vx,vy,vz):
    return dt*GravityForceZ(x+0.5*b_rX(x,y,z,vx,vy,vz), y+0.5*b_rY(x,y,z,vx,vy,vz), z+0.5*b_rZ(x,y,z,vx,vy,vz), vx+0.5*b_vX(x,y,z,vx,vy,vz), vy+0.5*b_vY(x,y,z,vx,vy,vz), vz+0.5*b_vZ(x,y,z,vx,vy,vz))

def d_rX(x,y,z,vx,vy,vz):
    return (vx+c_vX(x,y,z,vx,vy,vz))*dt
def d_rY(x,y,z,vx,vy,vz):
    return (vy+c_vY(x,y,z,vx,vy,vz))*dt
def d_rZ(x,y,z,vx,vy,vz):
    return (vz+c_vZ(x,y,z,vx,vy,vz))*dt

def d_vX(x,y,z,vx,vy,vz):
    return dt*GravityForceX(x+c_rX(x,y,z,vx,vy,vz), y+c_rY(x,y,z,vx,vy,vz), z+c_rZ(x,y,z,vx,vy,vz), vx+c_vX(x,y,z,vx,vy,vz), vy+c_vY(x,y,z,vx,vy,vz), vz+c_vZ(x,y,z,vx,vy,vz))
def d_vY(x,y,z,vx,vy,vz):
    return dt*GravityForceY(x+c_rX(x,y,z,vx,vy,vz), y+c_rY(x,y,z,vx,vy,vz), z+c_rZ(x,y,z,vx,vy,vz), vx+c_vX(x,y,z,vx,vy,vz), vy+c_vY(x,y,z,vx,vy,vz), vz+c_vZ(x,y,z,vx,vy,vz))
def d_vZ(x,y,z,vx,vy,vz):
    return dt*GravityForceZ(x+c_rX(x,y,z,vx,vy,vz), y+c_rY(x,y,z,vx,vy,vz), z+c_rZ(x,y,z,vx,vy,vz), vx+c_vX(x,y,z,vx,vy,vz), vy+c_vY(x,y,z,vx,vy,vz), vz+c_vZ(x,y,z,vx,vy,vz))


filename2 = 'Object_orbit_Classic_Theory_of_Gravitaty.dat' 
outfile2 = open(filename2, 'w')

#saves and output

while t0<T:
    outfile2.write ('%4.9f %4.9f %4.9f %4.3f \n' % (r2[0], r2[1], r2[2],  t0))

    l2=r2[0]
    m2=r2[1]
    k2=r2[2]
    v_1=v2[0]
    v_2=v2[1]
    v_3=v2[2]
    r2[0]=r2[0]+(a_rX(l2,m2,k2,v_1,v_2,v_3)+2*b_rX(l2,m2,k2,v_1,v_2,v_3)+2*c_rX(l2,m2,k2,v_1,v_2,v_3)+d_rX(l2,m2,k2,v_1,v_2,v_3))/6
    r2[1]=r2[1]+(a_rY(l2,m2,k2,v_1,v_2,v_3)+2*b_rY(l2,m2,k2,v_1,v_2,v_3)+2*c_rY(l2,m2,k2,v_1,v_2,v_3)+d_rY(l2,m2,k2,v_1,v_2,v_3))/6
    r2[2]=r2[2]+(a_rZ(l2,m2,k2,v_1,v_2,v_3)+2*b_rZ(l2,m2,k2,v_1,v_2,v_3)+2*c_rZ(l2,m2,k2,v_1,v_2,v_3)+d_rZ(l2,m2,k2,v_1,v_2,v_3))/6

    v2[0]=v2[0]+(a_vX(l2,m2,k2,v_1,v_2,v_3)+2*b_vX(l2,m2,k2,v_1,v_2,v_3)+2*c_vX(l2,m2,k2,v_1,v_2,v_3)+d_vX(l2,m2,k2,v_1,v_2,v_3))/6
    v2[1]=v2[1]+(a_vY(l2,m2,k2,v_1,v_2,v_3)+2*b_vY(l2,m2,k2,v_1,v_2,v_3)+2*c_vY(l2,m2,k2,v_1,v_2,v_3)+d_vY(l2,m2,k2,v_1,v_2,v_3))/6
    v2[2]=v2[2]+(a_vZ(l2,m2,k2,v_1,v_2,v_3)+2*b_vZ(l2,m2,k2,v_1,v_2,v_3)+2*c_vZ(l2,m2,k2,v_1,v_2,v_3)+d_vZ(l2,m2,k2,v_1,v_2,v_3))/6
    
    t0=t0+dt


outfile2.close()

x2=[]
y2=[]
z2=[]

f=open('Object_orbit_Classic_Theory_of_Gravitaty.dat', 'r')
for line in f: 
    words = line.split() # splitting the string into words  
    x2.append(float(words[0]))
    y2.append(float(words[1]))
    z2.append(float(words[2]))
f.close()

ax.scatter(0, 0, 0)
ax.plot(x2,y2,z2)
ax.grid()
plt.show()