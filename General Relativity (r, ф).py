import math
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import integrate

# known values (revolving around TON 618 - supermassive black hole {hereinafter BH})

G = 6.67430 * 10 ** (-11) #Gravitational constant in MS
M = (6.6 * 10 ** 10) * 1.98847 * 10 ** 30 #BH mass in kg
m0 = 555 #body's mass in kg
a_u = 149.597870700 * 10 ** 9 #astronomical unit in m
c = 2.99792458 * 10 ** 8 #light speed in m/sec
T1 = 2500 * 24 * 3600 #BH revolving period (it revolves around itself)
T2 = 15000 * 24 * 3600 #body's revolving period (body revolves around the BH)
R = 2600 * a_u #BH radius in m
pie = 3.141592653589793 #number π
r_btw = 10 * R #distance to BH

# data input and calculations

print('enter orbit eccentricity')
e = float(input()) #orbit eccentricity
print('...')

omega = 2 * pie / T1 #angular velocity
a = 0.4 * R ** 2 * omega / c #metric's parameter for revolving
r_g = 2 * G * M / (c**2) #Schwarzschield Radius 
L_star = M * a * r_g # BH angular momentum

E1 = - m0 * (pie * G * M/(2 * T2)) ** (2/3) # satellite energy in the classical limit
E0 = E1 + m0 * c ** 2 # full energy
L = ((e**2-1) * (G*M)**2 * m0**3 / (2 * E1)) ** 0.5 # satellite's angular momentum

#Nondimensionalisation parameteres (on R)

alpha = L / (m0 * c * r_g) # nondimensionalised angular momentum
beta = E0 / (m0 * c ** 2)  # nondimensionalised full energy
a_ = a / R                 # nondimensionalised а
r_g_ = r_g / R             # nondimensionalised r_g
print('Nondimensionalisation parameteres')
print('α =', alpha,'| β =', beta, ' | nondimensionalised а =', a_, ' | nondimensionalised Schwarzschield Radius =', r_g_)
print('...')

# functions input

def r_n(ri):      # Nondimensionalisation of the variable
  return(ri / R)
def delta(r):     # delta
  return((r ** 2) - (r * r_g) + (a ** 2))
def delta_(r):    # Nondimensionalised delta
  return((r ** 2) - (r * r_g_) + (a_ ** 2))
def y(r):         # derrivate r of satellite's time
  return((1/r ** 4) * ((r ** 2 + a ** 2) * E0 - a * c * L) ** 2 - (delta(r)/r ** 4) * ((a * E0 - c * L) ** 2 + c ** 4 * m0 **2 * r **2))
def y_(r):        # transfirmed y(r)
  return((((r ** 2 + a ** 2) * E0 - a * c * L) ** 2) / (delta(r) * ((a * E0 - c * L) ** 2 + c ** 4 * m0 **2 * r **2)) - 1)

#finding minimum and maximum turning points

r_i = G * m0 * M / (2 * abs(E1)) * (1 - e)
solution = optimize.root(y_, r_i, method = 'lm', options = {'xtol':1e-15, 'ftol':1e-40, 'maxiter':200000})
r_min = solution.x[0]
print('minimum turning point')
print(r_min)
print('...')
print('solution for r_min')
print(solution)
print('...')
r_min_ = r_min / R #nondimensionalisation of r_min

r_i = G * m0 * M / (2 * abs(E1)) * (1 + e)
solution = optimize.root(y_, r_i, method = 'lm', options = {'xtol':1e-15, 'ftol':1e-40, 'maxiter':200000})
r_max = solution.x[0]
print('maximum turning point')
print(r_max)
print('...')
print('solution of r_max')
print(solution)
print('...')
r_max_ = r_max / R #nondimensionalisation of r_max

#functions of dimensionalised integration and inaccuracy (error) check

def dPhi(r): #integration of polar angle
  return(((r ** 2) * ((c * L / delta(r)) * (1 - (r_g / r)) + E0 * (r_g * a/(r * delta(r))))) / (((E0 * (r ** 2 + a ** 2) - a * c * L) ** 2 - delta(r) * ((a * E0 - c * L) ** 2 + c ** 4 * m0 ** 2 * r ** 2)) ** 0.5))

def dPhi_(r): #Nondimensionalised integration of polar angle
  return ((r ** 2) * (alpha * r_g_ * (1 - r_g_ / r) / delta_(r) + a_ * beta * r_g_ / (r * delta_(r))) / (((r ** 2 + a_ ** 2) * beta - a_ * alpha * r_g_) ** 2 - delta_(r) * ((a_ * beta - alpha * r_g_) ** 2 + r ** 2)) ** 0.5)

def c_dt(r): #integration of time
  return(((r ** 2) * (-(r_g * a / (r * delta(r))) * L * c + (E0 / delta(r)) * (r ** 2 + a ** 2 + (r_g * a ** 2 / r)))) / (((E0 * (r ** 2 + a ** 2) - a * c * L) ** 2 - delta(r) * ((a * E0 - c * L) ** 2 + c ** 4 * m0 ** 2 * r ** 2)) ** 0.5))

phi_i = integrate.quadrature(dPhi, r_min, r_max)[0]
print('integrating via scipy.integrate')
print('φ =', phi_i)
print('φ - π =', phi_i - pie)
print('...')

#more accurate nondimensionalised integration by Gauss-Chebyshev quadrature method

def f(x):
  r_ = 0.5 * x * (r_max_ - r_min_) + 0.5 * (r_max_ + r_min_)
  l = (alpha * r_g_ * (r_ ** 2 - r_g_ * r_) / delta_(r_) + a_ * beta * r_g_ * r_ / (delta_(r_))) / (((r_ ** 2 + a_ ** 2) * beta - a_ * alpha * r_g_) ** 2 - delta_(r_) * ((a_ * beta - alpha * r_g_) ** 2 + r_ ** 2)) ** 0.5
  return((1 - x ** 2) ** 0.5 * l)

W = 110
support_seq_1 = [f(math.cos((i - 0.5) * pie / W)) for i in range (1, W + 1)]
integ = 0.5 * (r_max_ - r_min_) * pie / W * sum(support_seq_1)

print('integration by Gauss-Chebyshev quadrature method')
print('φ =', integ) 
print('φ - π =', integ - pie)
support_seq_1.clear()
print('...')

#integration and data saves

r = []
phi_st = []
t = []
N = 110 #number of intervals
K = 10 #number of half-turns

r.append(r_max_)  #r_max
phi_st.append(0)  #initial phase

for i in range(1, K * N + 1):
    r.append(0)
    phi_st.append(0)

q = 1
for j in range(1,K + 1):
  if j % 2 != 0:
    for i in range(1, N + 1):
      r[q] = r[q] = 0.5 * math.cos((i - 0.5) * pie / N) * (r_max_ - r_min_) + 0.5 * (r_max_ + r_min_)
      support_seq_1 = [f(math.cos((l - 0.5) * pie / N)) for l in range (1, i + 1)]
      integ = 0.5 * (r_max_ - r_min_) * pie / N * sum(support_seq_1)
      support_seq_1.clear()
      phi_st[q] = integ + phi_st[(j - 1) * N]
      q = q + 1
  if j % 2 ==0:
    for i in range(N, 0, -1):
      r[q] = r[q] = 0.5 * math.cos((i - 0.5) * pie / N) * (r_max_ - r_min_) + 0.5 * (r_max_ + r_min_)
      support_seq_2 = [f(math.cos((l - 0.5) * pie / N)) for l in range (N, i - 1, -1)]
      integ = 0.5 * (r_max_ - r_min_) * pie / N * sum(support_seq_2)
      support_seq_2.clear()
      phi_st[q] = integ + phi_st[(j - 1) * N]
      q = q + 1

print('comparing the result from integration and true value of К * π')
print((phi_st[len(phi_st) - 1] - K * pie))
print('...')
print('difference with 2 * π for 1 turn')
print(3 * pie * r_g_ / ((0.5 * r_min_ + 0.5 * r_max_) * (1 - e ** 2)))

#output of the graph of dependence of polar angle of nondimensionalised r {r / R}

print('to look at the next graph, close this one')

plt.plot(phi_st, r, 'r', linewidth = 1)
plt.minorticks_on()
plt.xlabel(u'$\phi $ ',fontsize=19)
plt.ylabel(u'$r/R_{Sun} $',fontsize=19)
plt.tick_params(direction='in',length=6, width=1,bottom=1, top=1, left=1, right=1,labelsize=16)
plt.tick_params(which='minor',direction='in',length=3, width=1,bottom=1, top=1, left=1, right=1)
plt.grid()
plt.show()

#output of the graph in two-dimensional coordinate system (x, y)

print('number of values of spherical coordinates')
print(len(r), len(phi_st))

x_coor = []
y_coor = []

for i in range(0, len(r)):       #changing coordinate system (r, φ) into coordinate system (x, y)
  x_coor.append(r[i] * math.cos(phi_st[i]) * R)
  y_coor.append(r[i] * math.sin(phi_st[i]) * R)

plt.scatter(0,0)
plt.plot(x_coor, y_coor, 'r', linewidth = 1)
plt.xlabel(u'$x [m] $ ',fontsize=19)
plt.ylabel(u'$y [m] $',fontsize=19) 
plt.tick_params(direction='in',length=6, width=1,bottom=1, top=1, left=1, right=1,labelsize=16)
plt.tick_params(which='minor',direction='in',length=3, width=1,bottom=1, top=1, left=1, right=1)
plt.grid()
plt.show()

# clearing integration data

r.clear()
phi_st.clear()
t.clear()
x_coor.clear()
y_coor.clear()