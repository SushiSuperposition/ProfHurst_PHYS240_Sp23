# Three-body/ Kozai Project

'''
The purpose of this code is to explore the famous three-body problem in Classical Mechanics.
The three-body problem is the problem of taking the intial positions and velocities (or momenta)
of three point masses and solving for their subsequen motion according to Newton's laws
of motion and Newton's law of universal gravitation. The system is chaotic for most intial conditions,
and numerical methods are generally required.
'''

# Equations of Motion
'''
Newton’s Law of Gravitation says that any two point masses have an attractive force between them 
(called the gravitational force), the magnitude of which is directly proportional to the product 
of their masses and inversely proportional to the square of the distance between them. The equation below represents 
this law in vector form:

F = G * Mm/r^2 (Unit vector points away from body M to m) (Equation 1)

According to Newton's second law of motion, the net force on an object produces
a net change in momentum of the object (i.e. F = MA). So, applying the above equation
to the body having mass M, we get the following differential equation of motion for the body:

M d^2r_1/dt^2 = G * Mm/r^3 (second-order differential equation) (Equation 2)

The acceleration of an object is the change in velocity of the object with time so the second order
differential of position can be replaced with a first order differential of velocity. Similarly,
the velocity can be expressed as a first order differential of the position:

m_i * dv_i/dt = G * m_im_j/r^3_ij (Equation 3)

dr_i/dt = v_i (Equation 4)
'''


# Functions used in project

# Force of Gravity of Sun on Earth
def force_sun_earth(r):
    """
    Calculate the gravitational force vector between Earth and the Sun based on the given position vector.

    Parameters:
        r: numpy.array, shape (3,), position vector representing the position of Earth with respect to the Sun in three dimensions.

    Returns:
        F: numpy.array, shape (3,), gravitational force vector between Earth and the Sun in three dimensions.
    """

    F = np.zeros(3)  # Force array for x,y,z force values
    F_mag = G * Ms * Me / np.linalg.norm(r) ** 2  # Magnitude of gravitational force

    theta = math.atan2(np.abs(r[1]), np.abs(r[0]))  # Calculate the angle theta in the x-y plane
    phi = math.atan2(np.abs(r[2]), np.sqrt(r[0] ** 2 + r[1] ** 2))  # Calculate the angle phi with respect to the z-axis

    F[0] = F_mag * np.cos(theta)  # Calculate the x-component of the force
    F[1] = F_mag * np.sin(theta)  # Calculate the y-component of the force
    F[2] = F_mag * np.sin(phi)  # Calculate the z-component of the force

    if r[0] > 0:
        F[0] = -F[0]  # If the x-coordinate is positive, reverse the sign of the x-component of the force

    if r[1] > 0:
        F[1] = -F[1]  # If the y-coordinate is positive, reverse the sign of the y-component of the force

    if r[2] > 0:
        F[2] = -F[2]  # If the z-coordinate is positive, reverse the sign of the z-component of the force

    return F  # Return the force vector F


# Force of Gravity of Sun on Moon
def force_sun_moon(r):
    """
    Calculate the gravitational force vector between Moon and the Sun based on the given position vector.

    Parameters:
        r: numpy.array, shape (3,), position vector representing the position of Moon with respect to the Sun in three dimensions.

    Returns:
        F: numpy.array, shape (3,), gravitational force vector between Moon and the Sun in three dimensions.
    """

    F = np.zeros(3)  # Force array for x,y,z force values
    F_mag = G * Ms * Mm / np.linalg.norm(r) ** 2  # Magnitude of gravitational force

    theta = math.atan2(np.abs(r[1]), np.abs(r[0]))  # Calculate the angle theta in the x-y plane
    phi = math.atan2(np.abs(r[2]), np.sqrt(r[0] ** 2 + r[1] ** 2))  # Calculate the angle phi with respect to the z-axis

    F[0] = F_mag * np.cos(theta)  # Calculate the x-component of the force
    F[1] = F_mag * np.sin(theta)  # Calculate the y-component of the force
    F[2] = F_mag * np.sin(phi)  # Calculate the z-component of the force

    if r[0] > 0:
        F[0] = -F[0]  # If the x-coordinate is positive, reverse the sign of the x-component of the force

    if r[1] > 0:
        F[1] = -F[1]  # If the y-coordinate is positive, reverse the sign of the y-component of the force

    if r[2] > 0:
        F[2] = -F[2]  # If the z-coordinate is positive, reverse the sign of the z-component of the force

    return F  # Return the force vector F


# Force of Gravity of Earth on Moon
def force_earth_moon(re, rm):
    """
    Calculate the gravitational force vector between Moon and the Earth based on the given position vector.

    Parameters:
        re: numpy.array, shape (3,), position vector representing the position of Earth in three dimensions.
        rm: numpy.array, shape (3,), position vector representing the position of Moon in three dimensions.

    Returns:
        F: numpy.array, shape (3,), gravitational force vector between Moon and the Sun in three dimensions.
    """
    r = np.zeros(3)  # Position array for x,y,z position values
    F = np.zeros(3)  # Force array for x,y,z force values

    r[0] = re[0] - rm[0]
    r[1] = re[1] - rm[1]
    r[2] = re[2] - rm[2]
    r_mag = np.linalg.norm(r)

    F_mag = G * Me * Mm / r_mag ** 2  # Magnitude of gravitational force

    theta = math.atan2(np.abs(r[1]), np.abs(r[0]))  # Calculate the angle theta in the x-y plane
    phi = math.atan2(np.abs(r[2]), np.sqrt(r[0] ** 2 + r[1] ** 2))  # Calculate the angle phi with respect to the z-axis

    F[0] = F_mag * np.cos(theta)  # Calculate the x-component of the force
    F[1] = F_mag * np.sin(theta)  # Calculate the y-component of the force
    F[2] = F_mag * np.sin(phi)  # Calculate the z-component of the force

    if r[0] > 0:
        F[0] = -F[0]  # If the x-coordinate is positive, reverse the sign of the x-component of the force

    if r[1] > 0:
        F[1] = -F[1]  # If the y-coordinate is positive, reverse the sign of the y-component of the force

    if r[2] > 0:
        F[2] = -F[2]  # If the z-coordinate is positive, reverse the sign of the z-component of the force

    return F  # Return the force vector F


# Calculate the net force for celestial objects
def force(r, obj, ro, vo):
    """
    Calculate the net force acting on an object.

    Parameters:
        r (float): Distance between the object and theSsun.
        obj (str): Name of the object ('earth' or 'moon').
        ro (float): Distance of the second object from the Sun.
        vo (float): Velocity of the second object.

    Returns:
        float: Net force acting on the object.
    """
    if obj == 'earth':
        return force_sun_earth(r) + force_earth_moon(r, ro)
    if obj == 'moon':
        return force_sun_moon(r) - force_earth_moon(r, ro)


# Velocity of celestial object
def velocity(t, r, v, obj, ro, vo):
    """
    Calculate the velocity of an object.

    Parameters:
        t (float): Time.
        r (float): Distance of the object from the Sun.
        v (float): Current velocity of the object.
        obj (str): Name of the object ('earth' or 'moon').
        ro (float): Distance of the second object from the Sun.
        vo (float): Velocity of the second object.

    Returns:
        float: Updated velocity of the object.
    """

    return v


# Acceleration of celestial object
def acceleration(t, r, v, obj, ro, vo):
    """
    Calculate the acceleration of a object.

    Parameters:
        t (float): Time.
        r (float): Distance of the object from the Sun.
        v (float): Velocity of the object.
        obj (str): Name of the object ('earth' or 'moon').
        ro (float): Distance of the second object from the Sun.
        vo (float): Velocity of the second object.

    Returns:
        float: Acceleration of the object.
    """
    F = force(r, obj, ro, vo)
    if obj == 'earth':
        a = F / Me  # Newton's second Law (F=MA)
    if obj == 'moon':
        a = F / Mm  # Newton's second Law (F=MA)
    return a


# Kinetic energy of the Moon
def KineticEnergy(v):
    """
    Calculate the kinetic energy of the Moon.

    Parameters:
        v (array-like): Velocity array of the Moon.

    Returns:
        float: Kinetic energy of the Moon.
    """
    v_mag = np.linalg.norm(v)  # Calculate the norm (magnitude) of the velocity array
    return 0.5 * Mm * v_mag ** 2


# Potential energy of the Moon
def PotentialEnergy(r):
    """
    Calculate the potential energy of the Moon.

    Parameters:
        r (array-like): Radius array of the Moon.

    Returns:
        float: Potential energy of the Moon.
    """
    f_mag = np.linalg.norm(force_sun_moon(r))  # Calculate the norm (magnitude) of the force array
    r_mag = np.linalg.norm(r)  # Calculate the norm (magnitude) of the radius array
    return -f_mag * r_mag


# Angular Momentum of the Moon
def AngularMomentum(r, v, theta):
    """
    Calculate the angular momentum of the Moon.

    Parameters:
        r (array-like): Radius array of the Moon.
        v (array-like): Velocity array of the Moon.
        theta (float): Angle between the radius vector and velocity vector.

    Returns:
        float: Angular momentum of the Moon.
    """
    r_mag = np.linalg.norm(r)  # Calculate the norm (magnitude) of the radius.
    v_mag = np.linalg.norm(v)  # Calculate the norm (magnitude) of the velocity.
    return Mm * r_mag * v_mag * np.sin(theta)  # Angular momentum equation: L = Mm * r * v * sin(theta)


# Function for eccentricity
def Eccentricity(r):
    """
    Calculate the eccentricity and related parameters of the orbit.

    Parameters:
        r (array-like): Array of radius vectors.

    Returns:
        tuple: Tuple containing eccentricity (e), semi-major axis (major),
               semi-minor axis (minor), apogee (rmax), and perigee (rmin).
    """
    r_mag = np.linalg.norm(r, axis=1)

    rmax = np.max(np.abs(r_mag))  # Calculate the apogee (km)
    rmin = np.min(np.abs(r_mag))  # Calculate the perigee (km)

    major = (rmax + rmin) / 2  # Calculate semi-major axis
    minor = np.sqrt(rmax * rmin)  # Calculate the semi-minor axis

    e = np.sqrt(1 - (minor ** 2 / major ** 2))  # Eccentricity equation based on Keplarian Law's

    return e, major, minor, rmax, rmin


# Runge-Kutta 4 method
def RK4(t, r, v, tStep, obj, ro, vo):
    """
    Perform one iteration of the Runge-Kutta 4 method for numerical integration.

    Parameters:
        t (float): Current time.
        r (array-like): Array of position vectors.
        v (array-like): Array of velocity vectors.
        tStep (float): Time step size.
        obj (str): Object identifier ('earth' or 'moon').
        ro (array-like): Reference position vector.
        vo (array-like): Reference velocity vector.

    Returns:
        array-like: Array containing the updated position vector (r_new) and
                    velocity vector (v_new).
    """
    k11 = velocity(t, r, v, obj, ro, vo)
    k21 = acceleration(t, r, v, obj, ro, vo)

    k12 = velocity(t + 0.5 * tStep, r + 0.5 * tStep * k11, v + 0.5 * tStep * k21, obj, ro, vo)
    k22 = acceleration(t + 0.5 * tStep, r + 0.5 * tStep * k11, v + 0.5 * tStep * k21, obj, ro, vo)

    k13 = velocity(t + 0.5 * tStep, r + 0.5 * tStep * k12, v + 0.5 * tStep * k22, obj, ro, vo)
    k23 = acceleration(t + 0.5 * tStep, r + 0.5 * tStep * k12, v + 0.5 * tStep * k22, obj, ro, vo)

    k14 = velocity(t + tStep, r + tStep * k13, v + tStep * k23, obj, ro, vo)
    k24 = acceleration(t + tStep, r + tStep * k13, v + tStep * k23, obj, ro, vo)

    r_new = r + tStep * (k11 + 2. * k12 + 2. * k13 + k14) / 6.
    v_new = v + tStep * (k21 + 2. * k22 + 2. * k23 + k24) / 6.

    z = np.zeros([2, 2])
    z = [r_new, v_new]
    return z


# Verlet/Leapfrog method
def verlet(t, r, v, tStep, obj, ro, vo):
    """
    Perform one iteration of the Verlet/Leapfrog method for numerical integration.

    Parameters:
        t (float): Current time.
        r (array-like): Array of position vectors.
        v (array-like): Array of velocity vectors.
        tStep (float): Time step size.
        obj (str): Object identifier ('earth' or 'moon').
        ro (array-like): Reference position vector.
        vo (array-like): Reference velocity vector.

    Returns:
        array-like: Array containing the updated position vector (r_new) and
                    velocity vector (v_new).
    """
    a = acceleration(t, r, v, obj, ro, vo)

    r_new = r + v * tStep + 0.5 * a * tStep ** 2
    a_new = acceleration(t, r_new, v, obj, ro, vo)

    v_new = v + 0.5 * (a + a_new) * tStep

    z = np.zeros([2, 2])
    z = [r_new, v_new]
    return z

    return (line1, line2)


# Start by importing all required modules for the simulation.
import math
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm  # Provide a progress bar for the process
from mpl_toolkits.mplot3d import Axes3D  # Need in three dimensions to see change in inclination for Kozai mechanism

Mm = 7.346e22  # Mass of Moon (kg)
Me = 5.972e24  # Mass of Earth (kg)
Ms = 1.988e30  # Mass of Sun (kg)

# Time interval
ti = 0
tf = eval(input("Please provide the orbital time in years: "))  # Prompt user for orbital period
tf = tf * 3.1536e7  # Seconds in a year
nStep = 10000  # Itertions
t = np.linspace(ti, tf, nStep)
tStep = t[2] - t[1]  # Uniform spacing

# Numerical Method selection
NumericalMethod = eval(input("Choose numerical method:\n1) Runge-Kutta\n2) Verlet\nEnter your choice: "))

# Initialization
Ke = np.zeros(nStep)  # Kinetic energy
Pe = np.zeros(nStep)  # Potential energy
Te = np.zeros(nStep)  # Total Energy
Lz = np.zeros(nStep)  # Angular Momentum needed for Kozai mechanism

# Define constants and reference quantities
G = 6.6743e-20  # Gravitational Constant (km^3 kg^-1 yr^-2)

rm0 = [1.525e8, 0, 0]  # Average distance from the Moon to the Sun (Universe Today)
re0 = [1.5e8, 0, 0]  # Average distance from the Earth to the Sun (1 AU)
rm_re0 = 1.503e8 - 1.5e8  # Average distance from Moon to Earth

rm = np.zeros([nStep, 3])  # position vector of Moon (3D)
re = np.zeros([nStep, 3])  # position vector of Earth (3D)

vm = np.zeros([nStep, 3])  # velocity vector of Moon (3D)
ve = np.zeros([nStep, 3])  # velocity vector of Earth (3D)

vm_mag = np.sqrt(Ms * G / rm0[0])  # Derived from Eqn. 1 remembering centripetal acceleration = v^2/r
ve_mag = np.sqrt(Ms * G / re0[0])  # Derived from Eqn. 1 remembering centripetal acceleration = v^2/r

# Inclination angle
inclination = eval(input("Set in inclination angle: "))  # Prompt user for inclination angle of Moon
inclination = np.radians(inclination)

vm0 = [0, vm_mag * np.cos(inclination), vm_mag * np.sin(inclination)]
ve0 = [0, ve_mag, 0]

# Initialize arrays with initial values
t[0] = ti
rm[0, :] = rm0
re[0, :] = re0
vm[0, :] = vm0
ve[0, :] = ve0
Ke[0] = KineticEnergy(ve[0, :])
Pe[0] = PotentialEnergy(re[0, :])
Te[0] = Ke[0] + Pe[0]

# Iterate over time interval
for iStep in tqdm(range(0, nStep - 1)):

    # Calculate new position and velocity using the desired method

    if NumericalMethod == 1:  # RK4
        [re[iStep + 1, :], ve[iStep + 1, :]] = RK4(t[iStep], re[iStep, :], ve[iStep, :], tStep, 'earth', rm[iStep, :],
                                                   vm[iStep, :])
        [rm[iStep + 1, :], vm[iStep + 1, :]] = RK4(t[iStep], rm[iStep, :], vm[iStep, :], tStep, 'moon', re[iStep, :],
                                                   ve[iStep, :])

        Ke[iStep + 1] = KineticEnergy(vm[iStep + 1, :])
        Pe[iStep + 1] = PotentialEnergy(rm[iStep + 1, :])
        Lz[iStep + 1] = AngularMomentum(rm[iStep + 1, :], vm[iStep + 1, :], inclination)

    if NumericalMethod == 2:  # Verlet
        [re[iStep + 1, :], ve[iStep + 1, :]] = verlet(t[iStep], re[iStep, :], ve[iStep, :], tStep, 'earth',
                                                      rm[iStep, :], vm[iStep, :])
        [rm[iStep + 1, :], vm[iStep + 1, :]] = verlet(t[iStep], rm[iStep, :], vm[iStep, :], tStep, 'moon', re[iStep, :],
                                                      ve[iStep, :])

        Ke[iStep + 1] = KineticEnergy(vm[iStep + 1, :])
        Pe[iStep + 1] = PotentialEnergy(rm[iStep + 1, :])
        Lz[iStep + 1] = AngularMomentum(rm[iStep + 1, :], vm[iStep + 1, :], inclination)

    # Calculate total energy of the system
    Te[iStep + 1] = Ke[iStep + 1] + Pe[iStep + 1]

# Calculate Moon's eccentricity at given inclination angle
Ec, major, minor, apogee, perigee = Eccentricity(rm)
print(f"Eccentricity of the Moon's orbit at inclination angle {np.degrees(inclination)}: {Ec}")

# Plotting
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(0, 0, 0, s=200, c='orange')
end_xe = re[-1, 0]
end_ye = re[-1, 1]
end_ze = re[-1, 2]
end_xm = rm[-1, 0]
end_ym = rm[-1, 1]
end_zm = rm[-1, 2]
ax.plot3D(re[:, 0], re[:, 1], re[:, 2], c='blue', linewidth=1)
ax.plot3D(rm[:, 0], rm[:, 1], rm[:, 2], c='grey', linewidth=1)
ax.scatter(end_xe, end_ye, end_ze, s=40, c='blue')
ax.scatter(end_xm, end_ym, end_zm, s=10, c='grey')

limits = 2.5e8
ax.set_xlim([-limits, limits])  # Set x-axis limits
ax.set_ylim([-limits, limits])  # Set y-axis limits
ax.set_zlim([-limits, limits])  # Set z-axis limits

ax.set_xlabel(r'$x$ position ($1.0 \times 10^{8}$ km)', fontsize=12, labelpad=5)
ax.set_ylabel(r'$y$ position ($1.0 \times 10^{8}$ km)', fontsize=12, labelpad=5)
ax.set_zlabel(r'$z$ position ($1.0 \times 10^{8}$ km)', fontsize=12, labelpad=5)
ax.set_title(f'Orbital Motion for {tf / 3.1536e7} years, $\Theta$ = {np.degrees(inclination)}', fontsize=14)

# Set axis tick labels font size
ax.tick_params(labelsize=10)

# Rotate the view angle
prompt = input('Custom viewing angle? [y/n]').lower()
if prompt == 'y':
    azimuth_angle = eval(input("Set azimuth angle in degrees: "))
    elevation_angle = eval(input("Set elevation angle in degrees: "))
    ax.view_init(elevation_angle, azimuth_angle)
else:
    print('Set to standard angle: ')

plt.show()

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

# Energy plot
ax1 = axes[0]
ax1.plot(t, Ke, ls='-.', label='Kinetic')
ax1.plot(t, Pe, ls='--', label='Potential')
ax1.plot(t, Te, ls='-', label='Total')
ax1.set_xlabel('Time (sec)')
ax1.set_ylabel(r'Energy ($M~km^3/sec^2$)')
ax1.legend()

# Angular Momentum plot
ax2 = axes[1]
ax2.plot(t[1:], Lz[1:], ls='-.', label='Angular Momentum')
ax2.set_xlabel('Time (sec)')
ax2.set_ylabel(r'Angular Momentum Z ($M~km^2/sec$)')
ax2.legend()

plt.tight_layout()
plt.show()

# List of eccentricity values recording from previous
# section for i >= 40 degrees in order to replicate Figure 6 in Ford et al. (2000).
ec_list = [0.01190264754056412, 0.00556461585008428, 0.010889524966638186, 0.02299698627281592, 0.03757053960036706,
           0.05526165317880892, 0.07359735960766663, 0.08987221468665313, 0.10082015775087794, 0.1031270595680474,
           0.0944413504266101, 0.0744204041313068, 0.04510152485799893, 0.026823678138207345, 0.05208650455592897,
           0.08993193398703138, 0.12228750586081652, 0.14555421300224397, 0.1555626943616907]
i_list = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90]

# Plotting the points
plt.scatter(i_list, ec_list, color='black')

# Plotting the line
plt.plot(i_list, ec_list, color='black')

# Adding labels and title with LaTeX formatting
plt.xlabel('i (deg)')
plt.ylabel('$e_{max}$')
plt.title('Maximum eccentricity of a single oscillation')

# Displaying the legend
plt.legend()

# Displaying the plot
plt.show()

# Attempt to do some curve fitting to analytical solution for a Sun-Earth system
from scipy.optimize import curve_fit


# Define the function to fit
def quadratic_func(angle, a, b, c):
    """
    Evaluate the quadratic function.

    Parameters:
        angle (array-like): Input values.
        a (float): Coefficient of the quadratic term.
        b (float): Coefficient of the linear term.
        c (float): Constant term.

    Returns:
        array-like: Output values evaluated by the quadratic function.
    """
    return a * angle ** 2 + b * angle + c


# Generate data
theta = np.linspace(0, 2 * np.pi, 10000)  # Range of theta from 0 to 2pi
l = 1.5e8  # 1 AU
e = 0.01671  # Eccentricity of Earth (NASA)
r = l / (1 + e * np.cos(theta))
re = np.linalg.norm(re, axis=1)
# Perform the curve fit
popt, pcov = curve_fit(quadratic_func, theta, re)

# Extract the optimized parameters
a_fit, b_fit, c_fit = popt

# Calculate the fitted curve
r_fit = quadratic_func(theta, a_fit, b_fit, c_fit)

# Calculate the root mean squared error (RMSE)
rmse = np.sqrt(np.mean((r - re) ** 2))

# Plot the data and the RMSE
fig, ax = plt.subplots()
ax.plot(r * np.sin(theta), r * np.cos(theta), 'b--')
ax.plot(re * np.sin(theta), re * np.cos(theta), 'r--')
ax.plot(r_fit * np.sin(theta), r_fit * np.cos(theta), 'g--')

# Add labels for the legend
ax.plot([], [], 'b--', label="Approx. Analytical Earth.")
ax.plot([], [], 'r--', label="Method Earth.")
ax.plot([], [], 'g--', label="Best Fit Earth.")

ax.legend()
ax.set_title(f"Eccentricity of Earth, RMSE: {rmse:.2f}")

plt.show()

