
###############################################################################
'''Given a companion's sky plane position at two epochs, this program produces
contour plots of the companion's possible orbital elements as functions of its
unknown line of sight (z, vz) coordinates at the first observation epoch. For
further details see Pearce, Wyatt & Kennedy 2015. This version of the code uses
dimensionless parameters. Hence the the companion's sky plane coordinates are
inputted as B and phi (see the paper for details), and the code considers the
dimensionless line of sight coordinates rho and nu (where rho = z/R, nu = vz/V,
where R and V are the companion's relative sky plane position and velocity
respectively. The code also outputs dimensionless orbital elements (i.e. the
semimajor axis is given in units of the binary sky separation). The user should
not have to edit anything beyond line 43. Note: the longitude of ascending node
is defined relative to the primary - companion separation vector at the first
epoch of observation. The default inputs are from Paul Kalas' 2013 paper on
Fomalhaut b, and should reproduce Figure 2 in Pearce, Wyatt & Kennedy 2015. The
code produces a Python plot and, if toggle_save_data = 1, outputs six .txt
files, one file for each orbital element. These files have the rho and nu
values along the edges, and a grid of the corresponding orbital element
values .'''
###############################################################################
# User inputs

# Companion's B and phi (defined in Pearce & Wyatt 2015). B < 1 for the orbit
# to potentially be bound. Phi in degrees, and defined as 0 <= phi <= 180 deg.
# Example: Fomalhaut b has B = 0.53 and phi = 21.
B = 0.53
phi = 21.

# Number of rho and nu points in the grid (integers)
N_rho = 100
N_nu = 100

# Save plot data as .txt files? (1 = yes, anything else = no)
toggle_save_data = 1

# Lists of contour levels for each orbital element. Angles in degrees
contour_levels = {'ap': [1.2,2,5], \
                  'e': [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95], \
                  'i': [20,40,60,75,80,90,100,105,110,120,140,160], \
                  'O': [0,15,25,90,180,195,205,270], \
                  'w': [0,45,90,135,180,225,270], \
                  'f': [0,40,80,120,150,180,225,270,305]}

###############################################################################
# Import libraries
import numpy as np                  # Numerical functions
import matplotlib.pyplot as plt     # Plotting functions
from matplotlib import gridspec     # Subplot gridding
###############################################################################
# Convert phi to radians
phi *= np.pi/180.

print
print 'B =', B
print 'phi =', phi/np.pi*180, 'deg'
print

###############################################################################
# Functions

def get_nu_max_rho(rho):
    '''For a given value of rho, where |rho| < rho_max, find maximum value of
    |nu| resulting in a bound orbit. Derived using v^2 r < 2 mu'''

    nu_max_rho = (B**-1*(1+rho**2)**-.5 - 1)**.5

    return nu_max_rho

#------------------------------------------------------------------------------
def get_rho_nu_data():
    '''Calculates the maximum values of rho and nu resulting in a bound orbit,
    and gets lists of tested rho and nu values'''

    # Calculate max values of |rho|, |nu| resulting in a bound orbit (Equation
    # 4)
    rho_max = (B**-2 - 1.)**.5
    nu_max = (B**-1 - 1.)**.5

    # Generate list of tested rho and nu values
    rho_list = np.linspace(-rho_max, rho_max, N_rho)
    nu_list = np.linspace(-nu_max, nu_max, N_nu)

    # Generate line separating bound and unbound orbits
    bound_rho_line = list(rho_list[:]) + list(reversed(rho_list[:]))
    bound_nu_line = []

    for rho_ind in range(2*N_rho):
        rho = bound_rho_line[rho_ind]

        nu_max_rho = get_nu_max_rho(rho)

        if rho_ind <= N_rho: bound_nu_line += [nu_max_rho]
        else: bound_nu_line += [-nu_max_rho]

    # Return dictionary of lists and values
    rho_nu_data = {'rho_list': rho_list, \
                 'nu_list': nu_list, \
                 'rho_list': rho_list, \
                 'nu_list': nu_list, \
                 'bound_rho_line': bound_rho_line, \
                 'bound_nu_line': bound_nu_line}

    return rho_nu_data

#------------------------------------------------------------------------------
def calc_elements(rho, nu):
    '''Derives dimensionless orbital elements from position and velocity using
    the method of Murray and Durmott 1999 (equations 2.126 - 2.139). Phi is in
    radians'''

    # -------------------------- Calculate elements ---------------------------

    # Dimensionless semimajor axis: ap = a / R
    ap = (.5*((1+rho**2)**-.5 - B*(1+nu**2))**-1)
    if ap < 0: ap = 1.e9    # If unbound, set a high to tidy subplot

    hp = (rho**2 - 2*rho*nu*np.cos(phi) + nu**2 + np.sin(phi)**2)**.5
    hxp = - rho*np.sin(phi)
    hyp = rho*np.cos(phi) - nu

    # Eccentricity
    e = (1 - 2*B*hp**2/ap)**.5

    # Inclination
    i = np.arccos(np.sin(phi) / hp)

    Si = np.sin(i)

    # Theta (w+f) and longitude of ascending node
    if i == 0: O, theta = 0., 0.        # Theta = 0 because r = x
    else:

        O = np.arctan2(hxp/(hp*Si), -hyp/(hp*Si))

        theta = np.arctan2(rho/(1+rho**2)**.5/Si, \
            (np.cos(O)*(1+rho**2)**.5)**-1*(1+rho*np.sin(O)*np.cos(i)/Si))

        while theta < 0: theta += 2*np.pi
        while theta >= 2*np.pi: theta -= 2*np.pi

    sgn = np.sign(np.cos(phi) + rho*nu)

    # True anomaly
    f = np.arctan2(ap*(1-e**2)/(hp*e)*(1+nu**2-hp**2/(1+rho**2))**.5*sgn, \
        (ap*(1-e**2)/(1+rho**2)**.5 - 1)/e)

    # Argument of pericentre
    w = theta - f

    # Convert angles to degrees, and define to lie between 0 and 360 deg:
    i *= 180./np.pi
    O *= 180./np.pi
    w *= 180./np.pi
    f *= 180./np.pi

    while O < 0: O += 360.
    while O >= 360: O -= 360.
    while f < 0: f += 360.
    while f >= 360: f -= 360.
    while w < 0: w += 360.
    while w >= 360: w -= 360.

    # Add elements to dictionary
    elements = {'ap': ap, \
                'e': e, \
                'i': i, \
                'O': O, \
                'w': w, \
                'f': f}

    return elements

#------------------------------------------------------------------------------
def get_element_grids(rho_nu_data):
    '''Cycles through rho and nu values. For each combination resulting in a
    bound orbit, calculates the corresponding orbital elements. Outputs grids
    of orbital elements for contour plotting.'''

    print 'Calculating elements...'

    # Unpack lists of rho and nu values
    rho_list = rho_nu_data['rho_list']
    nu_list = rho_nu_data['nu_list']

    # Initiate element grids, to be filled with orbital elements corresponding
    # to rho and nu values
    ap_mat = np.zeros((N_nu, N_rho))
    e_mat = np.zeros((N_nu, N_rho))
    i_mat = np.zeros((N_nu, N_rho))
    O_mat = np.zeros((N_nu, N_rho))
    w_mat = np.zeros((N_nu, N_rho))
    f_mat = np.zeros((N_nu, N_rho))

    # Cycle through rho and nu values, and derive orbital elements
    for rho_ind in range(len(rho_list)):
        rho = rho_list[rho_ind]

        for nu_ind in range(len(nu_list)):
            nu = nu_list[nu_ind]

            # Calculate elements corresponding to these rho, nu coordinates
            elements = calc_elements(rho, nu)

            # Add elements to matrices
            ap_mat[nu_ind][rho_ind] = elements['ap']
            e_mat[nu_ind][rho_ind] = elements['e']
            i_mat[nu_ind][rho_ind] = elements['i']
            O_mat[nu_ind][rho_ind] = elements['O']
            w_mat[nu_ind][rho_ind] = elements['w']
            f_mat[nu_ind][rho_ind] = elements['f']

    # Output matrices
    element_matrices = {'ap': ap_mat, \
                         'e': e_mat, \
                         'i': i_mat, \
                         'O': O_mat, \
                         'w': w_mat, \
                         'f': f_mat}

    return element_matrices

#------------------------------------------------------------------------------
def save_data(rho_nu_data, element_matrices):
    '''Saves the element grids as .txt files, with the rho and nu values along
    the grid edges. Also save the bound / unbound divide line.'''

    print 'Saving data...'

    # Unpack v and nu lists and bounding lines
    rho_list = rho_nu_data['rho_list']
    nu_list = rho_nu_data['nu_list']
    bound_rho_line = rho_nu_data['bound_rho_line']
    bound_nu_line = rho_nu_data['bound_nu_line']

    # Cycle through all six elements, outputting files for each
    for elmnt_str in ['ap','e','i','O','w','f']:

        # Get element grid
        mat = element_matrices[elmnt_str]

		# Flip grid in up / down direction (so y axis is positive at top and
		# negative at bottom of output file)
        flipped_mat = np.flipud(mat)

        # Define output file
        grid_file = file('%s_grid.txt' % elmnt_str, 'w')

        # Write out rho values along top of grid
        line = ''
        for rho in rho_list: line += ' %s' % rho
        line += '\n'
        grid_file.write(line)

        # Cycle through rho and nu values, adding the corresponding elements to
        # the grid. Also add nu values down the left hand side of the grid
        for nu_ind in range(len(nu_list)):

			# Read out nu_list in reverse order (so y axis is positive at top
			# and negative at bottom of output file)
            nu = nu_list[-1-nu_ind]

            line = '%s' % nu

            for rho_ind in range(len(rho_list)):
                rho = rho_list[rho_ind]        
            
                # Get element matrix entry and add to line
                elmnt = flipped_mat[nu_ind][rho_ind]
                line += ' %s' % elmnt

            line += '\n'
            grid_file.write(line)

        grid_file.close()

    # Write bound/unbound divide file
    bound_line_file = file('bound_line.txt', 'w')
    bound_line_file.write('rho    nu\n')

    for ind in range(2*N_rho):
        rho, nu = bound_rho_line[ind], bound_nu_line[ind]
        bound_line_file.write('%s %s\n' % (rho, nu))

    bound_line_file.close()

#------------------------------------------------------------------------------
def make_individual_cntr_plt(fig, gs, elmnt_str, rho_nu_data, element_matrices,
        subplot_pars, contour_levels):
    '''For the given orbital element (elmnt_str), construct the corresponding
    subplot of element contours vs rho and nu'''

    # Unpack all required values
    rho_list = rho_nu_data['rho_list']
    nu_list = rho_nu_data['nu_list']
    bound_rho_line = rho_nu_data['bound_rho_line']
    bound_nu_line = rho_nu_data['bound_nu_line']
    mat = element_matrices[elmnt_str]
    subplot_title = subplot_pars[elmnt_str]['title']
    subplot_number = subplot_pars[elmnt_str]['number']

    # Make subplot
    ax = fig.add_subplot(gs[subplot_number])

    # Plot and label contours, with appropriate numbers of decimal places
    CS = plt.contour(rho_list, nu_list, mat, contour_levels[elmnt_str])
    if elmnt_str in ['ap', 'e']:
        ax.clabel(CS, inline=1, fontsize=10, fmt = '%0.1f')
    else: ax.clabel(CS, inline=1, fontsize=10, fmt = '%1.0f')

    # Add bound / unbound line
    ax.plot(bound_rho_line, bound_nu_line, 'k--')

    # Set subplot axis labels, if necessary
    if subplot_number in [0,1,2]:
        ax.xaxis.tick_top()
        ax.xaxis.set_ticks_position('both')
        ax.xaxis.set_label_position('top')

    if subplot_number in [2,5]:
        ax.yaxis.tick_right()
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_label_position('right')
        ax.set_ylabel(r'$\dot{z} / V$', labelpad=20, rotation=270,\
            fontsize = 16, fontname="Times New Roman")

    if subplot_number in [0,3]: ax.set_ylabel(r'$\dot{z} / V$',\
        fontsize = 16, fontname="Times New Roman")

    if subplot_number in [1,4]:
        ax.tick_params(labelleft='off')

    ax.set_xlabel(r'$z$ / R', fontsize = 16, fontname="Times New Roman")

    # Add subplot title
    ax.text(.05,.95, subplot_title, transform=ax.transAxes, ha='left', \
        va='top', fontsize = 20, fontname="Times New Roman", \
        bbox=dict(facecolor='white', edgecolor='white', pad=1), zorder=4)

#------------------------------------------------------------------------------
def make_contour_plots(rho_nu_data, element_matrices):
    '''Plots contours for all six elements as functions of rho and nu.'''

    print 'Making contour plots...'

    # Initialise figure
    fig = plt.figure()

    # Set subplot ratios and white space widths
    gs = gridspec.GridSpec(2, 3)
    gs.update(hspace=0., wspace = 0.)

    subplot_pars = {'ap': {'title': '$a / R$', 'number': 0}, \
                    'e': {'title': '$e$', 'number': 1}, \
                    'i': {'title': '$i$ / $^\circ$', 'number': 2}, \
                    'O': {'title': '$\Omega$ / $^\circ$', 'number': 3}, \
                    'w': {'title': '$\omega$ / $^\circ$', 'number': 4}, \
                    'f': {'title': '$f$ / $^\circ$', 'number': 5}}    

    elmnt_strs = ['ap', 'e', 'i', 'O', 'w', 'f']

    # Generate the contour plot for each orbital element
    for elmnt_str in elmnt_strs:
        make_individual_cntr_plt(fig, gs, elmnt_str, rho_nu_data, \
            element_matrices, subplot_pars, contour_levels)

    # Display figure
    plt.show()

###############################################################################
# Program

# Define tested region of rho, nu space
rho_nu_data = get_rho_nu_data()

# Cycle through rho, nu values, and derive orbital elements at each set of
# values
element_matrices = get_element_grids(rho_nu_data)

# Save data if necessary
if toggle_save_data == 1: save_data(rho_nu_data, element_matrices)

# Make contour plots
make_contour_plots(rho_nu_data, element_matrices)

print 'Program complete'
print
###############################################################################

