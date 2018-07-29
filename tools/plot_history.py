# Input: history

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

def main(name, type, plot, x, y, z):
    # Load data from CSV
    dat = np.genfromtxt(name, delimiter='\t',skip_header=1)

    X_dat = dat[:,x]
    Y_dat = dat[:,y]
    Z_dat = dat[:,z]

    h = 3
    h_dat = dat[:,h]

    # Convert from pandas dataframes to numpy arrays
    X, Y, Z, = np.array([]), np.array([]), np.array([])
    H = np.array([])
    cap_z = 400.0
    for i in range(len(X_dat)):
        #if Z_dat[i] < cap_z and X_dat[i] > 0 and Y_dat[i] > 0 and X_dat[i] < 20 and dat[:,2][i] == 200:
        #if True:       # You can add restrictions on the data here
        if (not np.isnan(Z_dat[i]) and Z_dat[i] < cap_z) and dat[:,2][i] == 150:# and dat[:,3][i] == 0.9: 
            X = np.append(X,X_dat[i])
            Y = np.append(Y,Y_dat[i])
            Z = np.append(Z,Z_dat[i])
            H = np.append(H,h_dat[i])
    
    # create x-y points to be used in heatmap
    xi = np.linspace(X.min(),X.max(),1000)
    yi = np.linspace(Y.min(),Y.max(),1000)
    # Z is a matrix of x-y values
    zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method=type)
    

    if plot == 0:
        # Create a heatmap
        CS = plt.contourf(xi, yi, zi, 300, cmap=plt.cm.rainbow)
    elif plot == 1:
        #Create a scatter plot
        CS = plt.scatter(X, Y, c=Z, s=50, edgecolor='')
    elif plot == 2:
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	p3d = ax.scatter(X, Y, H, c=Z)

	plt.colorbar(p3d)

    # Plot the minimum value ontop of the heatmap 
    min=10000
    with open(name, 'r') as r:
            for num, line in enumerate(r,1):
                list_line = line.split()
                if "Warning" not in line and 'error' not in line:
                    if list_line[-1] != 'None' and float(list_line[-1]) < min:
                        min = float(list_line[-1])
                        min_line = line
    
    min_r = min_line.split()[0]
    min_k = min_line.split()[1]
    
    print min_line
    
    #plt.plot([min_r], [min_k], 'go')
    
    #plt.colorbar()  
    plt.show()
    

if __name__ == "__main__":

    helpstring = """ Input: fort.99, name Output: *png """

    if len(sys.argv) != 3:
        print helpstring
        sys.exit()
  
    name = sys.argv[1]
    type = sys.argv[2]

    # Plot = 0 heatmap, 1 scatterplot
    plot = 2

    # The x, y, and z columns to plot (starting counting at 0)
    x = 0
    y = 1
    z = 5

    main(name, type, plot, x, y, z)


