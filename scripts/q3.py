import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy.stats import binned_statistic_2d
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines

# File path for fits file
fits_file_path = '/home/el1as/github/astro3/astro_comp_assignment3/data/nihao_uhd_simulation_g8.26e11_xyz_positions_and_oxygen_ao.fits'

# Open the fits file to discover that useful data is in Second HDU
hdul = fits.open(fits_file_path)
data = hdul[1].data

# Access columns (x, y, z, A(O))
x = data['x']
y = data['y']
z = data['z']
A_o = data['A_O']

# Define galactocentric radii
R_gal = np.sqrt(x**2 + y**2 + z**2)

# Rescaled R_Gal to its standard deviation (Hyperparameter 1)
mean_rgal = np.mean(R_gal)
std_rgal = np.std(R_gal)
R_gal_scaled = (R_gal - mean_rgal) / std_rgal

# Initial guess (Hyperparameter 2)
initial_slope_guess = -0.1
initial_intercept_guess = 9

# Define a linear function for curve_fit 
def linear_func(x, m, c):
    return m * x + c

# Perform the fit using curve_fit with an initial guess on the scaled R_gal
popt, pcov = curve_fit(linear_func, R_gal_scaled, A_o, p0=[initial_slope_guess, initial_intercept_guess])
scaled_slope, scaled_intercept = popt
slope_uncertainty, intercept_uncertainty = np.sqrt(np.diag(pcov))

# Transform slope and intercept back to original scale
slope = scaled_slope / std_rgal
intercept = scaled_intercept - (scaled_slope * mean_rgal / std_rgal)

# Calculate fitted values and residuals using the unscaled slope and intercept
A_o_fit = linear_func(R_gal, slope, intercept)
residuals = A_o - A_o_fit

# Print fit results with uncertainties
print(f"Slope: {slope:.6f} ± {slope_uncertainty:.6f}")
print(f"Intercept: {intercept:.6f} ± {intercept_uncertainty:.6f}")

# Calculate goodness-of-fit metrics
ss_res = np.sum((A_o - A_o_fit) ** 2)
ss_tot = np.sum((A_o - np.mean(A_o)) ** 2)
r_squared = 1 - (ss_res / ss_tot)
rms = np.sqrt(np.mean(residuals**2))

print(f"R-squared: {r_squared:.6f}")
print(f"RMS: {rms:.6f}")

# Define the number of bins 
num_bins = 5  
bins = np.linspace(min(R_gal), max(R_gal), num_bins + 1)

# Arrays to store RMS and bin centers
rms_per_bin = []
bin_centers = (bins[:-1] + bins[1:]) / 2

# Loop over each bin
for i in range(num_bins):
    # Get the data points within the current bin
    bin_mask = (R_gal >= bins[i]) & (R_gal < bins[i + 1])
    R_gal_bin = R_gal[bin_mask]
    A_o_bin = A_o[bin_mask]
    A_o_fit_bin = A_o_fit[bin_mask]
    
    # Calculate RMS for the current bin
    residuals_bin = A_o_bin - A_o_fit_bin
    rms_bin = np.sqrt(np.mean(residuals_bin ** 2))
    rms_per_bin.append(rms_bin)

# Print RMS and R^2 values for each bin
for i in range(num_bins):
    print(f"R_range: {bins[i]} to {bins[i + 1]}")
    print(f"  RMS = {rms_per_bin[i]:.6f}")
    

# Plot
fig, axes = plt.subplots(1, 2, figsize=(10, 5))

# Panel 1: Logarithmic density plot of R_gal vs A_o with linear fit
im0 = axes[0].hexbin(R_gal, A_o, gridsize=80, cmap='magma', bins='log')
fit_line, = axes[0].plot(R_gal, A_o_fit, color='BLACK', label='Linear Fit', linewidth=1)  # Assign fit line to a variable

# Manually create a legend entry for the hexbin plot
hexbin_legend = mlines.Line2D([], [], color='orange', marker='h', markersize=10, linestyle='None', label='Data')

# Add colorbar
fig.colorbar(im0, ax=axes[0])

# Set axis labels and title
axes[0].set_xlabel('R$Gal$ [kpc]')
axes[0].set_ylabel('A(O)')
axes[0].set_title('Logarithmic density plot of R$Gal$ vs. A(O)')

# Combine both the custom hexbin legend and the linear fit legend
axes[0].legend(handles=[hexbin_legend, fit_line])


# Panel 2: Scatter plot of residuals
axes[1].scatter(R_gal, residuals, s=.1, color='red')
axes[1].set_xlabel('R$Gal$ [kpc]')
axes[1].set_ylabel(r'$\Delta A(O)$')
axes[1].axhline(0, color='black', linestyle='--')  # Add a horizontal line at zero
axes[1].set_title('Residuals of the linear fit function')

plt.tight_layout()
fig.savefig('/home/el1as/github/astro3/astro_comp_assignment3/figures/rg_vs_ao_fit_and_residuals.png')
plt.close()




































# Define the bins for the x and y axes
bins_x = np.linspace(min(x), max(x), 100)  
bins_y = np.linspace(min(y), max(y), 100)  

# Calculate the 2D histograms for the median simulated A(O), fitted A(O), and residuals
stat_simulated, _, _, _ = binned_statistic_2d(x, y, A_o, statistic='median', bins=[bins_x, bins_y])
stat_fitted, _, _, _ = binned_statistic_2d(x, y, A_o_fit, statistic='median', bins=[bins_x, bins_y])
stat_residuals, _, _, _ = binned_statistic_2d(x, y, residuals, statistic='median', bins=[bins_x, bins_y])

# Plot the 3-panel figure
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel (a): Simulated A(O)
im1 = axes[0].imshow(stat_simulated.T, origin='lower', cmap='magma', extent=[min(x), max(x), min(y), max(y)])
axes[0].set_xlabel('x')
axes[0].set_ylabel('y')
axes[0].set_title(' 2D-histogram of the median simulated A(O)')
fig.colorbar(im1, ax=axes[0])

# Panel (b): Fitted A(O)
im2 = axes[1].imshow(stat_fitted.T, origin='lower', cmap='magma', extent=[min(x), max(x), min(y), max(y)])
axes[1].set_xlabel('x')
axes[1].set_ylabel('y')
axes[1].set_title('2D-histogram of the median fitted A(O)')
fig.colorbar(im2, ax=axes[1])

# Panel (c): Residuals ∆A(O)
im3 = axes[2].imshow(stat_residuals.T, origin='lower', cmap='magma', extent=[min(x), max(x), min(y), max(y)])
axes[2].set_xlabel('x')
axes[2].set_ylabel('y')
axes[2].set_title(' 2D-histogram of the median residuals ∆A(O)')
fig.colorbar(im3, ax=axes[2])

# Adjust layout and save the figure
plt.tight_layout()
fig.savefig('/home/el1as/github/astro3/astro_comp_assignment3/figures/median_AO_3panel.png')
plt.close()