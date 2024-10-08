import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy.stats import binned_statistic_2d

# File path for fits file
fits_file_path = '/home/el1as/github/astro3/astro_comp_assignment3/data/nihao_uhd_simulation_g8.26e11_xyz_positions_and_oxygen_ao.fits'

# open the fits file path to discover that useful data is in Second HDU
hdul = fits.open(fits_file_path)
data = hdul[1].data

# note that data is in (x, y, z, A(O)), so we can access the columns by their names
x = data['x']
y = data['y']
z = data['z']
A_o = data['A_O']

# Define galactocentric radii
R_gal = np.sqrt(x**2 + y**2 + z**2)






# Define a linear function for curve_fit
def linear_func(x, m, c):
    return m * x + c

# Perform the fit using curve_fit
popt, pcov = curve_fit(linear_func, R_gal, A_o)
slope, intercept = popt
slope_uncertainty, intercept_uncertainty = np.sqrt(np.diag(pcov))

# Calculate fitted values and residuals
A_o_fit = linear_func(R_gal, slope, intercept)
residuals = A_o - A_o_fit

# Print fit results with uncertainties
print(f"Slope: {slope:.4f} ± {slope_uncertainty:.4f}")
print(f"Intercept: {intercept:.4f} ± {intercept_uncertainty:.4f}")






fig, axes = plt.subplots(1, 2, figsize=(10, 5))

# Panel 1: Logarithmic density plot of R_gal vs A_o with linear fit
axes[0].hexbin(R_gal, A_o, gridsize=50, cmap='magma', bins='log')
axes[0].plot(R_gal, A_o_fit, color='blue', label=f'Fit: slope={slope:.2f}', linewidth=2)  # Overplot linear fit
axes[0].set_xlabel('R$Gal$ ')
axes[0].set_ylabel('A(O)')
axes[0].set_title('Logarithmic density plot with linear fit')
axes[0].legend()

# Panel 2: Scatter plot of residuals
axes[1].scatter(R_gal, residuals, s=.1, color='red')
axes[1].set_xlabel('R$Gal$ ')
axes[1].set_ylabel(r'$\Delta A(O)$ (Residuals)')
axes[1].axhline(0, color='black', linestyle='--')  # Add a horizontal line at zero
axes[1].set_title('Residuals of the fit')

# # Panel 2: Logarithmic density plot of residuals
# axes[1].hexbin(R_gal, residuals, gridsize=50, cmap='magma', bins='log')
# axes[1].set_xlabel('R$Gal$ ')
# axes[1].set_ylabel(r'$\Delta A(O)$ (Residuals)')
# axes[1].axhline(0, color='black', linestyle='--')  # Add a horizontal line at zero
# axes[1].set_title('Residuals of the fit')

plt.tight_layout()
fig.savefig('/home/el1as/github/astro3/astro_comp_assignment3/figures/rg_vs_ao_fit_and_residuals.png')
plt.close()


# Define the bin edges for the x and y axes
bins_x = np.linspace(min(x), max(x), 80)  
bins_y = np.linspace(min(y), max(y), 80)  

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
axes[0].set_title('Median Simulated A(O)')
fig.colorbar(im1, ax=axes[0])

# Panel (b): Fitted A(O)
im2 = axes[1].imshow(stat_fitted.T, origin='lower', cmap='magma', extent=[min(x), max(x), min(y), max(y)])
axes[1].set_xlabel('x')
axes[1].set_ylabel('y')
axes[1].set_title('Median Fitted A(O)')
fig.colorbar(im2, ax=axes[1])

# Panel (c): Residuals ∆A(O)
im3 = axes[2].imshow(stat_residuals.T, origin='lower', cmap='magma', extent=[min(x), max(x), min(y), max(y)])
axes[2].set_xlabel('x')
axes[2].set_ylabel('y')
axes[2].set_title('Median Residuals ∆A(O)')
fig.colorbar(im3, ax=axes[2])

# Adjust layout and save the figure
plt.tight_layout()
fig.savefig('/home/el1as/github/astro3/astro_comp_assignment3/figures/median_AO_3panel.png')
plt.close()