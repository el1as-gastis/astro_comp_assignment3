from astroquery.gaia import Gaia
import matplotlib.pyplot as plt
import numpy as np

query = """
SELECT *
FROM gaiadr3.gaia_source AS gaia
JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS xmatch USING (source_id)
JOIN gaiadr3.tmass_psc_xsc_join AS xjoin USING (clean_tmass_psc_xsc_oid)
JOIN gaiadr1.tmass_original_valid AS tmass ON
   xjoin.original_psc_source_id = tmass.designation
WHERE
gaia.phot_g_mean_mag < 14
AND 1=CONTAINS(
    POINT('ICRS', gaia.ra, gaia.dec),
    CIRCLE('ICRS', 132.825, 11.8, 1.0)
)
"""

# Launch the query
job = Gaia.launch_job(query)
# Get the results
result = job.get_results()

# Apply quality cut: filter out rows where 'ph_qual' is not 'AAA'
filtered_1 = result[result['ph_qual'] == 'AAA']

# Filter to identify stars with negative or non-positive parallaxes
filtered_2 = filtered_1[filtered_1['parallax'] > 0]


# INFO FOR PLOT 1
# Calculate distance in parsecs
distance = 1000 / filtered_2['parallax']  # parallax is in milliarcseconds
# Calculate absolute G magnitude
absolute_g_mag = filtered_2['phot_g_mean_mag'] - 5 * np.log10(distance) + 5
# Calculate Gaia BP-RP color (assuming BP and RP magnitudes are available in the data)
bp_rp_color = filtered_2['phot_bp_mean_mag'] - filtered_2['phot_rp_mean_mag']



# INFO FOR PLOT 2
# Calculate J-Ks color index
j_ks_color = filtered_2['j_m'] - filtered_2['ks_m']
# Ks-band magnitude (apparent magnitude)
ks_mag = filtered_2['ks_m']



# Create a figure with two panels
fig, axes = plt.subplots(1, 2, figsize=(10, 5))

# Panel 1: Color-Magnitude Diagram (CMD)
axes[0].scatter(bp_rp_color, absolute_g_mag, s=1, color='blue')
axes[0].invert_yaxis()  # Invert y-axis for absolute magnitude
axes[0].set_xlabel('BP - RP [mag]')
axes[0].set_ylabel('Absolute G [mag]')
axes[0].set_title('BP - RP vs Absolute G CMD of M67')

# Panel 2: 
axes[1].scatter(j_ks_color, ks_mag, s=1, color='red')
axes[1].invert_yaxis()  # Invert y-axis for absolute magnitude
axes[1].set_xlabel('J - $K_s$ [mag]')
axes[1].set_ylabel('Apparent $K_s$ [mag]')
axes[1].set_title('J - $K_s$ vs Apparent $K_s$ CMD of M67')

# Adjust layout and show plot
plt.tight_layout()
plt.savefig('/home/el1as/github/astro3/astro_comp_assignment3/figures/cmds M67.png', dpi=200)