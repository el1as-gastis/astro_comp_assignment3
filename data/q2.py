from astroquery.gaia import Gaia

query = """
SELECT gaia.source_id, gaia.ra AS g_ra, gaia.dec AS g_dec, gaia.phot_g_mean_mag, gaia.parallax, gaia.pmra, gaia.pmdec, tmass.*
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
print(result)


# Apply quality cut: filter out rows where 'ph_qual' is not 'AAA'
filtered_result = result[result['ph_qual'] == 'AAA']
print(filtered_result)

# Filter to identify stars with negative or non-positive parallaxes
positive_parallax = filtered_result[filtered_result['parallax'] > 0]
print(positive_parallax)