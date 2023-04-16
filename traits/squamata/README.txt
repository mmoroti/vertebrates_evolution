The data used for the analyses is presented in the 'RData' files, named according to the analysis type they were used for (see details below and in the related article).
Additionally, full species-level data, as well as 'legend' and 'references' sheets are presented in the excel file "species_level_data.xlsx".

See the code for data manipulation, analysis, and plotting in the file "Squamate_viviparity_code.Rmd".
The above code also uses parts of the "Code for Calculating a Corrected SE using Moran's I in SEM" by Jarrett E.K. Byrnes (see link below):
  https://github.com/jebyrnes/spatial_correction_lavaan

The two 'csv' files are summaries of model estimates, used for creating plots via the R script above.


#### General legend for all the 'RData' files:

The environment of each 'RData' file includes eight datasets, where "dat" is the dataset including the grid cells used for the global analysis, and the other seven are subsets of it by biogeographical realms (see details below).
'dat': Global dataset
'dataAA': Australasia
'dataAT': Afrotropics
'dataIM': Indomalaya
'dataNA': Nearctic
'dataNT': Neotropics
'dataOC': Oceania
'dataPA': Palearctic


### Assemblage-level datasets:

Each dataframe (e.g., "dataAA") in the environment includes 14 variables, where each row is a grid cell, with variables representing characteristics within it:
"WorldID": ID of the grid cell
"ovi.richness": Richness of oviparous species
"vivi.richness": Richness of viviparous species
"proportion.viviparity": Proportion of viviparous species out of all species (viviparous+oviparous) within the grid cell
"long": longitude centroid
"lat": latitude centroid
"bio1": mean annual temperature (°C)
"bio4": temperature seasonality (standard deviation×100)
"bio15": precipitation seasonality (coefficient of variation×100)
"inter_temp_var": interannual variation in mean temperature (coefficient of variation×100)
"inter_prec_var": interannual variation in total precipitation (coefficient of variation×100)
"elevation": elevation (meters above sea level)
"realm_id": Biogeographical realm (according to the acronyms in the "general legend" above)
"land_type": the land type in the grid cell - 'mainland', 'large islands', or 'small islands'.

## The six datasets below differ in the distributional range type of the species included, based on elevational data for each species. Species' EOO (extent of occurrence) range maps were filtered by the elevational ranges if these were available [thus creating ESH (extent of suitable habitat) maps], or left as is (i.e., EOO) when elevational data was lacking.

# The next three files are based on both literature data on reproductive mode and on phylogenetically imputed values for species lacking such data.
'assemblage_level_data_imputed_obs_elev_only.RData': data based only on elevational data from the literature (i.e., excluding EOO ranges and ESH ranges that are based on imputed elevation values).
'assemblage_level_data_imputed_ESH_only.RData': data based both on elevational data from the literature and imputed elevation values (i.e., excluding EOO ranges only).
'assemblage_level_data_imputed.RData': including all ESH and EOO range maps.

# The next three files have the same description as above, but are based only on species with literature data on reproduction (i.e., where reproductive mode was not imputed).
'assemblage_level_data_non-imputed.RData'
'assemblage_level_data_non-imputed_ESH_only.RData'
'assemblage_level_data_non-imputed_obs_elev_only.RData'


### Species-level datasets:

The file is based on the "species_level_data.xlsx", and differs from it in the following way:
1) The dataframes in the RData file exclude 18 species with mixed reproductive mode.
2) The dataframes in the RData file exclude the following columns, which are present only in the "species_level_data.xlsx" file: ‘Major.clade’, ‘family’, ‘family.etc’, ‘source.reproductive.mode’, ‘phy_imp_support’, ‘elevMin_m’, ‘elevMin_m’, ‘elevation_verbatim’, ‘source.elev’, ‘lat’.
3) Similar to the assemblage-level RData file, the species-level RData file includes separate dataframes for each realm.
For the full legend for this dataset, see the excel file.

For further information, see the related article: https://doi.org/10.1111/geb.13598
