# Bushinsky-and-Cerovecki-SAMW-formation-properties-2022

Processing and plotting script associated with: "Subantarctic Mode Water Biogeochemical Formation Properties and Interannual Variability", Seth M. Bushinsky<sup>1</sup> and Ivana Cerovecki<sup>2</sup>

<sup>1</sup>Department of Oceanography, School of Ocean and Earth Science and Technology, University of Hawai'i at Manoa, Honolulu, HI  
<sup>2</sup>Scripps Institution of Oceanography, University of California San Diego, La Jolla, CA

Corresponding author: Seth Bushinsky (seth.bushinsky@hawaii.edu)

___

Processing script: Bushinsky_and_Cerovecki_SAMW_2022_R1.m

Data files needed (see manuscript for detailed citations):
- SO_calc_09-Mar-2022_w_calc_param_pco2_W14_wCalcs.mat - May 2021 SOCCOM Snapshot (doi.org/10.6075/J0T43SZG ) plus Drucker and Riser 2016 UW Oxygen dataset (https://argo.ucsd.edu/data/data-from-gdacs/ ). Secondary QC performed. Calculations of potential density, pot. temperature, MLD. Other calculated parameters are included but not used in this analysis. https://www.dropbox.com/s/5wasxkdf7zgvdyk/SO_calc_09-Mar-2022_w_calc_param_pco2_W14_wCalcs.mat?dl=0 
- Glodapv2.2020_Merged_Master_File.mat - https://www.glodap.info/index.php/merged-and-adjusted-data-product/ 
- SOCATv2021 Southern Oceans- https://www.socat.info/index.php/data-access/, read and processed using code (Read_SOCATv3_v2021.m) available from the same page.  SOCAT data quality flags ABCD and WOCE QC flag 2. 
- co2_GHGreference.901469012_surface.txt - https://gml.noaa.gov/ccgg/mbl/data.php. Surface data from 1979-01–01 to 2020-01-01.
- MLD_av_2005_OCT_2021_SouthOc.mat - Monthly 1x1 gridded MLDs from core Argo data; included in Data/ 
- MW_Vol_RG_Argo20[05-21]_All_Densities_updated.mat - monthly MW volumes binned by density on a 1x1 deg grid; included in Data/


____

Plotting script: Bushinsky_and_Cerovecki_SAMW_2022_figures_R1.m

Data files needed: 

- regional_boundaries_5zone.mat - SO Front boundaries as defined/used in Gray et al. 2018 and Bushinsky et al. 2019; included in Data/
- sam_nino_2005_2021.mat - Monthly SAM and ENSO indices; included in Data/
- BSOSE output for figures (included in Data/): 
  - MLD_AnAv_SAMW_range8_9_SEPac_200_m_025_ST_iter135_45_64_246_290_se_pacific.mat
  - BSOSE_iter135_[DIC,NO3,O2]_2013_2019_SouthOc_MLD_150m_AugSep_BSOSE_mask_v2.mat

Colormaps: 
- ColorBrewer - https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps
- Other colormaps from: https://www.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps . Ander Biguri (2022). Perceptually uniform colormaps (https://www.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps), MATLAB Central File Exchange. 

____

Other functions called (in Subfunctions/):
- ncep_pressure_matching.m - loads annual files of 4x daily surface pressure from NCEP Reanalysis: https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/surface/ , finds nearest lat/lon and interpolates to time.
- CO2SYSSOCCOM.m - CO2sys version modified by Nancy Williams for SOCCOM-relevant constants
- GGo2_units.m - calculates oxygen saturation concentration from Garcia and Gordon 1992
- SW properties toolbox - https://researchdata.edu.au/csiro-marine-research-library-2006/692093 
- fit_cosine.m, fit_harmonics.m, correlate.m- from from Dr. Kathie Kelly's Data Analysis course at UW ~ 2014
- pCO2_from_fCO2.m - for calculating pCO2 from fCO2 and temperature
- mld_dbm.m - De Boyer Montegut et al. 2004 method for calculating MLD, written by Dr. Noel Pelland
- importfile_surf_CO2_2019_02_22.m - for reading in NOAA ESRL surface CO2 files\
- ph2osat_smb.m - calculate water vapor pressure as a function of T and S

Matlab toolboxes required:
- Matlab v9.7 (R2019b)
- Statistics and Machine Learning v11.6
- Mapping Toolbox v4.9
- Curve Fitting Toolbox v3.5.10
- Image Processing Toolbox
