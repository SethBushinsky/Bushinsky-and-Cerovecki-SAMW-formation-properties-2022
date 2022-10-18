% Processing script associated with: "Subantarctic Mode Water
% Biogeochemical Formation Properties and Interannual Variability", Seth M. Bushinsky1 and Ivana Cerove?ki2
%
% 1Department of Oceanography, School of Ocean and Earth Science and Technology, University of Hawai?i at M?noa, Honolulu, HI
% 2Scripps Institution of Oceanography, University of California San Diego, La Jolla, CA
%
% Corresponding author: Seth Bushinsky (seth.bushinsky@hawaii.edu)

%% Load data:

% Load float data
float_data_dir = [data_dir 'ARGO_O2_Floats/Global/SOCCOM/2021_05_05_Snapshot_LoRes_LIAR/'];
float_data_filename = 'SO_calc_09-Mar-2022_w_calc_param_pco2_W14_wCalcs.mat'; % 0.03 MLD cutoff - % updated with fully QC'd dataset
load([float_data_dir float_data_filename])

cerovecki_dir = [home_dir 'Work/Projects/2020_02 SO SAMW Variability_Cerovecki/'];


% set geographic limits for Pacific and Indian SAMW formation regions
geo_lims.pac_lims = [170 290];
geo_lims.ind_lims = [68 170];
geo_lims.lat_lims_pac = [-64 -45];
geo_lims.lat_lims_ind = [-55 -30];

month_lims_for_mean = [8 9];
MLD_cutoff = 200;% only consider float data with MLs greater than this depth
clear SOCCOM_float_directory float_data_dir float_data_filename

% Load Glodap data

gdap = load([data_dir 'Data_Products/GLODAP/GLODAPv2.2020_Merged_Master_File.mat']);
disp('loading glodap')

% filter glodap to Southern Ocean, calculate O2 sauration concentration using Garcia and Gordon and
% delta O2
disp('starting glodap calculations')
clear gdap_SO
gdap.PDENS = sw_pden(gdap.G2salinity,gdap.G2temperature,gdap.G2pressure, 0);
gdap.GMT_Matlab = datenum(gdap.G2year, gdap.G2month, gdap.G2day, gdap.G2hour, gdap.G2minute, zeros(size(gdap.G2minute)));

gdap_fields_no_flags = {'G2longitude', 'G2latitude', 'GMT_Matlab', 'G2year', 'G2month', 'G2cruise', 'G2station', 'G2pressure', 'G2temperature', 'G2theta', 'PDENS'};
gdap_fields_w_flags = {'G2salinity','G2oxygen', 'G2nitrate', 'G2tco2', 'G2talk', 'G2fco2'};

gdap.G2longitude(gdap.G2longitude<0) = gdap.G2longitude(gdap.G2longitude<0)+360;

SO_gdap_index = gdap.G2latitude<-30 & ~isnan(gdap.GMT_Matlab);

for g = 1:length(gdap_fields_no_flags)
    gdap_SO.(gdap_fields_no_flags{g}) = gdap.(gdap_fields_no_flags{g})(SO_gdap_index);
end

disp('processing glodap data');
% set data to nan if it doesn't have a flag of 2
for g = 1:length(gdap_fields_w_flags)
    temp_data = gdap.(gdap_fields_w_flags{g});
    temp_data(gdap.([gdap_fields_w_flags{g} 'f'])~=2) = nan;
    
    gdap_SO.(gdap_fields_w_flags{g}) =temp_data(SO_gdap_index);
end

gdap_SO.G2o2sat = GGo2_units(gdap_SO.G2temperature , gdap_SO.G2salinity , 'umol');
gdap_SO.G2deltaO2 = gdap_SO.G2oxygen - gdap_SO.G2o2sat;

clear gdap temp_T temp_T_0 temp_S temp_S_0 temp_P temp_P_0 s c g temp_data

% Calculate pCO2 from GLODAP tco2 and alk
[DATA,~,~]= CO2SYSSOCCOM(gdap_SO.G2talk, gdap_SO.G2tco2 , ...
    1,2,gdap_SO.G2salinity, gdap_SO.G2temperature, ...
    gdap_SO.G2temperature,...
    gdap_SO.G2pressure,gdap_SO.G2pressure,20,1.8,1,10,3);

gdap_SO.G2pCO2 = (DATA(:,19));

clear SO_gdap_index g DATA
% Load SOCAT dataset

% load socat
disp('loading socat')
% socat = load([home_dir 'Data/Data_Products/SOCAT/v2020/SOCATv2020_SouthernOceans_ABCD.mat']);
socat = load([data_dir 'Data_Products/SOCAT/v2021/SOCATv2021_SouthernOceans_ABCD_WOCE_2.mat']);

%% Calculate SLP climatology and calculate delta pCO2 for all datasets

disp('Calculating 10-yr SLP climatology')

% Assemble a list of lat/lon position for all float and shipboard data
% rounded to the nearest 2.5 deg. This will be used for the SLP adjustment
total_lat_lon_list = [];

for f = 1:length(SO_SNs)
    
    temp_date_vec = datevec(Argo.(SO_SNs{f}).GMT_Matlab);
    geo_index = (Argo.(SO_SNs{f}).Lat>geo_lims.lat_lims_pac(1) & Argo.(SO_SNs{f}).Lat<geo_lims.lat_lims_pac(2) & ...
        Argo.(SO_SNs{f}).Lon>=geo_lims.pac_lims(1) & Argo.(SO_SNs{f}).Lon<geo_lims.pac_lims(2)) | ...
        (Argo.(SO_SNs{f}).Lat>geo_lims.lat_lims_ind(1) & Argo.(SO_SNs{f}).Lat<geo_lims.lat_lims_ind(2) & ...
        Argo.(SO_SNs{f}).Lon>=geo_lims.ind_lims(1) & Argo.(SO_SNs{f}).Lon<geo_lims.ind_lims(2)) & ...
        temp_date_vec(:,2)>=month_lims_for_mean(1) & temp_date_vec(:,2)<=month_lims_for_mean(end);
    
    temp_lats = Argo.(SO_SNs{f}).Lat(geo_index);
    temp_lats_nearest_pt5 = round(temp_lats.*.4)./.4;
    
    temp_lons = Argo.(SO_SNs{f}).Lon(geo_index);
    temp_lons_nearest_pt5 = round(temp_lons.*.4)./.4;
    
    
    temp_lat_lon_list = unique([temp_lats_nearest_pt5 temp_lons_nearest_pt5], 'rows');
    
    total_lat_lon_list = [total_lat_lon_list ; temp_lat_lon_list];
    
    clear temp_lat_lon_list temp_lons temp_lats_nearest_pt5 temp_lons_nearest_pt5 geo_index temp_date_vec temp+lons temp_lats
end

% only save unique points
total_lat_lon_list = unique(total_lat_lon_list, 'rows');

temp_date_vec = datevec(gdap_SO.GMT_Matlab);

% find unique gdap lat/lon combos
geo_index = (gdap_SO.G2latitude>geo_lims.lat_lims_pac(1) & gdap_SO.G2latitude<geo_lims.lat_lims_pac(2) & ...
    gdap_SO.G2longitude>=geo_lims.pac_lims(1) & gdap_SO.G2longitude<geo_lims.pac_lims(2)) | ...
    (gdap_SO.G2latitude>geo_lims.lat_lims_ind(1) & gdap_SO.G2latitude<geo_lims.lat_lims_ind(2) & ...
    gdap_SO.G2longitude>=geo_lims.ind_lims(1) & gdap_SO.G2longitude<geo_lims.ind_lims(2)) & ...
    temp_date_vec(:,2)>=month_lims_for_mean(1) & temp_date_vec(:,2)<=month_lims_for_mean(end);

temp_lats = gdap_SO.G2latitude(geo_index);
temp_lats_nearest_pt5 = round(temp_lats.*.4)./.4;

temp_lons = gdap_SO.G2longitude(geo_index);
temp_lons_nearest_pt5 = round(temp_lons.*.4)./.4;

temp_lat_lon_list = unique([temp_lats_nearest_pt5 temp_lons_nearest_pt5], 'rows');

% concatenate w/ float list
total_lat_lon_list = [total_lat_lon_list ; temp_lat_lon_list];

clear temp_lat_lon_list temp_lons temp_lats_nearest_pt5 temp_lons_nearest_pt5 geo_index temp_date_vec temp+lons temp_lats


% find unique socat lat/lon combos
geo_index = (socat.latitude>geo_lims.lat_lims_pac(1) & socat.latitude<geo_lims.lat_lims_pac(2) & ...
    socat.longitude>=geo_lims.pac_lims(1) & socat.longitude<geo_lims.pac_lims(2)) | ...
    (socat.latitude>geo_lims.lat_lims_ind(1) & socat.latitude<geo_lims.lat_lims_ind(2) & ...
    socat.longitude>=geo_lims.ind_lims(1) & socat.longitude<geo_lims.ind_lims(2)) & ...
    socat.mon>=month_lims_for_mean(1) & socat.mon<=month_lims_for_mean(end);

temp_lats = socat.latitude(geo_index);
temp_lats_nearest_pt5 = round(temp_lats.*.4)./.4;

temp_lons = socat.longitude(geo_index);
temp_lons_nearest_pt5 = round(temp_lons.*.4)./.4;

temp_lat_lon_list = unique([temp_lats_nearest_pt5 temp_lons_nearest_pt5], 'rows');

% concatenate with float and ship list
total_lat_lon_list = [total_lat_lon_list ; temp_lat_lon_list];

total_lat_lon_list = unique(total_lat_lon_list, 'rows');

clear temp_lat_lon_list temp_lons temp_lats_nearest_pt5 temp_lons_nearest_pt5 geo_index temp_date_vec temp+lons temp_lats f
% At each location, calculate a mean seasonal cycle of pressure based on 10 years of NCEP reanalysis data.
% fit a harmonic mean to each location, use this for your pressure
% adjustment. This allows comparison of delta pCO2 between different years
% so that differences between ship and float delta pCO2 is actually due to
% differences in pCO2 and not passing storms, or other anomalies.

NCEP_pressures = [];
select_seasonal_press = NaN(length(total_lat_lon_list), 365);

press_year_day = 1:365;
press_dates = datenum(2011,1,press_year_day.*10);

for q = 1:length(total_lat_lon_list)
    [press_int, NCEP_pressures, ~, ~] = ncep_pressure_matching(press_dates, total_lat_lon_list(q,1), total_lat_lon_list(q,2), NCEP_pressures);
    % this seems to work for now.
    nharm=2;cutoff=.4;L=365.25;
    
    [var_amp2,var_phase2,~,var_offset2,~]= ...
        fit_harmonics(press_int', press_dates', nharm, L, cutoff);
    
    yt = var_offset2*ones(length(press_dates'),1);
    for j=1:nharm
        yt=yt+var_amp2(j)*cos(2*pi*j*press_dates'/L + var_phase2(j));
    end
    seasonal_var = yt;
    
    % % code for checking whether the harmonic fit was working properly
    % var_anomaly = press_int' - seasonal_var;
    % clf
    % subplot(2,1,1)
    % plot(press_dates, press_int, 'k-')
    % hold on
    % plot(press_dates, seasonal_var, 'r-', 'linewidth', 2)
    % subplot(2,1,2)
    % plot(press_dates, var_anomaly, 'k-')
    % title([num2str((nansum(var_anomaly.^2))^.5)])
    select_seasonal_press(q,:) = seasonal_var(1:365)';
end
clear yt j nharm var_amp2 var_phase2 var_frac2 var_offset2 var_residual2 press_int NCEP_pressures L seasonal_var cutoff

% calculate delta pCO2 for each dataset - start with floats
%  use the same calculation as for all three datasets and 10-yr avg pressure to make sure everything is uniform

noaa_co2_dir = [data_dir 'Data_Products/CO2_atmos/NOAA_ESRL/'];
noaa_co2_filename = 'co2_GHGreference.901469012_surface.txt';

noaa_co2 = importfile_surf_CO2_2019_02_22([noaa_co2_dir noaa_co2_filename]);

disp('Calculating delta pCO2 for Argo')
for f = 1:length(SO_SNs)
    if ~isfield(Argo.(SO_SNs{f}), 'pCO2_LIAR') % checks for presence of pCO2 estimates
        continue
    end
    disp(f)
    
    % calculate delta pCO2 -
    % loop through datapoints, matching lat/time
    Argo.(SO_SNs{f}).pC_dry_atmosphere = NaN(length(Argo.(SO_SNs{f}).GMT_Matlab),1);
    for d = 1:length(Argo.(SO_SNs{f}).GMT_Matlab)
        % find the closest matching noaa_co2 time. If time is greater than
        % 30 days away, sets pCO2 dry atmosphere (pC_dry_atmosphere) to nan
        date_index = min(abs(noaa_co2.Matlab_time - Argo.(SO_SNs{f}).GMT_Matlab(d)))==abs(noaa_co2.Matlab_time - Argo.(SO_SNs{f}).GMT_Matlab(d));
        if isnan(Argo.(SO_SNs{f}).Lat(d)) || Argo.(SO_SNs{f}).GMT_Matlab(d)> noaa_co2.Matlab_time(end)+30 % should check that this works! Seems to.
            Argo.(SO_SNs{f}).pC_dry_atmosphere(d) = nan;
        else
            lat_index = min(abs(noaa_co2.lat_co2 - Argo.(SO_SNs{f}).Lat(d)))==abs(noaa_co2.lat_co2 - Argo.(SO_SNs{f}).Lat(d));
            
            % save dry atmospheric CO2 into each individual float
            Argo.(SO_SNs{f}).pC_dry_atmosphere(d) = noaa_co2.co2_atmos(lat_index, date_index);
        end
    end
    
    
    % calculate water vapor pressure from T and S
    pH2O = ph2osat_smb(Argo.(SO_SNs{f}).Temp_C(:,1), Argo.(SO_SNs{f}).Sal(:,1));
    
    % find the closest SLP using nearest 2.5 deg lat/lon
    temp_lats = Argo.(SO_SNs{f}).Lat;
    temp_lats_nearest_pt5 = round(temp_lats.*.4)./.4;
    
    temp_lons = Argo.(SO_SNs{f}).Lon;
    temp_lons_nearest_pt5 = round(temp_lons.*.4)./.4;
    
    % finding SLP for each pt:
    temp_date_vec = datevec(Argo.(SO_SNs{f}).GMT_Matlab);
    Argo.(SO_SNs{f}).NCEP_SLP_CLIM = NaN(length(Argo.(SO_SNs{f}).GMT_Matlab),1);
    
    for d = 1:length(Argo.(SO_SNs{f}).GMT_Matlab)
        
        lat_lon_index = find(total_lat_lon_list(:,1)==temp_lats_nearest_pt5(d) & total_lat_lon_list(:,2)==temp_lons_nearest_pt5(d));
        if isempty(lat_lon_index)
            continue
        end
        
        lat_lon_index = lat_lon_index(1); % temp cludge b/c there are duplicates in the lat/lon list
        
        year_day = Argo.(SO_SNs{f}).GMT_Matlab(d) - datenum(temp_date_vec(d,1), 1, 0);
        
        day_index = min(abs(press_year_day - year_day))==abs(press_year_day - year_day);
        
        % saves out climatological SLP
        Argo.(SO_SNs{f}).NCEP_SLP_CLIM(d) = nanmean(select_seasonal_press(lat_lon_index, day_index));
    end
    
    Argo.(SO_SNs{f}).pC_moist_atmosphere = Argo.(SO_SNs{f}).pC_dry_atmosphere.*(Argo.(SO_SNs{f}).NCEP_SLP_CLIM./1013.25-pH2O); % uatm wet air
    
    % "unified" here refers to the fact that this is calculated using a
    % climatological SLP, not the insitu SLP so that I can compare delta
    % pCO2 between different years and not have short-term changes in SLP
    % changing things
    Argo.(SO_SNs{f}).delta_pCO2_unified = Argo.(SO_SNs{f}).pCO2_LIAR(:,1) - Argo.(SO_SNs{f}).pC_moist_atmosphere;
    
    clear lat_lon_index year_day day_index d temp_date_vec temp_lons temp_lons_nearest_pt5 temp_lats temp_lats_nearest_pt5 pH2O
end

clear f d date_index lat_index
% glodap delta pCO2 calculation

disp('calculating delta pCO2 for glodap')
% loop through datapoints, matching lat/time
gdap_SO.pC_dry_atmosphere = NaN(length(gdap_SO.GMT_Matlab),1);
decile_vec = round(length(gdap_SO.GMT_Matlab)./10).*[1:10];
decile_count = 1;
for d = 1:length(gdap_SO.GMT_Matlab)
    date_index = min(abs(noaa_co2.Matlab_time - gdap_SO.GMT_Matlab(d)))==abs(noaa_co2.Matlab_time - gdap_SO.GMT_Matlab(d));
    if isnan(gdap_SO.G2latitude(d)) || gdap_SO.GMT_Matlab(d)>noaa_co2.Matlab_time(end) + 30 % double check that this is working
        gdap_SO.pC_dry_atmosphere(d) = nan;
    else
        lat_index = min(abs(noaa_co2.lat_co2 - gdap_SO.G2latitude(d)))==abs(noaa_co2.lat_co2 - gdap_SO.G2latitude(d));
        
        gdap_SO.pC_dry_atmosphere(d) = nanmean(reshape(noaa_co2.co2_atmos(lat_index, date_index),[],1));
    end
    %     if d/length(gdap_SO.G2latitude)*10==round(d/length(gdap_SO.G2latitude)*10)
    %         disp([num2str((d/length(gdap_SO.G2latitude))*100,4) ' % done'])
    %     end
    if d>decile_vec(decile_count)
        disp([num2str((d/length(gdap_SO.G2latitude))*100,4) ' % done'])
        decile_count = decile_count+1;
    end
end

%     S_CO2 = CO2_Sol(Temperature,Salinity); % % mmol m-3 atm-1

pH2O = ph2osat_smb(gdap_SO.G2temperature, gdap_SO.G2salinity);

temp_date_vec = datevec(gdap_SO.GMT_Matlab);

% limit the search to a specific geographic range - can remove and just add
% a filter to check if lat_lon_index is empty - not sure which is quicker
geo_index = find((gdap_SO.G2latitude>geo_lims.lat_lims_pac(1) & gdap_SO.G2latitude<geo_lims.lat_lims_pac(2) & ...
    gdap_SO.G2longitude>=geo_lims.pac_lims(1) & gdap_SO.G2longitude<geo_lims.pac_lims(2)) | ...
    (gdap_SO.G2latitude>geo_lims.lat_lims_ind(1) & gdap_SO.G2latitude<geo_lims.lat_lims_ind(2) & ...
    gdap_SO.G2longitude>=geo_lims.ind_lims(1) & gdap_SO.G2longitude<geo_lims.ind_lims(2)) & ...
    temp_date_vec(:,2)>=month_lims_for_mean(1) & temp_date_vec(:,2)<=month_lims_for_mean(end));

temp_lats = gdap_SO.G2latitude(geo_index);
temp_lats_nearest_pt5 = round(temp_lats.*.4)./.4;

temp_lons = gdap_SO.G2longitude(geo_index);
temp_lons_nearest_pt5 = round(temp_lons.*.4)./.4;

gdap_SO.NCEP_SLP_CLIM = NaN(length(gdap_SO.GMT_Matlab),1);

for d = 1:length(geo_index)
    lat_lon_index = find(total_lat_lon_list(:,1)==temp_lats_nearest_pt5(d) & total_lat_lon_list(:,2)==temp_lons_nearest_pt5(d));
    lat_lon_index = lat_lon_index(1); % temp cludge b/c there are duplicates in the lat/lon list
    
    year_day = gdap_SO.GMT_Matlab(geo_index(d)) - datenum(temp_date_vec(geo_index(d),1), 1, 0);
    
    day_index = min(abs(press_year_day - year_day))==abs(press_year_day - year_day);
    
    gdap_SO.NCEP_SLP_CLIM(geo_index(d)) = nanmean(select_seasonal_press(lat_lon_index, day_index));
    
end

% gdap_SO.SLP = ncep_pressure_matching(gdap_SO.GMT_Matlab, gdap_SO.G2latitude, gdap_SO.G2longitude, []);
pC_moist_atmosphere = gdap_SO.pC_dry_atmosphere.*(gdap_SO.NCEP_SLP_CLIM./1013.25-pH2O); % uatm wet air

gdap_SO.pC_moist_atmosphere = pC_moist_atmosphere;
gdap_SO.delta_pCO2_unified = gdap_SO.G2pCO2 - gdap_SO.pC_moist_atmosphere;

clear pH2O pC_moist_atmosphere pC_dry_atmosphere decile count d decile_vec lat_index date_index year_day day_index lat_lon_index temp_lats temp_lons temp_lats_nearest_pt5 temp_lons_nearest_pt5
clear d q temp_date_vec geo_index decile_count
% SOCAT delta pCO2 calculation
disp('Calculating SOCAT delta pCO2')
socat.calc_pCO2 = pCO2_from_fCO2(socat.fCO2rec, socat.SST);
socat.GMT_Matlab = datenum(socat.yr, socat.mon, socat.day, socat.hh, socat.mm, socat.ss);
socat.pdens = sw_pden(socat.sal, socat.SST, 0, 0);

% loop through datapoints, matching lat/time
pC_dry_atmosphere = NaN(length(socat.GMT_Matlab),1);
decile_vec = round(length(socat.GMT_Matlab)./10).*[1:10];
decile_count = 1;

for d = 1:length(socat.GMT_Matlab)
    date_index = min(abs(noaa_co2.Matlab_time - socat.GMT_Matlab(d)))==abs(noaa_co2.Matlab_time - socat.GMT_Matlab(d));
    %     disp(' check that this code is working!!!!')
    if isnan(socat.latitude(d)) || socat.GMT_Matlab(d)>noaa_co2.Matlab_time(end) + 30 || socat.GMT_Matlab(d)<noaa_co2.Matlab_time(1) - 30% double check that this is working
        pC_dry_atmosphere(d) = nan;
    else
        lat_index = min(abs(noaa_co2.lat_co2 - socat.latitude(d)))==abs(noaa_co2.lat_co2 - socat.latitude(d));
        
        pC_dry_atmosphere(d) = nanmean(reshape(noaa_co2.co2_atmos(lat_index, date_index),[],1));
    end
    if d>decile_vec(decile_count)
        disp([num2str((d/length(socat.longitude))*100,4) ' % done'])
        decile_count = decile_count+1;
    end
end
%
%     S_CO2 = CO2_Sol(Temperature,Salinity); % % mmol m-3 atm-1

pH2O = ph2osat_smb(socat.SST, socat.sal);

temp_lats = socat.latitude;
temp_lats_nearest_pt5 = round(temp_lats.*.4)./.4;

temp_lons = socat.longitude;
temp_lons_nearest_pt5 = round(temp_lons.*.4)./.4;

temp_date_vec = datevec(socat.GMT_Matlab);

socat.NCEP_SLP_CLIM = NaN(length(socat.GMT_Matlab),1);

for d = 1:length(socat.GMT_Matlab)
    lat_lon_index = find(total_lat_lon_list(:,1)==temp_lats_nearest_pt5(d) & total_lat_lon_list(:,2)==temp_lons_nearest_pt5(d));
    if isempty(lat_lon_index)
        continue
    end
    
    lat_lon_index = lat_lon_index(1); % temp cludge b/c there are duplicates in the lat/lon list
    
    year_day = socat.GMT_Matlab(d) - datenum(temp_date_vec(d,1), 1, 0);
    
    day_index = min(abs(press_year_day - year_day))==abs(press_year_day - year_day);
    
    socat.NCEP_SLP_CLIM(d) = nanmean(select_seasonal_press(lat_lon_index, day_index));
    
end
clear year_day lat_lon_index day_index d


pC_moist_atmosphere = pC_dry_atmosphere.*(socat.NCEP_SLP_CLIM./1013.25-pH2O); % uatm wet air

socat.pC_dry_atmosphere = pC_dry_atmosphere;
socat.pC_moist_atmosphere = pC_moist_atmosphere;
socat.delta_pCO2_unified = socat.calc_pCO2 - socat.pC_moist_atmosphere;

clear pH2O pC_moist_atmosphere pC_dry_atmosphere decile_count date_index lat_index temp_date_vec temp_lats temp_lons temp_lats_nearest_pt5 temp_lons_nearest_pt5
clear noaa_co2_dir noaa_co2_filename noaa_co2 decile_vec gdap_fields_w_flags gdap_fields_no_flags press_dates press_year_day select_seasonal_press total_lat_lon_list

%%  Get all float data into two bins - Pacific and Indian
% only criteria for separation is longitude
% eventually will create "store_argo", "store_gdap", "store_socat", which
% have all property observations in vectors/arrays (for 2d data) separated
% into each region. Used for later calculation of binned formation
% properties and in Figure 1
disp('Storing data by region')
clear store_argo
regions = {'pacific', 'indian'};

parameters_1D = {'GMT_Matlab', 'Lat', 'Lon', 'MLD', 'ml_o2', 'ml_DIC', 'ml_pCO2', 'ml_delta_o2', 'ml_gamma', 'ml_no3', 'Oxygen_flux', 'ml_pdens', 'ml_Temp_C', 'ml_Sal', 'delta_pCO2_unified', 'float_SN', 'CO2_flux', 'ml_PTMP', 'SLP'};
parameters_2D = {'Press_db', 'Temp_C', 'Sal', 'PDENS', 'o2_umol_kg', 'Nitrate', 'pCO2', 'DIC', 'AOU'};
clear region_obs
daily_flag = 0; % using original profiles, not interpolated daily values
disp([' daily flag ' num2str(daily_flag)]) 

if daily_flag==0
    parameter_list = parameters_1D;
    parameter_list_2D = parameters_2D;

else
    parameter_list = parameters_1D_daily;
    parameter_list_2D = parameters_2D;
end

% first calculate a mixed layer potential density

for f = 1:length(SO_SNs)
    if ~isfield(Argo.(SO_SNs{f}), 'MLD')
        continue
    end
    Argo.(SO_SNs{f}).ml_pdens = NaN(length(Argo.(SO_SNs{f}).GMT_Matlab),1);
    Argo.(SO_SNs{f}).ml_Temp_C = NaN(length(Argo.(SO_SNs{f}).GMT_Matlab),1);
    Argo.(SO_SNs{f}).ml_Sal = NaN(length(Argo.(SO_SNs{f}).GMT_Matlab),1);
    Argo.(SO_SNs{f}).ml_PTMP = NaN(length(Argo.(SO_SNs{f}).GMT_Matlab),1);

    Argo.(SO_SNs{f}).float_SN = NaN(length(Argo.(SO_SNs{f}).GMT_Matlab),1);

    Argo.(SO_SNs{f}).float_SN(:) = str2double((SO_SNs{f}(2:end)));
    for p=1:length(Argo.(SO_SNs{f}).GMT_Matlab)
        Argo.(SO_SNs{f}).ml_pdens(p) = nanmean(Argo.(SO_SNs{f}).PDENS(p, Argo.(SO_SNs{f}).Press_db(p,:)<=Argo.(SO_SNs{f}).MLD(p)));
        Argo.(SO_SNs{f}).ml_Temp_C(p) = nanmean(Argo.(SO_SNs{f}).Temp_C(p, Argo.(SO_SNs{f}).Press_db(p,:)<=Argo.(SO_SNs{f}).MLD(p)));
        Argo.(SO_SNs{f}).ml_PTMP(p) = nanmean(Argo.(SO_SNs{f}).PTMP(p, Argo.(SO_SNs{f}).Press_db(p,:)<=Argo.(SO_SNs{f}).MLD(p)));

        Argo.(SO_SNs{f}).ml_Sal(p) = nanmean(Argo.(SO_SNs{f}).Sal(p, Argo.(SO_SNs{f}).Press_db(p,:)<=Argo.(SO_SNs{f}).MLD(p)));
        
    end
end

% region_lims = 
SO_lat_lim = -35;
daily_flag=0;

clear store_argo
            
% initialize store_argo arrays for pacific and indian regions
for r=1:length(regions)
    for p1 = 1:length(parameter_list)
        store_argo.(regions{r}).(parameter_list{p1}) = [];
    end
    for p2 = 1:length(parameters_2D)
        store_argo.(regions{r}).(parameters_2D{p2}) = [];
    end
end

clear r p1
bad_flag=0;
% for each float, create an index of profiles in each region, then append
% the data in store_argo vectors for each parameter
for f = 1:size(SO_SNs,1)

    clear region_index

    region_index.pacific = find(Argo.(SO_SNs{f}).Lat>geo_lims.lat_lims_pac(1) & Argo.(SO_SNs{f}).Lat<geo_lims.lat_lims_pac(2) & ...
        Argo.(SO_SNs{f}).Lon>=geo_lims.pac_lims(1) & Argo.(SO_SNs{f}).Lon<geo_lims.pac_lims(2));
    
    region_index.indian = find(Argo.(SO_SNs{f}).Lat>geo_lims.lat_lims_ind(1) & Argo.(SO_SNs{f}).Lat<geo_lims.lat_lims_ind(2) & ...
        Argo.(SO_SNs{f}).Lon>=geo_lims.ind_lims(1) & Argo.(SO_SNs{f}).Lon<geo_lims.ind_lims(2));
    
    
    for r=1:length(regions)
        if ~isempty(region_index.(regions{r}))
            disp(f)
            for p1 = 1:length(parameter_list)
                if isfield(Argo.(SO_SNs{f}), parameter_list{p1})
                    if ~isempty(Argo.(SO_SNs{f}).(parameter_list{p1}))
                        store_argo.(regions{r}).(parameter_list{p1}) = [store_argo.(regions{r}).(parameter_list{p1}); Argo.(SO_SNs{f}).(parameter_list{p1})(region_index.(regions{r}),1)];
                    end

                else
                    
                    temp_vec = NaN(length(region_index.(regions{r})),1);
                    store_argo.(regions{r}).(parameter_list{p1}) = [store_argo.(regions{r}).(parameter_list{p1}); temp_vec];
                    
                end
                clear temp_vec
            end
            
            for p2 = 1:length(parameters_2D)
                if isfield(Argo.(SO_SNs{f}), parameters_2D{p2})
                    % While SOCCOM data are all the same size, the UW ARGO
                    % data are not.  This gives every profile a length of
                    % 1000 - wastes space, but an OK fix for now
                    temp_data = Argo.(SO_SNs{f}).(parameters_2D{p2})(region_index.(regions{r}),:);
                    temp_vec = NaN( size(temp_data,1), 1000);
                    temp_vec(:,1:size(temp_data,2)) = temp_data;
                    
                    store_argo.(regions{r}).(parameters_2D{p2}) = [store_argo.(regions{r}).(parameters_2D{p2}); temp_vec];
                    
                else % if the field does not exist, need to add that much NaN space
                    temp_vec = NaN(length(region_index.(regions{r})),1000);
                    store_argo.(regions{r}).(parameters_2D{p2}) = [store_argo.(regions{r}).(parameters_2D{p2}); temp_vec];
                    
                end
                
                if size(store_argo.(regions{r}).(parameters_2D{p2}),1) ~= size(store_argo.(regions{r}).GMT_Matlab,1)
                    disp(f)
                    bad_flag =1;
                    break
                end
                clear temp_vec temp_data
            end
        end
        if bad_flag==1
            break
        end
    end
    clear region_index
    if bad_flag==1
        break
    end
end
%
clear f p1 p

for r=1:length(regions)
    
    temp_vec = datevec(store_argo.(regions{r}).GMT_Matlab);
    store_argo.(regions{r}).date_vec = temp_vec;
    store_argo.(regions{r}).year = temp_vec(:,1);
    
    store_argo.(regions{r}).month = temp_vec(:,2);
    clear temp_vec
    
    store_argo.(regions{r}).ml_no3(store_argo.(regions{r}).ml_no3<-10) = nan;
end

clear r daily_flag parameters_1D parameters_2D p2 bad_flag

%% ML calculations for shipboard data
disp('socat ML calculations')
% load gridded Argo-derived MLDs and perform calculations on Glodap and Socat
% used for Figure 1
argo_MLD = load([home_dir 'Work/Projects/2020_02 SO SAMW Variability_Cerovecki/Data_from_Ivana/2021_12_10/MLD_av_2005_OCT_2021_SouthOc.mat']);

argo_MLD.GMT_Matlab = datenum(2005,1:size(argo_MLD.mld,3),15); % first date is Jan 2005 and data is monthly so make a time vector according to the size of the mld array
argo_MLD.datevec = datevec(argo_MLD.GMT_Matlab);

% Match SOCAT points to gridded Argo MLDs 
socat.MLD = NaN(length(socat.longitude),1);
argo_MLD.latvec = argo_MLD.LAT(1,:);
argo_MLD.lonvec = argo_MLD.LON(:,1);
argo_MLD.wint_mean = nanmean(argo_MLD.mld(:,:,argo_MLD.datevec(:,2)>=month_lims_for_mean(1) & argo_MLD.datevec(:,2)<=month_lims_for_mean(2) & ...
    argo_MLD.datevec(:,1)<=2020),3);

decile_vec = round(length(socat.GMT_Matlab)./10).*[1:10];
decile_count = 1;
for s = 1:length(socat.longitude)
    
    date_index = socat.yr(s)==argo_MLD.datevec(:,1) & socat.mon(s)==argo_MLD.datevec(:,2);
    if sum(date_index)==0
        continue
    end
    
    lat_index = min(abs(argo_MLD.latvec - socat.latitude(s)))==abs(argo_MLD.latvec - socat.latitude(s));
    
    if (abs(socat.latitude(s) - argo_MLD.latvec(lat_index))>1) % if the match is too far off or if the latitude exceeds 28S, skip
        continue
    end
    lon_index = min(abs(argo_MLD.lonvec - socat.longitude(s)))==abs(argo_MLD.lonvec - socat.longitude(s));
    
    temp_MLD = argo_MLD.mld(lon_index, lat_index, date_index);
    socat.MLD(s) = nanmean(reshape(temp_MLD,1,[]));
    
    clear date_index lat_index lon_index temp_MLD
    %     if s==round(length(socat.longitude)/10)
    %         disp([num2str((s/length(socat.longitude))*100,4) ' % done'])
    %     end
    if s>decile_vec(decile_count)
        disp([num2str((s/length(socat.longitude))*100,4) ' % done'])
        decile_count = decile_count+1;
    end
end

clear decile_count decile_vec s 

disp('socat ML calcs done')
% calculate ML means for GDAP
gdap_para = {'GMT_Matlab'; 'G2latitude';'G2longitude';'G2MLD';'G2oxygen';'G2tco2';'G2pCO2';'G2deltaO2';...
    '';'G2nitrate';'';'PDENS';'G2temperature';'G2salinity';'delta_pCO2_unified';'';'';'G2theta'; 'G2cruise'};

for p1 = 1:length(gdap_para)
    if ~isempty(gdap_para{p1})
        
        gdap_SO.ML_avg.(gdap_para{p1}) = [];
    end
end

for p1 = 1:length(gdap_para)
    if ~isempty(gdap_para{p1})
        
        gdap_SO.Surf_avg.(gdap_para{p1}) = [];
    end
end


% Gdap mixed layer calcs - calculate MLD for each cruise, then calculate ML
% averages and Surface (<= 25m) averages
disp('Performing mixed layer calculations for glodap')
cruise_numbers = unique(gdap_SO.G2cruise);
stations_per_cruise = NaN(size(cruise_numbers));
for c = 1:length(cruise_numbers)
    cruise_index = find(gdap_SO.G2cruise==cruise_numbers(c));
    
    station_numbers =  unique(gdap_SO.G2station(cruise_index));
    
    stations_per_cruise(c) = numel(station_numbers);
    for s=1:length(station_numbers)
        station_index = find(gdap_SO.G2cruise==cruise_numbers(c) & gdap_SO.G2station==station_numbers(s));
        
        temp_T_0 = gdap_SO.G2temperature(station_index);
        temp_S_0 = gdap_SO.G2salinity(station_index);
        temp_P_0 = gdap_SO.G2pressure(station_index);
        
        temp_T = temp_T_0(~isnan(temp_T_0) & ~isnan(temp_S_0) & ~isnan(temp_P_0));
        temp_S = temp_S_0(~isnan(temp_T_0) & ~isnan(temp_S_0) & ~isnan(temp_P_0));
        temp_P = temp_P_0(~isnan(temp_T_0) & ~isnan(temp_S_0) & ~isnan(temp_P_0));
        
        temp_MLD = mld_dbm(temp_T, temp_S, ...
            temp_P, 0);
        
        % MLD averages
        for p1 = 1:length(gdap_para)
            if ~isempty(gdap_para{p1})
                if ~strcmp(gdap_para{p1}, 'G2MLD')
                    temp_var = gdap_SO.(gdap_para{p1})(station_index);
                    
                    gdap_SO.ML_avg.(gdap_para{p1}) = [gdap_SO.ML_avg.(gdap_para{p1}); nanmean(temp_var(temp_P_0<=temp_MLD))];
                else
                    gdap_SO.ML_avg.(gdap_para{p1}) = [gdap_SO.ML_avg.(gdap_para{p1}); temp_MLD];
                end
            end
        end
        
        % upper 25 m averages
        for p1 = 1:length(gdap_para)
            if ~isempty(gdap_para{p1})
                if ~strcmp(gdap_para{p1}, 'G2MLD')
                    temp_var = gdap_SO.(gdap_para{p1})(station_index);
                    
                    gdap_SO.Surf_avg.(gdap_para{p1}) = [gdap_SO.Surf_avg.(gdap_para{p1}); nanmean(temp_var(temp_P_0<=25))];
                else
                    gdap_SO.Surf_avg.(gdap_para{p1}) = [gdap_SO.Surf_avg.(gdap_para{p1}); temp_MLD]; % still store temp_MLD b/c that's the only way you have of knowing whether the data are part of a deep ML
                end
            end
        end
        
        % % code for plotting ML calculations
        %         if c==2
        %         clf
        %         plot(gdap_SO.G2oxygen(station_index), gdap_SO.G2pressure(station_index), 'x-'); set(gca, 'ydir', 'reverse')
        %         hold on
        %         plot(gdap_SO.ML_avg.G2oxygen(end), gdap_SO.ML_avg.G2MLD(end), 'ro')
        %         set(gca, 'ylim', [0 500])
        %         pause
        %         end
    end
    
end
disp('Gdap ML calcs done')

clear temp_T_0 temp_S_0 temp_P_0 temp_T temp_S temp_P temp_MLD temp_var p1 c s cruise_index station_numbers cruise_numbers press_dates press_year_day station_index stations_per_cruise
%% put shipboard data into regions
disp('Separating shipboard data into regions')
% only criteria for separation is longitude

clear store_socat store_gdap

socat_para = {'GMT_Matlab';'latitude';'longitude';'MLD';'';'';'calc_pCO2'};
socat_para(12) = {'pdens'};
socat_para(13) = {'SST'};
socat_para(14) = {'sal'};
socat_para(15) = {'delta_pCO2_unified'};
socat_para(16) = {''};
socat_para(17) = {''};
socat_para(18) = {'SST'};
socat_para(19) = {'NCEP_SLP'};

% initialize "store" structures
for r=1:length(regions)
    for p1 = 1:length(parameter_list)
        store_socat.(regions{r}).(parameter_list{p1}) = [];
        store_gdap.(regions{r}).(parameter_list{p1}) = [];
        
        store_gdap.surf.(regions{r}).(parameter_list{p1}) = [];
        
    end
end

clear r p1

% find matching GLODAP data
region_index.pacific = find(gdap_SO.ML_avg.G2latitude>geo_lims.lat_lims_pac(1) & gdap_SO.ML_avg.G2latitude<geo_lims.lat_lims_pac(2) & ...
    gdap_SO.ML_avg.G2longitude>=geo_lims.pac_lims(1) & gdap_SO.ML_avg.G2longitude<geo_lims.pac_lims(2));

region_index.indian = find(gdap_SO.ML_avg.G2latitude>geo_lims.lat_lims_ind(1) & gdap_SO.ML_avg.G2latitude<geo_lims.lat_lims_ind(2) & ...
    gdap_SO.ML_avg.G2longitude>=geo_lims.ind_lims(1) & gdap_SO.ML_avg.G2longitude<geo_lims.ind_lims(2));


for r=1:length(regions)
    if isempty(region_index.(regions{r}))
        continue
    end
    for p1 = 1:length(parameter_list)
        if ~isempty(gdap_para{p1})
            store_gdap.(regions{r}).(parameter_list{p1}) = gdap_SO.ML_avg.(gdap_para{p1})(region_index.(regions{r}));
            store_gdap.surf.(regions{r}).(parameter_list{p1}) = gdap_SO.Surf_avg.(gdap_para{p1})(region_index.(regions{r}));
            
        end
    end
    store_gdap.(regions{r}).G2cruise = gdap_SO.ML_avg.G2cruise(region_index.(regions{r}));

end
clear region_index r

clear p1 p
% find matching SOCAT data
region_index.pacific = find(socat.latitude>geo_lims.lat_lims_pac(1) & socat.latitude<geo_lims.lat_lims_pac(2) & ...
    socat.longitude>=geo_lims.pac_lims(1) & socat.longitude<geo_lims.pac_lims(2));

region_index.indian = find(socat.latitude>geo_lims.lat_lims_ind(1) & socat.latitude<geo_lims.lat_lims_ind(2) & ...
    socat.longitude>=geo_lims.ind_lims(1) & socat.longitude<geo_lims.ind_lims(2));

for r=1:length(regions)
    if ~isempty(region_index.(regions{r}))
        for p1 = 1:length(parameter_list)
            if ~isempty(socat_para{p1})
                store_socat.(regions{r}).(parameter_list{p1}) = socat.(socat_para{p1})(region_index.(regions{r}));
                
            end
        end
        
    end
    
    store_socat.(regions{r}).expocode = socat.Expocode(region_index.(regions{r}));
end
clear region_index r

clear p1 p

% add some date vectors, years, and months to the store_gdap and
% store_socat structures
for r=1:length(regions)
    
    temp_vec = datevec(store_gdap.(regions{r}).GMT_Matlab);
    store_gdap.(regions{r}).date_vec = temp_vec;
    store_gdap.(regions{r}).year = temp_vec(:,1);
    store_gdap.(regions{r}).month = temp_vec(:,2);
    clear temp_vec
    
    temp_vec = datevec(store_gdap.surf.(regions{r}).GMT_Matlab);
    store_gdap.surf.(regions{r}).date_vec = temp_vec;
    store_gdap.surf.(regions{r}).year = temp_vec(:,1);
    store_gdap.surf.(regions{r}).month = temp_vec(:,2);
    clear temp_vec
    
    temp_vec = datevec(store_socat.(regions{r}).GMT_Matlab);
    store_socat.(regions{r}).date_vec = temp_vec;
    store_socat.(regions{r}).year = temp_vec(:,1);
    store_socat.(regions{r}).month = temp_vec(:,2);
    clear temp_vec
end
clear r

clear parameter_list parameter_list_2D socat_para gdap_para
disp('Done')

%% load mode water volumes
% volumes are stored in 0.05 bins 

past_final_year = '2021'; 
clear mode_vol
disp('Loading density-binned MW volumes')
% load annual MW volume files calculated from core Argo
mode_vol_dir = [cerovecki_dir 'Data_from_Ivana/MW_Vol_RG_Argo/'];
mode_files = dir([mode_vol_dir 'MW_Vol*.mat']);

% mode_vol array created - each mode_vol year has volume binned in 0.05
% density bins and size(360 - lon, 65 - lat, 25 - dens, 12 - month)
for m = 1:length(mode_files)
    mode_year = mode_files(m).name(15:18);
    
    % only want data through to 2020 for now to be consistent with
    % everything else
    if strcmp(mode_year, past_final_year)
        clear mode_year

        break
    end
    mode_vol.(['y' mode_year]) = load([mode_vol_dir mode_files(m).name]);
    
    % set first and last densities to nans - these are not real numbers
    % according to Ivana
    mode_vol.(['y' mode_year]).vol_vs_sigth(:,:,1,:) = nan;
    mode_vol.(['y' mode_year]).vol_vs_sigth(:,:,end,:) = nan;

    clear mode_year
    
    
end
clear m mode_files mode_vol_dir
%

mode_years = fieldnames(mode_vol);

dens_interval = mean(diff(mode_vol.(mode_years{1}).sig_th));

% calculate a mean MW volume for plotting (only selected months)
temp_vol_array = NaN(length(mode_years), length(mode_vol.(mode_years{1}).X), length(mode_vol.(mode_years{1}).Y));

for my = 1:length(mode_years)
% first average over selected months, then sum along all densities   
    temp_vol_array(my,:,:) = nansum(nanmean(mode_vol.(mode_years{my}).vol_vs_sigth(:,:,:,month_lims_for_mean),4),3);
end
%
% then average along all years:
mode_vol.avg = squeeze(nanmean(temp_vol_array,1));

% volume of mode water formed and mean concentration
properties = {'ml_o2' 'ml_delta_o2' 'ml_DIC' 'ml_no3'  'ml_pCO2', 'ml_PTMP', 'ml_Sal', 'MLD', 'ml_pdens', 'delta_pCO2_unified', 'Oxygen_flux', 'CO2_flux'};

% latitude index does not change for each region
clear mw_prop temp_vol_array

% initialize mw_prop structures for each region
for r = 1:length(regions)
    
    mw_prop.(regions{r}).vol = NaN(length(mode_years), 12, length(mode_vol.y2005.sig_th),1); % (year, month, density)
    
    for p = 1:length(properties)
        mw_prop.(regions{r}).(properties{p}) = NaN(length(mode_years), 12, length(mode_vol.y2005.sig_th),2);
    end
end

% Store summed volume and mean/std Argo BGC properties for each density bin
% across entire lat/lon range
% will create ml_[property] arrays in pacific and indian regions. Each
% ml_[property] has size (16 - year, 12 - month, 25 - density bin, 2 - mean/std
% dev). 
% Selection criteria applied: density, month, year, ML cutoff
for r = 1:length(regions)
    
    if strncmp((regions{r}), 'pacific', 7)
        lon_lims = geo_lims.pac_lims;
        lat_lims = geo_lims.lat_lims_pac;
    elseif strncmp((regions{r}), 'indian',6)
        lon_lims = geo_lims.ind_lims;
        lat_lims = geo_lims.lat_lims_ind;
    end
    
    % longitude index for mode water volume and argo property storage
    lon_index = mode_vol.y2005.X>=lon_lims(1) & mode_vol.y2005.X < lon_lims(2);
    lat_index = mode_vol.y2005.Y>=lat_lims(1) & mode_vol.y2005.Y<=lat_lims(2);

    % loop through each year
    for my = 1:length(mode_years)
        
        m_year = str2double(mode_years{my}(2:end));
        
        for month = 1:12
            
            % steps through the density bins loaded from the mode water
            % volumes that Ivana sent - %2021_12 update, 0.05 bins 
            for s = 1:length(mode_vol.(mode_years{my}).sig_th)-1
                % sums volume across lon and lat, while stepping through
                % year, month, lon, lat
                mw_prop.(regions{r}).vol(my, month, s) = nansum(reshape(mode_vol.(mode_years{my}).vol_vs_sigth(lon_index,lat_index,s,month),1,[]));
                
                % finds argo data that 
                % fall into the original density bins of the loaded SAMW
                % volume for each month and year - this data is already
                % stored by region and does not need lat/lon criteria at
                % this point.
                % Ivana's volumes correspond to >= sigth(s) & < sigth(s+1)

                argo_index = ... %store_argo.(regions{r}).Lon>=lon_lims(1) & store_argo.(regions{r}).Lon<lon_lims(2) & store_argo.(regions{r}).Lat>=lat_lims(1) & store_argo.(regions{r}).Lat<lat_lims(2) ...
                    store_argo.(regions{r}).ml_pdens-1000>=(mode_vol.(mode_years{my}).sig_th(s))  & store_argo.(regions{r}).ml_pdens-1000<mode_vol.(mode_years{my}).sig_th(s+1)  & ...
                    store_argo.(regions{r}).month == month & store_argo.(regions{r}).year == m_year  & store_argo.(regions{r}).MLD>MLD_cutoff;
                
                for p = 1:length(properties) % stores the argo mixed layer data
                    mw_prop.(regions{r}).(properties{p})(my, month, s,1) = nanmean(store_argo.(regions{r}).(properties{p})(argo_index));
                    mw_prop.(regions{r}).(properties{p})(my, month, s,2) = nanstd(store_argo.(regions{r}).(properties{p})(argo_index));
                end
                clear argo_index
            end
            
        end
        clear m_year
    end
    clear lon_index lon_lims lat_lims
end
clear my month r s p lat_index 


% Similar to above, but now calculating mean bgc properties for all months
% and all years and for winter months only. 

% initialize "all" structures - this combines all years, all months across
% entire geographic range, means/std per density bin
% Initialize "wint" structures, which is the same but only for winter
% months defined at the start of this script (month_lims_for_mean)
% "wint structures are used for vol/prop vs. density figures (Fig. 3) and
% associated table 1
% properties in each structure has size (25-density bins, 2 - mean/std)
for r = 1:length(regions)
    mw_prop.(regions{r}).all.vol = NaN(12, length(mode_vol.y2005.sig_th),2); % (month, density, mean/std)
        
    for p = 1:length(properties)
        mw_prop.(regions{r}).all.(properties{p}) = NaN(length(mode_vol.y2005.sig_th),2); % (Density, mean/std)
        mw_prop.(regions{r}).wint.(properties{p}) = NaN(length(mode_vol.y2005.sig_th),2);

    end
    mw_prop.(regions{r}).wint.vol = NaN(length(mode_vol.y2005.sig_th),2); % (Density, mean/std)

end

for r = 1:length(regions)
    
    % applying latitude and longitude criteria for MW volume
    if strncmp((regions{r}), 'pacific', 7)
        lon_lims = geo_lims.pac_lims;
        lat_lims = geo_lims.lat_lims_pac;

        disp('pac')
    elseif strncmp((regions{r}), 'indian',6)
        lon_lims = geo_lims.ind_lims;
        lat_lims = geo_lims.lat_lims_ind;

        disp('ind')
    end
    lon_index = mode_vol.y2005.X>=lon_lims(1) & mode_vol.y2005.X < lon_lims(2);
    lat_index = mode_vol.y2005.Y>=lat_lims(1) & mode_vol.y2005.Y<=lat_lims(2);

    clear lat_lims lon_lims
    % loops through densities from initial mode water volume data
    for s = 1:length(mode_vol.(mode_years{1}).sig_th)-1 
        % create a monthly climatology of volume in density bins
        for month = 1:12
            
            temp_vols = [];
            
            for my = 1:length(mode_years) % adding up all volumes for a given density and month across years
                temp_vols = [temp_vols ; nansum(reshape(mode_vol.(mode_years{my}).vol_vs_sigth(lon_index, lat_index, s, month),[],1))];
            end
            clear my
            mw_prop.(regions{r}).all.vol(month, s,1) = nanmean(temp_vols);
            mw_prop.(regions{r}).all.vol(month, s,2) = nanstd(temp_vols);
            clear temp_vols
            
        end
        clear month
               
        % create an annual mean of bgc properties for each density
        % Ivana's volumes correspond to >= sigth(s) & < sigth(s+1)
        argo_index = store_argo.(regions{r}).ml_pdens-1000>=(mode_vol.(mode_years{1}).sig_th(s))  & store_argo.(regions{r}).ml_pdens-1000<mode_vol.(mode_years{1}).sig_th(s+1) ...
            & store_argo.(regions{r}).MLD>MLD_cutoff;

        for p = 1:length(properties)
            mw_prop.(regions{r}).all.(properties{p})(s,1) = nanmean(store_argo.(regions{r}).(properties{p})(argo_index));
            mw_prop.(regions{r}).all.(properties{p})(s,2) = nanstd(store_argo.(regions{r}).(properties{p})(argo_index));
        end
        clear p argo_index
        
        % create a wintertime mean of bgc properties for each density using
        % the month_lim criteria - filtered using density bin, winter
        % months, MLD deeper than cutoff
        argo_index = store_argo.(regions{r}).ml_pdens-1000>=(mode_vol.(mode_years{1}).sig_th(s)) & store_argo.(regions{r}).ml_pdens-1000<mode_vol.(mode_years{1}).sig_th(s+1) ... 
            & store_argo.(regions{r}).month>=month_lims_for_mean(1) & store_argo.(regions{r}).month<=month_lims_for_mean(2) & store_argo.(regions{r}).MLD>MLD_cutoff;
       
        for p = 1:length(properties)
            mw_prop.(regions{r}).wint.(properties{p})(s,1) = nanmean(store_argo.(regions{r}).(properties{p})(argo_index));
            mw_prop.(regions{r}).wint.(properties{p})(s,2) = nanstd(store_argo.(regions{r}).(properties{p})(argo_index));
            
        end
        clear p argo_index
        
    end
    clear s lon_index
    
  
    % wintertime mean volume using month_lim criteria
    vol_mean_temp = NaN(length(mode_years), size(mw_prop.(regions{r}).vol,3));
    for my = 1:length(mode_years)
        vol_mean_temp(my,:) = nanmean(squeeze(mw_prop.(regions{r}).vol(my,  month_lims_for_mean(1):month_lims_for_mean(2),:)),1);
        
    end
    mw_prop.(regions{r}).wint.vol(:,1) = nanmean(vol_mean_temp,1);
    mw_prop.(regions{r}).wint.vol(:,2) = nanstd(vol_mean_temp,1);

    clear vol_mean_temp my
end
clear r
clear lat_index
%%

% finding binned mode water volume percentages from the wintertime months used in this
% study
% vol_5_per_index is the index of volume bins that fall within the >=5%
% criteria
% sig_thet_range is the range of densities that satisfy the >=5% criteria
for r = 1:length(regions)
    mw_prop.(regions{r}).all.vol_per = NaN(12,length(mode_vol.(mode_years{1}).sig_th)); % initialize array: (month, density)

    % for each winter month, calculate the percentage of total SAMW volume in each
    % density bin
    for month = month_lims_for_mean(1):month_lims_for_mean(end)
        mw_prop.(regions{r}).all.vol_per(month,:) = mw_prop.(regions{r}).all.vol(month,:,1)./nansum(mw_prop.(regions{r}).all.vol(month,:,1));
    end
    
    % take the average across all months
    mw_prop.(regions{r}).all.vol_per_mean = nanmean(mw_prop.(regions{r}).all.vol_per,1); %(density)
    
    % find those density bins with more than 5% of MW volume, rounding to
    % the nearest percent so that I don't cut anything off that is close
    mw_prop.(regions{r}).all.per_5_index = mw_prop.(regions{r}).all.vol_per_mean>=.05;
    
    per_5_temp = find(mw_prop.(regions{r}).all.per_5_index);
    
    mw_prop.(regions{r}).all.sig_thet_range = [mode_vol.y2005.sig_th(mw_prop.(regions{r}).all.per_5_index) mode_vol.y2005.sig_th(per_5_temp(end)+1)];
%         mw_prop.(regions{r}).all.sig_thet_vals = mode_vol.y2005.sig_th(mw_prop.(regions{r}).all.per_5_index);

clear per_5_temp
end
clear r month 
%%
% save the mean value in mode water volumes that are greater than 5 %
% currently "hard coded" values for month range and density bins
for r = 1:length(regions)
    
    % initialize array for means over several months, both for unweighted
    % means and volume weighted means
    % search criteria are: density range, month limits (winter months), MLD
    % cutoff
    mw_prop.(regions{r}).all.mean_prop = [];
    mw_prop.(regions{r}).all.mean_prop_weighted = [];

    %initialize array for annual means
    % search criteria are: density range, month limits (winter months),
    % year, MLD cutoff
    mw_prop.(regions{r}).all.prop_yr_mean = NaN(length(properties), 15);
    mw_prop.(regions{r}).all.prop_yr_std = NaN(length(properties), 15);
    
    % save out the float SNs that fall within each region
    mw_prop.(regions{r}).all.winter_SNs = [];
    
    % for each property, go back to original argo dataset and save out
    % properties.
    for p = 1:length(properties)
        
        % overall mean for select months
        argo_index = store_argo.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(1) & ...
            store_argo.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) ...
            & store_argo.(regions{r}).month>=month_lims_for_mean(1) & store_argo.(regions{r}).month<=month_lims_for_mean(2) & store_argo.(regions{r}).MLD>MLD_cutoff;
        
        mw_prop.(regions{r}).all.mean_prop.(properties{p})(1) = nanmean(store_argo.(regions{r}).(properties{p})(argo_index));
        mw_prop.(regions{r}).all.mean_prop.(properties{p})(2) = nanstd(store_argo.(regions{r}).(properties{p})(argo_index));
        
        % save the float SNs that fit these criteria
        mw_prop.(regions{r}).all.winter_SNs = unique([mw_prop.(regions{r}).all.winter_SNs ; unique(store_argo.(regions{r}).float_SN(argo_index))]);
        
        % create a volume weighted mean
        temp_prop = NaN(length(mw_prop.(regions{r}).all.sig_thet_range)-1,2);
        for s = 1:length(mw_prop.(regions{r}).all.sig_thet_range)-1
            argo_index = store_argo.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(s) &  ...
                store_argo.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(s+1) ...
                & store_argo.(regions{r}).month>=month_lims_for_mean(1) & store_argo.(regions{r}).month<=month_lims_for_mean(2) & store_argo.(regions{r}).MLD>MLD_cutoff;
            
            temp_prop(s,1) = nanmean(store_argo.(regions{r}).(properties{p})(argo_index));
            temp_prop(s,2) = nanstd(store_argo.(regions{r}).(properties{p})(argo_index));

        end
        prop_weights = mw_prop.(regions{r}).all.vol_per_mean(mw_prop.(regions{r}).all.per_5_index)./sum(mw_prop.(regions{r}).all.vol_per_mean(mw_prop.(regions{r}).all.per_5_index));
        
        mw_prop.(regions{r}).all.mean_prop_weighted.(properties{p})(1) = sum(temp_prop(:,1).*prop_weights'); % weighted mean
        mw_prop.(regions{r}).all.mean_prop_weighted.(properties{p})(2) = sum(temp_prop(:,2).*prop_weights'); % std

        clear argo_index prop_weights temp_prop 
        % annual means
        for y = 1:length(mode_years)

            m_year = str2double(mode_years{y}(2:end));

            argo_index = store_argo.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(1) &  ...
                store_argo.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) ...
                & store_argo.(regions{r}).month>=month_lims_for_mean(1) & store_argo.(regions{r}).month<=month_lims_for_mean(2) & store_argo.(regions{r}).year==m_year  & store_argo.(regions{r}).MLD>MLD_cutoff;
            mw_prop.(regions{r}).all.prop_yr_mean(p,y) = nanmean(store_argo.(regions{r}).(properties{p})(argo_index));
            mw_prop.(regions{r}).all.prop_yr_std(p,y) = nanstd(store_argo.(regions{r}).(properties{p})(argo_index));
            
            clear argo_index m_year

        end
        
    end
    
    clear sig_thet_range p y month
end
clear r
%%

% % % annual climatology of mw properties:
% % for r = 1:2
% %     for p=1:length(properties)
% %         mw_prop.(regions{r}).annual_climatology.(properties{p}) = NaN(12,2);
% %     end
% % end
% % 
% % for r = 1:2
% % %     sig_thet_range = mode_vol.y2005.sig_th(mw_prop.(regions{r}).all.per_5_index);
% %     
% %     
% %     for month=1:12
% %         argo_index =store_argo.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(1) & ...
% %             store_argo.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) ...
% %             & store_argo.(regions{r}).month==month & store_argo.(regions{r}).MLD>MLD_cutoff;
% %         
% %         for p = 1:length(properties)
% %             mw_prop.(regions{r}).annual_climatology.(properties{p})(month, 1) = nanmean(store_argo.(regions{r}).(properties{p})(argo_index));
% %             mw_prop.(regions{r}).annual_climatology.(properties{p})(month, 2) = nanstd(store_argo.(regions{r}).(properties{p})(argo_index));
% %             
% %         end
% %     end
% % end
% % 
% % clear p r s month lat_lims lon_lims argo_index

%% comparison to Carter et al. 2021 properties
% carter = load([cerovecki_dir 'Carter preformed properties/PreformedPropertiesDefault_v1.mat']);

%% Carter et al. 2014 comparison:
% Use nitrate difference and redfield to estimate equivalent DIC and Alk
% differences and calculate a new pCO2

carter.no3_diff = mean([23.7 - 20.32 ; 24.27 - 21.93]);
carter.sst_diff = mean([4.6 - 5.247; 4.17 - 4.68]);

carter.extra_OM_deg = carter.no3_diff.*106./16;

carter.dic_change = carter.extra_OM_deg;
carter.alk_change = -.15*carter.extra_OM_deg;


[DATA,~,~]= CO2SYSSOCCOM([2273.1 2273.5], [2122.2 2129.0] , ...
    1,2,[34.181 34.129], [5.247 4.680], ...
    [5.247 4.680],...
    [0 0],[0 0],[4.87 6.91], [1.44 1.54],1,10,3);

carter.carter_pCO2_orig = (DATA(:,19));

[DATA,~,~]= CO2SYSSOCCOM([2273.1 2273.5] + carter.alk_change, [2122.2 2129.0]+carter.dic_change , ...
    1,2,[34.181 34.129], [5.247 4.680], ...
    [5.247 4.680]+ carter.sst_diff,...
    [0 0],[0 0],[4.87 6.91], [1.44 1.54],1,10,3);

carter.carter_pCO2_new = (DATA(:,19));

