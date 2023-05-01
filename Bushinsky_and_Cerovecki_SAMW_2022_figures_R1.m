% figures for "Subantarctic Mode Water
% Biogeochemical Formation Properties and Interannual Variability", Seth M. Bushinsky1 and Ivana Cerove?ki2
%
% 1Department of Oceanography, School of Ocean and Earth Science and Technology, University of Hawai?i at M?noa, Honolulu, HI
% 2Scripps Institution of Oceanography, University of California San Diego, La Jolla, CA
%
% Corresponding author: Seth Bushinsky (seth.bushinsky@hawaii.edu)
%% First run "Bushinsky_and_Ceroveecki_SAMW_2022_R1.m"

plot_dir = [home_dir 'Work/Manuscripts/2020_05 SAMW formation properties/figures/R1/'];


%% Fig 1 Maps of mean MLD, all winter obs, winter obs w/in density range of SAMW and >200 m MLDs

% load SO front criteria used in Gray et al. 2018 and Bushinsky et al. 2019
% for plotting
so_fronts = load([data_dir 'ARGO_O2_Floats/Front_definitions/Gray_5_regions/regional_boundaries_5zone.mat']);

plot_filename = ['Fig_1_map_Months ' num2str(month_lims_for_mean(1)) ' to ' num2str(month_lims_for_mean(2)) '_v6'];
coast = load('coast');

colormap = brewermap(10, 'Paired');
boundary_width=1.5;
coast_gray = [.8 .8 .8];

marker_size = 2;

clf
set(gcf, 'units', 'inches')
paper_w = 7; paper_h =8;
set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]); clear paper_w paper_h

ph = NaN(3,1);


ph(1) = subplot(3,1,1); % top plot, mean winter ML


colormap = flipud(inferno);
set(gcf, 'colormap', colormap)

lat_map = double(argo_MLD.LAT(1,:))';
lon_map = double(argo_MLD.LON(:,1));

R = georasterref('RasterSize', [length(lat_map) length(lon_map)], 'RasterInterpretation', 'Cells', ...
    'LatitudeLimits', [min(lat_map)-.5 max(lat_map)+.5], 'LongitudeLimits', [min(lon_map)-.5 max(lon_map)+.5]);


data1 = argo_MLD.wint_mean';
temp_output = data1;

% Creates a polar stereographic plot with the origin and bounds below
axesm('lambertstd', 'MapLatLimit',[-80 -37],'MapParallels',[-75 -15],'MapLonLimit',[68 -70])
axis off; framem on; gridm on; mlabel on; plabel on;
setm(gca,'MLabelParallel',-20, 'fontsize', 14)
hold on
meshm(temp_output, R);


geoshow(coast.lat,coast.long,'DisplayType','polygon', 'facecolor', coast_gray)
clear temp_output
geoshow(gca, so_fronts.lat_pf, so_fronts.lon_pf, 'displaytype', 'line', 'linewidth', boundary_width, 'color', 'k') % polar front
geoshow(gca, so_fronts.lat_saf, so_fronts.lon_saf, 'displaytype', 'line', 'linewidth', boundary_width, 'color', 'k') % Subantarctic Front

setm(gca, 'mlabellocation', [0 ], 'mlabelparallel', 0)
setm(gca, 'plabellocation', [-180 -180], 'plabelmeridian', -150)
geoshow([-90 -30], [geo_lims.pac_lims(1) geo_lims.pac_lims(1)], 'displaytype', 'line', 'color', 'k', 'linestyle', '--', 'linewidth', 2)


ph(2) = subplot(3,1,2); % plot all winter data

colormap2 = brewermap(10, 'Paired');

hold on
coast = load('coast');
axesm('lambertstd', 'MapLatLimit',[-80 -37],'MapParallels',[-75 -15],'MapLonLimit',[68 -70])
axis off; framem on; gridm on; mlabel on; plabel on;

hold on

% Filter Socat data using winter month criteria and only looking at data after
% 1990
socat_wint_index = (socat.mon>=month_lims_for_mean(1) &  socat.mon<=month_lims_for_mean(length(month_lims_for_mean))) & socat.yr>1990;

s_h = geoshow(socat.latitude(socat_wint_index), socat.longitude(socat_wint_index), ...
    'displaytype', 'multipoint', 'marker', 'd', 'markerfacecolor', colormap2(3,:), 'markeredgecolor', colormap2(4,:), 'markerfacecolor', colormap2(4,:) ,'markersize', marker_size-1);


% plot all winter Argo data, regardless of whether there is pCO2 available or not
for f = 1:length(SO_SNs)
    
    argo_vec = datevec(Argo.(SO_SNs{f}).GMT_Matlab);
    wint_index = argo_vec(:,2)>=month_lims_for_mean(1) & argo_vec(:,2)<=month_lims_for_mean(length(month_lims_for_mean));
    a_h = geoshow(Argo.(SO_SNs{f}).Lat(wint_index), Argo.(SO_SNs{f}).Lon(wint_index), 'displaytype', 'multipoint', 'marker', 'o', 'markeredgecolor', colormap2(2,:),'markerfacecolor', colormap2(2,:), 'markersize', marker_size) ;
    
end

% Filter Glodap data using winter month criteria and only looking at data after
% 1990
gdap_wint_index = (gdap_SO.G2month>=month_lims_for_mean(1) & gdap_SO.G2month<=month_lims_for_mean(length(month_lims_for_mean))) & gdap_SO.G2year>1990;

g_h = geoshow(gdap_SO.G2latitude(gdap_wint_index), gdap_SO.G2longitude(gdap_wint_index), ...
    'displaytype', 'multipoint', 'marker', 's', 'markerfacecolor', colormap2(9,:), 'markeredgecolor', colormap2(10,:), 'markerfacecolor', colormap2(10,:) ,'markersize', marker_size+1) ;

hold on
geoshow(gca, so_fronts.lat_pf, so_fronts.lon_pf, 'displaytype', 'line', 'linewidth', boundary_width, 'color', 'k') % polar front
geoshow(gca, so_fronts.lat_saf, so_fronts.lon_saf, 'displaytype', 'line', 'linewidth', boundary_width, 'color', 'k') % Subantarctic Front

% plot longitude boundary between pacific and indian
geoshow(coast.lat,coast.long,'DisplayType','polygon', 'facecolor', coast_gray)
geoshow([-90 -30], [geo_lims.pac_lims(1) geo_lims.pac_lims(1)], 'displaytype', 'line', 'color', 'k', 'linestyle', '--', 'linewidth', 2)


clear wint_index
pco2 = 0; % 0 = map, 1 = delta pCO2, 2 = time

ph(3) = subplot(3,1,3); % wintertime observations w/in SAMW density range
hold on
% coast = load('coast');
axesm('lambertstd', 'MapLatLimit',[-80 -37],'MapParallels',[-75 -15],'MapLonLimit',[68 -70])
axis off; framem on; gridm on; mlabel on; plabel on;
setm(gca,'MLabelParallel',-20, 'fontsize', 14)

hold on
coast_gray = [.8 .8 .8];

% Criteria for plot SAMW density range, winter months, MLD cutoff, after
% 1990
for r=1:2
    
    socat_index = store_socat.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(1) &  ...
        store_socat.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) & ...
        store_socat.(regions{r}).month>=month_lims_for_mean(1) & store_socat.(regions{r}).month<=month_lims_for_mean(2) & store_socat.(regions{r}).MLD>=MLD_cutoff & ...
        store_socat.(regions{r}).year>1990;
    
    
    gdap_index = store_gdap.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(1) &  ...
        store_gdap.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) & ...
        store_gdap.(regions{r}).month>=month_lims_for_mean(1) & store_gdap.(regions{r}).month<=month_lims_for_mean(2) & store_gdap.(regions{r}).MLD>=MLD_cutoff &...
        store_gdap.(regions{r}).year>1990;
    
    
    %
    argo_index = store_argo.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(1) &  ...
        store_argo.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) & ...
        store_argo.(regions{r}).month>=month_lims_for_mean(1) & store_argo.(regions{r}).month<=month_lims_for_mean(2) & store_argo.(regions{r}).MLD>=MLD_cutoff;% & ...
    %         ;
    
    if pco2==1
        scatterm(store_socat.(regions{r}).Lat(socat_index), store_socat.(regions{r}).Lon(socat_index), 20, store_socat.(regions{r}).delta_pCO2(socat_index), 'filled')
        
        scatterm(store_gdap.(regions{r}).Lat(gdap_index), store_gdap.(regions{r}).Lon(gdap_index), 20, store_gdap.(regions{r}).delta_pCO2(gdap_index), 'filled')
        
        scatterm(store_argo.(regions{r}).Lat(argo_index), store_argo.(regions{r}).Lon(argo_index), 20, store_argo.(regions{r}).delta_pCO2(argo_index), 'filled', 'marker', 's')
        
    elseif pco2==2
        
        scatterm(store_socat.(regions{r}).Lat(socat_index), store_socat.(regions{r}).Lon(socat_index), 20, store_socat.(regions{r}).year(socat_index) + store_socat.(regions{r}).month(socat_index)./12, 'filled')
        
        scatterm(store_gdap.(regions{r}).Lat(gdap_index), store_gdap.(regions{r}).Lon(gdap_index), 20, store_gdap.(regions{r}).year(gdap_index)+ store_gdap.(regions{r}).month(gdap_index)./12, 'filled')
        
        scatterm(store_argo.(regions{r}).Lat(argo_index), store_argo.(regions{r}).Lon(argo_index), 20, store_argo.(regions{r}).year(argo_index)+ store_argo.(regions{r}).month(argo_index)./12, 'filled', 'marker', 's')
        
        
        
    else
        % for general map, don't filter based on pCO2
        argo_index = store_argo.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(1) &  ...
            store_argo.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) & ...
            store_argo.(regions{r}).month>=month_lims_for_mean(1) & store_argo.(regions{r}).month<=month_lims_for_mean(2) & store_argo.(regions{r}).MLD>=MLD_cutoff;
        
        s_h = geoshow(store_socat.(regions{r}).Lat(socat_index), store_socat.(regions{r}).Lon(socat_index), ...
            'displaytype', 'multipoint', 'marker', 'd', 'markerfacecolor', colormap2(1,:), 'markeredgecolor', colormap2(4,:), 'markerfacecolor', colormap2(4,:) ,'markersize', marker_size+1);
        
        
        a_h = geoshow(store_argo.(regions{r}).Lat(argo_index), store_argo.(regions{r}).Lon(argo_index), 'displaytype', 'multipoint', 'marker', 'o', 'markeredgecolor', colormap2(2,:), 'markerfacecolor', colormap2(2,:),'markersize', marker_size) ;
        
        g_h = geoshow(store_gdap.(regions{r}).Lat(gdap_index), store_gdap.(regions{r}).Lon(gdap_index), ...
            'displaytype', 'multipoint', 'marker', 's', 'markerfacecolor', colormap2(3,:), 'markeredgecolor', colormap2(10,:), 'markerfacecolor', colormap2(10,:) ,'markersize', marker_size+1) ;
        
    end
    
end
if pco2==1
    set(gcf, 'colormap2', flipud(brewermap(30, 'RdYlBu')))
    c1 = colorbar;
    
    caxis([-50 50])
    ylabel(c1, '\DeltapCO_2')
elseif pco2==2
    c_map = (brewermap(30, 'GnBu'));
    
    set(gcf, 'colormap2', c_map(1:end,:))
    c1 = colorbar;
    
    caxis([1990 2020])
    ylabel(c1, 'Year')
    
    
end

hold on
geoshow(gca, so_fronts.lat_pf, so_fronts.lon_pf, 'displaytype', 'line', 'linewidth', boundary_width, 'color', 'k') % polar front
geoshow(gca, so_fronts.lat_saf, so_fronts.lon_saf, 'displaytype', 'line', 'linewidth', boundary_width, 'color', 'k') % Subantarctic Front

% plot winter 200 m MLD contour
contourm(double(argo_MLD.LAT), double(argo_MLD.LON), argo_MLD.wint_mean, 'levellist', 200, 'color', [.6 .6 .6],'linestyle', '-', 'linewidth', 2)

geoshow(coast.lat,coast.long,'DisplayType','polygon', 'facecolor', coast_gray)
geoshow([-90 -30], [geo_lims.pac_lims(1) geo_lims.pac_lims(1)], 'displaytype', 'line', 'color', 'k', 'linestyle', '--', 'linewidth', 2)

setm(gca, 'mlabellocation', [0 ], 'mlabelparallel', 0)
setm(gca, 'plabellocation', [-180 -180], 'plabelmeridian', -150)

all_pos = NaN(1,4);
for q = 1:3
    all_pos(q,:) = get(ph(q), 'position');
end

m_scale = 0.10;
for q = 1:3
    set(ph(q), 'position', all_pos(q,:)+[-0.1 -0.1 m_scale m_scale])
end

clear temp_prop temp_lat temp_lon


top_pos = get(ph(1), 'position');
set(ph(1), 'position', top_pos+[0 0 0 0])

c1 = colorbar(ph(1), 'location', 'eastoutside');
ylabel(c1, 'Mean MLD (m)', 'interpreter', 'none', 'fontsize', 13)
set(ph(1), 'position', top_pos+[0 0 0 0])

print(gcf, '-dpdf', '-r800', [plot_dir plot_filename '_w_colorbar.pdf' ])
% print(gcf, '-dpng', '-r400', [plot_dir plot_filename '_w_colorbar.png' ])

clear wint_index r s_h g_h a_h argo_vec top_pos ph q R s socat_index socat_wint_index m_scale gdap_wint_index f c1
clear all_pos argo_index colormap colormap2 data1 gdap_index lat_map lon_map marker_size pco2 plot_filename 

%% Figure 2
% Use existing lat/lon bounds but increase the northward boundary of latitude to show re-ventilation sites:
font_size = 15;

r=1;

n_lon_lim = -30;

pac_SNs = {};

for f = 1:length(SO_SNs)
    if r==1
        geo_index = (Argo.(SO_SNs{f}).Lat>geo_lims.lat_lims_pac(1) & Argo.(SO_SNs{f}).Lat<n_lon_lim & ...
            Argo.(SO_SNs{f}).Lon>=geo_lims.pac_lims(1) & Argo.(SO_SNs{f}).Lon<geo_lims.pac_lims(2));
    else
        geo_index = (Argo.(SO_SNs{f}).Lat>geo_lims.lat_lims_ind(1) & Argo.(SO_SNs{f}).Lat<n_lon_lim & ...
            Argo.(SO_SNs{f}).Lon>=geo_lims.ind_lims(1) & Argo.(SO_SNs{f}).Lon<geo_lims.ind_lims(2));
    end
    if sum(geo_index)>0
        pac_SNs{end+1,1} = SO_SNs{f};
    end
end


% using float 5904695 as an example section
f = 32;
if r==1
    plot_filename = ['Pacific_' pac_SNs{f}];
else
    plot_filename = ['Indian_' pac_SNs{f}];
    
end
plot_filename = ['Fig_2_' plot_filename];

% map to go w float section plot: Fig 2a
set(gcf, 'colormap', plasma)

clf
set(gcf, 'units', 'inches')
paper_w = 5; paper_h =4;
set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]); clear paper_w paper_h

subplot(1,1,1)
hold on
coast = load('coast');
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -37])
axis off; framem on; gridm on; mlabel on; plabel on;
setm(gca,'MLabelParallel',-37, 'PLabelMeridian', 0, 'fontsize', font_size-5)

hold on
date_vec = datevec(Argo.(pac_SNs{f}).GMT_Matlab);
%     scatter(Argo.(pac_SNs{f}).Lon, Argo.(pac_SNs{f}).Lat,50, date_vec(:,1), 'filled')
scatterm(Argo.(pac_SNs{f}).Lat,Argo.(pac_SNs{f}).Lon, 80, date_vec(:,1)+ date_vec(:,2)./12, 'filled')
geoshow(coast.lat,coast.long,'DisplayType','polygon', 'facecolor', coast_gray)

pos_map = get(gca, 'position');
geoshow(gca, so_fronts.lat_saf, so_fronts.lon_saf, 'displaytype', 'line', 'linewidth', 1.5, 'color', 'k')% Subantarctic Front
geoshow(gca, so_fronts.lat_pf, so_fronts.lon_pf, 'displaytype', 'line', 'linewidth', 1.5, 'color', 'k')% Polar Front

if r==1
    geoshow([-90 -90 -30 -30 -90],[geo_lims.pac_lims(1) geo_lims.pac_lims(2) geo_lims.pac_lims(2) geo_lims.pac_lims(1) geo_lims.pac_lims(1)], 'displaytype', 'line', 'linewidth', 1.5, 'color', 'r')
else
    geoshow([-90 -90 -30 -30 -90],[geo_lims.ind_lims(1) geo_lims.ind_lims(2) geo_lims.ind_lims(2) geo_lims.ind_lims(1) geo_lims.ind_lims(1)], 'displaytype', 'line', 'linewidth', 1.5, 'color', 'b')
end

c1 = colorbar;
ylabel(c1, 'Year')
pos_cbar = get(c1, 'position');

%
set(gca, 'position', (pos_map+[-.1 -.04 0 0]).*[1 1 1 1]);

set(c1, 'position', pos_cbar+[.005 -.04 .013 -.05], 'fontsize', 13.5, 'ticks', 2016:2021);
set(gca, 'fontsize', font_size)

print(gcf, '-dpdf', '-r800', [plot_dir plot_filename '_map_R1.pdf' ])

clear c1 pos_cbar pos_map
%% plot subplot b
coast_gray = [.8 .8 .8];
set(gcf, 'colormap', turbo)
date_grid = repmat(Argo.(pac_SNs{f}).GMT_Matlab, 1, size(Argo.(pac_SNs{f}).Press_db,2));

clf
set(gcf, 'units', 'inches')
paper_w = 10; paper_h =6;
set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]); clear paper_w paper_h

subplot(1,1,1)

contourf(date_grid', Argo.(pac_SNs{f}).Press_db', Argo.(pac_SNs{f}).o2_umol_kg', 'levellist', [180:5:350], 'linestyle', 'none'); set(gca, 'ylim', [0 600])

hold on
% add MLD, SAMW density boundaries to the plot
set(gca, 'ydir', 'reverse', 'xtick', datenum(2017:2021, 1, 1));
plot(Argo.(pac_SNs{f}).GMT_Matlab, Argo.(pac_SNs{f}).MLD, 'color', [.5 .5 .5], 'linewidth', 3);

if r==1
    contour(date_grid', Argo.(pac_SNs{f}).Press_db', Argo.(pac_SNs{f}).PDENS'-1000, [26.8 26.8], 'linestyle', '-', 'color','b', 'linewidth', 3)
    contour(date_grid', Argo.(pac_SNs{f}).Press_db', Argo.(pac_SNs{f}).PDENS'-1000, [27.05 27.05], 'linestyle', '-', 'color', 'k', 'linewidth', 3)
else
    contour(date_grid', Argo.(pac_SNs{f}).Press_db', Argo.(pac_SNs{f}).PDENS'-1000, [26.6 26.6], 'linestyle', '-', 'color','b', 'linewidth', 3)
    contour(date_grid', Argo.(pac_SNs{f}).Press_db', Argo.(pac_SNs{f}).PDENS'-1000, [26.85 26.85], 'linestyle', '-', 'color', 'k', 'linewidth', 3)
end


c1 = colorbar;
ylabel(c1, '[O_2] (\mumol kg^-^1)', 'fontsize', font_size)
ylabel('Pressure (db)')
caxis([220 320])
if r==1
    %         title([num2str(f) ' ' pac_SNs{f} ' - Pacific focus'])
    title([ 'WMO: ' pac_SNs{f}(2:end)])
    
else
    title([num2str(f) ' ' pac_SNs{f} ' - Indian focus'])
    
end
% show when the floats are in different regions
if r==1
    pac_index =  Argo.(pac_SNs{f}).Lon>=geo_lims.pac_lims(1) & Argo.(pac_SNs{f}).Lon<geo_lims.pac_lims(2);
    
    temp_vec = ones(length(pac_index),1);
    plot(Argo.(pac_SNs{f}).GMT_Matlab(pac_index), temp_vec(pac_index), 'or-', 'markerface', 'r');
else
    pac_index =  Argo.(pac_SNs{f}).Lon>=geo_lims.ind_lims(1) & Argo.(pac_SNs{f}).Lon<geo_lims.ind_lims(2);
    
    temp_vec = ones(length(pac_index),1);
    plot(Argo.(pac_SNs{f}).GMT_Matlab(pac_index), temp_vec(pac_index), 'ob-', 'markerface', 'b');
    
end

datetick('x', 'YYYY-mm', 'keepticks')
set(gca, 'fontsize', font_size)


%     print(gcf, '-dpng', [plot_dir plot_filename '.png' ])
print(gcf, '-dpdf', '-r800', [plot_dir plot_filename '_R1.pdf' ])

% pause
clear temp_o2 r pos_cbar pos_map pac_SNs pac_index n_lon_lim geo_index f date_vec date_grid c1 temp_vec plot_filename


%% figure 3 - SAMW formation properties and volume vs. density
% also save out data for Table 1
% comparison between MW Volume and BGC/T/S props vs. density

carter_overlay = 1; % overlay points from Carter et al. 2014

plot_filename = ['Fig_3A_Mode Water formation properties- weighted ' ...
    ' Months ' num2str(month_lims_for_mean(1)) ' to ' num2str(month_lims_for_mean(2)) ' Dens. int for mean ' num2str(dens_interval) ' MLD cut ' num2str(MLD_cutoff) '_v16_.03MLD'];


% temp value to make the plots until I have the 200m data:
MLD_cutoff_temp = 200; %disp('Replace with new bSOSE data');
% sose.indian = load([cerovecki_dir 'Data_from_Ivana/2021_04_08/MLD_BIN_average_1_511_Ind_' num2str(MLD_cutoff_temp) '_m_05_ST_winter_iter135_lat_30_55_68_170.mat']);
% sose.pacific = load([cerovecki_dir 'Data_from_Ivana/2021_04_08/MLD_BIN_average_1_511_Pac_' num2str(MLD_cutoff_temp) '_m_05_ST_winter_iter135_lat_45_64_170_290']);

sose.indian = load([cerovecki_dir 'Data_from_Ivana/2021_08/BIN_MLD_CLIM_iter135_Ind_200_m_05_ST_lat_30_55_lon_68_170_regional_indian.mat']);
sose.pacific = load([cerovecki_dir 'Data_from_Ivana/2021_08/BIN_MLD_CLIM_iter135_Pac_200_m_05_ST_lat_45_64_lon_170_290_regional_pacific']);


out_pac_ind_mean = NaN(2,6);
out_pac_ind_std = NaN(2,6);
plot_colors = brewermap(10, 'Paired');

prop_lims = [240 340;...
    -20 20;...
    2090 2200;...
    0 35;...
    365 480;...
    0 14;...
    33.5 36];

prop_names = {'[O_2] (\mumol kg^-^1)'; ...
    ' '; ...
    '[DIC] (\mumol kg^-^1)' ; ...
    '[NO_3^-] (\mumol kg^-^1)'; ...
    '{\itp}CO_2 (\muatm)'; ...
    '\theta (\circC)'; ...
    'Sal. (PSS-78)'} ;

plot_map = [1 4 3 5];

clf

set(gcf, 'units', 'inches')
paper_w = 12; paper_h =12;
set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]); clear paper_w paper_h

d = NaN(8,1);
clear fig_3_data

for r=1:2
    
    axis_font_size = 12;
    
    
    for pl = 1:4%length(properties)
        if r ==1
            plot_number = pl*2-r;
        else
            plot_number = pl*2;
        end
        d(plot_number) = subplot(4,2,plot_number);
        hold on
        
        ylabel('MW Vol (m^3)', 'fontsize', axis_font_size)
        if pl<4
            xlabel('\sigma_\theta', 'fontsize', axis_font_size)
        else
            xlabel('\sigma_\theta (kg m^-^3)', 'fontsize', axis_font_size)
        end
        grid on
    end
    
    for pl = 1:length(plot_map)
        if r ==1
            plot_number = pl*2-r;
        else
            plot_number = pl*2;
        end
        
        axes(d(plot_number))
        
        
        fig_3_data.vol_sig_th = mode_vol.(mode_years{1}).sig_th';
        fig_3_data.vol = mw_prop.(regions{r}).wint.vol;
                
        dens_index = find(mw_prop.(regions{r}).all.per_5_index);
        
        intermediate_value_first = (mw_prop.(regions{r}).wint.vol(dens_index(1)) + mw_prop.(regions{r}).wint.vol(dens_index(1)-1))./2;
        intermediate_value_last = (mw_prop.(regions{r}).wint.vol(dens_index(end)) + mw_prop.(regions{r}).wint.vol(dens_index(end)+1))./2;
        
        p1 = patch([mw_prop.(regions{r}).all.sig_thet_range(1) mw_prop.(regions{r}).all.sig_thet_range(1) mw_prop.(regions{r}).all.sig_thet_range(1:end-1)+.025  ...
            mw_prop.(regions{r}).all.sig_thet_range(end) mw_prop.(regions{r}).all.sig_thet_range(end)], ...
            [0 intermediate_value_first mw_prop.(regions{r}).wint.vol(mw_prop.(regions{r}).all.per_5_index,1)' intermediate_value_last 0], [.9 .9 .9]);
        p1.LineStyle = 'none';
        
        errorbar(mode_vol.(mode_years{1}).sig_th+.025, mw_prop.(regions{r}).wint.vol(:,1), mw_prop.(regions{r}).wint.vol(:,2), 'color', plot_colors(2,:), 'linestyle', '-', 'linewidth', 2)
        
       
        set(gca, 'ylim', [0 17e14])
    end
    
    if r==1 && carter_overlay==1
        % overlay carter et al. 2014 samples on these plots
        axes(d(1))
        yyaxis right
        plot(26.997, 297.8, 'x', 'color', [.5 .5 .5], 'markersize', 8, 'linewidth', 2)
        plot(27.021, 296.8, 'o', 'color', [.5 .5 .5], 'markersize', 8, 'linewidth', 2)
        
        axes(d(3)); yyaxis right
        plot(26.997, 20.32, 'x', 'color', [.5 .5 .5], 'markersize', 8, 'linewidth', 2)
        plot(27.021, 21.93, 'o', 'color', [.5 .5 .5], 'markersize', 8, 'linewidth', 2)
        
        axes(d(5)); yyaxis right
        plot(26.997, 2122.2, 'x', 'color', [.5 .5 .5], 'markersize', 8, 'linewidth', 2)
        plot(27.021, 2129, 'o', 'color', [.5 .5 .5], 'markersize', 8, 'linewidth', 2)
        
        plot(26.997, 2122.2+carter.dic_change, 'kx', 'markerfacecolor', 'k', 'markersize', 8, 'linewidth', 2)
        plot(27.021, 2129+carter.dic_change, 'ko', 'markersize', 8, 'linewidth', 2)
        
        axes(d(7)); yyaxis right
        plot(26.997, carter.carter_pCO2_orig(1), 'x', 'color', [.5 .5 .5], 'markersize', 8, 'linewidth', 2)
        plot(27.021, carter.carter_pCO2_orig(2) , 'o', 'color', [.5 .5 .5], 'markersize', 8, 'linewidth', 2)
        
        plot(26.997, carter.carter_pCO2_new(1), 'kx', 'markerfacecolor', 'k', 'markersize', 8, 'linewidth', 2)
        plot(27.021, carter.carter_pCO2_new(2), 'ko', 'markersize', 8, 'linewidth', 2)
        
    end
    
    for pl = 1:length(plot_map)
        if r ==1
            plot_number = pl*2-r;
        else
            plot_number = pl*2;
        end
        
        axes(d(plot_number))
        
        p = plot_map(pl);
        if p==1
            p_color = plot_colors(10,:);
            b_sose_color = plot_colors(9,:);
            sose_prop = 'o2_mld_clim';
        elseif p==5
            p_color = plot_colors(4,:);
            b_sose_color = plot_colors(3,:);
            
            sose_prop = 'pco2_mld_clim';
        elseif p==3
            p_color = plot_colors(6,:);
            b_sose_color = plot_colors(5,:);
            
            sose_prop = 'dic_mld_clim';
        elseif p==4
            p_color = plot_colors(8,:);
            b_sose_color = plot_colors(7,:);
            sose_prop = 'no3_mld_clim';
        elseif p==6
            p_color = plot_colors(8,:);
            sose_prop = 'T_mld_clim';
        elseif p==7
            p_color = plot_colors(p,:);
            sose_prop = 'S_mld_clim';
        elseif p==2
            p_color = plot_colors(p,:);
            sose_prop = [];
        end
        
        yyaxis right
        hold on
        errorbar( mode_vol.(mode_years{1}).sig_th+.025, squeeze(mw_prop.(regions{r}).wint.(properties{p})(:,1)), squeeze(mw_prop.(regions{r}).wint.(properties{p})(:,2)), '-', 'linewidth', 2, 'color', p_color)
        
        fig_3_data.prop_sig_th = mode_vol.(mode_years{1}).sig_th';
        
        fig_3_data.(regions{r}).(properties{p}) = mw_prop.(regions{r}).wint.(properties{p});
        
        if ~isempty(sose_prop)
            if r==1
                  plot(nanmean(sose.(regions{r}).R_mld_clim(5:end,month_lims_for_mean(1):month_lims_for_mean(end)),2), ...
                    nanmean(sose.(regions{r}).(sose_prop)(5:end,month_lims_for_mean(1):month_lims_for_mean(end)),2), 'color', b_sose_color, 'linestyle', '--', 'linewidth', 3, 'marker', 'none')
                 
            else
                  plot(nanmean(sose.(regions{r}).R_mld_clim(3:17,month_lims_for_mean(1):month_lims_for_mean(end)),2), ...
                    nanmean(sose.(regions{r}).(sose_prop)(3:17,month_lims_for_mean(1):month_lims_for_mean(end)),2), 'color', b_sose_color, 'linestyle', '--', 'linewidth', 3, 'marker', 'none')
            end
        end
        % title contains the mean properties across the month range specified
        % earlier - should include in text somewhere (first title)
        if pl==1
%             title(upper(regions{r}))
            if r==1
                title('Pacific')
            elseif r==2
                title('Indian')
            end
        else
            
        end
        
        out_pac_ind_mean(r,pl+2) = mw_prop.(regions{r}).all.mean_prop_weighted.(properties{p})(1);
        out_pac_ind_std(r,pl+2) = mw_prop.(regions{r}).all.mean_prop_weighted.(properties{p})(2);
        
        %         ylabel(properties{p}, 'interpreter', 'none')
        
        ylabel(prop_names{p})
        set(gca, 'ylim', prop_lims(p,:))
        set(gca, 'ycolor', p_color)
        %         if r==1
        set(gca, 'xlim', [26.4 27.3]);
        %         else
        %             set(gca, 'xlim', [26.4 27.2]);
        %         end
        yyaxis left
        set(gca, 'ycolor', plot_colors(2,:))
        
        %     set(gca, 'ylim', [0 11e14])
    end
    %
    
    
    if strncmp((regions{r}), 'pacific', 7)
        lon_lims = geo_lims.pac_lims;
        lat_lims = geo_lims.lat_lims_pac;
    elseif strncmp((regions{r}), 'indian',6)
        lon_lims = geo_lims.ind_lims;
        lat_lims = geo_lims.lat_lims_ind;
    end
    
    
end

% text(27.9,8e14, ['Lat lims: ' num2str(lat_lims) ' Lon lims: ' num2str(lon_lims) ' weighted by vol.'], 'fontsize', 15)
% text(27.9,2e14, [' Months: ' num2str(month_lims_for_mean(1)) ' to ' num2str(month_lims_for_mean(2)) ' Dens. int for mean: ' num2str(dens_interval) ' MLD cut ' num2str(MLD_cutoff)], 'fontsize', 15)

print(gcf, '-dpdf', '-r800', [plot_dir plot_filename '.pdf' ])
% print(gcf, '-dpng', '-r500', [plot_dir plot_filename '.png' ])

%% figure 3B - SAMW formation T, S and volume vs. density

plot_filename = ['Fig_3B_Mode Water formation properties- weighted ' ...
    ' Months ' num2str(month_lims_for_mean(1)) ' to ' num2str(month_lims_for_mean(2)) ' Dens. int for mean ' num2str(dens_interval) ' MLD cut ' num2str(MLD_cutoff) '_v16_.03MLD'];

% comparison to Ivana's values
% temp value to make the plots until I have the 200m data:
% MLD_cutoff_temp = 200; %disp('Replace with new bSOSE data');

% sose.indian = load([cerovecki_dir 'Data_from_Ivana/2021_04_08/' 'MLD_BIN_average_1_511_Ind_' num2str(MLD_cutoff_temp) '_m_05_ST_winter_iter135_lat_30_55_68_170.mat']);
% sose.pacific = load([cerovecki_dir 'Data_from_Ivana/2021_04_08/' 'MLD_BIN_average_1_511_Pac_' num2str(MLD_cutoff_temp) '_m_05_ST_winter_iter135_lat_45_64_170_290']);

plot_colors = brewermap(9, 'Set1');
bsose_color_sal = brewermap(10, 'Set2');
bsose_color_temp = brewermap(10, 'Set3');

plot_map = [6 7];

clf
set(gcf, 'units', 'inches')
paper_w = 12; paper_h =12;
set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]); clear paper_w paper_h


d = NaN(8,1);

for r=1:2
    
    axis_font_size = 12;
    
    
    
    for pl = 1:2%length(properties)
        if r ==1
            plot_number = pl*2-r;
        else
            plot_number = pl*2;
        end
        d(plot_number) = subplot(4,2,plot_number);
        hold on
        
        ylabel('MW Vol (m^3)', 'fontsize', axis_font_size)
        if pl==1
            xlabel('\sigma_\theta', 'fontsize', axis_font_size)
        elseif pl==2
            xlabel('\sigma_\theta (kg m^-^3)', 'fontsize', axis_font_size)
        end
        grid on
    end
    
    for pl = 1:length(plot_map)
        if r ==1
            plot_number = pl*2-r;
        else
            plot_number = pl*2;
        end
        
        axes(d(plot_number))
       
        dens_index = find(mw_prop.(regions{r}).all.per_5_index);
        
        intermediate_value_first = (mw_prop.(regions{r}).wint.vol(dens_index(1)) + mw_prop.(regions{r}).wint.vol(dens_index(1)-1))./2;
        intermediate_value_last = (mw_prop.(regions{r}).wint.vol(dens_index(end)) + mw_prop.(regions{r}).wint.vol(dens_index(end)+1))./2;
        
        
        
        p1 = patch([mw_prop.(regions{r}).all.sig_thet_range(1) mw_prop.(regions{r}).all.sig_thet_range(1) mw_prop.(regions{r}).all.sig_thet_range(1:end-1)+.025  ...
            mw_prop.(regions{r}).all.sig_thet_range(end) mw_prop.(regions{r}).all.sig_thet_range(end)], ...
            [0 intermediate_value_first mw_prop.(regions{r}).wint.vol(mw_prop.(regions{r}).all.per_5_index,1)' intermediate_value_last 0], [.9 .9 .9]);
        p1.LineStyle = 'none';
        
        
        errorbar(mode_vol.(mode_years{1}).sig_th+.025, mw_prop.(regions{r}).wint.vol(:,1), mw_prop.(regions{r}).wint.vol(:,2), 'color', plot_colors(2,:), ...
            'linestyle', '-', 'linewidth', 2)
       
        set(gca, 'ylim', [0 17e14])
    end
    
    if r==1 && carter_overlay==1
        % overlay carter et al. 2014 samples on these plots
        axes(d(1))
        yyaxis right
        plot(26.997, 5.247, 'x', 'color', [.5 .5 .5], 'markersize', 8, 'linewidth', 2)
        plot(27.021, 4.68, 'o', 'color', [.5 .5 .5], 'markersize', 8, 'linewidth', 2)
        
        axes(d(3)); yyaxis right
        plot(26.997, 34.181, 'x', 'color', [.5 .5 .5], 'markersize', 8, 'linewidth', 2)
        plot(27.021, 34.129, 'o', 'color', [.5 .5 .5], 'markersize', 8, 'linewidth', 2)
    end
    
    for pl = 1:length(plot_map)
        if r ==1
            plot_number = pl*2-r;
        else
            plot_number = pl*2;
        end
        
        axes(d(plot_number))
        
        p = plot_map(pl);
        if p==1
            p_color = plot_colors(10,:);
            b_sose_color = plot_colors(9,:);
            sose_prop = 'o2_win';
        elseif p==5
            p_color = plot_colors(4,:);
            b_sose_color = plot_colors(3,:);
            
            sose_prop = 'pco2_win';
        elseif p==3
            p_color = plot_colors(6,:);
            b_sose_color = plot_colors(5,:);
            
            sose_prop = 'dic_win';
        elseif p==4
            p_color = plot_colors(8,:);
            b_sose_color = plot_colors(7,:);
            sose_prop = 'no3_win';
        elseif p==6
            p_color = plot_colors(8,:);
            b_sose_color = bsose_color_temp(8,:);
            
            sose_prop = 'T_mld_clim';
        elseif p==7
            p_color = plot_colors(p,:);
            b_sose_color = bsose_color_sal(7,:);
            
            sose_prop = 'S_mld_clim';
        elseif p==2
            p_color = plot_colors(p,:);
            sose_prop = [];
        end
        
        yyaxis right
        hold on
        errorbar( mode_vol.(mode_years{1}).sig_th+.025, squeeze(mw_prop.(regions{r}).wint.(properties{p})(:,1)), squeeze(mw_prop.(regions{r}).wint.(properties{p})(:,2)), '-',...
            'linewidth', 2, 'color', p_color)
        
        fig_3_data.(regions{r}).(properties{p}) = mw_prop.(regions{r}).wint.(properties{p});
        
        if ~isempty(sose_prop)
            if r==1

                plot(nanmean(sose.(regions{r}).R_mld_clim(5:end,month_lims_for_mean(1):month_lims_for_mean(end)),2), ...
                    nanmean(sose.(regions{r}).(sose_prop)(5:end,month_lims_for_mean(1):month_lims_for_mean(end)),2), 'color', b_sose_color, 'linestyle', '--', 'linewidth', 3)
                
            else
                plot(nanmean(sose.(regions{r}).R_mld_clim(3:17,month_lims_for_mean(1):month_lims_for_mean(end)),2), ...
                    nanmean(sose.(regions{r}).(sose_prop)(3:17,month_lims_for_mean(1):month_lims_for_mean(end)),2), 'color', b_sose_color, 'linestyle', '--', 'linewidth', 3)
                
            end
        end
        % title contains the mean properties across the month range specified
        % earlier - should include in text somewhere (first title)
        if pl==1
            %             title(upper(regions{r}))
            if r==1
                title('Pacific')
            elseif r==2
                title('Indian')
            end
        else
            
        end
        
        out_pac_ind_mean(r,pl) = mw_prop.(regions{r}).all.mean_prop_weighted.(properties{p})(1);
        out_pac_ind_std(r,pl) = mw_prop.(regions{r}).all.mean_prop_weighted.(properties{p})(2);
        
        
        ylabel(prop_names{p})
        
        set(gca, 'ylim', prop_lims(p,:))
        set(gca, 'ycolor', p_color)
        %         if r==1
        set(gca, 'xlim', [26.4 27.3]);
        %         else
        %             set(gca, 'xlim', [26.4 27.2]);
        %         end
        yyaxis left
        set(gca, 'ycolor', plot_colors(2,:))
        
        %     set(gca, 'ylim', [0 11e14])
    end
    %
    
    
    if strncmp((regions{r}), 'pacific', 7)
        lon_lims = geo_lims.pac_lims;
        lat_lims = geo_lims.lat_lims_pac;
    elseif strncmp((regions{r}), 'indian',6)
        lon_lims = geo_lims.ind_lims;
        lat_lims = geo_lims.lat_lims_ind;
    end
    
    
end


save([plot_dir plot_filename '.mat' ], 'fig_3_data')

print(gcf, '-dpdf', '-r800', [plot_dir plot_filename '.pdf' ])

clear p_color p p1 pl r sose_prop plot_number plot_map plot_colors intermediate_value_first intermediate_value_last
clear dens_index d carter_overlay b_sose_color
%% Create monthly means using the same filtering criteria as the rest of figure 4 - Run prior to Fig. 4
clear temp_filt

data_sources = {'argo', 'gdap', 'socat'};

year_range = 1990:2020;
temp_filt.GMT_Matlab = NaN(length(year_range).*12,1);

props_out = [2 10 6 7 9 13 14 1 15];

properties = {'ml_o2' 'ml_delta_o2' 'ml_DIC' 'ml_no3'  'ml_pCO2', 'ml_PTMP', 'ml_Sal', 'MLD', 'ml_pdens', 'delta_pCO2_unified', 'Oxygen_flux', 'CO2_flux', 'Lat', 'Lon', 'SLP'};


for r=1:2
    
    for p=props_out
        
        for ds = 1:length(data_sources)
            temp_filt.(regions{r}).(data_sources{ds}).(properties{p}) = NaN(length(year_range).*12,2);
            temp_filt.(regions{r}).(data_sources{ds}).(properties{p}) = NaN(length(year_range).*12,2);
        end
    end
end

clear r p ds
SNs_out = [];
for yy = 1:length(year_range)
    for mon = 1:12
        temp_filt.GMT_Matlab((yy-1)*12+mon) = datenum(year_range(yy), mon, 15);
        
        for p=props_out
            for r=1:2
                %                 sig_thet_range = mode_vol.y2005.sig_th(mw_prop.(regions{r}).all.per_5_index);
                
                if p==10 % only filter for delta pCO2 if you are looking at that (make sure this doesn't change any calculations...) - changed to compare ml ptmp to ml o2 and possibly add to supplemental figure of O2 vs. T IAV
                    argo_index = store_argo.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(1) &  ...
                        store_argo.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) & ...
                        store_argo.(regions{r}).month>=month_lims_for_mean(1) & store_argo.(regions{r}).month<=month_lims_for_mean(2) & store_argo.(regions{r}).MLD>=MLD_cutoff & ...
                        store_argo.(regions{r}).month==mon & store_argo.(regions{r}).year==year_range(yy) ...
                        & ~isnan(store_argo.(regions{r}).delta_pCO2_unified);
                else
                    
                    %                     p==2 ||  p==1 % don't filter for delta pCO2 if you are looking at delta O2 or ml o2, otherwise do
                    argo_index = store_argo.(regions{r}).ml_pdens-1000>mw_prop.(regions{r}).all.sig_thet_range(1) &  ...
                        store_argo.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) & ...
                        store_argo.(regions{r}).month>=month_lims_for_mean(1) & store_argo.(regions{r}).month<=month_lims_for_mean(2) & store_argo.(regions{r}).MLD>=MLD_cutoff & ...
                        store_argo.(regions{r}).month==mon & store_argo.(regions{r}).year==year_range(yy);
                    
                end
                temp_filt.(regions{r}).argo.(properties{p})((yy-1)*12+mon,1) = nanmean(store_argo.(regions{r}).(properties{p})(argo_index));
                temp_filt.(regions{r}).argo.(properties{p})((yy-1)*12+mon,2) = nanstd(store_argo.(regions{r}).(properties{p})(argo_index));
                
                if p==2
                    SNs_out = [SNs_out; store_argo.(regions{r}).float_SN(argo_index)]; % list of SNs used in delta O2 calculation 
                end
                
                if p~=2 && p~=1 && p~=15
                    socat_index = store_socat.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(1) &  ...
                        store_socat.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) & ...
                        store_socat.(regions{r}).month>=month_lims_for_mean(1) & store_socat.(regions{r}).month<=month_lims_for_mean(2) & store_socat.(regions{r}).MLD>=MLD_cutoff & ...
                        store_socat.(regions{r}).month==mon & store_socat.(regions{r}).year==year_range(yy)...
                        & ~isnan(store_socat.(regions{r}).delta_pCO2_unified);
                    
                    
                    temp_filt.(regions{r}).socat.(properties{p})((yy-1)*12+mon,1) = nanmean(store_socat.(regions{r}).(properties{p})(socat_index));
                    temp_filt.(regions{r}).socat.(properties{p})((yy-1)*12+mon,2) = nanstd(store_socat.(regions{r}).(properties{p})(socat_index));
                    
                end
                
                if p~=15
                    
                if p==2 ||  p==1
                    gdap_index = store_gdap.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(1) &  ...
                        store_gdap.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) & ...
                        store_gdap.(regions{r}).month>=month_lims_for_mean(1) & store_gdap.(regions{r}).month<=month_lims_for_mean(2) & store_gdap.(regions{r}).MLD>=MLD_cutoff & ...
                        store_gdap.(regions{r}).month==mon & store_gdap.(regions{r}).year==year_range(yy);
                    
                    temp_filt.(regions{r}).gdap.(properties{p})((yy-1)*12+mon,1) = nanmean(store_gdap.(regions{r}).(properties{p})(gdap_index));
                    temp_filt.(regions{r}).gdap.(properties{p})((yy-1)*12+mon,2) = nanstd(store_gdap.(regions{r}).(properties{p})(gdap_index));
                    
                else
                    gdap_index = store_gdap.surf.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(1) &  ...
                        store_gdap.surf.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) & ...
                        store_gdap.surf.(regions{r}).month>=month_lims_for_mean(1) & store_gdap.surf.(regions{r}).month<=month_lims_for_mean(2) & store_gdap.surf.(regions{r}).MLD>=MLD_cutoff & ...
                        store_gdap.surf.(regions{r}).month==mon & store_gdap.surf.(regions{r}).year==year_range(yy) ...
                        & ~isnan(store_gdap.surf.(regions{r}).delta_pCO2_unified);
                    
                    temp_filt.(regions{r}).gdap.(properties{p})((yy-1)*12+mon,1) = nanmean(store_gdap.surf.(regions{r}).(properties{p})(gdap_index));
                    temp_filt.(regions{r}).gdap.(properties{p})((yy-1)*12+mon,2) = nanstd(store_gdap.surf.(regions{r}).(properties{p})(gdap_index));
                    
                end
                end
                
            end
        end
        
    end
end

%
out_delt_o2_pCO2_mean = NaN(2,6);
out_delt_o2_pCO2_std = NaN(2,6);

data_sources = {'argo' 'gdap' 'socat'};

for r = 1:2
    
    for q = 1:length(props_out)
        for dd = 1:3
            
            
            out_delt_o2_pCO2_mean(r, dd + 3*(q-1)) = nanmean(temp_filt.(regions{r}).(data_sources{dd}).(properties{props_out(q)})(:,1));
            out_delt_o2_pCO2_std(r, dd + 3*(q-1)) = nanstd(temp_filt.(regions{r}).(data_sources{dd}).(properties{props_out(q)})(:,1));
            
        end
    end
end

clear r q dd p gdap_index socat_index argo_index yy mon

%% Table 2
% combined average O2 (Argo and GDAP) and delta pCO2 (Argo and SOCAT)
% "out_delta_o2_pCO2_[mean/std]" are used for table 2
for r = 1:2
    
    q = 1; % ml delta O2
    disp(properties{props_out(q)})
    out_delt_o2_pCO2_mean(r, 3) = nanmean([temp_filt.(regions{r}).(data_sources{1}).(properties{props_out(q)})(:,1) ; ...
        temp_filt.(regions{r}).(data_sources{2}).(properties{props_out(q)})(:,1)]);
    out_delt_o2_pCO2_std(r, 3) = nanstd([temp_filt.(regions{r}).(data_sources{1}).(properties{props_out(q)})(:,1) ; ...
        temp_filt.(regions{r}).(data_sources{2}).(properties{props_out(q)})(:,1)]);
    
    q = 2; % delta pCO2
    disp(properties{props_out(q)})
    
    out_delt_o2_pCO2_mean(r, 25) = nanmean([temp_filt.(regions{r}).(data_sources{1}).(properties{props_out(q)})(:,1); ...
        temp_filt.(regions{r}).(data_sources{3}).(properties{props_out(q)})(:,1)]);
    
    out_delt_o2_pCO2_std(r, 25) = nanstd([temp_filt.(regions{r}).(data_sources{1}).(properties{props_out(q)})(:,1); ...
        temp_filt.(regions{r}).(data_sources{3}).(properties{props_out(q)})(:,1)]);
    
    
end

clear r q

%% saving out o2 and ptmp for supplemental figure
temp_date_vec = datevec(temp_filt.GMT_Matlab);
o2_iav_argo.ptmp = NaN(length(year_range),1);
o2_iav_argo.ml_o2 = NaN(length(year_range),1);

for y = 1:length(year_range)
    
    y_index = temp_date_vec(:,1) == year_range(y);
    o2_iav_argo.ptmp(y) = nanmean(temp_filt.pacific.argo.ml_PTMP(y_index,1));
    o2_iav_argo.ml_o2(y) = nanmean(temp_filt.pacific.argo.ml_o2(y_index,1));
    o2_iav_argo.SLP(y) = nanmean(temp_filt.pacific.argo.SLP(y_index,1));

end


%% Figure 4 - delta pCO2 and O2
c_map = brewermap(10, 'Paired');
marker_size = 3;

clf

set(gcf, 'units', 'inches')
paper_w = 8; paper_h =6;
set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]); clear paper_w paper_h

plot_count = 0;
d = NaN(4,1);


for p=[2 10]
    
    for r=1:2
        plot_count = plot_count+1;
        d(plot_count) = subplot(2,2,plot_count);
        hold on; grid on
%         title(upper(regions{r}))
        if r==1
            title('Pacific')
        elseif r==2
            title('Indian')
        end
        
        plot([datenum(1991,1,1) datenum(2021,1,1)], [0 0], 'k-')
        
        if ~isempty(store_gdap.(regions{r}).(properties{p}))
            
            gdap_index = store_gdap.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(1) &  ...
                store_gdap.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) & ...
                store_gdap.(regions{r}).month>=month_lims_for_mean(1) & store_gdap.(regions{r}).month<=month_lims_for_mean(2) & store_gdap.(regions{r}).MLD>=MLD_cutoff;
            pG = plot(store_gdap.(regions{r}).GMT_Matlab(gdap_index), store_gdap.(regions{r}).(properties{p})(gdap_index), 's', 'color', c_map(9,:),'markerfacecolor', c_map(9,:), 'markersize', marker_size);
            
            disp([regions{r} ' ' properties{p} ' gdap ' num2str(nanmean(store_gdap.(regions{r}).(properties{p})(gdap_index))) ' std ' num2str(nanstd(store_gdap.(regions{r}).(properties{p})(gdap_index)))])
            
            
            errorbar(temp_filt.GMT_Matlab, temp_filt.(regions{r}).gdap.(properties{p})(:,1), temp_filt.(regions{r}).gdap.(properties{p})(:,2), 'color', c_map(10,:), 'marker', 's', 'linewidth', 1.5, 'markersize', marker_size-1)
            
            disp([regions{r} ' ' properties{p} ' gdap month ' num2str(nanmean(temp_filt.(regions{r}).gdap.(properties{p})(:,1))) ...
                ' std ' num2str(nanstd(temp_filt.(regions{r}).gdap.(properties{p})(:,1)))])
            
            
        end
        if ~isempty(store_socat.(regions{r}).(properties{p}))
            
            socat_index = store_socat.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(1) &  ...
                store_socat.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) & ...
                store_socat.(regions{r}).month>=month_lims_for_mean(1) & store_socat.(regions{r}).month<=month_lims_for_mean(2) & store_socat.(regions{r}).MLD>=MLD_cutoff ;
            pS = plot(store_socat.(regions{r}).GMT_Matlab(socat_index), store_socat.(regions{r}).(properties{p})(socat_index), 'd', 'color', c_map(3,:), 'markerfacecolor', c_map(3,:), 'markersize', marker_size);
            
            disp([regions{r} ' ' properties{p} ' socat ' num2str(nanmean(store_socat.(regions{r}).(properties{p})(socat_index))) ...
                ' std ' num2str(nanstd(store_socat.(regions{r}).(properties{p})(socat_index)))])
            
            errorbar(temp_filt.GMT_Matlab, temp_filt.(regions{r}).socat.(properties{p})(:,1), temp_filt.(regions{r}).socat.(properties{p})(:,2), 'color', c_map(4,:), 'marker', 'd', 'linewidth', 1.5, 'markersize', marker_size-1)
            
            disp([regions{r} ' ' properties{p} ' socat month ' num2str(nanmean(temp_filt.(regions{r}).socat.(properties{p})(:,1))) ...
                ' std ' num2str(nanstd(temp_filt.(regions{r}).socat.(properties{p})(:,1)))])
            
        end
        
        argo_index = store_argo.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(1) &  ...
            store_argo.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) & ...
            store_argo.(regions{r}).month>=month_lims_for_mean(1) & store_argo.(regions{r}).month<=month_lims_for_mean(2) & store_argo.(regions{r}).MLD>=MLD_cutoff;
        pA = plot(store_argo.(regions{r}).GMT_Matlab(argo_index), store_argo.(regions{r}).(properties{p})(argo_index), 'o', 'markerfacecolor', c_map(1,:), 'color', c_map(1,:), 'markersize', marker_size);
        
        disp([regions{r} ' ' properties{p} ' argo ' num2str(nanmean(store_argo.(regions{r}).(properties{p})(argo_index))) ...
            ' std ' num2str(nanstd(store_argo.(regions{r}).(properties{p})(argo_index)))])
        
        eA = errorbar(temp_filt.GMT_Matlab, temp_filt.(regions{r}).argo.(properties{p})(:,1), temp_filt.(regions{r}).argo.(properties{p})(:,2), 'color', c_map(2,:), 'marker', 'o', 'linewidth', 1.5, 'markersize', marker_size-1);
        
        disp([regions{r} ' ' properties{p} ' argo month ' num2str(nanmean(temp_filt.(regions{r}).argo.(properties{p})(:,1))) ...
            ' std ' num2str(nanstd(temp_filt.(regions{r}).argo.(properties{p})(:,1)))])
        
        
        if p==2
            ylabel('\Delta[O_2] (\mumol kg^-^1)')
        elseif p==10
            ylabel('\Delta{\itp}CO_2 (\muatm)')
        end
        
        %         ylabel(properties{p})
        if plot_count==3
            if ~isempty(store_socat.(regions{r}).(properties{p}))
                
                legend([pA pG pS], 'Argo', 'GLODAP', 'SOCAT', 'location', 'northwest')
            else
                legend([pA pG], 'Argo', 'Glodap', 'location', 'northwest')
            end
        end
        set(gca, 'xtick', datenum(1990:5:2020, 1, 1))
        datetick('x', 'YYYY', 'keepticks')
        
        clear eA pA argo_index socat_index pS pG gdap_index
    end
    
    disp(' ')
end

set([d(3) d(4)], 'ylim', [-55 55], 'ytick', -60:20:60)
set([d(1) d(2)], 'ylim', [-30 5])
x1 = xlabel('Year');

x1_pos = get(x1, 'position');
%
set(x1, 'position', x1_pos+[-365*20 -5 0])

plot_filename = ['Fig_4_Time Series Ind Pac- Pac Lat ' num2str(geo_lims.lat_lims_pac) ...
    ' Lon ' num2str(geo_lims.pac_lims) ' Ind Lat ' num2str(geo_lims.lat_lims_ind) ' Lon ' num2str(geo_lims.ind_lims) ...
    ' Months ' num2str(month_lims_for_mean(1)) ' to ' num2str(month_lims_for_mean(2))  ...
    ' MLD cut ' num2str(MLD_cutoff) ' ' properties{2} ' '  properties{10} '_GDAP_SOCAT_v9'];
print(gcf, '-dpdf', '-r800',  [plot_dir plot_filename '.pdf'])

clear x1 x1_pos d p eA pA argo_index plot_count r temp_filt temp_date_vec y y_index

%% Figure 5 - Temp and density histograms

prop_names = {'[O_2] (\mumol kg^-^1)'; ...
    ' '; ...
    '[DIC] (\mumol kg^-^1)' ; ...
    '[NO_3^-] (\mumol kg^-^1)'; ...
    'pCO_2 (\muatm)'; ...
    '\theta (\circC)'; ...
    'Sal'; ...
    ' ' ; ...
    '\sigma_\theta (kg m^-^3)' ; ...
    '\DeltapCO_2 (\muatm)'; ...
    ' ' ; ...
    ' ' ; ...
    ' ' ;...
    'Lon. (\circE)'} ;



c_map = brewermap(10, 'Paired');
marker_size = 3;

plot_glodap = 0;

clf

set(gcf, 'units', 'inches')
paper_w = 8; paper_h =8;
set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]); clear paper_w paper_h

plot_count = 0;
d = NaN(4,1);


for p=[10 14 9 6]
    if p==10
        BW = 5; % Bin Width
    elseif p==14
        BW = 20;
    elseif p==9
        BW = 0.025;
    elseif p==6
        BW = 0.5;
    end
    
    for r=1:2
        plot_count = plot_count+1;
        d(plot_count) = subplot(4,2,plot_count);
        hold on; grid on
        
        if plot_count<3
            title(upper(regions{r}))
            if r==1
                title('Pacific')
            elseif r==2
                title('Indian')
            end
        end
        %         sig_thet_range = mode_vol.y2005.sig_th(mw_prop.(regions{r}).all.per_5_index);
        
        %         plot([datenum(1991,1,1) datenum(2021,1,1)], [0 0], 'k-')
        
        if p==9
            prop_adjust=-1000;
        else
            prop_adjust = 0;
        end
        if plot_glodap==1
            
            if ~isempty(store_gdap.(regions{r}).(properties{p}))
                gdap_index = store_gdap.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(1) &  ...
                    store_gdap.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) & ...
                    store_gdap.(regions{r}).month>=month_lims_for_mean(1) & store_gdap.(regions{r}).month<=month_lims_for_mean(2) & store_gdap.(regions{r}).MLD>=MLD_cutoff ...
                    & ~isnan(store_gdap.(regions{r}).delta_pCO2);
                %             pG = plot(store_gdap.(regions{r}).GMT_Matlab(gdap_index), store_gdap.(regions{r}).(properties{p})(gdap_index), 's', 'color', c_map(9,:),'markerfacecolor', c_map(9,:), 'markersize', marker_size);
                
                pG = histogram(store_gdap.(regions{r}).(properties{p})(gdap_index)+ prop_adjust, 10, 'Normalization', 'probability');
                pG.FaceColor = c_map(10,:);
                
                %             disp([regions{r} ' ' properties{p} ' gdap ' num2str(nanmean(store_gdap.(regions{r}).(properties{p})(gdap_index))) ' std ' num2str(nanstd(store_gdap.(regions{r}).(properties{p})(gdap_index)))])
                
                
                %             errorbar(temp_filt.GMT_Matlab, temp_filt.(regions{r}).gdap.(properties{p})(:,1), temp_filt.(regions{r}).gdap.(properties{p})(:,2), 'color', c_map(10,:), 'marker', 's', 'linewidth', 1.5, 'markersize', marker_size-1)
                
                %             disp([regions{r} ' ' properties{p} ' gdap month ' num2str(nanmean(temp_filt.(regions{r}).gdap.(properties{p})(:,1))) ...
                %                 ' std ' num2str(nanstd(temp_filt.(regions{r}).gdap.(properties{p})(:,1)))])
                
                
                
            end
        end
        argo_index = store_argo.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(1) &  ...
            store_argo.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) & ...
            store_argo.(regions{r}).month>=month_lims_for_mean(1) & store_argo.(regions{r}).month<=month_lims_for_mean(2) & store_argo.(regions{r}).MLD>=MLD_cutoff ...
            & ~isnan(store_argo.(regions{r}).delta_pCO2_unified);
        %         pA = plot(store_argo.(regions{r}).GMT_Matlab(argo_index), store_argo.(regions{r}).(properties{p})(argo_index), 'o', 'markerfacecolor', c_map(1,:), 'color', c_map(1,:), 'markersize', marker_size);
        pA = histogram(store_argo.(regions{r}).(properties{p})(argo_index)+ prop_adjust, 10,'Normalization', 'probability');
        pA.FaceColor = c_map(2,:);
        
        if ~isempty(store_socat.(regions{r}).(properties{p}))
            
            socat_index = store_socat.(regions{r}).ml_pdens-1000>=mw_prop.(regions{r}).all.sig_thet_range(1) &  ...
                store_socat.(regions{r}).ml_pdens-1000<mw_prop.(regions{r}).all.sig_thet_range(end) & ...
                store_socat.(regions{r}).month>=month_lims_for_mean(1) & store_socat.(regions{r}).month<=month_lims_for_mean(2) & store_socat.(regions{r}).MLD>=MLD_cutoff ...
                & ~isnan(store_socat.(regions{r}).delta_pCO2_unified);
            %             pS = plot(store_socat.(regions{r}).GMT_Matlab(socat_index), store_socat.(regions{r}).(properties{p})(socat_index), 'd', 'color', c_map(3,:), 'markerfacecolor', c_map(3,:), 'markersize', marker_size);
            pS = histogram(store_socat.(regions{r}).(properties{p})(socat_index)+ prop_adjust, 10,'Normalization', 'probability');
            pS.FaceColor = c_map(4,:);
            %             disp([regions{r} ' ' properties{p} ' socat ' num2str(nanmean(store_socat.(regions{r}).(properties{p})(socat_index))) ...
            %                 ' std ' num2str(nanstd(store_socat.(regions{r}).(properties{p})(socat_index)))])
            %
            %             errorbar(temp_filt.GMT_Matlab, temp_filt.(regions{r}).socat.(properties{p})(:,1), temp_filt.(regions{r}).socat.(properties{p})(:,2), 'color', c_map(4,:), 'marker', 'd', 'linewidth', 1.5, 'markersize', marker_size-1)
            %
            %             disp([regions{r} ' ' properties{p} ' socat month ' num2str(nanmean(temp_filt.(regions{r}).socat.(properties{p})(:,1))) ...
            %                 ' std ' num2str(nanstd(temp_filt.(regions{r}).socat.(properties{p})(:,1)))])
            
        end
        
        
        %         disp([regions{r} ' ' properties{p} ' argo ' num2str(nanmean(store_argo.(regions{r}).(properties{p})(argo_index))) ...
        %             ' std ' num2str(nanstd(store_argo.(regions{r}).(properties{p})(argo_index)))])
        %
        %         eA = errorbar(temp_filt.GMT_Matlab, temp_filt.(regions{r}).argo.(properties{p})(:,1), temp_filt.(regions{r}).argo.(properties{p})(:,2), 'color', c_map(2,:), 'marker', 'o', 'linewidth', 1.5, 'markersize', marker_size-1);
        %
        %         disp([regions{r} ' ' properties{p} ' argo month ' num2str(nanmean(temp_filt.(regions{r}).argo.(properties{p})(:,1))) ...
        %             ' std ' num2str(nanstd(temp_filt.(regions{r}).argo.(properties{p})(:,1)))])
        
        
        xlabel(prop_names{p})
        %         ylabel(properties{p})
        
        if plot_glodap==1
            pA.BinLimits = [min([pA.BinLimits pG.BinLimits pS.BinLimits]) max([pA.BinLimits pG.BinLimits pS.BinLimits])];
            
            pG.BinLimits = [min([pA.BinLimits pG.BinLimits pS.BinLimits]) max([pA.BinLimits pG.BinLimits pS.BinLimits])];
            pS.BinLimits = [min([pA.BinLimits pG.BinLimits pS.BinLimits]) max([pA.BinLimits pG.BinLimits pS.BinLimits])];
            pA.NumBins = 10; pS.NumBins = 10; pG.NumBins = 10;
            
        else
            %                         pA.BinLimits = [min([pA.BinLimits pS.BinLimits]) max([pA.BinLimits pS.BinLimits])];
            %                         pS.BinLimits = [min([pA.BinLimits pS.BinLimits]) max([pA.BinLimits pS.BinLimits])];
            %                         pA.NumBins = 10; pS.NumBins = 10;
            
            set([pA pS], 'BinWidth', BW)
        end
        if plot_count==2
            if plot_glodap==1
                legend([pA, pS, pG], 'Argo', 'SOCAT', 'GLODAP')
            else
                legend([pA, pS], 'Argo', 'SOCAT')
            end
        end
        % pS.BinEdges = pA.BinEdges;
        % pG.BinEdges = pA.BinEdges;
    end
    
    %     disp(' ')
end
% set([d(3) d(4)], 'ylim', [33.5 35.5])
% set([d(1) d(2)], 'ylim', [2 14])
y1 = ylabel(d(5), 'Relative Frequency');
y1_pos = get(y1, 'Position');
%
set(y1, 'position', y1_pos+[-.03 .33 0])

d7_pos = get(d(7), 'position');

set(d(7), 'position', d7_pos+[0 -0.005 0 0])
set(d(8), 'position', get(d(8), 'position')+[0 -0.005 0 0])


plot_filename = ['Fig_5_Histogram Props- Pac Lat ' num2str(geo_lims.lat_lims_pac) ...
    ' Lon ' num2str(geo_lims.pac_lims) ' Ind Lat ' num2str(geo_lims.lat_lims_ind) ' Lon ' num2str(geo_lims.ind_lims) ...
    ' Months ' num2str(month_lims_for_mean(1)) ' to ' num2str(month_lims_for_mean(2))  ...
    ' MLD cut ' num2str(MLD_cutoff) ' ' properties{2} ' '  properties{10} '_SOCAT_v11'];
print(gcf, '-dpdf', '-r800',  [plot_dir plot_filename '.pdf'])

clear p pA plot_count plot_glodap prop_adjust y1 y1_pos socat_index r pS d7_pos d c_map bsose_color_sal bsose_color_temp argo_index 
%% Figure 6 - bSOSE results
sose_data_dir = [home_dir 'Work/Projects/2020_02 SO SAMW Variability_Cerovecki/Data_from_Ivana/2021_12_29/'];

bsose = load([sose_data_dir 'MLDanom_Timav_mo_8_9_yrs_2013_2019_iter135_SouthOc_mask_200_35_64.mat']);

%% Fig 6 Maps

front_bounds = load([data_dir 'ARGO_O2_Floats/Front_definitions/Gray_5_regions/regional_boundaries_5zone.mat']);

bsose_years = 2013:2019;
bsose_props = {'mld_h_an' 'mld_dic_an' 'mld_no3_an' 'mld_o2_an'};
bsose_prop_names = {'MLD' 'DIC' 'Nitrate' 'Oxygen'};
color_map_names = {'PRGn' 'PuOr' 'PiYg' 'RdBu'};
bsose_prop_units = {'(m)' '(\mumol kg^-^1)' '(\mumol kg^-^1)' '(\mumol kg^-^1)'};

c_lims = [-100 100 ;  -10 10 ; -2 2; -15 15];

for pp = 1:length(bsose_props)
    plot_filename = ['Fig_6_bsose_annual_anomalies_' bsose_props{pp} '_v2' ];
    
    clf
    set(gcf, 'units', 'inches')
    paper_w = 9; paper_h =9.7;
    
    coast = load('coast');
    % colormap = inferno;
    colormap = brewermap(30, color_map_names{pp});
    if pp==2 || pp==4
        colormap = flipud(colormap);
    end
    set(gcf, 'colormap', colormap)
    
    set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]);
    lat_map = double(bsose.yc(1,:))';
    lon_map = double(bsose.xc(:,1));
    %
    R = georasterref('RasterSize', [length(lat_map) length(lon_map)], 'RasterInterpretation', 'Cells', ...
        'LatitudeLimits', [min(lat_map)-.5 max(lat_map)+.5], 'LongitudeLimits', [min(lon_map)-.5 max(lon_map)+.5]);
    boundary_width=1.5;
    
    ph = NaN(length(bsose_years),1);
    plot_index = 0;
    
    for yy=1:length(bsose_years)
        plot_index = plot_index +1;
        
        
        
        ph(yy) = subplot(7,1,plot_index);
        
        
        data1 = squeeze(bsose.(bsose_props{pp})(:,:,yy))';
        % data1(isnan(data1))=0;
        temp_output = data1;
        
        % Creates a polar stereographic plot with the origin and bounds below
        % axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -34], 'FFaceColor', [.9 .9 .9])
        axesm('lambertstd', 'MapLatLimit',[-60 -42],'MapParallels',[-75 -15],'MapLonLimit',[170 -70])
        
        axis off; framem on; gridm on; mlabel on; plabel on;
        setm(gca,'MLabelParallel',-10, 'fontsize', 14)
        hold on
        hm = meshm(temp_output, R);
        
        set(hm, 'facealpha', 'texturemap', 'alphadata', double(~isnan(temp_output)))
        set(gcf, 'renderer', 'opengl')
        
        % surfm(lon_map, lat_map, temp_output)
        
        geoshow(coast.lat,coast.long,'DisplayType','polygon', 'facecolor', [.5 .5 .5])
        clear temp_output
        % geoshow(gca, lat_siz, lon_siz, 'displaytype', 'line', 'linewidth', boundary_width, 'color', 'k') % Seasonal Ice Zone
        geoshow(gca, front_bounds.lat_pf, front_bounds.lon_pf, 'displaytype', 'line', 'linewidth', boundary_width, 'color', 'k') % polar front
        geoshow(gca, front_bounds.lat_saf, front_bounds.lon_saf, 'displaytype', 'line', 'linewidth', boundary_width, 'color', 'k') % Subantarctic Front
        % geoshow(gca, lat_stf, lon_stf, 'displaytype', 'line', 'linewidth', boundary_width, 'color', 'k') % Subtropical Front
        setm(gca, 'mlabellocation', [0 90 180 ], 'mlabelparallel', -40)
        setm(gca, 'plabellocation', [-35 -90], 'plabelmeridian', -150)
        
        caxis(c_lims(pp,:))
        %     set(ph(1), 'clim', map_prop_lims(p,:))
        %
        set(gca, 'color', [.5 .5 .5]);
        %
        %
        
        %         c1 = colorbar;
        % ylabel(c1, 'Mean MLD 2005-2020 (m)', 'interpreter', 'none')
        % set(c1, 'fontsize', 18)
        % col_pos = c1.Position;
        % c1.Position = col_pos.*[1 1 1 .9]+[0.03 .04 0.015 0];
        if yy==1
            title(bsose_prop_names{pp},'fontsize', 16)
        end
        text(-.65,-.83, num2str(bsose_years(yy)), 'fontsize', 15)
        % text(-1.1, 1.1, [' Months: ' num2str(month_lims_for_mean(1)) ' to ' num2str(month_lims_for_mean(2))], 'fontsize', 15)
    end
    %
    all_pos = NaN(length(bsose_years),4);
    for yy = 1:length(bsose_years)
        all_pos(yy,:) = get(ph(yy), 'position');
    end
    %
    
    y_offset = -.035;
    m_scale = 0.08;
    for yy = 1:length(bsose_years)
        set(ph(yy), 'position', all_pos(yy,:)+[-0.1 y_offset m_scale m_scale])
    end
    
    
%     print(gcf, '-dpdf', '-r300',  [plot_dir plot_filename '.pdf'])
    %
    c1 = colorbar('southoutside');
    ylabel(c1, bsose_prop_units{pp})
    set(ph(yy), 'position', all_pos(yy,:)+[-0.1 y_offset m_scale m_scale])

    c1_pos = get(c1, 'position');
    
    set(c1, 'position', c1_pos+[0.1 0.008 -0.2 0], 'fontsize', 12)
    print(gcf, '-dpdf', '-r800',  [plot_dir plot_filename 'w_colorbar.pdf'])
    
end


clear temp_prop temp_lat temp_lon y_offset m_scale all_pos c1 c1_pos R temp_output hm 

%% Fig 6 top panel - climate indices time series
c_index_data_dir = [home_dir 'Work/Projects/2020_02 SO SAMW Variability_Cerovecki/Data_from_Ivana/2021_12_10/'];
c_indices_IC = load([c_index_data_dir 'sam_nino_2005_2021.mat']);

c_map = brewermap(4, 'Paired');
% load  sam_nino_2005_2019.mat

NUM_SMOOTH=1;
sam_sm=smooth(c_indices_IC.sam,NUM_SMOOTH); % <<---- want to smooth? This is 2-month smoothing
enso_sm=smooth(c_indices_IC.nino34,NUM_SMOOTH);
ntt=length(c_indices_IC.sam);
time=[1:ntt];

figure(1); clf

set(gcf, 'units', 'inches')
paper_w = 9; paper_h =6;

set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]);

subplot(3,1,1)
% set(gca, 'ColorOrder', [0 0 1; 1 0 0; 0 1 0; 0.5 0.5 0.5], 'NextPlot', 'replacechildren');

y_ind = 96;

% [hax,hLine1,hLine2]=plotyy(time(y_ind+1:ntt),sam_sm(y_ind+1:ntt),time(y_ind+1:ntt),c_indices_IC.nino34(y_ind+1:ntt)); hold on; grid on;
hLine1=plot(time(y_ind+1:ntt),sam_sm(y_ind+1:ntt), '-'); hold on; grid on;
hLine3 = plot(time(1:y_ind+1),sam_sm(1:y_ind+1), '-');



yyaxis right

hLine2=plot(time(y_ind+1:ntt),enso_sm(y_ind+1:ntt), '-'); hold on; grid on;
hLine4 = plot(time(1:y_ind+1),enso_sm(1:y_ind+1), '-');
hax = gca;

set(hLine2,'linewidth',2); set(hLine1,'linewidth',2);

set(hLine1, 'Color', c_map(2,:))
set(hLine2, 'Color', c_map(4,:))

% [~,hLine3,hLine4]=plotyy(time(1:y_ind+1),sam_sm(1:y_ind+1),time(1:y_ind+1),c_indices_IC.nino34(1:y_ind+1)); hold on; grid on;
set(hLine3,'linewidth',2); set(hLine4,'linewidth',2);
%
set(hLine3, 'Color', c_map(1,:))
set(hLine4, 'Color', c_map(3,:))
%
% set(hax(1), 'yColor', c_map(2,:))
% set(hax(2), 'yColor', c_map(4,:))

yyaxis left
set(gca,'xTick',[1:12:ntt]);
set(gca, 'yColor', c_map(2,:))

set(hax(1),'XTick',[1:12:ntt]);box on ;
set(hax(1),'XTickLabel',{' ',' ',' ',' ',' ',' '});
set(gca, 'xlim', [1 ntt-21])

set(hax(1),'XTickLabel',{'2005',' ','2007',' ','2009',' ','2011',' ','2013',' ','2015',' ','2017',' ','2019', ' '});
axis(hax(1),[1 ntt -5 5]);
set(hax(1),'YTick',[-4:2:4]); ylabel(hax(1),'SAM' ,'fontsize',10);
%
yyaxis right
set(gca, 'yColor', c_map(4,:))

% set(hax(2),'XTick',[1:12:ntt]);
set(gca,'ytick',[-3:1:3]);
axis(gca,[1 ntt -3 3]);
%
set(gca, 'xlim', [1 ntt-21])

set(gca,'fontsize',11); grid on
ylabel(gca,'NINO3.4 (^{\circ}C)','fontsize',10);
plot(time,zeros(size(time)),'linewidth',1,'color',[0.7 0.7 0.7]); grid on
% text(3.8,1.95,'a)','fontsize',12);
%
plot_filename = 'Fig_6 climate indices_v2';
%
print(gcf, '-dpdf', '-r800',  [plot_dir plot_filename '.pdf'])

clear all_pos c1 c1_pos c_lims c_map color_map_names colormap data1 enso_sm hax hLine1 hLine2 hLine3 hLine4 hm ntt NUM_SMOOTH paper_h paper_w ph 
clear plot_filename plot_index pp R sam_sm y_ind y_offset yy lat_map lon_map
clear time temp_time temp reg_index reg_names cmap

%% Fig 7: bSOSE IAV
bSOSE_anav = load([cerovecki_dir '/Data_from_Ivana/2021_12_10/MLD_AnAv_SAMW_range8_9_SEPac_200_m_025_ST_iter135_45_64_246_290_se_pacific.mat']);
bsose_fields = fieldnames(bSOSE_anav);
reg_names = {'samw_c_pacific' 'samw_se_pacific' 'samw_pacific'};

plot_filename = ['Fig_7_BSOSE_IAV_v3'];
x = 2013:2019;
plot_symbols = {'x', 'o', 'd'};
plot_colors = brewermap(3, 'Set1');
plot_colors (3,:) = [0 0 0];

clf
set(gcf, 'units', 'inches')
paper_w = 7; paper_h =9;
set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]); clear paper_w paper_h
plot_count = 0;
for s=[7 5 4 1 2 3]
    
    if s==2
        %         sose_prop = 'no3_mld_clim';
        prop_name = '[NO_3^-]';
        prop_unit = '\mumol kg^-^1';
        %         anomaly_text = '[NO_3^-] anomaly';
    elseif s==1
        %         sose_prop = 'dic_mld_clim';
        prop_name = '[DIC]' ;
        prop_unit = '\mumol kg^-^1';
        
        %         anomaly_text = '[DIC] anomaly';
        
    elseif s==3
        %         sose_prop = 'o2_mld_clim';
        prop_name = '[O_2]' ;
        prop_unit = '\mumol kg^-^1';
        
        %         anomaly_text = '[O_2] anomaly';
        
    elseif s==4
        prop_name = 'Sal.' ;
        prop_unit = 'PSS-78';

    elseif s==5
        %         sose_prop = 'T_mld_clim';
        prop_name =  '\theta';
        prop_unit = '\circC';
        
        %         anomaly_text = '\theta anomaly';
        
    elseif s==7
        prop_name = 'Vol.' ;
        prop_unit = 'm^3';
        
    end
    
    plot_count = plot_count+1;
    subplot(3,2,plot_count)
    hold on
    
    l_n = NaN(size(reg_names,2),1);
    
    for r = 1:3
        
        
        reg_index = find(cellfun(@isempty,strfind(bsose_fields, reg_names{r}))==0);
        
        y = bSOSE_anav.(bsose_fields{reg_index(s)});
        
        l_n(r) =  plot(x, y, plot_symbols{r}, 'linestyle', '-', 'linewidth',2, 'color', plot_colors(r,:));
        
    end
    ylabel(prop_unit)
    title(prop_name)
    if plot_count==1
        l1 = legend(l_n, 'C Pac', 'SE Pac', 'Pac', 'location', 'northeast', 'fontsize', 11);
        %                 l1 = legend(l_n, reg_names, 'location', 'northeast', 'fontsize', 11);
        
    end
    
end
x1=xlabel('Year');

set(x1, 'position', get(x1, 'position')+[-5.2 -1.5 0]);
l1_pos = get(l1, 'position');

set(l1, 'position', l1_pos+[.01 .01 0 0])

print(gcf, '-dpdf', '-r800',  [plot_dir plot_filename '.pdf'])

clear l1 l1_pos x1 prop_name l_n prop_unit plot_count r s plot_colors x  
%% Supplemental Figures

%% Fig S1 - SAMW formation properties and volume vs. density

plot_filename = ['Fig_S1_Mode Water formation properties overalaid- weighted ' ...
    ' Months ' num2str(month_lims_for_mean(1)) ' to ' num2str(month_lims_for_mean(2)) ' Dens. int for mean ' num2str(dens_interval) ' MLD cut ' num2str(MLD_cutoff) '_v14_.03MLD'];

plot_colors = brewermap(10, 'Dark2');

prop_lims = [240 340;...
    -20 20;...
    2090 2200;...
    0 35;...
    365 470;...
    0 14;...
    33.5 36];

prop_names = {'[O_2] (\mumol kg^-^1)'; ...
    ' '; ...
    '[DIC] (\mumol kg^-^1)' ; ...
    '[NO_3^-] (\mumol kg^-^1)'; ...
    '{\itp}CO_2 (\muatm)'; ...
    '\theta (\circC)'; ...
    'PSS-78'} ;

plot_map = [1 4 3 5 6 7];

clf

set(gcf, 'units', 'inches')
paper_w = 5; paper_h =8;
set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]); clear paper_w paper_h

d = NaN(8,1);
axis_font_size = 12;

for r=1:2
    
    
    for pl = 1:length(plot_map)
        hold on
        xlabel('\sigma_\theta', 'fontsize', axis_font_size)
        grid on
        
     
        plot_number = pl;
        %         axes(d(plot_number))
        d(plot_number) = subplot(6,1,plot_number);
        
        p = plot_map(pl);
        
        if r==1
            pcol = 1;
        else
            pcol = 3;
        end
        errorbar( mode_vol.(mode_years{1}).sig_th, squeeze(mw_prop.(regions{r}).wint.(properties{p})(:,1)), squeeze(mw_prop.(regions{r}).wint.(properties{p})(:,2)), 'm', 'linewidth', 2, 'color', plot_colors(pcol,:))
        
        fig_2_data.prop_sig_th = mode_vol.(mode_years{1}).sig_th';
        
        fig_2_data.(regions{r}).(properties{p}) = mw_prop.(regions{r}).wint.(properties{p});
        
        
        ylabel(prop_names{p}, 'fontsize', 10)
        set(gca, 'ylim', prop_lims(p,:))
        set(gca, 'ycolor', 'k')
        set(gca, 'xlim', [26.4 27.3]);
       
    end
   
    
    if pl==6 && r==2
        legend('Pacific', 'Indian')
    end
end

set(d(1), 'ytick', 240:40:340)
set(d(4), 'ytick', 370:30:460)

print(gcf, '-dpdf', '-r800', [plot_dir plot_filename '.pdf' ])

clear pcol d plot_number pl lon_lims lat_lims paper_w paper_h plot_lims plot_map
clear r 

%% Fig S2 - BSOSE float bgc difference evaluation

clear bsose_comp

BSOSE_dir = [data_dir 'Model_Output/SOSE/2013-2019_ITER135_1_6deg/'];

viz_files = dir([BSOSE_dir 'VIZ*.nc']);
for f = 1:length(viz_files)
    file_info = ncinfo([BSOSE_dir viz_files(f).name]);
    
    for v = 1:length(file_info.Variables)
        bsose_comp(f).year.(file_info.Variables(v).Name) = ncread([BSOSE_dir viz_files(f).name], file_info.Variables(v).Name);
    end
end

clear f v
%%
clear var_diff
comp_vars = {'prof_T'; 'prof_S'; 'prof_O2'; 'prof_NO3'; 'prof_DIC'};

for cv = 1:length(comp_vars)
    var_diff.(comp_vars{cv}) = [];
end
var_diff.prof_date = [];
var_diff.MLD = [];
var_diff.prof_lat = [];
var_diff.prof_lon = [];

% calculate differences for each year, concatenate into a single array for
% each variable
for f = 1:length(viz_files)
    for cv = 1:length(comp_vars)
        temp_var = bsose_comp(f).year.(comp_vars{cv}) - bsose_comp(f).year.([comp_vars{cv} 'estim']);
        var_diff.(comp_vars{cv}) = cat(2, var_diff.(comp_vars{cv}), temp_var);
        
    end
    % add datetime vectors
    var_diff.prof_date = [var_diff.prof_date ; bsose_comp(f).year.prof_date];
    var_diff.prof_lat = [var_diff.prof_lat ; bsose_comp(f).year.prof_lat];
    var_diff.prof_lon = [var_diff.prof_lon ; bsose_comp(f).year.prof_lon];
    
    temp_MLD = NaN(length(bsose_comp(f).year.prof_date),1);
    % calculate MLDs
    for p = 1:size(bsose_comp(f).year.(comp_vars{cv}),2)
        temp_MLD(p) = ...
            mld_dbm(bsose_comp(f).year.prof_Testim(:,p), ...
            bsose_comp(f).year.prof_Sestim(:,p), ...
            bsose_comp(f).year.prof_depth, 1);
    end
    var_diff.MLD = [var_diff.MLD ; temp_MLD];


end

clear cv f temp_var

var_diff.date_vec = datevec(var_diff.prof_date);
var_diff.prof_lon(var_diff.prof_lon<0) = var_diff.prof_lon(var_diff.prof_lon<0)+360;

%% calculate binned means

% split by region:
for r = 1:3
    
    
    if r==1
        var_date_ind = (var_diff.prof_lat>geo_lims.lat_lims_pac(1) & var_diff.prof_lat<geo_lims.lat_lims_pac(2) & ...
            var_diff.prof_lon>=geo_lims.pac_lims(1) & var_diff.prof_lon<geo_lims.pac_lims(2)) & ...
            sum(var_diff.date_vec(:,2)==month_lims_for_mean,2)==1 & var_diff.MLD>=200;
        region_name = regions{r};
        
    elseif r==2
        var_date_ind = (var_diff.prof_lat>geo_lims.lat_lims_ind(1) & var_diff.prof_lat<geo_lims.lat_lims_ind(2) & ...
            var_diff.prof_lon>=geo_lims.ind_lims(1) & var_diff.prof_lon<geo_lims.ind_lims(2)) & ...
            sum(var_diff.date_vec(:,2)==month_lims_for_mean,2)==1 & var_diff.MLD>=200;
        region_name = regions{r};
    else
        var_date_ind = sum(var_diff.date_vec(:,2)==month_lims_for_mean,2)==1 & var_diff.MLD>=200;
        region_name = 'all';
        
    end
    
    
    
    % 
    
    % vertical binning
    bin_size = 25;
    var_diff.press_bins = 0:bin_size:2000;
    var_diff.press_bins = var_diff.press_bins';
    
    for cv = 1:length(comp_vars)
        var_diff.bin.(region_name).(comp_vars{cv}) = NaN(length(var_diff.press_bins),2);
        
        if strcmp(comp_vars{cv}, 'prof_DIC') ...
                || strcmp(comp_vars{cv}, 'prof_NO3') ...
                || strcmp(comp_vars{cv}, 'prof_O2')
            scale = 1000/1.0245;
        else
            scale=1;
        end
        
        for d = 1:length(var_diff.press_bins)-1
            var_depth_ind = bsose_comp(1).year.prof_depth> var_diff.press_bins(d) ...
                & bsose_comp(1).year.prof_depth<= var_diff.press_bins(d+1);
            
            temp_var = var_diff.(comp_vars{cv})(var_depth_ind,var_date_ind).*scale;
            var_diff.bin.(region_name).(comp_vars{cv})(d,1) = nanmean(reshape(temp_var,[],1));
            var_diff.bin.(region_name).(comp_vars{cv})(d,2) = (nanmean((reshape(temp_var,[],1)).^2)).^0.5; %rmse
            
        end
    end
end

%% plot
% var_filt = sum(var_diff.date_vec(:,2)==month_lims_for_mean,2)==1;

for r = 1:2
    
    if        r<3
        region_name = regions{r};
    else
        region_name = 'all';
    end
    clf
    set(gcf, 'units', 'inches')
    paper_w = 13; paper_h =6;
    set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]);
    
    for cv = 1:length(comp_vars)
        subplot(1,5,cv)
        hold on
        
        
        if strcmp(comp_vars{cv}, 'prof_DIC') ...
                || strcmp(comp_vars{cv}, 'prof_NO3') ...
                || strcmp(comp_vars{cv}, 'prof_O2')
            scale = 1000;
            units = '\mumol kg^-^1';
            
        elseif strcmp(comp_vars{cv}, 'prof_T')
            scale=1;
            units = '\circC';
        elseif strcmp(comp_vars{cv}, 'prof_S')
            scale=1;
            
            units = 'PSS';
        end
        
        if cv==1
            x_lims = [-3 2];
        elseif cv==2
            x_lims = [-0.4 .3];
        elseif cv==3
            x_lims = [-45 30];
        elseif cv==4
            x_lims = [-5 5];
        elseif cv==5
            x_lims = [-20 20];
        end
        
        
        %     var_mean = nanmean(var_diff.(comp_vars{cv})(:,var_filt).*scale,2);
        %     var_std = nanstd((var_diff.(comp_vars{cv})(:,var_filt).*scale)')';
        %
        %     boundedline( var_mean,  bsose_comp(1).year.prof_depth, var_std, 'orientation', 'horiz')
        set(gca, 'ydir', 'reverse', 'ylim', [0 200], 'xlim', x_lims)
        plot([0 0], [0 2000], '--k', 'linewidth', 1)
        
        errorbar(var_diff.bin.(region_name).(comp_vars{cv})(:,1), var_diff.press_bins+bin_size./2, ...
            var_diff.bin.(region_name).(comp_vars{cv})(:,2), 'horizontal', 'linewidth', 2)
        
        xlabel(units)
        if cv==1
        ylabel('Depth')
        end
        %calculate each profile average error in the top 200m
        
        %     all_prof = nanmean(var_diff.(comp_vars{cv})(:,var_filt).*scale,1);
        prof_avg = nanmean(var_diff.bin.(region_name).(comp_vars{cv})(var_diff.press_bins<=200,1));
        
        %     prof_rmse = (nanmean(var_diff.bin.(comp_vars{cv})(var_diff.press_bins<=200,1).^2)).^0.5
        prof_rmse = nanmean(var_diff.bin.(region_name).(comp_vars{cv})(var_diff.press_bins<=200,2));
        
        title([comp_vars{cv}(6:end) ': ' num2str(prof_avg,2) ' \pm ' num2str(prof_rmse,3)])
        set(gca, 'fontsize', 14)
        
    end
    
    
    plot_filename = ['Fig_S2_BSOSE_VIZ_prof_diff_plot_' region_name];
    
    print(gcf, '-dpdf', '-r800',  [plot_dir plot_filename '_v2.pdf'])
end
%% Crossover comparisons for Figure S3

%% pre-calculate SOCAT distance
Argo_socat_dist = [];
day_range = 1;
distance_range = 25; % km
explore_plots=0;

for f =  1:length(SO_SNs)
    
    if ~isfield(Argo.(SO_SNs{f}), 'pCO2_LIAR')
        continue
    end
    if sum(sum(~isnan(Argo.(SO_SNs{f}).pCO2_LIAR(:,1:5))))==0
        %         disp('blue')
        continue
    end
    %     Argo_socat_dist.(SO_SNs{f}).socat_dist_all = NaN(length(Argo.(SO_SNs{f}).GMT_Matlab), length(socat.latitude));
    Argo_socat_dist.(SO_SNs{f}).profile_time_match = [];
    %     Argo_socat_dist.(SO_SNs{f}).socat_dist_all = [];
    
    Argo_socat_dist.(SO_SNs{f}).mean_dist = [];
    Argo_socat_dist.(SO_SNs{f}).mean_time = [];
    Argo_socat_dist.(SO_SNs{f}).mean_pCO2 = [];
    Argo_socat_dist.(SO_SNs{f}).mean_sst = [];
    Argo_socat_dist.(SO_SNs{f}).mean_sal = [];
    Argo_socat_dist.(SO_SNs{f}).mean_lat = [];
    Argo_socat_dist.(SO_SNs{f}).mean_lon = [];
    Argo_socat_dist.(SO_SNs{f}).mean_pCO2_sstcorr = [];
    Argo_socat_dist.(SO_SNs{f}).expocode = {};
    
    for p=1:length(Argo.(SO_SNs{f}).GMT_Matlab)
        
        time_match = (socat.GMT_Matlab>= Argo.(SO_SNs{f}).GMT_Matlab(p)-day_range)  & (socat.GMT_Matlab<=  Argo.(SO_SNs{f}).GMT_Matlab(p)+day_range);
        if sum(time_match)>0
            lon_temp_argo = Argo.(SO_SNs{f}).Lon(p);
            lon_temp_argo(lon_temp_argo>180) = lon_temp_argo(lon_temp_argo>180)-360;
            
            lon_temp_socat = socat.longitude(time_match);
            lon_temp_socat(lon_temp_socat>180) = lon_temp_socat(lon_temp_socat>180)-360;
            
            socat_dist = distance('gc',Argo.(SO_SNs{f}).Lat(p), lon_temp_argo, socat.latitude(time_match), lon_temp_socat).*111.19; % distance in km
            
            if sum(socat_dist<=distance_range)
%                 disp(p)
                
                
                Argo_socat_dist.(SO_SNs{f}).mean_dist(end+1,1) = nanmean(socat_dist(socat_dist<distance_range));
                Argo_socat_dist.(SO_SNs{f}).mean_dist(end,2) = nanstd(socat_dist(socat_dist<distance_range));
                
                temp_time = socat.GMT_Matlab(time_match);
                Argo_socat_dist.(SO_SNs{f}).mean_time(end+1,1) = nanmean(temp_time(socat_dist<distance_range));
                Argo_socat_dist.(SO_SNs{f}).mean_time(end,2) = nanstd(temp_time(socat_dist<distance_range));
                
                temp_pCO2 = socat.calc_pCO2(time_match);
                Argo_socat_dist.(SO_SNs{f}).mean_pCO2(end+1,1) = nanmean(temp_pCO2(socat_dist<distance_range));
                Argo_socat_dist.(SO_SNs{f}).mean_pCO2(end,2) = nanstd(temp_pCO2(socat_dist<distance_range));
                
                temp_sst = socat.SST(time_match);
                Argo_socat_dist.(SO_SNs{f}).mean_sst(end+1,1) = nanmean(temp_sst(socat_dist<distance_range));
                Argo_socat_dist.(SO_SNs{f}).mean_sst(end,2) = nanstd(temp_sst(socat_dist<distance_range));
                
                temp_sal = socat.sal(time_match);
                Argo_socat_dist.(SO_SNs{f}).mean_sal(end+1,1) = nanmean(temp_sal(socat_dist<distance_range));
                Argo_socat_dist.(SO_SNs{f}).mean_sal(end,2) = nanstd(temp_sal(socat_dist<distance_range));
                
                temp_lat = socat.latitude(time_match);
                Argo_socat_dist.(SO_SNs{f}).mean_lat(end+1,1) = nanmean(temp_lat(socat_dist<distance_range));
                Argo_socat_dist.(SO_SNs{f}).mean_lat(end,2) = nanstd(temp_lat(socat_dist<distance_range));
                
                temp_lon = socat.longitude(time_match);
                filt_temp_lon = temp_lon(socat_dist<distance_range);
                
                if (max(filt_temp_lon) - min(filt_temp_lon)) > 300
                    filt_temp_lon(filt_temp_lon<100) = filt_temp_lon(filt_temp_lon<100) + 360;
                end
                filt_temp_lon_mean = nanmean(filt_temp_lon);
                if filt_temp_lon_mean>360
                    filt_temp_lon_mean = filt_temp_lon_mean-360;
                end
                Argo_socat_dist.(SO_SNs{f}).mean_lon(end+1,1) = nanmean(filt_temp_lon_mean);
                Argo_socat_dist.(SO_SNs{f}).mean_lon(end,2) = nanstd(filt_temp_lon);
                
                
                % recalculate the SOCAT pCO2 to account for the temperature
                % difference with the float
                ArgALK_CO2sys = NaN(length(temp_pCO2(socat_dist<distance_range)),1);
                ArgSST_CO2sys = NaN(length(temp_pCO2(socat_dist<distance_range)),1);
                Press_CO2sys = NaN(length(temp_pCO2(socat_dist<distance_range)),1);
                Press_CO2sys(:) = 5;
                
                ArgALK_CO2sys(:) = Argo.(SO_SNs{f}).TALK_LIAR(p,1);
                ArgSST_CO2sys(:) = Argo.(SO_SNs{f}).Temp_C(p,1);
                
                [DATA,~,~]= CO2SYSSOCCOM(ArgALK_CO2sys, temp_pCO2(socat_dist<distance_range) , ...
                    1,4,temp_sal(socat_dist<distance_range), temp_sst(socat_dist<distance_range), ...
                    ArgSST_CO2sys,...
                    Press_CO2sys,Press_CO2sys,20,1.8,1,10,3);
                
                Argo_socat_dist.(SO_SNs{f}).mean_pCO2_sstcorr(end+1,1) = nanmean(DATA(:,19));
                Argo_socat_dist.(SO_SNs{f}).mean_pCO2_sstcorr(end,2) = nanstd(DATA(:,19));
                
                Argo_socat_dist.(SO_SNs{f}).profile_time_match(end+1) = p;
                
                temp_expocode = socat.Expocode(time_match);
                temp_expocode = temp_expocode(socat_dist<distance_range);
                
                Argo_socat_dist.(SO_SNs{f}).expocode{end+1,1} =  unique(temp_expocode);
                
                %                 Argo_socat_dist.(SO_SNs{f}).socat_dist_all(end+1,:) = NaN(1, length(socat.latitude));
                %                 Argo_socat_dist.(SO_SNs{f}).socat_dist_all(end,time_match) = socat_dist;
                if explore_plots==1
                    clf
                    subplot(4,2,1); plot(socat_dist(socat_dist<distance_range), 'x'); hold on
                    plot(Argo_socat_dist.(SO_SNs{f}).mean_dist(end,1), 'ro');
                    
                    subplot(4,2,2); plot(temp_time(socat_dist<distance_range), 'x'); hold on
                    plot(Argo_socat_dist.(SO_SNs{f}).mean_time(end,1), 'ro');
                    
                    plot(Argo.(SO_SNs{f}).GMT_Matlab(p), 'm^')
                    
                    subplot(4,2,3); plot(temp_pCO2(socat_dist<distance_range+50), 'gs'); hold on
                    plot(DATA(:,19), 'x'); hold on;
                    plot(Argo_socat_dist.(SO_SNs{f}).mean_pCO2(end,1), 'ro');
                    plot(Argo.(SO_SNs{f}).pCO2_LIAR(p, Argo.(match_SNs{f}).Press_db(p,:)<=10), 'm^')
                    
                    subplot(4,2,4); plot(temp_sst(socat_dist<distance_range), 'x'); hold on
                    plot(Argo_socat_dist.(SO_SNs{f}).mean_sst(end,1), 'ro');
                    plot(Argo.(SO_SNs{f}).Temp_C(p, Argo.(match_SNs{f}).Press_db(p,:)<=10), 'm^')
                    
                    subplot(4,2,5); plot(temp_sal(socat_dist<distance_range), 'x'); hold on
                    plot(Argo_socat_dist.(SO_SNs{f}).mean_sal(end,1), 'ro');
                    plot(Argo.(SO_SNs{f}).Sal(p, Argo.(match_SNs{f}).Press_db(p,:)<=10), 'm^')
                    
                    subplot(4,2,6); plot(temp_lat(socat_dist<distance_range), 'x'); hold on
                    plot(Argo_socat_dist.(SO_SNs{f}).mean_lat(end,1), 'ro');
                    
                    plot(Argo.(SO_SNs{f}).Lat(p), 'm^')
                    
                    subplot(4,2,7); plot(temp_lon(socat_dist<distance_range), 'x'); hold on
                    plot(Argo_socat_dist.(SO_SNs{f}).mean_lon(end,1), 'ro');
                    
                    plot(Argo.(SO_SNs{f}).Lon(p), 'm^')
                    
                    pause
                end
                
            end
        end
    end
    

    
    disp(f)
    %         save([data_dir 'Argo_SOCAT_dist_' num2str(f) '.mat'], 'Argo_socat_dist')
    
end

clear explore_plots temp_expocode DATA ArgALK_CO2sys Press_CO2sys temp_pCO2 temp_sal temp_sst
clear temp_lon temp_lat filt_temp_lon filt_temp_lon_mean lon_temp_argo lon_temp_socat
clear clear socat_dist f time_match ArgSST_CO2sys 
%% creating the difference table
match_SNs = fieldnames(Argo_socat_dist);
diff_table = NaN(150,17);
diff_index = 0;

pCO2_type = 'pCO2_LIAR';

for f=1:length(match_SNs)
    for q = 1:length(Argo_socat_dist.(match_SNs{f}).profile_time_match)
        diff_index = diff_index+1;
        matched_profiles = Argo_socat_dist.(match_SNs{f}).profile_time_match;
        %         disp(f)
        %         disp(q)
        argo_match_depths = Argo.(match_SNs{f}).Press_db(matched_profiles(q),:)<=10;
        %                 argo_match_depths = 1:4;
        
        %         plot(Argo_socat_dist.(match_SNs{f}).mean_pCO2(q,1), nanmean(Argo.(match_SNs{f}).pCO2_LIAR(matched_profiles(q),argo_match_depths)), 's')
        %         diff_table(diff_index,1) = Argo_socat_dist.(match_SNs{f}).mean_pCO2(q,1);
        diff_table(diff_index,1) = Argo_socat_dist.(match_SNs{f}).mean_pCO2_sstcorr(q,1);
        
        diff_table(diff_index,2) =  nanmean(Argo.(match_SNs{f}).(pCO2_type)(matched_profiles(q),argo_match_depths));
        
        diff_table(diff_index,3) = Argo_socat_dist.(match_SNs{f}).mean_sst(q,1);
        diff_table(diff_index,4) = nanmean(Argo.(match_SNs{f}).Temp_C(matched_profiles(q),argo_match_depths));
        
        diff_table(diff_index,5) = Argo_socat_dist.(match_SNs{f}).mean_sal(q,1);
        diff_table(diff_index,6) = nanmean(Argo.(match_SNs{f}).Sal(matched_profiles(q),argo_match_depths));
        
        diff_table(diff_index,7) = Argo_socat_dist.(match_SNs{f}).mean_dist(q,1);
        
        diff_table(diff_index,8) = Argo_socat_dist.(match_SNs{f}).mean_time(q,1);
        diff_table(diff_index,9) = (Argo.(match_SNs{f}).GMT_Matlab(matched_profiles(q),1));
        
        diff_table(diff_index,10) = Argo_socat_dist.(match_SNs{f}).mean_lat(q,1);
        diff_table(diff_index,11) = (Argo.(match_SNs{f}).Lat(matched_profiles(q),1));
        
        diff_table(diff_index,12) = Argo_socat_dist.(match_SNs{f}).mean_lon(q,1);
        diff_table(diff_index,13) = (Argo.(match_SNs{f}).Lon(matched_profiles(q),1));
        
        if abs(diff_table(diff_index,12) - diff_table(diff_index,13))>10  && abs(diff_table(diff_index,12) - diff_table(diff_index,13))<355
            disp(f); disp(diff_index); disp(' ');
        end
        diff_table(diff_index,14) = Argo_socat_dist.(match_SNs{f}).mean_time(q,1) - Argo.(match_SNs{f}).GMT_Matlab(matched_profiles(q),1);
        
        diff_table(diff_index,15) = Argo_socat_dist.(match_SNs{f}).mean_pCO2(q,2);
        diff_table(diff_index,16) = nanstd(Argo.(match_SNs{f}).(pCO2_type)(matched_profiles(q),argo_match_depths));
        
        % float deployment
        diff_table(diff_index,16) = Argo.(match_SNs{f}).GMT_Matlab(1);
        diff_table(diff_index,17) = str2double(match_SNs{f}(2:end));
        
    end
end

clear mathed_profiles argo_match_depths f q 
%% Figure S3 
plot_filename = ['Fig_S3_SOCAT_Float_1_to_1_SOCATv2020_AB ' num2str(day_range) ' days and '  num2str(distance_range) ' km SST Corr ' pCO2_type];
plasma_map = plasma;

clf
set(gcf, 'units', 'inches')
paper_w = 12; paper_h =5;
set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]);

cmap = brewermap(8, 'Paired');
set(gcf, 'colormap', plasma_map)
d1 = subplot(1,2,1);
hold on
grid on

socat_minus_argo = diff_table(:,1) - diff_table(:,2);


plot([250 440], [250 440], '-k')
plot([250 440], [240 430], '--k')
plot([250 440], [260 450], '--k')

xlabel('SOCAT pCO_2 (\muatm)')
ylabel('Float pCO_2 (\muatm)')
set(gca, 'ylim', [250 440], 'xlim', [250 440])
d2 = subplot(1,2,2);

hold on
plot([0 0], get(gca, 'ylim'), '--k')


socat_dens = sw_dens(diff_table(:,5), diff_table(:,3), 1);
float_dens = sw_dens(diff_table(:,6), diff_table(:,4), 1);
dens_diff = abs(socat_dens - float_dens);
%
density_threshold = 0.03;
plot(d1, diff_table(dens_diff<=density_threshold,1), diff_table(dens_diff<=density_threshold,2), ...
    'marker', 's','markeredgecolor', 'k', 'markerfacecolor', 'b', 'linestyle', 'none')
hold on
title(d1, 'SOCCOM and SOCAT Matchup comparison')
% plot(d1, diff_table(dens_diff>density_threshold,1), diff_table(dens_diff>density_threshold,2), 'ro', 'markersize', 15)
diff_filtered_dens_difference = diff_table(:,1) - diff_table(:,2);
diff_filtered_dens_difference(dens_diff>density_threshold)=nan;
%
cla(gca)

hold on
h2 = histogram(diff_filtered_dens_difference);


h2.BinWidth = 5;
h2.FaceColor = cmap(4,:);
% set(gca, 'ylim', [0 20])
plot([0 0], get(gca, 'ylim'), '--k')
orig_x_lim = get(d2, 'xlim');
grid on
set(d2, 'xlim', [-max(abs(orig_x_lim)) max(abs(orig_x_lim))])
xlabel('SOCAT minus Argo (\muatm)')
title('Differences')
disp(['Green: remove > ' num2str(density_threshold) 'kg/m3 diff. Mean filtered: ' ...
    num2str(nanmean(diff_filtered_dens_difference),3) ', std: ' num2str(nanstd(diff_filtered_dens_difference),3)])
disp( ['original n: ' num2str(sum(~isnan(socat_minus_argo))) ' ' ...
    'filtered n: ' num2str(sum(~isnan(diff_filtered_dens_difference)))])

% print(gcf,'-dpdf', '-r800', [plot_dir plot_filename 'density_comp ' num2str(density_threshold) '.pdf'])


clear d2 h2 orig_x_lim d1
%% Figure S4 - bSOSE O2 vs. T w/ o2 solubility overlaid
% updating w/ new data to see if relationship still holds
% run Fig 2 and Fig 6 first to load needed data
% also fig 3 processing if you end up overlaying obs on this
plot_colors = brewermap(3, 'Set1');
plot_colors (3,:) = [0 0 0];


plot_filename = 'Fig_S4_bSOSE o2_T_v3';
clf

set(gcf, 'units', 'inches')
paper_w = 7; paper_h =7;
set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]); clear paper_w paper_h
subplot(1,1,1); hold on; grid on

x_temps = 4:.5:9;
calc_sol = GGo2_units(x_temps, 34.2, 'umol');

% l_n = NaN(size(var_names,2)+1,1);
l_n = NaN(size(reg_names,2),1);

for r = 1:size(reg_names,2)
    
    reg_index = find(cellfun(@isempty,strfind(bsose_fields, reg_names{r}))==0);
    
    y = bSOSE_anav.(bsose_fields{reg_index(3)});
    x = bSOSE_anav.(bsose_fields{reg_index(5)});
   
    
    l_n(r) = plot(x, y, plot_symbols{r}, 'linewidth',2, 'color', plot_colors(r,:));
    
    
end

l_n(r+1) = plot(o2_iav_argo.ptmp(:,1),  o2_iav_argo.ml_o2(:,1), 's', 'markerfacecolor', 'm', 'markersize', 6, 'markeredgecolor', 'k');


l_n(r+2) = plot(x_temps, calc_sol, '--k', 'linewidth', 2);
l_n(r+3) = plot(x_temps, calc_sol.*nanmean(o2_iav_argo.SLP)./1013.25, '-.', 'color', [.5 .5 .5], 'linewidth', 1.5);
set(gca, 'xlim', [4.5 8.7], 'fontsize', 13)%, 'ylim', [280 304])

xlabel('\theta (\circC)')
ylabel('[O_2] (\mumol kg^-^1)');
set(gcf, 'renderer', 'painters')
% l1 = legend(l_n, 'BSOSE C Pac', 'BSOSE SE Pac', 'BSOSE Pac', 'Argo Pac', '[O_2]_s_a_t',  ['[O_2]_s_a_t x' '$\frac{SLP}{1013.25}$'], 'location', 'northeast', 'interpreter', 'latex');
% legend('$\frac{SLP}{1013.25}$')
l1 = legend(l_n, 'BSOSE C Pac', 'BSOSE SE Pac', 'BSOSE Pac', 'Argo Pac', '$[O_2]_{sat}$', '$[O_2]_{sat}\times\frac{\overline{SLP}}{1013.25}$','Interpreter','latex', 'fontsize', 14);
print(gcf, '-dpdf', '-r800',  [plot_dir plot_filename '.pdf'])

clear l1 l_n r x y reg_index x_temps calc_sol
%% Figure S5 - Hovmouller plots

hov_props = {'DIC' 'NO3' 'O2'};
hov_prop_titles = {'[DIC] anomaly' '[NO_3^-] anomaly' '[O_2] anomaly'};
for q=1:length(hov_props)
    hov_data.(hov_props{q}) = load([cerovecki_dir 'Data_from_Ivana/2021_12_28 Data for Hovmullers/BSOSE_iter135_' hov_props{q} '_2013_2019_SouthOc_MLD_150m_AugSep_BSOSE_mask_v2.mat']);
end

% bsose_prop_names = {'MLD' 'DIC' 'Nitrate' 'Oxygen'};
color_map_names = {'PuOr' 'PiYg' 'RdBu'};
c_lims = [10 ; 2; 15];

for q=1:3

plot_filename = ['Fig_S5_bsose_hovmuller_anomalies_' hov_props{q} '_v1' ];

clf
set(gcf, 'units', 'inches')
paper_w = 7; paper_h =4;
set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]);

colormap = brewermap(30, color_map_names{q});
if q==1 || q==3
    colormap = flipud(colormap);
end

set(gcf, 'colormap', colormap)
unit='\mumol kg^-^1';

data1 = hov_data.(hov_props{q}).([hov_props{q} '_yav']);
[nx nt]=size(data1);
time=[1:nt];
X=hov_data.(hov_props{q}).XX;

nyr=floor(nt/12);
yplot=[1:12:nt];
yhov=[time(1):0.1:time(end)];

[hh]=subplot(1,1,1);
AA(q,:)=get(hh,'position');

% AA3=0.7*AA(1,3); ddd=0.02;
% AA4=AA(1,4)+ddd; dd=0.04; DDCB=0.125;
% Aleft=AA(1,1);


% figure(2); clf; LW=2.2;
FST=13; FSA=12; FSCB=11;
% subplot('position',[Aleft AA(1,2)+ddd  AA3 AA4]);
lim1=c_lims(q);
% stp=10; vcont=[-lim1:stp:lim1];
pcolor(X,time,data1'); shading interp; hold on ;

% title(['a)   DIC                                 '],'fontsize',FST)
title(hov_prop_titles{q}, 'fontsize', FST);
[cl]=colorbar; hold on ; [CC]=get(cl,'position');

% set(cl,'position',[CC(1)+DDCB CC(2) 0.6*CC(3) CC(4)])
hhll=ylabel(cl,unit,'fontsize',FSCB);
% YLH=get(hhll,'position');
% YLH2=YLH(2); YLH3=YLH(3); %set(hhll,'position',[YLH(1)-0.05 YLH2 YLH3 ])
set(gca,'YTick',[1:12:nt]); set(gca,'YTickLabel',{'2013','2014','2015','2016','2017','2018','2019'});
%
% cd(dir_my_scripts); colormap('redblue');
axis([170 290 1 nt]); ylabel([' Year'],'fontsize',FSA);
%
for itt=1:12:nt
    plot(X,squeeze(time(itt))*ones(size(X)),'linewidth',1.5,'color','k');
end; set(gca,'tickdir','out'); caxis([-1 1]*lim1) ;
set(gca,'xtick',[150:30:300]);set(gca,'xticklabel',{'150E';' 180';' 150W';'120W';'90W'}); axx=axis; set(gca,'fontsize',FSA)

print(gcf, '-dpdf', '-r800',  [plot_dir plot_filename 'w_colorbar.pdf'])
end

%% Comparison to Carter et al. 2021 properties
carter = load([cerovecki_dir 'Carter preformed properties/PreformedPropertiesDefault_v1.mat']);

depth_index = carter.Depth==200;

[carter.lon_grid, carter.lat_grid] = meshgrid(carter.Longitude, carter.Latitude);

temp_carter_o2 = carter.pO(:,:, depth_index);
temp_carter_no3 = carter.pN(:,:, depth_index);

temp_lon_carter = carter.lon_grid';
temp_lat_carter = carter.lat_grid';

% carter_lat_index = temp_lat_carter<=max(argo_MLD.LAT(1,:));
carter_lat_index_1 = temp_lat_carter<=max(argo_MLD.LAT(1,:));
carter_lat_index_2 = temp_lat_carter>=min(argo_MLD.LAT(1,:));
carter_lat_index = carter_lat_index_1 & carter_lat_index_2;

temp_lat_carter = temp_lat_carter(:,carter_lat_index(1,:));
temp_lon_carter = temp_lon_carter(:,carter_lat_index(1,:));
temp_carter_o2 = temp_carter_o2(:,carter_lat_index(1,:));
temp_carter_no3 = temp_carter_no3(:,carter_lat_index(1,:));

temp_lon_carter_shifted = NaN(size(temp_lon_carter));
temp_lat_carter_shifted = NaN(size(temp_lon_carter));
temp_o2_carter_shifted = NaN(size(temp_lon_carter));
temp_no3_carter_shifted = NaN(size(temp_lon_carter));

temp_lon_carter_shifted(1:20,:) = temp_lon_carter(341:end,:)-360;
temp_lon_carter_shifted(21:end,:) = temp_lon_carter(1:340,:);

temp_lat_carter_shifted(1:20,:) = temp_lat_carter(341:end,:);
temp_lat_carter_shifted(21:end,:) = temp_lat_carter(1:340,:);

temp_o2_carter_shifted(1:20,:) = temp_carter_o2(341:end,:);
temp_o2_carter_shifted(21:end,:) = temp_carter_o2(1:340,:);

temp_no3_carter_shifted(1:20,:) = temp_carter_no3(341:end,:);
temp_no3_carter_shifted(21:end,:) = temp_carter_no3(1:340,:);
%
set(gcf, 'colormap', turbo(50))
clf
subplot(3,2,1)
pcolor(carter.lon_grid, carter.lat_grid, carter.pO(:,:,depth_index)'); shading flat; colorbar



subplot(3,2,2)
pcolor(argo_MLD.LON, argo_MLD.LAT, argo_MLD.wint_mean); shading flat; colorbar

subplot(3,2,3)
pcolor(temp_lon_carter, temp_lat_carter, temp_carter_o2); shading flat; colorbar

subplot(3,2,4)
pcolor(temp_lon_carter_shifted, temp_lat_carter_shifted, temp_o2_carter_shifted); shading flat; colorbar


samw_index_pac = argo_MLD.wint_mean>=MLD_cutoff & argo_MLD.LON>=geo_lims.pac_lims(1) & argo_MLD.LON<=geo_lims.pac_lims(2) ...
    & argo_MLD.LAT>=geo_lims.lat_lims_pac(1) & argo_MLD.LAT<=geo_lims.lat_lims_pac(2);

temp_o2_carter_pac = temp_o2_carter_shifted;
temp_o2_carter_pac(~samw_index_pac)=nan;

temp_no3_carter_pac = temp_no3_carter_shifted;
temp_no3_carter_pac(~samw_index_pac)=nan;


subplot(3,2,5)
pcolor(temp_lon_carter_shifted, temp_lat_carter_shifted, temp_o2_carter_pac); shading flat; colorbar


samw_index_ind = argo_MLD.wint_mean>=MLD_cutoff & argo_MLD.LON>=geo_lims.ind_lims(1) & argo_MLD.LON<=geo_lims.ind_lims(2) ...
    & argo_MLD.LAT>=geo_lims.lat_lims_ind(1) & argo_MLD.LAT<=geo_lims.lat_lims_ind(2);

temp_o2_carter_ind = temp_o2_carter_shifted;
temp_o2_carter_ind(~samw_index_ind)=nan;

temp_no3_carter_ind = temp_no3_carter_shifted;
temp_no3_carter_ind(~samw_index_ind)=nan;

subplot(3,2,6)
pcolor(temp_lon_carter_shifted, temp_lat_carter_shifted, temp_o2_carter_ind); shading flat; colorbar

disp('depth')
disp(carter.Depth(depth_index))
disp('O2 Pac')
disp(nanmean(reshape(temp_o2_carter_pac,[],1)))
disp(nanstd(reshape(temp_o2_carter_pac,[],1)))

disp('O2 Ind')
disp(nanmean(reshape(temp_o2_carter_ind,[],1)))
disp(nanstd(reshape(temp_o2_carter_ind,[],1)))

disp('no3 Pac')
disp(nanmean(reshape(temp_no3_carter_pac,[],1)))
disp(nanstd(reshape(temp_no3_carter_pac,[],1)))

disp('no3 Ind')
disp(nanmean(reshape(temp_no3_carter_ind,[],1)))
disp(nanstd(reshape(temp_no3_carter_ind,[],1)))%%

