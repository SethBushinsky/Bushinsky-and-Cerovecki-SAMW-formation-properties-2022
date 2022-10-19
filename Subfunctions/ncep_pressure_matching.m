function [press_int, NCEP_pressures, pressure_time_matlab, local_pressure] = ncep_pressure_matching(GMT_matlab_dates, lat, lon, NCEP_pressures)
% 2019_02_25 - make sure that input lons are matched to 0-360 to match the
% NCEP pressure grid
% v2016, corrected to allow lon to cross the international date line.  
% Seth Bushinsky, smb4@hawaii.edu
if isempty(NCEP_pressures), %#ok<*NOCOL>
    [NCEP_pressures] = ncep_pressure_load;
end

if~isnumeric(lat),
    lat = str2double(lat);
    lon = str2double(lon);
end

% NCEP pressure is 0-360.  Change input lons to match.
lon_neg_index = find(lon<0);
lon(lon_neg_index) = lon(lon_neg_index) + 360;


if (size(GMT_matlab_dates,1)>size(lat,1)) || (size(GMT_matlab_dates,2)>size(lat,2)),
    lat_matched = zeros(size(GMT_matlab_dates));
    lat_matched(:) = lat;
    
    lon_matched = zeros(size(GMT_matlab_dates));
    lon_matched(:) = lon;
else
    lat_matched = lat;
    lon_matched = lon;
end

press_int = zeros(size(GMT_matlab_dates));
press_int(:)=nan;
pressure_time_matlab = NCEP_pressures.time;

for z = 1:length(GMT_matlab_dates),
    lat_index = find(NCEP_pressures.lat>lat_matched(z)-2 & NCEP_pressures.lat<lat_matched(z)+2);
    if ~isempty(lat_index)
        
        
        if size(lat_index,1)>1,
            a = lat_matched(z) - NCEP_pressures.lat(lat_index);
            if eq(abs(a(1)), abs(a(2)))  % for that odd case when the lat is exactly in between grid points
                lat_index = lat_index(1);
            else
                b = find(min(abs(a))==abs(a));
                lat_index = lat_index(b);
            end
        end
        
        if lon_matched(z)<0,
            lon_matched(z) = lon_matched(z)+360;
        end
        
        lon_index = find(NCEP_pressures.lon>=lon_matched(z)-2 & NCEP_pressures.lon<lon_matched(z)+2);
        
        if size(lon_index,1)>1,
            c = lon_matched(z) - NCEP_pressures.lon(lon_index);
            if eq(abs(c(1)), abs(c(2)))  % for that odd case when the lon is exactly in between grid points
                lon_index = lon_index(1);
            else
                d = find(min(abs(c))==abs(c));
                lon_index = lon_index(d);
            end
        elseif lon_matched(z)>=357.5,
            c = lon_matched(z) - [NCEP_pressures.lon(144); 360];
            d = find(min(abs(c))==abs(c));
            if d==1,
                lon_index = 144;
            elseif d==2,
                lon_index = 1;
            end
        end
        % disp(z)
        local_pressure = (squeeze(NCEP_pressures.pressure(lon_index, lat_index, :)))./100; %converts to hpa
        press_int(z) = interp1(pressure_time_matlab, local_pressure, GMT_matlab_dates(z));
    end
end
if ~exist('local_pressure', 'var'),
    local_pressure=[];
end
end

function [NCEP_pressures] =  ncep_pressure_load

pressure_directory = [data_dir 'Data_Products/SLP/NCEPNCAR/'];
try
    load([pressure_directory 'ncep_pressure.mat']);
catch %err
    %     disp(err)
    NCEP_pressures.time_stamp = 0;
end
if NCEP_pressures.time_stamp<now-14,
    
    pressure = [];
    pressure_time = [];
    pressure_lat = [];
    pressure_lon = [];
    pressure_temp = [];
    
    files = dir([pressure_directory 'pres.sfc.*']);
    for f = 1:length(files),
        if f==1, % only need to make lat/lon once
            pressure_lat = ncread([pressure_directory files(f).name], 'lat');
            pressure_lon = ncread([pressure_directory files(f).name], 'lon');
        end
        pressure_temp{f} = ncread([pressure_directory files(f).name], 'pres');
        pressure_time_temp{f} = ncread([pressure_directory files(f).name], 'time');
    end
    for f = 1:length(files),
        pressure = cat(3, pressure, pressure_temp{f});
        pressure_time = cat(1, pressure_time, pressure_time_temp{f});
    end
    
    time_stamp = now;
    
    pressure_time_matlab = pressure_time./24+datenum('1-1-1800 00:00:00');
    
    NCEP_pressures.pressure = pressure;
    NCEP_pressures.time = pressure_time_matlab;
    NCEP_pressures.lat = pressure_lat;
    NCEP_pressures.lon = pressure_lon;
    NCEP_pressures.time_stamp = time_stamp;
    
    clear pressure_temp pressure_time_temp
    save([pressure_directory 'ncep_pressure'], 'NCEP_pressures')
end
end