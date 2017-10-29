%% Total2netCDF_v30.m
% This function converts hourly total files from Codar format in netCDF-4
% format according to the standard struvture defined within INCREASE and
% JERICO-NEXT projects.
% Total data are set on a regular grid before to be converted in
% gridded netCDF. U_std, V_std, UVcov and GDOP are set as well.

% Current velocity measurement unit is set to m/s.

% ************
%   WARNING
% ************
% Only water U amd water V velocitites are set to m/s, errors are still
% measured in cm/s, as from the original radial Codar files.

% This version is designed for HFR_Combiner_TirLig_v21 and next releases.

% This version adds to the netCDF structure the variables std_u, std_V,
% cov, site_lat, site_lon and sites_code.

% The release 3.0 implements the data and metadata structure defined in
% INCREASE project compliant to CMEMS needs.

% INPUT:
%         mat_tot: structure containing radial file in Codar format
%         grid: two columns matrix containing the longitude values (first
%         column) and the latitude values (second column) of the
%         geographical grid where totals have been combined.
%         map_lon_lim: map longitude range (vector of two values: upper bound and lower bound).
%         map_lat_lim: map latitude range (vector of two values: upper bound and lower bound).
%         search_radius: radius where radials are collected for the
%         combination of each total.
%         Total_QC_params: structure containing parameters for total QC tests
%         mask: masking map (not used).
%         dest_nc_tot : destination folder where to save the netCDF file
%         servTH_nc_tot: THREDDS server address.
%         servRD_nc_tot: ISMAR SP disk address.

% OUTPUT:
%         T2C_err: error flag (0 = correct, 1 = error)


% Author: Lorenzo Corgnati
% Date: November 11, 2016

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [T2C_err] = Total2netCDF_v30(mat_tot, grid, map_lon_lim, map_lat_lim, search_radius, Total_QC_params, sitesLat, sitesLon, sitesCodes, mask, dest_nc_tot, servTH_nc_tot, servRD_nc_tot)

display(['[' datestr(now) '] - - ' 'Total2netCDF_v30.m started.']);

T2C_err = 0;

%% Prepares data
% Set netcdf format
ncfmt = 'netcdf4_classic';

% Sets total data on a regular grid.
try
    lonGrid = unique(mat_tot.LonLat(:,1));
    latGrid = unique(mat_tot.LonLat(:,2));
    depth = 0;
catch err
    display(['[' datestr(now) '] - - ' err.message]);
    T2C_err = 1;
end

if (T2C_err == 0)
    % Prepare variables
   
    U_grid = netcdf.getConstant('NC_FILL_DOUBLE').*ones(length(lonGrid),length(latGrid),1);
    V_grid = netcdf.getConstant('NC_FILL_DOUBLE').*ones(length(lonGrid),length(latGrid),1);
    
    U_std = netcdf.getConstant('NC_FILL_DOUBLE').*ones(length(lonGrid),length(latGrid),1);
    V_std = netcdf.getConstant('NC_FILL_DOUBLE').*ones(length(lonGrid),length(latGrid),1);
    
    covariance = netcdf.getConstant('NC_FILL_DOUBLE').*ones(length(lonGrid),length(latGrid),1);
    
    GDOP = netcdf.getConstant('NC_FILL_DOUBLE').*ones(length(lonGrid),length(latGrid),1);
    
end

if (T2C_err == 0)
    % Populate variables
    try
        for i=1:length(mat_tot.LonLat(:,1))
            lonGrid_idx = find(lonGrid==mat_tot.LonLat(i,1));
            latGrid_idx = find(latGrid==mat_tot.LonLat(i,2));
            % U and V components of current velocity
            if (not(isnan(mat_tot.U(i))))
                U_grid(lonGrid_idx,latGrid_idx,1) = mat_tot.U(i)*0.01;
            end
            if (not(isnan(mat_tot.V(i))))
                V_grid(lonGrid_idx,latGrid_idx,1) = mat_tot.V(i)*0.01;
            end
            % U and V standard errors
            if (not(isnan(mat_tot.ErrorEstimates(1,1).Uerr(i))))
                U_std(lonGrid_idx,latGrid_idx,1) = sqrt(mat_tot.ErrorEstimates(1,1).Uerr(i))*0.01;
            end
            if (not(isnan(mat_tot.ErrorEstimates(1,1).Verr(i))))
                V_std(lonGrid_idx,latGrid_idx,1) = sqrt(mat_tot.ErrorEstimates(1,1).Verr(i))*0.01;
            end
            % UV covariance
            if (not(isnan(mat_tot.ErrorEstimates(1,1).UVCovariance(i))))
                covariance(lonGrid_idx,latGrid_idx,1) = mat_tot.ErrorEstimates(1,1).UVCovariance(i)*0.0001;
            end
            % GDOP
            if (not(isnan(mat_tot.ErrorEstimates(1,1).TotalErrors(i))))
                GDOP(lonGrid_idx,latGrid_idx,1) = mat_tot.ErrorEstimates(1,1).TotalErrors(i);
            end
        end
    catch err
        display(['[' datestr(now) '] - - ' err.message]);
        T2C_err = 1;
    end
end


% Gets dimensions
if (T2C_err == 0)
    try
%        time_dim = size(mat_tot.TimeStamp,1);
        time_dim = netcdf.getConstant('unlimited');
        lat_dim = size(latGrid,1);
        lon_dim = size(lonGrid,1);
        nSites_dim = length(sitesLat);
        depth_dim = 1;
    catch err
        display(['[' datestr(now) '] - - ' err.message]);
        T2C_err = 1;
    end
end

% Sets reference time
if (T2C_err == 0)
    timeref = datenum(1970,1,1);
    time_units = ['seconds since ' datestr(timeref, 'yyyy-mm-dd') 'T' datestr(timeref, 'HH:MM:SS') 'Z'];
    [year,mon,day,hr,minutes,sec] = datevec(timeref);
%    base_stamp = sprintf('%.4d, %.2d, %.2d, %.2d',year,mon,day,hr);
%    time_units = sprintf('seconds since %.4d-%.2d-%.2d %.2d:%.2d:%.2d',year,mon,day,hr,minutes,sec);
end

% Sets data creation time and date and start and stop time and date of the
% data time window.
if (T2C_err == 0)
    try
        creation = datestr(mat_tot.TimeStamp);
        creationTime = [creation(length(creation)-7:length(creation)) ' UTC'];
        stopTime = creationTime;
        creation = datevec(mat_tot.TimeStamp);
        creationDate = [num2str(creation(1)) '-' num2str(creation(2)) '-' num2str(creation(3))];
        stopDate = creationDate;
        creation = datestr(mat_tot.TimeStamp-1/24);
        if (length(creation) == 11)
            startTime = '00:00:00 UTC';
        else
            startTime = [creation(length(creation)-7:length(creation)) ' UTC'];
        end
        creation = datevec(mat_tot.TimeStamp-1/24);
        startDate = [num2str(creation(1)) '-' num2str(creation(2)) '-' num2str(creation(3))];
    catch err
        display(['[' datestr(now) '] - - ' err.message]);
        T2C = 1;
    end
end

% Sets ADCC compliant data creation, coverage times and spatial resolution.
if (T2C_err == 0)
    try
        % File creation datetime
        dateCreated = [datestr(now, 'yyyy-mm-dd') 'T' datestr(now, 'HH:MM:SS') 'Z'];
        % Data coverage period
        coverageStart = addtodate(mat_tot.TimeStamp, -30, 'minute');
        timeCoverageStart = [datestr(coverageStart, 'yyyy-mm-dd') 'T' datestr(coverageStart, 'HH:MM:SS') 'Z'];
        coverageEnd = addtodate(mat_tot.TimeStamp, 30, 'minute');
        timeCoverageEnd = [datestr(coverageEnd, 'yyyy-mm-dd') 'T' datestr(coverageEnd, 'HH:MM:SS') 'Z'];
        % Geospatial information
        latRes = diff(latGrid);
        lonRes = diff(lonGrid);
    catch err
        display(['[' datestr(now) '] - - ' err.message]);
        T2C = 1;
    end
end

% Sets nc output file name and citation string
if (T2C_err == 0)
    try
        ts = datevec(mat_tot.TimeStamp);
        time_str = sprintf('%.4d_%.2d_%.2d_%.2d%.2d',ts(1,1),ts(1,2),ts(1,3),ts(1,4),ts(1,5));
        ncfile = [dest_nc_tot 'TOTL_' time_str '.nc'];
        servTH_ncfile = [servTH_nc_tot 'TOTL_' time_str '.nc'];
        servRD_ncfile = [servRD_nc_tot 'TOTL_' time_str '.nc'];
        citation_str = ['Data collected and processed by CNR-ISMAR within RITMARE and Jerico-Next projects -  Year ' num2str(ts(1,1))];
    catch err
        display(['[' datestr(now) '] - - ' err.message]);
        T2C_err = 1;
    end
end

% Sets collection time
if (T2C_err == 0)
    ts = datevec(mat_tot.TimeStamp);
    time_coll = [datestr(ts, 'yyyy-mm-dd') 'T' datestr(ts, 'HH:MM:SS') 'Z'];
end

% Sets the data ID
dataID = ['TirLig_' strrep(time_str(1:10), '_', '-') '_' time_str(12:13) 'Z'];

% Deletes the eventually present netCDF file with the same name
try
    delete(ncfile);
catch err
    display(['[' datestr(now) '] - - ' err.message]);
    T2C_err = 1;
end

% Builds the array containing processing parameters
procParamsVec = [single(Total_QC_params.GDOPThr), single(Total_QC_params.VelThr/100), single(2), single(3), single(0), single(search_radius), single(min(grid(:,2))), single(max(grid(:,2))), single(min(grid(:,1))), single(max(grid(:,1))), single(1.5)];

%%

%% Perform QC tests
[overall_QCflag, varianceThreshold_QCflag, GDOP_QCflag, dataDensity_QCflag, radialBalance_QCflag, velocityThreshold_QCflag] = TotalQCtests_v10(mat_tot, Total_QC_params);

%%

%% Creates vars with their dimensions
if (T2C_err == 0)
    try
        nccreate(ncfile,'TIME',...
            'Dimensions',{'TIME',time_dim},...
            'Datatype','int32',...
            'Format',ncfmt);
        
        nccreate(ncfile,'LATITUDE',...
            'Dimensions',{'LATITUDE',lat_dim},...
            'Datatype','single',...
            'Format',ncfmt);
        
        nccreate(ncfile,'LONGITUDE',...
            'Dimensions',{'LONGITUDE',lon_dim},...
            'Datatype','single',...
            'Format',ncfmt);
        
        nccreate(ncfile,'DEPH',...
            'Dimensions',{'DEPH',depth_dim},...
            'Datatype','single',...
            'Format',ncfmt);
        
        nccreate(ncfile,'EWCT',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','double',...
            'FillValue', netcdf.getConstant('NC_FILL_DOUBLE'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'NSCT',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','double',...
            'FillValue',netcdf.getConstant('NC_FILL_DOUBLE'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'EWCS',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','double',...
            'FillValue',netcdf.getConstant('NC_FILL_DOUBLE'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'NSCS',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','double',...
            'FillValue',netcdf.getConstant('NC_FILL_DOUBLE'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'CCOV',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','double',...
            'FillValue',netcdf.getConstant('NC_FILL_DOUBLE'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'GDOP',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','double',...
            'FillValue',netcdf.getConstant('NC_FILL_DOUBLE'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'QCflag',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','int16',...
            'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'VAR_QC',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','int16',...
            'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'GDOP_QC',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','int16',...
            'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'Density_QC',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','int16',...
            'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'Balance_QC',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','int16',...
            'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'CSPD_QC',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','int16',...
            'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
            'Format',ncfmt);
             
        nccreate(ncfile,'SLAT',...
            'Dimensions',{'NSIT', nSites_dim},...
            'Datatype','single',...
            'Format',ncfmt);
        
        nccreate(ncfile,'SLON',...
            'Dimensions',{'NSIT', nSites_dim},...
            'Datatype','single',...
            'Format',ncfmt);
                
        nccreate(ncfile,'SCOD',...
            'Dimensions',{'SMXL', 9,'NSIT',nSites_dim},...
            'Datatype','char',...
            'Format',ncfmt);
        
        nccreate(ncfile,'LMSK',...
            'Dimensions',{'NLMK', 1},...
            'Datatype','int8',...
            'Format',ncfmt);
        
        nccreate(ncfile,'PRPC',...
            'Dimensions',{'NPRP', 11},...
            'Datatype','double',...
            'Format',ncfmt);
        
        %% Creates attributes for the variables
        ncwriteatt(ncfile,'TIME','long_name',char('Time of Measurement UTC'));
        ncwriteatt(ncfile,'TIME','standard_name',char('time'));
%         ncwriteatt(ncfile,'TIME','base_date',char(base_stamp));
        ncwriteatt(ncfile,'TIME','units',char(time_units));
        ncwriteatt(ncfile,'TIME','calendar',char('Gregorian'));
        ncwriteatt(ncfile,'TIME','axis',char('T'));
                
        ncwriteatt(ncfile,'LATITUDE','long_name',char('Latitude'));
        ncwriteatt(ncfile,'LATITUDE','standard_name',char('latitude'));
        ncwriteatt(ncfile,'LATITUDE','units',char('degrees_north'));
        ncwriteatt(ncfile,'LATITUDE','axis',char('Y'));
        
        ncwriteatt(ncfile,'LONGITUDE','long_name',char('Longitude'));
        ncwriteatt(ncfile,'LONGITUDE','standard_name',char('longitude'));
        ncwriteatt(ncfile,'LONGITUDE','units',char('degrees_east'));
        ncwriteatt(ncfile,'LONGITUDE','axis',char('X'));
        
        ncwriteatt(ncfile,'DEPH','long_name',char('Depth of measurement'));
        ncwriteatt(ncfile,'DEPH','standard_name',char('depth'));
        ncwriteatt(ncfile,'DEPH','units',char('m'));
        ncwriteatt(ncfile,'DEPH','axis',char('Z')); 
        ncwriteatt(ncfile,'DEPH','positive',char('down'));
        ncwriteatt(ncfile,'DEPH','reference',char('sea_level'));

        ncwriteatt(ncfile,'EWCT','long_name',char('Surface Eastward Sea Water Velocity'));
        ncwriteatt(ncfile,'EWCT','standard_name',char('surface_eastward_sea_water_velocity'));
        ncwriteatt(ncfile,'EWCT','units',char('m s-1'));
        ncwriteatt(ncfile,'EWCT','scale_factor',double(1));
        ncwriteatt(ncfile,'EWCT','add_offset',double(0));
        ncwriteatt(ncfile,'EWCT','ioos_category',char('Currents'));
        ncwriteatt(ncfile,'EWCT','coordsys',char('geographic'));
%        ncwriteatt(ncfile,'EWCT','cell_methods',char('time: mean over hours time'));
        ncwriteatt(ncfile,'EWCT','valid_range',[double(-3000.0),double(3000.0)]);
%         ncwriteatt(ncfile,'EWCT','valid_min',double(-3000.0));
%         ncwriteatt(ncfile,'EWCT','valid_max',double(3000.0));
        
        ncwriteatt(ncfile,'NSCT','long_name',char('Surface Northward Sea Water Velocity'));
        ncwriteatt(ncfile,'NSCT','standard_name',char('surface_northward_sea_water_velocity'));
        ncwriteatt(ncfile,'NSCT','units',char('m s-1'));
        ncwriteatt(ncfile,'NSCT','scale_factor',double(1));
        ncwriteatt(ncfile,'NSCT','add_offset',double(0));
        ncwriteatt(ncfile,'NSCT','ioos_category',char('Currents'));
        ncwriteatt(ncfile,'NSCT','coordsys',char('geographic'));
%        ncwriteatt(ncfile,'NSCT','cell_methods',char('time: mean over hours time'));
        ncwriteatt(ncfile,'NSCT','valid_range',[double(-3000.0),double(3000.0)]);
%         ncwriteatt(ncfile,'NSCT','valid_min',double(-3000.0));
%         ncwriteatt(ncfile,'NSCT','valid_max',double(3000.0));
        
        ncwriteatt(ncfile,'EWCS','long_name',char('Standard Deviation of Surface Eastward Sea Water Velocity'));
%        ncwriteatt(ncfile,'EWCS','standard_name',char('surface_eastward_sea_water_velocity_standard_error'));
        ncwriteatt(ncfile,'EWCS','units',char('m s-1'));
        ncwriteatt(ncfile,'EWCS','valid_range',[double(-3000.0),double(3000.0)]);
%         ncwriteatt(ncfile,'EWCS','valid_min',double(-3000.0));
%         ncwriteatt(ncfile,'EWCS','valid_max',double(3000.0));
        ncwriteatt(ncfile,'EWCS','scale_factor',double(1));
        ncwriteatt(ncfile,'EWCS','add_offset',double(0));
        
        ncwriteatt(ncfile,'NSCS','long_name',char('Standard Deviation of Surface Northward Sea Water Velocity'));
%        ncwriteatt(ncfile,'NSCS','standard_name',char('surface_northward_sea_water_velocity_standard_error'));
        ncwriteatt(ncfile,'NSCS','units',char('m s-1'));
        ncwriteatt(ncfile,'NSCS','valid_range',[double(-3000.0),double(3000.0)]);
%         ncwriteatt(ncfile,'NSCS','valid_min',double(-3000.0));
%         ncwriteatt(ncfile,'NSCS','valid_max',double(3000.0));
        ncwriteatt(ncfile,'NSCS','scale_factor',double(1));
        ncwriteatt(ncfile,'NSCS','add_offset',double(0));
        
        ncwriteatt(ncfile,'CCOV','long_name',char('Covariance of Surface Sea Water Velocity'));
%         ncwriteatt(ncfile,'CCOV','standard_name',char('surface_sea_water_velocity_covariance'));
        ncwriteatt(ncfile,'CCOV','units',char('m2 s-2'));
        ncwriteatt(ncfile,'CCOV','valid_range',[double(-3000.0),double(3000.0)]);
%         ncwriteatt(ncfile,'CCOV','valid_min',double(-3000.0));
%         ncwriteatt(ncfile,'CCOV','valid_max',double(3000.0));
        ncwriteatt(ncfile,'CCOV','scale_factor',double(1));
        ncwriteatt(ncfile,'CCOV','add_offset',double(0));
        
        ncwriteatt(ncfile,'GDOP','long_name',char('Geometrical Dilution of precision'));
%         ncwriteatt(ncfile,'GDOP','standard_name',char('gdop'));
        ncwriteatt(ncfile,'GDOP','units',char('1'));
        ncwriteatt(ncfile,'GDOP','valid_range',[double(-100.0),double(100.0)]);
%         ncwriteatt(ncfile,'GDOP','valid_min',double(-100.0));
%         ncwriteatt(ncfile,'GDOP','valid_max',double(100.0));
        ncwriteatt(ncfile,'GDOP','scale_factor',double(1));
        ncwriteatt(ncfile,'GDOP','add_offset',double(0));
        ncwriteatt(ncfile,'GDOP','comment',char(['The Geometric Dilution of Precision (GDOP) is the coefficient of the uncertainty, which relates the uncertainties in radial and velocity vectors.' ...
            ' The GDOP is a unit-less coefficient, which characterizes the effect that radar station geometry has on the measurement and position determination errors.' ...
            ' A low GDOP corresponds to an optimal geometric configuration of radar stations, and results in accurate surface current data. Essentially, GDOP is a quantitative way to relate the radial and velocity vector uncertainties.'...
            ' Setting a threshold on GDOP for total combination avoids the combination of radials with an intersection angle below a certain value.' ...
            ' GDOP is a useful metric for filtering errant velocities due to poor geometry.']));
        
        ncwriteatt(ncfile,'QCflag','long_name',char('Overall Quality Flags'));
        ncwriteatt(ncfile,'QCflag','units',char('1'));
        ncwriteatt(ncfile,'QCflag','valid_range',int16([0 9]));
        ncwriteatt(ncfile,'QCflag','flag_values',int16([0 1 2 3 4 7 8 9]));
        ncwriteatt(ncfile,'QCflag','flag_meanings',char('unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value'));
        ncwriteatt(ncfile,'QCflag','comment',char('OceanSITES quality flagging for all QC tests'));
        ncwriteatt(ncfile,'QCflag','scale_factor',int16(1));
        ncwriteatt(ncfile,'QCflag','add_offset',int16(0));
        
        ncwriteatt(ncfile,'VAR_QC','long_name',char('Variance Threshold Quality flags'));
        ncwriteatt(ncfile,'VAR_QC','units',char('1'));
        ncwriteatt(ncfile,'VAR_QC','valid_range',int16([0 9]));
        ncwriteatt(ncfile,'VAR_QC','flag_values',int16([0 1 2 3 4 7 8 9]));
        ncwriteatt(ncfile,'VAR_QC','flag_meanings',char('unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value'));
        ncwriteatt(ncfile,'VAR_QC','comment',char('OceanSITES quality flagging for variance threshold QC test. Threshold set to m2/s2'));
        ncwriteatt(ncfile,'VAR_QC','scale_factor',int16(1));
        ncwriteatt(ncfile,'VAR_QC','add_offset',int16(0));
        
        ncwriteatt(ncfile,'GDOP_QC','long_name',char('GDOP Threshold Quality flags'));
        ncwriteatt(ncfile,'GDOP_QC','units',char('1'));
        ncwriteatt(ncfile,'GDOP_QC','valid_range',int16([0 9]));
        ncwriteatt(ncfile,'GDOP_QC','flag_values',int16([0 1 2 3 4 7 8 9]));
        ncwriteatt(ncfile,'GDOP_QC','flag_meanings',char('unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value'));
        ncwriteatt(ncfile,'GDOP_QC','comment',char(['OceanSITES quality flagging for GDOP threshold QC test. Threshold set to ' Total_QC_params.GDOPThr '.']));
        ncwriteatt(ncfile,'GDOP_QC','scale_factor',int16(1));
        ncwriteatt(ncfile,'GDOP_QC','add_offset',int16(0));
        
        ncwriteatt(ncfile,'Density_QC','long_name',char('Data density Threshold Quality flags'));
        ncwriteatt(ncfile,'Density_QC','units',char('1'));
        ncwriteatt(ncfile,'Density_QC','valid_range',int16([0 9]));
        ncwriteatt(ncfile,'Density_QC','flag_values',int16([0 1 2 3 4 7 8 9]));
        ncwriteatt(ncfile,'Density_QC','flag_meanings',char('unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value'));
        ncwriteatt(ncfile,'Density_QC','comment',char('OceanSITES quality flagging for Data density threshold QC test. Threshold set to contributing radials.'));
        ncwriteatt(ncfile,'Density_QC','scale_factor',int16(1));
        ncwriteatt(ncfile,'Density_QC','add_offset',int16(0));
        
        ncwriteatt(ncfile,'Balance_QC','long_name',char('Balance of number of radials from each contributing site Quality flags'));
        ncwriteatt(ncfile,'Balance_QC','units',char('1'));
        ncwriteatt(ncfile,'Balance_QC','valid_range',int16([0 9]));
        ncwriteatt(ncfile,'Balance_QC','flag_values',int16([0 1 2 3 4 7 8 9]));
        ncwriteatt(ncfile,'Balance_QC','flag_meanings',char('unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value'));
        ncwriteatt(ncfile,'Balance_QC','comment',char('OceanSITES quality flagging for Balance of number of radials from each contributing site QC test. Balance threshold set to .'));
        ncwriteatt(ncfile,'Balance_QC','scale_factor',int16(1));
        ncwriteatt(ncfile,'Balance_QC','add_offset',int16(0));
        
        ncwriteatt(ncfile,'CSPD_QC','long_name',char('Velocity threshold Quality flags'));
        ncwriteatt(ncfile,'CSPD_QC','units',char('1'));
        ncwriteatt(ncfile,'CSPD_QC','valid_range',int16([0 9]));
        ncwriteatt(ncfile,'CSPD_QC','flag_values',int16([0 1 2 3 4 7 8 9]));
        ncwriteatt(ncfile,'CSPD_QC','flag_meanings',char('unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value'));
        ncwriteatt(ncfile,'CSPD_QC','comment',char('OceanSITES quality flagging for Velocity threshold QC test. Threshold set to m/s.'));
        ncwriteatt(ncfile,'CSPD_QC','scale_factor',int16(1));
        ncwriteatt(ncfile,'CSPD_QC','add_offset',int16(0));
        
        ncwriteatt(ncfile,'SLAT','long_name',char('Contributing radar site latitudes'));
        ncwriteatt(ncfile,'SLAT','standard_name',char('latitude'));
        ncwriteatt(ncfile,'SLAT','units','degrees_north');
        
        ncwriteatt(ncfile,'SLON','long_name',char('Contributing radar site longitudes'));
        ncwriteatt(ncfile,'SLON','standard_name',char('longitude'));
        ncwriteatt(ncfile,'SLON','units','degrees_east');
        
        ncwriteatt(ncfile,'SCOD','long_name',char('Contributing radar site codes'));
        ncwriteatt(ncfile,'SCOD','comment',char(['RDLI: Ideal Antenna Pattern is applied to the radial data.'...
            ' RDLM: Measured Antenna Pattern is applied to the radial data.']));
        
        ncwriteatt(ncfile,'LMSK','long_name',char('Flag for masking dry grid points'));
%         ncwriteatt(ncfile,'LMSK','standard_name',char('land_mask'));
        ncwriteatt(ncfile,'LMSK','units',char('1'));
        ncwriteatt(ncfile,'LMSK','comment',char(['Flag = 0: masking inactive' ...
            ' Flag = 1: masking active']));
        
        ncwriteatt(ncfile,'PRPC','long_name',char('Near Real Time data processing parameters'));
%         ncwriteatt(ncfile,'PRPC','standard_name',char('processing_parameters'));
        ncwriteatt(ncfile,'PRPC','units',char('1'));
        ncwriteatt(ncfile,'PRPC','comment',char(['01) Maximum GDOP threshold'...
            ' 02) Maximum velocity threshold [m/s]' ...
            ' 03) Minimum number of sites required'...
            ' 04) Minimum number of radials required' ...
            ' 05) Maximum angular gap to interpolate radials [degrees_true] (0 = no interpolation)' ...
            ' 06) Spatial search radius for radial solutions [km]' ...
            ' 07) Minimum latitude for total combination [degrees_north]' ...
            ' 08) Maximum latitude for total combination [degrees_north]' ...
            ' 09) Minimum longitude for total combination [degrees_east]' ...
            ' 10) Maximum longitude for total combination [degrees_east]' ...
            ' 11) Spatial Resolution [km]']));
        
        %% Writes values in variables
        ncwrite(ncfile,'TIME',int32((mat_tot.TimeStamp-timeref)*86400));
        ncwrite(ncfile,'LATITUDE',latGrid);
        ncwrite(ncfile,'LONGITUDE',lonGrid);
        ncwrite(ncfile,'DEPH',depth);
        ncwrite(ncfile,'EWCT',U_grid);
        ncwrite(ncfile,'NSCT',V_grid);
        ncwrite(ncfile,'EWCS',U_std);
        ncwrite(ncfile,'NSCS',V_std);
        ncwrite(ncfile,'CCOV',covariance);
        ncwrite(ncfile,'GDOP',GDOP);
        ncwrite(ncfile,'SLAT',sitesLat);
        ncwrite(ncfile,'SLON',sitesLon);
        ncwrite(ncfile,'SCOD',sitesCodes');
        ncwrite(ncfile,'LMSK',double(1));
        ncwrite(ncfile,'PRPC',procParamsVec);
        ncwrite(ncfile,'QCflag',overall_QCflag);
        ncwrite(ncfile,'VAR_QC',varianceThreshold_QCflag);
        ncwrite(ncfile,'GDOP_QC',GDOP_QCflag);
        ncwrite(ncfile,'Density_QC',dataDensity_QCflag);
        ncwrite(ncfile,'Balance_QC',radialBalance_QCflag);
        ncwrite(ncfile,'CSPD_QC',velocityThreshold_QCflag);
        
        %% Defines global attributes
        ncwriteatt(ncfile,'/','title',char('Near Real Time Surface Ocean Velocity'));
        ncwriteatt(ncfile,'/','institution',char('National Research Council - Institute of Marine Science, U.O.S. La Spezia'));
        ncwriteatt(ncfile,'/','Conventions',char('CF-1.6, Unidata, OceanSITES, ACDD, INSPIRE'));
        ncwriteatt(ncfile,'/','summary',char('The data set consists of maps of total velocity of the surface current in the North-Western Tyrrhenian Sea and Ligurian Sea averaged over a time interval of 1 hour around the cardinal hour. Surface ocean velocities estimated by HF Radar are representative of the upper 0.3-2.5 meters of the ocean. The main objective of near real time processing is to produce the best product from available data at the time of processing. Total velocities are derived using least square fit that maps radial velocities measured from individual sites onto a cartesian grid. The final product is a map of the horizontal components of the ocean currents on a regular grid in the area of overlap of two or more radar stations.'));
        ncwriteatt(ncfile,'/','source',char('coastal structure'));
        ncwriteatt(ncfile,'/','network',char('HF Radar Surface Observation'));
        ncwriteatt(ncfile,'/','keywords',char('OCEAN CURRENTS, SURFACE WATER, RADAR, SCR-HF'));
        ncwriteatt(ncfile,'/','keywords_vocabulary',char('GCMD Science Keywords'));
        ncwriteatt(ncfile,'/','history',char([time_coll ' data collected. ' dateCreated ' netCDF file created and sent to TAC']));
        ncwriteatt(ncfile,'/','data_language',char('eng'));
        ncwriteatt(ncfile,'/','data_character_set',char('utf8'));
        ncwriteatt(ncfile,'/','topic_category',char('oceans'));
        ncwriteatt(ncfile,'/','reference_system',char('EPSG:4806'));
        ncwriteatt(ncfile,'/','metadata_language',char('eng'));
        ncwriteatt(ncfile,'/','metadata_character_set',char('utf8'));
        ncwriteatt(ncfile,'/','metadata_contact',char('lorenzo.corgnati@sp.ismar.cnr.it'));
        ncwriteatt(ncfile,'/','metadata_date_stamp',char(dateCreated)); 
        ncwriteatt(ncfile,'/','abstract',char('The data set consists of maps of total velocity of the surface current in the North-Western Tyrrhenian Sea and Ligurian Sea averaged over a time interval of 1 hour around the cardinal hour. Surface ocean velocities estimated by HF Radar are representative of the upper 0.3-2.5 meters of the ocean. Total velocities are derived using least square fit that maps radial velocities measured from individual sites onto a cartesian grid. The final product is a map of the horizontal components of the ocean currents on a regular grid in the area of overlap of two or more radar stations.'));
        ncwriteatt(ncfile,'/','standard_name_vocabulary',char('NetCDF Climate and Forecast (CF) Metadata Convention Standard Name Table Version 1.6'));
        ncwriteatt(ncfile,'/','id',char(dataID));
        ncwriteatt(ncfile,'/','naming_authority',char('it.cnr.ismar'));
        ncwriteatt(ncfile,'/','cdm_data_type',char('Grid'));
        ncwriteatt(ncfile,'/','project',char('RITMARE and Jerico-Next'));
        ncwriteatt(ncfile, '/','time_coverage_start',char(timeCoverageStart));
        ncwriteatt(ncfile, '/','time_coverage_end',char(timeCoverageEnd));
        ncwriteatt(ncfile, '/','time_coverage_duration',char('PT1H'));
        ncwriteatt(ncfile, '/','time_coverage_resolution',char('PT1H'));
        ncwriteatt(ncfile,'/','geospatial_lat_min',char(num2str(min(map_lat_lim))));
        ncwriteatt(ncfile,'/','geospatial_lat_max',char(num2str(max(map_lat_lim))));      
        ncwriteatt(ncfile,'/','geospatial_lon_min',char(num2str(min(map_lon_lim))));
        ncwriteatt(ncfile,'/','geospatial_lon_max',char(num2str(max(map_lon_lim))));
        ncwriteatt(ncfile,'/','date_created',char(dateCreated));
        ncwriteatt(ncfile,'/','processing_level',char('3B'));
        ncwriteatt(ncfile,'/','license',char('HF radar sea surface current velocity dataset by CNR-ISMAR is licensed under a Creative Commons Attribution 4.0 International License. You should have received a copy of the license along with this work. If not, see http://creativecommons.org/licenses/by/4.0/.'));
        ncwriteatt(ncfile,'/','creator_name',char('Lorenzo Corgnati'));
        ncwriteatt(ncfile,'/','creator_url',char('http://radarhf.ismar.cnr.it'));
        ncwriteatt(ncfile,'/','creator_email',char('lorenzo.corgnati@sp.ismar.cnr.it'));
        ncwriteatt(ncfile,'/','acknowledgment',char('ISMAR HF Radar Network has been established within RITMARE and Jerico-Next projects. The network has been designed, implemented and managed through the efforts of ISMAR UOS La Spezia.'));
        ncwriteatt(ncfile,'/','comment',char('Total velocities are derived using least square fit that maps radial velocities measured from individual sites onto a cartesian grid. The final product is a map of the horizontal components of the ocean currents on a regular grid in the area of overlap of two or more radar stations.'));
        ncwriteatt(ncfile,'/','netcdf_version',char(netcdf.inqLibVers));
        ncwriteatt(ncfile,'/','netcdf_format',char(ncfmt));
        ncwriteatt(ncfile,'/','metadata_convention',char('Unidata Dataset Discovery v1.0 compliant. NOAA GNOME format compliant.'));
        ncwriteatt(ncfile,'/','platform',char('CNR-ISMAR HF Radar Network'));
        ncwriteatt(ncfile,'/','sensor',char('CODAR SeaSonde'));
        ncwriteatt(ncfile,'/','date_modified',char(dateCreated));
        ncwriteatt(ncfile,'/','grid_resolution',char('1.5 Km'));
        ncwriteatt(ncfile,'/','grid_mapping',char('Transverse Mercator'));
        ncwriteatt(ncfile,'/','citation',char(citation_str));        
        ncwriteatt(ncfile,'/','institution_reference',char('http://www.ismar.cnr.it'));
        ncwriteatt(ncfile,'/','operational_manager',char('Carlo Mantovani'));
        ncwriteatt(ncfile,'/','operational_manager_email',char('carlo.mantovani@cnr.it'));
        ncwriteatt(ncfile,'/','format_version',char('v1.0'));
        ncwriteatt(ncfile,'/','data_mode',char('R'));
        ncwriteatt(ncfile,'/','update_interval',char('void'));
        ncwriteatt(ncfile,'/','site_code',char(''));
        ncwriteatt(ncfile,'/','area',char('Mediterranean Sea'));
        ncwriteatt(ncfile,'/','regional_description',char('North-Western Tyrrhenian Sea and Ligurian Sea, Italy'));
        ncwriteatt(ncfile,'/','date_issued',char(dateCreated));
        ncwriteatt(ncfile,'/','software_name',char('HFR_Combiner_TirLig'));
        ncwriteatt(ncfile,'/','software_version',char('v2.1'));
        ncwriteatt(ncfile,'/','quality_control',char('Level-B: Velocity Threshold + GDOP Threshold'));
        ncwriteatt(ncfile,'/','file_quality_index',int32(0));
        ncwriteatt(ncfile,'/','interpolation',char(''));
        ncwriteatt(ncfile,'/','geospatial_lat_resolution',latRes(1));
        ncwriteatt(ncfile,'/','geospatial_lon_resolution',lonRes(1));
        ncwriteatt(ncfile,'/','geospatial_lat_units',char('degrees_north'));
        ncwriteatt(ncfile,'/','geospatial_lon_units',char('degrees_east'));
        ncwriteatt(ncfile,'/','geospatial_vertical_max', char('0'));
        ncwriteatt(ncfile,'/','geospatial_vertical_min', char('1'));
        ncwriteatt(ncfile,'/','geospatial_vertical_positive', char('down'));
        ncwriteatt(ncfile,'/','geospatial_vertical_resolution', char('1'));
        ncwriteatt(ncfile,'/','geospatial_vertical_units', char('m'));
        ncwriteatt(ncfile,'/','references',char('HFR_Progs Matlab Documentation - Copyright (C) 2006-7 David M. Kaplan; Otero,M. (2008).NETCDF DESCRIPTION FOR NEAR REAL-TIME SURFACE CURRENTS PRODUCED BY THE HF-RADAR NETWORK. https://cordc.ucsd.edu/projects/mapping/documents/HFRNet_RTV-NetCDF.pdf'));
        ncwriteatt(ncfile,'/','publisher_name',char('Lorenzo Corgnati'));
        ncwriteatt(ncfile,'/','publisher_url',char('http://radarhf.ismar.cnr.it'));
        ncwriteatt(ncfile,'/','publisher_email',char('lorenzo.corgnati@sp.ismar.cnr.it'));
        ncwriteatt(ncfile,'/','contributor_name',char('Vega Forneris, Cristina Tronconi'));
        ncwriteatt(ncfile,'/','contributor_role',char('THREDDS expert, metadata expert'));
      
        
%        ncwriteatt(ncfile,'/','summary',char('The data set consists of maps of total velocity of the surface current in the North-Western Tyrrhenian Sea and Ligurian Sea averaged over a time interval of 1 hour around the cardinal hour. Surface ocean velocities estimated by HF Radar are representative of the upper 0.3-2.5 meters of the ocean. The main objective of near real time processing is to produce the best product from available data at the time of processing. Total velocities are derived using least square fit that maps radial velocities measured from individual sites onto a cartesian grid. The final product is a map of the horizontal components of the ocean currents on a regular grid in the area of overlap of two or more radar stations.'));       
        
%        ncwriteatt(ncfile,'/','creation_date',char(creationDate));
%        ncwriteatt(ncfile,'/','creation_time',char(creationTime));
        
    catch err
        display(['[' datestr(now) '] - - ' err.message]);
        T2C_err = 1;
    end
end

% Copies the netCDF file to the THREDDS server
if (T2C_err == 0)
    try
        [cp_status] = copyfile(ncfile,servTH_ncfile,'f');
    catch err
        display(['[' datestr(now) '] - - ' err.message]);
        T2C_err = 1;
    end
end

% Copies the netCDF file to the RadarDisk
if (T2C_err == 0)
    try
        [cp_status] = copyfile(ncfile,servRD_ncfile,'f');
    catch err
        display(['[' datestr(now) '] - - ' err.message]);
        T2C_err = 1;
    end
end
return

