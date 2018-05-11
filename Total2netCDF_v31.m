%% Total2netCDF_v31.m
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

% The release 3.1 implements the data and metadata structure defined in
% INCREASE project and suited for CMEMS-INSTAC requirements.

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
%         sitesCodes: codes of the contributing radar sites.
%         mask: masking map (not used).
%         dest_nc_tot : destination folder where to save the netCDF file
%         servRD_nc_tot: ISMAR SP disk address.
%         lastPatternStr: string containing the last pattern measurement date

% OUTPUT:
%         T2C_err: error flag (0 = correct, 1 = error)


% Author: Lorenzo Corgnati
% Date: March 15, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [T2C_err] = Total2netCDF_v31(mat_tot, grid, map_lon_lim, map_lat_lim, search_radius, Total_QC_params, sitesLat, sitesLon, sitesCodes, mask, dest_nc_tot, servRD_nc_tot, lastPatternStr)

display(['[' datestr(now) '] - - ' 'Total2netCDF_v31.m started.']);

T2C_err = 0;

warning('off', 'all');

%% Prepare data
% Set netcdf format
ncfmt = 'netcdf4_classic';

% Set total data on a regular grid.
try
    lonGrid = unique(mat_tot.LonLat(:,1));
    latGrid = unique(mat_tot.LonLat(:,2));
    depth = 0;
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

if (T2C_err == 0)
    % Prepare variables
   
    mat_tot.U_grid = netcdf.getConstant('NC_FILL_DOUBLE').*ones(length(lonGrid),length(latGrid),1);
    mat_tot.V_grid = netcdf.getConstant('NC_FILL_DOUBLE').*ones(length(lonGrid),length(latGrid),1);
    
    mat_tot.U_std = netcdf.getConstant('NC_FILL_DOUBLE').*ones(length(lonGrid),length(latGrid),1);
    mat_tot.V_std = netcdf.getConstant('NC_FILL_DOUBLE').*ones(length(lonGrid),length(latGrid),1);
    
    mat_tot.covariance = netcdf.getConstant('NC_FILL_DOUBLE').*ones(length(lonGrid),length(latGrid),1);
    
    mat_tot.GDOP = netcdf.getConstant('NC_FILL_DOUBLE').*ones(length(lonGrid),length(latGrid),1);
    
end

if (T2C_err == 0)
    % Populate variables
    try
        for i=1:length(mat_tot.LonLat(:,1))
            lonGrid_idx = find(lonGrid==mat_tot.LonLat(i,1));
            latGrid_idx = find(latGrid==mat_tot.LonLat(i,2));
            % U and V components of current velocity
            if (not(isnan(mat_tot.U(i))))
                mat_tot.U_grid(lonGrid_idx,latGrid_idx,1) = mat_tot.U(i)*0.01;
            end
            if (not(isnan(mat_tot.V(i))))
                mat_tot.V_grid(lonGrid_idx,latGrid_idx,1) = mat_tot.V(i)*0.01;
            end
            % U and V standard errors
            if (not(isnan(mat_tot.ErrorEstimates(1,1).Uerr(i))))
                mat_tot.U_std(lonGrid_idx,latGrid_idx,1) = sqrt(mat_tot.ErrorEstimates(1,1).Uerr(i))*0.01;
            end
            if (not(isnan(mat_tot.ErrorEstimates(1,1).Verr(i))))
                mat_tot.V_std(lonGrid_idx,latGrid_idx,1) = sqrt(mat_tot.ErrorEstimates(1,1).Verr(i))*0.01;
            end
            % UV covariance
            if (not(isnan(mat_tot.ErrorEstimates(1,1).UVCovariance(i))))
                mat_tot.covariance(lonGrid_idx,latGrid_idx,1) = mat_tot.ErrorEstimates(1,1).UVCovariance(i)*0.0001;
            end
            % GDOP
            if (not(isnan(mat_tot.ErrorEstimates(1,1).TotalErrors(i))))
                mat_tot.GDOP(lonGrid_idx,latGrid_idx,1) = mat_tot.ErrorEstimates(1,1).TotalErrors(i);
            end
        end
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        T2C_err = 1;
    end
end


% Gets dimensions
if (T2C_err == 0)
    try
        time_dim = size(mat_tot.TimeStamp,1);
%         time_dim = netcdf.getConstant('unlimited');
        lat_dim = size(latGrid,1);
        lon_dim = size(lonGrid,1);
        depth_dim = 1;
        maxSite_dim = 50;
        string15_dim = 15;
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        T2C_err = 1;
    end
end

% Set reference time
% if (T2C_err == 0)
%     timeref = datenum(1970,1,1);
%     time_units = ['seconds since ' datestr(timeref, 'yyyy-mm-dd') 'T' datestr(timeref, 'HH:MM:SS') 'Z'];
%     [year,mon,day,hr,minutes,sec] = datevec(timeref);
% end

if (T2C_err == 0)
    timeref = datenum(1950,1,1);
    time_units = ['days since ' datestr(timeref, 'yyyy-mm-dd') 'T' datestr(timeref, 'HH:MM:SS') 'Z'];
    [year,mon,day,hr,minutes,sec] = datevec(timeref);
end

% Set data creation time and date and start and stop time and date of the
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
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        T2C = 1;
    end
end

% Set ADCC compliant data creation, coverage times and spatial resolution.
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
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        T2C = 1;
    end
end

% Set nc output file name, citation string and distribution string
if (T2C_err == 0)
    try
        ts = datevec(mat_tot.TimeStamp);
        time_str = sprintf('%.4d_%.2d_%.2d_%.2d%.2d',ts(1,1),ts(1,2),ts(1,3),ts(1,4),ts(1,5));
        ncfile = [dest_nc_tot 'TOTL_' time_str '.nc'];
        servRD_ncfile = [servRD_nc_tot 'TOTL_' time_str '.nc'];
        citation_str = ['These data were collected and made freely available by the Copernicus project and the programs that contribute to it. Data collected and processed by CNR-ISMAR within RITMARE and Jerico-Next projects -  Year ' num2str(ts(1,1))];
        distribution_str = 'These data follow Copernicus standards; they are public and free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.';
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        T2C_err = 1;
    end
end

% Sets collection time
if (T2C_err == 0)
    ts = datevec(mat_tot.TimeStamp);
    time_coll = [datestr(ts, 'yyyy-mm-dd') 'T' datestr(ts, 'HH:MM:SS') 'Z'];
end

% Sets the data ID
dataID = ['HFR_TirLig_Total_' strrep(time_str(1:10), '_', '-') '_' time_str(12:13) 'Z'];

% Deletes the eventually present netCDF file with the same name
try
    delete(ncfile);
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    T2C_err = 1;
end

% Builds the array containing processing parameters
procParamsVec = [single(Total_QC_params.GDOPThr), single(Total_QC_params.VelThr/100), single(2), single(3), single(0), single(search_radius), single(min(grid(:,2))), single(max(grid(:,2))), single(min(grid(:,1))), single(max(grid(:,1))), single(1.5)];

%%

%% Perform QC tests

% Build the names of the files of the previous two hours
[twoHoursBefore, oneHourBefore] = TimeStamp2TS(mat_tot.TimeStamp);
twoHoursBefore_rep1 = strrep(ncfile,time_str,twoHoursBefore);
Total_QC_params.TempDerThr.hour2 = strrep(twoHoursBefore_rep1,time_str(1:length(time_str)-5),twoHoursBefore(1:length(twoHoursBefore)-5));
oneHoursBefore_rep1 = strrep(ncfile,time_str,oneHourBefore);
Total_QC_params.TempDerThr.hour1 = strrep(oneHoursBefore_rep1,time_str(1:length(time_str)-5),oneHourBefore(1:length(oneHourBefore)-5));
oneHoursBefore_rep1_rd = strrep(servRD_ncfile,time_str,oneHourBefore);
oneHoursBefore_rep2_rd = strrep(oneHoursBefore_rep1_rd,time_str(1:length(time_str)-5),oneHourBefore(1:length(oneHourBefore)-5));
oneHoursBefore_rep3_rd = strrep(oneHoursBefore_rep2_rd,time_str(1:length(time_str)-8),oneHourBefore(1:length(oneHourBefore)-8));
Total_QC_params.TempDerThr.hour1_RD = strrep(oneHoursBefore_rep3_rd,time_str(1:length(time_str)-11),oneHourBefore(1:length(oneHourBefore)-11));

[overall_QCflag, varianceThreshold_QCflag, temporalDerivativeThreshold_QCflag, GDOP_QCflag, dataDensity_QCflag, velocityThreshold_QCflag] = TotalQCtests_v10(mat_tot, Total_QC_params);

%%

%% Set the time, position and depth quality flags
% Time quality flag
sdnTime_QCflag = 1;
% Position quality flag
sdnPosition_QCflag = netcdf.getConstant('NC_FILL_SHORT').*int16(ones(length(lonGrid),length(latGrid),1));
sdnPosition_QCflag(mat_tot.U_grid~=netcdf.getConstant('NC_FILL_DOUBLE')) = 1;

% Depth quality flag
sdnDepth_QCflag = 1;

%%

%% Creates vars with their dimensions
if (T2C_err == 0)
    try
        nccreate(ncfile,'TIME',...
            'Dimensions',{'TIME',time_dim},...
            'Datatype','single',...
            'Format',ncfmt);
        
        nccreate(ncfile,'LATITUDE',...
            'Dimensions',{'LATITUDE',lat_dim},...
            'Datatype','single',...
            'Format',ncfmt);
        
        nccreate(ncfile,'LONGITUDE',...
            'Dimensions',{'LONGITUDE',lon_dim},...
            'Datatype','single',...
            'Format',ncfmt);
        
        nccreate(ncfile,'crs',...
            'Dimensions',{'crs',time_dim},...
            'Datatype','int16',...
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
        
        nccreate(ncfile,'TIME_SEADATANET_QC',...
            'Dimensions',{'TIME',time_dim},...
            'Datatype','int16',...
            'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'POSITION_SEADATANET_QC',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','int16',...
            'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'DEPTH_SEADATANET_QC',...
            'Dimensions',{'TIME',time_dim},...
            'Datatype','int16',...
            'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'QCflag',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','int16',...
            'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'VART_QC',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','int16',...
            'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'GDOP_QC',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','int16',...
            'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'DDNS_QC',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','int16',...
            'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'CSPD_QC',...
            'Dimensions',{'LONGITUDE',lon_dim,'LATITUDE',lat_dim, 'DEPH', depth_dim, 'TIME',time_dim},...
            'Datatype','int16',...
            'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'NARX',...
            'Dimensions',{'TIME',time_dim},...
            'Datatype','int16',...
            'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'NATX',...
            'Dimensions',{'TIME',time_dim},...
            'Datatype','int16',...
            'FillValue',netcdf.getConstant('NC_FILL_SHORT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'SLTR',...
            'Dimensions',{'MAXSITE',maxSite_dim,'TIME',time_dim},...
            'Datatype','single',...
            'FillValue',netcdf.getConstant('NC_FILL_FLOAT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'SLNR',...
            'Dimensions',{'MAXSITE',maxSite_dim,'TIME',time_dim},...
            'Datatype','single',...
            'FillValue',netcdf.getConstant('NC_FILL_FLOAT'),...
            'Format',ncfmt);
             
        nccreate(ncfile,'SLTT',...
            'Dimensions',{'MAXSITE',maxSite_dim,'TIME',time_dim},...
            'Datatype','single',...
            'FillValue',netcdf.getConstant('NC_FILL_FLOAT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'SLNT',...
            'Dimensions',{'MAXSITE',maxSite_dim,'TIME',time_dim},...
            'Datatype','single',...
            'FillValue',netcdf.getConstant('NC_FILL_FLOAT'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'SCDR',...
            'Dimensions',{'STRING15',string15_dim,'MAXSITE',maxSite_dim,'TIME',time_dim},...
            'Datatype','char',...
            'FillValue',netcdf.getConstant('NC_FILL_CHAR'),...
            'Format',ncfmt);
                
        nccreate(ncfile,'SCDT',...
            'Dimensions',{'STRING15',string15_dim,'MAXSITE',maxSite_dim,'TIME',time_dim},...
            'Datatype','char',...
            'FillValue',netcdf.getConstant('NC_FILL_CHAR'),...
            'Format',ncfmt);
        
        nccreate(ncfile,'PRPC',...
            'Dimensions',{'NPRP', 11},...
            'Datatype','double',...
            'Format',ncfmt);
        
        %% Creates attributes for the variables
        ncwriteatt(ncfile,'TIME','long_name',char('Time of Measurement UTC'));
        ncwriteatt(ncfile,'TIME','standard_name',char('time'));
        ncwriteatt(ncfile,'TIME','units',char(time_units));
        ncwriteatt(ncfile,'TIME','calendar',char('Julian'));
        ncwriteatt(ncfile,'TIME','axis',char('T'));
        ncwriteatt(ncfile,'TIME','ancillary_variables',char('TIME_SEADATANET_QC'));
                
        ncwriteatt(ncfile,'LATITUDE','long_name',char('Latitude'));
        ncwriteatt(ncfile,'LATITUDE','standard_name',char('latitude'));
        ncwriteatt(ncfile,'LATITUDE','units',char('degrees_north'));
        ncwriteatt(ncfile,'LATITUDE','axis',char('Y'));
        ncwriteatt(ncfile,'LATITUDE','ancillary_variables',char('POSITION_SEADATANET_QC'));
        
        ncwriteatt(ncfile,'LONGITUDE','long_name',char('Longitude'));
        ncwriteatt(ncfile,'LONGITUDE','standard_name',char('longitude'));
        ncwriteatt(ncfile,'LONGITUDE','units',char('degrees_east'));
        ncwriteatt(ncfile,'LONGITUDE','axis',char('X'));
        ncwriteatt(ncfile,'LONGITUDE','ancillary_variables',char('POSITION_SEADATANET_QC'));
        
        ncwriteatt(ncfile,'crs','grid_mapping_name',char('latitude_longitude'));
        ncwriteatt(ncfile,'crs','epsg_code',char('EPSG:4326'));
        ncwriteatt(ncfile,'crs','semi_major_axis',6378137.0);
        ncwriteatt(ncfile,'crs','inverse_flattening',298.257223563);
        
        ncwriteatt(ncfile,'DEPH','long_name',char('Depth of measurement'));
        ncwriteatt(ncfile,'DEPH','standard_name',char('depth'));
        ncwriteatt(ncfile,'DEPH','units',char('m'));
        ncwriteatt(ncfile,'DEPH','axis',char('Z')); 
        ncwriteatt(ncfile,'DEPH','positive',char('down'));
        ncwriteatt(ncfile,'DEPH','reference',char('sea_level'));
        ncwriteatt(ncfile,'DEPH','ancillary_variables',char('DEPTH_SEADATANET_QC'));

        ncwriteatt(ncfile,'EWCT','long_name',char('Surface Eastward Sea Water Velocity'));
        ncwriteatt(ncfile,'EWCT','standard_name',char('surface_eastward_sea_water_velocity'));
        ncwriteatt(ncfile,'EWCT','units',char('m s-1'));
        ncwriteatt(ncfile,'EWCT','scale_factor',double(1));
        ncwriteatt(ncfile,'EWCT','add_offset',double(0));
        ncwriteatt(ncfile,'EWCT','ioos_category',char('Currents'));
        ncwriteatt(ncfile,'EWCT','coordsys',char('geographic'));
        ncwriteatt(ncfile,'EWCT','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
%        ncwriteatt(ncfile,'EWCT','cell_methods',char('time: mean over hours time'));
        ncwriteatt(ncfile,'EWCT','valid_range',[double(-10.0),double(10.0)]);
%         ncwriteatt(ncfile,'EWCT','valid_min',double(-10.0));
%         ncwriteatt(ncfile,'EWCT','valid_max',double(10.0));
        ncwriteatt(ncfile,'EWCT','ancillary_variables',char('QCflag, VART_QC, CSPD_QC, DDNS_QC, GDOP_QC'));
        
        ncwriteatt(ncfile,'NSCT','long_name',char('Surface Northward Sea Water Velocity'));
        ncwriteatt(ncfile,'NSCT','standard_name',char('surface_northward_sea_water_velocity'));
        ncwriteatt(ncfile,'NSCT','units',char('m s-1'));
        ncwriteatt(ncfile,'NSCT','scale_factor',double(1));
        ncwriteatt(ncfile,'NSCT','add_offset',double(0));
        ncwriteatt(ncfile,'NSCT','ioos_category',char('Currents'));
        ncwriteatt(ncfile,'NSCT','coordsys',char('geographic'));
        ncwriteatt(ncfile,'NSCT','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
%        ncwriteatt(ncfile,'NSCT','cell_methods',char('time: mean over hours time'));
        ncwriteatt(ncfile,'NSCT','valid_range',[double(-10.0),double(10.0)]);
%         ncwriteatt(ncfile,'NSCT','valid_min',double(-10.0));
%         ncwriteatt(ncfile,'NSCT','valid_max',double(10.0));
        ncwriteatt(ncfile,'NSCT','ancillary_variables',char('QCflag, VART_QC, CSPD_QC, DDNS_QC, GDOP_QC'));
        
        ncwriteatt(ncfile,'EWCS','long_name',char('Standard Deviation of Surface Eastward Sea Water Velocity'));
%        ncwriteatt(ncfile,'EWCS','standard_name',char('surface_eastward_sea_water_velocity_standard_error'));
        ncwriteatt(ncfile,'EWCS','units',char('m s-1'));
        ncwriteatt(ncfile,'EWCS','valid_range',[double(-10.0),double(10.0)]);
        ncwriteatt(ncfile,'EWCS','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
%         ncwriteatt(ncfile,'EWCS','valid_min',double(-10.0));
%         ncwriteatt(ncfile,'EWCS','valid_max',double(10.0));
        ncwriteatt(ncfile,'EWCS','scale_factor',double(1));
        ncwriteatt(ncfile,'EWCS','add_offset',double(0));
        ncwriteatt(ncfile,'EWCS','ancillary_variables',char('QCflag, VART_QC'));
        
        ncwriteatt(ncfile,'NSCS','long_name',char('Standard Deviation of Surface Northward Sea Water Velocity'));
%        ncwriteatt(ncfile,'NSCS','standard_name',char('surface_northward_sea_water_velocity_standard_error'));
        ncwriteatt(ncfile,'NSCS','units',char('m s-1'));
        ncwriteatt(ncfile,'NSCS','valid_range',[double(-10.0),double(10.0)]);
        ncwriteatt(ncfile,'NSCS','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
%         ncwriteatt(ncfile,'NSCS','valid_min',double(-10.0));
%         ncwriteatt(ncfile,'NSCS','valid_max',double(10.0));
        ncwriteatt(ncfile,'NSCS','scale_factor',double(1));
        ncwriteatt(ncfile,'NSCS','add_offset',double(0));
        ncwriteatt(ncfile,'NSCS','ancillary_variables',char('QCflag, VART_QC'));
        
        ncwriteatt(ncfile,'CCOV','long_name',char('Covariance of Surface Sea Water Velocity'));
%         ncwriteatt(ncfile,'CCOV','standard_name',char('surface_sea_water_velocity_covariance'));
        ncwriteatt(ncfile,'CCOV','units',char('m2 s-2'));
        ncwriteatt(ncfile,'CCOV','valid_range',[double(-10.0),double(10.0)]);
        ncwriteatt(ncfile,'CCOV','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
%         ncwriteatt(ncfile,'CCOV','valid_min',double(-10.0));
%         ncwriteatt(ncfile,'CCOV','valid_max',double(10.0));
        ncwriteatt(ncfile,'CCOV','scale_factor',double(1));
        ncwriteatt(ncfile,'CCOV','add_offset',double(0));
        ncwriteatt(ncfile,'CCOV','ancillary_variables',char('QCflag'));
        
        ncwriteatt(ncfile,'GDOP','long_name',char('Geometrical Dilution of precision'));
%         ncwriteatt(ncfile,'GDOP','standard_name',char('gdop'));
        ncwriteatt(ncfile,'GDOP','units',char('1'));
        ncwriteatt(ncfile,'GDOP','valid_range',[double(-20.0),double(20.0)]);
        ncwriteatt(ncfile,'GDOP','coordinates',char('TIME DEPH LATITUDE LONGITUDE'));
%         ncwriteatt(ncfile,'GDOP','valid_min',double(-20.0));
%         ncwriteatt(ncfile,'GDOP','valid_max',double(20.0));
        ncwriteatt(ncfile,'GDOP','scale_factor',double(1));
        ncwriteatt(ncfile,'GDOP','add_offset',double(0));
        ncwriteatt(ncfile,'GDOP','comment',char(['The Geometric Dilution of Precision (GDOP) is the coefficient of the uncertainty, which relates the uncertainties in radial and velocity vectors.' ...
            ' The GDOP is a unit-less coefficient, which characterizes the effect that radar station geometry has on the measurement and position determination errors.' ...
            ' A low GDOP corresponds to an optimal geometric configuration of radar stations, and results in accurate surface current data. Essentially, GDOP is a quantitative way to relate the radial and velocity vector uncertainties.'...
            ' Setting a threshold on GDOP for total combination avoids the combination of radials with an intersection angle below a certain value.' ...
            ' GDOP is a useful metric for filtering errant velocities due to poor geometry.']));
        ncwriteatt(ncfile,'GDOP','ancillary_variables',char('QCflag, GDOP_QC'));
        
        ncwriteatt(ncfile,'TIME_SEADATANET_QC','long_name',char('Time SeaDataNet quality flag'));
        ncwriteatt(ncfile,'TIME_SEADATANET_QC','units',char('1'));
        ncwriteatt(ncfile,'TIME_SEADATANET_QC','valid_range',int16([0 9]));
        ncwriteatt(ncfile,'TIME_SEADATANET_QC','flag_values',int16([0 1 2 3 4 7 8 9]));
        ncwriteatt(ncfile,'TIME_SEADATANET_QC','flag_meanings',char('unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value'));
        ncwriteatt(ncfile,'TIME_SEADATANET_QC','comment',char('OceanSITES quality flagging for temporal coordinate.'));
        ncwriteatt(ncfile,'TIME_SEADATANET_QC','scale_factor',int16(1));
        ncwriteatt(ncfile,'TIME_SEADATANET_QC','add_offset',int16(0));
        
        ncwriteatt(ncfile,'POSITION_SEADATANET_QC','long_name',char('Position SeaDataNet quality flag'));
        ncwriteatt(ncfile,'POSITION_SEADATANET_QC','units',char('1'));
        ncwriteatt(ncfile,'POSITION_SEADATANET_QC','valid_range',int16([0 9]));
        ncwriteatt(ncfile,'POSITION_SEADATANET_QC','flag_values',int16([0 1 2 3 4 7 8 9]));
        ncwriteatt(ncfile,'POSITION_SEADATANET_QC','flag_meanings',char('unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value'));
        ncwriteatt(ncfile,'POSITION_SEADATANET_QC','comment',char('OceanSITES quality flagging for position coordinates.'));
        ncwriteatt(ncfile,'POSITION_SEADATANET_QC','scale_factor',int16(1));
        ncwriteatt(ncfile,'POSITION_SEADATANET_QC','add_offset',int16(0));
        
        ncwriteatt(ncfile,'DEPTH_SEADATANET_QC','long_name',char('Time SeaDataNet quality flag'));
        ncwriteatt(ncfile,'DEPTH_SEADATANET_QC','units',char('1'));
        ncwriteatt(ncfile,'DEPTH_SEADATANET_QC','valid_range',int16([0 9]));
        ncwriteatt(ncfile,'DEPTH_SEADATANET_QC','flag_values',int16([0 1 2 3 4 7 8 9]));
        ncwriteatt(ncfile,'DEPTH_SEADATANET_QC','flag_meanings',char('unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value'));
        ncwriteatt(ncfile,'DEPTH_SEADATANET_QC','comment',char('OceanSITES quality flagging for depth coordinate.'));
        ncwriteatt(ncfile,'DEPTH_SEADATANET_QC','scale_factor',int16(1));
        ncwriteatt(ncfile,'DEPTH_SEADATANET_QC','add_offset',int16(0));
        
        ncwriteatt(ncfile,'QCflag','long_name',char('Overall Quality Flags'));
        ncwriteatt(ncfile,'QCflag','units',char('1'));
        ncwriteatt(ncfile,'QCflag','valid_range',int16([0 9]));
        ncwriteatt(ncfile,'QCflag','flag_values',int16([0 1 2 3 4 7 8 9]));
        ncwriteatt(ncfile,'QCflag','flag_meanings',char('unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value'));
        ncwriteatt(ncfile,'QCflag','comment',char('OceanSITES quality flagging for all QC tests.'));
        ncwriteatt(ncfile,'QCflag','scale_factor',int16(1));
        ncwriteatt(ncfile,'QCflag','add_offset',int16(0));
        
        ncwriteatt(ncfile,'VART_QC','long_name',char('Variance Threshold Quality flags'));
        ncwriteatt(ncfile,'VART_QC','units',char('1'));
        ncwriteatt(ncfile,'VART_QC','valid_range',int16([0 9]));
        ncwriteatt(ncfile,'VART_QC','flag_values',int16([0 1 2 3 4 7 8 9]));
        ncwriteatt(ncfile,'VART_QC','flag_meanings',char('unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value'));
        ncwriteatt(ncfile,'VART_QC','comment',char(['OceanSITES quality flagging for variance threshold QC test. ' ...
            'Test not applicable to Direction Finding systems. The Temporal Derivative test is applied.' ...
            'Threshold set to ' num2str(Total_QC_params.TempDerThr.threshold) ' m/s. ']));
        ncwriteatt(ncfile,'VART_QC','scale_factor',int16(1));
        ncwriteatt(ncfile,'VART_QC','add_offset',int16(0));
        
        ncwriteatt(ncfile,'GDOP_QC','long_name',char('GDOP Threshold Quality flags'));
        ncwriteatt(ncfile,'GDOP_QC','units',char('1'));
        ncwriteatt(ncfile,'GDOP_QC','valid_range',int16([0 9]));
        ncwriteatt(ncfile,'GDOP_QC','flag_values',int16([0 1 2 3 4 7 8 9]));
        ncwriteatt(ncfile,'GDOP_QC','flag_meanings',char('unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value'));
        ncwriteatt(ncfile,'GDOP_QC','comment',char(['OceanSITES quality flagging for GDOP threshold QC test. ' ...
            'Threshold set to ' num2str(Total_QC_params.GDOPThr) '.']));
        ncwriteatt(ncfile,'GDOP_QC','scale_factor',int16(1));
        ncwriteatt(ncfile,'GDOP_QC','add_offset',int16(0));
        
        ncwriteatt(ncfile,'DDNS_QC','long_name',char('Data density Threshold Quality flags'));
        ncwriteatt(ncfile,'DDNS_QC','units',char('1'));
        ncwriteatt(ncfile,'DDNS_QC','valid_range',int16([0 9]));
        ncwriteatt(ncfile,'DDNS_QC','flag_values',int16([0 1 2 3 4 7 8 9]));
        ncwriteatt(ncfile,'DDNS_QC','flag_meanings',char('unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value'));
        ncwriteatt(ncfile,'DDNS_QC','comment',char(['OceanSITES quality flagging for Data density threshold QC test. ' ...
            'Threshold set to ' num2str(Total_QC_params.DataDensityThr) ' radials.']));
        ncwriteatt(ncfile,'DDNS_QC','scale_factor',int16(1));
        ncwriteatt(ncfile,'DDNS_QC','add_offset',int16(0));
         
        ncwriteatt(ncfile,'CSPD_QC','long_name',char('Velocity threshold Quality flags'));
        ncwriteatt(ncfile,'CSPD_QC','units',char('1'));
        ncwriteatt(ncfile,'CSPD_QC','valid_range',int16([0 9]));
        ncwriteatt(ncfile,'CSPD_QC','flag_values',int16([0 1 2 3 4 7 8 9]));
        ncwriteatt(ncfile,'CSPD_QC','flag_meanings',char('unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value'));
        ncwriteatt(ncfile,'CSPD_QC','comment',char(['OceanSITES quality flagging for Velocity threshold QC test. ' ...
            'Threshold set to ' num2str(Total_QC_params.VelThr) ' m/s.']));
        ncwriteatt(ncfile,'CSPD_QC','scale_factor',int16(1));
        ncwriteatt(ncfile,'CSPD_QC','add_offset',int16(0));
        
        ncwriteatt(ncfile,'NARX','long_name',char('Number of Receive Antennas'));
        ncwriteatt(ncfile,'NARX','units',char('1'));
        ncwriteatt(ncfile,'NARX','valid_range',int16([0 maxSite_dim]));
%         ncwriteatt(ncfile,'NARX','coordinates',char('TIME'));
        ncwriteatt(ncfile,'NARX','scale_factor',int16(1));
        ncwriteatt(ncfile,'NARX','add_offset',int16(0));
        
        ncwriteatt(ncfile,'NATX','long_name',char('Number of Transmit Antennas'));
        ncwriteatt(ncfile,'NATX','units',char('1'));
        ncwriteatt(ncfile,'NATX','valid_range',int16([0 maxSite_dim]));
%         ncwriteatt(ncfile,'NATX','coordinates',char('TIME'));
        ncwriteatt(ncfile,'NATX','scale_factor',int16(1));
        ncwriteatt(ncfile,'NATX','add_offset',int16(0));
        
        ncwriteatt(ncfile,'SLTR','long_name',char('Receive Antennas Latitudes'));
        ncwriteatt(ncfile,'SLTR','standard_name',char('latitude'));
        ncwriteatt(ncfile,'SLTR','units','degrees_north');
        ncwriteatt(ncfile,'SLTR','valid_range',single([-90 90]));
        ncwriteatt(ncfile,'SLTR','coordinates',char('TIME MAXSITE'));
        ncwriteatt(ncfile,'SLTR','scale_factor',single(1));
        ncwriteatt(ncfile,'SLTR','add_offset',single(0));
        
        ncwriteatt(ncfile,'SLNR','long_name',char('Receive Antennas Longitudes'));
        ncwriteatt(ncfile,'SLNR','standard_name',char('longitude'));
        ncwriteatt(ncfile,'SLNR','units','degrees_east');
        ncwriteatt(ncfile,'SLNR','valid_range',single([-180 180]));
        ncwriteatt(ncfile,'SLNR','coordinates',char('TIME MAXSITE'));
        ncwriteatt(ncfile,'SLNR','scale_factor',single(1));
        ncwriteatt(ncfile,'SLNR','add_offset',single(0));
        
        ncwriteatt(ncfile,'SLTT','long_name',char('Transmit Antennas Latitudes'));
        ncwriteatt(ncfile,'SLTT','standard_name',char('latitude'));
        ncwriteatt(ncfile,'SLTT','units','degrees_north');
        ncwriteatt(ncfile,'SLTT','valid_range',single([-90 90]));
        ncwriteatt(ncfile,'SLTT','coordinates',char('TIME MAXSITE'));
        ncwriteatt(ncfile,'SLTT','scale_factor',single(1));
        ncwriteatt(ncfile,'SLTT','add_offset',single(0));
        
        ncwriteatt(ncfile,'SLNT','long_name',char('Transmit Antennas Longitudes'));
        ncwriteatt(ncfile,'SLNT','standard_name',char('longitude'));
        ncwriteatt(ncfile,'SLNT','units','degrees_east');
        ncwriteatt(ncfile,'SLNT','valid_range',single([-180 180]));
        ncwriteatt(ncfile,'SLNT','coordinates',char('TIME MAXSITE'));
        ncwriteatt(ncfile,'SLNT','scale_factor',single(1));
        ncwriteatt(ncfile,'SLNT','add_offset',single(0));       
        
        ncwriteatt(ncfile,'SCDR','long_name',char('Receive antenna Codes'));

        ncwriteatt(ncfile,'SCDT','long_name',char('Transmit antenna Codes'));
    
        ncwriteatt(ncfile,'PRPC','long_name',char('Near Real Time data processing parameters'));
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
%         ncwrite(ncfile,'TIME',int32((mat_tot.TimeStamp-timeref)*86400));
        ncwrite(ncfile,'TIME',single(mat_tot.TimeStamp-timeref));
        ncwrite(ncfile,'LATITUDE',latGrid);
        ncwrite(ncfile,'LONGITUDE',lonGrid);
        ncwrite(ncfile,'DEPH',depth);
        ncwrite(ncfile,'EWCT',mat_tot.U_grid);
        ncwrite(ncfile,'NSCT',mat_tot.V_grid);
        ncwrite(ncfile,'EWCS',mat_tot.U_std);
        ncwrite(ncfile,'NSCS',mat_tot.V_std);
        ncwrite(ncfile,'CCOV',mat_tot.covariance);
        ncwrite(ncfile,'GDOP',mat_tot.GDOP);
        ncwrite(ncfile,'NARX',length(sitesLat));
        ncwrite(ncfile,'NATX',length(sitesLat));       
        ncwrite(ncfile,'SLTR',sitesLat');
        ncwrite(ncfile,'SLNR',sitesLon');
        ncwrite(ncfile,'SLTT',sitesLat');
        ncwrite(ncfile,'SLNT',sitesLon');       
        ncwrite(ncfile,'SCDR',sitesCodes');
        ncwrite(ncfile,'SCDT',sitesCodes');
        ncwrite(ncfile,'PRPC',procParamsVec);
        ncwrite(ncfile,'TIME_SEADATANET_QC',sdnTime_QCflag);
        ncwrite(ncfile,'POSITION_SEADATANET_QC',sdnPosition_QCflag);
        ncwrite(ncfile,'DEPTH_SEADATANET_QC',sdnDepth_QCflag);
        ncwrite(ncfile,'QCflag',overall_QCflag);
        ncwrite(ncfile,'VART_QC',temporalDerivativeThreshold_QCflag);
        ncwrite(ncfile,'GDOP_QC',GDOP_QCflag);
        ncwrite(ncfile,'DDNS_QC',dataDensity_QCflag);
        ncwrite(ncfile,'CSPD_QC',velocityThreshold_QCflag);
        
        %% Define global attributes
        
        % MANDATORY ATTRIBUTES
        % Discovery and Identification
        ncwriteatt(ncfile,'/','site_code',char('HFR_TirLig'));
        ncwriteatt(ncfile,'/','platform_code',char('HFR_TirLig_Total'));
        ncwriteatt(ncfile,'/','data_mode',char('R'));
        ncwriteatt(ncfile,'/','DoA_estimation_method',char('Direction Finding'));
        ncwriteatt(ncfile,'/','calibration_type',char('APM'));
        ncwriteatt(ncfile,'/','last_calibration_date',char(lastPatternStr));
        ncwriteatt(ncfile,'/','calibration_link',char('carlo.mantovani@cnr.it'));
        ncwriteatt(ncfile,'/','title',char('Near Real Time Surface Ocean Velocity by HFR_TirLig'));
        ncwriteatt(ncfile,'/','summary',char('The data set consists of maps of total velocity of the surface current in the North-Western Tyrrhenian Sea and Ligurian Sea averaged over a time interval of 1 hour around the cardinal hour. Surface ocean velocities estimated by HF Radar are representative of the upper 0.3-2.5 meters of the ocean.'));
        ncwriteatt(ncfile,'/','source',char('coastal structure'));
        ncwriteatt(ncfile,'/','source_platform_category_code',char('17'));
        ncwriteatt(ncfile,'/','institution',char('National Research Council - Institute of Marine Science, S.S. Lerici'));
        ncwriteatt(ncfile,'/','institution_edmo_code',char('134'));
        ncwriteatt(ncfile,'/','data_assembly_center',char('European HFR Node'));
        ncwriteatt(ncfile,'/','id',char(dataID));
                
        % Geo-spatial-temporal
        ncwriteatt(ncfile,'/','data_type', char('HF radar total data'));
        ncwriteatt(ncfile,'/','feature_type',char('surface'));
        ncwriteatt(ncfile,'/','geospatial_lat_min',char(num2str(min(map_lat_lim))));
        ncwriteatt(ncfile,'/','geospatial_lat_max',char(num2str(max(map_lat_lim))));      
        ncwriteatt(ncfile,'/','geospatial_lon_min',char(num2str(min(map_lon_lim))));
        ncwriteatt(ncfile,'/','geospatial_lon_max',char(num2str(max(map_lon_lim))));
        ncwriteatt(ncfile,'/','geospatial_vertical_min', char('0'));
        ncwriteatt(ncfile,'/','geospatial_vertical_max', char('0.48'));
        ncwriteatt(ncfile, '/','time_coverage_start',char(timeCoverageStart));
        ncwriteatt(ncfile, '/','time_coverage_end',char(timeCoverageEnd));
        % Conventions used
        ncwriteatt(ncfile,'/','format_version',char('v2.1'));
        ncwriteatt(ncfile,'/','Conventions',char('CF-1.6, OceanSITES-Manual-1.2, Copernicus-InSituTAC-SRD-1.4, CopernicusInSituTAC-ParametersList-3.1.0, Unidata, ACDD, INSPIRE'));
        % Publication information
        ncwriteatt(ncfile,'/','update_interval',char('void'));
        ncwriteatt(ncfile,'/','citation',char(citation_str));  
        ncwriteatt(ncfile,'/','distribution_statement',char(distribution_str));
        ncwriteatt(ncfile,'/','publisher_name',char('Lorenzo Corgnati'));
        ncwriteatt(ncfile,'/','publisher_url',char('http://radarhf.ismar.cnr.it'));
        ncwriteatt(ncfile,'/','publisher_email',char('lorenzo.corgnati@sp.ismar.cnr.it'));
        ncwriteatt(ncfile,'/','license',char('HF radar sea surface current velocity dataset by CNR-ISMAR is licensed under a Creative Commons Attribution 4.0 International License. You should have received a copy of the license along with this work. If not, see http://creativecommons.org/licenses/by/4.0/.'));
        ncwriteatt(ncfile,'/','acknowledgment',char('ISMAR HF Radar Network has been established within RITMARE and Jerico-Next projects. The network has been designed, implemented and managed through the efforts of ISMAR S.S. Lerici.'));
        % Provenance
        ncwriteatt(ncfile,'/','date_created',char(dateCreated));      
        ncwriteatt(ncfile,'/','history',char([time_coll ' data collected. ' dateCreated ' netCDF file created and sent to European HFR Node']));
        ncwriteatt(ncfile,'/','date_modified',char(dateCreated));
        ncwriteatt(ncfile,'/','date_update',char(dateCreated));
        ncwriteatt(ncfile,'/','processing_level',char('3B'));
        ncwriteatt(ncfile,'/','contributor_name',char('Vega Forneris, Cristina Tronconi'));
        ncwriteatt(ncfile,'/','contributor_role',char('THREDDS expert, metadata expert'));
        ncwriteatt(ncfile,'/','contributor_email',char('vega.forneris@artov.isac.cnr.it, cristina.tronconi@artov.isac.cnr.it'));

        % RECOMMENDED ATTRIBUTES
        % Discovery and Identification        
        ncwriteatt(ncfile,'/','project',char('RITMARE and Jerico-Next'));
        ncwriteatt(ncfile,'/','naming_authority',char('it.cnr.ismar'));
        ncwriteatt(ncfile,'/','keywords',char('OCEAN CURRENTS, SURFACE WATER, RADAR, SCR-HF'));
        ncwriteatt(ncfile,'/','keywords_vocabulary',char('GCMD Science Keywords'));
        ncwriteatt(ncfile,'/','comment',char('Total velocities are derived using least square fit that maps radial velocities measured from individual sites onto a cartesian grid. The final product is a map of the horizontal components of the ocean currents on a regular grid in the area of overlap of two or more radar stations.'));
        ncwriteatt(ncfile,'/','data_language',char('eng'));
        ncwriteatt(ncfile,'/','data_character_set',char('utf8'));
        ncwriteatt(ncfile,'/','metadata_language',char('eng'));
        ncwriteatt(ncfile,'/','metadata_character_set',char('utf8'));
        ncwriteatt(ncfile,'/','topic_category',char('oceans'));
        ncwriteatt(ncfile,'/','network',char('ISMAR_HFR_TirLig'));
        % Geo-spatial-temporal
        ncwriteatt(ncfile,'/','area',char('Mediterranean Sea'));
        ncwriteatt(ncfile,'/','geospatial_lat_units',char('degrees_north'));
        ncwriteatt(ncfile,'/','geospatial_lon_units',char('degrees_east'));
        ncwriteatt(ncfile,'/','geospatial_lat_resolution',latRes(1));
        ncwriteatt(ncfile,'/','geospatial_lon_resolution',lonRes(1));
        ncwriteatt(ncfile,'/','geospatial_vertical_resolution', char('0.48'));
        ncwriteatt(ncfile,'/','geospatial_vertical_units', char('m'));
        ncwriteatt(ncfile,'/','geospatial_vertical_positive', char('down')); 
        ncwriteatt(ncfile, '/','time_coverage_duration',char('PT1H'));
        ncwriteatt(ncfile, '/','time_coverage_resolution',char('PT1H'));
        ncwriteatt(ncfile,'/','reference_system',char('EPSG:4806'));
        ncwriteatt(ncfile,'/','grid_resolution',char('1.5 Km'));
        ncwriteatt(ncfile,'/','cdm_data_type',char('Grid'));
        % Conventions used
        ncwriteatt(ncfile,'/','netcdf_version',char(netcdf.inqLibVers));
        ncwriteatt(ncfile,'/','netcdf_format',char(ncfmt));
        
        % OTHER ATTRIBUTES    
        ncwriteatt(ncfile,'/','metadata_contact',char('lorenzo.corgnati@sp.ismar.cnr.it'));
        ncwriteatt(ncfile,'/','metadata_date_stamp',char(dateCreated)); 
        ncwriteatt(ncfile,'/','abstract',char('The data set consists of maps of total velocity of the surface current in the North-Western Tyrrhenian Sea and Ligurian Sea averaged over a time interval of 1 hour around the cardinal hour. Surface ocean velocities estimated by HF Radar are representative of the upper 0.3-2.5 meters of the ocean. Total velocities are derived using least square fit that maps radial velocities measured from individual sites onto a cartesian grid. The final product is a map of the horizontal components of the ocean currents on a regular grid in the area of overlap of two or more radar stations.'));
        ncwriteatt(ncfile,'/','standard_name_vocabulary',char('NetCDF Climate and Forecast (CF) Metadata Convention Standard Name Table Version 1.6'));     
        ncwriteatt(ncfile,'/','creator_name',char('Lorenzo Corgnati'));
        ncwriteatt(ncfile,'/','creator_url',char('http://radarhf.ismar.cnr.it'));
        ncwriteatt(ncfile,'/','creator_email',char('lorenzo.corgnati@sp.ismar.cnr.it'));        
        ncwriteatt(ncfile,'/','sensor',char('CODAR SeaSonde'));
        ncwriteatt(ncfile,'/','grid_mapping',char('Transverse Mercator'));            
        ncwriteatt(ncfile,'/','institution_reference',char('http://www.ismar.cnr.it'));
        ncwriteatt(ncfile,'/','operational_manager',char('Carlo Mantovani'));
        ncwriteatt(ncfile,'/','operational_manager_email',char('carlo.mantovani@cnr.it'));       
        ncwriteatt(ncfile,'/','regional_description',char('North-Western Tyrrhenian Sea and Ligurian Sea, Italy'));
        ncwriteatt(ncfile,'/','date_issued',char(dateCreated));
        ncwriteatt(ncfile,'/','software_name',char('HFR_Combiner_TirLig'));
        ncwriteatt(ncfile,'/','software_version',char('v3.1'));       
        ncwriteatt(ncfile,'/','references',char('HFR_Progs Matlab Documentation - Copyright (C) 2006-7 David M. Kaplan; Otero,M. (2008).NETCDF DESCRIPTION FOR NEAR REAL-TIME SURFACE CURRENTS PRODUCED BY THE HF-RADAR NETWORK. https://cordc.ucsd.edu/projects/mapping/documents/HFRNet_RTV-NetCDF.pdf'));
        
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        T2C_err = 1;
    end
end

% Copy the netCDF file to the RadarDisk
if (T2C_err == 0)
    try
        [cp_status] = copyfile(ncfile,servRD_ncfile,'f');
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        T2C_err = 1;
    end
end
return

