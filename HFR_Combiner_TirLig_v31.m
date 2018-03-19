%% HFR_Combiner_TirLig_v31.m
% This application automatically reads Radial data from HF radar
% installation in Tyrrhenian and Ligurian Sea (MONT, TINO), generates totals
% and stores radial and total files.
% Grid: 1.5 km - Combination Radius: 3 km.

% The 1.1 release is able to recover old radial data used for uncomplete
% total maps generation in order to remake totals using new radial data
% sent with a delay by the communication system.

% The 1.1.x release works on shared folder radardisk and radarserver and
% plots total maps in an optimized way (reduced area and more visible
% arrows) and saves incremental copies of RadStruct.mat file.

% The 1.1.5 release reads the input parameters from an external text file,
% namely the folder paths, the latitude and longitudes values for maps
% representation and for total maps evaluation, the velocity threshold for
% totals cleaning, the spatial threshold for totals generataion and the
% time period limits for totals evaluation.

% The 1.2 release starts from the architecture of 1.1.5 version and
% implements cleaning and corrections on the radial data in order to
% operate a first level Quality control.

% The 1.3 release stores hourly total data in netCDF format.

% The 1.3.1 release stores hourly total data in netCDF format on
% the THREDDS server, into the RITMARE dedicated folders.

% The 1.3.2 release runs a function at its startup in order to delete from
% the RadarDisk RadialRealTime folder the radial files older than a week
% (Archivalist rsyncs radial files of the last week). This task speeds up
% the directory lisiting and RadStruct operations.

% The 1.3.4 release generates total velocity files with current
% measurement unit in m/s according to the EmodNet requirements.

% This release stores the tuv and the netCDF total files in a folder structure
% with the scheme yyyy/yyyy_mm/yyyy_mm_dd/ locally and on the RadarDisk,
% and with the structure Last/yyyy_mm_dd/ on the THREDDS server.

% The 1.3.5 release makes local copies of the radial files locally and
% processes them locally. Then the resulting total products are moved in
% the destination folders. The local copies are deleted after processing.

% The 2.0 release creates netCDF 4 files containing total velocities, U and V
% errors, UV covariance and GDOP. The file format is compliant with CF-1.6
% Unidata Data Discoverty, INSPIRE and ACDD conventions.

% The 2.1 release stores hourly radial and total data in netCDF format on
% the THREDDS server, into the RITMARE dedicated folders.

% The 2.2 release stores hourly radial and total data in netCDF format in a
% folder structure compliant to the ISMAR THREDDS requirements.

% The 3.0 release implements the data and metadata structure compliant to
% CMEMS needs and performs QC tests.

% The 3.1 release implements the data and metadata structure compliant to
% CMEMS-INSTAC requirements and performs QC tests.

% Author: Lorenzo Corgnati
% Date: March 15, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

warning('off', 'all');

clear all
close all
clc

warning('off','all');

setup_nctoolbox;

display(['[' datestr(now) '] - - ' 'HFR_Combiner_TirLig_v31 started.']);

tic;
%% Parameters
% Reads the text file containing the parameters
% param_file = textread('/home/radarcombine/HFR_Combiner_TirLig/Src/parameters_TirLig_Meas_v31.txt',  '%s', 'whitespace', '\n');
param_file = textread('/Users/reverendo/Documents/CNR/RADAR/Script/hfrprogs/Src/parameters_TirLig_test.txt',  '%s', 'whitespace', '\n');

% RadarDisk ftp account
RD_ftp_host = param_file(2);
RD_ftp_host = RD_ftp_host{1};
RD_ftp_username = param_file(4);
RD_ftp_username = RD_ftp_username{1};
RD_ftp_password = param_file(6);
RD_ftp_password = RD_ftp_password{1};

% Folder paths
work_folder = param_file(8);
work_folder = work_folder{1};
src_folder = param_file(10);
src_folder = src_folder{1};
incr_folder = param_file(12);
incr_folder = incr_folder{1};
radProc_folder = param_file(14);
radProc_folder = radProc_folder{1};
site1_folder = param_file(16);
site1_folder = site1_folder{1};
site2_folder = param_file(18);
site2_folder = site2_folder{1};
site3_folder = param_file(20);
site3_folder = site3_folder{1};
site4_folder = param_file(22);
site4_folder = site4_folder{1};
site5_folder = param_file(24);
site5_folder = site5_folder{1};
site6_folder = param_file(26);
site6_folder = site6_folder{1};
total_maps_folder = param_file(28);
total_maps_folder = total_maps_folder{1};
total_maps_server_folder = param_file(30);
total_maps_server_folder = total_maps_server_folder{1};
tuv_local_folder = param_file(32);
tuv_local_folder = tuv_local_folder{1};
tuv_serv_folder = param_file(34);
tuv_serv_folder = tuv_serv_folder{1};
netcdf_thredds_folder = param_file(36);
netcdf_thredds_folder = netcdf_thredds_folder{1};

% Longitude and latitude limit values for maps representation.
lon_lim = textscan(char(param_file(38)), '%f');
lon_lim = lon_lim{1}';
lat_lim = textscan(char(param_file(40)), '%f');
lat_lim = lat_lim{1}';

% Longitude and latitude limit values of the grid for totals generation.
lon_grid = textscan(char(param_file(42)), '%f');
lon_grid = lon_grid{1}';
lat_grid = textscan(char(param_file(44)), '%f');
lat_grid = lat_grid{1}';
lon_points = textscan(char(param_file(46)), '%d');
lon_points = lon_points{1};
lat_points = textscan(char(param_file(48)), '%d');
lat_points = lat_points{1};

% Radial QC tests thresholds
% Median Filter Threshold
medFilt = textscan(char(param_file(50)), '%f');
medFilt = medFilt{1}';

% Variance Threshold
var_thr_R = textscan(char(param_file(52)), '%f');
var_thr_R = var_thr_R{1}';

% Temporal Derivative Threshold
temp_der_thr_str_R = textscan(char(param_file(54)), '%f');
temp_der_thr_R.threshold = temp_der_thr_str_R{1}';

% Velocity Threshold
maxspd_R = textscan(char(param_file(56)), '%f');
maxspd_R = maxspd_R{1}';

% Average Radial Bearing Range
avgRadBear_str_R = textscan(char(param_file(58)), '%f');
avgRadBear_R.site(1).code = 'MONT';
avgRadBear_R.site(1).range = avgRadBear_str_R{1}';

avgRadBear_str_R = textscan(char(param_file(60)), '%f');
avgRadBear_R.site(2).code = 'TINO';
avgRadBear_R.site(2).range = avgRadBear_str_R{1}';

% Radial Count Threshold
radCount_R = textscan(char(param_file(62)), '%f');
radCount_R = radCount_R{1}';

% Total QC tests thresholds
% Velocity threshold
maxspd_T = textscan(char(param_file(64)), '%f');
maxspd_T = maxspd_T{1};

% GDOP threshold
GDOPthresh = textscan(char(param_file(66)), '%f');
GDOPthresh = GDOPthresh{1};

% Variance Threshold
var_thr_T = textscan(char(param_file(68)), '%f');
var_thr_T = var_thr_T{1}';

% Temporal Derivative Threshold
temp_der_thr_str_T = textscan(char(param_file(70)), '%f');
temp_der_thr_T.threshold = temp_der_thr_str_T{1}';

% Data Density Threshold
dataDens_T = textscan(char(param_file(72)), '%f');
dataDens_T = dataDens_T{1}';

% Spatial threshold for total generation [km].
spatthresh = textscan(char(param_file(74)), '%f');
spatthresh = spatthresh{1};

% Totals evaluation start and stop times
start_str = param_file(76);
start_str = start_str{1};
stop_str = param_file(78);
stop_str = stop_str{1};

% Default number of range cells for radial data
range_cells_number = textscan(char(param_file(80)), '%f');
range_cells_number = range_cells_number{1};

display(['[' datestr(now) '] - - ' 'Parameters successfully read.']);

%%

%% Build structures containing QC tests parameters
% Radial QC tests parameters
Radial_QC_params.MedFilt = medFilt;
Radial_QC_params.VarThr = var_thr_R;
Radial_QC_params.TempDerThr = temp_der_thr_R;
Radial_QC_params.RadCnt = radCount_R;
Radial_QC_params.VelThr = maxspd_R;
Radial_QC_params.AvgRadBear = avgRadBear_R;

% Total QC tests parameters
Total_QC_params.GDOPThr = GDOPthresh;
Total_QC_params.VarThr = var_thr_T;
Total_QC_params.TempDerThr = temp_der_thr_T;
Total_QC_params.VelThr = maxspd_T;
Total_QC_params.DataDensityThr = dataDens_T;

%%

%% Initialization of RadStruct variable
% If the RadStruct.mat file has been saved, it's loaded, otherwise
% it's initialized. If the RadStruct.mat file is corrupted, the last
% non-corrupted incremental copy of the file is loaded.
src_dir = dir(src_folder);
RadStruct_flag = 0;
for rs=3:length(src_dir)
    if (strcmp(src_dir(rs,1).name,'RadStruct.mat'))
        RadStruct_flag = 1;
    end
end
if (RadStruct_flag ==1)
    try
        load(strcat(src_folder, 'RadStruct.mat'));
        display(['[' datestr(now) '] - - ' 'RadStruct successfully loaded.']);
    catch err
        % Search for the last non-corrupted incremental copy of RadStruct.mat
        incr_dir = dir(incr_folder);
        dot_flag = 1;
        while (dot_flag == 1)
            if (incr_dir(1).name(1) == '.')
                incr_dir = incr_dir(2:size(incr_dir,1));
            else
                dot_flag = 0;
            end
        end
        incr_err = 1;
        icr = length(incr_dir);
        while (incr_err == 1 && icr > 0)
            try
                load(strcat(incr_folder, incr_dir(icr).name));
                display(['[' datestr(now) '] - - ' 'Incremental version of RadStruct successfully loaded.']);
                incr_err = 0;
            catch err
                icr = icr - 1;
                display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            end
        end
    end
else
    RadStruct(1).TS = '';
    RadStruct(1).site1 = '';
    RadStruct(1).site2 = '';
    RadStruct(1).site3 = '';
    RadStruct(1).site4 = '';
    RadStruct(1).processedFlag = 0;
    RadStruct(1).missingID = [1 1 1 1];
    RadStruct(1).corruptedFlag = 0;
    
    display(['[' datestr(now) '] - - ' 'RadStruct successfully initialized.']);
end

%%


% %% Coastline generation for area of interest
% % If the AoI.fig file has been saved, it's loaded, otherwise the
% % map is generated.
% AoI_flag = 0;
% for aoi=3:length(src_dir)
%     if (strcmp(src_dir(aoi,1).name,'AoI.fig'))
%         AoI_flag = 1;
%     end
% end
% if (AoI_flag == 0)
%     m_proj('transverse mercator','longitudes',lon_lim,'latitudes',lat_lim);
%     m_gshhs_f('patch',[.7 .7 .7],'edgecolor','none');
%     m_grid('box', 'fancy', 'tickdir', 'in', 'xlabeldir','end','fontsize',10);
%     
%     [X,Y]=m_ll2xy(9.6533333,44.1458333);
%     line(X,Y,'marker','square','markersize',4,'color','r');
%     text(X,Y,' MONT','vertical','top');
%     
%     [X,Y]=m_ll2xy(9.8492167,44.0263667);
%     line(X,Y,'marker','square','markersize',4,'color','r');
%     text(X,Y,' TINO','vertical','top');
%     
%     %     [X,Y]=m_ll2xy(16.1847500,41.8890333);
%     %     line(X,Y,'marker','square','markersize',4,'color','r');
%     %     text(X,Y,' VIES','vertical','top');
%     %
%     %     [X,Y]=m_ll2xy(16.1922000,41.7825667);
%     %     line(X,Y,'marker','square','markersize',4,'color','r');
%     %     text(X,Y,' PUGN','vertical','top');
%     %
%     %     [X,Y]=m_ll2xy(16.1162167,41.7310667);
%     %     line(X,Y,'marker','square','markersize',4,'color','r');
%     %     text(X,Y,' MATT','vertical','top');
%     %
%     %     [X,Y]=m_ll2xy(15.9253,41.6206667);
%     %     line(X,Y,'marker','square','markersize',4,'color','r');
%     %     text(X,Y,' MANF','vertical','top');
%     
%     h = gcf;
%     saveas(h, [src_folder, 'AoI.fig']);
%     close;
%     
%     display(['[' datestr(now) '] - - ' 'Coastline successfully generated.']);
% end
% %%

%% Grid generation for radial combination
lon = linspace(lon_grid(1),lon_grid(2),lon_points);
lat = linspace(lat_grid(1),lat_grid(2),lat_points);
length_lon=length(lon);
length_lat=length(lat);

for i=1:length_lon
    lonG(1+(i-1)*length_lat:(i-1)*length_lat+length_lat) = lon(i)*ones(1,length_lat);
    latG(1+(i-1)*length_lat:(i-1)*length_lat+length_lat) = lat;
end

Grid(:,1) = lonG';
Grid(:,2) = latG';

display(['[' datestr(now) '] - - ' 'Grid for radial combination successfully generated.']);
%%

%% Mask area generation.
fname = makeCoast(lon_lim,lat_lim,'transverse mercator',strcat(src_folder, 'MaskMap.mat'),5);
load(strcat(src_folder, 'MaskMap.mat'));

display(['[' datestr(now) '] - - ' 'Mask area successfully generated.']);
%%

setupTime = toc;
%% Processes data.
kk = 5;
while (kk > 0)
    % RadarDisk RadialRealTime folder cleaning.
    APT_err = 0;
    try
        [cl_status] = cleanRRTfolder(work_folder);
    catch err
        APT_err = 1;
    end
    
    tic;
    radialfiles = dir(work_folder);
    fileListingTime = toc;
    size_rad = size(radialfiles,1);
    if (size_rad > 0)
        dot_flag = 1;
        while (dot_flag == 1)
            if (radialfiles(1).name(1) == '.')
                radialfiles = radialfiles(2:size(radialfiles,1));
            else
                dot_flag = 0;
            end
        end
    end
    
    % Checks the correct loading of the work folder. If it has not been
    % loaded correctly, the application pauses before to retry.
    while (size_rad < 1)
        pause(60);
        radialfiles = dir(work_folder);
        size_rad = size(radialfiles,1);
        if (size_rad > 0)
            dot_flag = 1;
            while (dot_flag == 1)
                if (radialfiles(1).name(1) == '.')
                    radialfiles = radialfiles(2:size(radialfiles,1));
                else
                    dot_flag = 0;
                end
            end
        end
    end
    
    display(['[' datestr(now) '] - - ' 'Radial files successfully listed.']);
    
    % Builds the radials information structure RadStruct as:
    % - TS = TimeStamp (string);
    % - site1 = VIAR radial filename (string);
    % - site2 = TINO radial filename (string);
    % - site3 = PCOR radial filename (string);
    % - site4 = PFIN radial filename (string);
    % - processedFlag = processing status flag (int, 0 if not yet
    %   processed, 1 if already processed;
    % - missingID = IDs of the eventually not evaluated radials
    %   (int array of 4 cells - the value of each cell is 0 if the
    %   corresponding radial is not missing, 1 if the corresponding
    %   radial is missing);
    % - corruptedFlag = health status of the total (int, 0 if healthy,
    %   1 if any error occurred).
    tic;
    num_rad = size(radialfiles,1); % Redundant (it's the same as size_rad)
    if (RadStruct_flag == 0)
        rs_idx = 1;
    else
        rs_idx = length(RadStruct) + 1;
    end
    while (num_rad > 0)
        curr_TS = radialfiles(1,1).name(11:25);
        curr_IDstr = radialfiles(1,1).name(6);
        switch curr_IDstr
            case 'V'
                curr_ID = 1; % VIAR
            case 'T'
                curr_ID = 2; % TINO
            case 'M'
                curr_ID = 3; % PCOR
            otherwise
                curr_ID = 4; % PFIN
        end
        flagMatch = 0;
        for rss=1:length(RadStruct) % Sweeps RadStruct
            if (strcmp(curr_TS, RadStruct(rss).TS)) % Finds an element in RadSruct with the same TS
                flagMatch = 1;
                if (sum(RadStruct(rss).missingID) == 2) % Checks if it's complete
                    break;
                else
                    if (RadStruct(rss).missingID(curr_ID) == 1) % Checks if the current ID is missing
                        switch curr_ID % Inserts the current radial
                            case 1
                                RadStruct(rss).site1 = radialfiles(1,1).name;
                                RadStruct(rss).missingID(1) = 0; % Sets the current ID as not missing
                            case 2
                                RadStruct(rss).site2 = radialfiles(1,1).name;
                                RadStruct(rss).missingID(2) = 0; % Sets the current ID as not missing
                            case 3
                                RadStruct(rss).site3 = radialfiles(1,1).name;
                                RadStruct(rss).missingID(3) = 0; % Sets the current ID as not missing
                            otherwise
                                RadStruct(rss).site4 = radialfiles(1,1).name;
                                RadStruct(rss).missingID(4) = 0; % Sets the current ID as not missing
                        end
                        RadStruct(rss).processedFlag = 0; % Sets the current group of radials for being reprocessed
                    else
                        break;
                    end
                end
            end
            if (flagMatch == 1)
                break;
            end
        end
        if (flagMatch == 0) % No TS match found
            RadStruct(rs_idx).missingID = [1 1 1 1];
            switch curr_ID % Inserts the current radial
                case 1
                    RadStruct(rs_idx).site1 = radialfiles(1,1).name;
                case 2
                    RadStruct(rs_idx).site2 = radialfiles(1,1).name;
                case 3
                    RadStruct(rs_idx).site3 = radialfiles(1,1).name;
                otherwise
                    RadStruct(rs_idx).site4 = radialfiles(1,1).name;
            end
            RadStruct(rs_idx).missingID(curr_ID) = 0;
            RadStruct(rs_idx).TS = curr_TS;
            RadStruct(rs_idx).processedFlag = 0;
            RadStruct(rs_idx).corruptedFlag = 0;
            rs_idx = rs_idx + 1;
        end
        % Saves RadStruct and removes the current radial from the radialfiles list
        save(strcat(src_folder, 'RadStruct.mat'), 'RadStruct');
        % Removes the processed radial from radialfiles list
        if (size(radialfiles,1) > 1)
            radialfiles = radialfiles(2:size(radialfiles,1));
            num_rad = size(radialfiles,1);
        else
            num_rad = 0;
        end
    end
    
    display(['[' datestr(now) '] - - ' 'RadStruct successfully filled in.']);
    
    % Sorts RadStruct by date
    RadStruct = RSsort_TS(RadStruct);
    display(['[' datestr(now) '] - - ' 'RadStruct successfully sorted.']);
    
    % Periodically saves incremental copies of RadStruct
    now = clock;
    if (mod(now(4),4) == 0) && (now(5) < 16)
        incr_str = strcat((int2str(now(1))), '_', (int2str(now(2))), '_', (int2str(now(3))), '_', (int2str(now(4))), '_', (int2str(now(5))));
        save(strcat(incr_folder, 'RadStruct_', incr_str, '.mat'), 'RadStruct');
        display(['[' datestr(now) '] - - ' 'Incremental copy of RadStruct successfully saved.']);
    end
    RadStructBuildingTime = toc;
    
    % Processes radials
    % Checks if start and stop times are present
    if (start_str == '0')
        start_time = 1;
    else
        start_time = findTS(RadStruct, start_str);
    end
    if (stop_str == '0')
        stop_time = length(RadStruct);
    else
        stop_time = findTS(RadStruct, stop_str);
    end
    display(['[' datestr(now) '] - - ' 'Processing start and stop times successfully checked.']);
    for rss=start_time:stop_time % Sweeps RadStruct
        if (RadStruct(rss).processedFlag == 0) % Checks if the group of radials has to be processed
            % Make local copies of the radial files and insert the filenames strings in a cell array of strings.
            inp = 1;
            for mis=1:4
                if (RadStruct(rss).missingID(mis) == 0)
                    switch mis
                        case 1
                            inputfiles(inp) = {[work_folder RadStruct(rss).site1]};
                            inp = inp + 1;
                        case 2
                            inputfiles(inp) = {[work_folder RadStruct(rss).site2]};
                            inp = inp + 1;
                        case 3
                            inputfiles(inp) = {[work_folder RadStruct(rss).site3]};
                            inp = inp + 1;
                        otherwise
                            inputfiles(inp) = {[work_folder RadStruct(rss).site4]};
                            inp = inp + 1;
                    end
                end
            end
            
            % Makes the totals.
            tic;
            num_input = size(inputfiles,2);
            %             if (num_input > 1)
            RadStruct(rss).corruptedFlag = SingleTotal_v41(inputfiles, RadStruct(rss).TS, lon_lim, lat_lim, Grid, spatthresh, Radial_QC_params, Total_QC_params, range_cells_number, ncst, total_maps_folder, tuv_local_folder, tuv_serv_folder, src_folder, total_maps_server_folder, netcdf_thredds_folder);
            RadStruct(rss).processedFlag = 1 - RadStruct(rss).corruptedFlag; % Sets the current group of radials as processed (if not corrupted)
            save(strcat(src_folder, 'RadStruct.mat'), 'RadStruct');
            display(['[' datestr(now) '] - - ' 'RadStruct successfully saved.']);
            %             end
            
            totalsMakingTime = toc;
            
%             % Copies the used radial files into the destination folders
%             % (ordered by site)
%             tic
%             for i=1:num_input
%                 curr_file = char(inputfiles(i));
%                 switch curr_file(116)
%                     case 'V' % VIAR
%                         [cp_status] = copyfile(curr_file,[site1_folder curr_file(111:139)],'f');
%                     case 'T' % TINO
%                         [cp_status] = copyfile(curr_file,[site2_folder curr_file(111:139)],'f');
%                     case 'M' % MONT
%                         [cp_status] = copyfile(curr_file,[site3_folder curr_file(111:139)],'f');
%                     case 'x'
%                         [cp_status] = copyfile(curr_file,[site4_folder curr_file(111:139)],'f');
%                     case 'y'
%                         [cp_status] = copyfile(curr_file,[site5_folder curr_file(111:139)],'f');
%                     otherwise
%                         [cp_status] = copyfile(curr_file,[site6_folder curr_file(111:139)],'f');
%                 end
%             end
%             filesCopyTime = toc;
            
%             % Clean the radial local working folder
%             rmdir([radProc_folder 'dati-radar'],'s');
            
            clear inputfiles;
        end
    end
    % Periodically saves incremental copies of RadStruct
    now = clock;
    if (mod(now(4),4) == 0) && (now(5) < 16)
        incr_str = strcat((int2str(now(1))), '_', (int2str(now(2))), '_', (int2str(now(3))), '_', (int2str(now(4))), '_', (int2str(now(5))));
        save(strcat(incr_folder, 'RadStruct_', incr_str, '.mat'), 'RadStruct');
    end
    
    % Reset start and stop times for radial processing
    start_str = '0';
    stop_str = '0';
    
    % Pause
    pause(900);
    
end

%%

