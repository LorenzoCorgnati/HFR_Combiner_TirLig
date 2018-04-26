%% SingleTotal_v41.m
% This function generates the total maps by combining radial data.
% The 2 release operates cleaning on radial data in order to limit range
% and bearing intervals of validity and on total data in order to avoid the
% combination of radials with an intersection angle below a certain
% threshold.

% This version runs functions for storing hourly radials and totals in
% netCDF format using m/s as current measurement unit.

% This release stores the tuv and the netCDF files in a folder structure
% with the scheme yyyy/yyyy_mm/yyyy_mm_dd/.

% The release 3.0 sets U_std, V_std, UVcov and GDOP fields in the netCDF
% total data.

% The release 3.1 generates netCDF files for each input radial file.

% The release 3.2 saves radial and total netCDF files in a folder structure
% compliant to the ISMAR THREDDS requirements.

% The 4.0 release implements the data and metadata structure compliant to
% CMEMS needs and performs QC tests. In this release, a new graphic
% visualisation for total maps is implemented.

% The 4.1 release implements the data and metadata structure compliant to
% CMEMS-INSTAC requirements and performs QC tests.

% This version is designed for HFR_Combiner_TirLig_v21 and next releases.

% Author: Lorenzo Corgnati
% Date: March 15, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [sT_err] = SingleTotal_v41(ruv_files, when, lon_lim, lat_lim, Grid, spatthresh, Radial_QC_params, Total_QC_params, range_cells_number, mask_cst, dest_maps, dest_tuvs, dest_serv_tuvs, src_folder, serv_maps, serv_netcdf)

display(['[' datestr(now) '] - - ' 'SingleTotal_v41.m started.']);

sT_err = 0;

warning('off', 'all');

% Radial data load (it makes two attempts)
try
    display(['[' datestr(now) '] - - ' 'loadRDLfile loading ...']);
    RADIAL = loadRDLFile(ruv_files, 'false', 'warning');
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    sT_err = 1;
end

% try
%     RADIAL = loadRDLFile(ruv_files, 'false', 'warning');
%     sT_err = 0;
% catch err
%     sT_err = 1;
% end

szR = size(RADIAL, 1);

% Radial data bearing sector cleaning (if Ideal Pattern)
if (sT_err == 0)
    for Rds=1:szR
        if (strcmp(RADIAL(Rds,1).Type, 'RDLIdeal'))
            tmpR = RADIAL(Rds,1);
            RADIAL(Rds,1) = cleanBearingSector(tmpR);
        end
    end
end

% Retrieve the Year, Month and Day folder names for the current netCDF file
year = when(1:4);
month = when(6:7);
day = when(9:10);
yearFolder = year;
monthFolder = [year '/' year '_' month];
dayFolder = [monthFolder '/' year '_' month '_' day];

% Saves Radials in netCDF format
if (sT_err == 0)
    for rd_idx=1:length(RADIAL)
        % Check if the proper folders for radial data are existing.
        % If not, create them.
        %     if (exist([path_dest_rad RADIAL(rd_idx,1).SiteName '/' yearFolder], 'dir') ~= 7)
        %         mkdir([path_dest_rad RADIAL(rd_idx,1).SiteName '/' yearFolder]);
        %     end
        %     if (exist([path_dest_rad RADIAL(rd_idx,1).SiteName '/' monthFolder], 'dir') ~= 7)
        %         mkdir([path_dest_rad RADIAL(rd_idx,1).SiteName '/' monthFolder]);
        %     end
        %     if (exist([path_dest_rad RADIAL(rd_idx,1).SiteName '/' dayFolder], 'dir') ~= 7)
        %         mkdir([path_dest_rad RADIAL(rd_idx,1).SiteName '/' dayFolder]);
        %     end
        
        % Set the netCDF destination folder paths for radial data
        path_dest_rad = [dest_tuvs(1:length(dest_tuvs)-12) 'Radials/netCDF/'];
        servTH_dest_rad = [serv_netcdf 'TirLig/Radials/v1.0/'];
        servRD_dest_rad = [dest_serv_tuvs(1:length(dest_serv_tuvs)-12) 'Radials/netCDF/'];
        
        % Version 1.0
        if (exist([servRD_dest_rad 'v1.0/' RADIAL(rd_idx,1).SiteName '/' yearFolder], 'dir') ~= 7)
            mkdir([servRD_dest_rad 'v1.0/' RADIAL(rd_idx,1).SiteName '/' yearFolder]);
        end
        if (exist([servRD_dest_rad 'v1.0/' RADIAL(rd_idx,1).SiteName '/' monthFolder], 'dir') ~= 7)
            mkdir([servRD_dest_rad 'v1.0/' RADIAL(rd_idx,1).SiteName '/' monthFolder]);
        end
        if (exist([servRD_dest_rad 'v1.0/' RADIAL(rd_idx,1).SiteName '/' dayFolder], 'dir') ~= 7)
            mkdir([servRD_dest_rad 'v1.0/' RADIAL(rd_idx,1).SiteName '/' dayFolder]);
        end
        
        if (exist([servTH_dest_rad RADIAL(rd_idx,1).SiteName '/Last/' year '_' month '_' day], 'dir') ~= 7)
            mkdir([servTH_dest_rad RADIAL(rd_idx,1).SiteName '/Last/' year '_' month '_' day]);
        end
        
        if (exist([path_dest_rad 'v1.0/' RADIAL(rd_idx,1).SiteName '/Last/' year '_' month '_' day], 'dir') ~= 7)
            mkdir([path_dest_rad 'v1.0/' RADIAL(rd_idx,1).SiteName '/Last/' year '_' month '_' day]);
        end
        
        % Version 2.1
        if (exist([servRD_dest_rad 'v2.1/' RADIAL(rd_idx,1).SiteName '/' yearFolder], 'dir') ~= 7)
            mkdir([servRD_dest_rad 'v2.1/' RADIAL(rd_idx,1).SiteName '/' yearFolder]);
        end
        if (exist([servRD_dest_rad 'v2.1/' RADIAL(rd_idx,1).SiteName '/' monthFolder], 'dir') ~= 7)
            mkdir([servRD_dest_rad 'v2.1/' RADIAL(rd_idx,1).SiteName '/' monthFolder]);
        end
        if (exist([servRD_dest_rad 'v2.1/' RADIAL(rd_idx,1).SiteName '/' dayFolder], 'dir') ~= 7)
            mkdir([servRD_dest_rad 'v2.1/' RADIAL(rd_idx,1).SiteName '/' dayFolder]);
        end
        
        if (exist([path_dest_rad 'v2.1/' RADIAL(rd_idx,1).SiteName '/Last/' year '_' month '_' day], 'dir') ~= 7)
            mkdir([path_dest_rad 'v2.1/' RADIAL(rd_idx,1).SiteName '/Last/' year '_' month '_' day]);
        end
        
        % Create netCDF
        sT_err = Radial2netCDF_v11(RADIAL(rd_idx,1), range_cells_number, [path_dest_rad 'v1.0/' RADIAL(rd_idx,1).SiteName '/Last/' year '_' month '_' day '/'], [servTH_dest_rad RADIAL(rd_idx,1).SiteName '/Last/' year '_' month '_' day '/'], [servRD_dest_rad 'v1.0/' RADIAL(rd_idx,1).SiteName '/' dayFolder '/']);
        display(['[' datestr(now) '] - - ' RADIAL(rd_idx,1).SiteName '_' when ' radial netCDF v1.0 file successfully created and stored.']);
        [sT_err, sitePatternDate(rd_idx,:)] = Radial2netCDF_v21(RADIAL(rd_idx,1), range_cells_number, Radial_QC_params, [path_dest_rad 'v2.1/' RADIAL(rd_idx,1).SiteName '/Last/' year '_' month '_' day '/'], [servRD_dest_rad 'v2.1/' RADIAL(rd_idx,1).SiteName '/' dayFolder '/']);
        display(['[' datestr(now) '] - - ' RADIAL(rd_idx,1).SiteName '_' when ' radial netCDF v2.1 file successfully created and stored.']);
    end
end

% Evaluate last pattern date
if (sT_err == 0)
    for rd_idx=1:length(RADIAL)
        patternTS(rd_idx) = datenum([str2double(sitePatternDate(rd_idx, 1:5)) str2double(sitePatternDate(rd_idx, 6:8)) str2double(sitePatternDate(rd_idx, 9:11)) str2double(sitePatternDate(rd_idx, 13:15)) str2double(sitePatternDate(rd_idx, 16:18)) str2double(sitePatternDate(rd_idx, 19:20))]);
    end
    lastPatternDate = max(patternTS);
    lastPatternVec = datevec(lastPatternDate);
    lastPatternStr = [datestr(lastPatternVec, 'yyyy-mm-dd') 'T' datestr(lastPatternVec, 'HH:MM:SS') 'Z'];
end

% Totals generation
if (sT_err == 0)
    TS = RADIAL(1,1).TimeStamp;
end

if (not(exist('TS', 'var')))
    sT_err = 1;
end

if (not(exist('RADIAL', 'var')))
    sT_err = 1;
end

if ((sT_err == 0) && (isnan(TS)))
    sT_err = 1;
end

if (sT_err == 0)
    try
        display(['[' datestr(now) '] - - ' 'makeTotals combining ...']);
        [TUV,R] = makeTotals(RADIAL, 'Grid', Grid, 'TimeStamp', TS, 'spatthresh', spatthresh, 'tempthresh', 1/24);
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        sT_err = 1;
    end
end

if (sT_err == 0)
    % Totals setting on a regular grid
    [TUVgrid,DIM,I] = gridTotals( TUV, 'true', 'true');
    
    % Totals masking
    [TUVmask,I] = maskTotals(TUVgrid,mask_cst,0);
    
    % Totals cleaning for graphic plot
    gdop_thr = sqrt(4);
    maxspd_T = 120;
    [TUVclean,I] = cleanTotals(TUVmask,maxspd_T,{'GDOPMaxOrthog','TotalErrors',gdop_thr});
    
    % Contributing site coordinates and codes retrieval
    for Rds=1:szR
        contrSiteLon(Rds) = RADIAL(Rds,1).SiteOrigin(1);
        contrSiteLat(Rds) = RADIAL(Rds,1).SiteOrigin(2);
        %         contrSiteCode(Rds,:) = [RADIAL(Rds,1).SiteName '_' RADIAL(Rds,1).Type(1:4)];
        contrSiteCode(Rds,:) = RADIAL(Rds,1).SiteName;
    end
    
    % Set the netCDF destination folder paths for total data
    path_dest_tot_v1 = [dest_tuvs(1:length(dest_tuvs)-5) 'netCDF/v1.0/Last/'];
    path_dest_tot_v2 = [dest_tuvs(1:length(dest_tuvs)-5) 'netCDF/v2.1/Last/'];
    servTH_dest_tot = [serv_netcdf 'TirLig/Totals/v1.0/Last/'];
    servRD_dest_tot = [dest_serv_tuvs(1:length(dest_serv_tuvs)-5) 'netCDF/'];
    
    % TUV saving
    U_nan = size(find((isnan(TUVmask.U))==1),1);
    V_nan = size(find((isnan(TUVmask.V))==1),1);
    U_size = size(TUVmask.U,1);
    V_size = size(TUVmask.V,1);
    
    if ((U_nan < U_size) && (V_nan < V_size))
        try
            % Check if the proper folders for total data are existing.
            % If not, create them.
            if (exist([dest_tuvs yearFolder], 'dir') ~= 7)
                mkdir([dest_tuvs yearFolder]);
            end            
            
            if (exist([dest_tuvs monthFolder], 'dir') ~= 7)
                mkdir([dest_tuvs monthFolder]);
            end
            
            if (exist([dest_tuvs dayFolder], 'dir') ~= 7)
                mkdir([dest_tuvs dayFolder]);
            end
            
            if (exist([dest_serv_tuvs yearFolder], 'dir') ~= 7)
                mkdir([dest_serv_tuvs yearFolder]);
            end
            
            if (exist([dest_serv_tuvs monthFolder], 'dir') ~= 7)
                mkdir([dest_serv_tuvs monthFolder]);
            end              
            
            if (exist([dest_serv_tuvs dayFolder], 'dir') ~= 7)
                mkdir([dest_serv_tuvs dayFolder]);
            end
            
            % Version 1.0
            if (exist([servRD_dest_tot 'v1.0/' yearFolder], 'dir') ~= 7)
                mkdir([servRD_dest_tot 'v1.0/' yearFolder]);
            end
            
            if (exist([servRD_dest_tot 'v1.0/' monthFolder], 'dir') ~= 7)
                mkdir([servRD_dest_tot 'v1.0/' monthFolder]);
            end
            
            if (exist([servRD_dest_tot 'v1.0/' dayFolder], 'dir') ~= 7)
                mkdir([servRD_dest_tot 'v1.0/' dayFolder]);
            end
            
            % Version 2.1
            if (exist([servRD_dest_tot 'v2.1/' yearFolder], 'dir') ~= 7)
                mkdir([servRD_dest_tot 'v2.1/' yearFolder]);
            end
            
            if (exist([servRD_dest_tot 'v2.1/' monthFolder], 'dir') ~= 7)
                mkdir([servRD_dest_tot 'v2.1/' monthFolder]);
            end
            
            if (exist([servRD_dest_tot 'v2.1/' dayFolder], 'dir') ~= 7)
                mkdir([servRD_dest_tot 'v2.1/' dayFolder]);
            end            
            
            % Version 1.0
            if (exist([path_dest_tot_v1 year '_' month '_' day], 'dir') ~= 7)
                mkdir([path_dest_tot_v1 year '_' month '_' day]);
            end
            
            % Version 2.1
            if (exist([path_dest_tot_v2 year '_' month '_' day], 'dir') ~= 7)
                mkdir([path_dest_tot_v2 year '_' month '_' day]);
            end
            
            if (exist([servTH_dest_tot year '_' month '_' day], 'dir') ~= 7)
                mkdir([servTH_dest_tot year '_' month '_' day]);
            end
            
            save(strcat([dest_tuvs dayFolder '/'], sprintf('%s.mat', when)), 'TUVmask');
            display(['[' datestr(now) '] - - ' when ' mat file successfully saved locally.']);
            [cp_status] = copyfile((strcat([dest_tuvs dayFolder '/'], sprintf('%s.mat', when))),(strcat([dest_serv_tuvs dayFolder '/'], sprintf('%s.mat', when))),'f');
            display(['[' datestr(now) '] - - ' when ' mat file successfully saved on RadarDisk.']);
        catch err
            display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            sT_err = 1;
        end
        
        % Totals maps generation
        %         openfig(strcat(src_folder, 'AoI.fig'), 'new', 'invisible');
        m_proj('transverse mercator','longitudes',lon_lim,'latitudes',lat_lim);
        m_gshhs_f('patch',[.7 .7 .7],'edgecolor','none');
        m_grid('box', 'fancy', 'tickdir', 'in', 'xlabeldir','end','fontsize',10);
        
%         [X,Y]=m_ll2xy(9.6533333,44.1458333);
%         line(X,Y,'marker','square','markersize',4,'color','r');
%         text(X,Y,' MONT','vertical','top');
        
        [X,Y]=m_ll2xy(9.6593000,44.1435167);
        line(X,Y,'marker','square','markersize',4,'color','r');
        text(X,Y,' PCOR','vertical','top');
        
        [X,Y]=m_ll2xy(9.8492167,44.0263667);
        line(X,Y,'marker','square','markersize',4,'color','r');
        text(X,Y,' TINO','vertical','top');
        
        [handles,TimeIndex] = plotData(TUVclean,'m_vec',TUV.TimeStamp, 0.005, 'headangle', 25, 'headlength', 5.5, 'shaftwidth', 1);
        
        when_fig = datestr(TUV.TimeStamp, 'dd-mmm-yyyy HH:MM:SS');
        title(sprintf('%s',when_fig));
        
        pos2 = [0.3  0.035  0.4  0.015];
        %ha2 = axes('position', pos2);
        hc = colorbar;
        caxis([0 100]);
        hc.Location = 'south';
        hc.Position = pos2;
        title(hc,'[cm/s]')
        
        % Saves the map in the local folder and copies it to the servers
        % (it makes two attempts)
        local_maps = [dest_tuvs(1:length(dest_tuvs)-5) 'Total_maps/'];
        try
            saveas(gcf,[local_maps when '.jpg']);
            display(['[' datestr(now) '] - - ' when ' map successfully saved locally.']);
        catch err
            display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            sT_err = 1;
        end
        
        try
            [cp_status] = copyfile(strcat(local_maps, sprintf('%s.jpeg', when)),strcat(dest_maps, sprintf('%s.jpeg', when)),'f');
            display(['[' datestr(now) '] - - ' when ' map successfully saved on RadarDisk.']);
        catch err
            display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            sT_err = 1;
        end
        
        try
            %print(strcat(serv_maps, sprintf('%s.jpeg', when)),'-djpeg');
            [cp_status] = copyfile(strcat(local_maps, sprintf('%s.jpeg', when)),strcat(serv_maps, sprintf('%s.jpeg', when)),'f');
            display(['[' datestr(now) '] - - ' when ' map successfully saved on web server.']);
        catch err
            display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            sT_err = 1;
        end
        close;
        
        % Saves Total in netCDF format
        sT_err = Total2netCDF_v21(TUVmask, Grid, lon_lim, lat_lim, maxspd_T, spatthresh, gdop_thr, contrSiteLat, contrSiteLon, contrSiteCode, mask_cst, [path_dest_tot_v1 year '_' month '_' day '/'], [servTH_dest_tot year '_' month '_' day '/'], [servRD_dest_tot 'v1.0/' dayFolder '/']);
        display(['[' datestr(now) '] - - ' when ' total netCDF v1.0 file successfully created and stored.']);
        sT_err = Total2netCDF_v31(TUVmask, Grid, lon_lim, lat_lim, spatthresh, Total_QC_params, contrSiteLat, contrSiteLon, contrSiteCode, mask_cst, [path_dest_tot_v2 year '_' month '_' day '/'], [servRD_dest_tot 'v2.1/' dayFolder '/'], lastPatternStr);
        display(['[' datestr(now) '] - - ' when ' total netCDF v2.1 file successfully created and stored.']);
    end
end


return