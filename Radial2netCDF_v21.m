%% Radial2netCDF_v21.m
% This function converts hourly radial files from Codar format in netCDF
% format according to the standard struvture defined within INCREASE and
% JERICO-NEXT projects.
% This version is designed for HFR_Combiner_TirLig_v21 and next releases.

% The release 2.1 implements the data and metadata structure defined in
% INCREASE project and suited for CMEMS-INSTAC requirements.

% INPUT:
%         HFRP_RUV: structure containing the radial data to be loaded using the
%                           HFR Progs library function (loadRDLFile)
%         range_cells_number: default number of range cells
%         Radial_QC_params: structure containing parameters for radial QC
%         tests
%         dest_nc_local : local folder where to save the netCDF file
%         dest_nc_rad : RadarDisk folder where to save the netCDF file

% OUTPUT:
%         R2C_err: error flag (0 = correct, 1 = error)
%         PatternDate: string containing the pattern measurement date


% Author: Lorenzo Corgnati
% Date: May 15, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [R2C_err, PatternDate] = Radial2netCDF_v21(HFRP_RUV, range_cells_number, Radial_QC_params, dest_nc_local, dest_nc_rd)

display(['[' datestr(now) '] - - ' 'Radial2netCDF_v21.m started.']);

R2C_err = 0;

warning('off', 'all');

%% Define the maximum number of allowed sites and of xlink references

maxsite = 50;
refmax = 1;

%%

%% Set the input files

rfile = HFRP_RUV.FileName{1,1};

%%

%% Build the output file name and the citation string

% Output file name
fileType = rfile(length(rfile)-25);
siteCode = rfile(length(rfile)-23:length(rfile)-20);
fileTime = rfile(length(rfile)-18:length(rfile)-4);

ncfile_local = [dest_nc_local 'RDL_' fileType '_' siteCode '_' fileTime '.nc'];
ncfile_rd = [dest_nc_rd 'RDL_' fileType '_' siteCode '_' fileTime '.nc'];

% Citation and distribution string
citation_str = ['These data were collected and made freely available by the Copernicus project and the programs that contribute to it. Data collected and processed by CNR-ISMAR within RITMARE, Jerico-Next and IMPACT projects -  Year ' fileTime(1:4)];
distribution_str = 'These data follow Copernicus standards; they are public and free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.';

%%

%% Define EDIOS and EDMO codes, site code, platform code, id and metadata resources

EDIOS_Series_ID = 'HFR_TirLig';
EDIOS_Platform_ID = siteCode;
EDMO_code = 134;
site_code = EDIOS_Series_ID;
platform_code = [EDIOS_Series_ID '_' EDIOS_Platform_ID];
id = [EDIOS_Series_ID '_' EDIOS_Platform_ID '_' strrep(fileTime(1:10), '_', '-') '_' fileTime(12:13) 'Z'];
TDS_catalog = 'http://150.145.136.27:8080/thredds/HF_RADAR/TirLig/TirLig_catalog.html';
xlink = ['<sdn_reference xlink:href="' TDS_catalog '" xlink:role="" xlink:type="URL"/>'];

%%

%% Load radial data
try
    % Load data using native Matlab load function
    r = load(rfile);
    
    R = struct( ...
        'lond', r(:, 1), ...
        'latd', r(:, 2), ...
        'velu', r(:, 3), ...
        'velv', r(:, 4), ...
        'vflg', r(:, 5), ...
        'espc', r(:, 6), ...
        'etmp', r(:, 7), ...
        'maxv', r(:, 8), ...
        'minv', r(:, 9), ...
        'ersc', r(:,10), ...
        'ertc', r(:,11), ...
        'xdst', r(:,12), ...
        'ydst', r(:,13), ...
        'rnge', r(:,14), ...
        'bear', r(:,15), ...
        'velo', r(:,16), ...
        'head', r(:,17), ...
        'sprc', r(:,18) ...
        );
    clear r;
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    R2C_err = 1;
end

%%

%% Retrieve information from original CODAR fields
if (R2C_err == 0)
    try
        % Read the file header
        header = HFRP_RUV.OtherMetadata.Header;
        
        % Retrieve information from header
        TableType = '';
        TableColumns = '';
        TableColumnTypes = '';
        TableRows = '';
        for head_idx=1:length(header)
            splitLine = regexp(header{head_idx}, ' ', 'split');
            
            % Range Resolution
            if (strcmp(splitLine{1}, '%RangeResolutionKMeters:'))
                RangeResolutionKMeters = str2num(splitLine{2});
            end
            if(strcmp(splitLine{1}, '%AngularResolution:'))
                AngularResolution = str2num(splitLine{2});
            end
            if(strcmp(splitLine{1}, '%Origin:'))
                JN = 1;
                try
                    Origin = [str2num(splitLine{3}), str2num(splitLine{7})]; % for Jerico-Next
                catch err
                    JN = 0;
                end
                if(JN == 0)
                    Origin = [str2num(splitLine{3}), str2num(splitLine{6})]; % for SSD
                end
            end
            if(strcmp(splitLine{1}, '%CTF:'))
                CTF = strrep(header{head_idx}(length('%CTF:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%FileType:'))
                FileType = strrep(header{head_idx}(length('%FileType:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%LLUVSpec:'))
                LLUVSpec = strrep(header{head_idx}(length('%LLUVSpec:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%UUID:'))
                UUID = strrep(header{head_idx}(length('%UUID:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%Manufacturer:'))
                manufacturer = strrep(header{head_idx}(length('%Manufacturer:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%Site:'))
                Site = strrep(header{head_idx}(length('%Site:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%TimeZone:'))
                TimeZone = strrep(header{head_idx}(length('%TimeZone:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%TimeCoverage:'))
                TimeCoverage = strrep(header{head_idx}(length('%TimeCoverage:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%Origin:'))
                Origin_str = strrep(header{head_idx}(length('%Origin:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%GreatCircle:'))
                GreatCircle = strrep(header{head_idx}(length('%GreatCircle:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%GeodVersion:'))
                GeodVersion = strrep(header{head_idx}(length('%GeodVersion:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%LLUVTrustData:'))
                LLUVTrustData = strtok(strrep(header{head_idx}(length('%LLUVTrustData:')+2:length(header{head_idx})), '"', ''), '%%');
            end
            if(strcmp(splitLine{1}, '%RangeStart:'))
                RangeStart = strrep(header{head_idx}(length('%RangeStart:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%RangeEnd:'))
                RangeEnd = strrep(header{head_idx}(length('%RangeEnd:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%RangeResolutionKMeters:'))
                RangeResolutionKMeters_str = strrep(header{head_idx}(length('%RangeResolutionKMeters:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%AntennaBearing:'))
                AntennaBearing = strrep(header{head_idx}(length('%AntennaBearing:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%ReferenceBearing:'))
                ReferenceBearing = strrep(header{head_idx}(length('%ReferenceBearing:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%AngularResolution:'))
                AngularResolution_str = strrep(header{head_idx}(length('%AngularResolution:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%SpatialResolution:'))
                SpatialResolution = strrep(header{head_idx}(length('%SpatialResolution:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%PatternType:'))
                PatternType = strrep(header{head_idx}(length('%PatternType:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%PatternDate:'))
                PatternDate = strrep(header{head_idx}(length('%PatternDate:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%PatternResolution:'))
                PatternResolution = strrep(header{head_idx}(length('%PatternResolution:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%TransmitCenterFreqMHz:'))
                TransmitCenterFreqMHz = strrep(header{head_idx}(length('%TransmitCenterFreqMHz:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%DopplerResolutionHzPerBin:'))
                DopplerResolutionHzPerBin = strrep(header{head_idx}(length('%DopplerResolutionHzPerBin:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%FirstOrderMethod:'))
                FirstOrderMethod = strrep(header{head_idx}(length('%FirstOrderMethod:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%BraggSmoothingPoints:'))
                BraggSmoothingPoints = strrep(header{head_idx}(length('%BraggSmoothingPoints:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%CurrentVelocityLimit:'))
                CurrentVelocityLimit = strrep(header{head_idx}(length('%CurrentVelocityLimit:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%BraggHasSecondOrder:'))
                BraggHasSecondOrder = strrep(header{head_idx}(length('%BraggHasSecondOrder:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%RadialBraggPeakDropOff:'))
                RadialBraggPeakDropOff = strrep(header{head_idx}(length('%RadialBraggPeakDropOff:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%RadialBraggPeakNull:'))
                RadialBraggPeakNull = strrep(header{head_idx}(length('%RadialBraggPeakNull:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%RadialBraggNoiseThreshold:'))
                RadialBraggNoiseThreshold = strrep(header{head_idx}(length('%RadialBraggNoiseThreshold:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%PatternAmplitudeCorrections:'))
                PatternAmplitudeCorrections = strrep(header{head_idx}(length('%PatternAmplitudeCorrections:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%PatternPhaseCorrections:'))
                PatternPhaseCorrections = strrep(header{head_idx}(length('%PatternPhaseCorrections:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%PatternAmplitudeCalculations:'))
                PatternAmplitudeCalculations = strrep(header{head_idx}(length('%PatternAmplitudeCalculations:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%PatternPhaseCalculations:'))
                PatternPhaseCalculations = strrep(header{head_idx}(length('%PatternPhaseCalculations:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%RadialMusicParameters:'))
                RadialMusicParameters = strrep(header{head_idx}(length('%RadialMusicParameters:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%MergedCount:'))
                MergedCount = strrep(header{head_idx}(length('%MergedCount:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%RadialMinimumMergePoints:'))
                RadialMinimumMergePoints = strrep(header{head_idx}(length('%RadialMinimumMergePoints:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%FirstOrderCalc:'))
                FirstOrderCalc = strrep(header{head_idx}(length('%FirstOrderCalc:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%MergeMethod:'))
                MergeMethod = strrep(header{head_idx}(length('%MergeMethod:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%PatternMethod:'))
                PatternMethod = strrep(header{head_idx}(length('%PatternMethod:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%TransmitSweepRateHz:'))
                TransmitSweepRateHz = strrep(header{head_idx}(length('%TransmitSweepRateHz:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%TransmitBandwidthKHz:'))
                TransmitBandwidthKHz = strrep(header{head_idx}(length('%TransmitBandwidthKHz:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%SpectraRangeCells:'))
                SpectraRangeCells = strrep(header{head_idx}(length('%SpectraRangeCells:')+2:length(header{head_idx})), '"', '');
            end
            if(strcmp(splitLine{1}, '%SpectraDopplerCells:'))
                SpectraDopplerCells = strrep(header{head_idx}(length('%SpectraDopplerCells:')+2:length(header{head_idx})), '"', '');
            end
            if((strcmp(splitLine{1}, '%TableType:')) & (isempty(TableType)))
                TableType = strrep(header{head_idx}(length('%TableType:')+2:length(header{head_idx})), '"', '');
            end
            if((strcmp(splitLine{1}, '%TableColumns:')) & (isempty(TableColumns)))
                TableColumns = strrep(header{head_idx}(length('%TableColumns:')+2:length(header{head_idx})), '"', '');
            end
            if((strcmp(splitLine{1}, '%TableColumnTypes:')) & (isempty(TableColumnTypes)))
                TableColumnTypes = strrep(header{head_idx}(length('%TableColumnTypes:')+2:length(header{head_idx})), '"', '');
            end
            if((strcmp(splitLine{1}, '%TableRows:')) & (isempty(TableRows)))
                TableRows = strrep(header{head_idx}(length('%TableRows:')+2:length(header{head_idx})), '"', '');
            end
        end
        
        % Evaluate last pattern date
        patternTS = datenum([str2double(PatternDate(1:5)) str2double(PatternDate(6:8)) str2double(PatternDate(9:11)) str2double(PatternDate(13:15)) str2double(PatternDate(16:18)) str2double(PatternDate(19:20))]);
        lastPatternVec = datevec(patternTS);
        lastPatternStr = [datestr(lastPatternVec, 'yyyy-mm-dd') 'T' datestr(lastPatternVec, 'HH:MM:SS') 'Z'];
        
        % Time Stamp
        R.time = HFRP_RUV.TimeStamp;
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        R2C_err = 1;
    end
    
    %%
    
    %% Retrieve time coverage period and metadata date stamp
    try
        coverageStart = addtodate(HFRP_RUV.TimeStamp, -30, 'minute');
        timeCoverageStart = [datestr(coverageStart, 'yyyy-mm-dd') 'T' datestr(coverageStart, 'HH:MM:SS') 'Z'];
        coverageEnd = addtodate(HFRP_RUV.TimeStamp, 30, 'minute');
        timeCoverageEnd = [datestr(coverageEnd, 'yyyy-mm-dd') 'T' datestr(coverageEnd, 'HH:MM:SS') 'Z'];
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        R2C_err = 1;
    end
    
    %%
    
    %% Retrieve the file creation datetime (i.e. metadata date stamp)
    try
        dateCreated = [datestr(now, 'yyyy-mm-dd') 'T' datestr(now, 'HH:MM:SS') 'Z'];
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        R2C_err = 1;
    end
    
    %%
    
    %% Set reference time
    try
        timeref = datenum(1950,1,1);
        time_units = [datestr(timeref, 'yyyy-mm-dd') 'T' datestr(timeref, 'HH:MM:SS') 'Z'];
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        R2C_err = 1;
    end
    
    %%
    
    
    %% Set collection time
    try
        ts = HFRP_RUV.TimeStamp;
        time_coll = [datestr(ts, 'yyyy-mm-dd') 'T' datestr(ts, 'HH:MM:SS') 'Z'];
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        R2C_err = 1;
    end
    
    %%
    
    %% Generate grid domain
    try
        
        %     range_dim = min(R.rnge):RangeResolutionKMeters:max(R.rnge);
        
        % Force the range variable to have the same number of values, in order to have
        % all radial files with the same range dimension. It's crucial for the
        % THREDDS time aggregation.
        range_dim = 0:RangeResolutionKMeters:range_cells_number;
        
        bearing_dim = 0:AngularResolution:360;
        
        [bearing, range] = meshgrid(bearing_dim, range_dim);
        
        % % CHECK - see that our grid captures the sample space right
        % figure
        % plot( range, bearing, 'go' )
        % hold on
        % plot( R.rnge, R.bear, 'b.')
        
        % Generate lat/lon from bearing and range
        % This is the most accurate method of defining lat/lon values for given
        % range & bearing. Other methods of interpolating and extrapolating using
        % griddata or otherwise from range, bearing, lon and lat are inconsistent
        % and generally not reliable from file to file. The reason we need to do
        % this at all is because there cannot be missing values in auxillary
        % coordinate variables.
        [latd, lond] = reckon(Origin(1), Origin(2), km2deg(range), bearing);
        
        % % CHECK
        % figure
        % plot(lond, latd, 'ko')
        % hold on
        % plot(R.lond, R.latd, 'r.')
        % daspect([1.2 1 1])
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        R2C_err = 1;
    end
    
    %%
    
    %% Determine the mapping index from data to full grid
    try
        range_map_idx = nan(1, numel(R.rnge));
        bearing_map_idx = nan(1, numel(R.bear));
        for I = 1 : numel(R.rnge)
            [~, range_map_idx(I)] = min(abs(range_dim - R.rnge(I)));
            [~, bearing_map_idx(I)] = min(abs(bearing_dim - R.bear(I)));
        end
        clear I;
        rb_map_idx = sub2ind(size(range), range_map_idx, bearing_map_idx);
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        R2C_err = 1;
    end
    
    %%
    
    %% Populate gridded matrices
    try
        % Eastward Velocity
        velu = nan(size(range));
        velu(rb_map_idx) = R.velu./100;
        %     velu(rb_map_idx) = 0;
        
        % Northward Velocity
        velv = nan(size(range));
        velv(rb_map_idx) = R.velv./100;
        %    velv(rb_map_idx) = 0;
        
        % % CHECK
        % figure
        % quiver(R.rnge, R.bear, R.velu, R.velv, 0, 'k')
        % hold on
        % quiver(range, bearing, velu, velv, 0, 'r')
        %
        % figure
        % quiver(R.lond, R.latd, R.velu./100, R.velv./100, 0, 'k')
        % daspect([1.2 1 1])
        % hold on
        % quiver(lond, latd, velu./100, velv./100, 0, 'r')
        
        % Current Speed - flip sign so positive velocity is away from radar
        % according to CF standard name 'radial_sea_water_velocity_away_from_instrument'.
        % CODAR reports positive speeds as toward the radar.
        velo = nan(size(range));
        velo(rb_map_idx) = -R.velo./100;
        
        % Current Heading
        head = nan(size(range));
        head(rb_map_idx) = R.head;
        
        % Current Bearing
        bear = nan(size(range));
        bear(rb_map_idx) = R.bear;
        
        % Vector On Water Flag
        owtr = nan(size(range));
        owtr(rb_map_idx) = R.vflg;
        
        % Spatial Quality
        espc = nan(size(range));
        espc(rb_map_idx) = R.espc./100;
        espc( espc==9.99 ) = NaN; % native bad-value
        
        % Temporal Quality
        etmp = nan(size(range));
        etmp(rb_map_idx) = R.etmp./100;
        etmp( etmp==9.99 ) = NaN; % native bad-value
        
        % Velocity Maximum - flip sign so positive velocity is away from radar
        % according to CF standard name 'radial_sea_water_velocity_away_from_instrument'.
        % CODAR reports positive speeds as toward the radar.
        maxv = nan(size(range));
        maxv(rb_map_idx) = -R.maxv/100;
        maxv( maxv==-9.99 ) = NaN; % native bad-value
        
        % Velocity Minimum - flip sign so positive velocity is away from radar
        % according to CF standard name 'radial_sea_water_velocity_away_from_instrument'.
        % CODAR reports positive speeds as toward the radar.
        minv = nan(size(range));
        minv(rb_map_idx) = -R.minv/100;
        minv( minv==-9.99 ) = NaN; % native bad-value
        
        % Spatial Count
        ersc = nan(size(range));
        ersc(rb_map_idx) = R.ersc;
        ersc(ersc>127) = 127;
        ersc( ersc==999 ) = NaN; % native bad-value
        
        % Temporal Count
        ertc = nan(size(range));
        ertc(rb_map_idx) = R.ertc;
        ertc(ertc>127) = 127;
        ertc( ertc==999 ) = NaN; % native bad-value
        
        % X-Distance
        xdst = nan(size(range));
        xdst(rb_map_idx) = R.xdst;
        
        % Y-Distance
        ydst = nan(size(range));
        ydst(rb_map_idx) = R.ydst;
        
        % Spectra Range Cell
        sprc = nan(size(range));
        sprc(rb_map_idx) = R.sprc;
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        R2C_err = 1;
    end
    
    %%
    
    %% Perform QC tests
    
    % Build the names of the files of the previous two hours
    [twoHoursBefore, oneHourBefore] = TimeStamp2TS(R.time);
    twoHoursBefore_rep1 = strrep(ncfile_local,fileTime,twoHoursBefore);
    Radial_QC_params.TempDerThr.hour2 = strrep(twoHoursBefore_rep1,fileTime(1:length(fileTime)-5),twoHoursBefore(1:length(twoHoursBefore)-5));
    oneHoursBefore_rep1 = strrep(ncfile_local,fileTime,oneHourBefore);
    Radial_QC_params.TempDerThr.hour1 = strrep(oneHoursBefore_rep1,fileTime(1:length(fileTime)-5),oneHourBefore(1:length(oneHourBefore)-5));
    oneHoursBefore_rep1_rd = strrep(ncfile_rd,fileTime,oneHourBefore);
    oneHoursBefore_rep2_rd = strrep(oneHoursBefore_rep1_rd,fileTime(1:length(fileTime)-5),oneHourBefore(1:length(oneHourBefore)-5));
    oneHoursBefore_rep3_rd = strrep(oneHoursBefore_rep2_rd,fileTime(1:length(fileTime)-8),oneHourBefore(1:length(oneHourBefore)-8));
    Radial_QC_params.TempDerThr.hour1_RD = strrep(oneHoursBefore_rep3_rd,fileTime(1:length(fileTime)-11),oneHourBefore(1:length(oneHourBefore)-11));
    
    [overall_QCflag, overWater_QCflag, varianceThreshold_QCflag, temporalDerivativeThreshold_QCflag, velocityThreshold_QCflag, medianFilter_QCflag, averageRadialBearing_QC_flag, radialVelocityMedianFiltered, radialCount_QC_flag] = RadialQCtests_v10(bear, lond, latd, owtr, etmp, head, velo, Radial_QC_params, Site(1:4));
    
    %%
    
    %% Scaling
    
    % Bearing and heading are only reported to 10ths and can be
    % reported as short unsigned integers when scaled. However, Bearing is a
    % coordinate variable and cannot be scaled by CF metadata convention.
    head = round(head.*10);
    
    %%
    
    %% Fill Values
    try
        % Type byte fill values
        ersc(isnan(ersc)) = netcdf.getConstant('NC_FILL_BYTE');
        ertc(isnan(ertc)) = netcdf.getConstant('NC_FILL_BYTE');
        sprc(isnan(sprc)) = netcdf.getConstant('NC_FILL_BYTE');
        
        % Type short fill values
        head(isnan(head)) = netcdf.getConstant('NC_FILL_SHORT');
        
        % Type float fill values
        velo(isnan(velo)) = netcdf.getConstant('NC_FILL_FLOAT');
        radialVelocityMedianFiltered(isnan(radialVelocityMedianFiltered)) = netcdf.getConstant('NC_FILL_FLOAT');
        velu(isnan(velu)) = netcdf.getConstant('NC_FILL_FLOAT');
        velv(isnan(velv)) = netcdf.getConstant('NC_FILL_FLOAT');
        espc(isnan(espc)) = netcdf.getConstant('NC_FILL_FLOAT');
        etmp(isnan(etmp)) = netcdf.getConstant('NC_FILL_FLOAT');
        maxv(isnan(maxv)) = netcdf.getConstant('NC_FILL_FLOAT');
        minv(isnan(minv)) = netcdf.getConstant('NC_FILL_FLOAT');
        xdst(isnan(xdst)) = netcdf.getConstant('NC_FILL_FLOAT');
        ydst(isnan(ydst)) = netcdf.getConstant('NC_FILL_FLOAT');
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        R2C_err = 1;
    end
    
    %%
    
    %% Set the time, position and depth quality flags
    % Time quality flag
    sdnTime_QCflag = 1;
    % Position quality flag
    sdnPosition_QCflag = netcdf.getConstant('NC_FILL_SHORT').*int16(ones(size(velu,1),size(velu,2),1));
    sdnPosition_QCflag(velu~=netcdf.getConstant('NC_FILL_FLOAT')) = 1;
    
    % Depth quality flag
    sdnDepth_QCflag = 1;
    
    %%
    
    %% Populate site latitude, site longitude and site code variables
    siteLat = Origin(1);
    siteLon = Origin(2);
    numSites = 1;
    siteLat((length(siteLat)+1):maxsite) = netcdf.getConstant('NC_FILL_FLOAT');
    siteLon((length(siteLon)+1):maxsite) = netcdf.getConstant('NC_FILL_FLOAT');
    siteCode(1,:) = Site(1:length(Site)-1);
    for sC_idx=2:maxsite
        siteCode(sC_idx,:) = '    ';
    end
    %     if(PatternType == 'Measured')
    %         siteCode = [Site(1:length(Site)-1) '_RDLM'];
    %     else
    %         siteCode = [Site(1:length(Site)-1) '_RDLI'];
    %     end
    
    %%
    
    
    %% Define dimensions
    
    % Deletes the eventually present netCDF file with the same name
    try
        delete(ncfile_local);
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        R2C_err = 1;
    end
    
    try
        % Keep the classic data model to increase compatibility since we aren't using
        % features specific to the NetCDF-4.0 model.
        mode = netcdf.getConstant('NETCDF4');
        mode = bitor(mode, netcdf.getConstant('CLASSIC_MODEL'));
        ncid = netcdf.create(ncfile_local, mode);
        hist_create = [time_coll ' data collected. ' dateCreated ' netCDF file created and sent to European HFR Node'];
        
        % Add dimensions (in order of T, Z, Y, X for CF/COARDS compliance)
        
        %         dimid_t = netcdf.defDim(ncid, 'TIME', netcdf.getConstant('unlimited'));
        dimid_t = netcdf.defDim(ncid, 'TIME', 1);
        dimid_bearing = netcdf.defDim(ncid, 'HEAD', numel( bearing_dim));
        dimid_range = netcdf.defDim(ncid, 'RNGE', numel( range_dim));
        dimid_depth = netcdf.defDim(ncid, 'DEPH', 1);
        dimid_maxsite = netcdf.defDim(ncid, 'MAXSITE', maxsite);
        dimid_refmax = netcdf.defDim(ncid, 'REFMAX', refmax);
        dimid_string4 = netcdf.defDim(ncid, 'STRING4', 4);
        dimid_string_site_code = netcdf.defDim(ncid, ['STRING' num2str(length(site_code))], length(site_code));
        dimid_string_platform_code = netcdf.defDim(ncid, ['STRING' num2str(length(platform_code))], length(platform_code));
        dimid_string_sdn_local_cdi_id = netcdf.defDim(ncid, ['STRING' num2str(length(id))], length(id));
        dimid_string_sdn_references = netcdf.defDim(ncid, ['STRING' num2str(length(TDS_catalog))], length(TDS_catalog));
        dimid_string_sdn_xlink = netcdf.defDim(ncid, ['STRING' num2str(length(xlink))], length(xlink));
        
        %dimid_lat = netcdf.defDim(ncid, 'lat', length(LatU));
        %dimid_lon = netcdf.defDim(ncid, 'lon', length(LonU));
        
        %%
        
        %% Add coordinate variables and attributes
        
        % For time coordinate variable, when using the standard or gregorian
        % calendar, units should be in seconds and shifted forward beyond
        % 1582/10/15 to avoid Julian -> Gregorian crossover and potential
        % problems with the udunits package, if used. Year length should be
        % exactly 365.2425 days long for the Gregorian/standard and proleptic
        % Gregorian calendar.
        %
        % MATLAB's datenum uses the proleptic gregorian calendar as defined
        % by the CF standard. When using standard_name attribute for time
        % coordinate be sure to include the calendar and units attributes.
        %
        % Use int for time data type because need 8 significant figures to
        % represent number of seconds in 1 year (31556952 seconds), float
        % isn't enough (7 sig. figures). Would need to use double beyond
        % int. Besides, don't need sub-second accuracy. Unix time will
        % overflow int in 2038. New epoch will be needed then to keep int
        % representation. Unix time epoch is: 1970-01-01 00:00:00Z = 0 seconds, 86,400 sec/day
        varid_t = netcdf.defVar(ncid, 'TIME', 'float', dimid_t);
        netcdf.putAtt(ncid, varid_t, 'long_name', 'Time of Measurement UTC');
        netcdf.putAtt(ncid, varid_t, 'standard_name', 'time');
        netcdf.putAtt(ncid, varid_t, 'units', ['days since ' time_units]);
        netcdf.putAtt(ncid, varid_t, 'calendar', 'Julian');
        netcdf.putAtt(ncid, varid_t, 'axis', 'T');
        netcdf.putAtt(ncid, varid_t, 'sdn_parameter_name', 'Elapsed time (since 1950-01-01T00:00:00Z)');
        netcdf.putAtt(ncid, varid_t, 'sdn_parameter_urn', 'SDN:P01::ELTJLD01');
        netcdf.putAtt(ncid, varid_t, 'sdn_uom_name', 'Days');
        netcdf.putAtt(ncid, varid_t, 'sdn_uom_urn', 'SDN:P06::UTAA');
        netcdf.putAtt(ncid, varid_t, 'ancillary_variables', 'TIME_SEADATANET_QC');
        
        % Bearing (arbitrary 'y' dimension)
        varid_bearing = netcdf.defVar(ncid, 'HEAD', 'float', dimid_bearing);
        netcdf.putAtt(ncid, varid_bearing, 'axis', 'Y');
        netcdf.putAtt(ncid, varid_bearing, 'long_name', 'Bearing Away From Instrument');
        netcdf.putAtt(ncid, varid_bearing, 'units', 'degrees_true');
        netcdf.putAtt(ncid, varid_bearing, 'sdn_parameter_name', 'Orientation (horizontal relative to true north) of measurement device {heading}');
        netcdf.putAtt(ncid, varid_bearing, 'sdn_parameter_urn', 'SDN:P01::HEADCM01');
        netcdf.putAtt(ncid, varid_bearing, 'sdn_uom_name', 'Degrees true');
        netcdf.putAtt(ncid, varid_bearing, 'sdn_uom_urn', 'SDN:P06::UABB');
        netcdf.putAtt(ncid, varid_bearing, 'ancillary_variables', 'POSITION_SEADATANET_QC');
        
        % Range (arbitrary 'x' dimension)
        varid_range = netcdf.defVar(ncid, 'RNGE', 'float', dimid_range);
        netcdf.putAtt(ncid, varid_range, 'axis', 'X');
        netcdf.putAtt(ncid, varid_range, 'long_name', 'Range Away From Instrument');
        netcdf.putAtt(ncid, varid_range, 'units', 'km');
        netcdf.putAtt(ncid, varid_range, 'sdn_parameter_name', 'Range (from fixed reference point) by unspecified GPS system');
        netcdf.putAtt(ncid, varid_range, 'sdn_parameter_urn', 'SDN:P01::RIFNAX01');
        netcdf.putAtt(ncid, varid_range, 'sdn_uom_name', 'Kilometres');
        netcdf.putAtt(ncid, varid_range, 'sdn_uom_urn', 'SDN:P06::ULKM');
        netcdf.putAtt(ncid, varid_range, 'ancillary_variables', 'POSITION_SEADATANET_QC');
        
        % Depth (arbitrary 'z' dimension)
        varid_depth = netcdf.defVar(ncid, 'DEPH', 'float', dimid_depth);
        netcdf.putAtt(ncid, varid_depth, 'long_name', 'Depth of Measurement');
        netcdf.putAtt(ncid, varid_depth, 'standard_name', 'depth');
        netcdf.putAtt(ncid, varid_depth, 'units', 'm');
        netcdf.putAtt(ncid, varid_depth, 'axis', 'Z');
        netcdf.putAtt(ncid, varid_depth, 'positive', 'down');
        netcdf.putAtt(ncid, varid_depth, 'reference', 'sea_level');
        netcdf.putAtt(ncid, varid_depth, 'sdn_parameter_name', 'Depth below surface of the water body');
        netcdf.putAtt(ncid, varid_depth, 'sdn_parameter_urn', 'SDN:P01::ADEPZZ01');
        netcdf.putAtt(ncid, varid_depth, 'sdn_uom_name', 'Metres');
        netcdf.putAtt(ncid, varid_depth, 'sdn_uom_urn', 'SDN:P06::ULAA');
        netcdf.putAtt(ncid, varid_depth, 'ancillary_variables', 'DEPTH_SEADATANET_QC');
        
        
        %% Add auxillary coordinate variables to provide mapping from range and bearing to lat, lon.
        
        % A minimum of 4 significant figures to the right of the decimal
        % place is needed to keep resolution below 10's of meters. Using
        % float data type for lat yields at least 5 significant digits to
        % the right of the decimal giving ~1/2m resolution in latitude.
        % Nine significant figures (at least 7 to the right of the
        % decimal) could be achieved using int data type but need to
        % introduce a scale factor.
        
        % Latitude
        varid_lat = netcdf.defVar( ncid, 'LATITUDE', 'float', [dimid_range dimid_bearing] );
        netcdf.defVarDeflate(ncid, varid_lat, true, true, 6);
        netcdf.putAtt( ncid, varid_lat, 'standard_name', 'latitude' );
        netcdf.putAtt( ncid, varid_lat, 'long_name', 'Latitude' );
        netcdf.putAtt( ncid, varid_lat, 'units', 'degrees_north' );
        netcdf.putAtt(ncid, varid_lat, 'valid_range', single( [-90 90]));
        netcdf.putAtt( ncid, varid_lat, 'coordinates', 'HEAD RNGE' );
        netcdf.putAtt(ncid, varid_lat, 'FillValue', netcdf.getConstant('NC_FILL_FLOAT'));
        %     netcdf.putAtt(ncid, varid_lat, 'axis', 'Y');
        netcdf.putAtt(ncid, varid_lat, 'sdn_parameter_name', 'Latitude north');
        netcdf.putAtt(ncid, varid_lat, 'sdn_parameter_urn', 'SDN:P01::ALATZZ01');
        netcdf.putAtt(ncid, varid_lat, 'sdn_uom_name', 'Degrees north');
        netcdf.putAtt(ncid, varid_lat, 'sdn_uom_urn', 'SDN:P06::DEGN');
        netcdf.putAtt(ncid, varid_lat, 'ancillary_variables', 'POSITION_SEADATANET_QC');
        
        % Longitude
        varid_lon = netcdf.defVar( ncid, 'LONGITUDE', 'float', [dimid_range dimid_bearing] );
        netcdf.defVarDeflate(ncid, varid_lon, true, true, 6);
        netcdf.putAtt( ncid, varid_lon, 'standard_name', 'longitude' );
        netcdf.putAtt( ncid, varid_lon, 'long_name', 'Longitude' );
        netcdf.putAtt( ncid, varid_lon, 'units', 'degrees_east' );
        netcdf.putAtt(ncid, varid_lon, 'valid_range', single( [-180 180]));
        netcdf.putAtt( ncid, varid_lon, 'coordinates', 'HEAD RNGE' );
        netcdf.putAtt(ncid, varid_lon, 'FillValue', netcdf.getConstant('NC_FILL_FLOAT'));
        netcdf.putAtt(ncid, varid_lon, 'sdn_parameter_name', 'Longitude east');
        netcdf.putAtt(ncid, varid_lon, 'sdn_parameter_urn', 'SDN:P01::ALONZZ01');
        netcdf.putAtt(ncid, varid_lon, 'sdn_uom_name', 'Degrees east');
        netcdf.putAtt(ncid, varid_lon, 'sdn_uom_urn', 'SDN:P06::DEGE');
        netcdf.putAtt(ncid, varid_lon, 'ancillary_variables', 'POSITION_SEADATANET_QC');
        
        % crs
        varid_crs = netcdf.defVar( ncid, 'crs', 'short', dimid_t);
        netcdf.defVarDeflate(ncid, varid_crs, true, true, 6);
        netcdf.putAtt( ncid, varid_crs, 'grid_mapping_name', 'latitude_longitude' );
        netcdf.putAtt( ncid, varid_crs, 'epsg_code', 'EPSG:4326' );
        netcdf.putAtt( ncid, varid_crs, 'semi_major_axis', single(6378137.0) );
        netcdf.putAtt(ncid, varid_crs, 'inverse_flattening', single(298.257223563));
        
        %%
        
        %% Add SDN namespace variables
        
        % To enforce homogeneity in the codes and interoperability with SDC, the site_code has to be set equal 
        % to the EDIOS Series id of the HFR network and platform_code must include the EDIOS Platform id of the 
        % HFR site, i.e.:
        %           SDN_CRUISE=site_code=EDIOS-Series-id
        %           SDN_STATION=platform_code=EDIOS-Series-id_Total (for total current data files)
        %           SDN_STATION=platform_code=EDIOS-Series-id_ EDIOS-Platform-id (for radial current data files)
        
        % SDN_CRUISE
        varid_sdncruise = netcdf.defVar( ncid, 'SDN_CRUISE', 'char', [dimid_string_site_code dimid_t]);
        netcdf.defVarDeflate(ncid, varid_sdncruise, true, true, 6);
        netcdf.putAtt( ncid, varid_sdncruise, 'long_name', 'Grid grouping label');
        
        % SDN_STATION
        varid_sdnstation = netcdf.defVar( ncid, 'SDN_STATION', 'char', [dimid_string_platform_code dimid_t]);
        netcdf.defVarDeflate(ncid, varid_sdnstation, true, true, 6);
        netcdf.putAtt( ncid, varid_sdnstation, 'long_name', 'Grid label');
        
        % SDN_LOCAL_CDI_ID
        varid_sdnlocalcdiid = netcdf.defVar( ncid, 'SDN_LOCAL_CDI_ID', 'char', [dimid_string_sdn_local_cdi_id dimid_t]);
        netcdf.defVarDeflate(ncid, varid_sdnlocalcdiid, true, true, 6);
        netcdf.putAtt( ncid, varid_sdnlocalcdiid, 'long_name', 'SeaDataCloud CDI identifier');
        netcdf.putAtt( ncid, varid_sdnlocalcdiid, 'cf_role', 'grid_id');
        
        % SDN_EDMO_CODE
        varid_sdnedmocode = netcdf.defVar( ncid, 'SDN_EDMO_CODE', 'short', dimid_t);
        netcdf.defVarDeflate(ncid, varid_sdnedmocode, true, true, 6);
        netcdf.putAtt( ncid, varid_sdnedmocode, 'long_name', 'European Directory of Marine Organisations code for the CDI partner');
        
        % SDN_REFERENCES
        varid_sdnreferences = netcdf.defVar( ncid, 'SDN_REFERENCES', 'char', [dimid_string_sdn_references dimid_t]);
        netcdf.defVarDeflate(ncid, varid_sdnreferences, true, true, 6);
        netcdf.putAtt( ncid, varid_sdnreferences, 'long_name', 'Usage metadata reference');
        
        % SDN_XLINK
        varid_sdnxlink = netcdf.defVar(ncid, 'SDN_XLINK', 'char', [dimid_string_sdn_xlink dimid_refmax dimid_t]);
        netcdf.defVarDeflate(ncid, varid_sdnxlink, true, true, 6);
        netcdf.putAtt( ncid, varid_sdnxlink, 'long_name', 'External resource linkages');
        
        %%
               
        %% Add data variables
        
        % radial_sea_water_velocity_away_from_instrument:
        % A velocity is a vector quantity. Radial velocity away from instrument
        % means the component of the velocity along the line of sight of the
        % instrument where positive implies movement away from the instrument (i.e.
        % outward). The "instrument" (examples are radar and lidar) is the device
        % used to make an observation.
        varid_speed = netcdf.defVar(ncid, 'RDVA', 'float', [dimid_range dimid_bearing dimid_depth dimid_t]);
        netcdf.defVarDeflate(ncid, varid_speed, true, true, 6);
        netcdf.putAtt(ncid, varid_speed, 'valid_range', single( [-10 10]));
        netcdf.putAtt(ncid, varid_speed, 'standard_name', 'radial_sea_water_velocity_away_from_instrument');
        netcdf.putAtt(ncid, varid_speed, 'units', 'm s-1');
        netcdf.putAtt(ncid, varid_speed, 'standard_name', 'radial_sea_water_velocity_away_from_instrument');
        netcdf.putAtt(ncid, varid_speed, 'long_name', 'Radial Sea Water Velocity Away From Instrument');
        netcdf.putAtt(ncid, varid_speed, 'FillValue', netcdf.getConstant('NC_FILL_FLOAT'));
        netcdf.putAtt(ncid, varid_speed, 'scale_factor', single(1));
        netcdf.putAtt(ncid, varid_speed, 'add_offset', single(0));
        netcdf.putAtt(ncid, varid_speed, 'sdn_parameter_name', 'Current speed (Eulerian) in the water body by directional range-gated radar');
        netcdf.putAtt(ncid, varid_speed, 'sdn_parameter_urn', 'SDN:P01::LCSAWVRD');
        netcdf.putAtt(ncid, varid_speed, 'sdn_uom_name', 'Metres per second');
        netcdf.putAtt(ncid, varid_speed, 'sdn_uom_urn', 'SDN:P06::UVAA');
        netcdf.putAtt( ncid, varid_speed, 'coordinates', 'TIME DEPH HEAD RNGE' );
        netcdf.putAtt(ncid, varid_speed, 'ancillary_variables', 'QCflag, OWTR_QC, MDFL_QC, CSPD_QC, RDCT_QC');
        
        % (radial current) direction
        % direction_of_radial_vector_away_from_instrument:
        % The direction_of_radial_vector_away_from_instrument is the direction in
        % which the instrument itself is pointing. The direction is measured
        % positive clockwise from due north. The "instrument" (examples are radar
        % and lidar) is the device used to make an observation. "direction_of_X"
        % means direction of a vector, a bearing.
        varid_direction = netcdf.defVar(ncid, 'DRVA', 'short', [dimid_range dimid_bearing dimid_depth dimid_t]);
        netcdf.defVarDeflate(ncid, varid_direction, true, true, 6);
        netcdf.putAtt(ncid, varid_direction, 'valid_range', int16( [0 360]));
        netcdf.putAtt(ncid, varid_direction, 'standard_name', 'direction_of_radial_vector_away_from_instrument');
        netcdf.putAtt(ncid, varid_direction, 'long_name', 'Direction Of Radial Vector Away From Instrument');
        netcdf.putAtt(ncid, varid_direction, 'FillValue', netcdf.getConstant('NC_FILL_SHORT'));
        netcdf.putAtt(ncid, varid_direction, 'add_offset', single(0));
        netcdf.putAtt(ncid, varid_direction, 'units', 'degrees_true');
        netcdf.putAtt(ncid, varid_direction, 'scale_factor', single(0.1));
        netcdf.putAtt(ncid, varid_direction, 'sdn_parameter_name', 'Current direction (Eulerian) in the water body by directional range-gated radar');
        netcdf.putAtt(ncid, varid_direction, 'sdn_parameter_urn', 'SDN:P01::LCDAWVRD');
        netcdf.putAtt(ncid, varid_direction, 'sdn_uom_name', 'Degrees True');
        netcdf.putAtt(ncid, varid_direction, 'sdn_uom_urn', 'SDN:P06::UABB');
        netcdf.putAtt( ncid, varid_direction, 'coordinates', 'TIME DEPH HEAD RNGE' );        
        netcdf.putAtt(ncid, varid_direction, 'ancillary_variables', 'QCflag, OWTR_QC, MDFL_QC, AVRB_QC, RDCT_QC');
        
        % u
        varid_u = netcdf.defVar(ncid, 'EWCT', 'float', [dimid_range dimid_bearing dimid_depth dimid_t]);
        netcdf.defVarDeflate(ncid, varid_u, true, true, 6);
        netcdf.putAtt(ncid, varid_u, 'valid_range', single( [-10 10]));
        netcdf.putAtt(ncid, varid_u, 'standard_name', 'surface_eastward_sea_water_velocity');
        netcdf.putAtt(ncid, varid_u, 'long_name', 'Surface Eastward Sea Water Velocity');
        netcdf.putAtt(ncid, varid_u, 'FillValue', netcdf.getConstant('NC_FILL_FLOAT'));
        netcdf.putAtt(ncid, varid_u, 'scale_factor', single(1));
        netcdf.putAtt(ncid, varid_u, 'add_offset', single(0));
        netcdf.putAtt(ncid, varid_u, 'units', 'm s-1');
        netcdf.putAtt(ncid, varid_u, 'sdn_parameter_name', 'Eastward current velocity in the water body');
        netcdf.putAtt(ncid, varid_u, 'sdn_parameter_urn', 'SDN:P01::LCEWZZ01');
        netcdf.putAtt(ncid, varid_u, 'sdn_uom_name', 'Metres per second');
        netcdf.putAtt(ncid, varid_u, 'sdn_uom_urn', 'SDN:P06::UVAA');
        netcdf.putAtt(ncid, varid_u, 'coordinates', 'TIME DEPH HEAD RNGE' );
        netcdf.putAtt(ncid, varid_u, 'ancillary_variables', 'QCflag, OWTR_QC, MDFL_QC, CSPD_QC, VART_QC, AVRB_QC, RDCT_QC');
        
        % v
        varid_v = netcdf.defVar(ncid, 'NSCT', 'float', [dimid_range dimid_bearing dimid_depth dimid_t]);
        netcdf.defVarDeflate(ncid, varid_v, true, true, 6);
        netcdf.putAtt(ncid, varid_v, 'valid_range', single( [-10 10]));
        netcdf.putAtt(ncid, varid_v, 'standard_name', 'surface_northward_sea_water_velocity');
        netcdf.putAtt(ncid, varid_v, 'long_name', 'Surface Northward Sea Water Velocity');
        netcdf.putAtt(ncid, varid_v, 'FillValue', netcdf.getConstant('NC_FILL_FLOAT'));
        netcdf.putAtt(ncid, varid_v, 'scale_factor', single(1));
        netcdf.putAtt(ncid, varid_v, 'add_offset', single(0));
        netcdf.putAtt(ncid, varid_v, 'units', 'm s-1');
        netcdf.putAtt(ncid, varid_v, 'sdn_parameter_name', 'Northward current velocity in the water body');
        netcdf.putAtt(ncid, varid_v, 'sdn_parameter_urn', 'SDN:P01::LCNSZZ01');
        netcdf.putAtt(ncid, varid_v, 'sdn_uom_name', 'Metres per second');
        netcdf.putAtt(ncid, varid_v, 'sdn_uom_urn', 'SDN:P06::UVAA');
        netcdf.putAtt(ncid, varid_v, 'coordinates', 'TIME DEPH HEAD RNGE' );
        netcdf.putAtt(ncid, varid_v, 'ancillary_variables', 'QCflag, OWTR_QC, MDFL_QC, CSPD_QC, VART_QC, AVRB_QC, RDCT_QC');
        
        % Spatial Quality
        varid_espc = netcdf.defVar(ncid, 'ESPC', 'float', [dimid_range dimid_bearing dimid_depth dimid_t]);
        netcdf.defVarDeflate(ncid, varid_espc, true, true, 6);
        netcdf.putAtt(ncid, varid_espc, 'valid_range', single( [-1e3 1e3] ));
        netcdf.putAtt(ncid, varid_espc, 'long_name', 'Radial Standard Deviation of Current Velocity over the Scatter Patch');
        netcdf.putAtt(ncid, varid_espc, 'units', 'm s-1');
        netcdf.putAtt(ncid, varid_espc, 'coordinates', 'TIME DEPH HEAD RNGE' );
        netcdf.putAtt(ncid, varid_espc, 'FillValue', netcdf.getConstant('NC_FILL_FLOAT'));
        netcdf.putAtt(ncid, varid_espc, 'scale_factor', single(1));
        netcdf.putAtt(ncid, varid_espc, 'add_offset', single(0));
        netcdf.putAtt(ncid, varid_espc, 'sdn_parameter_name', '');
        netcdf.putAtt(ncid, varid_espc, 'sdn_parameter_urn', '');
        netcdf.putAtt(ncid, varid_espc, 'sdn_uom_name', 'Metres per second');
        netcdf.putAtt(ncid, varid_espc, 'sdn_uom_urn', 'SDN:P06::UVAA');
        netcdf.putAtt(ncid, varid_espc, 'ancillary_variables', 'QCflag, VART_QC');
        
        % Temporal Quality
        varid_etmp = netcdf.defVar(ncid, 'ETMP', 'float', [dimid_range dimid_bearing dimid_depth dimid_t]);
        netcdf.defVarDeflate(ncid, varid_etmp, true, true, 6);
        netcdf.putAtt(ncid, varid_etmp, 'valid_range', single( [-1e3 1e3] ));
        netcdf.putAtt(ncid, varid_etmp, 'long_name', 'Radial Standard Deviation of Current Velocity over Coverage Period');
        netcdf.putAtt(ncid, varid_etmp, 'units', 'm s-1');
        netcdf.putAtt(ncid, varid_etmp, 'coordinates', 'TIME DEPH HEAD RNGE' );
        netcdf.putAtt(ncid, varid_etmp, 'FillValue', netcdf.getConstant('NC_FILL_FLOAT'));
        netcdf.putAtt(ncid, varid_etmp, 'scale_factor', single(1));
        netcdf.putAtt(ncid, varid_etmp, 'add_offset', single(0));
        netcdf.putAtt(ncid, varid_etmp, 'sdn_parameter_name', '');
        netcdf.putAtt(ncid, varid_etmp, 'sdn_parameter_urn', '');
        netcdf.putAtt(ncid, varid_etmp, 'sdn_uom_name', 'Metres per second');
        netcdf.putAtt(ncid, varid_etmp, 'sdn_uom_urn', 'SDN:P06::UVAA');
        netcdf.putAtt(ncid, varid_etmp, 'ancillary_variables', 'QCflag, VART_QC');
        
        % Velocity Maximum
        varid_maxv = netcdf.defVar(ncid, 'MAXV', 'float', [dimid_range dimid_bearing dimid_depth dimid_t]);
        netcdf.defVarDeflate(ncid, varid_maxv, true, true, 6);
        netcdf.putAtt(ncid, varid_maxv, 'long_name', 'Radial Sea Water Velocity Away From Instrument Maximum');
        netcdf.putAtt(ncid, varid_maxv, 'valid_range', single( [-10 10] ));
        netcdf.putAtt(ncid, varid_maxv, 'FillValue', netcdf.getConstant('NC_FILL_FLOAT'));
        netcdf.putAtt(ncid, varid_maxv, 'scale_factor', single(1));
        netcdf.putAtt(ncid, varid_maxv, 'add_offset', single(0));
        netcdf.putAtt(ncid, varid_maxv, 'units', 'm s-1' );
        netcdf.putAtt(ncid, varid_maxv, 'sdn_parameter_name', 'Current speed (Eulerian) in the water body by directional range-gated radar');
        netcdf.putAtt(ncid, varid_maxv, 'sdn_parameter_urn', 'SDN:P01::LCSAWVRD');
        netcdf.putAtt(ncid, varid_maxv, 'sdn_uom_name', 'Metres per second');
        netcdf.putAtt(ncid, varid_maxv, 'sdn_uom_urn', 'SDN:P06::UVAA');
        netcdf.putAtt(ncid, varid_maxv, 'coordinates', 'TIME DEPH HEAD RNGE' );
        netcdf.putAtt(ncid, varid_maxv, 'ancillary_variables', 'QCflag, MDFL_QC, CSPD_QC, VART_QC');
        
        % Velocity Minimum
        varid_minv = netcdf.defVar(ncid, 'MINV', 'float', [dimid_range dimid_bearing dimid_depth dimid_t]);
        netcdf.defVarDeflate(ncid, varid_minv, true, true, 6);
        netcdf.putAtt( ncid, varid_minv, 'long_name','Radial Sea Water Velocity Away From Instrument Minimum');
        netcdf.putAtt(ncid, varid_minv, 'valid_range', single( [-10 10] ));
        netcdf.putAtt(ncid, varid_minv, 'FillValue', netcdf.getConstant('NC_FILL_FLOAT'));
        netcdf.putAtt(ncid, varid_minv, 'scale_factor', single(1));
        netcdf.putAtt(ncid, varid_minv, 'add_offset', single(0));
        netcdf.putAtt(ncid, varid_minv, 'units', 'm s-1');
        netcdf.putAtt(ncid, varid_minv, 'sdn_parameter_name', 'Current speed (Eulerian) in the water body by directional range-gated radar');
        netcdf.putAtt(ncid, varid_minv, 'sdn_parameter_urn', 'SDN:P01::LCSAWVRD');
        netcdf.putAtt(ncid, varid_minv, 'sdn_uom_name', 'Metres per second');
        netcdf.putAtt(ncid, varid_minv, 'sdn_uom_urn', 'SDN:P06::UVAA');
        netcdf.putAtt(ncid, varid_minv, 'coordinates', 'TIME DEPH HEAD RNGE' );
        netcdf.putAtt(ncid, varid_minv, 'ancillary_variables', 'QCflag, MDFL_QC, CSPD_QC, VART_QC');
        
        % Spatial Count
        varid_ersc = netcdf.defVar(ncid, 'ERSC', 'byte', [dimid_range dimid_bearing dimid_depth dimid_t]);
        netcdf.defVarDeflate(ncid, varid_ersc, true, true, 6);
        netcdf.putAtt(ncid, varid_ersc, 'long_name', 'Radial Sea Water Velocity Spatial Quality Count' );
        netcdf.putAtt(ncid, varid_ersc, 'valid_range', single( [0 127] ));
        netcdf.putAtt(ncid, varid_ersc, 'FillValue', netcdf.getConstant('NC_FILL_BYTE'));
        netcdf.putAtt(ncid, varid_ersc, 'scale_factor', single(1));
        netcdf.putAtt(ncid, varid_ersc, 'add_offset', single(0));
        netcdf.putAtt(ncid, varid_ersc, 'units', '1');
        netcdf.putAtt(ncid, varid_ersc, 'sdn_parameter_name', '');
        netcdf.putAtt(ncid, varid_ersc, 'sdn_parameter_urn', '');
        netcdf.putAtt(ncid, varid_ersc, 'sdn_uom_name', 'Dimensionless');
        netcdf.putAtt(ncid, varid_ersc, 'sdn_uom_urn', 'SDN:P06::UUUU');
        netcdf.putAtt(ncid, varid_ersc, 'coordinates', 'TIME DEPH HEAD RNGE' );
        netcdf.putAtt(ncid, varid_ersc, 'ancillary_variables', 'QCflag');
        
        % Temporal Count
        varid_ertc = netcdf.defVar(ncid, 'ERTC', 'byte', [dimid_range dimid_bearing dimid_depth dimid_t]);
        netcdf.defVarDeflate(ncid, varid_ertc, true, true, 6);
        netcdf.putAtt(ncid, varid_ertc, 'long_name', 'Radial Sea Water Velocity Temporal Quality Count');
        netcdf.putAtt(ncid, varid_ertc, 'valid_range', single( [0 127] ));
        netcdf.putAtt(ncid, varid_ertc, 'FillValue', netcdf.getConstant('NC_FILL_BYTE'));
        netcdf.putAtt(ncid, varid_ertc, 'scale_factor', single(1));
        netcdf.putAtt(ncid, varid_ertc, 'add_offset', single(0));
        netcdf.putAtt(ncid, varid_ertc, 'units', '1');
        netcdf.putAtt(ncid, varid_ertc, 'sdn_parameter_name', '');
        netcdf.putAtt(ncid, varid_ertc, 'sdn_parameter_urn', '');
        netcdf.putAtt(ncid, varid_ertc, 'sdn_uom_name', 'Dimensionless');
        netcdf.putAtt(ncid, varid_ertc, 'sdn_uom_urn', 'SDN:P06::UUUU');
        netcdf.putAtt(ncid, varid_ertc, 'coordinates', 'TIME DEPH HEAD RNGE' );
        netcdf.putAtt(ncid, varid_ertc, 'ancillary_variables', 'QCflag');
        
        % X-Distance
        varid_xdst = netcdf.defVar(ncid, 'XDST', 'float', [dimid_range dimid_bearing]);
        netcdf.defVarDeflate(ncid, varid_xdst, true, true, 6);
        netcdf.putAtt(ncid, varid_xdst, 'long_name', 'Eastward Distance From Instrument');
        netcdf.putAtt(ncid, varid_xdst, 'valid_range', single( [0 1e3] ));
        netcdf.putAtt(ncid, varid_xdst, 'FillValue', netcdf.getConstant('NC_FILL_FLOAT'));
        netcdf.putAtt(ncid, varid_xdst, 'scale_factor', single(1));
        netcdf.putAtt(ncid, varid_xdst, 'add_offset', single(0));
        netcdf.putAtt(ncid, varid_xdst, 'units', 'km');
        netcdf.putAtt(ncid, varid_xdst, 'sdn_parameter_name', '');
        netcdf.putAtt(ncid, varid_xdst, 'sdn_parameter_urn', '');
        netcdf.putAtt(ncid, varid_xdst, 'sdn_uom_name', 'Kilometres');
        netcdf.putAtt(ncid, varid_xdst, 'sdn_uom_urn', 'SDN:P06::ULKM');
        netcdf.putAtt(ncid, varid_xdst, 'coordinates', 'HEAD RNGE' );
        netcdf.putAtt(ncid, varid_xdst, 'ancillary_variables', 'QCflag, OWTR_QC, MDFL_QC, CSPD_QC, VART_QC');
        
        % Y-Distance
        varid_ydst = netcdf.defVar(ncid, 'YDST', 'float', [dimid_range dimid_bearing]);
        netcdf.defVarDeflate(ncid, varid_ydst, true, true, 6);
        netcdf.putAtt(ncid, varid_ydst, 'long_name', 'Northward Distance From Instrument');
        netcdf.putAtt(ncid, varid_ydst, 'valid_range', single( [0 1e3] ));
        netcdf.putAtt(ncid, varid_ydst, 'FillValue', netcdf.getConstant('NC_FILL_FLOAT'));
        netcdf.putAtt(ncid, varid_ydst, 'scale_factor', single(1));
        netcdf.putAtt(ncid, varid_ydst, 'add_offset', single(0));
        netcdf.putAtt(ncid, varid_ydst, 'units', 'km');
        netcdf.putAtt(ncid, varid_ydst, 'sdn_parameter_name', '');
        netcdf.putAtt(ncid, varid_ydst, 'sdn_parameter_urn', '');
        netcdf.putAtt(ncid, varid_ydst, 'sdn_uom_name', 'Kilometres');
        netcdf.putAtt(ncid, varid_ydst, 'sdn_uom_urn', 'SDN:P06::ULKM');
        netcdf.putAtt(ncid, varid_ydst, 'coordinates', 'HEAD RNGE' );
        netcdf.putAtt(ncid, varid_ydst, 'ancillary_variables', 'QCflag, OWTR_QC, MDFL_QC, CSPD_QC, VART_QC');
        
        % Spectra Range Cell
        varid_sprc = netcdf.defVar(ncid, 'SPRC', 'byte', [dimid_range dimid_bearing dimid_depth dimid_t]);
        netcdf.defVarDeflate(ncid, varid_sprc, true, true, 6);
        netcdf.putAtt(ncid, varid_sprc, 'long_name', 'Radial Sea Water Velocity Cross Spectral Range Cell');
        netcdf.putAtt(ncid, varid_sprc, 'valid_range', single( [-1e3 1e3] ));
        netcdf.putAtt(ncid, varid_sprc, 'FillValue', netcdf.getConstant('NC_FILL_BYTE'));
        netcdf.putAtt(ncid, varid_sprc, 'scale_factor', single(1));
        netcdf.putAtt(ncid, varid_sprc, 'add_offset', single(0));
        netcdf.putAtt(ncid, varid_sprc, 'units', '1');
        netcdf.putAtt(ncid, varid_sprc, 'sdn_parameter_name', '');
        netcdf.putAtt(ncid, varid_sprc, 'sdn_parameter_urn', '');
        netcdf.putAtt(ncid, varid_sprc, 'sdn_uom_name', 'Dimensionless');
        netcdf.putAtt(ncid, varid_sprc, 'sdn_uom_urn', 'SDN:P06::UUUU');
        netcdf.putAtt(ncid, varid_sprc, 'coordinates', 'TIME DEPH HEAD RNGE' );
        netcdf.putAtt(ncid, varid_sprc, 'ancillary_variables', 'QCflag, OWTR_QC, MDFL_QC, CSPD_QC, VART_QC');
        
        % Number of receive antennas
        varid_narx = netcdf.defVar(ncid, 'NARX', 'short', dimid_t);
        netcdf.defVarDeflate(ncid, varid_narx, true, true, 6);
        netcdf.putAtt(ncid, varid_narx, 'long_name', 'Number of Receive Antennas');
        netcdf.putAtt(ncid, varid_narx, 'valid_range', int16([0 maxsite] ));
        netcdf.putAtt(ncid, varid_narx, 'FillValue', netcdf.getConstant('NC_FILL_SHORT'));
        netcdf.putAtt(ncid, varid_narx, 'scale_factor', int16(1));
        netcdf.putAtt(ncid, varid_narx, 'add_offset', int16(0));
        netcdf.putAtt(ncid, varid_narx, 'units', char('1'));
        netcdf.putAtt(ncid, varid_narx, 'sdn_parameter_name', '');
        netcdf.putAtt(ncid, varid_narx, 'sdn_parameter_urn', '');
        netcdf.putAtt(ncid, varid_narx, 'sdn_uom_name', 'Dimensionless');
        netcdf.putAtt(ncid, varid_narx, 'sdn_uom_urn', 'SDN:P06::UUUU');
        %         netcdf.putAtt(ncid, varid_sprc, 'coordinates', 'TIME' );
        
        % Number of transmit antennas
        varid_natx = netcdf.defVar(ncid, 'NATX', 'short', dimid_t);
        netcdf.defVarDeflate(ncid, varid_natx, true, true, 6);
        netcdf.putAtt(ncid, varid_natx, 'long_name', 'Number of Transmit Antennas');
        netcdf.putAtt(ncid, varid_natx, 'valid_range', int16([0 maxsite] ));
        netcdf.putAtt(ncid, varid_natx, 'FillValue', netcdf.getConstant('NC_FILL_SHORT'));
        netcdf.putAtt(ncid, varid_natx, 'scale_factor', int16(1));
        netcdf.putAtt(ncid, varid_natx, 'add_offset', int16(0));
        netcdf.putAtt(ncid, varid_natx, 'units', char('1'));
        netcdf.putAtt(ncid, varid_natx, 'sdn_parameter_name', '');
        netcdf.putAtt(ncid, varid_natx, 'sdn_parameter_urn', '');
        netcdf.putAtt(ncid, varid_natx, 'sdn_uom_name', 'Dimensionless');
        netcdf.putAtt(ncid, varid_natx, 'sdn_uom_urn', 'SDN:P06::UUUU');
        %         netcdf.putAtt(ncid, varid_natx, 'coordinates', 'TIME' );
        
        % Receive antenna latitudes
        varid_sltr = netcdf.defVar(ncid, 'SLTR', 'float', [dimid_maxsite dimid_t]);
        netcdf.defVarDeflate(ncid, varid_sltr, true, true, 6);
        netcdf.putAtt(ncid, varid_sltr, 'long_name', 'Receive Antenna Latitudes');
        netcdf.putAtt( ncid, varid_sltr, 'standard_name', 'latitude' );
        netcdf.putAtt(ncid, varid_sltr, 'valid_range', single( [-180 180] ));
        netcdf.putAtt(ncid, varid_sltr, 'FillValue', netcdf.getConstant('NC_FILL_FLOAT'));
        netcdf.putAtt(ncid, varid_sltr, 'scale_factor', single(1));
        netcdf.putAtt(ncid, varid_sltr, 'add_offset', single(0));
        netcdf.putAtt(ncid, varid_sltr, 'units', 'degrees_north');
        netcdf.putAtt(ncid, varid_sltr, 'sdn_parameter_name', 'Latitude north');
        netcdf.putAtt(ncid, varid_sltr, 'sdn_parameter_urn', 'SDN:P01::ALATZZ01');
        netcdf.putAtt(ncid, varid_sltr, 'sdn_uom_name', 'Degrees north');
        netcdf.putAtt(ncid, varid_sltr, 'sdn_uom_urn', 'SDN:P06::DEGN');
        netcdf.putAtt(ncid, varid_sltr, 'coordinates', 'TIME MAXSITE' );
        
        % Receive antenna longitudes
        varid_slnr = netcdf.defVar(ncid, 'SLNR', 'float', [dimid_maxsite dimid_t]);
        netcdf.defVarDeflate(ncid, varid_slnr, true, true, 6);
        netcdf.putAtt(ncid, varid_slnr, 'long_name', 'Receive Antenna Longitudes');
        netcdf.putAtt( ncid, varid_slnr, 'standard_name', 'longitude' );
        netcdf.putAtt(ncid, varid_slnr, 'valid_range', single( [-90 90] ));
        netcdf.putAtt(ncid, varid_slnr, 'FillValue', netcdf.getConstant('NC_FILL_FLOAT'));
        netcdf.putAtt(ncid, varid_slnr, 'scale_factor', single(1));
        netcdf.putAtt(ncid, varid_slnr, 'add_offset', single(0));
        netcdf.putAtt(ncid, varid_slnr, 'units', 'degrees_east');
        netcdf.putAtt(ncid, varid_slnr, 'sdn_parameter_name', 'Longitude east');
        netcdf.putAtt(ncid, varid_slnr, 'sdn_parameter_urn', 'SDN:P01::ALONZZ01');
        netcdf.putAtt(ncid, varid_slnr, 'sdn_uom_name', 'Degrees east');
        netcdf.putAtt(ncid, varid_slnr, 'sdn_uom_urn', 'SDN:P06::DEGE');
        netcdf.putAtt(ncid, varid_slnr, 'coordinates', 'TIME MAXSITE' );
        
        % Transmit antenna latitudes
        varid_sltt = netcdf.defVar(ncid, 'SLTT', 'float', [dimid_maxsite dimid_t]);
        netcdf.defVarDeflate(ncid, varid_sltt, true, true, 6);
        netcdf.putAtt(ncid, varid_sltt, 'long_name', 'Transmit Antenna Latitudes');
        netcdf.putAtt( ncid, varid_sltt, 'standard_name', 'latitude' );
        netcdf.putAtt(ncid, varid_sltt, 'valid_range', single( [-180 180] ));
        netcdf.putAtt(ncid, varid_sltt, 'FillValue', netcdf.getConstant('NC_FILL_FLOAT'));
        netcdf.putAtt(ncid, varid_sltt, 'scale_factor', single(1));
        netcdf.putAtt(ncid, varid_sltt, 'add_offset', single(0));
        netcdf.putAtt(ncid, varid_sltt, 'units', 'degrees_north');
        netcdf.putAtt(ncid, varid_sltt, 'sdn_parameter_name', 'Latitude north');
        netcdf.putAtt(ncid, varid_sltt, 'sdn_parameter_urn', 'SDN:P01::ALATZZ01');
        netcdf.putAtt(ncid, varid_sltt, 'sdn_uom_name', 'Degrees north');
        netcdf.putAtt(ncid, varid_sltt, 'sdn_uom_urn', 'SDN:P06::DEGN');
        netcdf.putAtt(ncid, varid_sltt, 'coordinates', 'TIME MAXSITE' );
        
        % Transmit antenna longitudes
        varid_slnt = netcdf.defVar(ncid, 'SLNT', 'float', [dimid_maxsite dimid_t]);
        netcdf.defVarDeflate(ncid, varid_slnt, true, true, 6);
        netcdf.putAtt(ncid, varid_slnt, 'long_name', 'Transmit Antenna Longitudes');
        netcdf.putAtt( ncid, varid_slnt, 'standard_name', 'longitude' );
        netcdf.putAtt(ncid, varid_slnt, 'valid_range', single( [-90 90] ));
        netcdf.putAtt(ncid, varid_slnt, 'FillValue', netcdf.getConstant('NC_FILL_FLOAT'));
        netcdf.putAtt(ncid, varid_slnt, 'scale_factor', single(1));
        netcdf.putAtt(ncid, varid_slnt, 'add_offset', single(0));
        netcdf.putAtt(ncid, varid_slnt, 'units', 'degrees_east');
        netcdf.putAtt(ncid, varid_slnt, 'sdn_parameter_name', 'Longitude east');
        netcdf.putAtt(ncid, varid_slnt, 'sdn_parameter_urn', 'SDN:P01::ALONZZ01');
        netcdf.putAtt(ncid, varid_slnt, 'sdn_uom_name', 'Degrees east');
        netcdf.putAtt(ncid, varid_slnt, 'sdn_uom_urn', 'SDN:P06::DEGE');
        netcdf.putAtt(ncid, varid_slnt, 'coordinates', 'TIME MAXSITE' );
        
        % Receive antenna codes
        varid_scdr = netcdf.defVar(ncid, 'SCDR', 'char', [dimid_string4 dimid_maxsite dimid_t]);
        netcdf.defVarDeflate(ncid, varid_scdr, true, true, 6);
        netcdf.putAtt(ncid, varid_scdr, 'long_name', 'Receive Antenna Codes');
        netcdf.putAtt(ncid, varid_scdr, 'sdn_parameter_name', '');
        netcdf.putAtt(ncid, varid_scdr, 'sdn_parameter_urn', '');
        netcdf.putAtt(ncid, varid_scdr, 'sdn_uom_name', 'Dimensionless');
        netcdf.putAtt(ncid, varid_scdr, 'sdn_uom_urn', 'SDN:P06::UUUU');
        
        % Transmit antenna codes
        varid_scdt = netcdf.defVar(ncid, 'SCDT', 'char', [dimid_string4 dimid_maxsite dimid_t]);
        netcdf.defVarDeflate(ncid, varid_scdt, true, true, 6);
        netcdf.putAtt(ncid, varid_scdt, 'long_name', 'Transmit Antenna Codes');
        netcdf.putAtt(ncid, varid_scdt, 'sdn_parameter_name', '');
        netcdf.putAtt(ncid, varid_scdt, 'sdn_parameter_urn', '');
        netcdf.putAtt(ncid, varid_scdt, 'sdn_uom_name', 'Dimensionless');
        netcdf.putAtt(ncid, varid_scdt, 'sdn_uom_urn', 'SDN:P06::UUUU');
        
        %%
        
        %% Add QC variables
        
        % Time QC Flag
        varid_tqc = netcdf.defVar(ncid, 'TIME_SEADATANET_QC', 'short', dimid_t);
        netcdf.defVarDeflate(ncid, varid_tqc, true, true, 6);
        netcdf.putAtt(ncid, varid_tqc, 'long_name', 'Time SeaDataNet quality flag');
        netcdf.putAtt(ncid, varid_tqc, 'valid_range', int16( [0 9]));
        netcdf.putAtt(ncid, varid_tqc, 'flag_values', int16( [0 1 2 3 4 7 8 9]));
        netcdf.putAtt(ncid, varid_tqc, 'flag_meanings', 'unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value');
        netcdf.putAtt(ncid, varid_tqc, 'comment', 'OceanSITES quality flagging for temporal coordinate.');
        netcdf.putAtt(ncid, varid_tqc, 'FillValue', netcdf.getConstant('NC_FILL_SHORT'));
        netcdf.putAtt(ncid, varid_tqc, 'scale_factor', int16(1));
        netcdf.putAtt(ncid, varid_tqc, 'add_offset', int16(0));
        netcdf.putAtt(ncid, varid_tqc, 'units', '1');
        
        % Position QC Flag
        varid_posqc = netcdf.defVar(ncid, 'POSITION_SEADATANET_QC', 'short', [dimid_range dimid_bearing dimid_depth dimid_t]);
        netcdf.defVarDeflate(ncid, varid_posqc, true, true, 6);
        netcdf.putAtt(ncid, varid_posqc, 'long_name', 'Position SeaDataNet quality flag');
        netcdf.putAtt(ncid, varid_posqc, 'valid_range', int16( [0 9]));
        netcdf.putAtt(ncid, varid_posqc, 'flag_values', int16( [0 1 2 3 4 7 8 9]));
        netcdf.putAtt(ncid, varid_posqc, 'flag_meanings', 'unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value');
        netcdf.putAtt(ncid, varid_posqc, 'comment', 'OceanSITES quality flagging for position coordinates');
        netcdf.putAtt(ncid, varid_posqc, 'FillValue', netcdf.getConstant('NC_FILL_SHORT'));
        netcdf.putAtt(ncid, varid_posqc, 'scale_factor', int16(1));
        netcdf.putAtt(ncid, varid_posqc, 'add_offset', int16(0));
        netcdf.putAtt(ncid, varid_posqc, 'units', '1');
        
        % Depth QC Flag
        varid_dqc = netcdf.defVar(ncid, 'DEPTH_SEADATANET_QC', 'short', dimid_t);
        netcdf.defVarDeflate(ncid, varid_dqc, true, true, 6);
        netcdf.putAtt(ncid, varid_dqc, 'long_name', 'Depth SeaDataNet quality flag');
        netcdf.putAtt(ncid, varid_dqc, 'valid_range', int16( [0 9]));
        netcdf.putAtt(ncid, varid_dqc, 'flag_values', int16( [0 1 2 3 4 7 8 9]));
        netcdf.putAtt(ncid, varid_dqc, 'flag_meanings', 'unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value');
        netcdf.putAtt(ncid, varid_dqc, 'comment', 'OceanSITES quality flagging for depth coordinate.');
        netcdf.putAtt(ncid, varid_dqc, 'FillValue', netcdf.getConstant('NC_FILL_SHORT'));
        netcdf.putAtt(ncid, varid_dqc, 'scale_factor', int16(1));
        netcdf.putAtt(ncid, varid_dqc, 'add_offset', int16(0));
        netcdf.putAtt(ncid, varid_dqc, 'units', '1');
        
        % Overall QC Flag
        varid_ovqc = netcdf.defVar(ncid, 'QCflag', 'short', [dimid_range dimid_bearing dimid_depth dimid_t]);
        netcdf.defVarDeflate(ncid, varid_ovqc, true, true, 6);
        netcdf.putAtt(ncid, varid_ovqc, 'long_name', 'Overall Quality Flags');
        netcdf.putAtt(ncid, varid_ovqc, 'valid_range', int16( [0 9]));
        netcdf.putAtt(ncid, varid_ovqc, 'flag_values', int16( [0 1 2 3 4 7 8 9]));
        netcdf.putAtt(ncid, varid_ovqc, 'flag_meanings', 'unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value');
        netcdf.putAtt(ncid, varid_ovqc, 'comment', 'OceanSITES quality flagging for all QC tests.');
        netcdf.putAtt(ncid, varid_ovqc, 'FillValue', netcdf.getConstant('NC_FILL_SHORT'));
        netcdf.putAtt(ncid, varid_ovqc, 'scale_factor', int16(1));
        netcdf.putAtt(ncid, varid_ovqc, 'add_offset', int16(0));
        netcdf.putAtt(ncid, varid_ovqc, 'units', '1');
        
        % Over-water QC Flag
        varid_owtr = netcdf.defVar(ncid, 'OWTR_QC', 'short', [dimid_range dimid_bearing dimid_depth dimid_t]);
        netcdf.defVarDeflate(ncid, varid_owtr, true, true, 6);
        netcdf.putAtt(ncid, varid_owtr, 'long_name', 'Over-water Quality Flags');
        netcdf.putAtt(ncid, varid_owtr, 'valid_range', int16( [0 9]));
        netcdf.putAtt(ncid, varid_owtr, 'flag_values', int16( [0 1 2 3 4 7 8 9]));
        netcdf.putAtt(ncid, varid_owtr, 'flag_meanings', 'unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value');
        netcdf.putAtt(ncid, varid_owtr, 'comment', 'OceanSITES quality flagging for Over-water QC test.');
        netcdf.putAtt(ncid, varid_owtr, 'FillValue', netcdf.getConstant('NC_FILL_SHORT'));
        netcdf.putAtt(ncid, varid_owtr, 'scale_factor', int16(1));
        netcdf.putAtt(ncid, varid_owtr, 'add_offset', int16(0));
        netcdf.putAtt(ncid, varid_owtr, 'units', '1');
        
        % Median Filter QC Flag
        varid_mdfl = netcdf.defVar(ncid, 'MDFL_QC', 'short', [dimid_range dimid_bearing dimid_depth dimid_t]);
        netcdf.defVarDeflate(ncid, varid_mdfl, true, true, 6);
        netcdf.putAtt(ncid, varid_mdfl, 'long_name', 'Median Filter Quality Flags');
        netcdf.putAtt(ncid, varid_mdfl, 'valid_range', int16( [0 9]));
        netcdf.putAtt(ncid, varid_mdfl, 'flag_values', int16( [0 1 2 3 4 7 8 9]));
        netcdf.putAtt(ncid, varid_mdfl, 'flag_meanings', 'unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value');
        netcdf.putAtt(ncid, varid_mdfl, 'comment', ['OceanSITES quality flagging for Median Filter QC test. Threshold set to ' num2str(Radial_QC_params.MedFilt(1)) ' km, ' num2str(Radial_QC_params.MedFilt(2)) ' deg, ' num2str(Radial_QC_params.MedFilt(3)) ' m/s, ']);
        netcdf.putAtt(ncid, varid_mdfl, 'FillValue', netcdf.getConstant('NC_FILL_SHORT'));
        netcdf.putAtt(ncid, varid_mdfl, 'scale_factor', int16(1));
        netcdf.putAtt(ncid, varid_mdfl, 'add_offset', int16(0));
        netcdf.putAtt(ncid, varid_mdfl, 'units', '1');
        
        % Variance Threshold QC Flag
        varid_vart = netcdf.defVar(ncid, 'VART_QC', 'short', [dimid_range dimid_bearing dimid_depth dimid_t]);
        netcdf.defVarDeflate(ncid, varid_vart, true, true, 6);
        netcdf.putAtt(ncid, varid_vart, 'long_name', 'Variance Threshold Quality Flags');
        netcdf.putAtt(ncid, varid_vart, 'valid_range', int16( [0 9]));
        netcdf.putAtt(ncid, varid_vart, 'flag_values', int16( [0 1 2 3 4 7 8 9]));
        netcdf.putAtt(ncid, varid_vart, 'flag_meanings', 'unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value');
        netcdf.putAtt(ncid, varid_vart, 'comment', ['OceanSITES quality flagging for Variance Threshold QC test. Test not applicable to Direction Finding systems. The Temporal Derivative test is applied. Threshold set to ' num2str(Radial_QC_params.TempDerThr.threshold) ' m/s.']);
        netcdf.putAtt(ncid, varid_vart, 'FillValue', netcdf.getConstant('NC_FILL_SHORT'));
        netcdf.putAtt(ncid, varid_vart, 'scale_factor', int16(1));
        netcdf.putAtt(ncid, varid_vart, 'add_offset', int16(0));
        netcdf.putAtt(ncid, varid_vart, 'units', '1');
        
        % Velocity Threshold QC Flag
        varid_cspd = netcdf.defVar(ncid, 'CSPD_QC', 'short', [dimid_range dimid_bearing dimid_depth dimid_t]);
        netcdf.defVarDeflate(ncid, varid_cspd, true, true, 6);
        netcdf.putAtt(ncid, varid_cspd, 'long_name', 'Velocity Threshold Quality Flags');
        netcdf.putAtt(ncid, varid_cspd, 'valid_range', int16( [0 9]));
        netcdf.putAtt(ncid, varid_cspd, 'flag_values', int16( [0 1 2 3 4 7 8 9]));
        netcdf.putAtt(ncid, varid_cspd, 'flag_meanings', 'unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value');
        netcdf.putAtt(ncid, varid_cspd, 'comment', ['OceanSITES quality flagging for Velocity Threshold QC test. Threshold set to ' num2str(Radial_QC_params.VelThr) ' m/s.']);
        netcdf.putAtt(ncid, varid_cspd, 'FillValue', netcdf.getConstant('NC_FILL_SHORT'));
        netcdf.putAtt(ncid, varid_cspd, 'scale_factor', int16(1));
        netcdf.putAtt(ncid, varid_cspd, 'add_offset', int16(0));
        netcdf.putAtt(ncid, varid_cspd, 'units', '1');
        
        % Average Radial Bearing QC Flag
        for aRB_idx=1:length(Radial_QC_params.AvgRadBear.site)
            if(strcmp(Radial_QC_params.AvgRadBear.site(aRB_idx).code, Site(1:4)))
                aRB_range = Radial_QC_params.AvgRadBear.site(aRB_idx).range;
            end
        end
        varid_avrb = netcdf.defVar(ncid, 'AVRB_QC', 'short', dimid_t);
        netcdf.defVarDeflate(ncid, varid_avrb, true, true, 6);
        netcdf.putAtt(ncid, varid_avrb, 'long_name', 'Average Radial Bearing Quality Flag');
        netcdf.putAtt(ncid, varid_avrb, 'valid_range', int16( [0 9]));
        netcdf.putAtt(ncid, varid_avrb, 'flag_values', int16( [0 1 2 3 4 7 8 9]));
        netcdf.putAtt(ncid, varid_avrb, 'flag_meanings', 'unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value');
        netcdf.putAtt(ncid, varid_avrb, 'comment', ['OceanSITES quality flagging for Average Radial Bearing QC test. Thresholds set to [' num2str(aRB_range(1)) '-' num2str(aRB_range(2)) '] deg.']);
        netcdf.putAtt(ncid, varid_avrb, 'FillValue', netcdf.getConstant('NC_FILL_SHORT'));
        netcdf.putAtt(ncid, varid_avrb, 'scale_factor', int16(1));
        netcdf.putAtt(ncid, varid_avrb, 'add_offset', int16(0));
        netcdf.putAtt(ncid, varid_avrb, 'units', '1');
        
        % Radial Count QC Flag
        varid_rdct = netcdf.defVar(ncid, 'RDCT_QC', 'short', dimid_t);
        netcdf.defVarDeflate(ncid, varid_rdct, true, true, 6);
        netcdf.putAtt(ncid, varid_rdct, 'long_name', 'Radial Count Quality Flag');
        netcdf.putAtt(ncid, varid_rdct, 'valid_range', int16( [0 9]));
        netcdf.putAtt(ncid, varid_rdct, 'flag_values', int16( [0 1 2 3 4 7 8 9]));
        netcdf.putAtt(ncid, varid_rdct, 'flag_meanings', 'unknown good_data probably_good_data potentially_correctable_bad_data bad_data nominal_value interpolated_value missing_value');
        netcdf.putAtt(ncid, varid_rdct, 'comment', ['OceanSITES quality flagging for Radial Count QC test. Thresholds set to ' num2str(Radial_QC_params.RadCnt) ' vectors.']);
        netcdf.putAtt(ncid, varid_rdct, 'FillValue', netcdf.getConstant('NC_FILL_SHORT'));
        netcdf.putAtt(ncid, varid_rdct, 'scale_factor', int16(1));
        netcdf.putAtt(ncid, varid_rdct, 'add_offset', int16(0));
        netcdf.putAtt(ncid, varid_rdct, 'units', '1');
        
        %%
        
        %% Add global attributes
        
        varid_global = netcdf.getConstant('global');
        
        % MANDATORY ATTRIBUTES
        
        % Discovery and Identification
        netcdf.putAtt(ncid, varid_global, 'site_code', site_code);
        netcdf.putAtt(ncid, varid_global, 'platform_code', platform_code);
        netcdf.putAtt(ncid, varid_global, 'data_mode', 'R');
        netcdf.putAtt(ncid, varid_global, 'DoA_estimation_method', 'Direction Finding');
        netcdf.putAtt(ncid, varid_global, 'calibration_type', 'APM');
        netcdf.putAtt(ncid, varid_global, 'last_calibration_date', lastPatternStr);
        netcdf.putAtt(ncid, varid_global, 'calibration_link', 'carlo.mantovani@cnr.it');
        netcdf.putAtt(ncid, varid_global, 'title', ['Near Real Time Surface Ocean Radial Velocity by HFR_TirLig_' Site]);
        nc_summary = sprintf('%s %s %s', ...
            'The dataset consists of maps of radial velocity of the surface current in the North-Western Tyrrhenian Sea and Ligurian Sea averaged over a time interval of 1 hour around the cardinal hour.', ...
            'HF-RADAR measurements of ocean velocity are radial in', ...
            'direction relative to the radar location and representative ', ...
            'of the upper 0.3 - 2.5 meters of the ocean.');
        netcdf.putAtt(ncid, varid_global, 'summary', nc_summary);
        netcdf.putAtt(ncid, varid_global, 'source', 'coastal structure');
        netcdf.putAtt(ncid, varid_global, 'source_platform_category_code', '17');
        netcdf.putAtt(ncid, varid_global, 'institution', 'CNR-ISMAR: National Research Council - Institute of Marine Sciences, S.S. Lerici');
        netcdf.putAtt(ncid, varid_global, 'institution_edmo_code', num2str(EDMO_code));
        netcdf.putAtt(ncid, varid_global, 'data_assembly_center', 'European HFR Node');
        netcdf.putAtt(ncid, varid_global, 'id', id);
        % Geo-spatial-temporal
        netcdf.putAtt(ncid, varid_global, 'data_type', 'HF radar radial data');
        netcdf.putAtt(ncid, varid_global, 'feature_type', 'surface');
        netcdf.putAtt(ncid, varid_global, 'geospatial_lat_min', num2str(single(min(min(latd)))));
        netcdf.putAtt(ncid, varid_global, 'geospatial_lat_max', num2str(single(max(max(latd)))));
        netcdf.putAtt(ncid, varid_global, 'geospatial_lon_min', num2str(single(min(min(lond)))));
        netcdf.putAtt(ncid, varid_global, 'geospatial_lon_max', num2str(single(max(max(lond)))));
        netcdf.putAtt(ncid, varid_global, 'geospatial_vertical_max', '0.48');
        netcdf.putAtt(ncid, varid_global, 'geospatial_vertical_min', '0');
        netcdf.putAtt(ncid, varid_global, 'time_coverage_start', timeCoverageStart);
        netcdf.putAtt(ncid, varid_global, 'time_coverage_end', timeCoverageEnd);
        % Conventions used
        netcdf.putAtt(ncid, varid_global, 'format_version', 'v2.1');
        netcdf.putAtt(ncid, varid_global, 'Conventions', 'CF-1.6, OceanSITES-Manual-1.2, Copernicus-InSituTAC-SRD-1.4, CopernicusInSituTAC-ParametersList-3.1.0, Unidata, ACDD, INSPIRE');
        % Publication information
        netcdf.putAtt(ncid, varid_global, 'update_interval', 'void');
        netcdf.putAtt(ncid, varid_global, 'citation', citation_str);
        netcdf.putAtt(ncid, varid_global, 'distribution_statement', distribution_str);
        netcdf.putAtt(ncid, varid_global, 'publisher_name', 'Lorenzo Corgnati');
        netcdf.putAtt(ncid, varid_global, 'publisher_url', 'http://radarhf.ismar.cnr.it');
        netcdf.putAtt(ncid, varid_global, 'publisher_email', 'lorenzo.corgnati@sp.ismar.cnr.it');
        netcdf.putAtt(ncid, varid_global, 'license', 'HF radar sea surface current velocity dataset by CNR-ISMAR is licensed under a Creative Commons Attribution 4.0 International License. You should have received a copy of the license along with this work. If not, see http://creativecommons.org/licenses/by/4.0/.');
        netcdf.putAtt(ncid, varid_global, 'acknowledgment', 'ISMAR HF Radar Network has been established within RITMARE, Jerico-Next and IMPACT projects. The network has been designed, implemented and managed through the efforts of ISMAR S.S. Lerici.');
        % Provenance
        netcdf.putAtt(ncid, varid_global, 'date_created', dateCreated);
        netcdf.putAtt(ncid, varid_global, 'history', hist_create);
        netcdf.putAtt(ncid, varid_global, 'date_modified', dateCreated);
        netcdf.putAtt(ncid, varid_global, 'date_update', dateCreated);
        netcdf.putAtt(ncid, varid_global, 'processing_level', '2B');
        netcdf.putAtt(ncid, varid_global, 'contributor_name', 'Vega Forneris, Cristina Tronconi');
        netcdf.putAtt(ncid, varid_global, 'contributor_role', 'THREDDS expert, metadata expert');
        netcdf.putAtt(ncid, varid_global, 'contributor_email', 'vega.forneris@artov.isac.cnr.it, cristina.tronconi@artov.isac.cnr.it');
        
        % RECOMMENDED ATTRIBUTES
        
        % Discovery and Identification
        netcdf.putAtt(ncid, varid_global, 'project', 'RITMARE, Jerico-Next and IMPACT');
        netcdf.putAtt(ncid, varid_global, 'naming_authority', 'it.cnr.ismar');
        netcdf.putAtt(ncid, varid_global, 'keywords', 'OCEAN CURRENTS, SURFACE WATER, RADAR, SCR-HF');
        netcdf.putAtt(ncid, varid_global, 'keywords_vocabulary', 'GCMD Science Keywords');
        nc_comment = sprintf('%s %s %s', ...
            'HF-RADAR measurements of ocean velocity are radial in', ...
            'direction relative to the radar location and representative ', ...
            'of the upper 0.3 - 2.5 meters of the ocean.');
        netcdf.putAtt(ncid, varid_global, 'comment', nc_comment);
        netcdf.putAtt(ncid, varid_global, 'data_language', 'eng');
        netcdf.putAtt(ncid, varid_global, 'data_character_set', 'utf8');
        netcdf.putAtt(ncid, varid_global, 'metadata_language', 'eng');
        netcdf.putAtt(ncid, varid_global, 'metadata_character_set', 'utf8');
        netcdf.putAtt(ncid, varid_global, 'topic_category', 'oceans');
        netcdf.putAtt(ncid, varid_global, 'network', 'ISMAR_HFR_TirLig');
        % Geo-spatial-temporal
        netcdf.putAtt(ncid, varid_global, 'area', 'Mediterranean Sea');
        netcdf.putAtt(ncid, varid_global, 'geospatial_lat_units', 'degrees_north');
        netcdf.putAtt(ncid, varid_global, 'geospatial_lon_units', 'degrees_east');
        netcdf.putAtt(ncid, varid_global, 'geospatial_lat_resolution', '');
        netcdf.putAtt(ncid, varid_global, 'geospatial_lon_resolution', '');
        netcdf.putAtt(ncid, varid_global, 'geospatial_vertical_resolution', '0.48');
        netcdf.putAtt(ncid, varid_global, 'geospatial_vertical_units', 'm');
        netcdf.putAtt(ncid, varid_global, 'geospatial_vertical_positive', 'down');
        netcdf.putAtt(ncid, varid_global, 'time_coverage_duration', 'PT1H');
        netcdf.putAtt(ncid, varid_global, 'time_coverage_resolution', 'PT1H');
        netcdf.putAtt(ncid, varid_global, 'reference_system', 'EPSG:4806');
        netcdf.putAtt(ncid, varid_global, 'cdm_data_type', 'Grid');
        % Conventions used
        netcdf.putAtt(ncid, varid_global, 'netcdf_version', netcdf.inqLibVers);
        netcdf.putAtt(ncid, varid_global, 'netcdf_format', 'netcdf4_classic');
        
        % OTHER ATTRIBUTES
        netcdf.putAtt(ncid, varid_global, 'metadata_contact', 'lorenzo.corgnati@sp.ismar.cnr.it');
        netcdf.putAtt(ncid, varid_global, 'metadata_date_stamp', dateCreated);
        nc_abstract = sprintf('%s %s %s', ...
            'The dataset consists of maps of radial velocity of the surface current in the North-Western Tyrrhenian Sea and Ligurian Sea averaged over a time interval of 1 hour around the cardinal hour.');
        netcdf.putAtt(ncid, varid_global, 'abstract', nc_abstract);
        netcdf.putAtt(ncid, varid_global, 'standard_name_vocabulary', 'NetCDF Climate and Forecast (CF) Metadata Convention Standard Name Table Version 1.6');
        netcdf.putAtt(ncid, varid_global, 'creator_name', 'Lorenzo Corgnati');
        netcdf.putAtt(ncid, varid_global, 'creator_email', 'lorenzo.corgnati@sp.ismar.cnr.it');
        netcdf.putAtt(ncid, varid_global, 'creator_url', 'http://www.ismar.cnr.it/');
        netcdf.putAtt(ncid, varid_global, 'sensor', 'CODAR SeaSonde');
        netcdf.putAtt(ncid, varid_global, 'operational_manager', 'Carlo Mantovani');
        netcdf.putAtt(ncid, varid_global, 'operational_manager_email', 'carlo.mantovani@cnr.it');
        netcdf.putAtt(ncid, varid_global, 'references', 'HFR_Progs Matlab Documentation - Copyright (C) 2006-7 David M. Kaplan; Otero,M. (2013). ENCODING NETCDF RADIAL DATA IN THE HF-RADAR NETWORK. http://cordc.ucsd.edu/projects/mapping/documents/HFRNet_Radial_NetCDF.pdf');
        netcdf.putAtt(ncid, varid_global, 'software_name', 'HFR_Combiner_TirLig');
        netcdf.putAtt(ncid, varid_global, 'software_version', 'v3.1');
        
        % Globals sourced from radial file metadata
        if (exist('CTF', 'var') == 0)
            CTF = '';
        end
        netcdf.putAtt(ncid, varid_global, 'CTF', CTF);
        
        if (exist('FileType', 'var') == 0)
            FileType = '';
        end
        netcdf.putAtt(ncid, varid_global, 'FileType', FileType);
        
        if (exist('LLUVSpec', 'var') == 0)
            LLUVSpec = '';
        end
        netcdf.putAtt(ncid, varid_global, 'LLUVSpec', LLUVSpec);
        
        if (exist('UUID', 'var') == 0)
            UUID = '';
        end
        netcdf.putAtt(ncid, varid_global, 'UUID', UUID);
        
        if (exist('manufacturer', 'var') == 0)
            manufacturer = '';
        end
        netcdf.putAtt(ncid, varid_global, 'manufacturer', manufacturer);
        
        if (exist('Site', 'var') == 0)
            Site = '';
        end
        netcdf.putAtt(ncid, varid_global, 'Site', Site);
        
        if (exist('TimeZone', 'var') == 0)
            TimeZone = '';
        end
        netcdf.putAtt(ncid, varid_global, 'TimeZone', TimeZone);
        
        if (exist('TimeCoverage', 'var') == 0)
            TimeCoverage = '';
        end
        netcdf.putAtt(ncid, varid_global, 'TimeCoverage', TimeCoverage);
        
        if (exist('Origin_str', 'var') == 0)
            Origin_str = '';
        end
        netcdf.putAtt(ncid, varid_global, 'Origin', Origin_str);
        
        if (exist('GreatCircle', 'var') == 0)
            GreatCircle = '';
        end
        netcdf.putAtt(ncid, varid_global, 'GreatCircle', GreatCircle);
        
        if (exist('GeodVersion', 'var') == 0)
            GeodVersion = '';
        end
        netcdf.putAtt(ncid, varid_global, 'GeodVersion', GeodVersion);
        
        if (exist('LLUVTrustData', 'var') == 0)
            LLUVTrustData = '';
        end
        netcdf.putAtt(ncid, varid_global, 'LLUVTrustData', LLUVTrustData);
        
        if (exist('RangeStart', 'var') == 0)
            RangeStart = '';
        end
        netcdf.putAtt(ncid, varid_global, 'RangeStart', RangeStart);
        
        if (exist('RangeEnd', 'var') == 0)
            RangeEnd = '';
        end
        netcdf.putAtt(ncid, varid_global, 'RangeEnd', RangeEnd);
        
        if (exist('RangeResolutionKMeters_str', 'var') == 0)
            RangeResolutionKMeters_str = '';
        end
        netcdf.putAtt(ncid, varid_global, 'RangeResolutionKMeters', RangeResolutionKMeters_str);
        
        if (exist('AntennaBearing', 'var') == 0)
            AntennaBearing = '';
        end
        netcdf.putAtt(ncid, varid_global, 'AntennaBearing', AntennaBearing);
        
        if (exist('ReferenceBearing', 'var') == 0)
            ReferenceBearing = '';
        end
        netcdf.putAtt(ncid, varid_global, 'ReferenceBearing', ReferenceBearing);
        
        if (exist('AngularResolution_str', 'var') == 0)
            AngularResolution_str = '';
        end
        netcdf.putAtt(ncid, varid_global, 'AngularResolution', AngularResolution_str);
        
        if (exist('SpatialResolution', 'var') == 0)
            SpatialResolution = '';
        end
        netcdf.putAtt(ncid, varid_global, 'SpatialResolution', SpatialResolution);
        
        if (exist('PatternType', 'var') == 0)
            PatternType = '';
        end
        netcdf.putAtt(ncid, varid_global, 'PatternType', PatternType);
        
        if (exist('PatternDate', 'var') == 0)
            PatternDate = '';
        end
        netcdf.putAtt(ncid, varid_global, 'PatternDate', PatternDate);
        
        if (exist('PatternResolution', 'var') == 0)
            PatternResolution = '';
        end
        netcdf.putAtt(ncid, varid_global, 'PatternResolution', PatternResolution);
        
        if (exist('TransmitCenterFreqMHz', 'var') == 0)
            TransmitCenterFreqMHz = '';
        end
        netcdf.putAtt(ncid, varid_global, 'TransmitCenterFreqMHz', TransmitCenterFreqMHz);
        
        if (exist('DopplerResolutionHzPerBin', 'var') == 0)
            DopplerResolutionHzPerBin = '';
        end
        netcdf.putAtt(ncid, varid_global, 'DopplerResolutionHzPerBin', DopplerResolutionHzPerBin);
        
        if (exist('FirstOrderMethod', 'var') == 0)
            FirstOrderMethod = '';
        end
        netcdf.putAtt(ncid, varid_global, 'FirstOrderMethod', FirstOrderMethod);
        
        if (exist('BraggSmoothingPoints', 'var') == 0)
            BraggSmoothingPoints = '';
        end
        netcdf.putAtt(ncid, varid_global, 'BraggSmoothingPoints', BraggSmoothingPoints);
        
        if (exist('CurrentVelocityLimit', 'var') == 0)
            CurrentVelocityLimit = '';
        end
        netcdf.putAtt(ncid, varid_global, 'CurrentVelocityLimit', CurrentVelocityLimit);
        
        if (exist('BraggHasSecondOrder', 'var') == 0)
            BraggHasSecondOrder = '';
        end
        netcdf.putAtt(ncid, varid_global, 'BraggHasSecondOrder', BraggHasSecondOrder);
        
        if (exist('RadialBraggPeakDropOff', 'var') == 0)
            RadialBraggPeakDropOff = '';
        end
        netcdf.putAtt(ncid, varid_global, 'RadialBraggPeakDropOff', RadialBraggPeakDropOff);
        
        if (exist('RadialBraggPeakNull', 'var') == 0)
            RadialBraggPeakNull = '';
        end
        netcdf.putAtt(ncid, varid_global, 'RadialBraggPeakNull', RadialBraggPeakNull);
        
        if (exist('RadialBraggNoiseThreshold', 'var') == 0)
            RadialBraggNoiseThreshold = '';
        end
        netcdf.putAtt(ncid, varid_global, 'RadialBraggNoiseThreshold', RadialBraggNoiseThreshold);
        
        if (exist('PatternAmplitudeCorrections', 'var') == 0)
            PatternAmplitudeCorrections = '';
        end
        netcdf.putAtt(ncid, varid_global, 'PatternAmplitudeCorrections', PatternAmplitudeCorrections);
        
        if (exist('PatternPhaseCorrections', 'var') == 0)
            PatternPhaseCorrections = '';
        end
        netcdf.putAtt(ncid, varid_global, 'PatternPhaseCorrections', PatternPhaseCorrections);
        
        if (exist('PatternAmplitudeCalculations', 'var') == 0)
            PatternAmplitudeCalculations = '';
        end
        netcdf.putAtt(ncid, varid_global, 'PatternAmplitudeCalculations', PatternAmplitudeCalculations);
        
        if (exist('PatternPhaseCalculations', 'var') == 0)
            PatternPhaseCalculations = '';
        end
        netcdf.putAtt(ncid, varid_global, 'PatternPhaseCalculations', PatternPhaseCalculations);
        
        if (exist('RadialMusicParameters', 'var') == 0)
            RadialMusicParameters = '';
        end
        netcdf.putAtt(ncid, varid_global, 'RadialMusicParameters', RadialMusicParameters);
        
        if (exist('MergedCount', 'var') == 0)
            MergedCount = '';
        end
        netcdf.putAtt(ncid, varid_global, 'MergedCount', MergedCount);
        
        if (exist('RadialMinimumMergePoints', 'var') == 0)
            RadialMinimumMergePoints = '';
        end
        netcdf.putAtt(ncid, varid_global, 'RadialMinimumMergePoints', RadialMinimumMergePoints);
        
        if (exist('FirstOrderCalc', 'var') == 0)
            FirstOrderCalc = '';
        end
        netcdf.putAtt(ncid, varid_global, 'FirstOrderCalc', FirstOrderCalc);
        
        if (exist('MergeMethod', 'var') == 0)
            MergeMethod = '';
        end
        netcdf.putAtt(ncid, varid_global, 'MergeMethod', MergeMethod);
        
        if (exist('PatternMethod', 'var') == 0)
            PatternMethod = '';
        end
        netcdf.putAtt(ncid, varid_global, 'PatternMethod', PatternMethod);
        
        if (exist('TransmitSweepRateHz', 'var') == 0)
            TransmitSweepRateHz = '';
        end
        netcdf.putAtt(ncid, varid_global, 'TransmitSweepRateHz', TransmitSweepRateHz);
        
        if (exist('TransmitBandwidthKHz', 'var') == 0)
            TransmitBandwidthKHz = '';
        end
        netcdf.putAtt(ncid, varid_global, 'TransmitBandwidthKHz', TransmitBandwidthKHz);
        
        if (exist('SpectraRangeCells', 'var') == 0)
            SpectraRangeCells = '';
        end
        netcdf.putAtt(ncid, varid_global, 'SpectraRangeCells', SpectraRangeCells);
        
        if (exist('SpectraDopplerCells', 'var') == 0)
            SpectraDopplerCells = '';
        end
        netcdf.putAtt(ncid, varid_global, 'SpectraDopplerCells', SpectraDopplerCells);
        
        if (exist('TableType', 'var') == 0)
            TableType = '';
        end
        netcdf.putAtt(ncid, varid_global, 'TableType', TableType);
        
        if (exist('TableColumns', 'var') == 0)
            TableColumns = '';
        end
        netcdf.putAtt(ncid, varid_global, 'TableColumns', TableColumns);
        
        if (exist('TableColumnTypes', 'var') == 0)
            TableColumnTypes = '';
        end
        netcdf.putAtt(ncid, varid_global, 'TableColumnTypes', TableColumnTypes);
        
        if (exist('TableRows', 'var') == 0)
            TableRows = '';
        end
        netcdf.putAtt(ncid, varid_global, 'TableRows', TableRows);
        
        %%
        
        
        %% Exit definition mode
        netcdf.endDef(ncid);
        
        %%
        
        %% Write variable data
        netcdf.putVar(ncid, varid_bearing, bearing_dim);
        netcdf.putVar(ncid, varid_range, range_dim);
        netcdf.putVar(ncid, varid_t, 0, R.time - datenum(1950,1,1));
        netcdf.putVar(ncid, varid_depth, 0);
        netcdf.putVar(ncid, varid_lat, latd);
        netcdf.putVar(ncid, varid_lon, lond);
        netcdf.putVar(ncid, varid_crs, 0);
        netcdf.putVar(ncid, varid_sdncruise, site_code);
        netcdf.putVar(ncid, varid_sdnstation, platform_code); 
        netcdf.putVar(ncid, varid_sdnedmocode, EDMO_code);
        netcdf.putVar(ncid, varid_sdnlocalcdiid, id);
        netcdf.putVar(ncid, varid_sdnreferences, TDS_catalog);
        netcdf.putVar(ncid, varid_sdnxlink, xlink');        
        netcdf.putVar(ncid, varid_speed, radialVelocityMedianFiltered);
        netcdf.putVar(ncid, varid_direction, head);
        netcdf.putVar(ncid, varid_u, velu);
        netcdf.putVar(ncid, varid_v, velv);
        netcdf.putVar(ncid, varid_espc, espc);
        netcdf.putVar(ncid, varid_etmp, etmp);
        netcdf.putVar(ncid, varid_maxv, maxv);
        netcdf.putVar(ncid, varid_minv, minv);
        netcdf.putVar(ncid, varid_ersc, ersc);
        netcdf.putVar(ncid, varid_ertc, ertc);
        netcdf.putVar(ncid, varid_xdst, xdst);
        netcdf.putVar(ncid, varid_ydst, ydst);
        netcdf.putVar(ncid, varid_sprc, sprc);
        netcdf.putVar(ncid, varid_narx, numSites);
        netcdf.putVar(ncid, varid_natx, numSites);
        netcdf.putVar(ncid, varid_sltr, siteLat);
        netcdf.putVar(ncid, varid_slnr, siteLon);
        netcdf.putVar(ncid, varid_sltt, siteLat);
        netcdf.putVar(ncid, varid_slnt, siteLon);
        netcdf.putVar(ncid, varid_scdr, siteCode');
        netcdf.putVar(ncid, varid_scdt, siteCode');
        netcdf.putVar(ncid, varid_tqc, sdnTime_QCflag);
        netcdf.putVar(ncid, varid_posqc, sdnPosition_QCflag);
        netcdf.putVar(ncid, varid_dqc, sdnDepth_QCflag);
        netcdf.putVar(ncid, varid_ovqc, overall_QCflag);
        netcdf.putVar(ncid, varid_owtr, overWater_QCflag);
        netcdf.putVar(ncid, varid_mdfl, medianFilter_QCflag);
        netcdf.putVar(ncid, varid_vart, temporalDerivativeThreshold_QCflag);
        netcdf.putVar(ncid, varid_cspd, velocityThreshold_QCflag);
        netcdf.putVar(ncid, varid_avrb, averageRadialBearing_QC_flag);
        netcdf.putVar(ncid, varid_rdct, radialCount_QC_flag);
        
        %%
        
        %% Close file
        netcdf.close(ncid);
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        R2C_err = 1;
    end
    %%
    
    %% Copy the netCDF file to RadarDisk and THREDDS server
    % Copies the netCDF file to the RadarDisk
    if (R2C_err == 0)
        try
            [cp_status] = copyfile(ncfile_local,ncfile_rd,'f');
        catch err
            display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
            R2C_err = 1;
        end
    end
    
    %%
    
end
return

