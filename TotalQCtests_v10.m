%% TotalQCtests_v10.m
% This function performs the QC tests on total velocity data. The tests are
% the ones defined in INCREASE project compliant to CMEMS needs.
% In particular, the following tests are performed:
%       - Velocity threshold
%       - GDOP threshold
%       - Variance threshold
%       - Data Density threshold
%       - Balance of contributing radials from different sites

% INPUT:
%         mat_tot: structure containing total file in Codar format
%         Total_QC_params: structure containing parameters for total QC tests

% OUTPUT:
%         overall: overall quality flag (good data value is assigned
%                         if and only if all QC tests are passed
%         varThr: Variance threshold quality flags
%         tempDer: Temporal Derivative quality flags
%         GDOPThr: GDOP threshold quality flags
%         dataDens: Data Density quality flags
%         velThr: Velocity threshold quality flags

% Author: Lorenzo Corgnati
% Date: November 11, 2016

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [overall, varThr, tempDer, GDOPThr, dataDens, velThr] = TotalQCtests_v10(mat_tot, Total_QC_params)

display(['[' datestr(now) '] - - ' 'TotalQCtests_v10.m started.']);

TQC_err = 0;

%% Prepare QC flag variables

% Sets total data on a regular grid.
try
    lonGrid = unique(mat_tot.LonLat(:,1));
    latGrid = unique(mat_tot.LonLat(:,2));
    depth = 0;
catch err
    display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
    TQC_err = 1;
end

overall = netcdf.getConstant('NC_FILL_SHORT').*int16(ones(length(lonGrid),length(latGrid),1));
varThr = netcdf.getConstant('NC_FILL_SHORT').*int16(ones(length(lonGrid),length(latGrid),1));
tempDer = netcdf.getConstant('NC_FILL_SHORT').*int16(ones(length(lonGrid),length(latGrid),1));
GDOPThr = netcdf.getConstant('NC_FILL_SHORT').*int16(ones(length(lonGrid),length(latGrid),1));
dataDens = netcdf.getConstant('NC_FILL_SHORT').*int16(ones(length(lonGrid),length(latGrid),1,1));
velThr = netcdf.getConstant('NC_FILL_SHORT').*int16(ones(length(lonGrid),length(latGrid),1));

%%

%% Perform QC tests
% Variance Threshold QC test
varVec = sqrt((mat_tot.ErrorEstimates(1,1).Uerr).^2 + (mat_tot.ErrorEstimates(1,1).Verr).^2).*0.0001;
varVec(isnan(mat_tot.U)) = NaN; % exclude grid points where no velocity data is present

% Velocity Threshold and GDOP Threshold QC tests
maxspd_T = Total_QC_params.VelThr*100;
gdop_thr = Total_QC_params.GDOPThr;
[TUV_clean, I] = cleanTotals(mat_tot, maxspd_T, {'GDOPMaxOrthog','TotalErrors',gdop_thr});
I(isnan(mat_tot.U)) = NaN; % exclude grid points where no velocity data is present

% Data Density Threshold
numRads = mat_tot.OtherMatrixVars.makeTotals_TotalsNumRads;
numRads(isnan(mat_tot.U)) = NaN; % exclude grid points where no velocity data is present

% Temporal Derivative QC test
tempDer_Thr = Total_QC_params.TempDerThr.threshold;
% Check if the files of the previous two hours exist
if ((exist(Total_QC_params.TempDerThr.hour2) == 2) && (exist(Total_QC_params.TempDerThr.hour1) == 2))
    tD_go = true;
else
    tD_go = false;
end

%%

%% Populate QC variables
if (TQC_err == 0)
    try
        for i=1:length(mat_tot.LonLat(:,1))
            lonGrid_idx = find(lonGrid==mat_tot.LonLat(i,1));
            latGrid_idx = find(latGrid==mat_tot.LonLat(i,2));
            
            % Velocity Threshold quality flags
            if (I(i) == 1 || I(i) == 3)
                velThr(lonGrid_idx,latGrid_idx,1) = 4;
            elseif (I(i) == 0 || I(i) == 2)
                velThr(lonGrid_idx,latGrid_idx,1) = 1;
            end
            
            % GDOP Threshold quality flags
            if (I(i) == 2 || I(i) == 3)
                GDOPThr(lonGrid_idx,latGrid_idx,1) = 4;
            elseif (I(i) == 0 || I(i) == 1)
                GDOPThr(lonGrid_idx,latGrid_idx,1) = 1;
            end
            
            % Variance Threshold quality flags
            if (not(isnan(varVec(i))))
                if (varVec(i) > Total_QC_params.VarThr)
                    varThr(lonGrid_idx,latGrid_idx,1) = 4;
                else
                    varThr(lonGrid_idx,latGrid_idx,1) = 1;
                end
            end
            
            % Data Density Threshold quality flag
            if (not(isnan(numRads(i))))
                if (numRads(i) < Total_QC_params.DataDensityThr)
                    dataDens(lonGrid_idx,latGrid_idx,1) = 4;
                else
                    dataDens(lonGrid_idx,latGrid_idx,1) = 1;
                end
            end
        end
        
        % Temporal Derivative quality flags
        if (tD_go)
            tempDer1h = tempDer;
            % Extract the U and V velocity fields from the previous two hours files
            totU1h = ncread(Total_QC_params.TempDerThr.hour1,'EWCT');
            totV1h = ncread(Total_QC_params.TempDerThr.hour1,'NSCT');
            totU2h = ncread(Total_QC_params.TempDerThr.hour2,'EWCT');
            totV2h = ncread(Total_QC_params.TempDerThr.hour2,'NSCT');
            % Evaluate total velocities
            totVel = sqrt(((mat_tot.U_grid).^2) + ((mat_tot.V_grid).^2));
            totVel(find(mat_tot.U_grid==netcdf.getConstant('NC_FILL_FLOAT'))) = NaN;
            totVel1h = sqrt(((totU1h).^2) + ((totV1h).^2));
            totVel2h = sqrt(((totU2h).^2) + ((totV2h).^2));
            for tVr=1:size(totVel,1)
                for tVc=1:size(totVel,2)
                    if (~isnan(totVel1h(tVr,tVc)))
                        if ((isnan(totVel(tVr,tVc))) || (isnan(totVel2h(tVr,tVc))))
                            tempDer1h(tVr,tVc) = 1;
                        elseif ((abs(totVel(tVr,tVc) - totVel1h(tVr,tVc)) < tempDer_Thr) && (abs(totVel2h(tVr,tVc) - totVel1h(tVr,tVc)) < tempDer_Thr))
                            tempDer1h(tVr,tVc) = 1;
                        else
                            tempDer1h(tVr,tVc) = 4;
                        end
                    end
                end
            end
            
            % Modify the VART_QC variable of the nc file of the previous
            % hour (both locally and on RadarDisk)
            ncwrite(Total_QC_params.TempDerThr.hour1,'VART_QC',tempDer1h);
            [cp_status] = copyfile(Total_QC_params.TempDerThr.hour1,Total_QC_params.TempDerThr.hour1_RD,'f');
            % ncwrite(Total_QC_params.TempDerThr.hour1_RD,'VART_QC',tempDer1h);
        end
        % Set the QC flag for the current hour to 0 (no QC performed)
        tempDer(find(mat_tot.U_grid ~= netcdf.getConstant('NC_FILL_FLOAT'))) = 0;
    catch err
        display(['[' datestr(now) '] - - ERROR in ' mfilename ' -> ' err.message]);
        TQC_err = 1;
    end
end

%%

%% Populate the overall quality variable

for ii=1:size(overall,1)
    for jj = 1:size(overall,2)
        if(velThr(ii,jj) ~= netcdf.getConstant('NC_FILL_SHORT'))
            if((varThr(ii,jj) == 1) && (tempDer(ii,jj) == 1) && (velThr(ii,jj) == 1) && (GDOPThr(ii,jj) == 1) && (dataDens(ii,jj) == 1))
                overall(ii,jj) = 1;
            else
                overall(ii,jj) = 4;
            end
        end
    end
end
%%

return