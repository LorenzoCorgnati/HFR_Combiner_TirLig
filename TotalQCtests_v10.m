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
%         GDOPThr: GDOP threshold quality flags
%         dataDens: Data Density quality flags
%         radBal: Radial Balance quality flags
%         velThr: Velocity threshold quality flags

% Author: Lorenzo Corgnati
% Date: November 11, 2016

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [overall, varThr, GDOPThr, dataDens, radBal, velThr] = TotalQCtests_v10(mat_tot, Total_QC_params)

display(['[' datestr(now) '] - - ' 'TotaQCtests_v10.m started.']);

TQC_err = 0;

%% Prepare QC flag variables

% Sets total data on a regular grid.
try
    lonGrid = unique(mat_tot.LonLat(:,1));
    latGrid = unique(mat_tot.LonLat(:,2));
    depth = 0;
catch err
    display(['[' datestr(now) '] - - ' err.message]);
    TQC_err = 1;
end

overall = netcdf.getConstant('NC_FILL_SHORT').*int16(ones(length(lonGrid),length(latGrid),1));
varThr = netcdf.getConstant('NC_FILL_SHORT').*int16(ones(length(lonGrid),length(latGrid),1));
GDOPThr = netcdf.getConstant('NC_FILL_SHORT').*int16(ones(length(lonGrid),length(latGrid),1));
dataDens = netcdf.getConstant('NC_FILL_SHORT').*int16(ones(length(lonGrid),length(latGrid),1,1));
radBal = netcdf.getConstant('NC_FILL_SHORT').*int16(ones(length(lonGrid),length(latGrid),1));
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
    catch err
        display(['[' datestr(now) '] - - ' err.message]);
        TQC_err = 1;
    end
end

%%

%% Populate the overall quality variable

for ii=1:size(overall,1)
    for jj = 1:size(overall,2)
        if(velThr(ii,jj) ~= netcdf.getConstant('NC_FILL_SHORT'))
            if((varThr(ii,jj) == 1) && (velThr(ii,jj) == 1) && (GDOPThr(ii,jj) == 1) && (dataDens(ii,jj) == 1) && (radBal(ii,jj) == 1))
                overall(ii,jj) = 1;
            else
                overall(ii,jj) = 4;
            end
        end
    end
end
%%

return