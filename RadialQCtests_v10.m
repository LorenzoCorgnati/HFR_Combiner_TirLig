%% RadialQCtests_v10.m
% This function performs the QC tests on radial velocity data. The tests are
% the ones defined in INCREASE project compliant to CMEMS needs.
% In particular, the following tests are performed:
%       - Velocity threshold
%       - Variance threshold
%       - Median Filter
%       - Average radial bearing
%       - Vector Over Water (performed in Radial2netCDF_v20.m)

% INPUT:
%         bear: radial velocity bearing variable from radial data
%         lond: longitude coordinates of the grid for radial velocities
%         latd: latitude coordinates of the grid for radial velocities
%         owtr: Vector Over Water quality flags computed in
%               Radial2netCDF_v20.m
%         etmp: temporal quality variable from radial data
%         head: radial velocity heading variable from radial data
%         radVel: radial velocities from radial data
%         Radial_QC_params: structure containing parameters for radial QC tests

% OUTPUT:
%         overall: overall quality flag (good data value is assigned
%                         if and only if all QC tests are passed
%         varThr: Variance threshold quality flags
%         velThr: Velocity threshold quality flags
%         overWater: Vector Over Water quality flags
%         medFilt: Median Filter quality flags
%         avgRadBear: Average Radial Bearing quality flags
%         radVel: radial velocities after the applicatiojn of the median
%         filter


% Author: Lorenzo Corgnati
% Date: January 23, 2017

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [overall, overWater, varThr, velThr, medFilt, avgRadBear, radVelMedianFiltered] = RadialQCtests_v10(bear, lond, latd, owtr, etmp, head, radVel, Radial_QC_params)

display(['[' datestr(now) '] - - ' 'RadialQCtests_v10.m started.']);

RQC_err = 0;

%% Prepare QC flag variables

overall = netcdf.getConstant('NC_FILL_SHORT').*int16(ones(size(owtr,1), size(owtr,2)));
overWater = owtr;
overWater(isnan(overWater)) = netcdf.getConstant('NC_FILL_SHORT');
varThr = netcdf.getConstant('NC_FILL_SHORT').*int16(ones(size(owtr,1), size(owtr,2)));
velThr = netcdf.getConstant('NC_FILL_SHORT').*int16(ones(size(owtr,1), size(owtr,2)));
medFilt = netcdf.getConstant('NC_FILL_SHORT').*int16(ones(size(owtr,1), size(owtr,2)));
avgRadBear = netcdf.getConstant('NC_FILL_SHORT');

%%

%% Prepare QC tests
% Variance Threshold QC test
varVec = etmp.^2;

% Velocity Threshold QC test
maxspd_R = Radial_QC_params.VelThr;

% Average Radial Bearing QC test
avgBear_HEAD = mean(head(~isnan(head)));
avgBear = mean(head(~isnan(bear)));

% Median Filter QC test
radVelMedianFiltered = radVel;


%%

%% Populate QC variables
if (RQC_err == 0)
    try
        % Over Water quality flags
        overWater(overWater==0) = 1;
        overWater(overWater==128) = 4;
    
        % Velocity Threshold quality flags
        velThr((abs(radVel) <= maxspd_R)) = 1;
        velThr((abs(radVel) > maxspd_R)) = 4;
        
        % Variance Threshold quality flags
        varThr((varVec > Radial_QC_params.VarThr)) = 4;
        varThr((varVec <= Radial_QC_params.VarThr)) = 1;
        
        % Average Radial Bearing quality flag
        if ((avgBear >= Radial_QC_params.AvgRadBear(1)) && (avgBear <= Radial_QC_params.AvgRadBear(2)))
            avgRadBear = 1;
        else
            avgRadBear = 4;
        end
        
        % Median Filter quality flags
        for i=1:size(radVel,1)
            for j = 1:size(radVel,2)
                if(~isnan(radVel(i,j)))
                    % Check vector distances
                    jj = j + 1;
                    radiusFlag = 0; % flag saying if the cell is inside (0) or outside (1) the search radius
                    % Scan the grid horizontally toward right
                    while((jj <= size(radVel,2)) && radiusFlag == 0)
                        [cellDist, d2km] = lldistkm([latd(i,j) lond(i,j)], [latd(i,jj) lond(i,jj)]);
                        if(cellDist > Radial_QC_params.MedFilt(1))
                            radiusFlag = 1;
                        end
                        jj = jj + 1;
                    end
                    
                    % Check if the search radius has been found or not
                    if(radiusFlag == 0) % if not, scan the grid vertically downward
                        ii = i + 1;
                        while((ii <= size(radVel,1)) && radiusFlag == 0)
                            [cellDist, d2km] = lldistkm([latd(i,j) lond(i,j)], [latd(ii,j) lond(ii,j)]);
                            if(cellDist > Radial_QC_params.MedFilt(1))
                                radiusFlag = 1;
                            end
                            ii = ii + 1;
                        end
                    elseif(radiusFlag == 1) % build the window around the current cell
                        windowSpan = (jj -1) - j;
                        vel4bear = radVel(max(1,i-windowSpan):min(i+windowSpan,size(radVel,1)),max(1,j-windowSpan):min(j+windowSpan,size(radVel,2)));
                        radiusFlag = 2;
                    end
                    
                    % Check if the search radius has been found or not
                    if(radiusFlag == 0) % if not, scan the grid horizontally toward left
                        jj = j - 1;
                        while((jj > 0) && radiusFlag == 0)
                            [cellDist, d2km] = lldistkm([latd(i,j) lond(i,j)], [latd(i,jj) lond(i,jj)]);
                            if(cellDist > Radial_QC_params.MedFilt(1))
                                radiusFlag = 1;
                            end
                            jj = jj - 1;
                        end
                    elseif(radiusFlag == 1) % build the window around the current cell
                        windowSpan = (ii -1) - i;
                        vel4bear = radVel(max(1,i-windowSpan):min(i+windowSpan,size(radVel,1)),max(1,j-windowSpan):min(j+windowSpan,size(radVel,2)));
                        radiusFlag = 2;
                    end
                    
                    % Check if the search radius has been found or not
                    if(radiusFlag == 0) % if not, scan the grid vertically upward
                        ii = i - 1;
                        while((ii > 0) && radiusFlag == 0)
                            [cellDist, d2km] = lldistkm([latd(i,j) lond(i,j)], [latd(ii,j) lond(ii,j)]);
                            if(cellDist > Radial_QC_params.MedFilt(1))
                                radiusFlag = 1;
                            end
                            ii = ii - 1;
                        end
                    elseif(radiusFlag == 1) % build the window around the current cell
                        windowSpan = j - (jj + 1);
                        vel4bear = radVel(max(1,i-windowSpan):min(i+windowSpan,size(radVel,1)),max(1,j-windowSpan):min(j+windowSpan,size(radVel,2)));
                        radiusFlag = 2;
                    end
                    
                    % Check if the search radius has been found or not
                    if(radiusFlag == 0) % all the grid lies inside the search radius
                        % Use the whole grid as window for median filter
                        vel4bear = radVel;
                    elseif(radiusFlag == 1) % build the window around the current cell
                        windowSpan = i - (ii + 1);
                        vel4bear = radVel(max(1,i-windowSpan):min(i+windowSpan,size(radVel,1)),max(1,j-windowSpan):min(j+windowSpan,size(radVel,2)));
                    end
                    
                    % Check bearing distances within vel4filt indices matrix
                    refBear = bear(i,j);
                    bearDist = bear(max(1,i-windowSpan):min(i+windowSpan,size(radVel,1)),max(1,j-windowSpan):min(j+windowSpan,size(radVel,2)));
                    vel4filt = vel4bear(((abs(bearDist-refBear)) <= Radial_QC_params.MedFilt(2)));
                    
                    % Evaluate quality flag
                    medVal = median(vel4filt(:));
                    if(abs(radVel(i,j) - medVal) <= Radial_QC_params.MedFilt(3))
                        medFilt(i,j) = 1;
                    else
                        radVelMedianFiltered(i,j) = medVal;
                        medFilt(i,j) = 8;
                    end
                end
            end
        end
        
    catch err
        display(['[' datestr(now) '] - - ' err.message]);
        RQC_err = 1;
    end
end

%%

%% Populate the overall quality variable

for ii=1:size(overall,1)
    for jj = 1:size(overall,2)
        if(velThr(ii,jj) ~= netcdf.getConstant('NC_FILL_SHORT'))
            if((varThr(ii,jj) == 1) && (velThr(ii,jj) == 1) && (overWater(ii,jj) == 1) && (medFilt(ii,jj) == 1) && (avgRadBear == 1))
                overall(ii,jj) = 1;
            else
                overall(ii,jj) = 4;
            end
        end
    end
end
%%

return