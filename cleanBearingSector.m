%% cleanBearingSector.m
% This function cleans out from the input radial data the bearing sectors
% lying outside the real field of view of the radar site.

% Author: Lorenzo Corgnati
% Date: December 5, 2013

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [clRD] = cleanBearingSector(drtRD)
% Selects the real field of view per each radar site
switch drtRD.SiteName
    case 'VIES'
        idx_tbc = find(drtRD.RangeBearHead(:,2) >= 270 | (drtRD.RangeBearHead(:,2) >= 0 & drtRD.RangeBearHead(:,2) <= 120));
    case 'PUGN'
        idx_tbc = find(drtRD.RangeBearHead(:,2) <= 80 | (drtRD.RangeBearHead(:,2) >= 250 & drtRD.RangeBearHead(:,2) <= 360));
    case 'MATT'
        idx_tbc = find(drtRD.RangeBearHead(:,2) <= 30 | (drtRD.RangeBearHead(:,2) >= 225 & drtRD.RangeBearHead(:,2) <= 360));
    case 'MANF'
        idx_tbc = find(drtRD.RangeBearHead(:,2) <= 45 | (drtRD.RangeBearHead(:,2) >= 210 & drtRD.RangeBearHead(:,2) <= 360));
    case 'TINO'
        idx_tbc = find(drtRD.RangeBearHead(:,2) <= 70 | (drtRD.RangeBearHead(:,2) >= 125 & drtRD.RangeBearHead(:,2) <= 360));
    case 'MONT'
        idx_tbc = find((drtRD.RangeBearHead(:,2) >= 225 & drtRD.RangeBearHead(:,2) <= 320));
end

clRD = subsrefRADIAL(drtRD, idx_tbc);
return
