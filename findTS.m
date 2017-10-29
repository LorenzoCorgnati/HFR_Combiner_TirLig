%% findTS.m
% This function retieves the position of the input time stamp inside
% the RadStruct structure.

% Author: Lorenzo Corgnati
% Date: November 24, 2013

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [position] = findTS(RadStruct, TS_str)
for sts=1:length(RadStruct)
    if (strcmp(RadStruct(sts).TS,TS_str))
        position = sts;
    end
end
return
