%% RSsort_TS.m
% This function sorts the RadStruct by TS field.

% Author: Lorenzo Corgnati
% Date: November 24, 2013

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [RSsorted] = RSsort_TS(RadStruct)
% Converts the structure into a cell array
RSfields = fieldnames(RadStruct);
RScell = struct2cell(RadStruct);
sz = size(RScell);

% Converts the cell array to a matrix
RScell = reshape(RScell, sz(1), []);

% Makes each field a column
RScell = RScell';

% Sort by first field "name"
RScell = sortrows(RScell, 1);

% Puts back into original cell array format
RScell = reshape(RScell', sz);

% Converts to Struct
RSsorted = cell2struct(RScell, RSfields, 1);
return
