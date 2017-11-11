%% TimeStamp2TS.m
% This function builds the strings in the form YYYY_MM_DD_HHHH related to
% the time stamps of the two previous hours with respect to the input time stamp. 

% INPUT:
%         present: time stamp of the current hour

% OUTPUT:
%         past2: time stamp of 2 hours before the current hour
%         past1: time stamp of 1 hour before the current hour


% Author: Lorenzo Corgnati
% Date: November 11, 2017

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [past2, past1] = TimeStamp2TS(present)

past2_num = present - 2/24;
past1_num = present - 1/24;

past2_vec = datevec(past2_num);
past1_vec = datevec(past1_num);

past2 = [num2str(past2_vec(1)) '_' sprintf('%02d',past2_vec(2)) '_' sprintf('%02d',past2_vec(3)) '_' sprintf('%02d',past2_vec(4)) sprintf('%02d',past2_vec(5))];
past1 = [num2str(past1_vec(1)) '_' sprintf('%02d',past1_vec(2)) '_' sprintf('%02d',past1_vec(3)) '_' sprintf('%02d',past1_vec(4)) sprintf('%02d',past1_vec(5))];

return
