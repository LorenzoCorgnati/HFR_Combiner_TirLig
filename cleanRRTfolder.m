%% cleanRRTfolder.m
% This function removes from the RadarDisk RadialRealtime folder the radial
% files older than a week (radar sites Archivalist rsyncs radial files of
% the last week only).

% Author: Lorenzo Corgnati
% Date: March 24, 2014

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [cl_status] = cleanRRTfolder(RRTfolder)

cl_status = 0;

% Reads the RadarDisk RadialRealtime folder content
folderList = dir(RRTfolder);

% Eliminates the fake names starting with '.'
dot_flag = 1;
while (dot_flag == 1)
    if (folderList(1).name(1) == '.')
        folderList = folderList(2:size(folderList,1));
    else
        dot_flag = 0;
    end
end

time_limit = datenum(now)-10;
for fL_idx=1:length(folderList)
    curr_ts = datenum([str2double(folderList(fL_idx,1).name(11:14)), str2double(folderList(fL_idx,1).name(16:17)), str2double(folderList(fL_idx,1).name(19:20)), 0, 0, 0]);
    if (curr_ts < time_limit)
        try
            delete([RRTfolder, folderList(fL_idx).name]);
        catch err
            cl_status = cl_status + 1;
        end
    end
end

return
