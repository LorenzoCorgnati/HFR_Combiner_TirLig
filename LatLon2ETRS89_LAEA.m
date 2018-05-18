%% LatLon2ETRS89_LAEA.m
% This function converts WGS84 latitude and longitude coordinates into 
% Lambertian Azimuthal Equal Area ETRS89 LAEA Easting and Northing coordinates.
% The grid is based on the ETRS89 Lambert Azimuthal Equal Area (ETRS89-LAEA)
% coordinate reference system with the centre of the projection at the point
% 52? N, 10? E and false easting: x0 = 4321000 m, false northing: y0 = 3210000 m.
% The origin of the grid coincides with the false origin of the ETRS89-LAEA 
% coordinate reference system (x=0, y=0).
%
% INPUT:
%         longitudeArray: array of the longitude values to be converted,
%                           expressed in degrees
%         latitudeArray: array of the latitude values to be converted,
%                           expressed in degrees

% OUTPUT:
%         eastingArray: array of the converted easting values, expressed in
%                           metres
%         northingArray: array of the converted northing values, expressed in
%                           metres


% Author: Lorenzo Corgnati
% Date: May 14, 2018

% E-mail: lorenzo.corgnati@sp.ismar.cnr.it
%%

function [eastingArray, northingArray] = LatLon2ETRS89_LAEA(longitudeArray, latitudeArray)

display(['[' datestr(now) '] - - ' 'LatLon2ETRS89_LAEA.m started.']);

L2E_err = 0;

%% Set the parameters of the ETRS89-LAEA reference system

% Coordinates of the center of projection
lat1 = 52;
lon0 = 10;
lat1r = lat1*pi/180;
lon0r = lon0*pi/180;

% False easting and false northing of the center of projection
FE = 4321000;
FN = 3210000;

% Latitude of the pole
latP = 90;
latPr = latP*pi/180;

% GRS1980 ellipsoid semiaxis and eccentricity
a = 6378137.0;
e = 0.081819191;

%%

%% Convert input lat/lon coordinates from degrees to radians

lonr = longitudeArray.*(pi/180);
latr = latitudeArray.*(pi/180);

%%

%% Evaluate Easting and Northing

q = (1-e^2).*((sin(latr))./(1-((e^2).*(sin(latr)).^2)) - (1/(2*e)).*log((1-e.*sin(latr))./(1+e.*sin(latr))));
q1 = (1-e^2).*((sin(lat1r))./(1-((e^2).*(sin(lat1r)).^2)) - (1/(2*e)).*log((1-e.*sin(lat1r))./(1+e.*sin(lat1r))));
qP = (1-e^2).*((sin(latPr))./(1-((e^2).*(sin(latPr)).^2)) - (1/(2*e)).*log((1-e.*sin(latPr))./(1+e.*sin(latPr))));

beta = asin(q./qP);
beta1 = asin(q1./qP);

m1 = (cos(lat1r))./(sqrt(1-(e^2).*(sin(lat1r)).^2));

Rq = a.*(sqrt(qP./2));

D = (a.*m1)./(Rq.*cos(beta1));

B = Rq.*sqrt(2./(1+(sin(beta1)).*(sin(beta))+(cos(beta1)).*(cos(beta)).*(cos(lonr-lon0r))));

E = B.*D.*cos(beta).*sin(lonr-lon0r);
N = (B./D).*((cos(beta1)).*(sin(beta))-(sin(beta1)).*(cos(beta)).*(cos(lonr-lon0r)));

eastingArray = E + FE;
northingArray = N + FN;

%%

return   