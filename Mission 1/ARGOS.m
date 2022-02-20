clc
clear all
close all
 

%% variables utiles :
Re = 6378.137e3 % Valeur du demi grand axe en km
f = 1/298.257223563 %  Valeur de l'applatissement
Rp = Re*(1-f); % Valeur du demi-petit axe 

% informations balise 
hb=0; %altitude
latb=44.8059; %Latitude (Bordeaux)
longb=-0.605349;% Longitude (Bordeaux)

GEb=Re+h;
GPb=Rp+h;
phib = atan((GP/GE)*tan(lat));

% coordonnées balise
xb=GEb*cos(phib)*cos(longb); % x dans un repère cartésien 
yb=GEb*cos(phib)*sin(longb);
zb=GPb*sin(phib);
vb=1.5; % vitesse balise

% informations satelite : 
hbs=800e3; %altitude
lats=44.8059; %Latitude (Satelite)
longs=-0.605349;% Longitude (Satelite)

GEb=Re+h;
GPb=Rp+h;
phib = atan((GP/GE)*tan(lat));

% coordonnées cartesiennes satelite
xs=GEs*cos(phis)*cos(longs); 
ys=GEs*cos(phis)*sin(longs);
zs=GPs*sin(phis);

vs=[0 0 8];







