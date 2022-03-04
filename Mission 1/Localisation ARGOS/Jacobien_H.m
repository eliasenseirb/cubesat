function [J_H]=Jacobien_H(lambda,phi,h,ft)
%Jacobienne de la fonction d'observation Doppler H

tau = 120e-3;                                       % durée d'une mesure
RE = 6378.137e3;                                  % Taille du demi grand axe en m
f = 1/298.257223563;                              % Aplatissement de l'ellipsoide
RP = RE*(1-f);                                    % Valeur du demi-petit axe 
                                        
c=physconst('LightSpeed');                        % Célérité de la lumière en m/s
GE=RE+h; 
GP=RP+h;
theta = atan((GP/GE)*tan(phi));                   % Latitude paramétrique (ou réduite)

% Coordonnées cartésiennes
x=GE*cos(theta)*cos(lambda); 
y=GE*cos(theta)*sin(lambda);
z=GP*sin(theta);

% informations balise 
hb= h;                       %altitude
phib=phi;                   %Latitude (Bordeaux)
lambdab= lambda ;           % Longitude (Bordeaux)

GEb=RE+hb;
GPb=RP+hb;
thetab = atan((GP/GE)*tan(phib));

% coordonnées balise
xb=GEb*cos(thetab)*cos(lambdab); % x dans un repère cartésien 
yb=GEb*cos(thetab)*sin(lambdab);
zb=GPb*sin(thetab);

% informations satellite : 
hs=1500e3;               %altitude du satellite en basse orbite
phis=44.8059;           %Latitude (Satellite)
lambdas=-0.605349;      % Longitude (Satellite)

GEs=RE+hs;
GPs=RP+hs;
thetas = atan((GPs/GEs)*tan(phis));

% coordonnées cartesiennes satellite
xs=GEs*cos(thetas)*cos(lambdas); 
ys=GEs*cos(thetas)*sin(lambdas);
zs=GPs*sin(thetas);
vs=7e3;                         %vitesse du satellite

zsf = zs + vs*tau;              %coordonnée selon z du satellite en fin de comptage
xsf = xs;                       %coordonnées en x et y sont inchangées puisque le satellite ne bouge que selon z
ysf = ys;

Rf= sum([xsf ysf zsf].*[xb yb zb]); %distance satellite-balise en fin de comptage
Rd= sum([xs ys zs].*[xb yb zb]);    %distance satellite-balise en début de comptage
vr = (Rf-Rd)/tau;                   % approximation de la vitesse radiale su satellite

dH_dxyz = (ft/c*tau)*[((x-xsf)/Rf-(x-xs)/Rd) ((y-ysf)/Rf-(y-ys)/Rd) ((z-zsf)/Rf-(z-zs)/Rd)];

dH_dft = 1- vr/c;

dtheta_dphi = (GP/GE)/((1-GP^2/GE^2)*cos(phi)^2+GP^2/GE^2);

dxyz_dlambdaphih = [-GE*cos(theta)*sin(lambda) -GE*sin(theta)*cos(lambda)*dtheta_dphi cos(theta)*sin(lambda);GE*cos(theta)*cos(lambda) -GE*sin(theta)*sin(lambda)*dtheta_dphi cos(theta)*sin(lambda);0 GP*cos(theta)*dtheta_dphi sin(theta)];

J_H = [dH_dxyz*dxyz_dlambdaphih dH_dft];


end
