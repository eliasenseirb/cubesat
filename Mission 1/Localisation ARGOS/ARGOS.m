clc
clear;
close all;

addpath("..\..\satellite_trajectory_link_budget")

%Simulation satellite

format long;

warning('off','all');



%% Constants
R=6378;         %Earth radius
d2r = pi/180;   %degrees to radians conversion
r2d = 180/pi;   %radians to degrees conversion
u =3998600.4415;  %Gravitational parameter


%% Inputs
H=1600 ;% satellite attitude (km)
a=H+R ;       %Semimajor axis in Km
e=0;          %Eccentricity 0.4
i=90;           %Inclination
omega=90;       %Right ascension of the ascending node
w=90;          %Argument of perigee
theta=90;      %Initial True anomaly
eps=20;         %minimum elevation angle in degrees
UTC='08-sep-2020 17:46:16'; %Satellite launch time
latmin=-90;%-90;
latmax=90;%90;

T=(2*pi)*(sqrt((a^3)/u)); %The Satellite period around the Earth in seconds
dT=T/20; %Time step
N=1; % number of turns around earth
T_f=N*T;
time = 0:dT:T_f-dT;

%% Output


%Satellite Ground Track
disp('The Simulator is running satrack.m function');
[S_lat, S_long, Ecc, E_time] = satrackFoV(a, e, i, omega, w, theta, UTC,time,T);

%% Display
%World map
disp('The Simulator is plotting the Worldmap');
worldmap world
load coastlines
[latcells, loncells] = polysplit(coastlat, coastlon);
plotm(coastlat, coastlon, 'green')
title(sprintf("Satellite trajectory, Simulation Time: %d T",N))
hold on

disp('The Simulator is drawing Satellite trajectory on the Worldmap');
for k = 1:length(S_lat)
    if S_lat(k)>latmin && S_lat(k)<latmax
        plotm(S_lat(k),S_long(k),'rs');
    end
end

hold on





disp('The Simulator is drawing the Field of View along Satellite trajectory and stores its lat/long coordinates in LATC/LONGC matrices ');

for t=1:length(Ecc)
    if S_lat(t)>latmin && S_lat(t)<latmax
                
        rov=FoV(Ecc(t),a, e, w, eps);
        [latc,longc] = scircle1(S_lat(t),S_long(t),rov);
        h2=plotm(latc,longc,'b-');

    end
end
hold off;
 


%% Méthode des moindres carrés
    %Estimation des coordonnées initiales

RE = 6378.137e3;                                  % Taille du demi grand axe en m
f = 1/298.257223563;                              % Aplatissement de l'ellipsoide
RP = RE*(1-f);  
Vs=7e3;                                           % Vitesse du satellite en m/s
c=physconst('LightSpeed');                        % Célérité de la lumière en m/s
hs=1500e3;                                        % Altitude du satellite en basse orbite
GE=RE+hs; 
GP=RP+hs;

ft0 = 868e6;                   % Fréquence d'emission par la plateforme
fr1 = ft0 + 5e3;               % Fréquence reçue au début du premier
fr2 = ft0 - 7e3;
R = 6371e3;                    % Rayon de la terre
etha = 30*d2r;                 % Angle entre l'axe joignant le sommet du cône et le centre de la sphère avec l'axe Z
alpha = 40*d2r;                % Azimuth (angle) du Nord vers l'axe du cone
phis = 50*d2r;                 % Latitude au point sous-sommet du cone
lambdasat = S_long(1)*d2r;


[phi0,lambda0]=init_localisation(Vs,fr1,ft0,hs,R,etha,alpha,phis,lambdasat);

phi0deg = phi0*r2d;
lambda0deg = lambda0*r2d;


% h0=0;                                                                    %Balise en mer par exemple
% lambda0 = 44.833328*d2r;
% phi0 = -0.56667*d2r;
% x0=[lambda0;phi0;h0;ft0];


    % Raffinement itératif (Méthode Gauss-Newton)
mk= 4;                                                                                              % Nombre de mesures de fréquences sur un passage satellite (doit etre >=3 pour pouvoir avoir assez d'equations)
% n_pass_satellite = 4;
% Z = [5000 1000 -5000 -8000;3000 2000 -4000 -9000;1000 500 -1000 -9500;3000 2500 -4000 -7500];       % Matrice contenant sur chaque ligne mk mesures d'effets Doppler sur un passage satellite

S_lat_rad = S_lat*d2r;
S_long_rad = S_long*d2r;

z = [ft0 + 5000;ft0+ 1000;ft0-5000;ft0-8000];
date_mes = [];
g=zeros(mk,1);

for k=1:mk
    if(z(k)>=0)
        g(k,1) = Doppler_func(lambda0,phi0,h0,ft0,1);
    else
        g(k,1) = Doppler_func(lambda0,phi0,h0,ft0,-1);
    end
end

sigma2 = 1;                                    % Variance du bruit
R = sigma2*eye(mk);

J = zeros(mk,4);
for k=1:mk
    J(k,:) = Jacobien_H(lambda0,phi0,h0,ft0,S_long(154+(k-1)),S_lat(154 + (k-1)),H)';
end

delta_x = inv(J'*inv(R)*J)*J'*inv(R)*(z-g);
x1 = x0 + delta_x;


        


% for i=1:mk
%     sigma2 = 2;
%     v = sqrt(sigma2)*randn(1,4);
%     G(i,:) = Z(i,:)-v;
% end




% for j=1:mk
%     gk0(j,:)= H(lambdak0,phik0,h,H(lambdak0,phik0,h,ftk0,1)-zk(j),1);
% end

% J=zeros(mk,size(zk,1));
% for k=1:mk
%     J(k,:)=Jacobien_H(lambdak0,phik0,h,H(lambdak0,phik0,h,ftk0,1)-zk(k));
%     
% end

% Xk_MAT = zeros(mk+1,4);
% Xk_MAT(1,:) = xk0;
% 
% for i=1:mk
%     J= Jacobien_H(lambdak0,phik0,h,ftk0)';
% 
%     if(zk(i)>=0)                                             % Effet doppler positif donc le satellite se rapproche de la balise
%         gk0 = H(lambdak0,phik0,h,ftk0,1);
%     else
%         gk0 = H(lambdak0,phik0,h,ftk0,-1);
%     end
%     
%     dxk0=inv(J'*inv(Rk)*J)*J'*inv(Rk)*(zk(i)-gk0);              % calcul de la petite variation pour raffiner l'estimation des coord
%     
%     xk1=xk0+dxk0;
%     Xk_MAT(i+1,:)=xk1;
%     lambdak0 = xk1(1);
%     phik0 = xk1(2);
%     ftk0 =xk1(4);
%     xk0=[lambdak0 phik0 xk1(3) ftk0];
% end
% 
% Xk_MAT(:,1:2)=Xk_MAT(:,1:2)*180/pi;




