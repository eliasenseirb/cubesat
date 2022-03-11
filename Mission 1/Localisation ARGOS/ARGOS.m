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

ft0 = 868e6;                   %Fréquence d'emission par la plateforme


% fr1 = ftk0 + 5000;              %fréquence reçue au début d'un passage satellite
% fr2 = ftk0 + 5000*15*60;        %fréquence reçue à la fin du même passage satellite

% phik0 = atan((GE^2*c*(fr1/ftk0 -1))/(sqrt(abs(Vs^2-c^2*(fr1/ftk0 -1)))));   % Expression de la longitude en fonction de fr1 en radian

%disp(Vs^2-c^2*(fr1/ftk0 -1))

h=0;                                                                    %Balise en mer par exemple
lambda0 = 44.833328;
phi0 = -0.56667;
x0=[lambda0 phi0 h ft0];


    % Raffinement itératif (Méthode Gauss-Newton)
mk= 4;                                                                                              % Nombre de mesures de fréquences sur un passage satellite (doit etre >=3 pour pouvoir avoir assez d'equations)
% n_pass_satellite = 4;
% Z = [5000 1000 -5000 -8000;3000 2000 -4000 -9000;1000 500 -1000 -9500;3000 2500 -4000 -7500];       % Matrice contenant sur chaque ligne mk mesures d'effets Doppler sur un passage satellite

zk = [5000 1000 -5000 -8000];
date = [100 300 700 850];  %dates en s



% for i=1:mk
%     sigma2 = 2;
%     v = sqrt(sigma2)*randn(1,4);
%     G(i,:) = Z(i,:)-v;
% end


% sigma2k = 1;                                    % Variance du bruit
% Rk = sigma2k*eye(mk);


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








