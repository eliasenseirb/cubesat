clear all;
close all;
clc

% Paramètres
d2r = pi/180;                           % degrés vers radian
r2d = 180/pi;                           % radians vers degrés
c = physconst('lightspeed');
re = 6371e3;                            %  Rayon de la terre
h = 1600e3;                             % Altitude du satellite
r = re +h;

theta_max = 20*d2r;              % Angle d'élévation max
i = 60*d2r;                             % inclinaison du satellite
vit_ang_terre = 2*pi/(24*3600);         % vitesse angulaire de la rotation terrestre
G = 6.6743*10^(-11);                    % Constante de gravitation
M = 5.972*10^24;                        % Masse de la Terre
T = 2*pi*sqrt(r^3/(G*M));               % Période de rotation du satellite
vit_ang_sat = 2*pi/T;                   % vitesse angulaire du satellite
omega_f = vit_ang_sat - vit_ang_terre*cos(i);   % vitesse angulaire du satellite dans l'ECF
diff_psi = [90:-2:-90]*d2r;                              % distance angulaire reliant deux points mesurés sur la surface terrestre
t = 1/(omega_f/(2*d2r));
temps = 0:round(t):90*round(t);

num = -re.*r.*sin(diff_psi).*cos(acos(re./r.*cos(theta_max))-theta_max).*omega_f;
denom = c.*sqrt(re^2 + r^2-2.*re.*r.*cos(diff_psi).*cos(acos(re./r.*cos(theta_max))-theta_max));

Doppler = num./denom;
Doppler_rate = (-re.*r.*omega_f.*cos(diff_psi).*cos(acos(re./r.*cos(theta_max))-theta_max).*omega_f.*denom - num.*c.*(2.*re.*r.*omega_f.*sin(diff_psi).*cos(acos(re./r.*cos(theta_max))-theta_max))./(2.*sqrt(denom)))./denom.^2;



figure()
plot(diff_psi*r2d,Doppler)
xlabel("Angle difference for psi (°)"),ylabel("Doppler")
%title('Evolution du Doppler en fonction de l elevation maximale')

figure()
plot(diff_psi*r2d,Doppler_rate)
xlabel("Angle difference for psi (°)"),ylabel("Doppler rate")
%title("Doppler rate en fonction de l elevation max")

