clc;
clear ;
close all;
%% Parametres "constants"

fT=2.4e9;                                                                   % frequence de la porteuse
Vs=7;                                                                       % Vitesse du sattelite en Km/s par rapport à la surface de la terre
c=299792.458;                                                               % Vitesse du la lumière en Km/s 
H=500;                                                                      % Altitude du sattelite par rapport au centre de la terre (en KM)
R=6371;                                                                     % Rayon de la terre (en KM)

%% Parametres de la simulation

Latitude=90;                                                                % La latitude du site d'emission en degrée (angle positif (0 --> 90) pour hémisphère nord
                                                                            % et negatif (0 --> -90) pour l'hémisphère sud) 
Theta_i=0;                                                                  % La position initiale du sattelite
Tmax=500000;

Theta_M=90-Latitude;

%% simulation 

t=0:1:Tmax;
ThetaS=Vs*t/(R+H)+Theta_i;

%alpha = atand(cosd(ThetaS-Theta_M)./(1+(R*tand(ThetaS-Theta_M)./(cosd(Theta_M)*((R+H)./cosd(ThetaS)-R./(cosd(Theta_M)*cosd(ThetaS-Theta_M)))))));
%fR=fT*(1-Vs*cosd(alpha)/c);
%EffetDoppler=fR-fT;
%figure,
%plot(t,EffetDoppler/1000,"r");
%xlabel('temps (s)')
%ylabel('Décalage Doppler en module (khz)')
%title('Décalage Doppler en fontion du temps pour un site de latitude 90°')

alphaTest=atand(((R+H)*cosd(ThetaS)-R*cosd(Theta_M))./(cosd(Theta_M)*R));
alpha=180-alphaTest+Theta_M+ThetaS;
fR=fT*(1-Vs*cosd(alpha)/c);
EffetDoppler=fR-fT;
figure,
plot(t,EffetDoppler/1000,"r");
xlabel('temps (s)')
ylabel('Décalage Doppler en module (khz)')
title('Décalage Doppler en fontion du temps pour un site de latitude 0° : début pire cas')
