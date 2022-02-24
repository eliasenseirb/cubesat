clc
clear;
close all;
 

%% Méthode des moindres carrés
    %Estimation des coordonnées initiales xk0=(lambdak0,phik0,ftk0)
ftk0 = 868e6;                   %Fréquence d'emission par la plateforme
phik0=44.8378;                  %Latitude de Bordeaux
lambdak0 = -0.594 ;             %Longitude de Bordeaux

% deltafk0 = (rand()-0.5)*ftk0;         %Effet doppler initial (freq reçue - freq transmise)
% frk0 = deltafk0 + ftk0;
% 
% if deltafk0>=0
%     eloignement =1;
% else
%     eloignement=-1;
% end

    %Méthode Gauss-Newton (xk1_est = xk0_est + deltaxk0_est)
h=20;                           %Altitude moyen à Bordeaux
mk= 4;                      %Nombre de mesures de fréquences (doit etre >=3 pour pouvoir avoir assez d'equations)
sigma2k = 2;
zk = [3000 3300 2700 2950];     %Mesures de mk frequences à l'instant k
gk0= H(lambdak0,phik0,h,ftk0,1);
Rk = sigma2k*eye(mk);

J=Jacobien_H(lambdak0,phik0,h,ftk0);








