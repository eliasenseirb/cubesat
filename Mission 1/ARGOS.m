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
zk = [5e6 1e7 5e7 1e8];     %Mesures de mk frequences à l'instant k
for j=1:mk
    gk0(:,j)= H(lambdak0,phik0,h,ftk0,1);
end
Rk = sigma2k*eye(mk);

J=zeros(mk,size(zk,2));
for k=1:mk;
    J(k,:)=Jacobien_H(lambdak0,phik0,h,H(lambdak0,phik0,h,ftk0,1)-zk(k));
    
end
%Jbis=Jacobien_H(lambdak0,phik0,h,ftk0); %calcul de la Jacobienne pour xo;
dxo=inv(J'*inv(Rk)*J)*J'*inv(Rk)*(zk-gk0)'; % calcul de dxo, formule 1.11

xk0=[lambdak0 phik0 h ftk0]
xk1=xk0+dxo';









