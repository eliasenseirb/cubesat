clc
clear all
close all

%% Implementation marche aléatoire :

% Paramètre de diffusion :
Dlambda=5;
Dphi=1;

Vf=100;
Dt= 1000; %Temps depuis le dernier passage satellite en seconde
longorigine=-0.5;
latorigine=44;
horigine=0;
vlamborig=5; % vitesse selon lambda/s
vphiorig=5; %vitesse en phi/s
ft0=1e8;
Q = [Dlambda*Dt 0 0; 0 Dphi*Dt 0; 0 0 Vf];
xo=[longorigine latorigine ft0];
xobis=[longorigine latorigine vlamborig vphiorig ft0];

% Marche aléatoire :

xest=zeros(3,10);
xest(:,1)=xo';

for k=2:10
    xest(:,k)=xest(:,k-1)+sqrt(Q)*randn(3,1);
end

% Marche aléatoire correlée :

N=100
xestbis=zeros(5,N);
xestbis(:,1)=xobis';
Qbis=[0 0 0 0 0; 0 0 0 0 0;0 0 2*Dlambda*Dt 0 0; 0 0 0 2*Dphi*Dt 0; 0 0 0 0 Vf];
M=[1 0 Dt 0 0; 0 1 0 Dt 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1];
for k=2:N
    xestbis(:,k)=M*xestbis(:,k-1)+sqrt(Qbis)*randn(5,1);
end
% Equation marche aléatoire biaisée : p96

de = 110000; % distance en longitude à l'équateur
Nthree=100
xestthree=zeros(3,N);
xestthree(:,1)=xo';
Qthree=[2*Dlambda*Dt 0 0; 0 2*Dphi*Dt 0; 0 0 Vf];
M=[Dt/de*cos(longorigine) 0; 0 Dt/de; 0 0];
alpha=0.3;

v=zeros(2,Nthree); %Ensemble des vecteurs vitesse
vinit1= [vlamborig vphiorig]' %vecteur vitesse initiale pour k=0
v(:,1)=vinit1;
%vinit2= [vlamborig vphiorig]' %vecteur vitesse initiale pour k=1


xestthree(:,2)=xestthree(:,1)+M*vinit1+sqrt(Qthree)*randn(3,1);
% 
% vktil=(xestthree(1:2,2)-xestthree(1:2,1))/Dt;
% vtmp=alpha*vktil+(1-alpha)*v(:,1); %vecteur vitesse pour la prochaine itération
% v(:,2)=vtmp;
% xestthree(:,3)=xestthree(:,2)+M*v(:,2)+sqrt(Qthree)*randn(3,1);

for k=2:Nthree
    vktil=(xestthree(1:2,k)-xestthree(1:2,k-1))/Dt;
    vtmp=alpha*vktil+(1-alpha)*v(:,k-1); %vecteur vitesse pour la prochaine itération
    v(:,k)=vtmp;
    xestthree(:,k+1)=xestthree(:,k)+M*v(:,k)+sqrt(Qthree)*randn(3,1);
end


%% Plot :
figure;
plot(xestbis(1,:),xestbis(2,:),'r*')
title("Simulation trajectoire balise marche aléatoire corrélée")


figure;
plot(xest(1,:),xest(2,:),'r*')
title("Simulation trajectoire balise marche aléatoire")

figure;
plot(xestthree(1,:),xestthree(2,:),'r*')
title("Simulation trajectoire balise marche aléatoire biaisée")
