clear;
close all;
clc

%% Param�tres
SF = 7 ;            %Nombre de bits/symbole
M=2^SF;             
Nbbits = 21000;     %Nombre de bits g�n�r�s
Ds = 1e6;           %Debit symbole
Ts=1/Ds;            %Temps symbole
B=600e3;            % Largeur de bande
P= 14;              %Puissance du signal �mis (en Dbm)


f0 = 1;             %Frequence min d'un chirp
f1= 250;            %Frequence max d'un chirp
t= 0:Ts:1e-3;       % Dur�e d'un chirp
t1 = 1e-6;

sb = randi([0,1],1,Nbbits);

sbMAT = reshape(sb,SF,length(sb)/SF);           %Matrice dont les colonnes sont des sous-sequences de SF bits

Sp = bit2int(sbMAT,SF,true);                    %Convertit en decimal les sequences de SF bits avec bit de poids fort � gauche (en haut de la colonne)

gammap = Sp/B;

s= zeros(size(gammap));

for i=1:length(gammap)
    s(i) = exp(1j*phi_p(t-(i-1)*Ts,gammap(i),M,Ts,Sp(i)));
end




ss = chirp(t,f0,t1,f1);

%% Figures
figure()
plot(t,ss)
title('chirp')

figure()
pspectrum(ss,Ds,'spectrogram')

