clear;
close all;
clc

%% Param�tres
SF = 7 ;            %Nombre de bits/symbole
M=2^SF;             
Nbbits = 21000;     %Nombre de bits g�n�r�s
B=600e3;            % Largeur de bande
P= 14;              %Puissance du signal �mis (en Dbm)
Ts=M/B;            %Temps symbole
Ds = 1/Ts;           %Debit symbole
Te = Ts/100;        %P�riode d'�chantillonnage 

SNR_dB = 40;           %Rapport signal sur bruit au niveau du r�cepteur

f0 = 1;             %Frequence min d'un chirp
f1= 250;            %Frequence max d'un chirp
t= 0:Te:1e-3;       % Dur�e d'un chirp
t1 = 1e-6;

sb = randi([0,1],1,Nbbits);

sbMAT = reshape(sb,SF,length(sb)/SF);           %Matrice dont les colonnes sont des sous-sequences de SF bits

Sp = bit2int(sbMAT,SF,true);                    %Convertit en decimal les sequences de SF bits avec bit de poids fort � gauche (en haut de la colonne)

gammap = Sp/B;

s =zeros(size(t));

for k=1:length(t)
    for i=1:length(gammap)
        s(k)=s(k)+exp(1j*2*pi*(t(k)-i*Ts)*(fc(t(k)-i*Ts,gammap(i),B,Ts)));
    end
end

%% Canal
h=1;

y=filter(h,1,s);

%% R�cepteur

Py = mean(abs(y).^2); % Puissance instantann�e du signal re�u
Pbruit = Py/10^(SNR_dB/10);
b = sqrt(Pbruit/2) * (randn(size(y)) + 1i*randn(size(y))); % vecteur de bruit AWG de variance Pbruit

x = y + b;


%% Figures


figure,
subplot 211
plot(t,abs(s)),title("Module de s")
subplot 212
plot(t,angle(s)),title("Phase de s")

figure,
subplot 211
plot(t,abs(x)),title("Module de x")
subplot 212
plot(t,angle(x)),title("Phase de x")

