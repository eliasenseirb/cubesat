clear;
close all;
clc

%% Paramètres
SF = 7 ;            %Nombre de bits/symbole
M=2^SF;

Fse=10;
B=600e3;            % Largeur de bande
P= 14;              %Puissance du signal émis (en Dbm)
Ts=M/B;            %Temps symbole
Ds = 1/Ts;         %Debit symbole
Te = Ts/M;        %Période d'échantillonnage
Nb_preambule_up = 5;
Nb_preambule_down=1; % SFD
N_sw = 2; % synchro word
taille_preambule = Nb_preambule_down+Nb_preambule_up+N_sw;
Nb_Chirp = 10; % nombre de Chirp qu'on souhaite dans le signal
SNR_dB = 40;           %Rapport signal sur bruit au niveau du récepteur
Nbbits = SF*Nb_Chirp;     %Nombre de bits générés
time = -Ts/2:Te:Ts/2-Te;                % base de temps sur laquelle les chirps sont générés

%% Transmetteur
sb = randi([0,1],1,Nbbits);     % génération des bits aléatoires
chirp_up= exp(1j*2*pi.*time*B/Ts.*time);    % Chirp up
chirp_down= exp(-1j*2*pi.*time*B/Ts.*time);     %Chirp down
sbMAT = reshape(sb,SF,length(sb)/SF);           %Matrice dont les colonnes sont des sous-sequences de SF bits

Sp = bit2int(sbMAT,SF,true);                    %Convertit en decimal les sequences de SF bits avec bit de poids fort à gauche (en haut de la colonne)
Dp = zeros(size(Sp));
% Modulation DCSS
for k=1:length(Sp)
    if k~=1
        Dp(k) = mod(Dp(k-1)+Sp(k),M);
    else
        Dp(k) = mod(Sp(k),M);
    end
end
Dp=[0,Dp];
gammap = Dp/B;  

%% SURECHANTILLONNE dun facteur 10
%%

preambule=[repmat(chirp_up,1,Nb_preambule_up),repmat(chirp_up,1,N_sw) ,repmat(chirp_down,1,Nb_preambule_down)]; % Préambule 
s=[];
for k=1:length(gammap)
    s = [s exp(1j*2*pi.*time.*fc(time,gammap(k),B,Ts))]; % génération des chirps
end
s=[preambule s];
%su=upsample(s,Fse);
su=s;
%% Canal
h=1;

y=filter(h,1,su);

%% Décalage temporel
decalage_temporel = randi([0,M-1],1); % génération d'un décalage aléatoire
y= [zeros(1,decalage_temporel),y]; % décalage du signal : on rajoute des 0 devant
freq =40e3;
y = y.*exp(-1j*2*pi*freq*(1:length(y)));
%% Récepteur

Py = mean(abs(y).^2); % Puissance instantanée du signal reçu
Pbruit = Py/10^(SNR_dB/10); % Puissance du bruit
%Pbruit=0;
b = sqrt(Pbruit/2) * (randn(size(y)) + 1i*randn(size(y))); % vecteur de bruit AWG de variance Pbruit

x = y + b; %ajout du bruit au signal

%% Estimation du décalage temporel
K_estime = preambule_detect(chirp_up,Nb_preambule_up,N_sw,x,M,decalage_temporel); % estimation du décalage temporel
synchro_temporelle= time_synchro(K_estime-taille_preambule,Nb_preambule_up,N_sw,x,M); % synchronisation temporelle
fprintf("L'écart entre le décalage temporel théorique et celui trouvé est de %i \n",(decalage_temporel-synchro_temporelle))

%%
x = x(synchro_temporelle:end);
DR_esti = doppler_rate_esti(x,M,Nb_preambule_up,chirp_up,Ts); %estimation doppler rate
temp=floor(length(x)/M); % Durée d'un chirp
x=x(1:temp*M); % on redimensionne x pour le reshape

sig_reshaped=reshape(x,[M,temp]); % on met en colonne les chirps
for i=1:Nb_preambule_up
    rdc(:,i) = sig_reshaped(:,i).*exp(-1j*DR_esti*Ts^2*(0:M-1).^2)'; % Dr compensation
end
nu_est = frac_CFO(rdc,Nb_preambule_up); % cfo estimation
lambda_est = STO_esti(rdc,M,chirp_up,nu_est,Nb_preambule_up);
z=sig_reshaped.*chirp_up'; % multiplication par le chirp brut

[max_fft, symbolesEstLoRa]=max(abs(fft(z))); % argmax des FFT
symbolesEstLoRa = M-(symbolesEstLoRa(taille_preambule:end)-1) ;% symboles estimés sans le préambule

% Amélio concavité sert que pour estimer le décalage doppler. 
%[symbole,maxi]= concave(z(:,8:end),symbolesEstLoRa,M); % amélioration de la localisation des max
for k=1:length(symbolesEstLoRa)-1
    symboleEst(k) =mod(symbolesEstLoRa(k+1)-symbolesEstLoRa(k),M); % calcul des symboles Sp 
    %new_symb_est(k)=round(mod(symbole(k+1)-symbole(k),M)); % calcul des symboles Sp
end 

bit_est = int2bit(symboleEst,SF);%bits estimés sans le préambule méthode "classique"
%bit_est2 = int2bit(new_symb_est,SF);%bits estimés sans le préambule méthode algo concavité
BER = mean(abs(sb-bit_est(:)')); % BER avec méthode "classique" (argmax des fft)
%BER2=mean(abs(sb-bit_est2(:)')); % BER avec l'algo de concavité

%% Figures

figure,
subplot 211
plot(abs(s)),title("Module de s")
subplot 212
plot(angle(s)),title("Phase de s")

figure,
subplot 211
plot(abs(fft(x))),title("Module de la fft du signal bruité")
subplot 212
plot(angle(fft(x))),title("Phase de la fft du signal bruité")
