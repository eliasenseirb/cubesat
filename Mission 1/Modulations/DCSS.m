clear;
close all;
clc

%% Paramètres
SF = 8 ;            %Nombre de bits/symbole
M=2^SF;

Fse=10; % Facteur de sur-échantillonnage
B=125e3;            % Largeur de bande la plus commun pour transmission LoRa
%B=600e3;            % Largeur de bande du sujet
P= 14;              %Puissance du signal émis (en Dbm)
Ts=M/B;            %Temps symbole
Ds = 1/Ts;         %Debit symbole
Te = Ts/M;        %Période d'échantillonnage
Nb_preambule_up = 7; % Preambule
Nb_preambule_down=1; % SFD
N_sw = 2; % synchro word
val_sw = 10; % valeur du mot de synchro
taille_preambule = Nb_preambule_down+Nb_preambule_up+N_sw;
Nb_Chirp = 10; % nombre de Chirp qu'on souhaite dans le signal
SNR_dB = -10;           %Rapport signal sur bruit au niveau du récepteur
Nbbits = SF*Nb_Chirp;     %Nombre de bits générés
time_upsampled = -Ts/2:Te/Fse:Ts/2-Te/Fse;                % base de temps sur laquelle les chirps sont générés
time = -Ts/2:Te:Ts/2-Te;
eb_n0_dB = -15:-9; % Liste des Eb/N0 en dB
eb_n0 = 10.^(eb_n0_dB/10); % Liste des Eb/N0
%% Transmetteur
sb = randi([0,1],1,Nbbits);     % génération des bits aléatoires
chirp_up_upsampled= exp(1j*2*pi.*time_upsampled*B/Ts.*time_upsampled);    % Chirp up sur échantillonné
chirp_down_upsampled= exp(-1j*2*pi.*time_upsampled*B/Ts.*time_upsampled);     %Chirp down sur échantillonné
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

Symbole_sync = [exp(1j*2*pi.*time_upsampled.*fc(time_upsampled,val_sw/B,B,Ts)) exp(1j*2*pi.*time_upsampled.*fc(time_upsampled,val_sw/B,B,Ts))]; % génération des chirps

preambule=[repmat(chirp_up_upsampled,1,Nb_preambule_up),Symbole_sync,repmat(chirp_down_upsampled,1,Nb_preambule_down)]; % Préambule 
s=[];
for k=1:length(gammap)
    s = [s exp(1j*2*pi.*time_upsampled.*fc(time_upsampled,gammap(k),B,Ts))]; % génération des chirps
end
s=[preambule s];
su=s;
for i = 1:length(eb_n0)
    error_cnt=0;
    bit_cnt=0;
    while error_cnt < 100
        %% Canal
        h=1;
        
        y=filter(h,1,su);
        
        %% Décalage temporel
        %decalage_temporel = randi([0,M-1],1); % génération d'un décalage aléatoire
        %y= [zeros(1,decalage_temporel),y]; % décalage du signal : on rajoute des 0 devant
        %deca = [exp(1j*2*pi.*time_upsampled(1:decalage_temporel).*fc(time_upsampled(1:decalage_temporel),0,B,Ts))]; % génération des chirps
        %y=[deca,y];
        %% Récepteur
        
        Py = mean(abs(y).^2); % Puissance instantanée du signal reçu
        %Pbruit = Py/10^(SNR_dB/10); % Puissance du bruit
        Pbruit = Py/10^(eb_n0_dB(i)/10); % Puissance du bruit
        %Pbruit=0;
        b = sqrt(Pbruit/2) * (randn(size(y)) + 1i*randn(size(y))); % vecteur de bruit AWG de variance Pbruit
        
        x = y + b; %ajout du bruit au signal 

        % Ajout du Doppler Rate
        temp = floor(length(x)/(M*Fse)); %durée d'un chirp
        sig = x(1:M*Fse*temp); % on coupe notre signal pour pouvoir le mettre en mode matrice (un chirp = une colonne)
        sig_Mat = reshape(sig,M*Fse,temp); % on transforme la matrice pour avoir des chirps en colonnes
        for j=1:Nb_preambule_up
            sig_Mat_Detect(:,j)=sig_Mat(:,j).*chirp_up_upsampled'; % On multiplie chaque colonne par le chirp brut conjugué
        end
        [~,test1] = max(abs(fft(sig_Mat_Detect))); % calcul de la fft de chaque colonne

        Cr=280; % valeur du Doppler Rate en Hz/s
        t=((0:length(x)-1)*Te/Fse).^2;
        x=x.*exp(1j*pi*Cr*t);
        temp = floor(length(x)/(M*Fse)); %durée d'un chirp
        sig = x(1:M*Fse*temp); % on coupe notre signal pour pouvoir le mettre en mode matrice (un chirp = une colonne)
        sig_Mat = reshape(sig,M*Fse,temp); % on transforme la matrice pour avoir des chirps en colonnes
        for j=1:Nb_preambule_up
            sig_Mat_Detect(:,j)=sig_Mat(:,j).*chirp_up_upsampled'; % On multiplie chaque colonne par le chirp brut conjugué
        end
        [~,test2] = max(abs(fft(sig_Mat_Detect))); % calcul de la fft de chaque colonne
        [test3,~] = concave(sig_Mat_Detect,test2-1,M);
        2/(Ts^2)*(test3(2)-test3(1))
        cd_estime=0;
        Np=Nb_preambule_up;
        for p=1:Np-1 %calcul equation 2.52
            somme_droite = 0;
            for l=p+1:Np
                   somme_droite = somme_droite+2/(Ts^2)*((test3(l)-test3(p)))/(l-p);
            end
            cd_estime=cd_estime+somme_droite/(Np-p);
        end
        cd_estime = cd_estime/(Np-1);
        %% Estimation du décalage temporel
        x2=x;
        x=x(1:Fse:end);% on travail en mode sous échantillonné pour tous les traitements
        Fse=1;
        temp = floor(length(x)/(M*Fse)); %durée d'un chirp
        sig = x(1:M*Fse*temp); % on coupe notre signal pour pouvoir le mettre en mode matrice (un chirp = une colonne)
        sig_Mat = reshape(sig,M*Fse,temp); % on transforme la matrice pour avoir des chirps en colonnes
        for j=1:Nb_preambule_up
            sig_Mat_Detect2(:,j)=sig_Mat(:,j).*chirp_up'; % On multiplie chaque colonne par le chirp brut conjugué
        end
        [~,test2] = max(abs(fft(sig_Mat_Detect2))); % calcul de la fft de chaque colonne
        [test3,~] = concave(sig_Mat_Detect2,test2-1,M);
        cd_estime2=0;
        Np=Nb_preambule_up;
        for p=1:Np-1 %calcul equation 2.52
            somme_droite = 0;
            for l=p+1:Np
                   somme_droite = somme_droite+2/((Ts*10)^2)*((test3(l)-test3(p)))/(l-p+10);
            end
            cd_estime2=cd_estime2+somme_droite/(Np-p);
        end
        cd_estime2 = cd_estime2/(Np-1+10);
        %K_estime = preambule_detect(chirp_up,Nb_preambule_up,N_sw,x,M,1); % estimation du décalage temporel
        %decalage_temporel
        %synchro_temporelle= time_synchro(K_estime,Nb_preambule_up,N_sw,x,M,1,time_upsampled,B,Ts); % synchronisation temporelle
        %synchro_temporelle
        %fprintf("L'écart entre le décalage temporel théorique et celui trouvé est de %i \n",(decalage_temporel-synchro_temporelle))
        
        %%
        DR_esti = doppler_rate_esti(x,M,Nb_preambule_up,chirp_up,Ts); %estimation doppler rate
        temp=floor(length(x)/M); % Durée d'un chirp
        x=x(1:temp*M); % on redimensionne x pour le reshape
        
        sig_reshaped=reshape(x,[M,temp]); % on met en colonne les chirps
        for j=1:Nb_preambule_up %on compense le dr que sur les up chirps du preambule
            rdc(:,j) = sig_reshaped(:,j).*exp(-1j*pi*DR_esti*Ts^2*(0:M-1).^2).'; % Dr compensation
        end
        nu_est = frac_CFO(rdc,Nb_preambule_up,M); % cfo estimation
        lambda_est = STO_esti(rdc,M,chirp_up,nu_est,Nb_preambule_up); % sto estimation
        % Compensation du cfo, sto et dr dans le payload


        z=sig_reshaped.*chirp_up'; % multiplication par le chirp brut conjugué
        
        [max_fft, symbolesEstLoRa]=max(abs(fft(z))); % argmax des FFT
        symbolesEstLoRa = M-(symbolesEstLoRa(taille_preambule+1:end)-1) ;% symboles estimés sans le préambule
        
        % Amélio concavité sert que pour estimer le décalage doppler. 
        %[symbole,maxi]= concave(z(:,8:end),symbolesEstLoRa,M); % amélioration de la localisation des max
        for k=1:length(symbolesEstLoRa)-1
            symboleEst(k) =mod(symbolesEstLoRa(k+1)-symbolesEstLoRa(k),M); % calcul des symboles Sp 
            %new_symb_est(k)=round(mod(symbole(k+1)-symbole(k),M)); % calcul des symboles Sp
        end 
        
        bit_est = int2bit(symboleEst,SF);%bits estimés sans le préambule méthode "classique"
        %bit_est2 = int2bit(new_symb_est,SF);%bits estimés sans le préambule méthode algo concavité
        BER = mean(abs(sb-bit_est(:).')); % BER avec méthode "classique" (argmax des fft)
        error_cnt=error_cnt+sum(sb~=bit_est(:).');
        bit_cnt=bit_cnt + Nbbits; %attention au bit_cnt
    end
    TEB(i) = error_cnt/bit_cnt;
end
%BER2=mean(abs(sb-bit_est2(:)')); % BER avec l'algo de concavité

%% Figures

figure,
semilogy(eb_n0_dB,TEB); 
grid()

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


