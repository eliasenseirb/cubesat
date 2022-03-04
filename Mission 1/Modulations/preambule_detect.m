function [K_est] = preambule_detect(chirp,Np,sig,M,val_decalage)
    %detection du preambule
    
    temp = floor(length(sig)/M); %durée d'un chirp
    sig = sig(1:M*temp); % on coupe notre signal pour pouvoir le mettre en mode matrice (un chirp = une colonne)
    sig_Mat = reshape(sig,M,temp); % on transforme la matrice pour avoir des chirps en colonnes
    sig_Mat_Detect = sig_Mat.*chirp'; % On multiplie chaque colonne par le chirp brut conjugué
    
    z = fft(sig_Mat_Detect); % calcul de la fft de chaque colonne
    
    for k=1:temp-1
        moy_fft(:,k) = (abs(z(:,k))+abs(z(:,k+1)))/2; %moyenne des fft sur 2 blocs consécutifs
        %test(:,k)= filter(b,a,[abs(z(:,k));abs(z(:,k+1))]); % filtrage par le filtre RII, sur 2 échantillons consécutifs
    end
    [absci,symb] = max(moy_fft); % calcul du max des ffts moyennées
    symbole = M-(symb-1); % calcul des symboles
    compteur=0; % compteur pour savoir combien de max valables on a détecté
    last=0; % variable pour savoir à quel indice on commence
    thresh = 30; % valeur de seuil
    for i=1:length(symbole)
        if compteur <Np-1 && symbole(i)>thresh
            compteur=compteur+1;
            last = i;
            if compteur == Np-1
                presence = 1;
                start = last - compteur+1;
                break;
            end
        else
            last=0;
            compteur = 0;
        end
    end
    
    if presence % il y a un preambule dans le sig
        sfd_locate = start+Np;
        sig_mat_sfd = sig_Mat(:,sfd_locate:sfd_locate+1).*chirp'; % on applique des up chirps la ou doit se situer notre sfd
        [maxval,K_est] = max(abs(fft(sig_mat_sfd))); % calcul de l'approx du décalage temporel 
        K_est = K_est+start*M;
    else
        K_est=0;
    end
end

