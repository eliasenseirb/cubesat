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
    [~,symb] = max(moy_fft); % calcul du max des ffts moyennées
    symbole = M-(symb-1); % calcul des symboles
    compteur=0; % compteur pour savoir combien de max valables on a détecté
    last=0; % variable pour savoir à quel indice on commence
    for i=2:length(symbole)
        if compteur <Np-1 && symbole(i)==symbole(i-1) % on regarde si 2 max consécutifs ont la même abscisse
            compteur=compteur+1; 
            last = i;
            if compteur == Np-1
                start_ind = last - compteur; % on note l'indice du premier max consécutif
                break;
            end
        else % si 2 max consécutifs n'ont pas la même abscisse on reprend tout de 0
            last=0;
            compteur = 0;
        end
    end
    
    sfd_locate = start_ind+Np; % indice du premier symbole sfd
    sig_mat_sfd = sig_Mat(:,sfd_locate:sfd_locate+1).*chirp'; % on applique des up chirps la ou doit se situer notre sfd
    [maxval,~] = max(abs(fft(sig_mat_sfd))); % calcul de l'approx du décalage temporel 
    [~,T_section_location] = max(maxval); 
    K_est = T_section_location+sfd_locate-1;
end

