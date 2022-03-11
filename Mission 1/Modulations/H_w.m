function [H] = H_w(Np,Nsw,w,x,M)

    temp=floor(length(x)/M); % Durée d'un chirp
    x=x(1:temp*M); % on redimensionne x pour le reshape
    
    sig=reshape(x,[M,temp]); % on met en colonne les chirps
    Y=abs(fft(sig));
    Y=Y(:);
    H = sum(max(Y(w*M:w*M+Np+Nsw-1)));
end