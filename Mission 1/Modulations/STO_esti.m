function [lambda_p] = STO_esti(rcd,M,chirp,nu_est,Np)   
    for j=1:Np
        res(:,j) = (rcd(:,j).*conj(chirp)'.*exp(-1j*2*pi*nu_est*(0:M-1)/M)').*exp(-1j*2*pi*(0:M-1).^2/M)';
        
    end
    res = res/sqrt(M);
    [~,ip]=max(abs(res));
    lambda_p = ip-floor(ip);
end

