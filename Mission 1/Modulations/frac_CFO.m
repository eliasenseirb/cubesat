function [nu] = frac_CFO(r,Nc)

    nu = 0; 
    for p=0:Nc-2
        nu=nu+angle(sum(r(:,p+1).*conj(r(:,p+2))))/(2*pi);
    end
    nu=nu/(Nc-1);
end

