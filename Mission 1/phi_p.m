function phi = phi_p(t,gammap,M,Ts,Sp)
    %fonction calculant la phase associée à chaque chirp
    
    if t>=0 && t<T-gammap
        phi = 2*pi*M*((t/(2*Ts))^2+(Sp/M-0.5)*t/Ts);
    else
        phi = 2*pi*M*((t/(2*Ts))^2+(Sp/M-(3/2))*t/Ts);
    end
end

