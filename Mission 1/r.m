function R_w = r(z,w)
    R_w = abs(sum(z'.*exp(-1j*((0:length(z)-1)-1)*w)));
end

