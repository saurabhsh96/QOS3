% Field in media by planar antenna
function [Eth, Ephi] = FieldMedia(n, kd, th, phi, r_obs)
    %Multiplication term to both the components
    mult_term = ((cos(th)).^n).*exp(-1j*kd*r_obs)./r_obs;
    Eth = cos(phi).*mult_term;
    Ephi = -sin(phi).*mult_term;
end