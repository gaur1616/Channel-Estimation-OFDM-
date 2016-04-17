function [H_interpolated] = interpolate(H,pilot_loc,Nfft,method)
% Input: H = Channel estimate using pilot sequence
% pilot_loc = Location of pilot sequence
% Nfft = FFT size
% method = ’linear’/’spline’
% Output: H_interpolated = interpolated channel
if pilot_loc(1)>1
slope = (H(2)-H_est(1))/(pilot_loc(2)-pilot_loc(1));
H = [H(1)-slope*(pilot_loc(1)-1) H]; pilot_loc = [1 pilot_loc];
end
if pilot_loc(end) <Nfft
slope = (H(end)-H(end-1))/(pilot_loc(end)-pilot_loc(end-1));
H = [H H(end)+slope*(Nfft-pilot_loc(end))];
pilot_loc = [pilot_loc Nfft];
end
if lower(method(1))=='l'
    H_interpolated=interp1(pilot_loc,H,(1:Nfft));
else H_interpolated = interp1(pilot_loc,H,(1:Nfft),'spline');
end