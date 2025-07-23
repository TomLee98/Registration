function r = WaveLength2RGB(w)
%WAVELENGTH2RGB This function transforms light wave length (unit as nm) to
%RGB representation
% Input:
%   - w: 1-by-1 positive double in range [380, 780], unit as nm, indicates
%        the light wave length 
% Output:
%   - r: 1-by-3 positive double integer in range [0, 255], as the RGB
%        representation
% The function based on web page: https://www.johndcook.com/wavelength_to_RGB.html
% See also: Grassmann's law

arguments
    w   (1,1)   double  {mustBeInRange(w, 380, 780)}
end

% Grassmann color law
if w >= 380 && w < 440
    red = -(w - 440) / (440 - 380);
    green = 0.0;
    blue = 1.0;
elseif w >= 440 && w < 490
    red = 0.0;
    green = (w - 440) / (490 - 440);
    blue = 1.0;
elseif w >= 490 && w < 510
    red = 0.0;
    green = 1.0;
    blue = -(w - 510) / (510 - 490);
elseif w >= 510 && w < 580
    red = (w - 510) / (580 - 510);
    green = 1.0;
    blue = 0.0;
elseif w >= 580 && w < 645
    red = 1.0;
    green = -(w - 645) / (645 - 580);
    blue = 0.0;
elseif w >= 645 && w <= 780
    red = 1.0;
    green = 0.0;
    blue = 0.0;
end

% Let the intensity fall off near the vision limits
if w >= 380 && w < 420
    factor = 0.3 + 0.7*(w - 380) / (420 - 380);
elseif w >= 420 && w < 701
    factor = 1.0;
elseif  w >= 701 && w <= 780
    factor = 0.3 + 0.7*(780 - w) / (780 - 700);
end

gamma = 0.80;       % how to make it reasonable

r = [red, green, blue];
r = r.*(r>0);

r = ceil(255*(r*factor).^gamma);

end

