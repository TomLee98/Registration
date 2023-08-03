function [R,G,B] = spec2RGB(wave_length,gamma,intensityMax)

if ~exist('gamma','var')
    gamma = 0.8;
end
if ~exist('intensityMax','var')
    intensityMax = 255;
end



if wave_length >= 380 && wave_length <= 439
    Red = -(wave_length - 440) / (440 - 350);
    Green = 0.0;
    Blue = 1.0;
elseif wave_length >= 440 && wave_length <= 489
    Red = 0.0;
    Green = (wave_length - 440) / (490 - 440);
    Blue = 1.0;
elseif wave_length >= 490 && wave_length <= 509
    Red = 0.0;
    Green = 1.0;
    Blue = -(wave_length - 510) / (510 - 490);
elseif wave_length >= 510 && wave_length <= 579
    Red = (wave_length - 510) / (580 - 510);
    Green = 1.0;
    Blue = 0.0;
elseif wave_length >= 580 && wave_length <= 644
    Red = 1.0;
    Green = -(wave_length - 645) / (645 - 580);
    Blue = 0.0;
elseif wave_length >= 645 && wave_length <= 780
    Red = 1.0;
    Green = 0.0;
    Blue = 0.0;
else
    Red = 0.0;
    Green = 0.0;
    Blue = 0.0;
end
  
if wave_length >= 350 && wave_length <= 419
    factor = 0.3 + 0.7*(wave_length - 380)/(420 - 380);
elseif wave_length >= 420 && wave_length <= 700 
    factor = 1.0;
elseif wave_length >= 701 && wave_length <= 780
    factor = 0.3 + 0.7*(780 - wave_length)/(780 - 700);
else
    factor = 0.0;
end
r=intensityMax*(Red*factor)^gamma;
g=intensityMax*(Green*factor)^gamma;
b=intensityMax*(Blue*factor)^gamma;

R=round(r);
G=round(g);
B=round(b);
end