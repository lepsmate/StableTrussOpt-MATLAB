% Assignment of XY vectors to individual elements of the beam
%
% Input: 
%   beams   .nbeams         - number of beams
%           .disc           - number of elements per beam
%           .XY             - vector determining the XY plane of the beam
%
% Output:
%   XY                      - vector determining the XY plane of the element
%
% (c) S. Glanc, 2024

function [XY]=XYtoElementFn(beams)
for p=1:beams.nbeams        
    for s=1:beams.disc
    XY(s+beams.disc*p-beams.disc,1)=beams.XY(p,1);
    XY(s+beams.disc*p-beams.disc,2)=beams.XY(p,2);
    XY(s+beams.disc*p-beams.disc,3)=beams.XY(p,3);
    end
end
end