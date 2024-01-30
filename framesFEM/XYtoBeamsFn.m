% Assignment of XY vectors to individual beams
%
% Input: 
%   beams           .nbeams         - number of beams
%                   .vectorX        - directional vector for individual beams

%
% Output:
%   XY                              - XY vector for individual beams     
%
% (c) S. Glanc, 2024


function [XY]=XYtoBeamsFn(beams)
    for b = 1:beams.nbeams
        if beams.vectorX(b,3) == 0
            XY(b,:) = [0 0 1];
        else
            XY(b,:) = [0 1 0];
        end
    end
end