% Assignment of XY vectors to individual elements of the beam
%
% Input: 
%   beams           .nbeams         - number of beams
%                   .vectorX        - directional vector for individual beams
%   beams.disc                      - beam discretization
%   numberOfBeam                    - number of beams
%   beamVectorXY                    - XY vector for individual beams   
%
% Output:
%   elementVectorXY                 - XY vector for individual elements     
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