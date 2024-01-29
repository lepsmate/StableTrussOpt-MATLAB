
% Creating vector X for individual beams
%
% Input:
%       beams   .nbeams     - number of beams
%               .nodesHead  - initial node(id)
%               .nodesEnd   - end node(id)
%       nodes   .x          - node coordinates
%               .y          - node coordinates
%               .z          - node coordinates
%
% Output:
%       vector              - directional vector for individual beams
%
% (c) S. Glanc, 2022

function [vector]=beamVectorXFn(beams,nodes)
    vector = zeros(beams.nbeams,3);
    for i = 1:beams.nbeams
        vector(i,1) = nodes.x(beams.nodesEnd(i)) - nodes.x(beams.nodesHead(i));
        vector(i,2) = nodes.y(beams.nodesEnd(i)) - nodes.y(beams.nodesHead(i));
        vector(i,3) = nodes.z(beams.nodesEnd(i)) - nodes.z(beams.nodesHead(i));
    end
end