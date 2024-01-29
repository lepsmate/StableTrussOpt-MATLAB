% Creating code numbers
%
% Input:
%       nodes       .dofs       - number of unknown displacements
%       beams       .nbeams     - number of beams
%                   .nodesHead  - initial node(id)
%                   .nodesEnd   - end node(id)
%
% Output:
%       codes                   - code numbers for individual beams
%
% (c) S. Glanc, 2023

function [codes]=codeNumbersFn(beams,nodes)
   [s,k] = size(nodes.dofs);
    dofs = zeros(s,k);
    m = 0;
    for g = 1:s
        for j = 1:k
            if nodes.dofs(g,j) == 1
                dofs(g,j) = nodes.dofs(g,j) + m;
                m = m+1;
            end
        end
    end
    k1 = k  ;
    k2 = k + 1 ; 
    k3 = 2*k ;
    codes = zeros(beams.nbeams,k3);
    for i = 1:beams.nbeams
    codes(i,1:k1) = dofs(beams.nodesHead(i),:);
    codes(i,k2:k3) = dofs(beams.nodesEnd(i),:);
    end
end
