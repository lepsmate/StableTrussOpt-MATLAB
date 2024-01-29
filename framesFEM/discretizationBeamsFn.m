% Beam discretization - creation of vectors for individual beams and creation   
%                       of unknowns on intermediate elements
%
% Input: 
%   beams       .nbeams         - number of beams
%               .nodesHead      - initial node(id)
%               .nodesEnd       - end node(id)
%               .codeNumbers    - code numbers
%               .disc           - number of elements per beam 
%           	.vectorX        - directional vector for individual beams    
%   nodes       .dofs           - number of unknown displacements
%
% Output:
%   elements    .codeNumbers    - code numbers of elements
%               .vectorX        - directional vector for individual elements
%               .nelement       - number of elements
%
% (c) S. Glanc, 2022
function [elements]=discretizationBeamsFn(beams,nodes)
    [~,k] = size(nodes.dofs);
    k1 = k  ;
    k2 = k + 1 ; 
    k3 = 2*k ;
    cislonezname=max(max(beams.codeNumbers))+1;
for p=1:beams.nbeams
    c=beams.disc(p);
    for s=1:c
    elemVector(s+c*p-c,1)=beams.vectorX(p,1)/c;
    elemVector(s+c*p-c,2)=beams.vectorX(p,2)/c;
    elemVector(s+c*p-c,3)=beams.vectorX(p,3)/c;
    end
end
for p=1:beams.nbeams
% code numbers for the first beam element --------------------------------------------    
    for f=1:k1
        elementsCodeNumber(1+c*(p-1),f)=beams.codeNumbers(p,f);
    end
    if c>1
    for f=k2:k3
        elementsCodeNumber(1+c*(p-1),f)=cislonezname;
        cislonezname=cislonezname+1;
    end
    end
%code numbers for the midle elements --------------------------------------------
    if c>1
    for h=2:c-1
        for f=1:k1
            elementsCodeNumber(h+c*(p-1),f)= elementsCodeNumber(h-1+c*(p-1),f+k1);
        end
        for f=k2:k3
            elementsCodeNumber(h+c*(p-1),f)=cislonezname;
            cislonezname=cislonezname+1;
        end
    end
    end
% code numbers for the last beam element --------------------------------------------
    if c>1
    for f=1:k1
        elementsCodeNumber(c+c*(p-1),f)=elementsCodeNumber(c-1+c*(p-1),f+k1);
    end
    end
    for f=k2:k3
        elementsCodeNumber(c+c*(p-1),f)=beams.codeNumbers(p,f);
    end
%-------------------------------------------------------------------------
end
elements.codeNumbers = elementsCodeNumber;
elements.vectorX = elemVector;
[elements.nelement,~] = size(elemVector);
end
