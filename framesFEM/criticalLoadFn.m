% Calculation of critical loads from stiffness matrix and stress matrix
%
% Input:
%       geometricMatrix     .global     - global stress matrix
%       stiffnesMatrix      .global     - global stiffness matrix
%       n                               - number of desired modes
% Output:
%       Results             .values     - critical loads;
%                           .vectors    - eigenmodes of deformed structure;
%
% (c) S. Glanc, 2023
function [Results]=criticalLoadFn(stiffnesMatrix,geometricMatrix,n)
    [eigenVectors,eigeinValues]=eigs(stiffnesMatrix,-geometricMatrix,n,'smallestabs');;
    eigeinValues=diag(eigeinValues);
    Results.values = eigeinValues;
    Results.vectors = eigenVectors;
end




            
