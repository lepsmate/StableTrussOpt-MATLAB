% Sorting positive eigenvalues and eigenvectors
% 
% Input:
%   values = eigenvalues
%   vectors = eigenvectors
% Output:
%   sorted_values = sorted eigenvalues
%   sorted_vectors = sorted eigenvectors
%
% (c) S. Glanc, M. Leps 2024
function [sorted_values,sorted_vectors] = sortValuesVectorFn(values,vectors)
    [~, index] = sort(abs(values));
    sorted_values = values(values(index)>=0);
    sorted_vectors = vectors(:,values(index)>=0);
end