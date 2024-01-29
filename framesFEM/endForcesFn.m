% Calculation of end forces and displacements
%
% Input: 
%   stiffnesMatrix          .global         - global stiffness matrix 
%                           .local          - local stiffness matrix - cell
%   endForces               .global     	- global loads at nodes
%   transformationMatrix    .matrices       - transformation matrices for individual elements 
%   elements                .codeNumbers    - code numbers of elements
%                           .vectorX        - directional vector for individual elements
%                           .nelement       - number of elements
%                           .ndofs          - number of unknown displacements
% 
% Output:
%   localEndForces                          - local end forces on elements      
%   displacements           .local          - local displacements on elements
% 
% (c) S. Glanc, 2023

function [localEndForces,displacements]=endForcesFn(stiffnesMatrix,endForces,transformationMatrix,elements)
    localEndForces = zeros(12,elements.nelement);
    psv = 6;
    psv2 = psv*2;
% Beams displacemnets
    r_global=zeros(elements.ndofs,1);
    r_global(:,1)= stiffnesMatrix.global\endForces.global;
    displacements.global = r_global;

% Elements displacemnets
    r_local=zeros(psv2,elements.nelement);
    for j=1:elements.nelement
        kcisla=elements.codeNumbers(j,:);
            for i=1:psv2
                if kcisla(i)==0
                r_local(i,j)=0;
                else
                r_local(i,j)=r_global(kcisla(i));   
                end
            end
    end

    for i=1:elements.nelement
        T=transformationMatrix.matrices{i}; 
        r=r_local(:,i);
        r_local(:,i)=T*r;
    end
% Elements end forces
    for i=1:elements.nelement
        K_l=stiffnesMatrix.local{i};
        r=r_local(:,i);

        localEndForces(:,i)=K_l*r;
    end
    displacements.local = r_local;
end