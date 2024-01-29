% Creation of local and global matrices of initial stresses
%
% Input: 
%   elements                .codeNumbers    - code numbers of elements
%                           .vectorX        - directional vector for individual elements
%                           .nelement       - number of elements
%                           .ndofs          - number of unknown displacements
%                           .sections   .E  - Young's modulus of the element
%                                       .v  - Poisson's ratio of the element
%   transformationMatrix    .matrices       - transformation matrices for individual elements 
%                           .lengths        - lengths of elements
%   endForces               .local          - local end forces
% 
% Output:
%   geometricMatrix         .local          - local matrix of initial stresses - cell
%                           .global         - global matrix of initial stresses       
%
% (c) S. Glanc, M. Lepš 2024

function [geometricMatrix]=geometricMatrixFnV2(elements,transformationMatrix,endForces)
    F=endForces.local;
    globalGeometrixMatrix=zeros(elements.ndofs,elements.ndofs);
    localGeometricMatrix{elements.nelement} = 0;
    for cp=1:elements.nelement
        L=transformationMatrix.lengths(cp);%m  
        v = elements.sections.v(cp);
        E = elements.sections.E(cp);
        %End forces
        Fx2=F(7,cp);
        %Constitutive Matrix
        D=sparse(6,6);
        D(1,1) = 1-v;
        D(1,2) = v;
        D(1,3) = v;
        D(2,2) = 1-v;
        D(2,3) = v;
        D(3,3) = 1-v;
        D(4,4) = (1-2*v)/2;
        D(5,5) = (1-2*v)/2;
        D(6,6) = (1-2*v)/2;
        D1=D';
        D = triu(D1.',1) + tril(D1);
        D = D .* (E / ((1+v)*(1-2*v)));

        if ( D(3, 3) ~= 0 )
            kappay = 6 * D(5, 5) / ( D(3, 3) * L * L );
        else
            kappay = 0;
        end

        if ( D(2, 2) ~= 0 ) 
            kappaz = 6 * D(6, 6) / ( D(2, 2) * L * L );
        else
            kappaz = 0;
        end

        Kg=zeros(12,12);
        kappay2 = kappay * kappay;
        kappaz2 = kappaz * kappaz;
        denomy = ( 1. + 2. * kappay ) * ( 1. + 2. * kappay ); 
        denomz = ( 1. + 2. * kappaz ) * ( 1. + 2. * kappaz );
        
        %Stress initial matrix of element
        Kg(2, 2) = ( 4. * kappaz2 + 4. * kappaz + 6. / 5. ) / denomz;
        Kg(2, 6) = ( L / 10. ) / denomz;
        Kg(2, 8) = ( -4. * kappaz2 - 4. * kappaz - 6. / 5. ) / denomz;
        Kg(2, 12) = ( L / 10. ) / denomz;

        Kg(3, 3) = ( 4. * kappay2 + 4. * kappay + 6. / 5. ) / denomy;
        Kg(3, 5) = ( -L / 10. ) / denomy;
        Kg(3, 9) = ( -4. * kappay2 - 4. * kappay - 6. / 5. ) / denomy;
        Kg(3, 11) = ( -L / 10. ) / denomy;

        Kg(5, 5) = L * L * ( kappay2 / 3. + kappay / 3. + 2. / 15. ) / denomy;
        Kg(5, 9) = ( L / 10. ) / denomy;
        Kg(5, 11) = -L * L * ( kappay2 / 3. + kappay / 3. + 1. / 30. ) / denomy;

        Kg(6, 6) = L * L * ( kappaz2 / 3. + kappaz / 3. + 2. / 15. ) / denomz;
        Kg(6, 8) = ( -L / 10. ) / denomz;
        Kg(6, 12) = -L * L * ( kappaz2 / 3. + kappaz / 3. + 1. / 30. ) / denomz;

        Kg(8, 8) = ( 4. * kappaz2 + 4. * kappaz + 6. / 5. ) / denomz;
        Kg(8, 12) = ( -L / 10. ) / denomz;

        Kg(9, 9) = ( 4. * kappay2 + 4. * kappay + 6. / 5. ) / denomy;
        Kg(9, 11) = ( L / 10. ) / denomy;

        Kg(11, 11) = L * L * ( kappay2 / 3. + kappay / 3. + 2. / 15. ) / denomy;
        Kg(12, 12) = L * L * ( kappaz2 / 3. + kappaz / 3. + 2. / 15. ) / denomz;

        minVal = min([Kg(2,2),Kg(3,3),Kg(5,5),Kg(6,6)])/1000;
        Kg(1,1) = minVal;
        Kg(1,7) = -minVal;
        Kg(7,7) = minVal;
        Kg(4,4) = minVal;
        Kg(4,10) = -minVal;
        Kg(10,10) = minVal;
        Kg = Kg .* Fx2/L;



        B1=Kg';
        B = triu(B1.',1) + tril(B1);
        Ksigma=B;
        T=transformationMatrix.matrices{cp};
        T_t=T';
        Ksigma=T_t*Ksigma*T;
        localGeometricMatrix{cp}=Ksigma;
        %Assembly to global sress initial matrix
        kcisla=elements.codeNumbers(cp,:);
                index = kcisla>0 ;
                globalGeometrixMatrix(kcisla(index),kcisla(index)) = globalGeometrixMatrix(kcisla(index),kcisla(index))+ Ksigma(index,index) ;

        geometricMatrix.global =globalGeometrixMatrix;
        geometricMatrix.local = localGeometricMatrix;
    end
end