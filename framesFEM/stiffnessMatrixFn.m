% Creation of local and global stiffness matrices
%
% In:
% elements              .sections.A         - cross-sectional area
%                       .sections.Iy        - moment of inertia about Y axis
%                       .sections.Iz        - moment of inertia about Z axis
%                       .sections.Ix        - torsional moment of inertia
%                       .sections.v         - Poisson's ratio
%                       .sections.E         - Young's modulus
%                       .nelement           - number of elements
%                       .ndofs              - number of unknown internal forces/displacements
%                       .codeNumbers        - element code numbers
% transformationMatrix  .matrices           - transformation matrices for individual elements
%                       .lengths            - lengths of elements
%
%
% Out:
% stiffnesMatrix        .local              - local stiffness matrices - cell
%                   	.global             - global stiffness matrix
%
% (c) S. Glanc, M. Lepš, 2024

function [stiffnesMatrix]=stiffnessMatrixFn(elements,transformationMatrix)
    globalStiffnessMatrix=zeros(elements.ndofs,elements.ndofs,'like',sdpvar);
    localStiffnessMatrix=cell(1,elements.nelement);
    psv = 6;
    K0 = zeros(12,12,'like',sdpvar);
    Kzeros = zeros(psv,psv,'like',sdpvar);
    for cp=1:elements.nelement                  % Id of elements
        A_el=elements.sections.A(cp);
        Iy_el=elements.sections.Iy(cp);
        Iz_el=elements.sections.Iz(cp);
        J_el=elements.sections.Ix(cp);
        E_el=elements.sections.E(cp);
        v_el=elements.sections.v(cp);
        L=transformationMatrix.lengths(cp);
        T=transformationMatrix.matrices{cp};
        T_t=T';
        % Stifness matrix
        % matrix K11
        K11=Kzeros;
        K11(1,1)=A_el/L;
        K11(2,2)=12*Iz_el/(L^3);
        K11(2,6)=6*Iz_el/(L^2);
        K11(3,3)=12*Iy_el/(L^3);
        K11(3,5)=-6*Iy_el/(L^2);
        K11(4,4)=J_el/(2*(1+v_el)*L);
        K11(5,3)=-6*Iy_el/(L^2);
        K11(5,5)=4*Iy_el/L;
        K11(6,2)=6*Iz_el/(L^2);
        K11(6,6)=4*Iz_el/L;

        % matrix K22
        K22=Kzeros;
        K22(1,1)	=	K11(1,1)	;	%	K22(1,1)	=	A_el/L;
        K22(2,2)	=	K11(2,2)	;	%	K22(2,2)	=	12*Iz_el/(L^3);
        K22(2,6)	=	-K11(2,6)	;	%	K22(2,6)	=	-6*Iz_el/(L^2);
        K22(3,3)	=	K11(3,3)	;	%	K22(3,3)	=	12*Iy_el/(L^3);
        K22(3,5)	=	-K11(3,5)	;	%	K22(3,5)	=	6*Iy_el/(L^2);
        K22(4,4)	=	K11(4,4)	;	%	K22(4,4)	=	J_el/(2*(1+v_el)*L);
        K22(5,3)	=	-K11(5,3)	;	%	K22(5,3)	=	6*Iy_el/(L^2);
        K22(5,5)	=	K11(5,5)	;	%	K22(5,5)	=	4*Iy_el/L;
        K22(6,2)	=	-K11(6,2)	;	%	K22(6,2)	=	-6*Iz_el/(L^2);
        K22(6,6)	=	K11(6,6)	;	%	K22(6,6)	=	4*Iz_el/L;



        % matrix K12
        K12=Kzeros;
        K12(1,1)	=	-K22(1,1)	;	%	K12(1,1)	=	-A_el/L;
        K12(2,2)	=	-K22(2,2)	;	%	K12(2,2)	=	-12*Iz_el/(L^3);
        K12(2,6)	=	-K22(2,6)	;	%	K12(2,6)	=	6*Iz_el/(L^2);
        K12(3,3)	=	-K22(3,3)	;	%	K12(3,3)	=	-12*Iy_el/(L^3);
        K12(3,5)	=	-K22(3,5)	;	%	K12(3,5)	=	-6*Iy_el/(L^2);
        K12(4,4)	=	-K22(4,4)	;	%	K12(4,4)	=	-J_el/(2*(1+v_el)*L);
        K12(5,3)	=	K22(5,3)	;	%	K12(5,3)	=	6*Iy_el/(L^2);
        K12(5,5)	=	2*Iy_el/L;	;	%	K12(5,5)	=	2*Iy_el/L;
        K12(6,2)	=	K22(6,2)	;	%	K12(6,2)	=	-6*Iz_el/(L^2);
        K12(6,6)	=	2*Iz_el/L;	;	%	K12(6,6)	=	2*Iz_el/L;

        % matrix K21
        K21=Kzeros;
        K21(1,1)	=	K12(1,1)	;	%	K21(1,1)	=	-A_el/L;
        K21(2,2)	=	K12(2,2)	;	%	K21(2,2)	=	-12*Iz_el/(L^3);
        K21(2,6)	=	-K12(2,6)	;	%	K21(2,6)	=	-6*Iz_el/(L^2);
        K21(3,3)	=	K12(3,3)	;	%	K21(3,3)	=	-12*Iy_el/(L^3);
        K21(3,5)	=	-K12(3,5)	;	%	K21(3,5)	=	6*Iy_el/(L^2);
        K21(4,4)	=	K12(4,4)	;	%	K21(4,4)	=	-J_el/(2*(1+v_el)*L);
        K21(5,3)	=	-K12(5,3)	;	%	K21(5,3)	=	-6*Iy_el/(L^2);
        K21(5,5)	=	K12(5,5)	;	%	K21(5,5)	=	2*Iy_el/L;
        K21(6,2)	=	-K12(6,2)	;	%	K21(6,2)	=	6*Iz_el/(L^2);
        K21(6,6)	=	K12(6,6)	;	%	K21(6,6)	=	2*Iz_el/L;

        K_tuhost=K0 ;
        K_tuhost(1:psv,1:psv)=K11(:,:);
        K_tuhost(psv+1:end,psv+1:end)=K22(:,:);
        K_tuhost(1:psv,psv+1:end)=K12(:,:);
        K_tuhost(psv+1:end,1:psv)=K21(:,:);

        K_tuhost=E_el*K_tuhost;
        localStiffnessMatrix{cp}=K_tuhost;
        K_tuhost=T_t*K_tuhost*T;

        % Assembly
        kcisla=elements.codeNumbers(cp,:);
        index = kcisla>0 ;
        globalStiffnessMatrix(kcisla(index),kcisla(index)) = globalStiffnessMatrix(kcisla(index),kcisla(index))+ K_tuhost(index,index) ;

    end
    stiffnesMatrix.global = globalStiffnessMatrix;
    stiffnesMatrix.local = localStiffnessMatrix;
end