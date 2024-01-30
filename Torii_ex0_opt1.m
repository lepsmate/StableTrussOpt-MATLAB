% Optimization of 3D truss structure employing frame elements
%
% Using YALMIP (https://yalmip.github.io/) and 
% PENLAB (https://web.mat.bham.ac.uk/kocvara/penlab/) libraries
%
% Topology taken from Torii, A.J., Lopez, R.H. & Miguel, L.F.F. 
% Modeling of global and local stability in optimization of truss-like structures using frame elements. 
% Struct Multidisc Optim 51, 1187–1198 (2015). https://doi.org/10.1007/s00158-014-1203-y
%
% (c) S. Glanc, M. Leps 2024
%
clear
close all
clc
addpath 'framesFEM'
%% Sections
%
nSections = 6 ; % number of cross-sections

r0 = 10*10^-3 % minimal radius 10 mm ; Original settings from (Torii et al., 2015) is r0 = 40*10^-3 mm

startingA = (15*pi*(ones(nSections,1)*r0).^2)/64 ; % cross section area

outA = sdpvar(nSections,1);
assign(outA, startingA);

sections.A = outA ; 

for j = 1:nSections
    sections.Iy(j,1) = outA(j)^2*113/(60*pi)  ;     % moment of inertia as a function of an area
end

sections.Iz = sections.Iy;                          % moment of inertia 
sections.Ix = sections.Iy*2;                        % moment of inertia 
sections.E  = ones(nSections,1) * 210*10^9;         % Young moduli of individual elements of beam
sections.v  = ones(nSections,1) * 0.3;              % Young moduli of individual elements of beam

ndisc = 8;                                          % number of FEM elements per beam

%% Nodes
nodes.x = [0;2;0;2;4];                              % x coordinates of nodes
nodes.y = [0;0;0;0;0];                              % y coordinates of nodes
nodes.z = [0;0;2;2;2];                              % z coordinates of nodes
nodes.nnodes = numel(nodes.x);                      % number of nodes

%% Supports
kinematic.x.nodes = [1;3];                          % node indices with restricted x-direction displacements
kinematic.y.nodes = [1;3];                          % node indices with restricted y-direction displacements
kinematic.z.nodes = [1;3];                          % node indices with restricted z-direction displacements
kinematic.rx.nodes = [];                            % node indices with restricted x-direction rotations
kinematic.ry.nodes = [];                            % node indices with restricted y-direction rotations
kinematic.rz.nodes = [1;3];                         % node indices with restricted z-direction rotations

nodes.dofs = true(nodes.nnodes,6);                  % no kinematic boundary conditions
nodes.dofs(kinematic.x.nodes,1) = false;            % mark prevented movement in x-direction
nodes.dofs(kinematic.y.nodes,2) = false;            % mark prevented movement in y-direction
nodes.dofs(kinematic.z.nodes,3) = false;            % mark prevented movement in z-direction
nodes.dofs(kinematic.rx.nodes,4) = false;           % mark prevented rotation in x-direction
nodes.dofs(kinematic.ry.nodes,5) = false;           % mark prevented rotation in y-direction
nodes.dofs(kinematic.rz.nodes,6) = false;           % mark prevented rotation in z-direction

%% Beams
beams.nodesHead = [1;3;3;2;4;2];                    % elements starting nodes
beams.nodesEnd  = [2;2;4;4;5;5];                    % elements ending nodes
beams.nbeams = numel(beams.nodesHead);              % number of beams

% Beams sections
beams.disc      = ones(beams.nbeams,1)*ndisc;       % disretization of beams
beams.sections  = 1:nSections;              % id of sections to beams

%% Ground structure                                 % graph with nodes and beams
plot3([nodes.x(beams.nodesHead) nodes.x(beams.nodesEnd)]', ...
     [nodes.y(beams.nodesHead) nodes.y(beams.nodesEnd)]', ...
     [nodes.z(beams.nodesHead) nodes.z(beams.nodesEnd)]', ...
     'k','LineWidth',3);
hold on;
scatter3(nodes.x, nodes.y, nodes.z, 'black', 'filled', 'o');
axis equal;
xlim([min(nodes.x)-0.1,max(nodes.x)+0.1]);          % to avoid tight limits
ylim([min(nodes.y)-0.1,max(nodes.y)+0.1]);          % to avoid tight limits
zlim([min(nodes.z)-0.1,max(nodes.z)+0.1]);          % to avoid tight limits
grid on

%% Loads
loads.y.nodes = [];                                 % node indices with x-direction forces
loads.y.value = [];                                 % magnitude of the x-direction forces
loads.x.nodes = [];                                 % node indices with y-direction forces
loads.x.value = [];                                 % magnitude of the y-direction forces 
loads.z.nodes = [5];                                % node indices with z-direction forces
loads.z.value = [-1000];                            % magnitude of the z-direction forces 

loads.rx.nodes = [];                                % node indices with x-direction moments
loads.rx.value = [];                                % magnitude of the x-direction moments
loads.ry.nodes = [];                                % node indices with y-direction moments
loads.ry.value = [];                                % magnitude of the y-direction moments 
loads.rz.nodes = [];                                % node indices with z-direction moments
loads.rz.value = []; 

%% Force vector
forceVector = sparse([loads.x.nodes*6-5; loads.y.nodes*6-4; loads.z.nodes*6-3;loads.rx.nodes*3-2; loads.ry.nodes*3-1; loads.rz.nodes*3], ...
                     1, ...
                     [loads.x.value; loads.y.value; loads.z.value;loads.rx.value; loads.ry.value; loads.rz.value], ...
                     nodes.nnodes*6, 1);
f = forceVector(reshape(reshape(nodes.dofs.',[],1).', 1, [])');

%% FEM
nodes.ndofs = sum(sum(nodes.dofs));                 % number of unknown dofs

beams.vectorX = beamVertexFn(beams,nodes);          % X vectors for beams
beams.codeNumbers = codeNumbersFn(beams,nodes);     % beam code numbers
beams.XY = XYtoBeamsFn(beams);                      % define the XY plane for beams

elements = discretizationBeamsFn(beams,nodes);      % discretisation of beams into elements
elements.XY = XYtoElementFn(beams);                 % XY plane assignment for elemtns  
elements.sections = sectionToElementFn(sections,beams); % sections assignment for elemtns  
elements.ndofs = max(max(elements.codeNumbers));    % number of unknown dofs for elements  

% Linear analysis 
endForces.global = sparse(elements.ndofs,1);                                
endForces.global(1:max(max(beams.codeNumbers))) = f;    % assembly force vector
transformationMatrix = transformationMatrixFn(elements);    % local and global transformation matrix and lengths in struct
stiffnesMatrix = stiffnessMatrixFn(elements,transformationMatrix); % local and global stifness matrix in struct

assign(outA, startingA);
realMatrix.global = value(stiffnesMatrix.global) ;
for i=1:size(stiffnesMatrix.local,2)
    realMatrix.local{i} = value(stiffnesMatrix.local{i}) ;
end

[endForces.local, displ] = endForcesFn(realMatrix,endForces,transformationMatrix,elements); % Solving FEM -> u = K\f

%% Optimization
geometricMatrix = geometricMatrixFnV2(elements,transformationMatrix,endForces);             % Geometric stiffness matrix

nonnegativeCrossSections = outA>=0.00001; % Limit for minimal cross-sections

Lambda = 10.0 ;  % critical coeficient
LMIconstraint = stiffnesMatrix.global+Lambda*geometricMatrix.global >= 0;

volume = sum(elements.sections.A .* transformationMatrix.lengths);

ops = sdpsettings('solver','penlab','verbose',1,'debug',1);
optimize([nonnegativeCrossSections, LMIconstraint], ...
          volume, ops);

%% Post-procesing
fprintf('Optimized objective: %f\n', value(volume));
disp('= Section areas ========')
fprintf('A%d = %f\n', [(1:numel(outA))', value(outA)]');

realMatrix.global = value(stiffnesMatrix.global) ;
Results = criticalLoadFn(realMatrix.global,geometricMatrix.global,100);

[sortedValues,sortedVectors]= sortValuesVectorFn(Results.values,Results.vectors);

disp('= Five smallest critical loads =====')
fprintf('lambda%d = %f\n', [(1:numel(sortedValues(1:5)))', value(sortedValues(1:5))]');

for i = 1:beams.nbeams
    idR = 1:12;
    idC = (i-1)*ndisc + 1 : i*ndisc;
    endForces.beams{i} = endForces.local(idR, idC);
    N(1,i)=endForces.beams{i}(7,1);
end

N = full(N'/1000);
disp('= Normal forces ========')
fprintf('N%d = %f kN \n', [(1:numel(N))', value(N)]');
sigma = N./value(outA)/1000;
disp('= Normal stresses ======')
fprintf('sigma%d = %f MPa \n', [(1:numel(sigma))', value(sigma)]');
% Graph of deformed structure
graphValueId = 1;
graphScale = 3;
graph = deformationGraphFn(nodes,beams,sortedVectors(:,graphValueId),graphScale);