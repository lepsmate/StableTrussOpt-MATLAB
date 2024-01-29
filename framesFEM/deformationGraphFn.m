% Drawing deformed structure
%
% Input:
%       beams       .nbeams     - number of beams
%                   .nodesHead  - initial node(id)
%                   .disc       - number of elements per beam
%       nodes       .dofs       - unknown displacements
%                   .x          - coordinates X
%                   .y          - coordinates Y
%                   .z          - coordinates Z
%                   .nnodes     - number of nodes
%       eigenVector             - eigenmode
%       scale                   - deformation scale
%
% Output:
%       graph                   - deformed structure plot
%
% (c) S. Glanc, 2024
function graph = deformationGraphFn(nodes,beams,eigenVector,scale)
id = 1;
nodes.disc.x = nodes.x;
nodes.disc.y = nodes.y;
nodes.disc.z = nodes.z;
for i = 1:beams.nbeams 
    for d = 1:(beams.disc(i)-1)
        nodes.disc.x(nodes.nnodes + id)=nodes.x(beams.nodesHead(i)) + beams.vectorX(i,1)/beams.disc(i)*(d);
        nodes.disc.y(nodes.nnodes + id)=nodes.y(beams.nodesHead(i)) + beams.vectorX(i,2)/beams.disc(i)*(d);
        nodes.disc.z(nodes.nnodes + id)=nodes.z(beams.nodesHead(i)) + beams.vectorX(i,3)/beams.disc(i)*(d);
        id = id+1;
    end
end
nodes.disc.nnodes = numel(nodes.disc.x); 
nodes.disc.dofs = true(nodes.disc.nnodes,6);                % no kinematic boundary conditions
nodes.disc.dofs(1:nodes.nnodes,:) = nodes.dofs;
nodes.disc.codeNumbers = zeros(size(nodes.disc.dofs)) ;
index = 1;
for i = 1:size(nodes.disc.dofs, 1)
    for j = 1:size(nodes.disc.dofs, 2)
        if nodes.disc.dofs(i, j)
            nodes.disc.codeNumbers(i, j) = index;
            index = index+1;
        end
    end
end

for k = 1:nodes.disc.nnodes
        if nodes.disc.codeNumbers(k,1) ~= 0
            nodes.disc.displacements.x(k,1)  = nodes.disc.x(k) + eigenVector(nodes.disc.codeNumbers(k,1))*scale;
        else
            nodes.disc.displacements.y(k,1)  = nodes.disc.y(k);
        end
        if nodes.disc.codeNumbers(k,2) ~= 0
            nodes.disc.displacements.y(k,1)  = nodes.disc.y(k) + eigenVector(nodes.disc.codeNumbers(k,2))*scale;
        else
            nodes.disc.displacements.y(k,1)  = nodes.disc.y(k);
        end
        if nodes.disc.codeNumbers(k,3) ~= 0
            nodes.disc.displacements.z(k,1)  = nodes.disc.z(k) + eigenVector(nodes.disc.codeNumbers(k,3))*scale;
        else
            nodes.disc.displacements.z(k,1)  = nodes.disc.z(k);
        end
end

id = nodes.nnodes+1;

for n = 1:beams.nbeams
    heads(1,n) = beams.nodesHead(n);    
    for d = 1:beams.disc(n)-1     
        heads(d+1,n) = id;
        ends(d,n) = id;
        id = id+1;
    end
    ends(beams.disc(n),n) = beams.nodesEnd(n);
end

heads = heads(:);
ends = ends(:);

xlim([min(nodes.disc.displacements.x)-0.5,max(nodes.disc.displacements.x)+0.5]);  % to avoid tight limits
ylim([min(nodes.disc.displacements.y)-0.5,max(nodes.disc.displacements.y)+0.5]);  % to avoid tight limits
zlim([min(nodes.disc.displacements.z)-0.5,max(nodes.disc.displacements.z)+0.5]);  % to avoid tight limits
graph = plot3([nodes.disc.displacements.x(heads) nodes.disc.displacements.x(ends)]', ...
     [nodes.disc.displacements.y(heads) nodes.disc.displacements.y(ends)]', ...
     [nodes.disc.displacements.z(heads) nodes.disc.displacements.z(ends)]', ...
     'k','LineWidth',3, 'Color', 'red');
 view([50 15]);
end