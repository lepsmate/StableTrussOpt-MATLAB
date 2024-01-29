% Assignment of sectional characteristics to individual elements of each beam
%
% Input: 
%   beams           .sections           - section id from variable sections
%                   .disc               - number of elements per beam
%   sections        .A                  - cross-sectional area
%                   .Iy                 - moment of inertia about Y axis
%                   .Iz                 - moment of inertia about Z axis
%                   .Ix                 - moment of inertia about X axis
%                   .E                  - Young's modulus of the material
%                   .v                  - Poisson's ratio of the material
%
% Output:
%   elemSections 	.A                  - cross-sectional area
%                   .Iy                 - moment of inertia about Y axis
%                   .Iz                 - moment of inertia about Z axis
%                   .Ix                 - moment of inertia about X axis
%                   .E                  - Young's modulus of the material
%                   .v                  - Poisson's ratio of the material
%
% (c) S. Glanc, 2024

function [elemSections]=sectionToElementFn(sections,beams)
    for p=1:beams.nbeams
            for s=1:beams.disc
                elemSections.A(s+beams.disc*(p-1))=sections.A(beams.sections(p));
                elemSections.Iy(s+beams.disc*(p-1))=sections.Iy(beams.sections(p));
                elemSections.Iz(s+beams.disc*(p-1))=sections.Iz(beams.sections(p)); 
                elemSections.Ix(s+beams.disc*(p-1))=sections.Ix(beams.sections(p)); 
                elemSections.E(s+beams.disc*(p-1))=sections.E(beams.sections(p));
                elemSections.v(s+beams.disc*(p-1))=sections.v(beams.sections(p)); 
            end
    end
end
