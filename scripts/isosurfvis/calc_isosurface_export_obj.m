clear; clc; close all;

% --- load source data ---
load("tmp/isosurfdata.mat")

% --- prepare inside/outsideness ---
% flip what's inside/outside, if we want to; just a matter of perspective
flip = false;
if flip
    vals   = -vals;
    isoval = -isoval;
end

% --- isosurfaces and isocaps ---
ISO  = isosurface(xyz, xyz, xyz, vals, isoval);
CAPS = isocaps(   xyz, xyz, xyz, vals, isoval);
ISO  = reducepatch(ISO, 0.9);  % reduce the typically overdense mesh by 10%; e.g. gets rid of duplicate vertices

% --- unitcell in R basis ---
UC.vertices = [[0,0,0]; [1,0,0]; [1,1,0]; [0,1,0]; 
               [0,0,1]; [1,0,1]; [1,1,1]; [0,1,1]] - [1,1,1]/2;
UC.faces = [[1,2,3,4]; [5,6,7,8]; [1,2,6,5]; [2,3,7,6]; [3,4,8,7]; [4,1,5,8]];
UC.lines = [[1,2]; [1,4]; [1,5];
            [3,2]; [3,4]; [3,7]; 
            [6,5]; [6,7]; [6,2];
            [8,7]; [8,4]; [8,5]];

% --- isonormals ---
% Calculate Iso-Normals of the surface: when using in Blender, this
% requires that we subsequently turn OFF Object data->Normals->Auto Smooth
% if we want to use Cycles. Including the isonormals in the models gives
% better normal data than could be gotten from the faces alone (because we
% have all the data available here)
NISO = isonormals(xyz, xyz, xyz, vals, ISO.vertices);

% --- transform to cartesian basis ---
if any(Rs{1} ~= [1;0;0]) || any(Rs{2} ~= [0;1;0]) || any(Rs{3} ~= [0;0;1])
    iso_vs  = ISO.vertices; caps_vs = CAPS.vertices; uc_vs = UC.vertices; niso_vs = NISO;
    Rm = [Rs{1} Rs{2} Rs{3}];
    invRmT = transpose(inv(Rm));
    for i = 1:3
        ISO.vertices(:,i)  =  iso_vs(:,1)*Rs{1}(i) +  iso_vs(:,2)*Rs{2}(i) +  iso_vs(:,3)*Rs{3}(i);
        CAPS.vertices(:,i) = caps_vs(:,1)*Rs{1}(i) + caps_vs(:,2)*Rs{2}(i) + caps_vs(:,3)*Rs{3}(i);
        UC.vertices(:,i)   =   uc_vs(:,1)*Rs{1}(i) +   uc_vs(:,2)*Rs{2}(i) +   uc_vs(:,3)*Rs{3}(i);
        % Normal vectors transform contravariantly, i.e. via the inverse-transpose
        % of the transformation matrix of coordinate vectors, see e.g. slides 78++ of
        % https://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-837-computer-graphics-fall-2012/lecture-notes/MIT6_837F12_Lec03.pdf
        NISO(:,i)          = niso_vs(:,1)*invRmT(1,i) + niso_vs(:,2)*invRmT(2,i) + niso_vs(:,3)*invRmT(3,i);
    end
    clear iso_vs caps_vs uc_vs niso_vs
end
% normalize normal vectors to unit length
NISO = NISO./(sqrt(NISO(:,1).^2+NISO(:,2).^2+NISO(:,3).^2)+eps);

% the normals of iso(...) are flipped relative to the .obj format; fix it
ISO.faces  = [ISO.faces(:,3) ISO.faces(:,2) ISO.faces(:,1)];
CAPS.faces = [CAPS.faces(:,3) CAPS.faces(:,2) CAPS.faces(:,1)];

% --- write wavefront OBJ files ---
% surface
OBJ.vertices                 = ISO.vertices;
OBJ.objects(1).type          = 'f';
OBJ.objects(1).data.vertices = ISO.faces;
OBJ.vertices_normal          = NISO;
OBJ.objects(1).data.normal   = ISO.faces;
write_wobj(OBJ,'tmp/isosurface.obj');
clear OBJ
% caps
OBJ.vertices                 = CAPS.vertices;
OBJ.objects(1).type          = 'f';
OBJ.objects(1).data.vertices = CAPS.faces;
write_wobj(OBJ,'tmp/isocaps.obj');
clear OBJ
% unit cell (manual, since write_wobj assumes triangular elements)
fid = fopen('tmp/unitcell.obj','w');
for i = 1:size(UC.vertices,1)
    fprintf(fid, 'v %.5f %.5f %.5f\n', UC.vertices(i,1), UC.vertices(i,2), UC.vertices(i,3));
end
% fprintf(fid, '# %g faces\n', size(UC.faces,1));
% for i = 1:size(UC.faces,1)
%     fprintf(fid, 'f %d %d %d %d\n', UC.faces(i,1), UC.faces(i,2), UC.faces(i,3), UC.faces(i,4));
% end
fprintf(fid, '# %g lines\n', size(UC.lines,1)); % lines only
for i = 1:size(UC.lines,1)
    fprintf(fid, 'l %d %d\n', UC.lines(i,1), UC.lines(i,2));
end
fclose(fid);