clear; clc; close all;
cols=flatcolors(); isocol = cols{7};

% load the data
load("isosurfdata.mat")

% we can flip what's inside/outside, if we want to; just a matter of
% perspective
flip = true;
vals = smooth3(vals);
if flip
    vals = -vals;
    isoval = -isoval;
end

% plot the isosurfaces and isocaps
fh=figure; set(fh,'color','w')
subaxis(1,1,1,'M',.25)
[faces,verts]=isosurface(xyz,xyz,xyz,vals,isoval);
p=patch('Vertices',verts,'Faces',faces,'FaceVertexCData',ones(size(verts,1),1).*isocol,...
    'FaceColor','interp','EdgeColor','none');
[icfaces, icverts] = isocaps(xyz,xyz,xyz,vals,isoval);
hold on
icp=patch('Vertices',icverts,'Faces',icfaces,'FaceVertexCData',ones(size(icverts,1),1).*isocol,...
    'FaceColor','interp','EdgeColor','none');

% visualize unitcell in R basis (i.e., always cubic)
r0 =   [-1 -1 -1; -1 -1 -1; -1 -1 -1; 1 1 1 ; 1 1 1 ; 1 1 1 ; -1 1 1; -1 1 1; 1 -1 1; 1 -1 1; 1 1 -1; 1 1 -1]/2;
move = [1 0 0   ; 0 1 0   ; 0 0 1   ; -1 0 0; 0 -1 0; 0 0 -1; 0 -1 0; 0 0 -1; -1 0 0; 0 0 -1; -1 0 0; 0 -1 0];
for j=1:size(r0,1)
    plot3(r0(j,1) + [0,move(j,1)],r0(j,2) + [0,move(j,2)],r0(j,3) + [0,move(j,3)],'-k','linewidth',.5);
end
hold off

% tidy up plot
pad=.1;
xlim([-1,1]/2+[-1,1]*pad); ylim([-1,1]/2+[-1,1]*pad); zlim([-1,1]/2+[-1,1]*pad);
xlabel('x'); ylabel('y'); zlabel('z')
axis vis3d% off
view(3)
lighting flat
light
