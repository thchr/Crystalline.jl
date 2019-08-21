clear; clc; close all;
cols=flatcolors(); 
isocol = cols{7}; bndcol = cols{6}; uccol = cols{6};

% load the data
load("isosurfdata.mat")
%isoval=isoval*.95

% we can flip what's inside/outside, if we want to; just a matter of
% perspective
flip = false;
vals = smooth3(vals);
if flip
    vals = -vals;
    isoval = -isoval;
end
% show outline of "connections" between interior and unit cells
showbndlines = true;

% plot the isosurfaces and isocaps
fh=figure; set(fh,'color','w')
subaxis(1,1,1,'M',.25)
ISO=isosurface(xyz,xyz,xyz,vals,isoval);
p=patch('Vertices',ISO.vertices,'Faces',ISO.faces,...
        'FaceVertexCData',ones(size(ISO.vertices,1),1).*isocol,...
        'FaceColor','interp','EdgeColor','none');
CAPS = isocaps(xyz,xyz,xyz,vals,isoval);
hold on
icp=patch('Vertices',CAPS.vertices,'Faces',CAPS.faces,'FaceVertexCData',ones(size(CAPS.vertices,1),1).*isocol,...
    'FaceColor','interp','EdgeColor','none','FaceAlpha',.75);

% plot boundary of isocaps (a bit fragile; e.g. ought to do this "double",
% to ensure caps are always without holes (requires a flip and check)
% even so, can get bugged out if connections are too widespread; works best
% with isolated connections
if showbndlines && ~isempty(CAPS.vertices)
    G=graph(CAPS.faces(:,1),CAPS.faces(:,2)); % there's some chance this is insufficient; but usually, it ought to be good enough...
    Gc=conncomp(G); % find connections between graph
    for j = 1:max(Gc)
        capverts = CAPS.vertices(Gc==j,:);
        bnd = [];
        for s = {[1:2], [2:3], [1,3]} % sort of a bug in 'boundary'; cannot match anything if an entire column has the same value; work around it...
            temp = boundary(capverts(:,s{:}));
            if numel(bnd) < numel(temp)
                bnd = temp;
            end
        end
        if numel(bnd) > 1
            plot3(capverts(bnd,1), capverts(bnd,2), capverts(bnd,3),'-','linewidth',.5,'color',bndcol)
        end
    end
end

% visualize unitcell in R basis (i.e., always cubic)
r0 =   [-1 -1 -1; -1 -1 -1; -1 -1 -1; 1 1 1 ; 1 1 1 ; 1 1 1 ; -1 1 1; -1 1 1; 1 -1 1; 1 -1 1; 1 1 -1; 1 1 -1]/2;
move = [1 0 0   ; 0 1 0   ; 0 0 1   ; -1 0 0; 0 -1 0; 0 0 -1; 0 -1 0; 0 0 -1; -1 0 0; 0 0 -1; -1 0 0; 0 -1 0];
for j=1:size(r0,1)
    plot3(r0(j,1) + [0,move(j,1)],r0(j,2) + [0,move(j,2)],r0(j,3) + [0,move(j,3)],'-k','linewidth',.5,'color',uccol);
end
hold off

% tidy up plot
pad=.1;
xlim([-1,1]/2+[-1,1]*pad); ylim([-1,1]/2+[-1,1]*pad); zlim([-1,1]/2+[-1,1]*pad);
xlabel('\itx'); ylabel('\ity'); zlabel('\itz')
set(gca,'LineWidth',.5)
axis vis3d% off
view(3)
lighting flat
light


% write patch objects to wavefront object format, for use in blender

% Calculate Iso-Normals of the surface: when using in Blender, this
% requires that we subsequently turn OFF Object data->Normals->Auto Smooth
% if we want to use Cycles. Including the isonormals in the models gives
% better normal data than could be gotten from the faces alone (because we
% have all the data available here)
NISO=isonormals(xyz,xyz,xyz,vals,ISO.vertices); 
L=(sqrt(NISO(:,1).^2+NISO(:,2).^2+NISO(:,3).^2)+eps); 
NISO = NISO./L;
ISO.faces=[ISO.faces(:,3) ISO.faces(:,2) ISO.faces(:,1)]; % for some reason, the normals of iso(...) are flipped relative to obj format; fix that
% Make OBJ structure
OBJ.vertices = ISO.vertices;
OBJ.vertices_normal = NISO;
OBJ.objects(1).type='f';
OBJ.objects(1).data.vertices=ISO.faces;
OBJ.objects(1).data.normal=ISO.faces;
write_wobj(OBJ,'isosurface.obj');


% Calculate Iso-Normals of the surface
clear OBJ
CAPS.faces=[CAPS.faces(:,3) CAPS.faces(:,2) CAPS.faces(:,1)]; % for some reason, the normals of iso(...) are flipped relative to obj format; fix that
OBJ.vertices = CAPS.vertices;
OBJ.objects(1).type='f';
OBJ.objects(1).data.vertices=CAPS.faces;
write_wobj(OBJ,'isocaps.obj');
