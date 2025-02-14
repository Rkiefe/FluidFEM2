%{
	Creator Rodrigo Kiefe, 14 February 2025

	Benchmark lid-driven cavity problem
	Solves the steady state fluid flow for an incompressible viscous fluid, solving for both velocity and pressure
	Following Larson "The Finite Element Method: Theory, Implementation, and Applications" chapter 12.2
%}

close all
clear
clc

% Constants
nu = 0.1;		% viscosity
u_lid = [1,0];  % Velocity of the fluid on the top surface
box_size = [2,2]; % width and hight
position = [0,0];

% Model
model = createpde();
rect1 = [3;
		 4;
		 position(1)-box_size(1)/2;
		 position(1)+box_size(1)/2;
		 position(1)+box_size(1)/2;
		 position(1)-box_size(1)/2;
		 position(2)-box_size(2)/2;
		 position(2)-box_size(2)/2;
		 position(2)+box_size(2)/2;
		 position(2)+box_size(2)/2];

gd = rect1;

ns = char('rect1'); % Create names for the shapes
ns = ns';
sf = 'rect1'; % Set logic for combining shapes

% Create geometry
[dl,bt] = decsg(gd,sf,ns);
geometryFromEdges(model,dl);

% Create a mesh
generateMesh(model,GeometricOrder="quadratic",Hmax=0,Hmin=0); % 0 -> use default values

[~,VE] = area(model.Mesh); % area of each element
p = model.Mesh.Nodes;
t = model.Mesh.Elements;

np = size(p,2);
nt = size(t,2);
ne = numel(unique(t(4:end,:))); 						% Number of midpoints
% ne = max(max(t(4:6,:)))-min(min(t(4:6,:))) + 1;		% Number of midpoints

% Build the matrix blocks
A = sparse(ne,ne);
B1 = sparse(nt,ne);
B2 = sparse(nt,ne);
for k = 1:nt
	nds = t(1:3,k);
	midPoints = t([5,6,4],k)-(np-ne);

	b = zeros(3,1);
	c = zeros(3,1);
	for i = 1:3
		[~,b(i),c(i)] = abc(p,nds,nds(i));
	end

	Sx = [-b(1) + b(2) + b(3);...
		   b(1) - b(2) + b(3);...
		   b(1) + b(2) - b(3)];

	Sy = [-c(1) + c(2) + c(3);...
		   c(1) - c(2) + c(3);...
		   c(1) + c(2) - c(3)];

	A(midPoints,midPoints) = A(midPoints,midPoints) + (Sx*Sx' + Sy*Sy')*VE(k);
	B1(k,midPoints) = -Sx' *VE(k);
	B2(k,midPoints) = -Sy' *VE(k);
end

% Whole matrix
LHS = [nu*A , sparse(ne,ne) , B1';...
	   sparse(ne,ne) , nu*A , B2';...
	   B1 , B2 , sparse(nt,nt)];

% Load vector
rhs = zeros(2*ne+nt,1);

% >> Augment the matrix to introduce the zero mean pressure
last = [zeros(2*ne,1); VE']; % last row and column
LHS = [LHS last; last' 0];
rhs = [rhs; 0];

neq = 2*ne+nt+1;	% Size of the solution to the matrix equation (Ax=b, neq = numel(x))

% Boundary conditions | Find the surface edges and set the velocity value
cavity_nodes = findNodes(model.Mesh,"region","Edge",[1,2,4]);
cavity_nodes = cavity_nodes(cavity_nodes-(np-ne)>0); 				 % Filter for only the midpoints
cavity_nodes = cavity_nodes - (np-ne); % Index start at 1

n_cavity = numel(cavity_nodes);

% Lid surface nodes
lid_nodes = findNodes(model.Mesh,"region","Edge",3);	% All surface nodes
lid_nodes = lid_nodes(lid_nodes-(np-ne)>0); 			% Filter for only the surface nodes
lid_nodes = lid_nodes - (np-ne);

n_lid = numel(lid_nodes);

% 'fixed' uses edge indexing starting at 1 which is the first edge of the whole mesh
% 'fixed' must also point to the index of i + n edges (ne) because of the matrix equation
fixed = [[lid_nodes,cavity_nodes],[lid_nodes,cavity_nodes]+ne]; 

% the free nodes are all of the matrix equation degrees of freedom 
% (all the edges not in a boundary condition + all the triangles + lagrange multiplier)
free = setdiff([1:neq],fixed);

% The boundary condition itself
gvals = zeros(2*(n_lid+n_cavity),1);
gvals(1:n_lid) = u_lid(1);
gvals(n_lid+n_cavity+1:2*n_lid+n_cavity) = u_lid(2);

rhs = rhs(free)-LHS(free,fixed)*gvals;  % shrink vector
LHS = LHS(free,free); 					% shrink matrix

sol = zeros(neq,1);  % allocate solution
sol(fixed) = gvals;  % insert no-slip values
sol(free) = LHS\rhs; % solve linear system

% Solution
velocity = zeros(ne,2);
velocity(:,1) = sol(1:ne); 
velocity(:,2) = sol(1+ne:2*ne); 
pressure = sol(2*ne+1:2*ne+nt);


% Plots
% Elmement centroids for plot
pc = zeros(size(t,2),2);
for k = 1:size(t,2)
	nds = t(1:3,k);
	pc(k,:) = mean(p(:,nds),2)';
end

figure
fig = pdemesh(model);
fig(1).Color = [fig(1).Color,0.1];
hold on
scatter(pc(:,1),pc(:,2),[],pressure,'filled')
cbar = colorbar; cbar.Label.String = "Pressure";
title("Pressure")

figure
fig = pdemesh(model);
fig(1).Color = [fig(1).Color,0.1];
hold on
quiver(p(1,np-ne+1:np)',p(2,np-ne+1:np)',velocity(:,1),velocity(:,2),'b')
title("Velocity field")


function [a,b,c] = abc(p,nds,nd) % [a,b,c]
	nd1 = 0;
	nd2 = 0;
	% Find opposing nodes to nd
	for i = 1:length(nds)
		if nds(i) ~= nd && nd1 == 0
			nd1 = nds(i);
		end

		if nd1 ~= nds(i) && nds(i) ~= nd && nd2 == 0
			nd2 = nds(i);
		end
	end
	
	% Vectors in-plane
	v1 = [p(1,nd)-p(1,nd1), p(2,nd)-p(2,nd1),1];
	v2 = [p(1,nd)-p(1,nd2), p(2,nd)-p(2,nd2),1];

	n = cross(v1,v2);

	a = 1 + (n(1)*p(1,nd) + n(2)*p(2,nd))/n(3);
	b = -n(1)/n(3);
	c = -n(2)/n(3);
end