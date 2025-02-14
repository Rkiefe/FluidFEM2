%{
	Solves the steady state fluid flow for an incompressible viscous fluid 

	Solving for the velocity and pressure

	Eq.: grad p - mu laplacian u = rho f

	f is a load force

	Following Larson chapter 12.2
%}

close all
clear
clc

model = createpde();
box_size = [2,2]; % width and hight
position = [0,0];

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

% Create names for the shapes
ns = char('rect1');
ns = ns';

% Set logic for combining shapes
sf = 'rect1'; 

% Create geometry
[dl,bt] = decsg(gd,sf,ns);
geometryFromEdges(model,dl);

generateMesh(model,GeometricOrder="linear",Hmax=0,Hmin=0);

p = model.Mesh.Nodes;
t = model.Mesh.Elements;

tic
t2e = Tri2Edge(p,t); % triangle-to-edge adjacency

nt = size(t,2); % number of triangles
ne = max(t2e(:)); % number of edges

[A11,B1,B2,areas] = NCAssembler(p,t2e,t); % assemble

nu = 0.1; % viscosity parameter

LHS = [ nu*A11 sparse(ne,ne) B1';
	 	sparse(ne,ne) nu*A11 B2';
	 	B1 B2 sparse(nt,nt)]; % LHS matrix

rhs = zeros(2*ne+nt,1); % RHS vector

last = [zeros(2*ne,1); areas]; % last row and column
LHS = [LHS last; last' 0];
rhs = [rhs; 0];

[xmid,ymid,edges] = EdgeMidPoints(p,t2e,t);
fixed=[]; % fixed nodes
gvals=[]; % nodal values of g

for i=1:length(edges) % loop over edges
	n=edges(i); % edge (ie. node) number
	x=xmid(i); % x-coordinate of edge mid-point
	y=ymid(i); % y-
	if (x<-0.99 | x>0.99 | y<-0.99 | y>0.99) % boundary
		fixed=[fixed; n; n+ne]; % fix velocity nodes
		u=0; v=0; % bc values
		if (y>0.99), u=1; end % u=1,v=0 on lid
		gvals=[gvals; u; v];
	end
end

neq=2*ne+nt+1; % number of equations
free=setdiff([1:neq],fixed);
rhs=rhs(free)-LHS(free,fixed)*gvals; % shrink vector
LHS=LHS(free,free); % shrink matrix
sol=zeros(neq,1); % allocate solution
sol(fixed)=gvals; % insert no-slip values
sol(free)=LHS\rhs; % solve linear system

U=sol(1:ne); V=sol(1+ne:2*ne); P=sol(2*ne+1:2*ne+nt);
toc

% Plots
% Elmement centroids for plot
pc = zeros(size(t,2),2);
for k = 1:size(t,2)
	nds = t(1:3,k);
	pc(k,:) = mean(p(:,nds),2)';
end

% Velocity plot
figure
fig = pdemesh(model);
fig(1).Color = [fig(1).Color,0.1];
hold on
quiver(xmid,ymid,U',V','b')
title("Larson - velocity")

% Pressure plot
figure
fig = pdemesh(model);
fig(1).Color = [fig(1).Color,0.1];
hold on
scatter(pc(:,1),pc(:,2),[],P,'filled')
colorbar
title("Larson - pressure")

function [A11,B1,B2,areas] = NCAssembler(p,t2e,t)
	nt=size(t,2);
	ne=max(t2e(:));
	A11=sparse(ne,ne);
	B1=sparse(nt,ne);
	B2=sparse(nt,ne);
	areas=zeros(nt,1);
	for i=1:nt
		vertex=t(1:3,i);
		x=p(1,vertex);
		y=p(2,vertex);
		[area,Sx,Sy]=CRGradients(x,y);
		edges=t2e(i,:);
		A11(edges,edges)=A11(edges,edges)+(Sx*Sx'+Sy*Sy')*area;
		B1(i,edges)=-Sx'*area;
		B2(i,edges)=-Sy'*area;
		areas(i)=area;
	end
end

function [xmid,ymid,e] = EdgeMidPoints(p,t2e,t)
	i=t(1,:); j=t(2,:); k=t(3,:); % triangle vertices
	t2e=t2e(:); % all edges in a long row
	start=[j i i]; % start vertices of all edges
	stop =[k k j]; % stop
	xmid=(p(1,start)+p(1,stop))/2; % mid point x-coordinates
	ymid=(p(2,start)+p(2,stop))/2; % y-
	[e,idx]=unique(t2e); % remove duplicate edges
	xmid=xmid(idx); % unique edge x-coordinates
	ymid=ymid(idx); % y-
end

function edges = Tri2Edge(p,t)
	np=size(p,2); % number of vertices
	nt=size(t,2); % number of triangles
	i=t(1,:); % i=1st vertex within all elements
	j=t(2,:); % j=2nd
	k=t(3,:); % k=3rd
	A=sparse(j,k,-1,np,np); % 1st edge is between (j,k)
	A=A+sparse(i,k,-1,np,np); % 2nd	(i,k)
	A=A+sparse(i,j,-1,np,np); % 3rd	(i,j)
	A=-((A+A.')<0);
	A=triu(A); % extract upper triangle of A
	[r,c,v]=find(A); % rows, columns, and values(=-1)
	v=[1:length(v)]; % renumber values (ie. edges)
	A=sparse(r,c,v,np,np); % reassemble A
	A=A+A'; % expand A to a symmetric matrix
	edges=zeros(nt,3);
	for k=1:nt
		edges(k,:)=[A(t(2,k),t(3,k))
					A(t(1,k),t(3,k))
					A(t(1,k),t(2,k))]';
	end
end

function [area,Sx,Sy] = CRGradients(x,y)
	[area,b,c]=HatGradients(x,y);
	Sx=[-b(1)+b(2)+b(3); b(1)-b(2)+b(3); b(1)+b(2)-b(3)];
	Sy=[-c(1)+c(2)+c(3); c(1)-c(2)+c(3); c(1)+c(2)-c(3)];
end

function [area,b,c] = HatGradients(x,y)
	area=polyarea(x,y);
	b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area;
	c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area;
end







function mesh = sketch(model,options)
	arguments
		model;
		options.Hmin = 0;
	 	options.Hmax = 0;
	 	options.GeometricOrder = "linear";
	end

	generateMesh(model, ...
				 GeometricOrder,options.GeometricOrder, ...
				 "Hmin",options.Hmin, ...
				 "Hmax",options.Hmax, ...
				 "Hgrad",1.5);

	[~,AE] = area(model.Mesh); % area of each element

	% >> Process mesh to get a list of edges
	[p,edgeList,t] = meshToPet(model.Mesh);

	t = t(1:3,:); % only store the node indices

	% edge list: nd1, nd2, edge index
	edgeList = edgeList([1,2,5],:);

	nv = length(p); 		% Number of nodes
	nt = length(t); 		% number of elements
	ne = length(edgeList); 	% number of edges

	% Elements of object
	InsideElements = []; nInside = 0;
	try
		InsideElements = findElements(model.Mesh,"region",Face=2);
		nInside = numel(InsideElements);
	end

	mesh.p = p; clear p
	mesh.t = t; clear t
	mesh.VE = AE; clear AE
	mesh.edgeList = edgeList; clear edgeList
	mesh.nv = nv;
	mesh.nt = nt;
	mesh.ne = ne;
	mesh.InsideElements = InsideElements; clear InsideElements
	mesh.nInside = nInside;
end

