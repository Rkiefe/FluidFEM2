%{
	Solves the steady state fluid flow for an incompressible fluid
	with irrotational velocity (no vortex)

	This equates to
		v = -grad u
			and
		laplacian u = 0
	where v is the velocity and u is the potential (pressure)

	This code follows the logic from the 2013 book "The finite element method: Theory, implementation and applications"
	from Mats G. Larson.
%}

close all
clear
clc

% Settings
dotSize = 40; % Size of scatter plot dots

% Mesh Parameters
hMax = 0; 	% maximum element size
hMin = 0;	% minimum element size

% PDE model and mesh
model = geometry(L=[20,10],		...
				 diameter=6,	...
				 position=[0,0],...
				 logic="sub");

% >> Plot model
fig_model = figure;
pdegplot(model,"facelabels","on","edgelabels","on")

% >> Generate Mesh
disp("Generating mesh...")
mesh = sketch(model,Hmax=hMax,Hmin=hMin);

disp("Number of elements: "+ mesh.nt)
disp("Number of inside elements: "+ mesh.nInside)
% pause()

% >> Plot Mesh
fig = figure; hold on
pdemesh(model)
% pause()

% >> Border conditions
h = zeros(max(mesh.edgeList(3,:)),1);
h(2) = 1e6;

gD = zeros(max(mesh.edgeList(3,:)),1);
gD(2) = 1;

gN = zeros(max(mesh.edgeList(3,:)),1);
gN(4) = 1;


% >> Run

% Mass matrix
M = boundaryMatrix(mesh,h);

% Stiffness matrix and load vector
A = stiffnessMatrix(mesh);
q = boundaryVector(mesh,gN - h.*gD);

u = (A+M)\q;

% >> Calculate the fluid velocity v
v = zeros(mesh.nt,2);
for k = 1:mesh.nt
	% Nodes of the element
	nds = mesh.t(:,k);

	% X component of the potential
    b = 0; 
    % Y component of the potential
    c = 0;
    
    % Calculate the total potential of the current element
    for ind = 1:length(nds)
        nd = nds(ind);

        [~,bi,ci] = abc(mesh.p,nds,nd);

        v(k,1) = v(k,1) - bi*u(nd);
        v(k,2) = v(k,2) - ci*u(nd);
    end
end
v_norm = sqrt(sum(v.^2,2)); % |v|

% >> Center of each element
pc = zeros(mesh.nt,2);
for k = 1:mesh.nt
	nds = mesh.t(:,k);
	pc(k,:) = mean(mesh.p(1:2,nds),2)';
end

% Plot vector Field of the velocity

figure
pdegplot(model,"facelabels","off"); hold on
plt = triplot(mesh.t',mesh.p(1,:),mesh.p(2,:),'k');
plt.Color = [plt.Color 0.05];

q = quiver(pc(:,1),pc(:,2),v(:,1),v(:,2));

% >> Plot |v|
scatter(pc(:,1),pc(:,2),dotSize,v_norm,'filled');
cbar = colorbar;


% --- Methods ---
function A = stiffnessMatrix(mesh)
	A = zeros(mesh.nv);
	for k = 1:mesh.nt
		% Nodes of the element
		nds = mesh.t(:,k);

		% For each node
		for i = 1:length(nds)
			[~,bi,ci] = abc(mesh.p,nds,nds(i)); % basis function parameters

			for j = i:length(nds) % matrix is symetric
				[~,bj,cj] = abc(mesh.p,nds,nds(j));

				A(nds(j),nds(i)) = A(nds(j),nds(i)) + (bi*bj + ci*cj)*mesh.VE(k);
				A(nds(i),nds(j)) = A(nds(j),nds(i)); % A is symetric
			end
		end
	end
end

function q = boundaryVector(mesh,g)
	q = zeros(mesh.nv,1);
	for e = 1:mesh.ne
		nds = mesh.edgeList(1:2,e);

		% Length of the edge
		l = norm(mesh.p(1:2,nds(2))-mesh.p(1:2,nds(1)));

		q(nds) = q(nds) + 0.5*l*g(mesh.edgeList(3,e));
	end
end

function C = lagVector(mesh) % Unused
	C = zeros(1,mesh.nv);
	for k = 1:mesh.nt
	    nds = mesh.t(:,k);

	    % For each node of the element
	    for ind = 1:length(nds)
	        nd = nds(ind);

	        % a b c d parameters for that element and node
	        [a,b,c] = abc(mesh.p,nds,nd);

	        % Corrdinate of element center
	        pcenter = mean(mesh.p(1:2,nds),2);

	        C(nd) = C(nd) + (a + b*pcenter(1) + c*pcenter(2))*mesh.VE(k);
	    end
	end
end

function M = boundaryMatrix(mesh,h)
	% The 1D mass matrix, along the border of the domain
	M = zeros(mesh.nv);
	for e = 1:mesh.ne
		nds = mesh.edgeList(1:2,e);

		l = norm(mesh.p(1:2,nds(2))-mesh.p(1:2,nds(1)));
		for i = 1:length(nds)
			for j = 1:length(nds)
				M(nds(i),nds(j)) = M(nds(i),nds(j)) + h(mesh.edgeList(end,e))*l/6;
			end
		end
	end 
end


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

function mesh = sketch(model,options)
	arguments
		model;
		options.Hmin = 0;
	 	options.Hmax = 0;
	end

	generateMesh(model, ...
				 "GeometricOrder","linear", ...
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