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

% Geometry type
geo = "circle"; % "circle"

% Mesh Parameters
hMax = 1; 	% maximum element size, set to 0 if you want automatic size
hMin = 1;		% minimum element size, set to 0 if you want automatic size

% PDE model and mesh
if geo == "circle"
	model = geometry(L=[20,10],		...
					 diameter=6,	...
					 position=[0,0],...
					 logic="sub");
else
	model = geometry(shape="rectangle",...
					 L=[20,10],		...
					 L2=[5,2],	...
					 position=[0,0],...
					 logic="sub");
end

% >> Plot model
fig_model = figure;
pdegplot(model,"facelabels","on","edgelabels","on")

% >> Generate Mesh
disp("Generating mesh...")
mesh = sketch(model,Hmax=hMax,Hmin=hMin);

disp("Number of elements: "+ mesh.nt)
disp("Number of inside elements: "+ mesh.nInside)
pause()

% >> Plot Mesh
fig = figure; hold on
pdemesh(model)
pause()

% >> Run fluid simulation
[v_vec,v] = fluid(mesh,geo,hMax=hMax,hMin=hMin);

% >> Center of each element
pc = zeros(mesh.nt,2);
for k = 1:mesh.nt
	nds = mesh.t(:,k);
	pc(k,:) = mean(mesh.p(1:2,nds),2)';
end

% >> Plot vector Field of the velocity
figure
pdegplot(model,"facelabels","off"); hold on
plt = triplot(mesh.t',mesh.p(1,:),mesh.p(2,:),'k');
plt.Color = [plt.Color 0.05];

q = quiver(pc(:,1),pc(:,2),v_vec(:,1),v_vec(:,2));

% Plot |v|
scatter(pc(:,1),pc(:,2),[],v,'filled');
cbar = colorbar;

% --- Methods ---
function [v_vec,v] = fluid(mesh,geo,options)

	arguments
		mesh
		geo = "rectangle";
		options.hMax = 0.1;
		options.hMin = 0;
	end

	% >> Border conditions
	if geo == "rectangle"
		h = zeros(max(mesh.edgeList(3,:)),1);
		h(1) = 1e6;

		gD = zeros(max(mesh.edgeList(3,:)),1);
		% gD(1) = 1;

		gN = zeros(max(mesh.edgeList(3,:)),1);
		gN(6) = 1;
	else
		h = zeros(max(mesh.edgeList(3,:)),1);
		h(1) = 1e6;

		gD = zeros(max(mesh.edgeList(3,:)),1);
		% gD(1) = 1;

		gN = zeros(max(mesh.edgeList(3,:)),1);
		gN(3) = 1;
	end

	% >> Run

	% Mass matrix
	% M = boundaryMatrix(mesh,h);		% Using dense matrix
	M = sparseBoundaryMatrix(mesh,h);	% Using sparse matrix

	% Stiffness matrix
	% A = stiffnessMatrix(mesh);		% Using dense matrix
	A = sparseStiffnessMatrix(mesh);	% Using sparse matrix

	% "Load" vector
	q = boundaryVector(mesh,gN - h.*gD);

	% u = (A+M)\q;				% When using dense matrix
	u = mldivide(A+M,q);		% When using sparse matrix

	% >> Calculate the fluid velocity v
	v_vec = zeros(mesh.nt,2);
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

	        v_vec(k,1) = v_vec(k,1) - bi*u(nd);
	        v_vec(k,2) = v_vec(k,2) - ci*u(nd);
	    end
	end
	v = sqrt(sum(v_vec.^2,2)); % |v|
end % End of fluid simulation

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

function mat = sparseBoundaryMatrix(mesh,h)
	% The 1D mass matrix, along the border of the domain
	% Each edge has 2 nodes, 3 combinations
	% 1: 1+2
	% 3: 1+1 2+2

	%These will be used on the nds arrays for each element to calculate the
	%corresponding aij
	nds_i=[1,1:2];
	nds_j=[2,1:2];

	M = zeros(3,mesh.ne);

	% The actual sparse matrix
 	mat = sparse(mesh.nv,mesh.nv);

	for e = 1:mesh.ne
		nds = mesh.edgeList(1:2,e);

		l = norm(mesh.p(1:2,nds(2))-mesh.p(1:2,nds(1)));
		
		M(:,e) = h(mesh.edgeList(end,e))*l/6;
	end

	% Cross terms
	ind = 1;
    mat = mat + sparse(mesh.edgeList(nds_i(ind),:),mesh.edgeList(nds_j(ind),:),M(ind,:),mesh.nv,mesh.nv);

    % Add the other diagonal (mat is symmetric)
    mat=mat+mat';

    % Add the main diagonal 
    for ind = 2:3
        mat = mat + sparse(mesh.edgeList(nds_i(ind),:),mesh.edgeList(nds_j(ind),:),M(ind,:),mesh.nv,mesh.nv);
    end
end


function mat = sparseStiffnessMatrix(mesh)
    % Each triangle has 3 nodes, 6 combinations:
    % 3: 12 ; 13 ; 23
    % 3: 11 , 22, 33

    %These will be used on the nds arrays for each element to calculate the
    %corresponding aij
    nds_i=[1,1,2,1:3];
    nds_j=[2,3,3,1:3];

    % All nonzero contributions to Aij in the same order as mesh.t
    A = zeros(6,mesh.nt);

    % The actual sparse array
    mat = sparse(mesh.nv,mesh.nv); % to include the lagrange multiplier technique, add +1

    for k = 1:mesh.nt
        
        %k = k_list(ind);
        
        % phi of each node
        b=zeros(3,1);
        c=zeros(3,1);

        nds = mesh.t(:,k);   % Nodes of that element

        % Determine all of the abcd coefficients for this element
        for i = 1:length(nds) 
            [~,b(i),c(i)] = abc(mesh.p,nds,nds(i));
        end

        A(:,k) = mesh.VE(k)*(b(nds_i).*b(nds_j) + ...
                             c(nds_i).*c(nds_j));
    end

    % >> Sum all of the contributions into the sparse array
    
    % Cross terms
    for ind = 1:3
        mat = mat + sparse(mesh.t(nds_i(ind),:),mesh.t(nds_j(ind),:),A(ind,:),mesh.nv,mesh.nv);
    end

    % Add the other diagonal (mat is symmetric)
    mat=mat+mat';

    % Add the main diagonal 
    for ind = 4:6
        mat = mat + sparse(mesh.t(nds_i(ind),:),mesh.t(nds_j(ind),:),A(ind,:),mesh.nv,mesh.nv);
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
