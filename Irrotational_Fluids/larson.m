%{
	This is straight from the Mats G. Larson's book
	"The finite element method: Theory, implementation and applications", page 88
	
	2D wing wind simulation from Larson 2013 book
%}

[p,t,potential] = WingFlowSolver2D();

triplot(t(1:3,:)',p(1,:),p(2,:));
hold on


function z = Kappa1(x,y)
	z=0;
	if (x>29.99)
		z=1.e+6; 
	end
end

function z = gD1(x,y)
	z=0;
end

function z = gN1(x,y)
	z=0;
	if (x<-29.99)
		z=1; 
	end
end

function [p,t,phi] = WingFlowSolver2D()
	g = Airfoilg();
	[p,e,t] = initmesh(g,'hmax',0.5);
	A = StiffnessAssembler2D(p,t,inline('1','x','y'));
	[R,r] = RobinAssembler2D(p,e,@Kappa1,@gD1,@gN1);
	phi = (A+R)\r;
end

function [area,b,c] = HatGradients(x,y)
	area=polyarea(x,y);
	b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area;
	c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area;
end

function A = StiffnessAssembler2D(p,t,a)
	np = size(p,2);
	nt = size(t,2);
	A = sparse(np,np); % allocate stiffness matrix
	for K = 1:nt
		loc2glb = t(1:3,K); % local-to-global map
		x = p(1,loc2glb); % node x-coordinates
		y = p(2,loc2glb); % node y-
		[area,b,c] = HatGradients(x,y);
		xc = mean(x); yc = mean(y); % element centroid
		abar = a(xc,yc); % value of a(x,y) at centroid
		AK = abar*(b*b'...
		+c*c')*area; % element stiffness matrix
		A(loc2glb,loc2glb) = A(loc2glb,loc2glb) ...
		+ AK; % add element stiffnesses to A
	end

end

function R = RobinMassMatrix2D(p,e,kappa)
	np = size(p,2); % number of nodes
	ne = size(e,2); % number of boundary edges
	R = sparse(np,np); % allocate boundary matrix
	for E = 1:ne
		loc2glb = e(1:2,E); % boundary nodes
		x = p(1,loc2glb); % node x-coordinates
		y = p(2,loc2glb); % node y-
		len = sqrt((x(1)-x(2))^2+(y(1)-y(2))^2); % edge length
		xc = mean(x); yc = mean(y); % edge mid-point
		k = kappa(xc,yc); % value of kappa at mid-point
		RE = k/6*[2 1; 1 2]*len; % edge boundary matrix
		R(loc2glb,loc2glb) = R(loc2glb,loc2glb) + RE;
	end
end

function r = RobinLoadVector2D(p,e,kappa,gD,gN)
	np = size(p,2);
	ne = size(e,2);
	r = zeros(np,1);
	for E = 1:ne
		loc2glb = e(1:2,E);
		x = p(1,loc2glb);
		y = p(2,loc2glb);
		len = sqrt((x(1)-x(2))^2+(y(1)-y(2))^2);
		xc = mean(x); yc = mean(y);
		tmp = kappa(xc,yc)*gD(xc,yc)+gN(xc,yc);
		rE = tmp*[1; 1]*len/2;
		r(loc2glb) = r(loc2glb) + rE;
	end
end

function [R,r] = RobinAssembler2D(p,e,kappa,gD,gN)
	R = RobinMassMatrix2D(p,e,kappa);
	r = RobinLoadVector2D(p,e,kappa,gD,gN);
end

function g = Airfoilg()
	g=[ 2	17.7218	16.0116	1.5737	1.6675	1	0
		2 16.0116 9.0610 1.6675 1.3668 1 0
		2 9.0610 -0.5759 1.3668 -0.1102 1 0 
		2 -0.5759 -9.5198 -0.1102 -1.8942 1 0
		2 -9.5198 -15.6511 -1.8942 -2.5938 1 0
		2 -15.6511 -18.1571 -2.5938 -1.7234 1 0
		2	-18.1571	-16.9459	-1.7234	0.2051	1	0
		2 -16.9459 -12.4137 0.2051 2.2238 1 0
		2	-12.4137	-5.4090	2.2238	3.4543	1	0
		2	-5.4090	2.8155	3.4543	3.5046	1	0
		2	2.8155	10.6777	3.5046	2.6664	1	0
		2	10.6777	16.3037	2.6664	1.7834	1	0
		2	16.3037	17.7218	1.7834	1.5737	1	0
		2	-30.0000	30.0000	-15.0000	-15.0000	1	0
		2	30.0000	30.0000	-15.0000	15.0000	1	0
		2	30.0000	-30.0000	15.0000	15.0000	1	0
		2	-30.0000	-30.0000	15.0000	-15.0000	1	0]';
end