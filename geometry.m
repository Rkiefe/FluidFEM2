function model = geometry(options)

arguments
	options.shape = "circle"; 	% or square
	options.logic = "add"; 		% or subtract
	options.L = [15,15]; 		% container size
	options.position = [0,0];	% position of 2nd object
	
	options.diameter = 1; 		% circle diameter
	options.L2 = [2,4];			% second rectangle dimensions
end

model = createpde;

if isequal(options.shape,"circle")
	
	% box with Circle
	box_size = options.L; % width and hight
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

	% Magnet
	% inner bottom circle
	position = options.position;
	radious = options.diameter/2;
	circ2 = [1;
			position(1);
			position(2);
			radious];
	circ2 = [circ2;zeros(6,1)];

	% Combine the geometries
	gd = [rect1,circ2];

	% Create names for the shapes
	ns = char('rect1','circ2');
	ns = ns';

	% Set logic for combining shapes
	if isequal(options.logic,"add")
		sf = 'rect1 + circ2'; 
	else
		sf = 'rect1 - circ2'; 
	end

	% Create geometry
	[dl,bt] = decsg(gd,sf,ns);

	geometryFromEdges(model,dl);

else
	% Just two boxes
	
	% Container
	box_size = options.L; % width and hight
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

	% Magnet
	rectangle_dimensions = options.L2; % width and hight
	position = options.position;

	rect2 = [3;
		4;
		position(1)-rectangle_dimensions(1)/2;
		position(1)+rectangle_dimensions(1)/2;
		position(1)+rectangle_dimensions(1)/2;
		position(1)-rectangle_dimensions(1)/2;
		position(2)-rectangle_dimensions(2)/2;
		position(2)-rectangle_dimensions(2)/2;
		position(2)+rectangle_dimensions(2)/2;
		position(2)+rectangle_dimensions(2)/2];

	% Combine the geometries
	gd = [rect1,rect2];

	% Create names for the shapes
	ns = char('rect1','rect2');
	ns = ns';

	% Set logic for combining shapes
	if isequal(options.logic,"add")
		sf = 'rect1 + rect2'; 
	else
		sf = 'rect1 - rect2'; 
	end

	% Create geometry
	[dl,bt] = decsg(gd,sf,ns);

	geometryFromEdges(model,dl);

end


end
