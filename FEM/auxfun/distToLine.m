function dist = distToLine(p0,p1,z)
x1 = p0(1,:); y1 = p0(2,:);
x2 = p1(1,:); y2 = p1(2,:);
x3 = z(1); y3 = z(2);
	

num = abs((x2 - x1) .* (y1 - y3) - (x1 - x3) .* (y2 - y1));

den = sqrt((x2 - x1).^ 2 + (y2 - y1).^ 2);

dist = num ./ den;