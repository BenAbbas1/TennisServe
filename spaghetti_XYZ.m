function [ x,y,z ] = spaghetti_XYZ( verts, width, n )
% based on tubecoords from the streamlines function in matlab
    d1 = diff(verts);
    zindex = find(~any(d1,2));
    verts(zindex,:) = [];
    
    if size(verts,1)<2
        x = []; y = []; z = [];
        return;
    end
    
    d1 = diff(verts);
    
    numverts = size(verts,1);
    unitnormals = zeros(numverts,3);
    
    % Radius of the tube is different in each direction to model
    % both linguini and spaghetti shapes
    if (size(width,2) > 1)
        if (size(width,1) == 3 && size(width,2) > 3)
            width = width';
        end
    end
    R  = zeros(numverts,3);
    switch numel(width)
        case 1
            R(:,1)  = (width/2).*ones(numverts,3);
            R(:,2)  = (width/2).*ones(numverts,3);
            R(:,3)  = (width/2).*ones(numverts,3);
        case 3
            R(:,1)  = (width(1)/2).*ones(numverts,1);
            R(:,2)  = (width(2)/2).*ones(numverts,1);
            R(:,3)  = (width(3)/2).*ones(numverts,1);
        case numverts
            width(zindex) = [];
            R(:,1)  = width./2;
            R(:,2)  = width./2;
            R(:,3)  = width./2;
        case numverts*3
            width(zindex) = [];
            R  = width./2;
        otherwise
            error('Dimension of radius should be (1,1), (1,3), (numverts,1), (numverts,3)');
    end
    
    d1(end+1,:) = d1(end,:);
    
    x10 = verts(:,1)';
    x20 = verts(:,2)';
    x30 = verts(:,3)';
    x11 = d1(:,1)';
    x21 = d1(:,2)';
    x31 = d1(:,3)';
    
    a = verts(2,:) - verts(1,:);
    b = [0 0 1];
    c = crossSimple(a,b);
    if ~any(c)
        b = [1 0 0];
        c = crossSimple(a,b);
    end
    b = crossSimple(c,a);
    normb = norm(b); if normb~=0, b = b/norm(b); end
    %b = b*R(1);
    
    unitnormals(1,:) = b;
    
    for j = 1:numverts-1;
        
        a = verts(j+1,:)-verts(j,:);
        c = crossSimple(a,b);
        b = crossSimple(c,a);
        normb = norm(b); if normb~=0, b = b/norm(b); end
        %  b = b*R(j);
        unitnormals(j+1,:) = b;
        
    end
    
    unitnormal1 = unitnormals(:,1)';
    unitnormal2 = unitnormals(:,2)';
    unitnormal3 = unitnormals(:,3)';
    
    speed = sqrt(x11.^2 + x21.^2 + x31.^2);
    
    % And the binormal vector ( B = T x N )
    binormal1 = (x21.*unitnormal3 - x31.*unitnormal2) ./ speed;
    binormal2 = (x31.*unitnormal1 - x11.*unitnormal3) ./ speed;
    binormal3 = (x11.*unitnormal2 - x21.*unitnormal1) ./ speed;
    
    % s is the coordinate along the circular cross-sections of the tube:
    s = (0:n)';
    s = (2*pi/n)*s;
    
    % Finally, the parametric surface.
    % Each of x1, x2, x3 is an (m+1)x(n+1) matrix.
    % Rows represent coordinates along the tube.  Columns represent coordinates
    % sgcfin each (circular) cross-section of the tube.
    
    xa1 = ones(n+1,1)*x10;
    xb1 = (cos(s)*unitnormal1 + sin(s)*binormal1);
    xa2 = ones(n+1,1)*x20;
    xb2 = (cos(s)*unitnormal2 + sin(s)*binormal2);
    xa3 = ones(n+1,1)*x30;
    xb3 = (cos(s)*unitnormal3 + sin(s)*binormal3);
    
    x1 = xa1 + repmat(R(:,1)',[n+1,1]).*xb1;
    x2 = xa2 + repmat(R(:,2)',[n+1,1]).*xb2;
    x3 = xa3 + repmat(R(:,3)',[n+1,1]).*xb3;
    
    x = x1';
    y = x2';
    z = x3';   
end
% simple cross product
function c=crossSimple(a,b)
    c(1) = b(3)*a(2) - b(2)*a(3);
    c(2) = b(1)*a(3) - b(3)*a(1);
    c(3) = b(2)*a(1) - b(1)*a(2);
end
