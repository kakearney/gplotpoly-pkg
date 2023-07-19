function [xp,yp,cp] = line2poly(x, y, w)
%LINE2POLY Create polygon based on a widened line
%
% [xp,yp,cp] = line2poly(x, y, w)
%
% Creates a polygon based on a line with a prescribed width orthoganal to
% that line.
%
% Input variables:
%
%   x:  vector, x coorinates of line
%
%   y:  vector, y coordinates of line
%
%   w:  width of line at each point in the line.  Can be either a vector
%       the same size as x and y, or a scalar (constant width for entire
%       line).
%
% Output variables:
%
%   xp: x coordinates of polygon patch
%
%   yp: y coordinates of polygon patch
%
%   cp: distance along centerline of the line, where (x(1),y(1)) = 0 and
%       (x(end,y(end)) = 1.  This can be used as a color value to create a
%       color gradient along the patch polygon from one end to the other.

% Copyright 2017 Kelly Kearney

x = x(:);
y = y(:);
if isscalar(w)
    w = ones(size(x)) * w;
else
    w = w(:);
end

nx = length(x);
p = [x y]';
w = w/2;

% Unit normal for each line segment

v = diff(p, [], 2);
N = zeros(2,nx-1);
for iv = 1:nx-1
    N(:,iv) = [0 -1; 1 0]*(v(:,iv)./norm(v(:,iv))); % Unit normal right at each vector
end

% Check direction at each turn

c = zeros(3,nx);
for iv = 2:nx-1
    c(:,iv) = cross([p(:,iv+1) - p(:,iv-1); 0], [p(:,iv)-p(:,iv-1); 0]);
end

% Distance along line at each point

dr = sqrt(v(1,:).^2 + v(2,:).^2);
dr = [0 cumsum(dr)];
dr = dr./max(dr);

% Calculate new points (o = out segment, i = in segment)

lo = p(:,1:end-1) + bsxfun(@times, N, w(1:end-1)');  
ro = p(:,1:end-1) - bsxfun(@times, N, w(1:end-1)');  
li = p(:,2:end) + bsxfun(@times, N, w(2:end)');    
ri = p(:,2:end) - bsxfun(@times, N, w(2:end)');

li = [nan(2,1) li];
ri = [nan(2,1) ri];
lo = [lo nan(2,1)];
ro = [ro nan(2,1)];

% Calculate where intersections/mitering needs to be added

[xl,yl,xr,yr,cl, cr] = deal(cell(nx,1));
xl{1} = lo(1,1);
yl{1} = lo(2,1);
xr{1} = ro(1,1);
yr{1} = ro(2,1);
cl{1} = dr(1);
cr{1} = dr(1);

xl{end} = li(1,end);
yl{end} = li(2,end);
xr{end} = ri(1,end);
yr{end} = ri(2,end);
cl{end} = dr(end);
cr{end} = dr(end);

for ii = 2:nx-1
    
    if c(3,ii) < 1e-10    % Straight, for all intents and purposes
        xr{ii} = ri(1,ii);
        yr{ii} = ri(2,ii);
        xl{ii} = li(1,ii);
        yl{ii} = li(2,ii);
        cr{ii} = dr(ii);
        cl{ii} = dr(ii);
    else
        
        x1 = [lo(1,ii-1) li(1,ii) ri(1,ii) ro(1,ii-1) lo(1,ii-1)];
        y1 = [lo(2,ii-1) li(2,ii) ri(2,ii) ro(2,ii-1) lo(2,ii-1)];
        x2 = [lo(1,ii) li(1,ii+1) ri(1,ii+1) ro(1,ii) lo(1,ii)];
        y2 = [lo(2,ii) li(2,ii+1) ri(2,ii+1) ro(2,ii) lo(2,ii)];

        [xint, yint] = polyxpoly(x1, y1, x2, y2);
        if length(xint) == 1
            warning('Issue with polygon intersection... check this');
        end
        [~,imax] = max((xint - x(ii)).^2 + (yint - y(ii)).^2);
    
        if c(3,ii) < 0 % Turns left
       
            xl{ii} = xint(imax);
            yl{ii} = yint(imax);

            xr{ii} = [ri(1,ii); ro(1,ii)];
            yr{ii} = [ri(2,ii); ro(2,ii)];

            cl{ii} = dr(ii);
            cr{ii} = [dr(ii); dr(ii)];

        else           % Turns right

            xr{ii} = xint(imax);
            yr{ii} = yint(imax);

            xl{ii} = [li(1,ii); lo(1,ii)];
            yl{ii} = [li(2,ii); lo(2,ii)];

            cr{ii} = dr(ii);
            cl{ii} = [dr(ii); dr(ii)];

        end
    end
end

% Concatenate left and right

xl = cat(1, xl{:});
yl = cat(1, yl{:});
xr = cat(1, xr{:});
yr = cat(1, yr{:});
cl = cat(1, cl{:});
cr = cat(1, cr{:});

xp = [xl; xr(end:-1:1); xl(1)];
yp = [yl; yr(end:-1:1); yl(1)];
cp = [cl; cr(end:-1:1); cl(1)];
