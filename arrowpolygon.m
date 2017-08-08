function [xout, yout] = arrowpolygon(x, y, width, headwidth, headlength, plotflag, usefast)
%ARROWPOLYGON Creates an arrow polygon along a specified line
%
% [xout, yout] = arrowpolygon(x, y, width, headwidth, headlength)
% [xout, yout] = arrowpolygon(x, y, width, headwidth, headlength, plotflag, usefast)
%
% This function returns an arrow-shaped polygon along a specified input
% line, which can be easily plotted.  It uses coordinate space rather than
% pixel space for its calculations, so it looks best when used on an axis
% with a 1:1 aspect ratio.  
%
% Input variables:
%
%   x:          x coordinates of line, given from tail to head
%
%   y:          y coordinates of line, given from tail to head
%
%   width:      width of the arrow
%
%   headwidth:  width of the base of the triangular arrow head
%
%   headlength: height from base to point of the triangular arrow head
%
%   plotflag:   true to turn on debugging plots.  Optional, Default is
%               false
%
%   usefast:    true to use faster (and less safe) version of polybool.
%               Optional, default is false.
%
% Output variables:
%
%   xout:       x coordinates of arrow polygon
%
%   yout:       y corrdinated of arrow polygon

% Copyright 2007 Kelly Kearney

%----------------------------
% Check input
%----------------------------

narginchk(5,7);

if nargin == 5
    plotflag = false;
    usefast = false;
elseif nargin == 6
   usefast = false;
end

x = x(:);
y = y(:);

npoint = length(x);
nseg = npoint - 1;

%----------------------------
% Arrow shaft
%----------------------------

% Circles around vertices (except first and last)

ncirc = 100;
theta = linspace(0,2*pi,ncirc);
theta = fliplr(theta); % so clockwise for polybool
xcirc = width/2 .* cos(theta);
ycirc = width/2 .* sin(theta);

[xgrid, xcircgrid] = meshgrid(x(2:end-1), xcirc);
xvert = xgrid + xcircgrid;
[ygrid, ycircgrid] = meshgrid(y(2:end-1), ycirc);
yvert = ygrid + ycircgrid;

% Rectangles around line segments

x1 = x(1:end-1);
x2 = x(2:end);
y1 = y(1:end-1);
y2 = y(2:end);

heading = atan(abs(y2-y1)./abs(x2-x1));
quad2 = (x2-x1) <  0 & (y2-y1) >= 0;
quad3 = (x2-x1) <  0 & (y2-y1) <  0;
quad4 = (x2-x1) >= 0 & (y2-y1) <  0;
heading(quad2) = pi - heading(quad2);
heading(quad3) = pi + heading(quad3);
heading(quad4) = 2*pi - heading(quad4);

xBottom = x(end) + headlength * 1.01 * cos(heading(end)+pi);
yBottom = y(end) + headlength * 1.01 * sin(heading(end)+pi);


headPerp1 = heading + pi/2;
headPerp2 = heading - pi/2;
xcorner1 = x1 + width/2*cos(headPerp1);
xcorner2 = x2 + width/2*cos(headPerp1);
xcorner3 = x1 - width/2*cos(headPerp1);
xcorner4 = x2 - width/2*cos(headPerp1);
ycorner1 = y1 + width/2*sin(headPerp1);
ycorner2 = y2 + width/2*sin(headPerp1);
ycorner3 = y1 - width/2*sin(headPerp1);
ycorner4 = y2 - width/2*sin(headPerp1);

xrect = nan(5, nseg);
yrect = nan(5, nseg);
for iseg = 1:nseg
    xrect(:,iseg) = [xcorner1(iseg) xcorner2(iseg) xcorner4(iseg) xcorner3(iseg) xcorner1(iseg)];
    yrect(:,iseg) = [ycorner1(iseg) ycorner2(iseg) ycorner4(iseg) ycorner3(iseg) ycorner1(iseg)];
end

%----------------------------
% Arrow head
%----------------------------

xBot1 = xBottom + headwidth/2 * cos(headPerp1(end));
xBot2 = xBottom + headwidth/2 * cos(headPerp2(end));
yBot1 = yBottom + headwidth/2 * sin(headPerp1(end));
yBot2 = yBottom + headwidth/2 * sin(headPerp2(end));

xhead = [xBot1 xBot2 x(end) xBot1];
yhead = [yBot1 yBot2 y(end) yBot1];

%----------------------------
% Erasure boxes
%----------------------------

% Near head

xherase1 = x(end) + headlength * cos(heading(end)+pi) * .999;
xherase2 = x(end) + width      * cos(heading(end));
yherase1 = y(end) + headlength * sin(heading(end)+pi) * .999;
yherase2 = y(end) + width      * sin(heading(end));

xhecorner1 = xherase1 + headwidth/2*1.00*cos(headPerp1(end));
xhecorner2 = xherase2 + headwidth/2*1.00*cos(headPerp1(end));
xhecorner3 = xherase1 - headwidth/2*1.00*cos(headPerp1(end));
xhecorner4 = xherase2 - headwidth/2*1.00*cos(headPerp1(end));
yhecorner1 = yherase1 + headwidth/2*1.00*sin(headPerp1(end));
yhecorner2 = yherase2 + headwidth/2*1.00*sin(headPerp1(end));
yhecorner3 = yherase1 - headwidth/2*1.00*sin(headPerp1(end));
yhecorner4 = yherase2 - headwidth/2*1.00*sin(headPerp1(end));

xherase = [xhecorner1 xhecorner2 xhecorner4 xhecorner3 xhecorner1];
yherase = [yhecorner1 yhecorner2 yhecorner4 yhecorner3 yhecorner1];

% Near tail

xterase1 = x(1) + width * cos(heading(1)+pi);
xterase2 = x(1);
yterase1 = y(1) + width * sin(heading(1)+pi);
yterase2 = y(1);

xtecorner1 = xterase1 + headwidth/2*1.00*cos(headPerp1(1));
xtecorner2 = xterase2 + headwidth/2*1.00*cos(headPerp1(1));
xtecorner3 = xterase1 - headwidth/2*1.00*cos(headPerp1(1));
xtecorner4 = xterase2 - headwidth/2*1.00*cos(headPerp1(1));
ytecorner1 = yterase1 + headwidth/2*1.00*sin(headPerp1(1));
ytecorner2 = yterase2 + headwidth/2*1.00*sin(headPerp1(1));
ytecorner3 = yterase1 - headwidth/2*1.00*sin(headPerp1(1));
ytecorner4 = yterase2 - headwidth/2*1.00*sin(headPerp1(1));

xterase = [xtecorner1 xtecorner2 xtecorner4 xtecorner3 xtecorner1];
yterase = [ytecorner1 ytecorner2 ytecorner4 ytecorner3 ytecorner1];

%----------------------------
% Debug plotting stuff
%----------------------------

if plotflag
    figure;
    axes; hold on;
    plot(x, y, '-kx');
    plot(xrect, yrect, 'r');
    plot(xvert, yvert, 'b');
    plot(xhead, yhead, 'g');
    plot(xherase, yherase, 'y');
    plot(xterase, yterase, 'y');
    set(gca, 'DataAspectRatio', [1 1 1]);
end

%----------------------------
% Merge polygons
%----------------------------

% Convert circles and rectangles to clockwise cell arrays for polybool

xrect = mat2cell(xrect, 5, ones(1,size(xrect,2)));
yrect = mat2cell(yrect, 5, ones(1,size(yrect,2)));
if npoint > 2
    xvert = mat2cell(xvert, ncirc, ones(1,size(xvert,2)));
    yvert = mat2cell(yvert, ncirc, ones(1,size(yvert,2)));
end

[xhead, yhead] = poly2cw(xhead, yhead);
[xrect, yrect] = poly2cw(xrect, yrect);
[xvert, yvert] = poly2cw(xvert, yvert);

% Merge circles and rectangles to form arrow shaft

xout = NaN;
yout = NaN;

if npoint > 2
    for ipoly = 1:size(xvert,2)
        if ~usefast
            [xout, yout] = polybool('union', xout, yout, xvert{ipoly}, yvert{ipoly});
        else
            [xout, yout] = polyboolfast('union', xout, yout, xvert{ipoly}, yvert{ipoly});
        end
    end
end

for ipoly = 1:size(xrect,2)
    if ~usefast
        [xout, yout] = polybool('union', xout, yout, xrect{ipoly}, yrect{ipoly});
    else
        [xout, yout] = polyboolfast('union', xout, yout, xrect{ipoly}, yrect{ipoly});
    end
end

% Erase any part of the shaft that would interfere with arrowhead

if ~usefast
    [xout, yout] = polybool('-', xout, yout, xherase, yherase);
else
    [xout, yout] = polyboolfast('diff', xout, yout, xherase, yherase);
end
    
if isempty(xout)
    error('The dimesions of your arrow are causing problems; try a narrower width or smaller head length');
end

if ~usefast
    [xout, yout] = polybool('-', xout, yout, xterase, yterase);
else
    [xout, yout] = polyboolfast('diff', xout, yout, xterase, yterase);
end
    
if isempty(xout)
    error('The dimesions of your arrow are causing problems; try a narrower width or smaller head length');
end

if plotflag
    plot(xout, yout, 'Color', [.8 .8 .8], 'LineWidth', 2, 'LineStyle', '--');
end

% Merge arrowhead

if ~usefast
    [xout, yout] = polybool('union', xout, yout, xhead, yhead);
else
    [xout, yout] = polyboolfast('union', xout, yout, xhead, yhead);
end
    
if plotflag
    plot(xout, yout, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
end

%----------------------------
% Subfunction: polyboolfast
% polybool without error 
% checks and assuming all 
% clockwise contours
%----------------------------

function [xout, yout] = polyboolfast(flag, x1, y1, x2, y2)

mappath = fullfile(matlabroot, 'toolbox', 'map', 'map', 'private');
Gpc = dir(fullfile(mappath, 'gpcmex*'));
gpcmexpath = fullfile(mappath, Gpc.name);

gpcmex = function_handle(gpcmexpath);

if strcmp(flag, 'diff')
    op = '-';
else
    op = flag;
end

try

    [xout, yout, emp] = handleEmptyInputs(op, x1, y1, x2, y2);

    if ~emp
    
        [x1, y1] = polysplit(x1, y1);
        [x2, y2] = polysplit(x2, y2);


        P1 = struct('x', x1, 'y', y1, 'ishole', false);
        P2 = struct('x', x2, 'y', y2, 'ishole', false);

        Out = gpcmex(flag, P1, P2);

        [xout, yout] = polyjoin({Out.x}, {Out.y});   
    end
    
catch ME
    
    if strcmp(ME.identifier, 'map:gpcmex:invalidPolygonArray')
        [xout, yout] = polybool(op, x1, y1, x2, y2);
    else
        rethrow(ME);
    end

end

%----------------------------
% Subfunction: stolen from 
% polybool
%----------------------------

function [x3, y3, emptyInput] = handleEmptyInputs(operation, x1, y1, x2, y2)
% Assuming that x1 and y1, and x2 and y2, are pairs of inputs having
% consistent sizes, return the appropriate values for x3 and y3 in the
% event that x1 and/or x2 are empty (or contain only NaN), and set
% emptyInput to true. Otherwise, set x3 and y3 to empty and set
% emptyInput to false. Operation has been validated and equals one of
% the following strings: 'int','union','xor','diff'.

% NaN-only arrays should behave the same way as empty arrays, so filter
% them up-front. Be careful, because all(isnan([])) evaluates to true.
% Also, be careful to preserve shape: return 1-by-0 given a row
% vector of NaN and a 0-by-1 given a column vector of NaN.
if  all(isnan(x1)) && ~isempty(x1)
    x1(1:end) = [];
    y1(1:end) = [];
end

if all(isnan(x2)) && ~isempty(x2)
    x2(1:end) = [];
    y2(1:end) = [];
end

if isempty(x2)
    if strcmp(operation,'int')
        % Intersection is empty, but preserve shape
        % by using x2 and y2 rather than [].
        x3 = x2;
        y3 = y2;
    else
        % Union, exclusive or, or difference with
        % empty leaves x1 and y1 unaltered.
        x3 = x1;
        y3 = y1;
    end
    emptyInput = true;
elseif isempty(x1)
    if any(strcmp(operation,{'int','diff'}))
        % Intersection or difference is empty, but preserve
        % shape by using x1 and y1 rather than [].
        x3 = x1;
        y3 = y1;        
    else
        % Union or exclusive or with empty leaves x2 and y2 unaltered.
        x3 = x2;
        y3 = y2;
    end
    emptyInput = true;
else
    x3 = [];
    y3 = [];
    emptyInput = false;
end





