function h = gplotpoly(adj, xy, nval, varargin)
%GPLOTPOLY Graph with polygon arrows
%
% h = gplotpoly(adj, xy, nval, p1, v1)
%
% Input variables:
%
%   adj:    n x n adjacency matrix, with non-zero values corresponding to
%           relative width of arrow
%
%   xy:     n x 2 array of x and y coordinates for the nodes
%
%   nval:   n x 1 array of values used to scale the size of each node
%
% Optional input values:
%
%   edgetype:   indicates how to generate edge coordinates ['curve']
%               'line':     straight lines
%               'curve':    curves counterclockwise from source to sink, as
%                           in gplotdc
%               'manual':   edge coordinates are passed via 'edgecoord'
%
%   edgecoord:  nnz x 1 cell array of two-column arrays, x and y
%               coordinates corresponding to each edge, in order of
%               find(adj) []
%
%   edgecurve:  for 'curve' edge type, degree of curvature [0.04]
%
%   edgeval:    nnz x 1 array, values corresponding to non-zero adj values.
%               Can be useful to use this if you require 0-width edges []
%
%   wlim:       arrow widths that correspond to minimum and maximum edge
%               data values, in axis coordinates [0.01 0.05]
%
%   elim:       data values corresponding to narrowest and widest arrows.
%               By default, this will be the minimum and maximum values in
%               adj
%
%   rlim:       node radius values corresponding to node data limits
%               [0.01 0.05]
%
%   nlim:       data values corresponding to smallest and largest nodes.
%               By default, this will be the minimum and maximum values in
%               nval
%
%   ecdata:     values used to colormap the arrow edges. [adj(adj~=0)]
%
%   ncdata:     values used to colormap the node edges [nval]
%
%   axis:       axis handle to plot to [gca]

% Copyright 2014 Kelly Kearney

if ~isvector(nval)
    error('nval must be vector');
end
nval = reshape(nval, 1, []);

Opt.edgetype = 'curve';
Opt.edgecoord = [];
Opt.edgecurve = 0.04;
Opt.wlim = [0.01 0.05];
Opt.elim = minmax(adj(adj~=0));
Opt.rlim = [0.01 0.05];
Opt.nlim = minmax(nval);
Opt.ecdata = adj(adj~=0);
Opt.ncdata = nval;
Opt.axis = gca;
Opt.arrow = false;
Opt.edgeval = [];

Opt = parsepv(Opt, varargin);

if ~isvector(Opt.ncdata)
    error('ncdata must be vector');
end
Opt.ncdata = reshape(Opt.ncdata, 1, []);

if Opt.elim(1) == Opt.elim(2)
    Opt.elim = Opt.elim + Opt.elim(1)*0.1*[-1 1];
end
if Opt.nlim(1) == Opt.nlim(2)
    Opt.nlim = Opt.nlim + Opt.nlim(1)*0.1*[-1 1];
end

% Construct path for edges

[ii,jj] = find(adj);

switch Opt.edgetype
    case 'line'
        X = [ xy(ii,1) xy(jj,1)]'; % TODO: arrow needs more points
        Y = [ xy(ii,2) xy(jj,2)]';
        xl = num2cell(X,1)';
        yl = num2cell(Y,1)';
    case 'curve'
        X = [ xy(ii,1) xy(jj,1)]';
        Y = [ xy(ii,2) xy(jj,2)]';

        [X,Y] = makeCurved(X,Y,Opt.edgecurve);
        xl = num2cell(X,1)';
        yl = num2cell(Y,1)';
    case 'manual'
        xl = cellfun(@(x) x(:,1), Opt.edgecoord, 'uni', 0);
        yl = cellfun(@(x) x(:,2), Opt.edgecoord, 'uni', 0);
%         error('Work in progress');
        
end

% Arrow coordinates

if isempty(Opt.edgeval)
    erel = adj(adj ~= 0);
else
    erel = Opt.edgeval;
end
w = interp1(Opt.elim, Opt.wlim, erel);
w(erel > Opt.elim(2)) = Opt.wlim(2);
w(erel < Opt.elim(1)) = Opt.wlim(1);

hwid = w*1.5; 
hlen = w;

if Opt.arrow
    [xa, ya] = cellfun(@(xl,yl,w,hw,hl) arrowpolygon(xl,yl,w,hw,hl,false,true), ...
                       xl, yl, num2cell(w), num2cell(hwid), num2cell(hlen), 'uni', 0);
else
    [xa, ya, ca] = cellfun(@(xl,yl,w) line2poly(xl, yl, w), xl, yl, num2cell(w), 'uni', 0);
end


% len = cellfun(@length, xa);
% npt = max(len);
% [xa2, ya2] = deal(zeros(npt,length(xa)));
% for ie = 1:size(xa2,2)
%     xa2(1:len(ie),ie) = xa{ie};
%     xa2(len(ie)+1:npt) = xa{ie}(end);
%     
%     ya2(1:len(ie),ie) = ya{ie};
%     ya2(len(ie)+1:npt) = ya{ie}(end);
%    
% end

% Node coordinates

th = linspace(0,2*pi,21)';

r = interp1(Opt.nlim, Opt.rlim, nval);
r(nval > Opt.nlim(2)) = Opt.rlim(2);
r(nval < Opt.nlim(1)) = Opt.rlim(1);

xn = bsxfun(@plus, cos(th) * r, xy(:,1)');
yn = bsxfun(@plus, sin(th) * r, xy(:,2)');

h.arrow = cellfun(@(x,y,c) patch(x,y,c, 'parent', Opt.axis), xa, ya, num2cell(Opt.ecdata));
h.node = patch(xn, yn, Opt.ncdata, 'parent', Opt.axis);



% if isscalar(Opt.edgeline)
%     figure;
%     h = gplotcolor(adj, xy, zeros(nnz(adj),3));
%     xyline = get(h, {'xdata', 'ydata'});
%     close(gcf);
% end
% 
% % Plot
% 
% flim = minmax(cat(1, flxint{:}));
% wlim = [0.01 0.05];
% 
% w = interp1(log10(flim), wlim, log10(adj(adj>0)));
% hwid = w*1.5; 
% hlen = w;
% 
% [xa, ya] = cellfun(@(xl,yl,w,hw,hl) arrowpolygon(xl,yl,w,hw,hl), ...
%     xyline(:,1), xyline(:,2), num2cell(w), num2cell(hwid), num2cell(hlen), 'uni', 0);
% 
% th = linspace(0,2*pi,13)';
% blim = minmax(Bi(1).bio(:,imin,1:27));
% rlim = [0.01 0.05];
% 
% r = interp1(log10(blim), rlim, log10(bval));
% xn = bsxfun(@plus, cos(th) * r, xy(:,1)');
% yn = bsxfun(@plus, sin(th) * r, xy(:,2)');
% 
% h = plotgrid('setup', cell(1));
% hold on;
% c = find(adj);
% h.arrow = cellfun(@(x,y,c) patch(x,y,c), xa, ya, num2cell(w));
% axis equal;
% h.node = patch(xn, yn, r);


function [xc,yc] = makeCurved(x,y,pct)

nx = size(x,2);
[xc,yc] = deal(zeros(17,nx));

for ix = 1:nx
    xtmp = x(:,ix);
    ytmp = y(:,ix);
    pctn = pct;
    for k = 1:4
        n = length(xtmp);
        m = 2*n - 1;
        xtmp2 = zeros(1,m);
        ytmp2 = zeros(1,m);
        xtmp2(1:2:m) = xtmp;
        ytmp2(1:2:m) = ytmp;
        xtmp2(2:2:m-1) = 0.5*(xtmp(1:n-1)+xtmp(2:n))+pctn*diff(ytmp);
        ytmp2(2:2:m-1) = 0.5*(ytmp(1:n-1)+ytmp(2:n))-pctn*diff(xtmp);
        xtmp = xtmp2;
        ytmp = ytmp2;
        pctn = pctn * 0.5;
    end
    xc(:,ix) = xtmp;
    yc(:,ix) = ytmp;
    
end