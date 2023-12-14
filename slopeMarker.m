function slopeMarker(p, xid, xlen, yval, mref, k, close, slide, style, varargin)
  xref = mref^2 * k;
  xstep = 4;
  x1 = xstep^(xid-1) * xref;
  x2 = x1 * xstep^xlen;
  y1 = yval;
  % y = C*(x/k)^(-p/2)
  y2 = y1 * (x2/x1)^(-p/2);
  if close == 1
    x = [x1,x2,x2,x1];
    y = [y1,y2,y1,y1];
  elseif close == 2
    x = [x1,x2,x1,x1];
    y = [y1,y2,y2,y1];
  else
    x = [x1,x2];
    y = [y1,y2];
  end
  loglog(x, y, style, 'LineWidth', 1, varargin{:});
  xt = sqrt(x1^(1-slide)*x2^(1+slide));
  yt = sqrt(y1^(1-slide)*y2^(1+slide));
  % Nudge the text a little to the side, normal to the line
  %slope = log(y2/y1) / log(x2/x1)
  xL = xlim; yL = ylim;
  axpos = get(gca, 'Position');
  figpos = get(gcf, 'Position');
  W = axpos(3) * figpos(3);
  H = axpos(4) * figpos(4);
  A = W/log(xL(2)/xL(1));
  B = H/log(yL(2)/yL(1));
  a = @(x) A*log(x/xL(1));
  b = @(y) B*log(y/yL(1));
  at = a(xt); bt = b(yt);
  a1 = a(x1); b1 = b(y1);
  a2 = a(x2); b2 = b(y2);
  n = [b1-b2, a2-a1];
  n = n/norm(n);
  if close == 1 || close == 0
    dist = 25; % pixels
    at = at + dist*n(1); bt = bt + dist*n(2);
  elseif close == 2 || close == -1
    dist = 20; % pixels
    at = at - dist*n(1); bt = bt - dist*n(2);
  end
  ainv = @(a) xL(1)*exp(a/A);
  binv = @(b) yL(1)*exp(b/B);
  xt = ainv(at);
  yt = binv(bt);
  text(xt, yt, sprintf('O(h^{%d})', p), 'HorizontalAlignment', 'center', varargin{:});
end
