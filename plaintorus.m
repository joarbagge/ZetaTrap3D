function s = plaintorus(a)
% Geometry generation following Barnett's idea in BIE3D wobblytorus.
% Generate a plain torus by revolving a generating curve about z-axis.
%
% Define generating curve [x,z] = [R,0] + a*[cos(th),sin(th)], rotate about z-axis
%
% Input:
%   a = minor radius (tube radius)
% Output:
%   s = struct containing geom info
%   s.Z =  [x(u,v); y(u,v); z(u,v)] func handle for surface parameterization
%   s.Zu = [xu(u,v); yu(u,v); zu(u,v)] func handle for the partials w.r.t. u
%   s.Zv = [xv(u,v); yv(u,v); zv(u,v)] func handle for the partials w.r.t. v
%       All of above assume u and v are row vectors.
% 
% Bowei Wu 4/28/20
% Added function "get_x_xp_xpp_xppp" computing sample points and
% derivatives more efficiently

if nargin == 0 && nargout == 0, test_wobblytorus; return; end
if nargin < 1, a = 0.2; end
m = 0; n = 0; k = 1;

% planar curve
ps = 0; % phase shift
ts = 0; % t shift
f = @(p,t) a * cos(m*t+n/k*p+ps); %generating curve
fp= @(p,t) -n/k*a*sin(m*t+n/k*p+ps); %partials
ft= @(p,t) -m*a*sin(m*t+n/k*p+ps);
fpp= @(p,t) -(n/k)^2*a*cos(m*t+n/k*p+ps); %2nd partials
ftt= @(p,t) -m^2*a*cos(m*t+n/k*p+ps);
fpt= @(p,t) -n/k*m*a*cos(m*t+n/k*p+ps);
fppp=@(p,t) n^3*a*sin(m*t+n*p+ps);  %3rd partials
fppt=@(p,t) n^2*m*a*sin(m*t+n*p+ps);
fptt=@(p,t) n*m^2*a*sin(m*t+n*p+ps);
fttt=@(p,t) m^3*a*sin(m*t+n*p+ps);
fpppp=@(p,t) n^4*a*cos(m*t+n*p+ps); %4th partials
fpppt=@(p,t) m*n^3*a*cos(m*t+n*p+ps);
fpptt=@(p,t) m^2*n^2*a*cos(m*t+n*p+ps);
fpttt=@(p,t) m^3*n*a*cos(m*t+n*p+ps);
ftttt=@(p,t) m^4*a*cos(m*t+n*p+ps);

%t = linspace(0,2*pi,100); plot(f(0,t).*cos(t),f(0,t).*sin(t)); axis equal %plot it

% toroidal surface
R = 1; r = 1;
x = @(p,t) (R + r * f(p,t) .* cos(t+ts)) .* cos(p/k);
y = @(p,t) (R + r * f(p,t) .* cos(t+ts)) .* sin(p/k);
z = @(p,t) r * f(p,t) .* sin(t+ts);

% derivatives
xp = @(p,t) -R/k * sin(p/k) + r * (fp(p,t) .* cos(p/k) - 1/k*f(p,t) .* sin(p/k)) .* cos(t+ts);
yp = @(p,t)  R/k * cos(p/k) + r * (fp(p,t) .* sin(p/k) + 1/k*f(p,t) .* cos(p/k)) .* cos(t+ts);
zp = @(p,t) r * fp(p,t) .* sin(t+ts);

xt = @(p,t) r * (ft(p,t) .* cos(t+ts) - f(p,t) .* sin(t+ts)) .* cos(p/k);
yt = @(p,t) r * (ft(p,t) .* cos(t+ts) - f(p,t) .* sin(t+ts)) .* sin(p/k);
zt = @(p,t) r * (ft(p,t) .* sin(t+ts) + f(p,t) .* cos(t+ts));

% 2nd derivatives
xpp = @(p,t) -R/k^2 * cos(p/k) + r * (fpp(p,t) .* cos(p/k) - 2/k * fp(p,t) .* sin(p/k) - 1/k^2*f(p,t) .* cos(p/k)) .* cos(t+ts);
ypp = @(p,t) -R/k^2 * sin(p/k) + r * (fpp(p,t) .* sin(p/k) + 2/k * fp(p,t) .* cos(p/k) - 1/k^2*f(p,t) .* sin(p/k)) .* cos(t+ts);
zpp = @(p,t) r * fpp(p,t) .* sin(t+ts);

xtt = @(p,t) r * (ftt(p,t) .* cos(t+ts) - 2 * ft(p,t) .* sin(t+ts) - f(p,t) .* cos(t+ts)) .* cos(p/k);
ytt = @(p,t) r * (ftt(p,t) .* cos(t+ts) - 2 * ft(p,t) .* sin(t+ts) - f(p,t) .* cos(t+ts)) .* sin(p/k);
ztt = @(p,t) r * (ftt(p,t) .* sin(t+ts) + 2 * ft(p,t) .* cos(t+ts) - f(p,t) .* sin(t+ts));

xpt = @(p,t) r * ( (fpt(p,t) .* cos(p/k) - 1/k*ft(p,t) .* sin(p/k)) .* cos(t+ts) - (fp(p,t) .* cos(p/k) - 1/k*f(p,t) .* sin(p/k)) .* sin(t+ts) );
ypt = @(p,t) r * ( (fpt(p,t) .* sin(p/k) + 1/k*ft(p,t) .* cos(p/k)) .* cos(t+ts) - (fp(p,t) .* sin(p/k) + 1/k*f(p,t) .* cos(p/k)) .* sin(t+ts) );
zpt = @(p,t) r * ( fpt(p,t) .* sin(t+ts) + fp(p,t) .* cos(t+ts) );

% output
s.Z  = @(u,v) [x(u,v); y(u,v); z(u,v)];
s.Zu = @(u,v) [xp(u,v); yp(u,v); zp(u,v)];
s.Zv = @(u,v) [xt(u,v); yt(u,v); zt(u,v)];
s.Zuu= @(u,v) [xpp(u,v);ypp(u,v);zpp(u,v)];
s.Zuv= @(u,v) [xpt(u,v);ypt(u,v);zpt(u,v)];
s.Zvv= @(u,v) [xtt(u,v);ytt(u,v);ztt(u,v)];
s.topo = 't';
s.type = 'periodic'; s.L = 0;

% faster output
s.f = f; %generating curve
s.fu = fp; %partials
s.fv = ft;
s.fuu = fpp; %second partials
s.fuv = fpt;
s.fvv = ftt;
s.fuuu = fppp; %third partials
s.fuuv = fppt;
s.fuvv = fptt;
s.fvvv = fttt;
s.fuuuu = fpppp; %fourth partials
s.fuuuv = fpppt;
s.fuuvv = fpptt;
s.fuvvv = fpttt;
s.fvvvv = ftttt;
s.fparam = [R,r,ts,m,n];
s.get_x_xp_xpp_xppp = @(u,v) get_x_xp_xpp_xppp(s,u,v); % compute X, X', X'' efficiently

function [X,Xu,Xv,Xuu,Xuv,Xvv,Xuuu,Xuuv,Xuvv,Xvvv,...
    Xuuuu,Xuuuv,Xuuvv,Xuvvv,Xvvvv] = get_x_xp_xpp_xppp(s,u,v)
% function computing X, X', X'', X''', X'''' efficiently

% get parameters
R = s.fparam(1); r = s.fparam(2); vs = s.fparam(3);
u = u(:).'; v = v(:).';

% prepare functions
f = s.f(u,v);
su = sin(u); sv = sin(v+vs);
cu = cos(u); cv = cos(v+vs);

% toroidal surface
rr = (R + r * f .* cv);
X = [rr.* cu;
     rr.* su;
      r * f .* sv];

% derivatives
if nargin > 1
    fv = s.fv(u,v);
    fu = fv;
    ru = r * fu .* cv;
    Xu = [-rr .* su + ru .* cu;
           rr .* cu + ru .* su;
           r * fu .* sv];
    rv = r * (fv .* cv - f .* sv);
    Xv = [rv .* cu;
          rv .* su;
          r * (fv .* sv + f .* cv)];
end

% 2nd derivatives
if nargout > 3
    fvv = s.fvv(u,v);
    fuv = fvv;
    fuu = fuv;
    ruu = r * fuu .* cv;
    Xuu = [-rr .* cu - 2 * ru .* su + ruu .* cu;
           -rr .* su + 2 * ru .* cu + ruu .* su;
            r * fuu .* sv];
    ruv = r * (fuv .* cv - fu .* sv);
    Xuv = [ruv .* cu - rv .* su;
           ruv .* su + rv .* cu;
           r * ( fuv .* sv + fu .* cv )];
    rvv = r * (fvv .* cv - 2 * fv .* sv - f .* cv);
    Xvv = [rvv .* cu;
           rvv .* su;
           r * (fvv .* sv + 2 * fv .* cv - f .* sv)];
end

% 3rd derivatives
if nargout > 6
    fvvv = s.fvvv(u,v);
    fuvv = fvvv;
    fuuv = fuvv;
    fuuu = fuuv;
    ruuu = r * fuuu .* cv;
    Xuuu = [ rr .* su - 3 * ru .* cu - 3 * ruu .* su + ruuu .* cu;
            -rr .* cu - 3 * ru .* su + 3 * ruu .* cu + ruuu .* su;
             r * fuuu .* sv];
    ruuv = r * (fuuv .* cv - fuu .*sv);
    Xuuv = [ruuv .* cu - 2 * ruv .* su - rv .* cu;
            ruuv .* su + 2 * ruv .* cu - rv .* su;
            r * (fuuv .* sv + fuu .* cv)];
    ruvv = r * (fuvv .* cv - 2 * fuv .* sv - fu .* cv);
    Xuvv = [ruvv .* cu - rvv .* su;
            ruvv .* su + rvv .* cu;
            r * (fuvv .* sv + 2 * fuv .* cv - fu .* sv)];
    rvvv = r * (fvvv .* cv - 3 * fvv .* sv - 3 * fv .* cv + f .* sv);
    Xvvv = [rvvv .* cu;
            rvvv .* su;
            r * (fvvv .* sv + 3 * fvv .* cv - 3 * fv .* sv - f .* cv)];
end

% 4th derivatives
if nargout > 10
    fvvvv = s.fvvvv(u,v);
    fuvvv = fvvvv;
    fuuvv = fuvvv;
    fuuuv = fuuvv;
    fuuuu = fuuuv;
    ruuuu = r * fuuuu .* cv;
    Xuuuu = [rr .* cu + 4 * ru .* su - 6 * ruu .* cu - 4 * ruuu .* su + ruuuu .* cu;
             rr .* su - 4 * ru .* cu - 6 * ruu .* su + 4 * ruuu .* cu + ruuuu .* su;
             r * fuuuu .* sv];
    ruuuv = r * (fuuuv .* cv - fuuu .*sv);
    Xuuuv = [ rv .* su - 3 * ruv .* cu - 3 * ruuv .* su + ruuuv .* cu;
             -rv .* cu - 3 * ruv .* su + 3 * ruuv .* cu + ruuuv .* su;
             r * (fuuuv .* sv + fuuu .* cv)];
    ruuvv = r * (fuuvv .* cv - 2 * fuuv .* sv - fuu .* cv);
    Xuuvv = [ruuvv .* cu - 2 * ruvv .* su - rvv .* cu;
             ruuvv .* su + 2 * ruvv .* cu - rvv .* su;
             r * (fuuvv .* sv + 2 * fuuv .* cv - fuu .* sv)];
    ruvvv = r * (fuvvv .* cv - 3 * fuvv .* sv - 3 * fuv .* cv + fu .* sv);
    Xuvvv = [ruvvv .* cu - rvvv .* su;
             ruvvv .* su + rvvv .* cu;
             r * (fuvvv .* sv + 3 * fuvv .* cv - 3 * fuv .* sv - fu .* cv)];   
    rvvvv = r * (fvvvv .* cv - 4 * fvvv .* sv - 6 * fvv .* cv + 4 * fv .* sv + f .* cv);
    Xvvvv = [rvvvv .* cu;
             rvvvv .* su;
             r * (fvvvv .* sv + 4 * fvvv .* cv - 6 * fvv .* sv - 4 * fv .* cv + f .* sv)];
end

function test_wobblytorus

k = 3;
s = wobblytorus2(3,2,0.2,k);

n = 30;
p = n*k+1; q = n+1;
[u,v] = meshgrid(linspace(0,2*k*pi,p),linspace(0,2*pi,q));
u = u(:)'; v = v(:)'; % turn into row vectors

% 3-by-n arrays of points and partials and normals
X = s.Z(u,v);
Xu= s.Zu(u,v);
Xv= s.Zv(u,v);

NX = cross(Xu,Xv);
J = vecnorm(NX); %length of normal (jacobian)
NX = NX./J; % unit normal

% plot it
x = reshape(X(1,:),q,p);
y = reshape(X(2,:),q,p);
z = reshape(X(3,:),q,p);
nx = reshape(NX(1,:),q,p);
ny = reshape(NX(2,:),q,p);
nz = reshape(NX(3,:),q,p);
mesh(x,y,z,'FaceAlpha',0); axis equal; hold on 
quiver3(x,y,z,nx,ny,nz); axis equal; hold off 

