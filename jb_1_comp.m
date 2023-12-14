% Compute potential at reference grid points, for convergence test
% Supports: Laplace, wobbly/plain torus, exact derivatives only

orders = [3, 5];
geometries = {'plain', 'wobbly'};
densities = {'const', 'trig'};
Nvs = 8*[1 2 4 8];
plotting = true;
recompute = false;
version = 'v1';

for j1=1:numel(orders)
    ord = orders(j1);
    for j2=1:numel(geometries)
        geom = geometries{j2};
        for j3=1:numel(densities)
            dens = densities{j3};
            for j4=1:numel(Nvs)
                Nv = Nvs(j4);
                fname = sprintf('data/jb_1_comp_%s_%d_%s_%s_%d.mat', version, ord, geom, dens, Nv);
                fsmall = sprintf('data/jb_1_comp_%s_%d_%s_%s_%d_small.mat', ...
                                 version, ord, geom, dens, Nv);
                if exist(fname, 'file') && ~recompute
                    fprintf('(%s exists, not recomputing)\n', fname);
                    continue;
                end
                % Computation core
                core_comp(fname, fsmall, ord, geom, dens, Nv, plotting, version);
            end % for j1
        end % for j2
    end % for j3
end % for j4

function core_comp(fname, fsmall, ord, geom, dens, Nv, plotting, version)

    fprintf('Comp %s ...', fname);

    % Define torus-like surface
    if strcmp(geom, 'plain')
        a = 0.5;
        s = plaintorus(a);
    elseif strcmp(geom, 'wobbly')
        rng(2021);
        m = 3; % petal number of generating curve
        n = 2; % twist number along toroidal direction
        a = 0.25;
        s = wobblytorus(m,n,a);
    end
    k = 2; % aspect ratio

    % Define density on surface
    if strcmp(dens, 'const')
        f = @(u,v) ones(size(u)).*ones(size(v));
    elseif strcmp(dens, 'trig')
        f = @(u,v) cos(2*u).*sin(v+pi/4);
    end

    if plotting
        % Plot density on surface
        s = quadr_doubleptr(s, [24*k,24]);
        sig = f(s.u, s.v');
        sig = sig(:);
        sfigure(1); clf; jb_plot_torus(s, sig);
        drawnow
    end

    % Set up quadrature and precompute local correction weights
    s = quadr_doubleptr(s, [k*Nv,Nv]);
    As = Lap3dLocCorrectionOnly(s, ord);

    % Evaluate density
    sig = f(s.u, s.v');
    sig = sig(:);

    % Construct coarse reference grid
    Nv_ref = 8;
    sref = quadr_doubleptr(s, [k*Nv_ref,Nv_ref]);
    pts = s.x';
    pts_ref = sref.x';
    [I, D] = knnsearch(pts, pts_ref);
    assert(max(abs(D)) < 10*eps, 'Grid mismatch!');

    % Evaluate Laplace SL potential
    % First use punctured trapezoidal rule
    pot = laplace_sl_mex(pts, pts, sig .* s.w');
    % Then add local correction
    pot = pot + As*sig;

    if plotting
        % Plot potential
        sfigure(2); clf; jb_plot_torus(s, pot); shading flat
        drawnow
    end

    % Pick out potential at coarse grid
    pot_ref = pot(I,:);

    % Save results to disk
    fprintf(' saving ...');
    save(fname, 's', 'pot');
    s = sref; pot = pot_ref;
    save(fsmall, 's', 'pot');
    fprintf(' done.\n');

end

% vim:set shiftwidth=4 softtabstop=4:
