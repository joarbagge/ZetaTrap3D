% Compute error and plot, for convergence test
% Uses potential values computed by jb_1_comp and stored in data/

orders = [3, 5];
geometries = {'plain', 'wobbly'};
densities = {'const', 'trig'};
Nvs = 8*[1 2 4 8];
plotting = true;
plot_avg = true;
version = 'v1';
err_against_ref = false;
ref = {version, 5, 64};

for j2=1:numel(geometries)
    geom = geometries{j2};
    for j3=1:numel(densities)
        dens = densities{j3};
        if plotting
            sfigure(j2*10 + j3); clf;
            maxV = -Inf; minV = Inf;
            labels = {};
            xV = [];
        end
        if err_against_ref
            % Start by loading reference solution
            fname = sprintf('data/jb_1_comp_%s_%d_%s_%s_%d_small.mat', ...
                            ref{1}, ref{2}, geom, dens, ref{3});
            [sref, pot_ref] = load_data(fname);
        end
        for j1=1:numel(orders)
            % Estimate pointwise error
            ord = orders(j1);
            potdiff = {};
            ss = {};
            I = [];
            for j4=1:numel(Nvs)
                Nv = Nvs(j4);
                fname = sprintf('data/jb_1_comp_%s_%d_%s_%s_%d_small.mat', ...
                                version, ord, geom, dens, Nv);
                if ~exist(fname, 'file')
                    continue;
                end
                [s, pot] = load_data(fname);
                if err_against_ref
                    potdiff{end+1} = abs(pot - pot_ref);
                    ss{end+1} = sref;
                    I(end+1) = j4;
                elseif j4 > 1
                    potdiff{end+1} = abs(pot - pot_old);
                    ss{end+1} = sold;
                    I(end+1) = j4-1;
                end
                sold = s;
                pot_old = pot;
            end % for j4
            % Compute average and maximum error
            avgErr = [];
            maxErr = [];
            for ii=1:numel(I)
                avgErr(end+1) = ss{ii}.w * potdiff{ii} / sum(ss{ii}.w);
                maxErr(end+1) = max(potdiff{ii});
            end
            % Compute average and maximum value (for normalization)
            if err_against_ref
                avgVal = sref.w * abs(pot_ref) / sum(sref.w);
                maxVal = max(abs(pot_ref));
            else
                avgVal = sold.w * abs(pot_old) / sum(sold.w);
                maxVal = max(abs(pot_old));
            end
            % Plot results
            if plotting
                [lab, k, maxv, minv, xv] = core_plot(ss, Nvs(I), avgErr, maxErr, ...
                                                     avgVal, maxVal, ord, plot_avg);
                maxV = max(maxV, maxv);
                minV = min(minV, minv);
                for jj=1:numel(lab)
                    labels{end+1} = lab{jj};
                end
                xV = unique([xV, xv]);
            end
        end % for j1
        if plotting
            title(sprintf('Laplace, Zeta v1A, %s torus, %s density', geom, dens));
            xlabel('Tot. num. grid points')
            ylabel('Est. rel. error')
            ylim([1e-8, 1e1])
            xlim([1e2, 48e3])
            set(gca, 'XTick', xV);
            grid on
            slopeMarker(3, 1, 4, 2*maxV, Nvs(I(1)), k(1), 0, 0, '--', 'Color', 'k');
            %loglog(Ntot, 1.5*maxErrRel(end)*(N/N(end)).^-ord, 'k:')
            %loglog(Ntot, (1/1.5)*avgErrRel(end)*(N/N(end)).^-ord, 'k:')
            legend(labels, 'Location', 'best');
        end
    end % for j2
end % for j3

function [s, pot] = load_data(fname)
    fd = load(fname);
    s = fd.s;
    pot = fd.pot;
end

function [labels, k, maxv, minv, xv] = core_plot(s, Nv, avgErr, maxErr, avgVal, maxVal, order, ...
                                                 plot_avg)
    k = [];
    for i=1:numel(s)
        k(end+1) = s{i}.Nu/s{i}.Nv;
    end
    Ntot = Nv.*(k.*Nv);
    %avgErrRel = avgErr / avgVal;
    avgErrRel = avgErr / maxVal;
    maxErrRel = maxErr / maxVal;
    % Start plotting
    labels = {};
    loglog(Ntot, maxErrRel, '^-', 'LineWidth', 1)
    labels{end+1} = sprintf('Z1A.%d (max)', order);
    hold on
    maxv = max(maxErrRel);
    if plot_avg
        set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')-1);
        loglog(Ntot, avgErrRel, '--', 'LineWidth', 1);
        labels{end+1} = sprintf('(avg)');
        minv = min(avgErrRel);
    else
        minv = min(maxErrRel);
    end
    xv = Ntot;
end

% vim:set shiftwidth=4 softtabstop=4:
