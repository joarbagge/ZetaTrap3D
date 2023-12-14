function jb_plot_torus(s, val)
% Based on showsurf.m
    Nr = [s.Nu,s.Nv]; % reshaping sizes
    x = reshape_wrapAround(s.x(1,:),Nr);
    y = reshape_wrapAround(s.x(2,:),Nr);
    z = reshape_wrapAround(s.x(3,:),Nr);
    if nargin >= 2
        val = reshape_wrapAround(val,Nr);
        surf(x,y,z,val);
    else
        surf(x,y,z);
    end
    axis equal; colorbar
end

function x = reshape_wrapAround(x0,Nr)
% reshape and wrap around to complete the torus-like surface for plotting
    x = reshape(x0,Nr);
    x = [x; x(1,:)]; x = [x, x(:,1)];
end

% vim:set shiftwidth=4 softtabstop=4:
