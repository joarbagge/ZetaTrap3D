% Run this in MATLAB to compile the MEX files

use_openmp = true;

if use_openmp
    % -fopenmp flag valid for GCC, other compilers may use a different flag
    mex CFLAGS='$CFLAGS -fopenmp' -lgomp laplace_sl_mex.c;
    mex CFLAGS='$CFLAGS -fopenmp' -lgomp stokes_sl_mex.c;
    %mex CFLAGS='$CFLAGS -fopenmp' -L/opt/MATLAB/R2022a/sys/os/glnxa64 -liomp5 laplace_sl_mex.c;
    %mex CFLAGS='$CFLAGS -fopenmp' -L/opt/MATLAB/R2022a/sys/os/glnxa64 -liomp5 stokes_sl_mex.c;
else
    mex laplace_sl_mex.c;
    mex stokes_sl_mex.c;
end
