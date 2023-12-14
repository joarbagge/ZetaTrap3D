#include <mex.h>
#include <math.h>
#define INV_SQRT_PI 0.564189583547756286948079 // = 1/sqrt(pi)
#define INV_4_PI 0.0795774715459476678844419 // = 1/(4*pi)

#define TRG_IN prhs[0] // targets
#define SRC_IN prhs[1] // sources
#define DEN_IN prhs[2] // density
#define DEL_IN prhs[3] // delta
#define POT_OUT plhs[0] // potential

/* For documentation, see matlab/src/integrate_laplace_sl.m. */

double smoothfun1(const double x);
void core_punctured_trapz(double* restrict pot, const double* restrict targets,
                          const double* restrict sources, const double* restrict density,
                          const double* restrict delta, size_t Nt, size_t Ns,
                          size_t numel_delta);
void core_beale_regularization(double* restrict pot, const double* restrict targets,
                               const double* restrict sources, const double* restrict density,
                               const double* restrict delta, size_t Nt, size_t Ns,
                               size_t numel_delta);

// Entry point
// Call signature: pot = laplace_sl_mex(targets, sources, density, delta)
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // Input
    if (nrhs < 3) {
        mexErrMsgIdAndTxt("CapsMex:NotEnoughInputs",
                          "Not enough input arguments (at least three required).");
    }
    const double* restrict targets = mxGetPr(TRG_IN);
    const double* restrict sources = mxGetPr(SRC_IN);
    const double* restrict density = mxGetPr(DEN_IN);
    const double* restrict delta;
    size_t numel_delta = 0;
    bool punctured_trapz = false;
    if (nrhs < 4) {
        punctured_trapz = true;
    } else {
        delta = mxGetPr(DEL_IN);
        numel_delta = mxGetNumberOfElements(DEL_IN);
        if (numel_delta == 1 && delta[0] == 0) {
            punctured_trapz = true;
        }
    }
    size_t Nt = mxGetM(TRG_IN);
    size_t Ns = mxGetM(SRC_IN);
    if (mxGetN(TRG_IN) != 3 || mxGetN(SRC_IN) != 3) {
        mexErrMsgIdAndTxt("CapsMex:IncorrectDimensions",
                          "Incorrect number of columns in input (should be 3).");
    }
    if (mxGetN(DEN_IN) != 1) {
        mexErrMsgIdAndTxt("CapsMex:IncorrectDimensions",
                          "Incorrect number of columns in density (should be 1).");
    }
    if (mxGetM(DEN_IN) != Ns) {
        mexErrMsgIdAndTxt("CapsMex:IncorrectDimensions",
                          "Incorrect number of rows in density (should match sources).");
    }
    if (numel_delta > 1 && numel_delta != Ns) {
        mexErrMsgIdAndTxt("CapsMex:IncorrectDimensions",
                          "Incorrect number of delta values given (should match sources).");
    }

    // Allocate output
    POT_OUT = mxCreateDoubleMatrix(Nt, 1, mxREAL);
    double* restrict pot = mxGetPr(POT_OUT);

    if (punctured_trapz) {
        core_punctured_trapz(pot, targets, sources, density, delta,
                             Nt, Ns, numel_delta);
    } else {
        core_beale_regularization(pot, targets, sources, density, delta,
                                  Nt, Ns, numel_delta);
    }
}

double smoothfun1(const double x)
{
    double x2 = x*x;
    double b = 2.0*x2 - 5.0;
    double c = 2.0/3.0;
    return erf(x) - c*INV_SQRT_PI*x*b*exp(-x2);
}

void core_punctured_trapz(double* restrict pot, const double* restrict targets,
                          const double* restrict sources, const double* restrict density,
                          const double* restrict delta, size_t Nt, size_t Ns,
                          size_t numel_delta)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t j=0; j<Nt; j++) { // loop over targets
        double r[3];
        double norm_r, tmp;
        double target[] = {targets[j], targets[j+Nt], targets[j+Nt*2]};
        double sum = 0;
        for (size_t i=0; i<Ns; i++) { // loop over sources
            double source[] = {sources[i], sources[i+Ns], sources[i+Ns*2]};
            double dens = density[i];
            for (size_t d=0; d<3; d++) {
                r[d] = target[d] - source[d];
            }
            norm_r = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

            // Compute potential (q/||r||)
            tmp = 1 / norm_r;
            // Remove the point where r==0
            if (norm_r < 1e-14) {
                tmp = 0;
            }
            sum += dens * tmp;
        }
        pot[j] = sum * INV_4_PI;
    }
}

void core_beale_regularization(double* restrict pot, const double* restrict targets,
                               const double* restrict sources, const double* restrict density,
                               const double* restrict delta, size_t Nt, size_t Ns,
                               size_t numel_delta)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t j=0; j<Nt; j++) { // loop over targets
        double r[3];
        double delt, norm_r, norm_r_delta, tmp;
        double target[] = {targets[j], targets[j+Nt], targets[j+Nt*2]};
        double sum = 0;
        for (size_t i=0; i<Ns; i++) { // loop over sources
            double source[] = {sources[i], sources[i+Ns], sources[i+Ns*2]};
            double dens = density[i];
            if (numel_delta == 1) {
                delt = delta[0];
            } else {
                delt = delta[i];
            }
            for (size_t d=0; d<3; d++) {
                r[d] = target[d] - source[d];
            }
            norm_r = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
            norm_r_delta = norm_r / delt;

            // Compute potential (q/||r||)
            tmp = smoothfun1(norm_r_delta) / norm_r;
            // Known limit for r -> 0
            if (norm_r_delta < 2.11e-8) {
                tmp = 16*INV_SQRT_PI/(3*delt);
            }
            sum += dens * tmp;
        }
        pot[j] = sum * INV_4_PI;
    }
}

// vim:set shiftwidth=4 softtabstop=4:
