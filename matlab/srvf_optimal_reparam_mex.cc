/*
 * libsrvf
 * =======
 *
 * A shape analysis library using the square root velocity framework.
 *
 * Copyright (C) 2012   FSU Statistical Shape Analysis and Modeling Group
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */
#include <srvf/plf.h>
#include <srvf/srvf.h>
#include <srvf/qmap.h>
#include <srvf/reparam.h>
#include <srvf/pointset.h>
#include <srvf/util.h>

#include <mex.h>
#include <vector>


void do_usage()
{
  mexPrintf(
    "USAGE: [G,T,s] = %s(X1,T1,X2,T2,tv1,tv2,nseeds)\n"
    "Inputs:\n"
    "\tX1, T1 = sample points and parameters of first curve\n"
    "\tX2, T2 = sample points and parameters of second curve\n"
    "\ttv1 = DP matching grid column parameters (default: use T1)\n"
    "\ttv2 = DP matching grid row parameters (default: use T2)\n"
    "\tnseeds = number of seed points to try on second curve (default is 1)\n"
    "Outputs:\n"
    "\tG,T = piecewise-linear reparametrization for second curve\n"
    "\ts : new seed should be tv2(s)\n",
    mexFunctionName()
  );
}

extern "C"
{
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Check number of arguments
  if (nlhs < 2 || nlhs > 3 || nrhs < 4 || nrhs > 7)
  { 
    mexPrintf("error: incorrect number of arguments\n");
    do_usage(); 
    return; 
  }

  // Check for incompatible SRVF dimensions
  size_t dim = mxGetM(prhs[0]);
  if (mxGetM(prhs[2]) != dim)
  { 
    mexPrintf("error: X1 and X2 have different dimensions\n"); 
    do_usage();
    return;
  }

  size_t npts1 = mxGetNumberOfElements(prhs[1]);
  size_t npts2 = mxGetNumberOfElements(prhs[3]);

  // Check for correct length on T1 and T2
  if ((npts1 != mxGetN(prhs[0])) || (npts2 != mxGetN(prhs[2])))
  {
    mexPrintf("error: incorrect length for T1 and/or T2\n");
    do_usage();
    return;
  }

  // Create sample pointsets.
  // Caller gives us the sample points in column-major, point-per-column 
  // ordering, which is the same thing as row-major, point-per-row ordering.
  srvf::Pointset samps1(dim, npts1, mxGetPr(prhs[0]), 
    srvf::Pointset::POINT_PER_ROW);
  srvf::Pointset samps2(dim, npts2, mxGetPr(prhs[2]), 
    srvf::Pointset::POINT_PER_ROW);

  // Create the parameter vectors T1 and T2
  double *T1_data = mxGetPr(prhs[1]);
  double *T2_data = mxGetPr(prhs[3]);
  std::vector<double> T1(&T1_data[0], &T1_data[npts1]);
  std::vector<double> T2(&T2_data[0], &T2_data[npts2]);

  // Create the grid parameter vectors tv1 and tv2
  std::vector<double> tv1, tv2;
  if (nrhs > 4){
    double *tv1_data = mxGetPr(prhs[4]);
    size_t ntv1 = mxGetNumberOfElements(prhs[4]);
    tv1 = std::vector<double>(&(tv1_data[0]), &(tv1_data[ntv1]));
  } else tv1 = T1;
  if (nrhs > 5){
    double *tv2_data = mxGetPr(prhs[5]);
    size_t ntv2 = mxGetNumberOfElements(prhs[5]);
    tv2 = std::vector<double>(&(tv2_data[0]), &(tv2_data[ntv2]));
  } else tv2 = T2;

  // Get the number of seedpoints to test
  size_t nseeds = 1;
  if (nrhs > 6)
  {
    nseeds = (size_t)mxGetScalar(prhs[6]);
  }

  // Find the optimal reparametrization using libsrvf
  srvf::Plf F1(samps1, T1);
  srvf::Plf F2(samps2, T2);
  srvf::Srvf Q1 = srvf::plf_to_srvf(F1);
  srvf::Srvf Q2 = srvf::plf_to_srvf(F2);
  srvf::Plf G = srvf::opencurves::optimal_reparam(Q1, Q2, tv1, tv2);

  // Allocate output variables
  plhs[0] = mxCreateDoubleMatrix(1, G.ncp(), mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1, G.ncp(), mxREAL);
  if (nlhs > 2) plhs[2] = mxCreateDoubleScalar(0.0);

  if (!plhs[0] || !plhs[1])
  {
    if (plhs[0]) mxDestroyArray(plhs[0]);
    if (plhs[1]) mxDestroyArray(plhs[1]);

    mexErrMsgIdAndTxt(mexFunctionName(), "mxCreateDoubleMatrix() failed.");
  }

  double *G_data = mxGetPr(plhs[0]);
  double *T_data = mxGetPr(plhs[1]);

  for (size_t i=0; i<G.ncp(); ++i)
  {
    G_data[i] = G.samps()[i][0];
    T_data[i] = G.params()[i];
  }
}
} // extern "C"
