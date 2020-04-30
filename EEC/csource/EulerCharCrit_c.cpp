//  array_calculation.cpp
//  
//
//  Created by wendy lin on 6/14/19.
//  Edited and improved by F. Telschow 7/11/19
//

#include "ECcrits_c.cpp"
#include <math.h>
#include <matrix.h>
#include <mex.h>

/* This is the gateway function connecting C++ and Matlab */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /******************************************************************
   * define required variables, pointers and constants
   ******************************************************************/
  double *z, *out;   // initialize pointers for input (z) and output (out)
  int cc;            // initialize conncetivity as an integer
  //
  mwSize i, j, k, n, n_nobound;
  // get pointer array with dimensions of input matrix
  const mwSize* DimMatrix = mxGetDimensions(prhs[0]);
    
  /******************************************************************
   * Check number of input and output arguments
   ******************************************************************/
  if (nrhs != 2) { 
    mexErrMsgIdAndTxt( "MATLAB:timestwoalt:invalidNumInputs",
            "Two input arguments required."); 
  } else if (nlhs > 1) {
    mexErrMsgIdAndTxt( "MATLAB:timestwoalt:maxlhs",
            "Too many output arguments."); 
  }
    
  /******************************************************************
   * save pointers/variables for input
   ******************************************************************/
  z  = (double *)mxGetPr(prhs[0]); // get pointer to input (z)
  cc = mxGetScalar(prhs[1]);       // get scalar input cc
  
  if( mxGetNumberOfDimensions(prhs[0]) == 3 ){
        // get dimensions of the matrix
        i = DimMatrix[0];
        j = DimMatrix[1];
        k = DimMatrix[2];

        n_nobound = (i-2) * (j-2) * (k-2); // number of elements in input array without padding

        // create an output matrix of size 2 x nvoxel of not padded image
        plhs[0] = mxCreateDoubleMatrix(2, n_nobound, mxREAL);

        /******************************************************************
        * create a C pointer to a copy of the output matrix 
        ******************************************************************/
        out = mxGetPr(plhs[0]);

        /******************************************************************
        * call the C subroutine
        ******************************************************************/
         ECcrit3D(z, cc, out, i, j, k);
        
  } else if( mxGetNumberOfDimensions(prhs[0]) == 2 && DimMatrix[1] != 1 ) {
          // get dimensions of the matrix
          i = DimMatrix[0];
          j = DimMatrix[1];

          n_nobound = (i-2) * (j-2); // number of elements in input array without padding

          // create an output matrix of size 2 x nvoxel of not padded image
          plhs[0] = mxCreateDoubleMatrix( 2, n_nobound, mxREAL );

          /******************************************************************
           * create a C pointer to a copy of the output matrix 
           ******************************************************************/
           out = mxGetPr( plhs[0] );

          /******************************************************************
           * call the C subroutine
           ******************************************************************/  
           ECcrit2D(z, cc, out, i, j);
  } else if( mxGetNumberOfDimensions(prhs[0]) == 2 && DimMatrix[1] == 1 ) {
          // get dimensions of the matrix
          i = DimMatrix[0];

          n_nobound = (i-2); // number of elements in input array without padding

          // create an output matrix of size 2 x nvoxel of not padded image
          plhs[0] = mxCreateDoubleMatrix(2, n_nobound, mxREAL);

          /******************************************************************
           * create a C pointer to a copy of the output matrix 
           ******************************************************************/
           out = mxGetPr(plhs[0]);

          /******************************************************************
           * call the C subroutine
           ******************************************************************/  
            ECcrit1D(z, out, i);
  } else {
            mexErrMsgIdAndTxt( "MATLAB:missing implementation",
                                "dimension must be smaller than 3."); 
  }
}
