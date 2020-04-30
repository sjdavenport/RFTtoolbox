/***********************************************
 *  This file contains c implementations of the
 *  Euler characteristic for different dimensions
 *  using the lower star/critical value trick
 * 
***********************************************/

#include <math.h>
#include <matrix.h>

/* This is the C++ subroutine computing the change in Euler characterisitic in 2D
 * Input:
 *   - *z: pointer to the input array 
 *   - cc: integer giving the connectivity (currently, only 8 is supported)
 *   - ADD all inputs with short discribtion
 */
void ECcrit1D(double *z, double *out, mwSize i)
{
    /*****************************************************************
     * Initialize variables and constants
     *****************************************************************/
    int ii, s_ind, ss, ind_ss; // variables for loops
    int dec = 0;   // initialize the change in Euler characteristic as being 1, since 1 voxel is added
    int nn  = 0;   // counter for accessing the values in dec while looping
    double x[9];  // initialize array for values of neighbourhood of a voxel
    int ind_t[9]; // initialize array for indices of neighbourhood of a voxel
    
    /*****************************************************************
     * loop over each voxel of the input z. Note that the original image
     * was zero paded by with one -Inf. Hence we start each dimension at
     * 1 and end at end-1.
     *****************************************************************/
    for(ii=1;ii<(i-1);ii++)
    {
            /*
             * getting the indices of the 3x3x3 neighbourhood of
             * the voxel (ii,jj,kk) and saving them in the ind_t array
             */
            ind_t[0] = ii-1;
            ind_t[1] = ii;
            ind_t[2] = ii+1;

            /* filling the x array using the indices of the
             * neighbourhood ind_t and the input array z
             */
            for(ss = 0; ss < 2; ss++)
            {
                s_ind = ind_t[ss];
                x[ss] = *(z+s_ind);
            }

            /*****************************************************************
             * let the even entries of out pointer point to the center voxel (ii,jj)
             *****************************************************************/
            *(out + nn*2) = x[1];

            /******************************************************************
             * compute the change in dec by checking which nD-faces are created
             ******************************************************************/
            // minima increase EC by one
            if(x[1] > x[0] && x[1] < x[2])
                dec = 1;
            // maxima decrease EC by one
            if(x[0] < x[1] && x[1] > x[2])
                dec = -1;
            
            // save the dec in the odd out locations
            *(out + nn*2 + 1) = dec;
            // increase counter moving through dec saving locations
            nn += 1;
            // reset dec to be 1
            dec = 0;
        }
}


/* This is the C++ subroutine computing the change in Euler characterisitic in 2D
 * Input:
 *   - *z: pointer to the input array 
 *   - cc: integer giving the connectivity (currently, only 8 is supported)
 *   - ADD all inputs with short discribtion
 */
void ECcrit2D(double *z, double cc, double *out, mwSize i, mwSize j)
{
    /*****************************************************************
     * Initialize variables and constants
     *****************************************************************/
    int ii, jj, s_ind, ss, ind_ss; // variables for loops
    int dec = 1;   // initialize the change in Euler characteristic as being 1, since 1 voxel is added
    int nn  = 0;   // counter for accessing the values in dec while looping
    double x[9];   // initialize array for values of neighbourhood of a voxel
    int ind_t[9];  // initialize array for indices of neighbourhood of a voxel
    
    /*****************************************************************
     * loop over each voxel of the input z. Note that the original image
     * was zero paded by with one -Inf. Hence we start each dimension at
     * 1 and end at end-1.
     *****************************************************************/
    for(ii=1;ii<(i-1);ii++)
    {
        for(jj=1;jj<(j-1);jj++)
        {
            /*
             * getting the indices of the 3x3 neighbourhood of
             * the voxel (ii,jj) and saving them in the ind_t array
             */
            ind_t[0] = (ii-1) + (jj-1)*i;
            ind_t[1] = (ii-1) + jj*i;
            ind_t[2] = (ii-1) + (jj+1)*i;
            ind_t[3] = ii     + (jj-1)*i;
            ind_t[4] = ii     + jj*i;
            ind_t[5] = ii     + (jj+1)*i;
            ind_t[6] = (ii+1) + (jj-1)*i;
            ind_t[7] = (ii+1) + jj*i;
            ind_t[8] = (ii+1) + (jj+1)*i;
            
            /* filling the x array using the indices of the
             * neighbourhood ind_t and the input array z
             */
            for(ss = 0; ss < 9; ss++)
            {
                s_ind = ind_t[ss];
                x[ss] = *(z+s_ind);
            }

            /*****************************************************************
             * let the even entries of out pointer point to the center voxel (ii,jj)
             *****************************************************************/
            *(out + nn*2) = x[4];

            /******************************************************************
             * compute the change in dec by checking which nD-faces are created
             ******************************************************************/
            // subtract number of edges
            double ind1[4] = {1, 3, 5, 7};
            double count_t = 0;
            for(ss=0; ss < cc; ss++){
                ind_ss = ind1[ss];
                // if new edge appears increase dec by 1
                if(x[ind_ss] > x[4])
                {
                    count_t += 1;
                }
            }

            dec = dec - count_t;

            // add number of faces
            if(x[0] > x[4] && x[1] > x[4] && x[3] > x[4])
              dec += 1;
            if(x[1] > x[4] && x[2] > x[4] && x[5] > x[4])
              dec += 1;
            if(x[3] > x[4] && x[6] > x[4] && x[7] > x[4])
              dec += 1;
            if(x[5] > x[4] && x[7] > x[4] && x[8] > x[4])
              dec += 1;

            // save the negative dec in the odd out locations (minus
            // sign convention used to recreate the ec curve)
            *(out + nn*2 + 1) = -dec;
            // increase counter moving through dec saving locations
            nn += 1;
            // reset dec to be 1
            dec = 1;
        }
    }
}


/* This is the C++ subroutine computing the change in Euler characterisitic in 3D
 * Input:
 *   - *z: pointer to the input array 
 *   - cc: integer giving the connectivity (currently, only 8 is supported)
 *   - ADD all inputs with short discribtion
 */
void ECcrit3D(double *z, double cc, double *out, mwSize i, mwSize j, mwSize k)
{
    /*****************************************************************
     * Initialize variables and constants
     *****************************************************************/
    int ii, jj, kk, s_ind, ss, ind_ss; // variables for loops
    int dec = 1;   // initialize the change in Euler characteristic as being 1, since 1 voxel is added
    int nn  = 0;   // counter for accessing the values in dec while looping
    double x[27];  // initialize array for values of neighbourhood of a voxel
    int ind_t[27]; // initialize array for indices of neighbourhood of a voxel
    
    /*****************************************************************
     * loop over each voxel of the input z. Note that the original image
     * was zero paded by with one -Inf. Hence we start each dimension at
     * 1 and end at end-1.
     *****************************************************************/
    for(ii=1;ii<(i-1);ii++)
    {
        for(jj=1;jj<(j-1);jj++)
        {
            for(kk=1;kk<(k-1);kk++)
            {
                /*
                 * getting the indices of the 3x3x3 neighbourhood of
                 * the voxel (ii,jj,kk) and saving them in the ind_t array
                 */
                ind_t[0] = (ii-1)+(jj-1)*i+(kk-1)*i*j;
                ind_t[1] = (ii-1)+(jj-1)*i+kk*i*j;
                ind_t[2] = (ii-1)+(jj-1)*i+(kk+1)*i*j;
                ind_t[3] = (ii-1)+jj*i+(kk-1)*i*j;
                ind_t[4] = (ii-1)+jj*i+kk*i*j;
                ind_t[5] = (ii-1)+jj*i+(kk+1)*i*j;
                ind_t[6] = (ii-1)+(jj+1)*i+(kk-1)*i*j;
                ind_t[7] = (ii-1)+(jj+1)*i+kk*i*j;
                ind_t[8] = (ii-1)+(jj+1)*i+(kk+1)*i*j;
                ind_t[9] = ii+(jj-1)*i+(kk-1)*i*j;
                ind_t[10] = ii+(jj-1)*i+kk*i*j;
                ind_t[11] = ii+(jj-1)*i+(kk+1)*i*j;
                ind_t[12] = ii+jj*i+(kk-1)*i*j;
                ind_t[13] = ii+jj*i+kk*i*j;
                ind_t[14] = ii+jj*i+(kk+1)*i*j;
                ind_t[15] = ii+(jj+1)*i+(kk-1)*i*j;
                ind_t[16] = ii+(jj+1)*i+kk*i*j;
                ind_t[17] = ii+(jj+1)*i+(kk+1)*i*j;
                ind_t[18] = (ii+1)+(jj-1)*i+(kk-1)*i*j;
                ind_t[19] = (ii+1)+(jj-1)*i+kk*i*j;
                ind_t[20] = (ii+1)+(jj-1)*i+(kk+1)*i*j;
                ind_t[21] = (ii+1)+jj*i+(kk-1)*i*j;
                ind_t[22] = (ii+1)+jj*i+kk*i*j;
                ind_t[23] = (ii+1)+jj*i+(kk+1)*i*j;
                ind_t[24] = (ii+1)+(jj+1)*i+(kk-1)*i*j;
                ind_t[25] = (ii+1)+(jj+1)*i+kk*i*j;
                ind_t[26] = (ii+1)+(jj+1)*i+(kk+1)*i*j;
                
                /* filling the x array using the indices of the
                 * neighbourhood ind_t and the input array z
                 */
                for(ss = 0; ss<27; ss++)
                {
                    s_ind = ind_t[ss];
                    x[ss] = *(z+s_ind);
                }
                
                /*****************************************************************
                 * let the even entries of out pointer point to the center voxel (ii,jj,kk)
                 *****************************************************************/
                *(out + nn*2) = x[13];
                
                /******************************************************************
                 * compute the change in dec by checking which nD-faces are created
                 ******************************************************************/
                // subtract number of edges (1D faces)
                double ind1[6] = {4, 10, 12, 14, 16, 22};
                double count_t = 0;
                // find number of newly created 1D faces
                for(ss=0;ss<6;ss++){
                    ind_ss = ind1[ss];
                        // if new edge appears increase dec by 1
                        if(x[ind_ss]>x[13])
                        {
                            count_t += 1;
                        }
                }
                // subtract number of 1D faces from dec
                dec = dec - count_t;

                // add number of newly created faces
                if(x[9] > x[13] && x[10] > x[13] && x[12] > x[13])
                  dec += 1;
                if(x[10] > x[13] && x[11] > x[13] && x[14] > x[13])
                  dec += 1;
                if(x[14] > x[13] && x[16] > x[13] && x[17] > x[13])
                  dec += 1;
                if(x[12] > x[13] && x[15] > x[13] && x[16] > x[13])
                  dec += 1;
                if(x[3] > x[13] && x[4] > x[13] && x[12] > x[13])
                  dec += 1;
                if(x[4] > x[13] && x[5] > x[13] && x[14] > x[13])
                  dec += 1;
                if(x[12] > x[13] && x[21] > x[13] && x[22] > x[13])
                  dec += 1;
                if(x[14] > x[13] && x[22] > x[13] && x[23] > x[13])
                  dec += 1;
                if(x[1] > x[13] && x[4] > x[13] && x[10] > x[13])
                  dec += 1;
                if(x[4] > x[13] && x[7] > x[13] && x[16] > x[13])
                  dec += 1;
                if(x[10] > x[13] && x[22] > x[13] && x[19] > x[13])
                  dec += 1;
                if(x[16] > x[13] && x[22] > x[13] && x[25] > x[13])
                  dec += 1;

                // subtract number of newly created 3D faces
                if(x[0] >= x[13] && x[1] >= x[13] && x[3] >= x[13] && x[4] >= x[13] && x[9] >= x[13] && x[10] >= x[13] && x[12] >= x[13] && x[13] >= x[13])
                  dec -= 1;
                if(x[1] >= x[13] && x[2] >= x[13] && x[4] >= x[13] && x[5] >= x[13] && x[10] >= x[13] && x[11] >= x[13] && x[13] >= x[13] && x[14] >= x[13])
                  dec -= 1;
                if(x[4] >= x[13] && x[5] >= x[13] && x[7] >= x[13] && x[8] >= x[13] && x[13] >= x[13] && x[14] >= x[13] && x[16] >= x[13] && x[17] >= x[13])
                  dec -= 1;
                if(x[3] >= x[13] && x[4] >= x[13] && x[6] >= x[13] && x[7] >= x[13] && x[12] >= x[13] && x[13] >= x[13] && x[15] >= x[13] && x[16] >= x[13])
                  dec -= 1;
                if(x[9] >= x[13] && x[10] >= x[13] && x[12] >= x[13] && x[13] >= x[13] && x[18] >= x[13] && x[19] >= x[13] && x[21] >= x[13] && x[22] >= x[13])
                  dec -= 1;
                if(x[10] >= x[13] && x[11] >= x[13] && x[13] >= x[13] && x[14] >= x[13] && x[19] >= x[13] && x[20] >= x[13] && x[22] >= x[13] && x[23] >= x[13])
                  dec -= 1;
                if(x[13] >= x[13] && x[14] >= x[13] && x[16] >= x[13] && x[17] >= x[13] && x[22] >= x[13] && x[23] >= x[13] && x[25] >= x[13] && x[26] >= x[13])
                  dec -= 1;
                if(x[12] >= x[13] && x[13] >= x[13] && x[15] >= x[13] && x[16] >= x[13] && x[21] >= x[13] && x[22] >= x[13] && x[24] >= x[13] && x[25] >= x[13])
                  dec -= 1;
                                
                // save the negative dec in the odd out locations (minus
                // sign convention used to recreate the ec curve)
                *(out + nn*2 + 1) = -dec;
                // increase counter moving through dec saving locations
                nn += 1;
                // reset dec to be 1
                dec = 1;
            }
        }
    }
}