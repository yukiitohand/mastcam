/* =====================================================================
 * proj_mastcam2MSLDEM_v4_mex.c
 * Read MSLDEM image data
 * Perform projection of mastcam pixels onto MSLDEM data
 * 
 * INPUTS:
 * 0 msldem_imgpath        char*
 * 1 msldem_header         struct
 * 2 msldemc_header        struct
 * 3 msldemc_northing      Double array [L_demc]
 * 4 msldemc_easting       Double array [S_demc]
 * 5 msldemc_imx           Doubles [S_demc x L_demc]
 * 6 msldemc_imy           Doubles [S_demc x L_demc]
 * 7 msldemc_imFOVmask     Boolean [S_demc x L_demc]
 * 8 msldemc_imFOVmaskd    Boolean [S_demc x L_demc]
 * 9 S_im                  int
 * 10 L_im                  int
 * 11 cam_C                 Double [3 x 1]
 * 12 cam_A                Double [3 x 1]
 * 13 PmC                  Double [3 x (L_im*S_im)]
 * 
 * 
 * OUTPUTS:
 * 0  im_north       [L_im x S_im]   Double / Float
 * 1  im_east        [L_im x S_im]
 * 2  im_elev        [L_im x S_im]
 * 3  im_refx        [L_im x S_im]   Integer
 * 4  im_refy        [L_im x S_im]   Integer
 * 5  im_refs        [L_im x S_im]   Integer
 * 6  im_range       [L_im x S_im]       Double / Float
 * 7  im_nnx         [L_im x S_im]   Integer
 * 8  im_nny         [L_im x S_im]   Integer
 * 
 *
 *
 * This is a MEX file for MATLAB.
 * ===================================================================== */
#include <stdint.h>
#include "io64.h"
#include "mex.h"
#include "math.h"
#include "matrix.h"
#include <string.h>
#include <stdio.h>

#include <stdlib.h>
#include "envi.h"

double get_sqr_dst(double ux,double uy,double uz,double vx,double vy,double vz)
{
    double dst;
    dst = pow(ux-vx,2) + pow(uy-vy,2) + pow(uz-vz,2);
    return dst;
}

/* main computation routine */
void proj_mastcam2MSLDEM(char *msldem_imgpath, EnviHeader msldem_header, 
        int msldemc_imxy_sample_offset, int msldemc_imxy_line_offset,
        int msldemc_samples, int msldemc_lines,
        double *msldemc_northing, double *msldemc_easting,
        double **msldemc_imx, double **msldemc_imy,
        bool **msldemc_imFOVmask, bool **msldemc_imFOVmaskd, int S_im, int L_im, double *cam_C, double *cam_A,
        double **PmC, double **im_north, double **im_east, double **im_elev,
        int32_T **im_refx, int32_T **im_refy, int32_T **im_refs,
        double **im_range, int32_T **im_nnx, int32_T **im_nny)
{
    int32_T c,l;
    int16_T sxy;
    long skip_pri;
    long skip_l, skip_r;
    float *elevl,*elevlp1;
    long ncpy;
    const int sz=sizeof(float);
    FILE *fid;
    double ppv1x,ppv1y,ppv2x,ppv2y,ppv3x,ppv3y;
    double ppv1gx,ppv1gy,ppv1gz,ppv2gx,ppv2gy,ppv2gz,ppv3gx,ppv3gy,ppv3gz,ppv4gx,ppv4gy,ppv4gz;
    double pdv1x,pdv1y,pdv1z,pdv2x,pdv2y,pdv2z;
    double detM;
    double Minv[2][2];
    int32_T x_min,y_min,x_max,y_max;
    int32_T xi,yi;
    double pipvx,pipvy;
    double pipvgx,pipvgy,pipvgz;
    double rnge;
    double pprm_sd,pprm_td; /* plane parameter for projected image plane */
    double dv1,dv2,dv3;
    double pd1,pd2,pd3,pp_sum,pd_dem;
    double pprm_s,pprm_t,pprm_1st;
    int32_T ppvx[4];
    int32_T ppvy[4];
    double dst_ppv[4];
    double dst_nn;
    int32_T di,di_dst_min;
    bool isinFOVd,isinFOV;
    
    
    
    fid = fopen(msldem_imgpath,"rb");
    
    /* skip lines */
    skip_pri = (long) msldem_header.samples * (long) msldemc_imxy_line_offset * (long) sz;
    // printf("%d*%d*%d=%ld\n",msldem_header.samples,msldemc_imxy_line_offset,s,skip_pri);
    fseek(fid,skip_pri,SEEK_CUR);
    
    /* read the data */
    elevl = (float*) malloc(sz*msldemc_samples);
    elevlp1 = (float*) malloc(sz*msldemc_samples);
    skip_l = (long) sz * (long) msldemc_imxy_sample_offset;
    skip_r = ((long) msldem_header.samples - (long) msldemc_samples)* (long) sz - skip_l;
    ncpy = (long) msldemc_samples* (long)sz;
    fseek(fid,skip_l,SEEK_CUR);
    fread(elevlp1,sz,msldemc_samples,fid);
    fseek(fid,skip_r,SEEK_CUR);
    
    // printf("%d,%d,%d\n",skip_l,msldemc_samples*s,skip_r);
    for(l=0;l < (msldemc_lines-1);l++){
        memcpy(elevl,elevlp1,ncpy);
        fseek(fid,skip_l,SEEK_CUR);
        fread(elevlp1,sz,msldemc_samples,fid);
        fseek(fid,skip_r,SEEK_CUR);
        //printf("l=%d\n",l);
        // decide the first and last indexes to be assessed.
        for(c=0;c<(msldemc_samples-1);c++){
            // process if 
            //printf("c=%d,mask_lc = %d,mask_lp1c = %d\n",c,msldemc_imFOVmask[c][l],msldemc_imFOVmask[c][l+1]);
            if(msldemc_imFOVmask[c][l] || msldemc_imFOVmask[c][l+1]){
                //printf("c=%d\n",c);
                for(sxy=0;sxy<2;sxy++){
                    if(sxy==0){
                        ppv1x = msldemc_imx[c][l];
                        ppv1y = msldemc_imy[c][l];
                        ppv2x = msldemc_imx[c+1][l];
                        ppv2y = msldemc_imy[c+1][l];
                        ppv3x = msldemc_imx[c][l+1];
                        ppv3y = msldemc_imy[c][l+1];
                        ppv1gx = msldemc_northing[l];
                        ppv1gy = msldemc_easting[c];
                        ppv1gz = (double) -elevl[c];
                        ppv2gx = msldemc_northing[l];
                        ppv2gy = msldemc_easting[c+1];
                        ppv2gz = (double) -elevl[c+1];
                        ppv3gx = msldemc_northing[l+1];
                        ppv3gy = msldemc_easting[c];
                        ppv3gz = (double) -elevlp1[c];
                        ppv4gx = msldemc_northing[l+1];
                        ppv4gy = msldemc_easting[c+1];
                        ppv4gz = (double) -elevlp1[c+1];
                        ppvx[0] = c; ppvx[1] = c+1; ppvx[2] = c; ppvx[3] = c+1;
                        ppvy[0] = l; ppvy[1] = l; ppvy[2] = l+1; ppvy[3] = l+1;
                        isinFOVd = (msldemc_imFOVmaskd[c][l] && msldemc_imFOVmaskd[c+1][l] && msldemc_imFOVmaskd[c][l+1]);
                        isinFOV = (msldemc_imFOVmask[c][l] && msldemc_imFOVmask[c+1][l] && msldemc_imFOVmask[c][l+1]);
                        
                    }
                    else{
                        ppv1x = msldemc_imx[c+1][l];
                        ppv1y = msldemc_imy[c+1][l];
                        ppv2x = msldemc_imx[c+1][l+1];
                        ppv2y = msldemc_imy[c+1][l+1];
                        ppv3x = msldemc_imx[c][l+1];
                        ppv3y = msldemc_imy[c][l+1];
                        ppv1gx = msldemc_northing[l];
                        ppv1gy = msldemc_easting[c+1];
                        ppv1gz = (double) -elevl[c+1];
                        ppv2gx = msldemc_northing[l+1];
                        ppv2gy = msldemc_easting[c+1];
                        ppv2gz = (double) -elevlp1[c+1];
                        ppv3gx = msldemc_northing[l+1];
                        ppv3gy = msldemc_easting[c];
                        ppv3gz = (double) -elevlp1[c];
                        ppv4gx = msldemc_northing[l];
                        ppv4gy = msldemc_easting[c];
                        ppv4gz = (double) -elevl[c];
                        ppvx[0] = c+1; ppvx[1] = c+1; ppvx[2] = c; ppvx[3] = c;
                        ppvy[0] = l; ppvy[1] = l+1; ppvy[2] = l+1; ppvy[3] = l;
                        isinFOVd = (msldemc_imFOVmaskd[c+1][l] && msldemc_imFOVmaskd[c+1][l+1] && msldemc_imFOVmaskd[c][l+1]);
                        isinFOV = (msldemc_imFOVmask[c][l] && msldemc_imFOVmask[c+1][l] && msldemc_imFOVmask[c][l+1]);
                        
                    }
                    //printf("sxy=%d\n",sxy);
                    //printf("isinFOVd=%d\n",isinFOVd);
                    if(isinFOVd)
                    {
                        // define some plane parameters
                        pdv1x = ppv2x - ppv1x; pdv1y = ppv2y - ppv1y;
                        pdv2x = ppv3x - ppv1x; pdv2y = ppv3y - ppv1y;
                        detM = pdv1x*pdv2y - pdv1y*pdv2x;
                        Minv[0][0] = pdv2y/detM;
                        Minv[0][1] = -pdv2x/detM;
                        Minv[1][0] = -pdv1y/detM;
                        Minv[1][1] = pdv1x/detM;
                        
                        
                        
                        // parameters for convert pprm in the image 
                        // plane into those for polygonal surface.
                        dv1 = cam_A[0]*(ppv1gx-cam_C[0]) + cam_A[1]*(ppv1gy-cam_C[1]) + cam_A[2]*(ppv1gz-cam_C[2]);
                        dv2 = cam_A[0]*(ppv2gx-cam_C[0]) + cam_A[1]*(ppv2gy-cam_C[1]) + cam_A[2]*(ppv2gz-cam_C[2]);
                        dv3 = cam_A[0]*(ppv3gx-cam_C[0]) + cam_A[1]*(ppv3gy-cam_C[1]) + cam_A[2]*(ppv3gz-cam_C[2]);
                        
                        
                        
                        x_min = ceil(fmin(fmin(ppv1x,ppv2x),ppv3x));
                        y_min = ceil(fmin(fmin(ppv1y,ppv2y),ppv3y));
                        x_max = ceil(fmax(fmax(ppv1x,ppv2x),ppv3x));
                        y_max = ceil(fmax(fmax(ppv1y,ppv2y),ppv3y));
                        
                        if(x_min<0)
                            x_min = 0;
                        if(x_min>S_im)
                            x_min = S_im;
                        if(x_max<0)
                            x_max = 0;
                        if(x_max>S_im)
                            x_max = S_im;
                        if(y_min<0)
                            y_min = 0;
                        if(y_min>L_im)
                            y_min = L_im;
                        if(y_max<0)
                            y_max = 0;
                        if(y_max>L_im)
                            y_max = L_im;
                        
                        
                        for(xi=x_min;xi<x_max;xi++){
                            for(yi=y_min;yi<y_max;yi++){
                                pipvx = (double) xi - ppv1x; pipvy = (double) yi - ppv1y; 
                                pprm_sd = Minv[0][0]*pipvx+Minv[0][1]*pipvy;
                                pprm_td = Minv[1][0]*pipvx+Minv[1][1]*pipvy;
                                pp_sum = pprm_sd+pprm_td;
                                if(pprm_sd>=0 && pprm_td>=0 && pp_sum<=1)
                                {
                                    // Conversion from plane parameters in the image plane
                                    // into plane parameters on the polygonal surface.
                                    pd2 = pprm_sd / dv2;
                                    pd3 = pprm_td / dv3;
                                    pd1 = (1-pp_sum) / dv1;
                                    pd_dem = pd1+pd2+pd3;
                                    pprm_s = pd2 / pd_dem;
                                    pprm_t = pd3 / pd_dem;
                                    pprm_1st = 1 - pprm_s - pprm_t;
                                    //printf("isinFOVd=%d\n",isinFOVd);
                                    
                                    // Evaluate distance
                                    //pipvgx = (1-pp_sum)*ppv1gx+pprm_px*ppv2gx+pprm_py*ppv3gx;
                                    //pipvgy = (1-pp_sum)*ppv1gy+pprm_px*ppv2gy+pprm_py*ppv3gy;
                                    //pipvgz = (1-pp_sum)*ppv1gz+pprm_px*ppv2gz+pprm_py*ppv3gz;
                                    pipvgx = pprm_1st*ppv1gx+pprm_s*ppv2gx+pprm_t*ppv3gx;
                                    pipvgy = pprm_1st*ppv1gy+pprm_s*ppv2gy+pprm_t*ppv3gy;
                                    pipvgz = pprm_1st*ppv1gz+pprm_s*ppv2gz+pprm_t*ppv3gz;
                                    //printf("isinFOVd=%d\n",isinFOVd);

                                    rnge = pow(pipvgx-cam_C[0],2) + pow(pipvgy-cam_C[1],2) + pow(pipvgz-cam_C[2],2);
                                    // printf("isinFOVd=%d\n",isinFOVd);
                                    // printf("rnge=%f\n",rnge);
                                    // printf("im_range[%d][%d]=%f\n",xi,yi,im_range[xi][yi]);
                                    

                                    if(rnge < im_range[xi][yi])
                                    {
                                        //printf("isinFOVd=%d\n",isinFOVd);
                                        im_range[xi][yi] = rnge;
                                        im_north[xi][yi] = pipvgx;
                                        im_east[xi][yi]  = pipvgy;
                                        im_elev[xi][yi]  = -pipvgz;
                                        im_refx[xi][yi]  = (int32_T) c;
                                        im_refy[xi][yi]  = (int32_T) l;
                                        im_refs[xi][yi]  = (int32_T) sxy;
                                        //printf("isinFOVd=%d\n",isinFOVd);

                                        // evaluate nearest neighbor
                                        //dst_ppv[0] = pow(pipvgx-ppv1gx,2)+pow(pipvgy-ppv1gy,2)+pow(pipvgz-ppv1gz,2);
                                        //dst_ppv[1] = pow(pipvgx-ppv2gx,2)+pow(pipvgy-ppv2gy,2)+pow(pipvgz-ppv2gz,2);
                                        //dst_ppv[2] = pow(pipvgx-ppv3gx,2)+pow(pipvgy-ppv3gy,2)+pow(pipvgz-ppv3gz,2);
                                        //dst_ppv[3] = pow(pipvgx-ppv4gx,2)+pow(pipvgy-ppv4gy,2)+pow(pipvgz-ppv4gz,2);
                                        dst_ppv[0] = get_sqr_dst(pipvgx,pipvgy,pipvgz,ppv1gx,ppv1gy,ppv1gz);
                                        dst_ppv[1] = get_sqr_dst(pipvgx,pipvgy,pipvgz,ppv2gx,ppv2gy,ppv2gz);
                                        dst_ppv[2] = get_sqr_dst(pipvgx,pipvgy,pipvgz,ppv3gx,ppv3gy,ppv3gz);
                                        dst_ppv[3] = get_sqr_dst(pipvgx,pipvgy,pipvgz,ppv4gx,ppv4gy,ppv4gz);
                                        //printf("isinFOVd=%d\n",isinFOVd);
                                        dst_nn = dst_ppv[0]; di_dst_min = 0;
                                        for(di=1;di<4;di++){
                                            if(dst_nn > dst_ppv[di]){
                                                dst_nn = dst_ppv[di];
                                                di_dst_min = di;
                                            }
                                        }
                                        im_nnx[xi][yi] = ppvx[di_dst_min];
                                        im_nny[xi][yi] = ppvy[di_dst_min];

                                    }
                                }
                                
                            }
                        }
                    }
                    else // if(isinFOV)
                    {
                        // If isinFOV but not in FOVd
                        for(xi=x_min;xi<S_im;xi++){
                            for(yi=y_min;yi<L_im;yi++){
                                // line plane intersect
                                
                                // evaluate direction
                                
                                // plane intersection pointing vector
                                
                        
                            }
                        }
                        
                    }
                        
                }
            
            }
            
        }
        
        
        
        

    }
    free(elevl);
    free(elevlp1);
    fclose(fid);
}

double** set_mxDoubleMatrix(const mxArray *pmi){
    mwSize M,N;
    mwIndex j;
    double **pm;
    M = mxGetM(pmi); N = mxGetN(pmi);
    pm = (double **) mxMalloc(N*sizeof(double*));
    pm[0] = mxGetDoubles(pmi);
    for(j=1;j<N;j++){
        pm[j] = pm[j-1]+M;
    }
    return pm;
}

bool** set_mxLogicalMatrix(const mxArray *pmi){
    mwSize M,N;
    mwIndex j;
    bool **pm;
    M = mxGetM(pmi); N = mxGetN(pmi);
    pm = (bool **) mxMalloc(N*sizeof(bool*));
    pm[0] = mxGetLogicals(pmi);
    for(j=1;j<N;j++){
        pm[j] = pm[j-1]+M;
    }
    return pm;
}

int32_T** set_mxInt32Matrix(const mxArray *pmi){
    mwSize M,N;
    mwIndex j;
    int32_T **pm;
    M = mxGetM(pmi); N = mxGetN(pmi);
    pm = (int32_T **) mxMalloc(N*sizeof(int32_T*));
    pm[0] = mxGetInt32s(pmi);
    for(j=1;j<N;j++){
        pm[j] = pm[j-1]+M;
    }
    return pm;
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *msldem_imgpath;
    EnviHeader msldem_header;
    mwSize msldemc_imxy_sample_offset,msldemc_imxy_line_offset;
    double *msldemc_northing;
    double *msldemc_easting;
    double **msldemc_imx;
    double **msldemc_imy;
    bool **msldemc_imFOVmask;
    bool **msldemc_imFOVmaskd;
    mwSize S_im,L_im;
    double *cam_C;
    double *cam_A;
    double **PmC;
    
    double **im_north;
    double **im_east;
    double **im_elev;
    int32_T **im_refx;
    int32_T **im_refy;
    int32_T **im_refs;
    double **im_range;
    int32_T **im_nnx;
    int32_T **im_nny;
    bool **msldem_nh;
    
    mwIndex si,li;
    mxLogical *pm;

    
    mwSize msldemc_samples, msldemc_lines;

    /* -----------------------------------------------------------------
     * CHECK PROPER NUMBER OF INPUTS AND OUTPUTS
     * ----------------------------------------------------------------- */
    /* if(nrhs!=11) {
        mexErrMsgIdAndTxt("proj_mastcam2MSLDEM_v4_mex:nrhs","Eleven inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("proj_mastcam2MSLDEM_v4_mex:nlhs","Five outputs required.");
    }
    */
    /* make sure the first input argument is scalar */
    /*
    if( !mxIsChar(prhs[0]) ) {
        mexErrMsgIdAndTxt("proj_mastcam2MSLDEM_v4_mex:notChar","Input 0 needs to be a character vector.");
    }
    */
    /* -----------------------------------------------------------------
     * I/O SETUPs
     * ----------------------------------------------------------------- */
    
    /* INPUT 0 msldem_imgpath */
    msldem_imgpath = mxArrayToString(prhs[0]);
    
    /* INPUT 1 msldem_header */
    msldem_header.samples = (int) mxGetScalar(mxGetField(prhs[1],0,"samples"));
    msldem_header.lines = (int) mxGetScalar(mxGetField(prhs[1],0,"lines"));
    msldem_header.bands = (int) mxGetScalar(mxGetField(prhs[1],0,"bands"));
    msldem_header.data_type = (int) mxGetScalar(mxGetField(prhs[1],0,"data_type"));
    msldem_header.byte_order = (int) mxGetScalar(mxGetField(prhs[1],0,"byte_order"));
    msldem_header.header_offset = (int) mxGetScalar(mxGetField(prhs[1],0,"header_offset"));
    
    /* INPUT 2 msldemc_sheader*/
    msldemc_imxy_sample_offset = (mwSize) mxGetScalar(mxGetField(prhs[2],0,"sample_offset"));
    msldemc_imxy_line_offset = (mwSize) mxGetScalar(mxGetField(prhs[2],0,"line_offset"));
    msldemc_samples = (mwSize) mxGetScalar(mxGetField(prhs[2],0,"samples"));
    msldemc_lines = (mwSize) mxGetScalar(mxGetField(prhs[2],0,"lines"));
    
    /* INPUT 3 msldem northing easting */
    msldemc_northing = mxGetDoubles(prhs[3]);
    msldemc_easting = mxGetDoubles(prhs[4]);
    
    /* INPUT 5/6 msldem imxy */
    msldemc_imx = set_mxDoubleMatrix(prhs[5]);    
    msldemc_imy = set_mxDoubleMatrix(prhs[6]);
    
    /* INPUT 7 msldem imFOV */
    msldemc_imFOVmask = set_mxLogicalMatrix(prhs[7]);

    /* INPUT 8 msldem imFOV */
    msldemc_imFOVmaskd = set_mxLogicalMatrix(prhs[8]);
    
    /* INPUT 9/10 image S_im, L_im */
    S_im = (mwSize) mxGetScalar(prhs[9]);
    L_im = (mwSize) mxGetScalar(prhs[10]);
    
    /* INPUT 11/12 camera model */
    cam_C = mxGetDoubles(prhs[11]);
    cam_A = mxGetDoubles(prhs[12]);
    
    /* INPUT 13 PmC */
    PmC = set_mxDoubleMatrix(prhs[13]);
    
    /* OUTPUT 0/1/2 north-east-elevation */
    plhs[0] = mxCreateDoubleMatrix(L_im,S_im,mxREAL);
    im_north = set_mxDoubleMatrix(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(L_im,S_im,mxREAL);
    im_east = set_mxDoubleMatrix(plhs[1]);
    
    plhs[2] = mxCreateDoubleMatrix(L_im,S_im,mxREAL);
    im_elev = set_mxDoubleMatrix(plhs[2]);
    
    /* OUTPUT 3/4/5 im_ref x/y/zs */
    plhs[3] = mxCreateNumericMatrix(L_im,S_im,mxINT32_CLASS,mxREAL);
    im_refx = set_mxInt32Matrix(plhs[3]);
    
    plhs[4] = mxCreateNumericMatrix(L_im,S_im,mxINT32_CLASS,mxREAL);
    im_refy = set_mxInt32Matrix(plhs[4]);
    
    plhs[5] = mxCreateNumericMatrix(L_im,S_im,mxINT32_CLASS,mxREAL);
    im_refs = set_mxInt32Matrix(plhs[5]);
    
    /* OUTPUT 6 im_range */
    plhs[6] = mxCreateDoubleMatrix(L_im,S_im,mxREAL);
    im_range = set_mxDoubleMatrix(plhs[6]);
    
    /* OUTPUT 7/8 im_nn x/y */
    plhs[7] = mxCreateNumericMatrix(L_im,S_im,mxINT32_CLASS,mxREAL);
    im_nnx = set_mxInt32Matrix(plhs[7]);
    
    plhs[8] = mxCreateNumericMatrix(L_im,S_im,mxINT32_CLASS,mxREAL);
    im_nny = set_mxInt32Matrix(plhs[8]);
    
    // Initialize matrices
    for(si=0;si<S_im;si++){
        for(li=0;li<L_im;li++){
            im_north[si][li] = NAN;
            im_east[si][li] = NAN;
            im_elev[si][li] = NAN;
            im_range[si][li] = INFINITY;
            im_refx[si][li] = -1;
            im_refy[si][li] = -1;
            im_refs[si][li] = -1;
            im_nnx[si][li] = -1;
            im_nny[si][li] = -1;
            
        }
    }
    
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    proj_mastcam2MSLDEM(msldem_imgpath, msldem_header, 
       msldemc_imxy_sample_offset, msldemc_imxy_line_offset,
       msldemc_samples, msldemc_lines,
       msldemc_northing, msldemc_easting,
       msldemc_imx, msldemc_imy,
       msldemc_imFOVmask,msldemc_imFOVmaskd, (int32_T) S_im, (int32_T) L_im, cam_C, cam_A,
       PmC, im_north, im_east, im_elev,
       im_refx, im_refy, im_refs,
       im_range, im_nnx, im_nny);
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(msldemc_imx);
    mxFree(msldemc_imy);
    mxFree(msldemc_imFOVmask);
    mxFree(msldemc_imFOVmaskd);
    mxFree(PmC);
    mxFree(im_north);
    mxFree(im_east);
    mxFree(im_elev);
    mxFree(im_refx);
    mxFree(im_refy);
    mxFree(im_refs);
    mxFree(im_range);
    mxFree(im_nnx);
    mxFree(im_nny);
    
}
