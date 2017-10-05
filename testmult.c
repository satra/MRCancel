#include "mex.h"
#include <math.h>
#include <float.h>

#define MYPOW(x)    (x*x)

static int initialized = 0;

static mxArray *corr_out     = NULL;
static double  *corr_out_ptr = NULL;

void cleanup(void) {
	mexPrintf("MEX-file is terminating, destroying array\n");
	mxDestroyArray(corr_out);
}

void mexFunction(int nlhs,
				 mxArray *plhs[],
				 int nrhs,
				 const mxArray *prhs[])
{
	int i0,j0;
	double *buf_ptr,*template_ptr;
    int buf_start, framelen, templen;
    double frame_sum = 0.0, frame_norm=0.0, temp_norm=0.0;
    
	buf_ptr = (double*)mxGetPr(prhs[0]);
	template_ptr = (double*)mxGetPr(prhs[1]);
    buf_start = (int)(*mxGetPr(prhs[2]));
    framelen  = (int)(*mxGetPr(prhs[3]));
    templen   = mxGetN(prhs[1]);
    
	if(!initialized) {
		mexPrintf("MEX-file initializing, creating arrays\n");

		/* Create persistent array and register its cleanup. */
		corr_out = mxCreateDoubleMatrix(1,framelen, mxREAL);
		mexMakeArrayPersistent(corr_out);
		mexAtExit(cleanup);
		initialized = 1;

		corr_out_ptr = mxGetPr(corr_out);
	}
    
    // norm of template
    temp_norm = 0.0;
    for(j0=0;j0<templen;j0++)
    {
        temp_norm += MYPOW(template_ptr[j0]);
    }
    temp_norm = sqrt(temp_norm);
    //printf("%f\n",temp_norm);
    //printf("%d\n",buf_start);
    
    // cross correlation
    for(i0=0;i0<framelen;i0++)
    {
        frame_sum = 0.0;frame_norm=0.0;
        for(j0=0;j0<templen;j0++)
        {
            frame_sum += template_ptr[j0]*buf_ptr[buf_start-1+i0+j0];
            frame_norm += MYPOW(buf_ptr[buf_start-1+i0+j0]);
        }
        corr_out_ptr[i0] = frame_sum/(sqrt(frame_norm)*temp_norm); 
    }
    
    plhs[0] = corr_out;
    return;
}