#include "mex.h"
#include <math.h>
#include <float.h>

#define MYPOW(x)    (x*x)
static int initialized = 0;
static double oldmaxc;
static int oldidx;

void mexFunction(int nlhs,
				 mxArray *plhs[],
				 int nrhs,
				 const mxArray *prhs[])
{
	int i0, j0, buffer_len, maxidx, win_start, win_end, framelen, extraspace;
	float maxc,norm1,norm2,rms, c;
	double *buffer_ptr,*template_ptr = NULL, *window_ptr, *index_ptr;
    double *buf2_ptr, rms_val, cval1, cval2;
    
	buffer_ptr  = mxGetPr(prhs[0]);
	buffer_len  = mxGetM(prhs[0]);
	buf2_ptr    = mxGetPr(prhs[1]);
    window_ptr  = mxGetPr(prhs[2]);
    rms_val     = *mxGetPr(prhs[3]);
    cval1       = *mxGetPr(prhs[4]);
    cval2       = *mxGetPr(prhs[5]);
    framelen    = (int)(*mxGetPr(prhs[6]));
    extraspace  =  (int)(*mxGetPr(prhs[7]));
    
    if (!initialized)
    {
        initialized = 1;
        oldmaxc = 0;
        oldidx = 1;
    }
	//mexPrintf("Length of input data: [%d] Sampling Rate[%d]\n",framelen,sampling_rate);
    //buf2 = mxCalloc(buffer_len,sizeof(double));
    
    win_start  =  (int)window_ptr[0];
    win_end  =  (int)window_ptr[1];
	maxc = -1;
	maxidx = 0;
	rms = 0.0f;

	for(i0=buffer_len-2*win_end;i0<buffer_len;i0++)
		rms += buf2_ptr[i0];

	rms	= sqrt(rms/(2*win_end));
	// mexPrintf("rms[%f]ws[%d]we[%d]len[%d]\n",rms,win_start,win_end,buffer_len);		

	if (rms>rms_val)
	{
		for(i0=win_start;i0<=win_end;i0++)
		{
			norm1 = 0.0f;norm2 = 0.0f;c=0.0f;
			for(j0=0;j0<i0;j0++)
			{
				norm1 += buf2_ptr[buffer_len-2*i0+j0];
				norm2 += buf2_ptr[buffer_len-i0+j0];
				c     += buffer_ptr[buffer_len-2*i0+j0]*buffer_ptr[buffer_len-i0+j0];
			}
			if ((fabs(norm1)>FLT_EPSILON) && (fabs(norm2)>FLT_EPSILON))
				c /= sqrt(norm1)*sqrt(norm2);
			else
				c = -1;
			if (c>maxc){
				maxc = c;
				maxidx = i0;
			}
		}
		mexPrintf("maxc[%f]\n",maxc);		
        if (maxc>cval1)
        {
			plhs[0] = mxCreateDoubleMatrix(maxidx+extraspace, 1, mxREAL);
			template_ptr = mxGetPr(plhs[0]);
			for(i0=0;i0<(oldidx+extraspace);i0++)
                template_ptr[i0] = buffer_ptr[buffer_len-2*maxidx+i0];
            if (nlhs==2)
            {
                plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
                index_ptr = mxGetPr(plhs[1]);
                index_ptr[0] = buffer_len-2*maxidx+i0;
            }
			return;
        }
        else if ((maxc<oldmaxc) & (oldmaxc>cval2))
		{
			plhs[0] = mxCreateDoubleMatrix(oldidx+extraspace, 1, mxREAL);
			template_ptr = mxGetPr(plhs[0]);
            for(i0=0;i0<(oldidx+extraspace);i0++)
                template_ptr[i0] = buffer_ptr[buffer_len-2*oldidx+i0-framelen];
            if (nlhs==2)
            {
                plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
                index_ptr = mxGetPr(plhs[1]);
                index_ptr[0] = buffer_len-2*oldidx-framelen;
            }
            
            return;
		}
        oldmaxc = maxc;
        oldidx = maxidx;
	}
	if (template_ptr==NULL)
		plhs[0] = mxCreateDoubleMatrix(0, 1, mxREAL);
    if (nlhs==2)
        plhs[1] = mxCreateDoubleMatrix(0, 1, mxREAL);
}