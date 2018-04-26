#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double absx(double x);

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double dt;
    double *V, *s, *g, *gL, *C, *Vthr, *Vres, *VL, *Vsyn, *tauR, *tauD;
    double *Idc, *It, *LFP;
    int n, TT, t, i, j, *nspk;
    
    double iisyn;
    
    double *Io, *Isyn, *IL, *Isti, *sr, *ss, *VV, *gg, *IIt, *spk;
    double *Vt;
    mxArray *sspk, *spktmp;
    double *spkptr;
    
    srand(clock());
    
    n = (int)(mxGetScalar(prhs[0]));
    dt = mxGetScalar(prhs[1]);
    TT = (int)(mxGetScalar(prhs[2]));
    
    V = mxGetPr(prhs[3]);
    s = mxGetPr(prhs[4]);
    C = mxGetPr(prhs[5]);
    VL = mxGetPr(prhs[6]);
    Vthr = mxGetPr(prhs[7]);
    Vres = mxGetPr(prhs[8]);
    Vsyn = mxGetPr(prhs[9]);
    gL = mxGetPr(prhs[10]);
    g  = mxGetPr(prhs[11]);
    tauR = mxGetPr(prhs[12]);
    tauD = mxGetPr(prhs[13]);
    Idc  = mxGetPr(prhs[14]);
    It   = mxGetPr(prhs[15]);
    
    
    plhs[1] = mxCreateDoubleMatrix(n,TT, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,TT, mxREAL);
    Vt = mxGetPr(plhs[1]);
    LFP = mxGetPr(plhs[2]);
    
    
    plhs[0] = mxCreateCellMatrix(1,n);
    sspk = plhs[0];
    
    Io = mxMalloc(n * sizeof(double));
    Isyn = mxMalloc(n* sizeof(double));
    IL = mxMalloc(n* sizeof(double));
    sr = mxMalloc(n* sizeof(double));
    ss = mxMalloc(n* sizeof(double));
    VV = mxMalloc(n* sizeof(double));
    Isti = mxMalloc(n* sizeof(double));
    spk = mxMalloc(n*TT* sizeof(double));
    nspk = mxMalloc(n* sizeof(int));
    
    
    for (i=0; i<n; i++){
        nspk[i] = 0;
        sr[i] = 0.;
    }
    
    IIt = It;
    
    
    for (t=0; t<TT; t++) {
        
        *LFP = 0;
        
        
        for (i=0; i<n; i++) {
            ss[i] = s[i] + dt* (sr[i]- s[i])/tauD[i];
                        sr[i] *=  1-dt/tauR[i];
        }
        
        gg = g;
        
        for (i=0; i<n; i++) {
            IL[i] = gL[i] * (VL[i] - V[i]);
            Isyn[i] = 0.;
            
            for (j=0; j<n; j++) {
                iisyn =  *gg  * ss[j] * (Vsyn[j]-V[i]);
                Isyn[i] += iisyn;
                if (Vsyn[i]>=0.)
                    *LFP -= absx(iisyn);
                gg++;
            }
            
            VV[i] = V[i] + dt*(IL[i] + Idc[i] + *IIt + Isyn[i]) / C[i];
            
            
            V[i] = VV[i];
            s[i] = ss[i];
            
            *Vt = V[i];
            Vt++;
            IIt++;
            
            if (V[i] > Vthr[i]) {
                V[i] = Vres[i];
                spk[TT*i+nspk[i]] = t*dt;
                nspk[i]++;
                sr[i] += 1.;
                
            }
            
        }
        
        LFP++;
        
        
    }
    
    
    for (i=0; i<n; i++) {
        
        
        spktmp = mxCreateDoubleMatrix(1,nspk[i], mxREAL);
        spkptr = mxGetPr(spktmp);
        for (j=0; j<nspk[i]; j++)
            spkptr[j] = spk[TT*i+j];
        
        mxSetCell(sspk, i, spktmp);
        
    }
    
    mxFree;
    
}


#define RANDOM_MAX					2.1475e+09

double absx(double x)
{
    if (x>0)
        return x;
    else return -x;
    }
