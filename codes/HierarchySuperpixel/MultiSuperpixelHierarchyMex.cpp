#include "mex.h"
#include "MultiSuperpixelHierarchyMex.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// model = superpixelHirarchyMex(image, edge, connect, iter);

	// que significa el 4?
	// if (nrhs < 2 || nrhs > 5)
		// mexErrMsgTxt("usage: model = MultisuperpixelHierarchyMex(image, edge, [connect], [iter], [selPerc], [sparseJoin])");
	// if (mxGetClassID(prhs[0]) != mxUINT8_CLASS)
	// 	mexErrMsgTxt("Input image must be uint8.");
	// if (mxGetClassID(prhs[1]) != mxUINT8_CLASS)
	// 	mexErrMsgTxt("Input edge map must be uint8.");

	unsigned short *image = (unsigned short *)mxGetData(prhs[0]);
	unsigned short *edge  = (unsigned short *)mxGetData(prhs[1]);
	const mwSize *dims;
	dims = mxGetDimensions(prhs[0]);

	int h = dims[0]; int w = dims[1];  
	int chnl = dims[2];
    int connect = mxGetScalar(prhs[2]);
	int iterSwitch = mxGetScalar(prhs[3]);
	int SelPerc = mxGetScalar(prhs[4]);
	int maxDistColor = mxGetScalar(prhs[5]);
	int maxSize = mxGetScalar(prhs[6]);
	bool euclidean = mxGetScalar(prhs[7]);
	bool meanColor = mxGetScalar(prhs[8]);

	//int connect = 24; int iterSwitch = 0;
    mexPrintf("Parameters Loaded");
    mexEvalString("drawnow;");
	
	MultiSuperpixelHierarchy SH;
	SH.init(h,w,chnl,connect,iterSwitch,SelPerc,maxDistColor,maxSize,euclidean,meanColor);
	SH.buildTree(image,edge);

	const char *fieldnames[] = {"parent","label","treeu","treev","nvertex","nregion","JoinedPixels"};
	mxArray *parent   = mxCreateNumericMatrix(h,w,mxINT32_CLASS,mxREAL);
	mxArray *label    = mxCreateNumericMatrix(h,w,mxINT32_CLASS,mxREAL);
	mxArray *treeu    = mxCreateNumericMatrix(h*w-1,1,mxINT32_CLASS,mxREAL);
	mxArray *treev    = mxCreateNumericMatrix(h*w-1,1,mxINT32_CLASS,mxREAL);
	mxArray *nvertex  = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
	mxArray *nregion  = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
	mxArray *JoinedPixels  = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
	memcpy(mxGetData(parent), SH.getParent(), sizeof(int)*h*w);
	memcpy(mxGetData(label),  SH.getLabel(),  sizeof(int)*h*w);
	memcpy(mxGetData(treeu),  SH.getTreeU(),  sizeof(int)*(h*w-1));
	memcpy(mxGetData(treev),  SH.getTreeV(),  sizeof(int)*(h*w-1));
	*(int*)mxGetData(nvertex) = h*w; *(int*)mxGetData(nregion) = SH.getRegion(); *(int*)mxGetData(JoinedPixels) = SH.getJoinedPixels();
	plhs[0] = mxCreateStructMatrix(1,1,7,fieldnames);
	mxSetFieldByNumber(plhs[0],0,0,parent);
	mxSetFieldByNumber(plhs[0],0,1,label);
	mxSetFieldByNumber(plhs[0],0,2,treeu);
	mxSetFieldByNumber(plhs[0],0,3,treev);
	mxSetFieldByNumber(plhs[0],0,4,nvertex);
	mxSetFieldByNumber(plhs[0],0,5,nregion);
	mxSetFieldByNumber(plhs[0],0,6,JoinedPixels);
}