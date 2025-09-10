#define _USE_MATH_DEFINES
#pragma warning(disable:4503)
// Generic mex includes
#include <mex.h>
#include <matrix.h>
#include <sstream>
#include <exception>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <ppl.h>

struct InputArgs
{
	const bool * ov_mask;
	size_t n_tforms;
	size_t n_c;
	const double * ov_loc;
	const double * c_bins;
	const double * ex_bins;
	size_t n_ex;
	const double * c_mag;
	const double * ex_mag;
	size_t n_bins;

	double nan = -1;
	double eps = -1;
};

struct OutputArgs
{
	double * c_hist;
	double * ex_hist;
};

void ParseInput(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], InputArgs& _inArgs, OutputArgs& _outArgs);
void CheckInputValidity(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void HandleArgs(const InputArgs& _inArgs, OutputArgs& _outArgs);
void CheckIdx(size_t _idx, const std::string& _idxName, size_t ulimit, const std::string& _uLimitName);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	try
	{
		InputArgs inArgs;
		OutputArgs outArgs;
		ParseInput(nlhs, plhs, nrhs, prhs, inArgs, outArgs);
		HandleArgs(inArgs, outArgs);
	}
	catch (std::exception ex)
	{
		mexErrMsgIdAndTxt("MATLAB:MexException", ex.what());
	}
}

void ParseInput(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], InputArgs& _inArgs, OutputArgs& _outArgs)
{
	CheckInputValidity(nlhs, plhs, nrhs, prhs);

	int argi = 0;
	_inArgs.ov_mask = (bool*)mxGetData(prhs[argi]);
	_inArgs.n_tforms = (size_t)mxGetN(prhs[argi]);
	_inArgs.n_c = (size_t)mxGetM(prhs[argi++]);
	_inArgs.ov_loc = (double*)mxGetData(prhs[argi]);
	size_t ov_locM = (size_t)mxGetM(prhs[argi]);
	size_t ov_locN = (size_t)mxGetN(prhs[argi++]);

	if (ov_locM != _inArgs.n_c || ov_locN != _inArgs.n_tforms)
	{
		std::stringstream ss;
		ss << "[ParseInput]: " << "ov_loc size should be the same as the size of ov_mask. Instead |ov_loc| = " << ov_locM<<" x "<< ov_locN
			<< ", |ov_mask| = " << _inArgs.n_c<<" x "<< _inArgs.n_tforms;
		throw std::runtime_error(ss.str().c_str());
	}
	_inArgs.c_bins = (double*)mxGetData(prhs[argi]);
	size_t c_bins_len = std::max({ (size_t)mxGetM(prhs[argi]),(size_t)mxGetN(prhs[argi]) });
	if (c_bins_len != _inArgs.n_c)
	{
		std::stringstream ss;
		ss << "[ParseInput]: " << "c_bins length should be the same as the number of contour points. Instead |c_bins| = " << c_bins_len << ", |n_c| = " << _inArgs.n_c;
		throw std::runtime_error(ss.str().c_str());
	}
	++argi;
	_inArgs.ex_bins = (double*)mxGetData(prhs[argi]);
	_inArgs.n_ex = std::max({ (size_t)mxGetM(prhs[argi]),(size_t)mxGetN(prhs[argi]) });
	++argi;
	_inArgs.c_mag = (double*)mxGetData(prhs[argi]);
	size_t c_mag_len = std::max({ (size_t)mxGetM(prhs[argi]),(size_t)mxGetN(prhs[argi]) });
	if (c_mag_len != _inArgs.n_c)
	{
		std::stringstream ss;
		ss << "[ParseInput]: " << "c_mag length should be the same as the number of contour points. Instead |c_bins| = " << c_mag_len << ", |n_c| = " << _inArgs.n_c;
		throw std::runtime_error(ss.str().c_str());
	}
	++argi;
	_inArgs.ex_mag = (double*)mxGetData(prhs[argi]);
	size_t ex_mag_len = std::max({ (size_t)mxGetM(prhs[argi]),(size_t)mxGetN(prhs[argi]) });
	if (ex_mag_len != _inArgs.n_ex)
	{
		std::stringstream ss;
		ss << "[ParseInput]: " << "ex_mag length should be the same as the number of ex points. Instead |ex_bins| = " << ex_mag_len << ", |n_ex| = " << _inArgs.n_ex;
		throw std::runtime_error(ss.str().c_str());
	}
	++argi;
	_inArgs.n_bins = (size_t)mxGetScalar(prhs[argi]);

	_inArgs.nan = mxGetNaN();
	_inArgs.eps = mxGetEps();
	plhs[0] = mxCreateDoubleMatrix(_inArgs.n_bins, _inArgs.n_tforms, mxREAL);
	_outArgs.c_hist = (double *)mxGetData(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(_inArgs.n_bins, _inArgs.n_tforms, mxREAL);
	_outArgs.ex_hist = (double *)mxGetData(plhs[1]);
}

void CheckInputValidity(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int expectedNrhs = 7;
	int expectedNlhs = 2;
	if (nrhs != expectedNrhs)
	{
		std::stringstream ss;
		ss << "[ParseInput]: " << expectedNrhs << " input arguments required, instead " << nrhs << " arguments were recieved.";
		mexErrMsgIdAndTxt("MATLAB:invalidNumInputs", ss.str().c_str());
	}
	if (nlhs != expectedNlhs)
	{
		std::stringstream ss;
		ss << "[ParseInput]: " << expectedNlhs << " output arguments required, instead " << nlhs << " arguments were recieved.";
		mexErrMsgIdAndTxt("MATLAB:invalidNumOutputs", ss.str().c_str());
	}

	std::vector<mxClassID> expectedInArgsClasses = { mxLOGICAL_CLASS,mxDOUBLE_CLASS,mxDOUBLE_CLASS,mxDOUBLE_CLASS,mxDOUBLE_CLASS,mxDOUBLE_CLASS,mxDOUBLE_CLASS };
	std::vector<size_t> expectedInArgsDims = { 2,2,2,2,2,2,2 };
	for (int i = 0; i < nrhs; ++i)
	{
		if (mxGetClassID(prhs[i]) != expectedInArgsClasses[i])
		{
			std::stringstream ssError;
			mxArray* tmp = mxCreateNumericArray(0, 0, expectedInArgsClasses[i], mxREAL);
			ssError << "[ParseInput]: Argument no. " << i + 1 << " should be of type '" << mxGetClassName(tmp) << "'. Instead was '" << mxGetClassName(prhs[i]) << "'";
			mexErrMsgIdAndTxt("MATLAB:invalidArgument", ssError.str().c_str());
		}
		size_t ndim = mxGetNumberOfDimensions(prhs[i]);
		if (ndim != expectedInArgsDims[i])
		{
			std::stringstream ssError;
			ssError << "[ParseInput]: Argument no. " << i + 1 << " is a " << ndim << "D matrix, instead of a " << expectedInArgsDims[i] << "D matrix";
			mexErrMsgIdAndTxt("MATLAB:invalidArgument", ssError.str().c_str());
		}
	}
}

void HandleArgs(const InputArgs& _inArgs, OutputArgs& _outArgs)
{

	Concurrency::parallel_for((size_t)0, _inArgs.n_tforms, [&](size_t ti)
	//for (size_t ti = 0; ti < _inArgs.n_tforms; ++ti)
	{
		const bool * tiov = &_inArgs.ov_mask[_inArgs.n_c*ti];
		const double * tiovloc = &_inArgs.ov_loc[_inArgs.n_c*ti];
		double * ti_c_hist = &_outArgs.c_hist[_inArgs.n_bins*ti];
		double * ti_ex_hist = &_outArgs.ex_hist[_inArgs.n_bins*ti];
		for (size_t ci = 0; ci < _inArgs.n_c; ++ci)
		{
			if (tiov[ci])
			{
				size_t c_binIdx = static_cast<size_t>(_inArgs.c_bins[ci] - 1);
				CheckIdx(c_binIdx, "binIdx", _inArgs.n_bins, "n_bins");
				ti_c_hist[c_binIdx] = ti_c_hist[c_binIdx] + _inArgs.c_mag[ci];

				size_t ex_px_locIdx = static_cast<size_t>(tiovloc[ci] - 1);
				CheckIdx(ex_px_locIdx, "ex_px_locIdx", _inArgs.n_ex, "n_ex");
				size_t ex_binIdx = static_cast<size_t>(_inArgs.ex_bins[ex_px_locIdx] - 1);
				CheckIdx(ex_binIdx, "ex_binIdx", _inArgs.n_bins, "n_bins");
				ti_ex_hist[ex_binIdx] = ti_ex_hist[ex_binIdx] + _inArgs.ex_mag[ex_px_locIdx];
			}
		}
	});
}

inline void CheckIdx(size_t _idx,const std::string& _idxName, size_t ulimit, const std::string& _uLimitName)
{
	if (_idx >= ulimit || _idx < 0)
	{
		std::stringstream ss;
		ss << "tried to access to non-existing bin idx. "<< _idxName<<" = " << _idx << ", "<< _uLimitName<<" = " << ulimit;
		throw std::runtime_error(ss.str().c_str());
	}
}