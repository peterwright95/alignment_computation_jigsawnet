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
	const double * Mi;
	double nan = -1;
	double eps = -1;
};

struct OutputArgs
{
	double * ov_limits;
};

void ParseInput(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], InputArgs& _inArgs, OutputArgs& _outArgs);
void CheckInputValidity(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void HandleArgs(const InputArgs& _inArgs, OutputArgs& _outArgs);
double MatMod(const double & a, const double & b);

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
	_inArgs.Mi = (double*)mxGetData(prhs[argi]);
	size_t n_Mi = std::max({ (size_t)mxGetM(prhs[argi]),(size_t)mxGetN(prhs[argi]) });
	if (n_Mi != _inArgs.n_c)
	{
		std::stringstream ss;
		ss << "[ParseInput]: " << "Mi length should be the same as the number of contour points. Instead |Mi| = " << n_Mi << ", |n_c| = " << _inArgs.n_c;
		throw std::runtime_error(ss.str().c_str());
	}
	_inArgs.nan = mxGetNaN();
	_inArgs.eps = mxGetEps();
	plhs[0] = mxCreateDoubleMatrix(2, _inArgs.n_tforms, mxREAL);
	_outArgs.ov_limits = (double *)mxGetData(plhs[0]);
}

void CheckInputValidity(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int expectedNrhs = 2;
	int expectedNlhs = 1;
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

	std::vector<mxClassID> expectedInArgsClasses = { mxLOGICAL_CLASS,mxDOUBLE_CLASS };
	std::vector<size_t> expectedInArgsDims = { 2,2 };
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
		std::vector<size_t> vi(_inArgs.n_c);
		std::iota(vi.begin(), vi.end(), 0);
		auto& ovEndItr = std::remove_if(vi.begin(), vi.end(), [&](size_t& a) {return !tiov[a]; });
		if (vi.begin() != ovEndItr)
		{
			std::sort(vi.begin(), ovEndItr, [&](const size_t & a, const size_t & b)
			{
				return _inArgs.Mi[a] < _inArgs.Mi[b];
			});
			ovEndItr=vi.insert(ovEndItr, vi[0]);

			auto first = vi.begin();
			auto largestItr = first;
			double largest = MatMod(_inArgs.Mi[*(first + 1)] - _inArgs.Mi[*first], 1 + _inArgs.eps);
			++first;
			for (; first != ovEndItr; ++first)
			{
				double firstMod = MatMod(_inArgs.Mi[*(first + 1)] - _inArgs.Mi[*first], 1 + _inArgs.eps);
				if (firstMod > largest)
				{
					largest = firstMod;
					largestItr = first;
				}
			}
			_outArgs.ov_limits[2 * ti] = _inArgs.Mi[*(largestItr + 1)];
			_outArgs.ov_limits[2 * ti+1] = _inArgs.Mi[*(largestItr)];
		}
		else
		{
			_outArgs.ov_limits[2 * ti] = _inArgs.nan;
			_outArgs.ov_limits[2 * ti + 1] = _inArgs.nan;
		}
	});
}

inline double MatMod(const double & a, const double & b)
{
	return std::fmod(b + std::fmod(a, b), b);
}
//std::numeric_limits<T>::epsilon()