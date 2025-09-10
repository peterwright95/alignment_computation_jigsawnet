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
	const double * p = nullptr;
	size_t n;
	const double * node = nullptr;
	const double * edge = nullptr;
	size_t nc;
	double TOL = 1.0e-12;
	double nan = -1;
	double eps = -1;
};

struct OutputArgs
{
	bool * in = nullptr;
	bool * on = nullptr;
};

void ParseInput(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], InputArgs& _inArgs, OutputArgs& _outArgs);
void CheckInputValidity(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void HandleArgs(const InputArgs& _inArgs, OutputArgs& _outArgs);
template <typename T>
void sort_indexes(T* _v, size_t _vsz, std::vector<size_t>& _idx /*out*/);

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
	_inArgs.p = (double*)mxGetData(prhs[argi]);
	size_t p_c = (size_t)mxGetN(prhs[argi]);
	_inArgs.n = (size_t)mxGetM(prhs[argi++]);
	if (p_c != 2)
	{
		std::stringstream ss;
		ss << "[ParseInput] - p must be an Nx2 array, instead was: " << _inArgs.n<< " x " << p_c;
		throw std::runtime_error(ss.str().c_str());
	}
	
	_inArgs.node = (double*)mxGetData(prhs[argi]);
	size_t node_c = (size_t)mxGetN(prhs[argi]);
	_inArgs.nc = (size_t)mxGetM(prhs[argi++]);
	if (node_c != 2)
	{
		std::stringstream ss;
		ss << "[ParseInput] - NODE must be an Mx2 array, instead was: " << _inArgs.nc << " x " << node_c;
		throw std::runtime_error(ss.str().c_str());
	}
	if (argi < nrhs)
	{
		_inArgs.edge = (double*)mxGetData(prhs[argi]);
		size_t edge_r = (size_t)mxGetM(prhs[argi]);
		size_t edge_c = (size_t)mxGetN(prhs[argi++]);
		if (edge_r != _inArgs.nc || edge_c != 2)
		{
			std::stringstream ss;
			ss << "[ParseInput] - EDGE must be an " << _inArgs.nc <<"x2 array. instead was: " << edge_r << " x " << edge_c;
			throw std::runtime_error(ss.str().c_str());
		}
		auto edges_minmax = std::minmax_element(_inArgs.edge, _inArgs.edge + _inArgs.nc * 2);
		if (*edges_minmax.first < 1 || *edges_minmax.second > _inArgs.nc)
		{
			std::stringstream ss;
			ss << "[ParseInput] - EDGE must contain indexes within the range [1, " << _inArgs.nc << "]. instead the range was [" << edges_minmax.first << ", " << edges_minmax.second <<"].";
			throw std::runtime_error(ss.str().c_str());
		}
		if (argi < nrhs)
		{
			_inArgs.TOL = (double)mxGetScalar(prhs[argi++]);
		}
	}
	
	_inArgs.nan = mxGetNaN();
	_inArgs.eps = mxGetEps();

	if (nlhs > 0)
	{
		plhs[0] = mxCreateLogicalMatrix(_inArgs.n, 1);
		_outArgs.in = (bool *)mxGetData(plhs[0]);
		if (nlhs > 1)
		{
			plhs[1] = mxCreateLogicalMatrix(_inArgs.n, 1);
			_outArgs.on = (bool *)mxGetData(plhs[1]);
		}
	}

	/*if (_outArgs.on == nullptr)
	{
		auto onArr = mxCreateLogicalMatrix(_inArgs.n, 1);
		_outArgs.on = (bool *)mxGetData(onArr);
	}*/
}

void CheckInputValidity(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int expectedNrhs = 2;
	int expectedNlhs = 2;
	if (nrhs < expectedNrhs)
	{
		std::stringstream ss;
		ss << "[ParseInput]: at least" << expectedNrhs << " input arguments required, instead " << nrhs << " arguments were recieved.";
		throw std::runtime_error(ss.str().c_str());
	}
	if (nlhs > expectedNlhs)
	{
		std::stringstream ss;
		ss << "[ParseInput]: at most" << expectedNlhs << " output arguments required, instead " << nlhs << " arguments were recieved.";
		throw std::runtime_error(ss.str().c_str());
	}

	std::vector<mxClassID> expectedInArgsClasses = { mxDOUBLE_CLASS,mxDOUBLE_CLASS,mxDOUBLE_CLASS,mxDOUBLE_CLASS };
	std::vector<size_t> expectedInArgsDims = { 2,2,2,2 };
	for (int i = 0; i < nrhs; ++i)
	{
		if (mxGetClassID(prhs[i]) != expectedInArgsClasses[i])
		{
			std::stringstream ssError;
			mxArray* tmp = mxCreateNumericArray(0, 0, expectedInArgsClasses[i], mxREAL);
			ssError << "[ParseInput]: Argument no. " << i + 1 << " should be of type '" << mxGetClassName(tmp) << "'. Instead was '" << mxGetClassName(prhs[i]) << "'";
			throw std::runtime_error(ssError.str().c_str());
		}
		size_t ndim = mxGetNumberOfDimensions(prhs[i]);
		if (ndim != expectedInArgsDims[i])
		{
			std::stringstream ssError;
			ssError << "[ParseInput]: Argument no. " << i + 1 << " is a " << ndim << "D matrix, instead of a " << expectedInArgsDims[i] << "D matrix";
			throw std::runtime_error(ssError.str().c_str());
		}
	}
}

void HandleArgs(const InputArgs& _inArgs, OutputArgs& _outArgs)
{
	// Choose the direction with the biggest range as the "y-coordinate" for the
	// test.This should ensure that the sorting is done along the best
	// direction for long and skinny problems wrt either the x or y axes.
	auto p_x_minmax = std::minmax_element(_inArgs.p, _inArgs.p + _inArgs.n);
	auto p_y_minmax = std::minmax_element(_inArgs.p + _inArgs.n, _inArgs.p + _inArgs.n*2);
	double dxy[2] = { *p_x_minmax.second - *p_x_minmax.first, *p_y_minmax.second - *p_y_minmax.first };
	double tol = _inArgs.TOL * std::min(dxy[0], dxy[1]);
	bool isflipped = dxy[0]>dxy[1];
	const double* px = isflipped ? (_inArgs.p + _inArgs.n): _inArgs.p;
	const double* py = isflipped ? _inArgs.p : (_inArgs.p + _inArgs.n);
	const double* nodex = isflipped ? (_inArgs.node + _inArgs.nc) : _inArgs.node;
	const double* nodey = isflipped ? _inArgs.node : (_inArgs.node + _inArgs.nc);
	std::vector<size_t> sp_i(_inArgs.n);
	sort_indexes(py, _inArgs.n, sp_i);
	
	for (size_t k = 0; k < _inArgs.nc; ++k)
	{
		//Nodes in current edge
		size_t n1 = k;
		size_t n2 = k + 1 == _inArgs.nc ? 0 : k + 1;
		if (_inArgs.edge != nullptr)
		{
			n1 = static_cast<size_t>(_inArgs.edge[k]) - 1;
			n2 = static_cast<size_t>(_inArgs.edge[k+_inArgs.nc]) - 1;
		}
		
		//Endpoints - sorted so that[x1, y1] & [x2, y2] has y1 <= y2
		//          - also get xmin = min(x1, x2), xmax = max(x1, x2)

		double y1 = nodey[n1];
		double y2 = nodey[n2];
		double x1 = nodex[n1];
		double x2 = nodex[n2];
		if (y1 >= y2)
		{
			std::swap(y1, y2);
			std::swap(x1, x2);
		}
		auto xminmax = std::minmax(x1, x2);
		double xmin = xminmax.first;
		double xmax = xminmax.second;

		// find the first point that lies beyond the lower y of the current edge
		auto start_itr = std::lower_bound(sp_i.begin(), sp_i.end(), y1, [&py](const size_t& _idx, const double & _val) { return py[_idx] < _val; });
		if (start_itr != sp_i.end())
		{
			//Loop through points
			for (auto j_itr = start_itr; j_itr != sp_i.end(); ++j_itr)
			{
				//Check the bounding-box for the edge before doing the intersection
				size_t j_idx = *j_itr;
				double Y = py[j_idx];
				if (Y <= y2)
				{
					double X = px[j_idx];
					if (X >= xmin)
					{
						if (X <= xmax)
						{
							if (_outArgs.on != nullptr)
							{
								// Check if we're "on" the edge
								_outArgs.on[j_idx] = _outArgs.on[j_idx] || std::abs((y2 - Y)*(x1 - X) - (y1 - Y)*(x2 - X)) < tol;
							}							
							
							/*if (_outArgs.on[j_idx])
							{
								_outArgs.in[j_idx] = true;
							}
							else*/
							{
								// Do the actual intersection test
								if ((Y < y2) && ((y2 - y1)*(X - x1) < (Y - y1)*(x2 - x1)))
								{
									if (_outArgs.in != nullptr)
									{
										_outArgs.in[j_idx] = !_outArgs.in[j_idx];
									}
								}
							}							
						}
					}
					else if (Y < y2) //&& !_outArgs.on[j_idx])) //Deal with points exactly at vertices
					{
						if (_outArgs.in != nullptr)
						{
							//Has to cross edge
							_outArgs.in[j_idx] = !_outArgs.in[j_idx];
						}
					}
				}
				else
				{
					// Due to the sorting, no points with >y
					// value need to be checked
					break;
				}
			}
		}
	}
}

template <typename T>
void sort_indexes(T* _v, size_t _vsz, std::vector<size_t>& _idx /*out*/) 
{
	_idx.resize(_vsz);
	std::iota(_idx.begin(), _idx.end(), 0);
	std::sort(_idx.begin(), _idx.end(),[&_v](const size_t& i1, const size_t& i2) { return _v[i1] < _v[i2]; });
}