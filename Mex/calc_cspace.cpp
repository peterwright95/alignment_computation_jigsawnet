#define _USE_MATH_DEFINES
#pragma warning(disable:4503)
// Generic mex includes
#include <mex.h>
#include <matrix.h>
#include <sstream>
#include <exception>
#include <vector>
#include <algorithm>

// Specific mex includes
#include <ppl.h>
#include <iostream>
#include <CGAL/Cartesian.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/Small_side_angle_bisector_decomposition_2.h>
#include <CGAL/Polygon_convex_decomposition_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                                   Point;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                Polygon_with_holes_2;
typedef CGAL::Aff_transformation_2<Kernel> Transformation;

struct InputArgs
{
	std::vector<std::vector<std::pair<double,double>>> P;
	std::vector<std::vector<std::pair<double,double>>> Q;
	const double * pairs = nullptr;
	size_t n_pairs = 0;
	const double * rots = nullptr;
	size_t n_rots = 0;
};

struct OutputArgs
{
	std::vector<std::vector<std::vector<double>>> cs;
};

void ParseInput(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], InputArgs& _inArgs, OutputArgs& _outArgs);
void CheckInputValidity(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void HandleArgs(const InputArgs& _inArgs, OutputArgs& _outArgs);
void CalcCS(Polygon_2 & _poly1, Polygon_2 & _poly2, double angle, std::vector<double> & _cs);
void mxArray2Poly(const mxArray * _array, std::vector<std::pair<double, double>> & _poly);
void GenerateOutput(int nlhs, mxArray *plhs[], OutputArgs& _outArgs);
void CheckPolygon(Polygon_2 & _poly, size_t _pi, const std::string& _polyArgStr);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	try
	{
		InputArgs inArgs;
		OutputArgs outArgs;
		ParseInput(nlhs, plhs, nrhs, prhs, inArgs, outArgs);
		HandleArgs(inArgs, outArgs);
		GenerateOutput(nlhs, plhs, outArgs);
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
	const mxArray * P = prhs[argi];
	if (mxIsEmpty(P))
	{
		mexErrMsgIdAndTxt("MATLAB:invalidInput", "First argument was empty");
	}
	size_t n_P = std::max<size_t>(mxGetM(P), mxGetN(P));
	++argi;
	const mxArray * Q = prhs[argi];
	size_t n_Q = 0;
	if (mxIsEmpty(Q))
	{
		Q = P;
		n_Q = n_P;
	}
	else
	{
		n_Q = std::max<size_t>(mxGetM(Q), mxGetN(Q));
	}
	++argi;
	_inArgs.pairs = (double*)mxGetData(prhs[argi]);
	size_t pairsM = mxGetM(prhs[argi]);
	_inArgs.n_pairs = mxGetN(prhs[argi]);
	if (pairsM != 2)
	{
		std::stringstream ss;
		ss << "[ParseInput]: " << "Third argument should be of size 2 x N. Instead was: "<< pairsM<< " x "<< _inArgs.n_pairs;
		mexErrMsgIdAndTxt("MATLAB:invalidInput", ss.str().c_str());
	}
	++argi;
	_inArgs.rots = (double*)mxGetData(prhs[argi]);
	_inArgs.n_rots = std::max<size_t>(mxGetM(prhs[argi]), mxGetN(prhs[argi]));

	_inArgs.Q.resize(n_Q);
	for (size_t qi = 0; qi < n_Q; ++qi)
	{
		auto& poly = _inArgs.Q[qi];
		mxArray2Poly(mxGetCell(Q, qi), poly);
	}
	
	_inArgs.P.resize(n_P);
	for (size_t pi = 0; pi < n_P; ++pi)
	{
		auto& poly = _inArgs.P[pi];
		mxArray2Poly(mxGetCell(P, pi), poly);
	}

	_outArgs.cs.resize(_inArgs.n_pairs, std::vector<std::vector<double>>(_inArgs.n_rots));
}

void CheckInputValidity(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int expectedNrhs = 4;
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

	std::vector<mxClassID> expectedInArgsClasses = { mxCELL_CLASS,mxCELL_CLASS,mxDOUBLE_CLASS, mxDOUBLE_CLASS };
	std::vector<size_t> expectedInArgsDims = { 2,2,2, 2 };
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
	//Concurrency::parallel_for((size_t)0, _inArgs.n_pairs, [&](size_t pi)
	for (size_t pi = 0; pi < _inArgs.n_pairs; ++pi)
	{
		size_t pair_i[] = { static_cast<size_t>(_inArgs.pairs[2 * pi] - 1), static_cast<size_t>(_inArgs.pairs[2 * pi + 1] - 1) };
		Concurrency::parallel_for((size_t)0, _inArgs.n_rots, [&](size_t ri)
		//for (size_t ri = 0; ri < _inArgs.n_rots; ++ri)
		{
			Polygon_2 poly1;
			for (auto& point : _inArgs.P[pair_i[0]])
			{
				poly1.push_back(Point(point.first, point.second));
			}
			CheckPolygon(poly1, pair_i[0], "P");
			Polygon_2 poly2;
			for (auto& point : _inArgs.Q[pair_i[1]])
			{
				poly2.push_back(Point(point.first, point.second));
			}
			CheckPolygon(poly2, pair_i[1],"Q");
			double rot = _inArgs.rots[ri];
			CalcCS(poly1, poly2, rot, _outArgs.cs[pi][ri]);
		});
	}//);
}

void CalcCS(Polygon_2 & _poly1, Polygon_2 & _poly2, double angle, std::vector<double> & _cs)
{
	_cs.clear();
	// for cs=P-RqQ; where P is the first polygon and Qis the second polygon
	// and Rq is the rotation of the second polygon. we add 180 to turn minkowski sum into minkowski diffrence.
	// and take -angle instead of angle because the polygon is upside-down
	// (from conversion of matlab image coordinates to regular coordinates)
	// so while imrotate rotates ccw we rotate cw

	double rad_angle = (180-angle) * M_PI / 180;
	Polygon_2 poly2R = CGAL::transform(Transformation(CGAL::ROTATION, sin(rad_angle), cos(rad_angle)), _poly2);
	Polygon_with_holes_2 sum;
	try
	{
		sum = CGAL::minkowski_sum_by_reduced_convolution_2(_poly1, poly2R);
	}
	catch (...)
	{
		return;
		//CGAL::Small_side_angle_bisector_decomposition_2<Kernel>  ssab_decomp;
		//sum = CGAL::minkowski_sum_2(_poly1, poly2R, ssab_decomp);
	}
	auto& ob = sum.outer_boundary();
	_cs.reserve(ob.size()*2);
	for (int i = 0; i<ob.size(); ++i)
	{
		_cs.push_back(CGAL::to_double(ob[i].x().exact()));
		_cs.push_back(CGAL::to_double(ob[i].y().exact()));
	}
}

void mxArray2Poly(const mxArray * _array, std::vector<std::pair<double, double>> & _poly)
{
	if (_array == nullptr)
	{
		mexErrMsgIdAndTxt("MATLAB:ArgumentNullException", "could not find polygon inside cell array");
	}
	size_t M = mxGetM(_array);
	size_t N = mxGetN(_array);
	if (M != 2)
	{
		std::stringstream ss;
		ss << "[mxArray2Poly]: " << "polygon inside cell should be of size 2 x N. Instead was: " << M << " x " << N;
		mexErrMsgIdAndTxt("MATLAB:invalidInput", ss.str().c_str());
	}
	const double * polydata = (double*)mxGetData(_array);
	_poly.clear();
	_poly.reserve(N);
	for (int i = 0; i < N; ++i)
	{
		_poly.emplace_back(polydata[2 * i], polydata[2 * i + 1]);
	}
}

void GenerateOutput(int nlhs, mxArray *plhs[], OutputArgs& _outArgs)
{
	auto& allcs = _outArgs.cs;
	size_t dims[] = { allcs[0].size(), allcs.size() };
	plhs[0] = mxCreateCellArray(2, dims);
	size_t csIdx = 0;
	for (auto& prvec : allcs)
	{
		for (auto& pqcs : prvec)
		{
			if (pqcs.size() != 0)
			{
				mxArray * pqCsArray = mxCreateDoubleMatrix(2, pqcs.size()/2, mxREAL);
				double * pr = mxGetPr(pqCsArray);
				std::copy(pqcs.begin(), pqcs.end(), pr);
				mxSetCell(plhs[0], csIdx, pqCsArray);
			}
			++csIdx;
		}
	}
}

void CheckPolygon(Polygon_2 & _poly, size_t _pi, const std::string& _polyArgStr)
{
	if (!_poly.is_simple())
	{
		std::stringstream ss;
		ss << "[HandleArgs]: " << "Polygon #" << _pi + 1 << " in " << _polyArgStr << " polygons set is not simple.";
		throw std::runtime_error(ss.str());
	}
	if (!_poly.is_counterclockwise_oriented())
	{
		std::stringstream ss;
		ss << "[HandleArgs]: " << "Polygon #" << _pi + 1 << " in " << _polyArgStr << " polygons set is not ccw oriented.";
		throw std::runtime_error(ss.str());
	}
}