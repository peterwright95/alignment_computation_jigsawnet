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
#include <boost/icl/interval_set.hpp>

struct InputArgs
{
	enum LogicalFilterMode { Vec, Mat, LB };
	const double * set1 = nullptr;
	size_t n_intervals1;
	const double * set2 = nullptr;
	size_t n_intervals2;
	const bool * passFilter = nullptr;
	const double * passFilter_i = nullptr;
	double lowerBound = -1;
	size_t n_passFilter_i = 0;
	LogicalFilterMode logicalMode = LogicalFilterMode::Mat;
	double nan = -1;
	double eps = -1;
};

struct OutputArgs
{
	double * iou = nullptr;
};

typedef boost::icl::interval_set<double> IntervalSetType;

void ParseInput(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], InputArgs& _inArgs, OutputArgs& _outArgs);
void CheckInputValidity(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void HandleArgs(const InputArgs& _inArgs, OutputArgs& _outArgs);
double IntervalSetSize(IntervalSetType&& _set);
void Double2Interval(const double * _interval, IntervalSetType& _set);
double CalcIOU(const double * _interval1, const double * _interval2, const double _nan);
double CalcIOU(const IntervalSetType& _set1, const IntervalSetType& _set2, const double _nan);
void filterIndMode(const InputArgs& _inArgs, OutputArgs& _outArgs);
void LogicalMatMode(const InputArgs& _inArgs, OutputArgs& _outArgs);
void LogicalMatNoFilterMode(const InputArgs& _inArgs, OutputArgs& _outArgs);
void LogicalVecMode(const InputArgs& _inArgs, OutputArgs& _outArgs);
void LogicalLBMode(const InputArgs& _inArgs, OutputArgs& _outArgs);
inline bool IsPassLowerBound(const double * _interval1, const double * _interval2, double _lb, double _eps);
inline bool IsIntersects(const double * _a, const double * _b);
inline double IntersectionSize(const double * _a, const double * _b);
inline double UnionSize(const double * _a, const double * _b);

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
	_inArgs.set1 = (double*)mxGetData(prhs[argi]);
	_inArgs.n_intervals1 = (size_t)mxGetN(prhs[argi]);
	size_t set1_r = (size_t)mxGetM(prhs[argi++]);
	if (set1_r != 2)
	{
		std::stringstream ss;
		ss << "[ParseInput] - first argument contain double pairs (2xN), instead was: " << set1_r << " x " << _inArgs.n_intervals1;
		throw std::runtime_error(ss.str().c_str());
	}

	_inArgs.set2 = (double*)mxGetData(prhs[argi]);
	_inArgs.n_intervals2 = (size_t)mxGetN(prhs[argi]);
	size_t set2_r = (size_t)mxGetM(prhs[argi++]);
	if (set2_r != 2)
	{
		std::stringstream ss;
		ss << "[ParseInput] - second argument contain double pairs (2xN), instead was: " << set2_r << " x " << _inArgs.n_intervals2;
		throw std::runtime_error(ss.str().c_str());
	}
	if (argi < nrhs)
	{
		size_t filter_r = (size_t)mxGetM(prhs[argi]);
		size_t filter_c = (size_t)mxGetN(prhs[argi]);
		size_t totalNIntervals = _inArgs.n_intervals1 * _inArgs.n_intervals2;
		if (mxGetClassID(prhs[argi]) == mxDOUBLE_CLASS)
		{
			if (filter_r == 1 && filter_c == 1)
			{
				_inArgs.lowerBound = (double)mxGetScalar(prhs[argi]);
				_inArgs.logicalMode = InputArgs::LogicalFilterMode::LB;
			}
			else
			{
				_inArgs.passFilter_i = (double*)mxGetData(prhs[argi]);
				_inArgs.n_passFilter_i = filter_r == 1 ? filter_c : filter_r;
			}			
		}
		else if (mxGetClassID(prhs[argi]) == mxLOGICAL_CLASS)
		{
			_inArgs.passFilter = (bool*)mxGetData(prhs[argi]);
			if ((filter_r == 1 && filter_c == totalNIntervals ) || (filter_c == 1 && filter_r == totalNIntervals))
			{
				_inArgs.logicalMode = InputArgs::LogicalFilterMode::Vec;
			}
			else if (filter_r == _inArgs.n_intervals1 && filter_c == _inArgs.n_intervals2)
			{
				_inArgs.logicalMode = InputArgs::LogicalFilterMode::Mat;
			}
			else
			{
				std::stringstream ss;
				ss << "[ParseInput] - filter argument expected size: " << _inArgs.n_intervals1 << " x " << _inArgs.n_intervals2 <<", or "<< totalNIntervals<<" x 1, or 1 x"<< totalNIntervals<<
					", instead was: " << filter_r << " x " << filter_c;
				throw std::runtime_error(ss.str().c_str());
			}
		}
		else
		{
			throw std::runtime_error("Unknown filter type. should be logical or double.");
		}
		++argi;
	}

	_inArgs.nan = mxGetNaN();
	_inArgs.eps = mxGetEps();
	
	if (nlhs == 1)
	{
		plhs[0] = mxCreateDoubleMatrix(_inArgs.n_intervals1, _inArgs.n_intervals2, mxREAL);
		_outArgs.iou = (double *)mxGetData(plhs[0]);
	}	
}

void CheckInputValidity(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int expectedNrhs = 2;
	int expectedNlhs = 1;
	if (nrhs < expectedNrhs)
	{
		std::stringstream ss;
		ss << "[ParseInput]: " << expectedNrhs << " input arguments required, instead " << nrhs << " arguments were recieved.";
		throw std::runtime_error(ss.str().c_str());
	}
	if (nlhs > expectedNlhs)
	{
		std::stringstream ss;
		ss << "[ParseInput]: " << expectedNlhs << " output arguments required, instead " << nlhs << " arguments were recieved.";
		throw std::runtime_error(ss.str().c_str());
	}

	std::vector<mxClassID> expectedInArgsClasses = { mxDOUBLE_CLASS,mxDOUBLE_CLASS};
	std::vector<size_t> expectedInArgsDims = { 2,2,2 };
	for (int i = 0; i < expectedNrhs; ++i)
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
	if (_inArgs.passFilter_i != nullptr)
	{
		filterIndMode(_inArgs, _outArgs);
	}
	else if (_inArgs.logicalMode == InputArgs::LogicalFilterMode::Mat)
	{
		if (_inArgs.passFilter == nullptr)
		{
			LogicalMatNoFilterMode(_inArgs, _outArgs);
		}
		else
		{
			LogicalMatMode(_inArgs, _outArgs);
		}
	}
	else if (_inArgs.logicalMode == InputArgs::LogicalFilterMode::LB)
	{
		LogicalLBMode(_inArgs, _outArgs);
	}
	else
	{
		LogicalVecMode(_inArgs, _outArgs);
	}
}
void filterIndMode(const InputArgs& _inArgs, OutputArgs& _outArgs)
{
	size_t n_pairs = static_cast<size_t>(_inArgs.n_intervals1*_inArgs.n_intervals2);
	Concurrency::parallel_for((size_t)0, _inArgs.n_passFilter_i, [&](size_t i)
	{
		size_t pi = static_cast<size_t>(_inArgs.passFilter_i[i]) - 1;
		if (pi >= 0 && pi < n_pairs)
		{
			size_t ai = pi % _inArgs.n_intervals1;
			size_t bi = pi / _inArgs.n_intervals1;
			if (ai >= 0 && ai < _inArgs.n_intervals1 && bi >= 0 && bi < _inArgs.n_intervals2)
			{
				double iou = CalcIOU(&_inArgs.set1[ai * 2], &_inArgs.set2[bi * 2], _inArgs.nan);
				if (_outArgs.iou != nullptr)
				{
					_outArgs.iou[pi] = iou;
				}
			}
		}
	});
}
void LogicalMatNoFilterMode(const InputArgs& _inArgs, OutputArgs& _outArgs)
{
	std::vector<IntervalSetType> intervals1(_inArgs.n_intervals1);
	Concurrency::parallel_for((size_t)0, _inArgs.n_intervals1, [&](size_t ai)
	{
		const double* a = &_inArgs.set1[ai * 2];
		Double2Interval(a, intervals1[ai]);
	});
	std::vector<IntervalSetType> intervals2(_inArgs.n_intervals2);
	Concurrency::parallel_for((size_t)0, _inArgs.n_intervals2, [&](size_t bi)
	{
		const double* b = &_inArgs.set2[bi * 2];
		Double2Interval(b, intervals2[bi]);
	});

	Concurrency::parallel_for((size_t)0, _inArgs.n_intervals1, [&](size_t ai)
		//for (size_t ai = 0; ai < _inArgs.n_intervals1; ++ai)
	{
		//Concurrency::parallel_for((size_t)0, _inArgs.n_intervals2, [&](size_t bi)
		for (size_t bi = 0; bi < _inArgs.n_intervals2; ++bi)
		{
			const double* interval1 = &_inArgs.set1[ai * 2];
			const double* interval2 = &_inArgs.set2[bi * 2];
			if (IsIntersects(interval1, interval2))
			{
				double iou = CalcIOU(intervals1[ai], intervals2[bi], _inArgs.nan);
				
				if (_outArgs.iou != nullptr)
				{
					size_t pairIdx = bi*_inArgs.n_intervals1 + ai;
					_outArgs.iou[pairIdx] = iou;
				}
			}
		}//);
	});
}
void LogicalMatMode(const InputArgs& _inArgs, OutputArgs& _outArgs)
{
	Concurrency::parallel_for((size_t)0, _inArgs.n_intervals1, [&](size_t ai)
	//for (size_t ai = 0; ai < _inArgs.n_intervals1; ++ai)
	{
		//Concurrency::parallel_for((size_t)0, _inArgs.n_intervals2, [&](size_t bi)
		for (size_t bi = 0; bi < _inArgs.n_intervals2; ++bi)
		{
			const double* interval1 = &_inArgs.set1[ai * 2];
			const double* interval2 = &_inArgs.set2[bi * 2];
			size_t pairIdx = bi*_inArgs.n_intervals1 + ai;
			if (_inArgs.passFilter == nullptr || _inArgs.passFilter[pairIdx])
			{
				IntervalSetType a;
				Double2Interval(interval1, a);
				IntervalSetType b;
				Double2Interval(interval2, b);
				double iou = CalcIOU(a, b, _inArgs.nan);
				if (_outArgs.iou != nullptr)
				{
					_outArgs.iou[pairIdx] = iou;
				}
			}
		}//);
	});
}
//void LogicalLBMode(const InputArgs& _inArgs, OutputArgs& _outArgs)
//{
//	std::vector<IntervalSetType> intervals1(_inArgs.n_intervals1);
//	Concurrency::parallel_for((size_t)0, _inArgs.n_intervals1, [&](size_t ai)
//	{
//		const double* a = &_inArgs.set1[ai * 2];
//		Double2Interval(a, intervals1[ai]);
//	});
//	std::vector<IntervalSetType> intervals2(_inArgs.n_intervals2);
//	Concurrency::parallel_for((size_t)0, _inArgs.n_intervals2, [&](size_t bi)
//	{
//		const double* b = &_inArgs.set2[bi * 2];
//		Double2Interval(b, intervals2[bi]);
//	});
//
//	Concurrency::parallel_for((size_t)0, _inArgs.n_intervals1, [&](size_t ai)
//	//for (size_t ai = 0; ai < _inArgs.n_intervals1; ++ai)
//	{
//		//Concurrency::parallel_for((size_t)0, _inArgs.n_intervals2, [&](size_t bi)
//		for (size_t bi = 0; bi < _inArgs.n_intervals2; ++bi)
//		{
//			const double* interval1 = &_inArgs.set1[ai * 2];
//			const double* interval2 = &_inArgs.set2[bi * 2];
//			/*bool isPassLowerBound = false;
//			double dif1 = interval1[1] - interval1[0];
//			if (dif1 < 0)
//			{
//				dif1 = 1 + _inArgs.eps + dif1;
//			}
//			double dif2 = interval2[1] - interval2[0];
//			if (dif2 < 0)
//			{
//				dif2 = 1 + _inArgs.eps + dif2;
//			}
//			if (dif1 < dif2)
//			{
//				isPassLowerBound = dif1 / dif2 >= _inArgs.lowerBound;
//			}
//			else
//			{
//				isPassLowerBound = dif2 / dif1 >= _inArgs.lowerBound;
//			}*/
//			if (IsPassLowerBound(interval1, interval2, _inArgs.lowerBound, _inArgs.eps) && IsIntersects(interval1, interval2))
//			{
//				double iou = CalcIOU(intervals1[ai], intervals2[bi], _inArgs.nan);
//				if (_outArgs.iou != nullptr)
//				{
//					size_t pairIdx = bi*_inArgs.n_intervals1 + ai;
//					_outArgs.iou[pairIdx] = iou;
//				}
//			}
//		}//);
//	});
//}
void LogicalLBMode(const InputArgs& _inArgs, OutputArgs& _outArgs)
{
	Concurrency::parallel_for((size_t)0, _inArgs.n_intervals1, [&](size_t ai)
		//for (size_t ai = 0; ai < _inArgs.n_intervals1; ++ai)
	{
		//Concurrency::parallel_for((size_t)0, _inArgs.n_intervals2, [&](size_t bi)
		for (size_t bi = 0; bi < _inArgs.n_intervals2; ++bi)
		{
			const double* interval1 = &_inArgs.set1[ai * 2];
			const double* interval2 = &_inArgs.set2[bi * 2];
			if (IsIntersects(interval1, interval2))
			{
				double iou = CalcIOU(interval1, interval2, _inArgs.nan);
				if (_outArgs.iou != nullptr)
				{
					size_t pairIdx = bi*_inArgs.n_intervals1 + ai;
					_outArgs.iou[pairIdx] = iou;
				}
			}
		}//);
	});
}
void LogicalVecMode(const InputArgs& _inArgs, OutputArgs& _outArgs)
{
	size_t n_pairs = _inArgs.n_intervals1*_inArgs.n_intervals2;
	Concurrency::parallel_for((size_t)0, n_pairs, [&](size_t pi) 
	{
		if (_inArgs.passFilter == nullptr || _inArgs.passFilter[pi])
		{
			size_t ai = pi%_inArgs.n_intervals1;
			size_t bi = pi/_inArgs.n_intervals1;
			if (ai >= 0 && ai < _inArgs.n_intervals1 && bi >= 0 && bi < _inArgs.n_intervals2)
			{
				double iou = CalcIOU(&_inArgs.set1[ai * 2], &_inArgs.set2[bi * 2], _inArgs.nan);
				if (_outArgs.iou != nullptr)
				{
					_outArgs.iou[pi] = iou;
				}
			}
		}		
	});
}
bool IsPassLowerBound(const double * _interval1, const double * _interval2, double _lb, double _eps)
{
	double dif1 = _interval1[1] - _interval1[0];
	if (dif1 < 0)
	{
		dif1 = 1 + _eps + dif1;
	}
	double dif2 = _interval2[1] - _interval2[0];
	if (dif2 < 0)
	{
		dif2 = 1 + _eps + dif2;
	}
	if (dif1 < dif2)
	{
		return dif1 / dif2 >= _lb;
	}
	return dif2 / dif1 >= _lb;
}

double CalcIOU(const double * _interval1, const double * _interval2, const double _nan)
{
	double intersectionSum = IntersectionSize(_interval1, _interval2);
	double iou = 0;
	if (intersectionSum != 0)
	{
		double unionSum = UnionSize(_interval1, _interval2);
		if (unionSum != 0)
		{
			iou = intersectionSum / unionSum;
		}
		else
		{
			iou = _nan;
		}
	}
	return iou;
}
double CalcIOU(const IntervalSetType& _set1, const IntervalSetType& _set2, const double _nan)
{
	double intersectionSum = IntervalSetSize(_set1 & _set2);
	double iou = 0;
	if (intersectionSum != 0)
	{
		double unionSum = IntervalSetSize(_set1 + _set2);
		if (unionSum != 0)
		{
			iou = intersectionSum / unionSum;
		}
		else
		{
			iou = _nan;
		}
	}
	return iou;
}
double IntervalSetSize(IntervalSetType&& _set)
{
	double sum = 0.0;
	for (auto& interval : _set)
	{
		sum += interval.upper() - interval.lower();
	}
	return sum;
}
void Double2Interval(const double * _interval, IntervalSetType& _set)
{
	if (_interval[0] > _interval[1])
	{
		_set.add(IntervalSetType::interval_type::closed(_interval[0], 1.0));
		_set.add(IntervalSetType::interval_type::closed(0.0, _interval[1]));
	}
	else
	{
		_set.add(IntervalSetType::interval_type::closed(_interval[0], _interval[1]));
	}
}

bool IsIntersects(const double * _a, const double * _b)
{
	double x1 = _a[0];
	double x2 = _a[1];
	double y1 = _b[0];
	double y2 = _b[1];
	//_a==[x1,x2], _b==[y1,y2]
	if (x1 <= x2)
	{
		if (y1 <= y2) // x1<=x2 && y1<=y2
		{
			return x1 <= y2 && y1 <= x2;
		}
		else // x1<=x2 && y1>y2
		{
			return y1 <= x2 || x1 <= y2;
		}
	}
	else
	{
		if (y1 <= y2) // x1>x2 && y1<=y2
		{
			return x1 <= y2 || y1 <= x2;
		}
		else // x1>x2 && y1>y2
		{
			return true;
		}
	}
}
double IntersectionSize(const double * _a, const double * _b)
{
	double x1 = _a[0];
	double x2 = _a[1];
	double y1 = _b[0];
	double y2 = _b[1];

	//_a==[x1,x2], _b==[y1,y2]
	if (x1 <= x2)
	{
		if (y1 <= y2) // x1<=x2 && y1<=y2
		{
			return std::min(x2, y2) - std::max(x1, y1);
		}
		else // x1<=x2 && y1>y2
		{
			return std::max(x2 - std::max(x1, y1), 0.0) + std::max(std::min(x2, y2) - x1, 0.0);
		}
	}
	else
	{
		if (y1 <= y2) // x1>x2 && y1<=y2
		{
			return std::max(y2 - std::max(x1, y1), 0.0) + std::max(std::min(x2, y2) - y1, 0.0);
		}
		else // x1>x2 && y1>y2
		{
			return std::max(1 - std::max(x1, y1), 0.0) + std::max(y2 - x1, 0.0) + std::max(x2 - y1, 0.0) + std::max(std::min(x2, y2), 0.0);
		}
	}
}
double UnionSize(const double * _a, const double * _b)
{
	double x1 = _a[0];
	double x2 = _a[1];
	double y1 = _b[0];
	double y2 = _b[1];

	//_a==[x1,x2], _b==[y1,y2]
	if (x1 <= x2)
	{
		if (y1 <= y2) // x1<=x2 && y1<=y2
		{
			return std::max(x2, y2) - std::min(x1, y1);
		}
		else // x1<=x2 && y1>y2
		{
			if (x1 <= y2)
			{
				if (y1 <= x2)
				{
					return 1;
				}
				else
				{
					return std::max(x2, y2) + 1 - y1;
				}
			}
			else
			{
				if (y1 <= x2)
				{
					return y2 + 1 - std::min(x1, y1);
				}
				else
				{
					return y2 + 1 - y1 + x2 - x1;
				}
			}
		}
	}
	else
	{
		if (y1 <= y2) // x1>x2 && y1<=y2
		{
			if (y1 <= x2)
			{
				if (x1 <= y2)
				{
					return 1;
				}
				else
				{
					return std::max(x2, y2) + 1 - x1;
				}
			}
			else
			{
				if (x1 <= y2)
				{
					return x2 + 1 - std::min(x1, y1);
				}
				else
				{
					return x2 + 1 - x1 + y2 - y1;
				}
			}
		}
		else // x1>x2 && y1>y2
		{
			if (y1 <= x2 || y2 >= x1)
			{
				return 1;
			}
			else
			{
				return 1 - std::min(x1, y1) + std::max(y2, x2);
			}
		}
	}
}