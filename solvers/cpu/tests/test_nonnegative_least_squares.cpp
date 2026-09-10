#include "NonnegativeLeastSquares.hpp"
#include <iostream>
#include <random>

namespace {
void Check(bool value,const char* message)
{
	if(!value)throw std::runtime_error(message);
}
template<class Function> void Reject(Function function)
{
	bool rejected=false;try{function();}catch(const std::exception&){rejected=true;}
	Check(rejected,"invalid NNLS operation was accepted");
}
void Kkt(const std::vector<double>& a,const std::vector<double>& b,const iga::NonnegativeLeastSquaresResult& fit)
{
	std::vector<long double> residual(b.begin(),b.end());
	for(std::size_t j=0;j<fit.coefficients.size();++j) {
		Check(fit.coefficients[j]>=0,"negative NNLS coefficient");
		for(std::size_t i=0;i<b.size();++i)residual[i]-=static_cast<long double>(a[j*b.size()+i])*fit.coefficients[j];
	}
	for(std::size_t j=0;j<fit.coefficients.size();++j) {
		long double dual=0;for(std::size_t i=0;i<b.size();++i)dual+=a[j*b.size()+i]*residual[i];
		Check(dual<2e-11L,"NNLS inactive KKT condition failed");
		if(fit.coefficients[j]>0)Check(std::abs(dual)<2e-11L,"NNLS passive KKT condition failed");
	}
}
}

int main()
{
	try {
		const std::vector<double> identity{1,0,0,1};
		const auto exact=iga::FitNonnegativeLeastSquares(identity,{2,3});
		Check(exact.within_tolerance&&std::abs(exact.coefficients[0]-2)<1e-14&&std::abs(exact.coefficients[1]-3)<1e-14,"identity fit failed");
		const auto constrained=iga::FitNonnegativeLeastSquares(identity,{-2,3});
		Check(!constrained.within_tolerance&&constrained.coefficients[0]==0&&std::abs(constrained.coefficients[1]-3)<1e-14&&std::abs(constrained.residual_norm-2)<1e-14L,"infeasible fit was misreported");
		Kkt(identity,{-2,3},constrained);
		const auto duplicate=iga::FitNonnegativeLeastSquares({1,0,1,0,0,0},{2,1});
		Check(!duplicate.within_tolerance&&std::abs(duplicate.coefficients[0]+duplicate.coefficients[1]-2)<1e-14,"rank-deficient fit failed");
		const auto zero=iga::FitNonnegativeLeastSquares(identity,{0,0});
		Check(zero.within_tolerance&&zero.positive_columns==0&&zero.iterations==0,"zero RHS failed");
		for(double scale:{1e-100,1.,1e100}) {
			const auto scaled=iga::FitNonnegativeLeastSquares({1e-100,0,0,1e100},{scale,2*scale});
			Check(scaled.within_tolerance&&scaled.relative_residual<1e-13L,"scaled NNLS failed");
		}
		// The exact positive solution is (1/2,1/2). After selecting the first
		// column the second reduced gradient is delta^2/2, while the residual
		// remains delta/2: an absolute epsilon floor must not claim stationarity.
		for(double delta:{1e-8,1e-9,1e-10}) {
			const std::vector<double> a{1,0,1,delta},b{1,delta/2};
			const auto fit=iga::FitNonnegativeLeastSquares(a,b);
			Check(fit.within_tolerance&&fit.relative_residual<=1e-13L,"small resolvable NNLS gradient stopped prematurely");
			Check(std::abs(fit.coefficients[0]-.5)<1e-6&&std::abs(fit.coefficients[1]-.5)<1e-6,"ill-conditioned positive solution differs");
			Kkt(a,b,fit);
		}
		std::mt19937 generator(49217);std::uniform_real_distribution<double> uniform(-1,1);
		std::size_t maximum_iterations=0,removals=0;
		for(unsigned trial=0;trial<200;++trial) {
			std::vector<double> a(3*7),b(3);for(auto& value:a)value=uniform(generator);for(auto& value:b)value=uniform(generator);
			const auto fit=iga::FitNonnegativeLeastSquares(a,b);Kkt(a,b,fit);
			maximum_iterations=std::max(maximum_iterations,fit.iterations);
			removals+=fit.removals;
		}
		Check(removals>0,"active-set removal was not exercised");
		Reject([&]{iga::FitNonnegativeLeastSquares({},{});});
		Reject([&]{iga::FitNonnegativeLeastSquares({1,2,3},{1,2});});
		Reject([&]{iga::FitNonnegativeLeastSquares({std::numeric_limits<double>::quiet_NaN()},{1});});
		Reject([&]{iga::FitNonnegativeLeastSquares({1},{std::numeric_limits<double>::infinity()});});
		Reject([&]{iga::NonnegativeLeastSquaresOptions o;o.max_workspace_bytes=1;iga::FitNonnegativeLeastSquares(identity,{1,1},o);});
		Reject([&]{iga::NonnegativeLeastSquaresOptions o;o.max_rows=1;iga::FitNonnegativeLeastSquares(identity,{1,1},o);});
		Reject([&]{iga::NonnegativeLeastSquaresOptions o;o.max_columns=1;iga::FitNonnegativeLeastSquares(identity,{1,1},o);});
		Reject([&]{iga::NonnegativeLeastSquaresOptions o;o.max_iterations=1;iga::FitNonnegativeLeastSquares(identity,{1,1},o);});
		Reject([&]{iga::FitNonnegativeLeastSquares({std::numeric_limits<double>::denorm_min()},{std::numeric_limits<double>::max()});});
		std::cout<<"nonnegative_least_squares_test: PASS analytic, scaling, 200 KKT, rank and cap gates; maximum_iterations="<<maximum_iterations<<'\n';
	} catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
	return 0;
}
