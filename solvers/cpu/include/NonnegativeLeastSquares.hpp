#ifndef IGA_NONNEGATIVE_LEAST_SQUARES_HPP
#define IGA_NONNEGATIVE_LEAST_SQUARES_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace iga {

struct NonnegativeLeastSquaresOptions {
	std::size_t max_rows=512,max_columns=8192,max_iterations=10000;
	std::size_t max_workspace_bytes=128u*1024u*1024u;
	double relative_tolerance=1e-13,absolute_tolerance=0;
};

struct NonnegativeLeastSquaresResult {
	std::vector<double> coefficients;
	long double residual_norm=0,relative_residual=0;
	std::size_t iterations=0,positive_columns=0,workspace_bound_bytes=0;
	std::size_t removals=0,rank_rejections=0;
	bool within_tolerance=false;
};

namespace nonnegative_least_squares_detail {
inline std::size_t Product(std::size_t a,std::size_t b)
{
	if(b&&a>std::numeric_limits<std::size_t>::max()/b) throw std::overflow_error("NNLS workspace size overflows");
	return a*b;
}
inline std::size_t Add(std::size_t a,std::size_t b)
{
	if(a>std::numeric_limits<std::size_t>::max()-b) throw std::overflow_error("NNLS workspace size overflows");
	return a+b;
}
inline long double Norm(const std::vector<long double>& vector)
{
	long double value=0;for(auto x:vector)value=std::hypot(value,x);return value;
}
}

// Column-major matrix, rows equal rhs.size(). Active-set NNLS with twice
// orthogonalized QR column insertion and Givens deletion. No normal equations
// or external linear-algebra dependency. An infeasible exact fit returns its
// residual with within_tolerance=false; callers must enforce their own gate.
inline NonnegativeLeastSquaresResult FitNonnegativeLeastSquares(
	const std::vector<double>& matrix,const std::vector<double>& rhs,
	NonnegativeLeastSquaresOptions options={})
{
	using namespace nonnegative_least_squares_detail;
	const std::size_t rows=rhs.size();
	if(!rows||matrix.empty()||matrix.size()%rows||!options.max_rows||!options.max_columns
		||!options.max_iterations||!options.max_workspace_bytes
		||!std::isfinite(options.relative_tolerance)||!(options.relative_tolerance>0)
		||!std::isfinite(options.absolute_tolerance)||options.absolute_tolerance<0)
		throw std::invalid_argument("invalid NNLS dimensions or options");
	const std::size_t columns=matrix.size()/rows;
	if(rows>options.max_rows||columns>options.max_columns) throw std::runtime_error("NNLS dimension cap reached");
	const auto square=Product(rows,rows);
	// Conservative bound includes result storage, temporaries and index/flag
	// vectors, but excludes the caller-owned matrix and right-hand side.
	const auto scalars=Add(Add(matrix.size(),Product(2,square)),Add(Product(8,columns),Product(8,rows)));
	const auto bytes=Add(Product(scalars,sizeof(long double)),Product(Product(8,columns),sizeof(std::size_t)));
	if(bytes>options.max_workspace_bytes) throw std::runtime_error("NNLS workspace cap reached");
	for(auto x:matrix)if(!std::isfinite(x))throw std::invalid_argument("NNLS matrix is not finite");
	for(auto x:rhs)if(!std::isfinite(x))throw std::invalid_argument("NNLS right-hand side is not finite");
	NonnegativeLeastSquaresResult result;result.workspace_bound_bytes=bytes;
	result.coefficients.assign(columns,0);
	std::vector<long double> b(rhs.begin(),rhs.end());
	const long double rhs_norm=Norm(b);
	if(!std::isfinite(rhs_norm))throw std::overflow_error("NNLS right-hand side norm is not finite");
	if(rhs_norm==0) { result.within_tolerance=true;return result; }
	for(auto& x:b)x/=rhs_norm;
	std::vector<long double> a(matrix.size()),norms(columns),x(columns),residual=b;
	for(std::size_t j=0;j<columns;++j) {
		long double norm=0;for(std::size_t i=0;i<rows;++i)norm=std::hypot(norm,static_cast<long double>(matrix[j*rows+i]));
		if(!std::isfinite(norm))throw std::overflow_error("NNLS column norm is not finite");
		norms[j]=norm;
		if(norm>0)for(std::size_t i=0;i<rows;++i)a[j*rows+i]=matrix[j*rows+i]/norm;
	}
	std::vector<long double> q(square),r(square),work(rows),z(rows);
	std::vector<std::size_t> passive;passive.reserve(rows);
	std::vector<bool> active(columns,false),blocked(columns,false);
	const long double epsilon=std::numeric_limits<long double>::epsilon();
	const auto tick=[&] {
		if(result.iterations==options.max_iterations)throw std::runtime_error("NNLS iteration cap reached");
		++result.iterations;
	};
	const auto remove=[&](std::size_t index) {
		++result.removals;
		const auto size=passive.size();active[passive[index]]=false;x[passive[index]]=0;
		for(std::size_t j=index;j+1<size;++j)for(std::size_t i=0;i<size;++i)r[i*rows+j]=r[i*rows+j+1];
		for(std::size_t j=index;j+1<size;++j) {
			const long double diagonal=r[j*rows+j],below=r[(j+1)*rows+j],length=std::hypot(diagonal,below);
			if(!(length>0))throw std::runtime_error("NNLS QR deletion lost rank");
			const long double cosine=diagonal/length,sine=below/length;
			for(std::size_t k=j;k+1<size;++k) {
				const auto top=r[j*rows+k],bottom=r[(j+1)*rows+k];
				r[j*rows+k]=cosine*top+sine*bottom;r[(j+1)*rows+k]=-sine*top+cosine*bottom;
			}
			for(std::size_t i=0;i<rows;++i) {
				const auto left=q[j*rows+i],right=q[(j+1)*rows+i];
				q[j*rows+i]=cosine*left+sine*right;q[(j+1)*rows+i]=-sine*left+cosine*right;
			}
		}
		passive.erase(passive.begin()+index);
	};
	while(true) {
		residual=b;
		for(auto j:passive)for(std::size_t i=0;i<rows;++i)residual[i]-=a[j*rows+i]*x[j];
		const auto residual_norm=Norm(residual);
		if(residual_norm<=options.relative_tolerance+options.absolute_tolerance/rhs_norm)break;
		// Columns have unit norm. Scale the reduced-gradient roundoff floor by
		// the current residual, otherwise a small but resolvable gradient can
		// stop a feasible ill-conditioned fit well before its residual gate.
		std::size_t selected=columns;long double maximum=32*epsilon*residual_norm;
		for(std::size_t j=0;j<columns;++j)if(!active[j]&&!blocked[j]&&norms[j]>0) {
			long double dual=0;for(std::size_t i=0;i<rows;++i)dual+=a[j*rows+i]*residual[i];
			if(dual>maximum) { maximum=dual;selected=j; }
		}
		if(selected==columns||passive.size()==rows)break;
		tick();const auto size=passive.size();
		for(std::size_t i=0;i<rows;++i)work[i]=a[selected*rows+i];
		for(std::size_t j=0;j<rows;++j)r[j*rows+size]=0;
		for(unsigned pass=0;pass<2;++pass)for(std::size_t j=0;j<size;++j) {
			long double dot=0;for(std::size_t i=0;i<rows;++i)dot+=q[j*rows+i]*work[i];
			r[j*rows+size]+=dot;for(std::size_t i=0;i<rows;++i)work[i]-=dot*q[j*rows+i];
		}
		const auto length=Norm(work);
		if(length<=128*epsilon*rows) { blocked[selected]=true;++result.rank_rejections;continue; }
		r[size*rows+size]=length;for(std::size_t i=0;i<rows;++i)q[size*rows+i]=work[i]/length;
		passive.push_back(selected);active[selected]=true;
		std::fill(blocked.begin(),blocked.end(),false);
		while(true) {
			for(std::size_t j=0;j<passive.size();++j) { z[j]=0;for(std::size_t i=0;i<rows;++i)z[j]+=q[j*rows+i]*b[i]; }
			for(std::size_t j=passive.size();j-- >0;) {
				for(std::size_t k=j+1;k<passive.size();++k)z[j]-=r[j*rows+k]*z[k];
				z[j]/=r[j*rows+j];if(!std::isfinite(z[j]))throw std::runtime_error("NNLS QR solution is not finite");
			}
			std::size_t hit=passive.size();long double alpha=1;
			for(std::size_t j=0;j<passive.size();++j)if(z[j]<=0) {
				const auto current=x[passive[j]];
				const long double step=current>0?current/(current-z[j]):0;
				if(hit==passive.size()||step<alpha) { alpha=step;hit=j; }
			}
			if(hit==passive.size()) { for(std::size_t j=0;j<passive.size();++j)x[passive[j]]=z[j];break; }
			for(std::size_t j=0;j<passive.size();++j)x[passive[j]]+=alpha*(z[j]-x[passive[j]]);
			x[passive[hit]]=0;
			for(std::size_t j=passive.size();j-- >0;)if(x[passive[j]]<=0) { tick();remove(j); }
			if(passive.empty())break;
		}
	}
	for(std::size_t j=0;j<columns;++j)if(x[j]>0) {
		const long double coefficient=x[j]*rhs_norm/norms[j];
		if(!std::isfinite(coefficient)||coefficient>std::numeric_limits<double>::max())throw std::overflow_error("NNLS coefficient is not representable");
		result.coefficients[j]=static_cast<double>(coefficient);
		if(result.coefficients[j]>0)++result.positive_columns;
	}
	// Audit the returned doubles against the original, unscaled problem.
	residual.assign(rhs.begin(),rhs.end());
	for(std::size_t j=0;j<columns;++j)if(result.coefficients[j]>0)
		for(std::size_t i=0;i<rows;++i)residual[i]-=static_cast<long double>(matrix[j*rows+i])*result.coefficients[j];
	result.residual_norm=Norm(residual);result.relative_residual=result.residual_norm/rhs_norm;
	if(!std::isfinite(result.residual_norm)||!std::isfinite(result.relative_residual))throw std::overflow_error("NNLS residual is not finite");
	result.within_tolerance=result.residual_norm<=options.absolute_tolerance+options.relative_tolerance*rhs_norm;
	return result;
}

}

#endif
