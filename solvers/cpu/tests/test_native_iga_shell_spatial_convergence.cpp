#include "NativeIgaShellSurface.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

std::vector<double> OpenUniformKnots(int elements,int degree)
{
	std::vector<double> knots(static_cast<std::size_t>(degree+1),0.0);
	for(int i=1;i<elements;++i) knots.push_back(static_cast<double>(i)/elements);
	knots.insert(knots.end(),static_cast<std::size_t>(degree+1),1.0);
	return knots;
}

std::vector<double> Solve(std::vector<double> matrix,std::vector<double> right)
{
	const std::size_t n=right.size();
	for(std::size_t column=0;column<n;++column) {
		std::size_t pivot=column;
		for(std::size_t row=column+1;row<n;++row)
			if(std::abs(matrix[row*n+column])>std::abs(matrix[pivot*n+column])) pivot=row;
		assert(std::abs(matrix[pivot*n+column])>1e-14);
		for(std::size_t j=0;j<n;++j) std::swap(matrix[column*n+j],matrix[pivot*n+j]);
		std::swap(right[column],right[pivot]);
		for(std::size_t row=column+1;row<n;++row) {
			const double factor=matrix[row*n+column]/matrix[column*n+column];
			for(std::size_t j=column;j<n;++j) matrix[row*n+j]-=factor*matrix[column*n+j];
			right[row]-=factor*right[column];
		}
	}
	std::vector<double> answer(n,0.0);
	for(std::size_t reverse=0;reverse<n;++reverse) {
		const std::size_t row=n-1-reverse;double value=right[row];
		for(std::size_t j=row+1;j<n;++j) value-=matrix[row*n+j]*answer[j];
		answer[row]=value/matrix[row*n+row];
	}
	return answer;
}

struct Observation { double error=0.0;int controls=0; };

Observation Measure(int elements)
{
	const int degree=3;
	const auto knots=OpenUniformKnots(elements,degree);
	const int controls=static_cast<int>(knots.size())-degree-1;
	std::vector<double> greville(static_cast<std::size_t>(controls));
	for(int i=0;i<controls;++i)
		for(int j=1;j<=degree;++j) greville[static_cast<std::size_t>(i)]+=knots[static_cast<std::size_t>(i+j)]/degree;
	std::vector<double> collocation(static_cast<std::size_t>(controls*controls)),values(static_cast<std::size_t>(controls));
	const double pi=std::acos(-1.0);
	for(int row=0;row<controls;++row) {
		const auto basis=iga::EvaluateNativeSplineBasis(knots,degree,greville[static_cast<std::size_t>(row)]);
		for(int column=0;column<controls;++column)
			collocation[static_cast<std::size_t>(row*controls+column)]=basis.value[static_cast<std::size_t>(column)];
		values[static_cast<std::size_t>(row)]=std::sin(pi*greville[static_cast<std::size_t>(row)]);
	}
	const auto coefficients=Solve(std::move(collocation),std::move(values));
	const std::array<double,4> x{{-0.8611363115940526,-0.3399810435848563,
		0.3399810435848563,0.8611363115940526}};
	const std::array<double,4> w{{0.3478548451374538,0.6521451548625461,
		0.6521451548625461,0.3478548451374538}};
	double error=0.0,reference=0.0;
	for(int ev=0;ev<elements;++ev) for(int eu=0;eu<elements;++eu)
		for(std::size_t gv=0;gv<4;++gv) for(std::size_t gu=0;gu<4;++gu) {
			const double u=(eu+0.5*(1.0+x[gu]))/elements;
			const double v=(ev+0.5*(1.0+x[gv]))/elements;
			const auto bu=iga::EvaluateNativeSplineBasis(knots,degree,u);
			const auto bv=iga::EvaluateNativeSplineBasis(knots,degree,v);
			double uu=0.0,uv=0.0,vv=0.0;
			for(int j=0;j<controls;++j) for(int i=0;i<controls;++i) {
				const double coefficient=coefficients[static_cast<std::size_t>(i)]
					*coefficients[static_cast<std::size_t>(j)];
				uu+=coefficient*bu.second[static_cast<std::size_t>(i)]*bv.value[static_cast<std::size_t>(j)];
				uv+=coefficient*bu.first[static_cast<std::size_t>(i)]*bv.first[static_cast<std::size_t>(j)];
				vv+=coefficient*bu.value[static_cast<std::size_t>(i)]*bv.second[static_cast<std::size_t>(j)];
			}
			const double exact_uu=-pi*pi*std::sin(pi*u)*std::sin(pi*v);
			const double exact_uv=pi*pi*std::cos(pi*u)*std::cos(pi*v);
			const double exact_vv=exact_uu;
			const double weight=w[gu]*w[gv]/(4.0*elements*elements);
			error+=weight*((uu-exact_uu)*(uu-exact_uu)+2.0*(uv-exact_uv)*(uv-exact_uv)
				+(vv-exact_vv)*(vv-exact_vv));
			reference+=weight*(exact_uu*exact_uu+2.0*exact_uv*exact_uv+exact_vv*exact_vv);
		}
	return {std::sqrt(error/reference),controls};
}

}

int main()
{
	const auto coarse=Measure(2),medium=Measure(4),fine=Measure(8);
	assert(coarse.error>medium.error&&medium.error>fine.error);
	const double order_1=std::log(coarse.error/medium.error)/std::log(2.0);
	const double order_2=std::log(medium.error/fine.error)/std::log(2.0);
	assert(order_1>1.5&&order_2>1.5);
	std::cout<<"native IGA shell curvature convergence passed errors="<<coarse.error<<','
		<<medium.error<<','<<fine.error<<" orders="<<order_1<<','<<order_2<<'\n';
	return 0;
}
