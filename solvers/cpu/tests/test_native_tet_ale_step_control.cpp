#include "NativeTetAleStepControl.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

int main()
{
	const std::vector<std::array<double,3>> points{{{0,0,0},{1,0,0},{0,1,0},{0,0,1}}};
	iga::NativeTetCell cell;cell.nodes={{0,1,2,3}};
	iga::NativeTetAleKinematics kinematics(points,{cell});
	const auto motion=[](double time) {
		std::vector<std::array<double,3>> displacement(4,{{0,0,0}});
		displacement[1][0]=-12.0*time;
		return displacement;
	};
	iga::NativeTetAleRetryPolicy policy;
	policy.maximum_attempts=4;policy.reduction_factor=0.5;policy.minimum_dt_s=0.01;
	const auto retried=iga::BeginNativeTetAleQualityControlledTrial(
		kinematics,0.2,motion,policy);
	assert(retried.attempts==3&&retried.rejected_attempts==2);
	assert(std::abs(retried.accepted_dt_s-0.05)<1e-15);
	assert(kinematics.HasTrial()&&kinematics.CommittedTime()==0.0);
	kinematics.CommitTrial();
	assert(std::abs(kinematics.CommittedTime()-0.05)<1e-15);

	iga::NativeTetAleKinematics failing(points,{cell});
	policy.maximum_attempts=3;
	bool rejected=false;
	try {
		(void)iga::BeginNativeTetAleQualityControlledTrial(failing,0.2,
			[](double) {
				std::vector<std::array<double,3>> displacement(4,{{0,0,0}});
				displacement[1][0]=-2.0;return displacement;
			},policy);
	} catch (const std::runtime_error&) { rejected=true; }
	assert(rejected&&!failing.HasTrial()&&failing.CommittedTime()==0.0);
	std::cout<<"native tetrahedral ALE quality retry tests passed\n";
	return 0;
}
