#include "PhaseProfile.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
double now = 0.0;
double Now() { return now; }
void Require(bool condition)
{
	if (!condition) throw std::runtime_error("phase accounting failed");
}
}

int main()
{
	using iga::ProfilePhase;
	iga::PhaseProfile profile(true, Now);
	{
		iga::PhaseScope outer(profile, ProfilePhase::Assembly);
		now = 2.0;
		{
			iga::PhaseScope child(profile, ProfilePhase::Communication);
			now = 5.0;
		}
		now = 7.0;
		outer.Stop();
		outer.Stop();
	}
	// Reentrant categories still partition exclusive time correctly.
	{
		iga::PhaseScope outer(profile, ProfilePhase::Assembly);
		{
			iga::PhaseScope child(profile, ProfilePhase::Assembly);
			now = 8.0;
		}
		now = 9.0;
	}
	try {
		iga::PhaseScope failing(profile, ProfilePhase::Output);
		now = 10.0;
		throw std::runtime_error("injected");
	} catch (const std::runtime_error&) {}
	Require(profile.Get(ProfilePhase::Assembly).inclusive == 10.0);
	Require(profile.Get(ProfilePhase::Assembly).exclusive == 6.0);
	Require(profile.Get(ProfilePhase::Assembly).calls == 3);
	Require(profile.Get(ProfilePhase::Communication).exclusive == 3.0);
	Require(profile.Get(ProfilePhase::Output).exclusive == 1.0);
	now = 12.0;
	std::ostringstream output;
	profile.Write(output, 1, 2, 1);
	Require(output.str().find("\"unscoped_s\":2") != std::string::npos);
	Require(output.str().find("\"rank\":1,\"ranks\":2,\"status\":1") != std::string::npos);
	iga::PhaseProfile disabled(false, Now);
	{
		iga::PhaseScope scope(disabled, ProfilePhase::Assembly);
		now = 15.0;
	}
	std::ostringstream empty;
	disabled.Write(empty, 0, 1, 0);
	Require(empty.str().empty() && disabled.Get(ProfilePhase::Assembly).calls == 0);
	std::cout << "phase_profile_test passed\n";
}
