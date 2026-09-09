#include "ElementAssemblyExecution.hpp"
#include <iostream>

namespace {
void Check(bool value, const char* message)
{
	if(!value) throw std::runtime_error(message);
}
template<class Work> void Reject(Work work, const char* diagnostic)
{
	bool rejected=false;
	try { work(); }
	catch(const std::exception& error) { rejected=std::string(error.what()).find(diagnostic)!=std::string::npos; }
	Check(rejected,"expected execution rejection was absent or incorrect");
}
}

int main(int argc, char** argv)
{
	const bool single=argc>1 && std::string(argv[1])=="--mpi-single";
	int provided=0;
	MPI_Init_thread(&argc,&argv,single ? MPI_THREAD_SINGLE : MPI_THREAD_FUNNELED,&provided);
	int status=0;
	try {
		::setenv("IGA_ASSEMBLY_THREADS","1",1);
		::unsetenv("IGA_ASSEMBLY_BATCH_SIZE");
		iga::ElementAssemblyExecution serial;
		Check(serial.Options().threads==1 && serial.Options().capacity==1,"serial memory default differs");
		serial.RequireCaller();
		for(const char* name:{"IGA_ASSEMBLY_THREADS","IGA_ASSEMBLY_BATCH_SIZE"}) {
			for(const char* value:{"0","-1","","2x","1,2","2147483648"}) {
				::setenv(name,value,1);
				Reject([] { iga::ElementAssemblyExecution invalid; },name);
			}
			::setenv(name,"1",1);
		}
		::setenv("IGA_ASSEMBLY_THREADS","2",1);
#ifdef _OPENMP
		if(provided<MPI_THREAD_FUNNELED) {
			Reject([] { iga::ElementAssemblyExecution invalid; },"MPI_THREAD_FUNNELED");
		} else {
			::setenv("IGA_ASSEMBLY_BATCH_SIZE","3",1);
			iga::ElementAssemblyExecution parallel;
			Check(parallel.Options().threads==2 && parallel.Options().capacity==3,"explicit resources differ");
			// Environment changes do not mutate a live runtime's policy.
			::setenv("IGA_ASSEMBLY_THREADS","1",1);
			Check(parallel.Options().threads==2,"execution settings were not frozen");
			bool rejected=false;
			std::thread other([&] { try { parallel.RequireCaller(); } catch(const std::logic_error&) { rejected=true; } });
			other.join(); Check(rejected,"foreign caller accepted");
		}
#else
		Reject([] { iga::ElementAssemblyExecution invalid; },"OpenMP build");
#endif
		std::cout << "element_assembly_execution mpi_thread=" << provided << " passed\n";
	} catch(const std::exception& error) { std::cerr << error.what() << '\n'; status=1; }
	MPI_Finalize();
	return status;
}
