#define main FullMovingRegressionMain
#include "test_moving_immersed_transient_flow.cpp"
#undef main

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);int status=0;
	try {
		RunRigidTranslation(iga::CutCellVolumeQuadratureStorageMode::Expanded);
		RunRigidTranslation(iga::CutCellVolumeQuadratureStorageMode::Compact);
		std::cout<<"moving_immersed_rigid_translation_test: PASS expanded/compact original trace and conservation gates"<<std::endl;
	} catch(const std::exception& error) { std::cerr<<error.what()<<std::endl;status=1; }
	PetscFinalize();return status;
}
