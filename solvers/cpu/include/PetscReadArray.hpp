#ifndef IGA_PETSC_READ_ARRAY_HPP
#define IGA_PETSC_READ_ARRAY_HPP

#include <petscvec.h>
#include <stdexcept>

namespace iga {

// Local borrowed view. Restore during exception unwinding before its Vec is
// destroyed; a zero-length view can have a null data pointer and still be held.
class PetscReadArray {
public:
	PetscReadArray() = default;
	PetscReadArray(const PetscReadArray&) = delete;
	PetscReadArray& operator=(const PetscReadArray&) = delete;
	~PetscReadArray()
	{
		if (vector_) VecRestoreArrayRead(vector_, &values_);
	}
	void Acquire(Vec vector)
	{
		if (vector_) throw std::logic_error("PETSc read view is already acquired");
		if (VecGetArrayRead(vector, &values_))
			throw std::runtime_error("VecGetArrayRead failed");
		vector_ = vector;
	}
	const PetscScalar* Data() const { return values_; }
	void Restore()
	{
		if (!vector_) return;
		if (VecRestoreArrayRead(vector_, &values_))
			throw std::runtime_error("VecRestoreArrayRead failed");
		vector_ = nullptr;
	}
private:
	Vec vector_ = nullptr;
	const PetscScalar* values_ = nullptr;
};

} // namespace iga

#endif
