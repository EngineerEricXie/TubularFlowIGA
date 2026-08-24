#ifndef GENERIC_TRANSPORT_ELEMENT_HPP
#define GENERIC_TRANSPORT_ELEMENT_HPP

#include "SimulationConfig.hpp"
#include "TransportElement.hpp"

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace iga {

struct GenericTransportMatrices {
	std::vector<PetscScalar> left;
	std::vector<PetscScalar> previous;
	std::vector<PetscScalar> source;
};

inline GenericTransportMatrices BuildGenericTransportElement(const Element& element,
	const std::vector<std::array<double, 3>>& nodal_velocity, const CompiledLinearSystem& system)
{
	if (system.fields.empty()) throw std::runtime_error("linear transport system has no fields");
	if (!(system.dt > 0.0)) throw std::runtime_error("linear transport time step must be positive");
	const auto fields = system.fields.size();
	const auto nen = element.connectivity.size();
	const auto ndof = fields * nen;
	GenericTransportMatrices matrices{
		std::vector<PetscScalar>(ndof * ndof, 0.0),
		std::vector<PetscScalar>(ndof * ndof, 0.0),
		std::vector<PetscScalar>(ndof, 0.0)};
	std::vector<int> supg(fields, 0);
	for (const auto& definition : system.stabilization) {
		const auto found = system.field_index.find(definition.equation);
		if (found == system.field_index.end()) throw std::runtime_error("stabilization references an unknown equation");
		if (definition.method != "supg" || definition.velocity != "prescribed")
			throw std::runtime_error("CPU linear transport currently supports SUPG with velocity source 'prescribed'");
		supg[found->second] = 1;
	}
	for (const auto& term : system.terms)
		if ((term.kind == TermKind::Advection) && term.velocity != "prescribed")
			throw std::runtime_error("CPU linear transport currently supports advection velocity source 'prescribed'");

	constexpr std::array<double, 4> points{{0.06943184420297371, 0.33000947820757187, 0.6699905217924281, 0.9305681557970262}};
	constexpr std::array<double, 4> weights{{0.3478548451374539, 0.6521451548625461, 0.6521451548625461, 0.3478548451374539}};
	for (std::size_t qz = 0; qz < 4; ++qz)
		for (std::size_t qy = 0; qy < 4; ++qy)
			for (std::size_t qx = 0; qx < 4; ++qx) {
				auto basis = EvaluateBasis(element, points[qx], points[qy], points[qz]);
				const auto measure = weights[qx] * weights[qy] * weights[qz] * basis.determinant;
				std::array<double, 3> velocity{};
				for (std::size_t a = 0; a < nen; ++a)
					for (int d = 0; d < 3; ++d)
						velocity[d] += basis.value[a] * nodal_velocity.at(static_cast<std::size_t>(element.connectivity[a]))[d];
				double inverse_length = 0.0;
				for (std::size_t a = 0; a < nen; ++a)
					inverse_length += std::abs(velocity[0]*basis.gradient[a][0]
						+ velocity[1]*basis.gradient[a][1] + velocity[2]*basis.gradient[a][2]);
				const auto tau_space = inverse_length > 0.0 ? 1.0 / inverse_length : 0.0;
				const auto tau_time = system.dt / 2.0;
				const auto tau = tau_space > 0.0
					? 1.0 / std::sqrt(1.0/(tau_space*tau_space) + 1.0/(tau_time*tau_time)) : 0.0;
				std::vector<double> streamline(nen);
				for (std::size_t a = 0; a < nen; ++a)
					streamline[a] = velocity[0]*basis.gradient[a][0]
						+ velocity[1]*basis.gradient[a][1] + velocity[2]*basis.gradient[a][2];

				for (std::size_t a = 0; a < nen; ++a) {
					for (const auto& term : system.terms) {
						const double test = basis.value[a] + (supg[term.equation] ? tau*streamline[a] : 0.0);
						if (term.kind == TermKind::VolumeSource) {
							matrices.source[(a*fields)+term.equation] += system.dt * term.coefficient * test * measure;
							continue;
						}
						for (std::size_t b = 0; b < nen; ++b) {
							const auto row = a*fields + term.equation;
							const auto column = b*fields + term.trial;
							const auto index = row*ndof + column;
							const auto gradient_dot = basis.gradient[a][0]*basis.gradient[b][0]
								+ basis.gradient[a][1]*basis.gradient[b][1] + basis.gradient[a][2]*basis.gradient[b][2];
							if (term.kind == TermKind::TimeDerivative) {
								const auto mass = term.coefficient * test * basis.value[b] * measure;
								matrices.left[index] += mass;
								matrices.previous[index] += mass;
							} else if (term.kind == TermKind::Diffusion) {
								const auto weak_gradient = supg[term.equation] ? test*gradient_dot : gradient_dot;
								matrices.left[index] += system.dt * term.coefficient * weak_gradient * measure;
							} else if (term.kind == TermKind::Advection) {
								matrices.left[index] += system.dt * term.coefficient * test
									* (velocity[0]*basis.gradient[b][0] + velocity[1]*basis.gradient[b][1]
										+ velocity[2]*basis.gradient[b][2]) * measure;
							} else if (term.kind == TermKind::LinearCoupling) {
								matrices.left[index] += system.dt * term.coefficient * test * basis.value[b] * measure;
							}
						}
					}
				}
			}
	return matrices;
}

} // namespace iga

#endif
