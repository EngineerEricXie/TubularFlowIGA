#ifndef IGA_TRANSPORT_BUDGET_HPP
#define IGA_TRANSPORT_BUDGET_HPP

#include "BoundaryFlow.hpp"
#include "GenericCaseInput.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <set>
#include <stdexcept>
#include <vector>

namespace iga {

struct TransportFieldBudget {
	long double previous_integral = 0.0L;
	long double current_integral = 0.0L;
	long double capacity_rate = 0.0L;
	long double current_capacity_rate_scale = 0.0L;
	long double previous_capacity_rate_scale = 0.0L;
	long double advection_volume = 0.0L;
	long double advection_outward = 0.0L;
	long double advection_absolute_flux = 0.0L;
	long double velocity_divergence = 0.0L;
	long double reaction = 0.0L;
	long double source = 0.0L;
	long double natural_outward = 0.0L;
	long double natural_absolute_flux = 0.0L;
	long double essential_inward = 0.0L;
	long double essential_absolute_exchange = 0.0L;
	long double volume_defect = 0.0L;
	long double surface_defect = 0.0L;
	long double volume_scale = 0.0L;
	long double surface_scale = 0.0L;
	long double free_residual_l2 = 0.0L;
	long double free_component_scale = 0.0L;
	long double essential_error_l2 = 0.0L;
	long double essential_reference_l2 = 0.0L;
	std::size_t free_rows = 0;
};

inline bool TransportBudgetWithin(long double defect, long double scale,
	double relative_tolerance, double absolute_tolerance)
{
	if (!std::isfinite(defect) || !std::isfinite(scale) || scale < 0.0L
		|| !std::isfinite(relative_tolerance) || relative_tolerance < 0.0
		|| !std::isfinite(absolute_tolerance) || absolute_tolerance < 0.0)
		throw std::invalid_argument("invalid transport budget metric or tolerance");
	return scale > 0.0L ? std::abs(defect)/scale <= relative_tolerance
		: std::abs(defect) <= absolute_tolerance;
}

// Independent field evaluation: interpolate c/grad(c) first, then integrate
// scalar terms and unconstrained test residuals. No production element matrix,
// PETSc matrix, or production residual is consumed. Geometry/basis primitives
// remain shared, so separate analytic geometry tests are still necessary.
class TransportBudgetAccumulator {
public:
	TransportBudgetAccumulator(std::size_t nodes, std::size_t elements,
		const CompiledLinearSystem& system, const SimulationConfiguration& configuration,
		const ResolvedScalarBoundaries& boundaries, const std::vector<double>& previous,
		const std::vector<double>& current, const std::vector<std::array<double, 3>>& velocity)
		: nodes_(nodes), system_(system), configuration_(configuration), boundaries_(boundaries),
		  previous_(previous), current_(current), velocity_(velocity), seen_(elements, false)
	{
		const auto fields = system_.fields.size();
		if (!nodes || !elements || !fields || nodes > std::numeric_limits<std::size_t>::max()/fields
			|| !(system_.dt > 0.0) || !std::isfinite(system_.dt))
			throw std::invalid_argument("invalid transport budget dimensions or dt");
		const auto rows = nodes*fields;
		if (previous_.size() != rows || current_.size() != rows || velocity_.size() != nodes
			|| boundaries_.constrained.size() != rows || boundaries_.value.size() != rows)
			throw std::invalid_argument("transport budget field dimensions differ");
		if (system_.velocity_source != "prescribed")
			throw std::invalid_argument("transport budget requires prescribed nodal velocity");
		for (const auto* values : {&previous_, &current_, &boundaries_.value})
			for (const auto value : *values) RequireFinite(value);
		for (const auto& value : velocity_) for (const auto component : value) RequireFinite(component);
		rows_.resize(rows);
		budgets_.resize(fields);
		supg_.assign(fields, false);
		for (const auto& stabilization : system_.stabilization) {
			const auto field = system_.field_index.find(stabilization.equation);
			if (field == system_.field_index.end() || field->second >= fields
				|| stabilization.method != "supg" || stabilization.velocity != "prescribed")
				throw std::invalid_argument("unsupported transport budget stabilization");
			supg_[field->second] = true;
		}
		for (const auto& term : system_.terms) {
			RequireFinite(term.coefficient);
			if (term.equation >= fields || term.trial >= fields)
				throw std::invalid_argument("transport budget term references an invalid field");
			if (term.kind == TermKind::Diffusion && supg_[term.equation])
				throw std::invalid_argument("transport budget does not support SUPG diffusion");
			if (term.kind == TermKind::Advection && term.velocity != "prescribed")
				throw std::invalid_argument("transport budget advection requires prescribed velocity");
		}
		for (const auto& boundary : configuration_.boundaries)
			for (const auto& condition : boundary.conditions) {
				if (!condition.waveform.empty())
					throw std::invalid_argument("transport budget requires static boundary values");
				for (const auto value : condition.value) RequireFinite(value);
				RequireFinite(condition.coefficient);
				RequireFinite(condition.exterior_value);
			}
	}

	void Add(const Element& element)
	{
		if (failed_) throw std::logic_error("transport budget accumulation previously failed");
		failed_ = true;
		if (element.id >= seen_.size() || seen_[element.id])
			throw std::invalid_argument("duplicate or out-of-range budget element");
		const auto nen = element.connectivity.size();
		if (!nen || element.extraction.size() != nen)
			throw std::invalid_argument("invalid budget element basis layout");
		std::set<std::int32_t> unique;
		for (const auto node : element.connectivity)
			if (node < 0 || static_cast<std::size_t>(node) >= nodes_ || !unique.insert(node).second)
				throw std::invalid_argument("invalid or duplicate budget connectivity");
		for (const auto& row : element.extraction) for (const auto value : row) RequireFinite(value);
		for (const auto& point : element.bezier_points) for (const auto value : point) RequireFinite(value);
		seen_[element.id] = true;
		const auto fields = system_.fields.size();
		FullCell4x4x4VolumeQuadratureProvider volume(element);
		for (const auto& point : volume.Rule().Points()) {
			const auto basis = EvaluateBoundaryBasis(element, point.parametric[0], point.parametric[1], point.parametric[2]);
			const long double measure = point.weight*basis.raw_determinant;
			std::vector<long double> current(fields), previous(fields);
			std::vector<std::array<long double, 3>> gradient(fields);
			std::array<long double, 3> velocity{};
			long double divergence = 0.0L;
			for (std::size_t a = 0; a < nen; ++a) {
				const auto node = static_cast<std::size_t>(element.connectivity[a]);
				for (std::size_t field = 0; field < fields; ++field) {
					current[field] += basis.value[a]*current_[node*fields+field];
					previous[field] += basis.value[a]*previous_[node*fields+field];
					for (int d = 0; d < 3; ++d)
						gradient[field][d] += basis.gradient[a][d]*current_[node*fields+field];
				}
				for (int d = 0; d < 3; ++d) {
					velocity[d] += basis.value[a]*velocity_[node][d];
					divergence += basis.gradient[a][d]*velocity_[node][d];
				}
			}
			std::vector<long double> streamline(nen);
			long double inverse_length = 0.0L;
			for (std::size_t a = 0; a < nen; ++a) {
				for (int d = 0; d < 3; ++d) streamline[a] += velocity[d]*basis.gradient[a][d];
				inverse_length += std::abs(streamline[a]);
			}
			const long double tau = inverse_length > 0.0L
				? 1.0L/std::hypot(inverse_length, 2.0L/system_.dt) : 0.0L;
			for (std::size_t field = 0; field < fields; ++field) {
				budgets_[field].previous_integral += measure*previous[field];
				budgets_[field].current_integral += measure*current[field];
			}
			for (const auto& term : system_.terms) {
				auto& budget = budgets_[term.equation];
				long double value = 0.0L;
				Component component = CapacityCurrent;
				switch (term.kind) {
				case TermKind::TimeDerivative:
					value = term.coefficient*(current[term.trial]-previous[term.trial])/system_.dt;
					budget.capacity_rate += measure*value;
					budget.current_capacity_rate_scale += measure*std::abs(term.coefficient*current[term.trial]/system_.dt);
					budget.previous_capacity_rate_scale += measure*std::abs(term.coefficient*previous[term.trial]/system_.dt);
					break;
				case TermKind::Advection:
					component = Advection;
					for (int d = 0; d < 3; ++d) value += term.coefficient*velocity[d]*gradient[term.trial][d];
					budget.advection_volume += measure*value;
					budget.velocity_divergence += measure*term.coefficient*current[term.trial]*divergence;
					break;
				case TermKind::LinearCoupling:
					component = Reaction;
					value = term.coefficient*current[term.trial];
					budget.reaction += measure*value;
					break;
				case TermKind::VolumeSource:
					component = Source;
					value = -term.coefficient;
					budget.source -= measure*value;
					break;
				case TermKind::Diffusion:
					component = Diffusion;
					break;
				default: throw std::invalid_argument("unknown transport budget term");
				}
				for (std::size_t a = 0; a < nen; ++a) {
					const auto row = static_cast<std::size_t>(element.connectivity[a])*fields+term.equation;
					const long double test = basis.value[a] + (supg_[term.equation] ? tau*streamline[a] : 0.0L);
					if (term.kind == TermKind::TimeDerivative) {
						rows_[row][CapacityCurrent] += measure*test*term.coefficient*current[term.trial]/system_.dt;
						rows_[row][CapacityPrevious] -= measure*test*term.coefficient*previous[term.trial]/system_.dt;
					} else if (component == Diffusion) {
						long double diffusive = 0.0L;
						for (int d = 0; d < 3; ++d) diffusive += basis.gradient[a][d]*gradient[term.trial][d];
						rows_[row][component] += measure*term.coefficient*diffusive;
					} else rows_[row][component] += measure*test*value;
				}
			}
		}
		BodyFittedSurface4x4QuadratureProvider surface(element);
		for (const auto& point : surface.Rule().Points()) {
			const auto boundary = std::find_if(configuration_.boundaries.begin(), configuration_.boundaries.end(),
				[&](const NamedBoundaryDefinition& value) { return value.label == point.boundary_id; });
			if (boundary == configuration_.boundaries.end())
				throw std::invalid_argument("budget surface has no configured boundary");
			const auto basis = EvaluateBoundaryBasis(element, point.parametric[0], point.parametric[1], point.parametric[2]);
			std::vector<long double> current(fields);
			long double normal_velocity = 0.0L;
			for (std::size_t a = 0; a < nen; ++a) {
				const auto node = static_cast<std::size_t>(element.connectivity[a]);
				for (std::size_t field = 0; field < fields; ++field)
					current[field] += basis.value[a]*current_[node*fields+field];
				for (int d = 0; d < 3; ++d)
					normal_velocity += basis.value[a]*velocity_[node][d]*point.normal[d];
			}
			for (const auto& term : system_.terms) if (term.kind == TermKind::Advection) {
				const auto flux = point.weight*term.coefficient*normal_velocity*current[term.trial];
				budgets_[term.equation].advection_outward += flux;
				budgets_[term.equation].advection_absolute_flux += std::abs(flux);
			}
			for (const auto& condition : boundary->conditions) {
				const auto field = system_.field_index.find(condition.field);
				if (field == system_.field_index.end()) continue;
				long double outward = 0.0L;
				if (condition.kind == FieldBoundaryKind::Flux) {
					if (condition.value.size() != 1) throw std::invalid_argument("invalid budget flux condition");
					outward = -condition.value[0];
				} else if (condition.kind == FieldBoundaryKind::Robin)
					outward = condition.coefficient*(current[field->second]-condition.exterior_value);
				else continue;
				budgets_[field->second].natural_outward += point.weight*outward;
				budgets_[field->second].natural_absolute_flux += point.weight*std::abs(outward);
				for (std::size_t a = 0; a < nen; ++a)
					rows_[static_cast<std::size_t>(element.connectivity[a])*fields+field->second][Natural]
						+= point.weight*basis.value[a]*outward;
			}
		}
		failed_ = false;
	}

	std::vector<TransportFieldBudget> Finish() const
	{
		if (failed_) throw std::logic_error("cannot publish a failed transport budget accumulation");
		if (std::find(seen_.begin(), seen_.end(), false) != seen_.end())
			throw std::invalid_argument("budget is missing database elements");
		auto result = budgets_;
		const auto fields = system_.fields.size();
		std::vector<std::array<long double, Components>> norms(fields);
		for (std::size_t row = 0; row < rows_.size(); ++row) {
			auto& budget = result[row%fields];
			long double residual = 0.0L;
			for (const auto value : rows_[row]) residual += value;
			if (boundaries_.constrained[row]) {
				budget.essential_inward += residual;
				budget.essential_absolute_exchange += std::abs(residual);
				budget.essential_error_l2 = std::hypot(budget.essential_error_l2,
					static_cast<long double>(current_[row])-boundaries_.value[row]);
				budget.essential_error_l2 = std::hypot(budget.essential_error_l2,
					static_cast<long double>(previous_[row])-boundaries_.value[row]);
				budget.essential_reference_l2 = std::hypot(budget.essential_reference_l2,
					static_cast<long double>(boundaries_.value[row]));
				budget.essential_reference_l2 = std::hypot(budget.essential_reference_l2,
					static_cast<long double>(boundaries_.value[row]));
			} else {
				++budget.free_rows;
				budget.free_residual_l2 = std::hypot(budget.free_residual_l2, residual);
				for (std::size_t component = 0; component < Components; ++component)
					norms[row%fields][component] = std::hypot(norms[row%fields][component], rows_[row][component]);
			}
		}
		for (std::size_t field = 0; field < fields; ++field) {
			auto& budget = result[field];
			for (const auto norm : norms[field]) budget.free_component_scale += norm;
			const auto common = budget.capacity_rate+budget.reaction-budget.source
				+budget.natural_outward-budget.essential_inward;
			const auto scale = budget.current_capacity_rate_scale+budget.previous_capacity_rate_scale
				+std::abs(budget.reaction)+std::abs(budget.source)
				+budget.natural_absolute_flux+budget.essential_absolute_exchange;
			budget.volume_defect = common+budget.advection_volume;
			budget.surface_defect = common+budget.advection_outward-budget.velocity_divergence;
			budget.volume_scale = scale+std::abs(budget.advection_volume);
			budget.surface_scale = scale+budget.advection_absolute_flux+std::abs(budget.velocity_divergence);
		}
		return result;
	}

private:
	enum Component : std::size_t { CapacityCurrent, CapacityPrevious, Advection, Reaction, Diffusion, Source, Natural, Components };
	static void RequireFinite(double value)
	{
		if (!std::isfinite(value)) throw std::invalid_argument("nonfinite transport budget input");
	}
	std::size_t nodes_;
	const CompiledLinearSystem& system_;
	const SimulationConfiguration& configuration_;
	const ResolvedScalarBoundaries& boundaries_;
	const std::vector<double>& previous_;
	const std::vector<double>& current_;
	const std::vector<std::array<double, 3>>& velocity_;
	std::vector<bool> seen_;
	std::vector<bool> supg_;
	std::vector<std::array<long double, Components>> rows_;
	std::vector<TransportFieldBudget> budgets_;
	bool failed_ = false;
};

} // namespace iga

#endif
