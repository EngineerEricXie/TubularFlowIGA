#ifndef NATIVE_TET_FEM_HPP
#define NATIVE_TET_FEM_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <istream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace iga {

struct NativeTetTriangle
{
	std::uint64_t id = 0;
	std::array<std::uint32_t, 3> nodes{};
	int boundary_label = -1;
};

struct NativeTetCell
{
	std::uint64_t id = 0;
	std::array<std::uint32_t, 4> nodes{};
};

struct NativeTetMesh
{
	std::vector<std::array<double, 3>> points;
	std::vector<NativeTetTriangle> boundary_triangles;
	std::vector<NativeTetCell> cells;
};

namespace native_tet_detail {

inline std::vector<std::string> Tokens(const std::string& line)
{
	std::istringstream input(line);
	std::vector<std::string> result;
	for (std::string value; input >> value;) result.push_back(value);
	return result;
}

inline std::uint64_t Unsigned(const std::string& token, const char* field)
{
	std::size_t used = 0;
	unsigned long long value = 0;
	try {
		value = std::stoull(token, &used);
	} catch (const std::exception&) {
		throw std::runtime_error(std::string("invalid Gmsh ")+field);
	}
	if (used != token.size()) throw std::runtime_error(std::string("invalid Gmsh ")+field);
	return static_cast<std::uint64_t>(value);
}

inline int Integer(const std::string& token, const char* field)
{
	std::size_t used = 0;
	long long value = 0;
	try {
		value = std::stoll(token, &used);
	} catch (const std::exception&) {
		throw std::runtime_error(std::string("invalid Gmsh ")+field);
	}
	if (used != token.size() || value < std::numeric_limits<int>::min()
		|| value > std::numeric_limits<int>::max())
		throw std::runtime_error(std::string("invalid Gmsh ")+field);
	return static_cast<int>(value);
}

inline double Real(const std::string& token, const char* field)
{
	std::size_t used = 0;
	double value = 0.0;
	try {
		value = std::stod(token, &used);
	} catch (const std::exception&) {
		throw std::runtime_error(std::string("invalid Gmsh ")+field);
	}
	if (used != token.size() || !std::isfinite(value))
		throw std::runtime_error(std::string("nonfinite or invalid Gmsh ")+field);
	return value;
}

inline std::string Line(std::istream& input, const char* context)
{
	std::string line;
	if (!std::getline(input, line))
		throw std::runtime_error(std::string("truncated Gmsh ")+context);
	return line;
}

inline std::size_t ElementNodeCount(int type)
{
	switch (type) {
	case 15: return 1;
	case 1: return 2;
	case 2: return 3;
	case 4: return 4;
	default: throw std::runtime_error("unsupported Gmsh element type "+std::to_string(type));
	}
}

inline std::array<std::uint32_t, 3> SortedFace(std::uint32_t a, std::uint32_t b,
	std::uint32_t c)
{
	std::array<std::uint32_t, 3> face{{a, b, c}};
	std::sort(face.begin(), face.end());
	return face;
}

} // namespace native_tet_detail

inline NativeTetMesh ReadNativeTetMeshGmsh41(std::istream& input)
{
	using namespace native_tet_detail;
	std::map<std::pair<int, int>, std::string> physical_names;
	std::map<std::pair<int, int>, int> entity_physical;
	std::unordered_map<std::uint64_t, std::array<double, 3>> tagged_points;
	struct TaggedElement { std::uint64_t id; int dimension; int entity; int type;
		std::vector<std::uint64_t> nodes; };
	std::vector<TaggedElement> tagged_elements;
	bool mesh_format = false;
	for (std::string section; std::getline(input, section);) {
		if (section.empty()) continue;
		if (section == "$MeshFormat") {
			const auto tokens = Tokens(Line(input, "mesh format"));
			if (tokens.size() != 3 || tokens[0] != "4.1" || tokens[1] != "0"
				|| tokens[2] != "8" || Line(input, "mesh format end") != "$EndMeshFormat")
				throw std::runtime_error("native FEM requires ASCII Gmsh 4.1 with 8-byte reals");
			mesh_format = true;
		} else if (section == "$PhysicalNames") {
			const auto count = Unsigned(Line(input, "physical-name count"), "physical-name count");
			for (std::uint64_t i = 0; i < count; ++i) {
				const std::string line = Line(input, "physical name");
				std::istringstream row(line);
				int dimension = -1, tag = -1;
				if (!(row >> dimension >> tag)) throw std::runtime_error("invalid Gmsh physical name");
				const auto first = line.find('"');
				const auto last = line.rfind('"');
				if (first == std::string::npos || last == first)
					throw std::runtime_error("Gmsh physical name must be quoted");
				if (!physical_names.emplace(std::make_pair(dimension, tag),
					line.substr(first+1, last-first-1)).second)
					throw std::runtime_error("duplicate Gmsh physical name");
			}
			if (Line(input, "physical-name end") != "$EndPhysicalNames")
				throw std::runtime_error("invalid Gmsh physical-name terminator");
		} else if (section == "$Entities") {
			const auto counts = Tokens(Line(input, "entity counts"));
			if (counts.size() != 4) throw std::runtime_error("invalid Gmsh entity counts");
			for (int dimension = 0; dimension < 4; ++dimension) {
				const auto count = Unsigned(counts[static_cast<std::size_t>(dimension)], "entity count");
				for (std::uint64_t i = 0; i < count; ++i) {
					const auto tokens = Tokens(Line(input, "entity"));
					const std::size_t physical_count_index = dimension == 0 ? 4 : 7;
					if (tokens.size() <= physical_count_index)
						throw std::runtime_error("invalid Gmsh entity record");
					const int entity = Integer(tokens[0], "entity tag");
					const auto physical_count = Unsigned(tokens[physical_count_index], "physical tag count");
					if (tokens.size() < physical_count_index+1+physical_count)
						throw std::runtime_error("truncated Gmsh entity physical tags");
					if (dimension >= 2) {
						if (physical_count != 1)
							throw std::runtime_error("every Gmsh surface and volume entity must have one physical tag");
						entity_physical[{dimension, entity}] = Integer(
							tokens[physical_count_index+1], "physical tag");
					}
				}
			}
			if (Line(input, "entity end") != "$EndEntities")
				throw std::runtime_error("invalid Gmsh entity terminator");
		} else if (section == "$Nodes") {
			const auto header = Tokens(Line(input, "node header"));
			if (header.size() != 4) throw std::runtime_error("invalid Gmsh node header");
			const auto blocks = Unsigned(header[0], "node block count");
			const auto total = Unsigned(header[1], "node count");
			for (std::uint64_t block = 0; block < blocks; ++block) {
				const auto block_header = Tokens(Line(input, "node block"));
				if (block_header.size() != 4) throw std::runtime_error("invalid Gmsh node block");
				const int dimension = Integer(block_header[0], "node entity dimension");
				const bool parametric = Integer(block_header[2], "parametric flag") != 0;
				const auto count = Unsigned(block_header[3], "node block size");
				std::vector<std::uint64_t> tags(static_cast<std::size_t>(count));
				for (auto& tag : tags) tag = Unsigned(Line(input, "node tag"), "node tag");
				for (const auto tag : tags) {
					const auto values = Tokens(Line(input, "node coordinates"));
					const std::size_t expected = 3+(parametric ? static_cast<std::size_t>(dimension) : 0);
					if (values.size() != expected) throw std::runtime_error("invalid Gmsh node coordinates");
					std::array<double, 3> point{{Real(values[0], "coordinate"),
						Real(values[1], "coordinate"), Real(values[2], "coordinate")}};
					if (!tagged_points.emplace(tag, point).second)
						throw std::runtime_error("duplicate Gmsh node tag");
				}
			}
			if (tagged_points.size() != total || Line(input, "node end") != "$EndNodes")
				throw std::runtime_error("Gmsh node count mismatch or invalid terminator");
		} else if (section == "$Elements") {
			const auto header = Tokens(Line(input, "element header"));
			if (header.size() != 4) throw std::runtime_error("invalid Gmsh element header");
			const auto blocks = Unsigned(header[0], "element block count");
			const auto total = Unsigned(header[1], "element count");
			for (std::uint64_t block = 0; block < blocks; ++block) {
				const auto block_header = Tokens(Line(input, "element block"));
				if (block_header.size() != 4) throw std::runtime_error("invalid Gmsh element block");
				const int dimension = Integer(block_header[0], "element entity dimension");
				const int entity = Integer(block_header[1], "element entity tag");
				const int type = Integer(block_header[2], "element type");
				const auto count = Unsigned(block_header[3], "element block size");
				const auto node_count = ElementNodeCount(type);
				for (std::uint64_t i = 0; i < count; ++i) {
					const auto values = Tokens(Line(input, "element"));
					if (values.size() != node_count+1) throw std::runtime_error("invalid Gmsh element record");
					TaggedElement element{Unsigned(values[0], "element tag"), dimension, entity, type, {}};
					for (std::size_t node = 0; node < node_count; ++node)
						element.nodes.push_back(Unsigned(values[node+1], "element node tag"));
					tagged_elements.push_back(std::move(element));
				}
			}
			if (tagged_elements.size() != total || Line(input, "element end") != "$EndElements")
				throw std::runtime_error("Gmsh element count mismatch or invalid terminator");
		} else if (!section.empty() && section.front() == '$') {
			throw std::runtime_error("unsupported Gmsh section "+section);
		} else {
			throw std::runtime_error("unexpected text outside a Gmsh section");
		}
	}
	if (!mesh_format || tagged_points.empty() || tagged_elements.empty())
		throw std::runtime_error("incomplete Gmsh mesh");
	NativeTetMesh mesh;
	std::unordered_map<std::uint64_t, std::uint32_t> node_index;
	std::vector<std::uint64_t> node_tags;
	node_tags.reserve(tagged_points.size());
	for (const auto& item : tagged_points) node_tags.push_back(item.first);
	std::sort(node_tags.begin(), node_tags.end());
	mesh.points.reserve(node_tags.size());
	for (const auto tag : node_tags) {
		if (mesh.points.size() > std::numeric_limits<std::uint32_t>::max())
			throw std::runtime_error("Gmsh mesh has too many nodes");
		node_index[tag] = static_cast<std::uint32_t>(mesh.points.size());
		mesh.points.push_back(tagged_points.at(tag));
	}
	for (const auto& tagged : tagged_elements) {
		if (tagged.type != 2 && tagged.type != 4) continue;
		const auto physical = entity_physical.find({tagged.dimension, tagged.entity});
		if (physical == entity_physical.end()) throw std::runtime_error("element entity has no physical tag");
		if (tagged.type == 2) {
			const auto name = physical_names.find({2, physical->second});
			const std::string prefix = "boundary_label_";
			if (name == physical_names.end() || name->second.compare(0, prefix.size(), prefix) != 0)
				throw std::runtime_error("surface physical name is not boundary_label_<integer>");
			const int label = Integer(name->second.substr(prefix.size()), "boundary label");
			NativeTetTriangle triangle;
			triangle.id = tagged.id;
			triangle.boundary_label = label;
			for (std::size_t i = 0; i < 3; ++i) triangle.nodes[i] = node_index.at(tagged.nodes[i]);
			mesh.boundary_triangles.push_back(triangle);
		} else {
			const auto name = physical_names.find({3, physical->second});
			if (name == physical_names.end() || name->second != "fluid")
				throw std::runtime_error("volume physical name must be fluid");
			NativeTetCell cell;
			cell.id = tagged.id;
			for (std::size_t i = 0; i < 4; ++i) cell.nodes[i] = node_index.at(tagged.nodes[i]);
			mesh.cells.push_back(cell);
		}
	}
	if (mesh.cells.empty() || mesh.boundary_triangles.empty())
		throw std::runtime_error("Gmsh mesh has no tetrahedra or labelled boundary triangles");
	std::map<std::array<std::uint32_t, 3>, unsigned> face_uses;
	for (const auto& cell : mesh.cells) {
		const auto& n = cell.nodes;
		for (const auto& face : {SortedFace(n[0],n[1],n[2]), SortedFace(n[0],n[1],n[3]),
			SortedFace(n[0],n[2],n[3]), SortedFace(n[1],n[2],n[3])}) ++face_uses[face];
	}
	std::set<std::array<std::uint32_t, 3>> labelled;
	for (const auto& triangle : mesh.boundary_triangles)
		if (!labelled.insert(SortedFace(triangle.nodes[0], triangle.nodes[1], triangle.nodes[2])).second)
			throw std::runtime_error("duplicate labelled Gmsh boundary triangle");
	std::set<std::array<std::uint32_t, 3>> boundary;
	for (const auto& item : face_uses) {
		if (item.second == 1) boundary.insert(item.first);
		else if (item.second != 2) throw std::runtime_error("non-manifold tetrahedral face");
	}
	if (boundary != labelled)
		throw std::runtime_error("labelled triangles do not equal tetrahedral boundary");
	return mesh;
}

struct NativeTaylorHoodTopology
{
	std::vector<std::array<std::uint32_t, 2>> edges;
	std::vector<std::array<std::uint32_t, 10>> cell_velocity_nodes;
	std::map<int, std::set<std::uint32_t>> boundary_velocity_nodes;
	std::map<int, std::set<std::uint32_t>> boundary_pressure_nodes;
};

inline NativeTaylorHoodTopology BuildNativeTaylorHoodTopology(const NativeTetMesh& mesh)
{
	NativeTaylorHoodTopology result;
	std::map<std::array<std::uint32_t, 2>, std::uint32_t> edge_indices;
	const std::array<std::array<int, 2>, 6> local_edges{{{{0,1}},{{0,2}},{{0,3}},
		{{1,2}},{{1,3}},{{2,3}}}};
	for (const auto& cell : mesh.cells) {
		std::array<std::uint32_t, 10> nodes{};
		std::copy(cell.nodes.begin(), cell.nodes.end(), nodes.begin());
		for (std::size_t local = 0; local < local_edges.size(); ++local) {
			std::array<std::uint32_t, 2> edge{{cell.nodes[local_edges[local][0]],
				cell.nodes[local_edges[local][1]]}};
			std::sort(edge.begin(), edge.end());
			auto inserted = edge_indices.emplace(edge, static_cast<std::uint32_t>(edge_indices.size()));
			nodes[4+local] = static_cast<std::uint32_t>(mesh.points.size())+inserted.first->second;
		}
		result.cell_velocity_nodes.push_back(nodes);
	}
	result.edges.resize(edge_indices.size());
	for (const auto& item : edge_indices) result.edges[item.second] = item.first;
	for (const auto& triangle : mesh.boundary_triangles) {
		auto& velocity = result.boundary_velocity_nodes[triangle.boundary_label];
		auto& pressure = result.boundary_pressure_nodes[triangle.boundary_label];
		for (const auto node : triangle.nodes) { velocity.insert(node); pressure.insert(node); }
		for (const auto pair : {std::array<int,2>{{0,1}}, std::array<int,2>{{0,2}},
			std::array<int,2>{{1,2}}}) {
			std::array<std::uint32_t, 2> edge{{triangle.nodes[pair[0]], triangle.nodes[pair[1]]}};
			std::sort(edge.begin(), edge.end());
			const auto found = edge_indices.find(edge);
			if (found == edge_indices.end()) throw std::runtime_error("boundary edge is absent from cells");
			velocity.insert(static_cast<std::uint32_t>(mesh.points.size())+found->second);
		}
	}
	return result;
}

struct NativeTaylorHoodBasis
{
	std::array<double, 10> velocity{};
	std::array<std::array<double, 3>, 10> velocity_gradient_reference{};
	std::array<double, 4> pressure{};
};

inline NativeTaylorHoodBasis EvaluateNativeTaylorHoodBasis(double r, double s, double t)
{
	const std::array<double, 4> lambda{{1.0-r-s-t, r, s, t}};
	const std::array<std::array<double, 3>, 4> gradient{{{{-1,-1,-1}},{{1,0,0}},
		{{0,1,0}},{{0,0,1}}}};
	const std::array<std::array<int, 2>, 6> edges{{{{0,1}},{{0,2}},{{0,3}},
		{{1,2}},{{1,3}},{{2,3}}}};
	NativeTaylorHoodBasis result;
	result.pressure = lambda;
	for (std::size_t i = 0; i < 4; ++i) {
		result.velocity[i] = lambda[i]*(2.0*lambda[i]-1.0);
		for (int component = 0; component < 3; ++component)
			result.velocity_gradient_reference[i][component] =
				(4.0*lambda[i]-1.0)*gradient[i][component];
	}
	for (std::size_t e = 0; e < edges.size(); ++e) {
		const auto first = static_cast<std::size_t>(edges[e][0]);
		const auto second = static_cast<std::size_t>(edges[e][1]);
		result.velocity[4+e] = 4.0*lambda[first]*lambda[second];
		for (int component = 0; component < 3; ++component)
			result.velocity_gradient_reference[4+e][component] = 4.0*(
				lambda[first]*gradient[second][component]
				+lambda[second]*gradient[first][component]);
	}
	return result;
}

struct NativeTetGeometry
{
	double determinant = 0.0;
	std::array<std::array<double, 3>, 4> barycentric_gradients{};
};

inline NativeTetGeometry EvaluateNativeTetGeometry(const NativeTetMesh& mesh,
	const NativeTetCell& cell)
{
	const auto& x0 = mesh.points.at(cell.nodes[0]);
	const auto& x1 = mesh.points.at(cell.nodes[1]);
	const auto& x2 = mesh.points.at(cell.nodes[2]);
	const auto& x3 = mesh.points.at(cell.nodes[3]);
	const std::array<std::array<double, 3>, 3> jacobian{{
		{{x1[0]-x0[0], x2[0]-x0[0], x3[0]-x0[0]}},
		{{x1[1]-x0[1], x2[1]-x0[1], x3[1]-x0[1]}},
		{{x1[2]-x0[2], x2[2]-x0[2], x3[2]-x0[2]}}}};
	const double determinant =
		jacobian[0][0]*(jacobian[1][1]*jacobian[2][2]-jacobian[1][2]*jacobian[2][1])
		-jacobian[0][1]*(jacobian[1][0]*jacobian[2][2]-jacobian[1][2]*jacobian[2][0])
		+jacobian[0][2]*(jacobian[1][0]*jacobian[2][1]-jacobian[1][1]*jacobian[2][0]);
	if (!(determinant > 0.0) || !std::isfinite(determinant))
		throw std::runtime_error("native FEM tetrahedron has non-positive Jacobian");
	std::array<std::array<double, 3>, 3> inverse{{
		{{(jacobian[1][1]*jacobian[2][2]-jacobian[1][2]*jacobian[2][1])/determinant,
			(jacobian[0][2]*jacobian[2][1]-jacobian[0][1]*jacobian[2][2])/determinant,
			(jacobian[0][1]*jacobian[1][2]-jacobian[0][2]*jacobian[1][1])/determinant}},
		{{(jacobian[1][2]*jacobian[2][0]-jacobian[1][0]*jacobian[2][2])/determinant,
			(jacobian[0][0]*jacobian[2][2]-jacobian[0][2]*jacobian[2][0])/determinant,
			(jacobian[0][2]*jacobian[1][0]-jacobian[0][0]*jacobian[1][2])/determinant}},
		{{(jacobian[1][0]*jacobian[2][1]-jacobian[1][1]*jacobian[2][0])/determinant,
			(jacobian[0][1]*jacobian[2][0]-jacobian[0][0]*jacobian[2][1])/determinant,
			(jacobian[0][0]*jacobian[1][1]-jacobian[0][1]*jacobian[1][0])/determinant}}}};
	NativeTetGeometry result;
	result.determinant = determinant;
	const std::array<std::array<double, 3>, 4> reference{{{{-1,-1,-1}},{{1,0,0}},
		{{0,1,0}},{{0,0,1}}}};
	for (std::size_t basis = 0; basis < 4; ++basis)
		for (int component = 0; component < 3; ++component)
			for (int reference_component = 0; reference_component < 3; ++reference_component)
				result.barycentric_gradients[basis][component] +=
					inverse[reference_component][component]*reference[basis][reference_component];
	return result;
}

struct NativeTetQuadraturePoint
{
	std::array<double, 3> reference{};
	double weight = 0.0;
};

inline std::vector<NativeTetQuadraturePoint> NativeTetDegreeFiveQuadrature()
{
	const std::array<double, 4> abscissa{{
		0.06943184420297371, 0.33000947820757187,
		0.6699905217924281, 0.9305681557970262}};
	const std::array<double, 4> weight{{
		0.17392742256872693, 0.32607257743127307,
		0.32607257743127307, 0.17392742256872693}};
	std::vector<NativeTetQuadraturePoint> result;
	result.reserve(64);
	for (std::size_t i = 0; i < 4; ++i)
		for (std::size_t j = 0; j < 4; ++j)
			for (std::size_t k = 0; k < 4; ++k) {
				const double a = abscissa[i], b = abscissa[j], c = abscissa[k];
				result.push_back({{{a, (1.0-a)*b, (1.0-a)*(1.0-b)*c}},
					weight[i]*weight[j]*weight[k]*(1.0-a)*(1.0-a)*(1.0-b)});
			}
	return result;
}

struct NativeTaylorHoodPhysicalBasis
{
	std::array<double, 10> velocity{};
	std::array<std::array<double, 3>, 10> velocity_gradients{};
	std::array<double, 4> pressure{};
};

inline NativeTaylorHoodPhysicalBasis EvaluateNativeTaylorHoodPhysicalBasis(
	const NativeTetGeometry& geometry, const std::array<double, 3>& reference)
{
	const auto basis = EvaluateNativeTaylorHoodBasis(reference[0], reference[1], reference[2]);
	NativeTaylorHoodPhysicalBasis result;
	result.velocity = basis.velocity;
	result.pressure = basis.pressure;
	for (std::size_t a = 0; a < 10; ++a)
		for (int component = 0; component < 3; ++component)
			for (int reference_component = 0; reference_component < 3; ++reference_component)
				result.velocity_gradients[a][component] +=
					basis.velocity_gradient_reference[a][reference_component]
					*geometry.barycentric_gradients[reference_component+1][component];
	return result;
}

struct NativeNavierStokesParameters
{
	double density = 0.0;
	double dynamic_viscosity = 0.0;
	// Per-component affine body acceleration in SI units. Coefficients are
	// ordered as constant, x, y, z and default to zero.
	std::array<std::array<double,4>,3> body_acceleration_m_s2{};
};

struct NativeTaylorHoodElementSystem
{
	static constexpr std::size_t dofs = 34;
	std::array<double, dofs> residual{};
	std::array<double, dofs*dofs> jacobian{};
};

inline NativeTaylorHoodElementSystem BuildNativeTaylorHoodAleNavierStokesElement(
	const NativeTetMesh& mesh, const NativeTetCell& cell,
	const std::array<double, NativeTaylorHoodElementSystem::dofs>& state,
	const std::array<std::array<double,3>,4>& mesh_velocity,
	const NativeNavierStokesParameters& parameters)
{
	if (!(parameters.density > 0.0) || !std::isfinite(parameters.density)
		|| !(parameters.dynamic_viscosity > 0.0)
		|| !std::isfinite(parameters.dynamic_viscosity))
		throw std::runtime_error("native FEM Navier-Stokes parameters must be finite and positive");
	for(const auto& component:parameters.body_acceleration_m_s2)for(double coefficient:component)
		if(!std::isfinite(coefficient))
			throw std::runtime_error("native FEM body acceleration must be finite");
	const auto geometry = EvaluateNativeTetGeometry(mesh, cell);
	NativeTaylorHoodElementSystem result;
	for (const auto& quadrature : NativeTetDegreeFiveQuadrature()) {
		const auto basis = EvaluateNativeTaylorHoodPhysicalBasis(geometry, quadrature.reference);
		std::array<double, 3> velocity{{0,0,0}};
		std::array<double, 3> grid_velocity{{0,0,0}};
		std::array<double,3> point{},body_acceleration{};
		std::array<std::array<double, 3>, 3> velocity_gradient{};
		double pressure = 0.0;
		for (std::size_t a = 0; a < 10; ++a)
			for (int component = 0; component < 3; ++component) {
				velocity[component] += basis.velocity[a]*state[3*a+component];
				for (int derivative = 0; derivative < 3; ++derivative)
					velocity_gradient[component][derivative] +=
						basis.velocity_gradients[a][derivative]*state[3*a+component];
			}
		for (std::size_t a = 0; a < 4; ++a) pressure += basis.pressure[a]*state[30+a];
		for(std::size_t a=0;a<4;++a)for(int component=0;component<3;++component)
			point[component]+=basis.pressure[a]*mesh.points[cell.nodes[a]][component];
		for(int component=0;component<3;++component)body_acceleration[component]=
			parameters.body_acceleration_m_s2[component][0]
			+parameters.body_acceleration_m_s2[component][1]*point[0]
			+parameters.body_acceleration_m_s2[component][2]*point[1]
			+parameters.body_acceleration_m_s2[component][3]*point[2];
		for (std::size_t a = 0; a < 4; ++a)
			for (int component = 0; component < 3; ++component)
				grid_velocity[component] += basis.pressure[a]*mesh_velocity[a][component];
		std::array<double,3> relative_velocity{};
		for(int component=0;component<3;++component)
			relative_velocity[component]=velocity[component]-grid_velocity[component];
		const double weighted_volume = quadrature.weight*geometry.determinant;
		for (std::size_t a = 0; a < 10; ++a)
			for (int component = 0; component < 3; ++component) {
				const std::size_t row = 3*a+component;
				double viscous = 0.0, convection = 0.0;
				for (int derivative = 0; derivative < 3; ++derivative) {
					viscous += parameters.dynamic_viscosity*(
						velocity_gradient[component][derivative]
						+velocity_gradient[derivative][component])*basis.velocity_gradients[a][derivative];
					convection += relative_velocity[derivative]*velocity_gradient[component][derivative];
				}
				result.residual[row] += weighted_volume*(viscous
					+parameters.density*basis.velocity[a]*convection
					-pressure*basis.velocity_gradients[a][component]
					-parameters.density*basis.velocity[a]*body_acceleration[component]);
				for (std::size_t c = 0; c < 10; ++c)
					for (int unknown_component = 0; unknown_component < 3; ++unknown_component) {
						const std::size_t column = 3*c+unknown_component;
						double tangent = parameters.dynamic_viscosity*(
							(component == unknown_component ?
								basis.velocity_gradients[c][0]*basis.velocity_gradients[a][0]
								+basis.velocity_gradients[c][1]*basis.velocity_gradients[a][1]
								+basis.velocity_gradients[c][2]*basis.velocity_gradients[a][2] : 0.0)
							+basis.velocity_gradients[c][component]
								*basis.velocity_gradients[a][unknown_component]);
						tangent += parameters.density*basis.velocity[a]*(
							basis.velocity[c]*velocity_gradient[component][unknown_component]
							+(component == unknown_component ?
								relative_velocity[0]*basis.velocity_gradients[c][0]
								+relative_velocity[1]*basis.velocity_gradients[c][1]
								+relative_velocity[2]*basis.velocity_gradients[c][2] : 0.0));
						result.jacobian[row*NativeTaylorHoodElementSystem::dofs+column] +=
							weighted_volume*tangent;
					}
				for (std::size_t c = 0; c < 4; ++c)
					result.jacobian[row*NativeTaylorHoodElementSystem::dofs+30+c] -=
						weighted_volume*basis.pressure[c]*basis.velocity_gradients[a][component];
			}
		for (std::size_t a = 0; a < 4; ++a) {
			const std::size_t row = 30+a;
			double divergence = 0.0;
			for (int component = 0; component < 3; ++component)
				divergence += velocity_gradient[component][component];
			result.residual[row] += weighted_volume*basis.pressure[a]*divergence;
			for (std::size_t c = 0; c < 10; ++c)
				for (int component = 0; component < 3; ++component)
					result.jacobian[row*NativeTaylorHoodElementSystem::dofs+3*c+component] +=
						weighted_volume*basis.pressure[a]*basis.velocity_gradients[c][component];
		}
	}
	return result;
}

inline NativeTaylorHoodElementSystem BuildNativeTaylorHoodNavierStokesElement(
	const NativeTetMesh& mesh, const NativeTetCell& cell,
	const std::array<double, NativeTaylorHoodElementSystem::dofs>& state,
	const NativeNavierStokesParameters& parameters)
{
	return BuildNativeTaylorHoodAleNavierStokesElement(mesh,cell,state,{},parameters);
}

} // namespace iga

#endif
