#include "CrossSectionTemplate.hpp"
#include "CheckedText.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <utility>

namespace tubular {

namespace {
constexpr double CoreSideRadius = 0.70;
constexpr double CoreBow = 0.15;
}

CrossSectionTemplates GenerateCircularTemplates(double target_size)
{
    if (!std::isfinite(target_size) || target_size < 1.0/128.0 || target_size > 1.0)
        throw std::runtime_error("cross-section target size must be in [1/128,1]");
    // A bowed square core plus four annular blocks avoids collapsed corners.
    // Side midpoints reach radius .70; core corners pull inward to (+/-.55,+/-.55).
    // This enlarges central cells and reduces the gap to the circular boundary.
    // Even side intervals put the shared diameter on mesh edges and give exact
    // reflection and quarter-turn permutations for the three junction arms.
    const double pi=std::acos(-1.0);
    const int n = 2*static_cast<int>(std::ceil(pi/(4.0*target_size)));
    const int rings=static_cast<int>(std::ceil(0.30/target_size));
    const auto node = [n](int i, int j) { return j*(n+1)+i; };
    CrossSectionTemplates t;
    t.generated = true;
    t.target_size = target_size;
    t.core_divisions = n;
    t.outer_layers = rings;
    std::vector<int> reflected;
    std::vector<TemplateFace> faces;
    for (int j=0; j<=n; ++j) for (int i=0; i<=n; ++i) {
        const double u=-1.0+2.0*double(i)/n;
        const double v=-1.0+2.0*double(j)/n;
        t.circle.push_back({u*(CoreSideRadius-CoreBow*v*v),v*(CoreSideRadius-CoreBow*u*u),0.0});
        reflected.push_back(node(i,n-j));
    }
    for (int j=0;j<n;++j) for(int i=0;i<n;++i)
        faces.push_back({node(i,j),node(i+1,j),node(i+1,j+1),node(i,j+1)});
    std::vector<int> inner;
    for(int i=0;i<n;++i) inner.push_back(node(i,0));
    for(int j=0;j<n;++j) inner.push_back(node(n,j));
    for(int i=n;i>0;--i) inner.push_back(node(i,n));
    for(int j=n;j>0;--j) inner.push_back(node(0,j));
    std::vector<int> previous=inner;
    for(int ring=1;ring<=rings;++ring) {
        const int offset=static_cast<int>(t.circle.size());
        const double alpha=double(ring)/rings;
        std::vector<int> current;
        for(int k=0;k<4*n;++k) {
            const double angle=-3.0*pi/4.0+k*pi/(2.0*n);
            Vec3 p=t.circle[inner[k]]*(1.0-alpha)+Vec3{std::cos(angle),std::sin(angle),0.0}*alpha;
            if(std::abs(p.x)<1e-14) p.x=0.0;
            if(std::abs(p.y)<1e-14) p.y=0.0;
            t.circle.push_back(p);
            reflected.push_back(offset+(3*n-k+4*n)%(4*n));
            current.push_back(offset+k);
        }
        for(int k=0;k<4*n;++k) {
            const int next=(k+1)%(4*n);
            faces.push_back({previous[k],current[k],current[next],previous[next]});
        }
        previous=std::move(current);
    }
    t.boundary_circle=previous;
    for(int turn=0;turn<4;++turn) {
        auto& rotation=t.quarter_turn_nodes[turn];rotation.resize(t.circle.size());
        for(int j=0;j<=n;++j) for(int i=0;i<=n;++i) {
            int x=i,y=j;
            for(int k=0;k<turn;++k) { const int old=x;x=n-y;y=old; }
            rotation[node(i,j)]=node(x,y);
        }
        const int core=(n+1)*(n+1);
        for(int ring=0;ring<rings;++ring) for(int k=0;k<4*n;++k)
            rotation[core+ring*4*n+k]=core+ring*4*n+(k+turn*n)%(4*n);
    }
    // Upper half first, then lower half: every quad lies wholly on one side.
    for(int sign:{1,-1}) for(const auto& face:faces) {
        double y=0.0;for(int k:face) y+=t.circle[k].y;
        if(sign*y>0.0) t.circle_faces.push_back(face);
    }
    const int count = static_cast<int>(t.circle.size());
    std::vector<int> third(count,-1);
    const double c = std::sqrt(3.0)/2.0;
    for (const auto& p : t.circle) t.merge.push_back({p.x,c*p.y,-std::abs(p.y)/2.0});
    for (int k=0;k<count;++k) {
        if (t.circle[k].y<0.0) continue;
        if (t.circle[k].y==0.0) third[k]=k;
        else {
            third[k]=static_cast<int>(t.merge.size());
            t.merge.push_back({t.circle[k].x,0.0,t.circle[k].y});
        }
    }
    std::vector<int> left(count), right(count);
    for (int k=0;k<count;++k) {
        left[k] = t.circle[k].y<0.0 ? third[reflected[k]] : k;
        right[k] = t.circle[k].y>0.0 ? third[k] : k;
    }
    t.branch_bottom=t.circle_faces;
    for (const auto& q : t.circle_faces) {
        TemplateFace a,b;
        for (int corner=0; corner<4; ++corner) { a[corner]=left[q[corner]]; b[corner]=right[q[corner]]; }
        t.branch_left.push_back(a); t.branch_right.push_back(b);
    }
    t.merge_faces=t.circle_faces;
    const int half=static_cast<int>(t.circle_faces.size()/2);
    for (int f=0; f<half; ++f) {
        TemplateFace q=t.circle_faces[f];
        for (auto& k:q) k=third[k];
        t.merge_faces.push_back(q);
    }
    t.branch_bottom_points.assign(t.merge.begin(),t.merge.begin()+count);
    for (int k=0; k<count; ++k) {
        t.branch_left_points.push_back(t.merge[left[k]]);
        t.branch_right_points.push_back(t.merge[right[k]]);
    }
    t.boundary_merge=t.boundary_circle;
    for (int k:t.boundary_circle) {
        if (t.circle[k].y>0.0) {
            t.left_boundary.push_back(k);
            t.right_boundary.push_back(reflected[k]);
            t.bottom_boundary.push_back(third[k]);
            t.boundary_merge.push_back(third[k]);
        }
    }
    for(int sign:{-1,1}) {
        int endpoint=-1,neighbor=-1;
        double end_x=-2.0,near_x=-2.0;
        for(int k:t.boundary_circle) {
            const double x=sign*t.circle[k].x;
            if(t.circle[k].y==0.0 && x>end_x) { endpoint=k;end_x=x; }
            if(t.circle[k].y>0.0 && x>near_x) { neighbor=k;near_x=x; }
        }
        if(endpoint<0 || neighbor<0) throw std::runtime_error("generated junction diameter is incomplete");
        t.axis_projections.push_back({endpoint,neighbor,reflected[neighbor],third[neighbor]});
    }
    return t;
}

namespace {
std::ofstream Output(const std::filesystem::path& path)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write cross-section template: "+path.string());
    out<<std::setprecision(17);
    return out;
}
void Preview(const std::vector<Vec3>& points, const std::vector<TemplateFace>& faces,
    const std::filesystem::path& path)
{
    auto out=Output(path);
    out<<"# vtk DataFile Version 2.0\nGenerated cross-section template\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "<<points.size()<<" double\n";
    for (const auto& p:points) out<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
    out<<"CELLS "<<faces.size()<<' '<<5*faces.size()<<'\n';
    for (const auto& q:faces) out<<"4 "<<q[0]<<' '<<q[1]<<' '<<q[2]<<' '<<q[3]<<'\n';
    out<<"CELL_TYPES "<<faces.size()<<'\n';
    for (std::size_t i=0;i<faces.size();++i) out<<"9\n";
    iga::FlushCheckedText(out);
}
}

void WriteCrossSectionTemplatePreviews(const CrossSectionTemplates& t, const std::filesystem::path& directory)
{
    if (!directory.empty()) std::filesystem::create_directories(directory);
    Preview(t.circle,t.circle_faces,directory/"cross_section_template.vtk");
    Preview(t.merge,t.merge_faces,directory/"merge_template.vtk");
}


} // namespace tubular
