#include <math.h>
#include <QFile>
#include <QTextStream>
#include <set>
#include <map>
#include <queue>
#include <unordered_set>
#include "data.h"
#include "mymath.h"
#include "color.h"
#include "linking.h"

#define D_TYPE float
using Vec3 = arma::Col<D_TYPE>::fixed<3>;
using Mat3 = arma::Mat<D_TYPE>::fixed<3,3>;

#define FACTOR_X   1.f 
#define FACTOR_Y   1.f
#define FACTOR_A  10.f

#define MIN_SUM_W  90.
#define MIN_SUM_W2  1.

#define EPSILON 1e-10

//=============================================================================
// Maths
//=============================================================================

inline int compute_distance_h(int h1, int h2)
{
    const int dh = abs(h1 - h2)%D_HEADING;
    return min(dh, D_HEADING - dh);
}

Vec3 project_point_segment_3d(Vec3 pt, Vec3 seg_A, Vec3 seg_B,
                              D_TYPE *t, D_TYPE *distance,
                              bool clip_t)
{
    const Vec3 AB = seg_B - seg_A;
    const Vec3 AP = pt    - seg_A;
    const D_TYPE AB_sq = arma::dot(AB, AB);
    
    *t = arma::dot(AP, AB) / AB_sq; // dans [0, 1] de A vers B
    if (clip_t)
        *t = std::max((D_TYPE)0., std::min((D_TYPE)1., *t));

    Vec3 proj = seg_A + (*t)*AB;

    if (distance != nullptr)
        *distance = arma::norm(proj - pt);
    
    return proj;
}

int compute_l1_distance_grid_coord(const Vec3 &xyz1, const Vec3 &xyz2)
{
    const arma::Col<D_TYPE>::fixed<3> un_scale({1./FACTOR_X, 1./FACTOR_Y, 1./FACTOR_A});
    return int(arma::sum(arma::abs(xyz1 - xyz2) % un_scale) + D_TYPE(0.5));
}

// intersection point between the 3D plan (O, vecN) and a segment AB
Vec3 project_point_segment_3d_perpendi(Vec3 O, Vec3 N, Vec3 A, Vec3 B, D_TYPE *t, D_TYPE *distance)//, bool tip)
{
    // TODO:
    // si t<0 ou t>1
    // on applique la partie ortho avec la distance mesuréer avec clip

    Vec3 proj({0., 0., 0.});
    D_TYPE denom = arma::dot(N, B - A);
    if (denom == D_TYPE(0.))
    {
        *distance = D_TYPE(1e20);
        return proj;
    }
    *t = arma::dot(N, O - A) / denom;
    
#if 0
    if (*t < 0. || 1. < *t)
    {
        *distance = D_TYPE(1e20);
        return proj;
    }
#endif

    proj = A + (*t)*(B - A);
    *t = std::max(D_TYPE(0.), std::min(D_TYPE(1.), *t));
    Vec3 proj_clipped = A + (*t)*(B - A);
    *distance = arma::norm(proj_clipped - O);
    //if (tip)
    //    return proj;
    //else
        return proj_clipped;
}

//=============================================================================
//
//=============================================================================

void save_Vec3(const Vec3 &p, std::ofstream &fs)
{
    fs.write(reinterpret_cast<const char*>(&p[0]), sizeof(p[0]));
    fs.write(reinterpret_cast<const char*>(&p[1]), sizeof(p[1]));
    fs.write(reinterpret_cast<const char*>(&p[2]), sizeof(p[2]));
}

Vec3 load_Vec3(std::ifstream &fs)
{
    Vec3 p;
    fs.read(reinterpret_cast<char*>(&p[0]), sizeof(p[0]));
    fs.read(reinterpret_cast<char*>(&p[1]), sizeof(p[1]));
    fs.read(reinterpret_cast<char*>(&p[2]), sizeof(p[2]));
    return p;
}

//=============================================================================
//
//=============================================================================

class Point
{
public:
    Vec3 xyz, xyz0;
    int h; // when only 1 heading, in the first steps
    bool headings[D_HEADING] = {false}; // when mixed headings, afterward
    double w;
    int line_id = -1;

    Point() {};
    Point(Vec3 xyz, int h, D_TYPE w, Vec3 xyz0={0,0,0})
    {
        this->xyz = xyz;
        this->h = h;
        this->headings[h] = true;
        this->w = w;
        this->xyz0 = xyz0;
    }
    // merge 2 points
    Point(Point *point1, Point *point2)
    {
        // pos
        this->xyz = D_TYPE(0.5)*(point1->xyz + point2->xyz);

        // heading
        this->h = -1;
        for (int ih=0; ih<D_HEADING; ih++)
            this->headings[ih] = point1->headings[ih] || point2->headings[ih];
        
        // weight
        this->w = point1->w + point2->w;
    }
    Point(std::ifstream &fs)
    {
        load(fs);
    }

    void save(std::ofstream &fs) const
    {
        save_Vec3(xyz, fs);
        save_Vec3(xyz0, fs);
        fs.write(reinterpret_cast<const char*>(&h), sizeof(h));
        for (int ih=0; ih<D_HEADING; ih++)
            fs.write(reinterpret_cast<const char*>(&headings[ih]), sizeof(headings[ih]));
        fs.write(reinterpret_cast<const char*>(&w), sizeof(w));
        fs.write(reinterpret_cast<const char*>(&line_id), sizeof(line_id));
    }
    
    void load(std::ifstream &fs)
    {
        xyz  = load_Vec3(fs);
        xyz0 = load_Vec3(fs);
        fs.read(reinterpret_cast<char*>(&h), sizeof(h));
        for (int ih=0; ih<D_HEADING; ih++)
            fs.read(reinterpret_cast<char*>(&headings[ih]), sizeof(headings[ih]));
        fs.read(reinterpret_cast<char*>(&w), sizeof(w));
        fs.read(reinterpret_cast<char*>(&line_id), sizeof(line_id));
    }
};


template <typename PointContainer>
tuple<Vec3, Vec3> compute_bb(PointContainer &points, bool xyz0=false)
{
    Vec3 min_xyz = { 1e30,  1e30,  1e30};
    Vec3 max_xyz = {-1e30, -1e30, -1e30};

    const size_t nb_points = points.size();
    for (size_t p=0; p<nb_points; p++)
    {
        if constexpr (std::is_pointer_v<typename PointContainer::value_type>) {
            // Dereference if the container holds pointers
            min_xyz = arma::min(min_xyz, xyz0 ? points[p]->xyz0 : points[p]->xyz);
            max_xyz = arma::max(max_xyz, xyz0 ? points[p]->xyz0 : points[p]->xyz);
        } else {
            // Use the object directly
            min_xyz = arma::min(min_xyz, xyz0 ? points[p].xyz0 : points[p].xyz);
            max_xyz = arma::max(max_xyz, xyz0 ? points[p].xyz0 : points[p].xyz);
        }
    }

    return {min_xyz, max_xyz};
}


double distance_point_bb(const Point &point, const tuple<Vec3,Vec3> &bb)
{
    const Vec3 min_xyz = std::get<0>(bb);
    const Vec3 max_xyz = std::get<1>(bb);

    Vec3 diff;

    for (int d=0; d<3; d++)
    {
        if (min_xyz[d] <= point.xyz[d] && point.xyz[d] <= max_xyz[d])
            diff[d] = 0.;
        else if (point.xyz[d] <= min_xyz[d])
            diff[d] = min_xyz[d] - point.xyz[d];
        else
            diff[d] = point.xyz[d] - max_xyz[d];
    }

    return arma::norm(diff);
}

D_TYPE distance_bb_bb(const tuple<Vec3,Vec3> &bb1, const tuple<Vec3,Vec3> &bb2)
{
    // Extract min and max points of the bounding boxes
    const Vec3 &min1 = std::get<0>(bb1);
    const Vec3 &max1 = std::get<1>(bb1);
    const Vec3 &min2 = std::get<0>(bb2);
    const Vec3 &max2 = std::get<1>(bb2);

    // Compute the distances along each axis
    D_TYPE dx = std::max(D_TYPE(0.0), std::max(min2[0] - max1[0], min1[0] - max2[0]));
    D_TYPE dy = std::max(D_TYPE(0.0), std::max(min2[1] - max1[1], min1[1] - max2[1]));
    D_TYPE dz = std::max(D_TYPE(0.0), std::max(min2[2] - max1[2], min1[2] - max2[2]));

    // Return the Euclidean distance
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

int distance_bb_bb_l1_grid(const tuple<Vec3,Vec3> &bb1, const tuple<Vec3,Vec3> &bb2)
{
    // Extract min and max points of the bounding boxes
    const Vec3 &min1 = std::get<0>(bb1);
    const Vec3 &max1 = std::get<1>(bb1);
    const Vec3 &min2 = std::get<0>(bb2);
    const Vec3 &max2 = std::get<1>(bb2);

    // Compute the distances along each axis
    D_TYPE dx = std::max(D_TYPE(0.0), std::max(min2[0] - max1[0], min1[0] - max2[0])) / FACTOR_X;
    D_TYPE dy = std::max(D_TYPE(0.0), std::max(min2[1] - max1[1], min1[1] - max2[1])) / FACTOR_Y;
    D_TYPE dz = std::max(D_TYPE(0.0), std::max(min2[2] - max1[2], min1[2] - max2[2])) / FACTOR_A;

    // Return the Euclidean distance
    return int(abs(dx) + abs(dy) + abs(dz) + 0.5);
}



class LinePoint: public Point
{
public:
    Vec3 tangent;
    vector<Point*> neighbors;

    LinePoint() {}
    LinePoint(const Point &point) : Point(point)
    {
    }
    LinePoint(const Point &point,
              const vector<Point*> &neighbors,
              const Vec3 &tangent) : Point(point)
    {
        this->neighbors = neighbors;
        this->tangent = tangent;
    }
    LinePoint(std::ifstream &fs)
    {
        load(fs);
    }
    LinePoint(LinePoint *point1, LinePoint *point2) : Point(point1, point2)
    {
        // merge neighbors
        set_union(point1->neighbors.begin(), point1->neighbors.end(),
                  point2->neighbors.begin(), point2->neighbors.end(), inserter(this->neighbors, this->neighbors.begin()));
    }

    void save(std::ofstream &fs) const
    {
        Point::save(fs);
        save_Vec3(tangent, fs);
        
        // neighbors
        const size_t nb_neighbors = neighbors.size();
        fs.write(reinterpret_cast<const char*>(&nb_neighbors), sizeof(nb_neighbors));
        for (size_t i=0; i<nb_neighbors; i++)
            neighbors[i]->save(fs);
    }

    void load(std::ifstream &fs)
    {
        Point::load(fs);
        tangent = load_Vec3(fs);

        // neighbors
        size_t nb_neighbors;
        fs.read(reinterpret_cast<char*>(&nb_neighbors), sizeof(nb_neighbors));
        neighbors.resize(nb_neighbors);
        for (size_t i=0; i<nb_neighbors; i++)
            neighbors[i] = new Point(fs);
    }
};

struct Line
{
    vector<LinePoint> points;
    Vec3 min_xyz;
    Vec3 max_xyz;
    
    size_t size() const
    {
        return points.size();
    }

    void update_bb()
    {
        tie(this->min_xyz, this->max_xyz) = compute_bb(points);
    }

    double distance_to_bb(const Point &point) const
    {
        return distance_point_bb(point, {this->min_xyz, this->max_xyz});
    }
};


class NearestNeighbors
{
    const int resolution[3] = {200, 200, 12};
    Vec3 min_xyz, max_xyz;
    vector<vector<vector<set<Point*>>>> grid;

public:
    inline tuple<int,int,int> get_cell(const Vec3 &xyz) const
    {
        const Vec3 normalized_xyz = (xyz - this->min_xyz)/(this->max_xyz - this->min_xyz);
        return tuple<int,int,int>(max(0, min(this->resolution[0]-1, (int)(normalized_xyz[0] * this->resolution[0]))),
                                  max(0, min(this->resolution[1]-1, (int)(normalized_xyz[1] * this->resolution[1]))),
                                  max(0, min(this->resolution[2]-1, (int)(normalized_xyz[2] * this->resolution[2]))));
    }

    NearestNeighbors(vector<Point> &points)
    {
        tie(this->min_xyz, this->max_xyz) = compute_bb(points);

        this->grid = vector<vector<vector<set<Point*>>>>(this->resolution[0], vector<vector<set<Point*>>>(
                                                         this->resolution[1], vector<set<Point*>>(
                                                         this->resolution[2])));

        const size_t nb_points = points.size();
        for (size_t p=0; p<nb_points; p++)
        {
            auto [cx, cy, cz] = get_cell(points[p].xyz);
            this->grid[cx][cy][cz].insert(&points[p]);
        }
    }

    void erase(Point *point)
    {
        const auto [cx, cy, cz] = get_cell(point->xyz);
        size_t avant = this->grid[cx][cy][cz].size();
        this->grid[cx][cy][cz].erase(point);
        size_t apres = this->grid[cx][cy][cz].size();
        assert(apres == avant - 1);
    }

    vector<Point*> get_neighbors(const Point& point, const D_TYPE distance, const int distance_h)
    {
        const auto [cx, cy, cz] = get_cell(point.xyz);
        const int dx = 1 + distance * resolution[0]/(max_xyz[0]-min_xyz[0]);
        const int dy = 1 + distance * resolution[1]/(max_xyz[1]-min_xyz[1]);
        const int dz = 1 + distance * resolution[2]/(max_xyz[2]-min_xyz[2]);

        vector<Point*> result;

        for (int x=clip(cx-dx, 0, resolution[0]); x<=clip(cx+dx, 0, resolution[0]-1); x++)
        for (int y=clip(cy-dy, 0, resolution[1]); y<=clip(cy+dy, 0, resolution[1]-1); y++)
        for (int z=clip(cz-dz, 0, resolution[2]); z<=clip(cz+dz, 0, resolution[2]-1); z++)
        {
            for (auto it=this->grid[x][y][z].begin(); it!=this->grid[x][y][z].end(); it++)
            {
                if (compute_distance_h((*it)->h, point.h) <= distance_h)
                {
                    if (arma::norm((*it)->xyz - point.xyz) <= distance)
                        result.push_back(*it);
                }
            }
        }

        return result;
    }
};


class ValueField
{
    Data *data;

public:
    ValueField(Data *data)
    {
        this->data = data;
    }

    void prepare_access(const Point& point, float p[4][4][4][4], float coords[4])
    {
        const float x = point.xyz[0] / FACTOR_X;
        const float y = point.xyz[1] / FACTOR_Y;
        const float a = point.xyz[2] / FACTOR_A;
        const float h = (float)point.h + .5f;

        data->prepare_interpolation(x, y, h, a, p, coords);
    }

    Vec3 get_tangent(const Point& point)
    {
        float p[4][4][4][4];
        float coords[4];
        prepare_access(point, p, coords);

        int derivative_order_dxdx[4] = {0, 2, 0, 0};
        int derivative_order_dxdy[4] = {1, 1, 0, 0};
        int derivative_order_dydy[4] = {2, 0, 0, 0};
        int derivative_order_dxda[4] = {0, 1, 0, 1};
        int derivative_order_dyda[4] = {1, 0, 0, 1};
        int derivative_order_dada[4] = {0, 0, 0, 2};

        float dxdx = nCubicInterpolate(4, (float*)p, coords, derivative_order_dxdx);
        float dxdy = nCubicInterpolate(4, (float*)p, coords, derivative_order_dxdy);
        float dydy = nCubicInterpolate(4, (float*)p, coords, derivative_order_dydy);
        float dxda = nCubicInterpolate(4, (float*)p, coords, derivative_order_dxda);
        float dyda = nCubicInterpolate(4, (float*)p, coords, derivative_order_dyda);
        float dada = nCubicInterpolate(4, (float*)p, coords, derivative_order_dada);

        Mat3 hessian = {{dxdx, dxdy, dxda},
                        {dxdy, dydy, dyda},
                        {dxda, dyda, dada}};
        Vec3 eigval;
        Mat3 eigvec;
        arma::eig_sym(eigval, eigvec, hessian); // The eigenvalues are in ascending order

        return arma::normalise(eigvec.col(2) % arma::Col<D_TYPE>::fixed<3>({FACTOR_X, FACTOR_Y, FACTOR_A}));
    }

    Vec3 get_gradient(const Point& point)
    {
        float p[4][4][4][4];
        float coords[4];
        prepare_access(point, p, coords);

        int derivative_order_dx[4] = {0, 1, 0, 0};
        int derivative_order_dy[4] = {1, 0, 0, 0};
        int derivative_order_da[4] = {0, 0, 0, 1};

        float dx = nCubicInterpolate(4, (float*)p, coords, derivative_order_dx);
        float dy = nCubicInterpolate(4, (float*)p, coords, derivative_order_dy);
        float da = nCubicInterpolate(4, (float*)p, coords, derivative_order_da);

        return Vec3({dx, dy, da});
    }

    D_TYPE get_value(const Point& point)
    {
        float p[4][4][4][4];
        float coords[4];
        prepare_access(point, p, coords);

        int derivative_order_val[4] = {0, 0, 0, 0};

        return nCubicInterpolate(4, (float*)p, coords, derivative_order_val);
    }
};




class Graph
{
public:
    
    Graph(vector<Line> &lines)
    {
        for (size_t l=0; l<lines.size(); l++)
        {
            for (size_t p=0; p<lines[l].size()-1; p++)
            {
                add_edge(&lines[l].points[p], &lines[l].points[p+1]);
            }
        }
    }

    Graph(string path)
    {
        std::ifstream fs(path, std::ios::binary);
        if (!fs)
        {
            std::cerr << "Failed to open file: " << path << std::endl;
            return;
        }

        // Step 1: Deserialize all LinePoints
        size_t num_points;
        fs.read(reinterpret_cast<char*>(&num_points), sizeof(num_points));

        std::vector<LinePoint*> id_to_point(num_points);
        for (size_t i = 0; i < num_points; ++i)
            id_to_point[i] = new LinePoint(fs);

        // Step 2: Deserialize adjacency list
        size_t num_edges;
        fs.read(reinterpret_cast<char*>(&num_edges), sizeof(num_edges));

        for (size_t i=0; i<num_edges; ++i)
        {
            size_t point_id;
            fs.read(reinterpret_cast<char*>(&point_id), sizeof(point_id));

            size_t num_neighbors;
            fs.read(reinterpret_cast<char*>(&num_neighbors), sizeof(num_neighbors));

            LinePoint* point = id_to_point[point_id];
            std::set<LinePoint*> neighbors;

            for (size_t j = 0; j < num_neighbors; ++j)
            {
                size_t neighbor_id;
                fs.read(reinterpret_cast<char*>(&neighbor_id), sizeof(neighbor_id));
                neighbors.insert(id_to_point[neighbor_id]);
            }

            adjList[point] = neighbors;
        }

        fs.close();
    }

    void save(const std::string &path) const
    {
        std::ofstream fs(path, std::ios::binary);
        if (!fs)
        {
            std::cerr << "Failed to open file: " << path << std::endl;
            return;
        }

        // Step 1: Serialize all unique LinePoints
        std::unordered_map<LinePoint*, size_t> point_to_id;
        std::vector<LinePoint*> id_to_point;

        size_t current_id = 0;
        for (const auto &[point, neighbors] : adjList)
        {
            if (point_to_id.find(point) == point_to_id.end())
            {
                point_to_id[point] = current_id++;
                id_to_point.push_back(point);
            }
            for (const auto &neighbor : neighbors)
            {
                if (point_to_id.find(neighbor) == point_to_id.end())
                {
                    point_to_id[neighbor] = current_id++;
                    id_to_point.push_back(neighbor);
                }
            }
        }

        // Write the number of points
        size_t num_points = id_to_point.size();
        fs.write(reinterpret_cast<const char*>(&num_points), sizeof(num_points));

        // Write each LinePoint
        for (const auto &point : id_to_point)
            point->save(fs);

        // Step 2: Serialize adjacency list
        size_t num_edges = adjList.size();
        fs.write(reinterpret_cast<const char*>(&num_edges), sizeof(num_edges));

        for (const auto &[point, neighbors] : adjList)
        {
            size_t point_id = point_to_id[point];
            fs.write(reinterpret_cast<const char*>(&point_id), sizeof(point_id));

            size_t num_neighbors = neighbors.size();
            fs.write(reinterpret_cast<const char*>(&num_neighbors), sizeof(num_neighbors));

            for (const auto &neighbor : neighbors)
            {
                size_t neighbor_id = point_to_id[neighbor];
                fs.write(reinterpret_cast<const char*>(&neighbor_id), sizeof(neighbor_id));
            }
        }

        fs.close();
    }

    void add_edge(LinePoint* src, LinePoint* dst)
    {
        adjList[src].insert(dst);
        adjList[dst].insert(src);
    }

    void erase_edge(LinePoint* src, LinePoint* dst)
    {
        adjList[src].erase(dst);
        adjList[dst].erase(src);
    }

    void collapse_edge(LinePoint* A, LinePoint* B, LinePoint* new_point)
    {
        adjList[A].insert(adjList[B].begin(), adjList[B].end());
        adjList[new_point] = adjList[A];
        adjList[new_point].erase(A);
        adjList[new_point].erase(B);
        adjList.erase(A);
        adjList.erase(B);

        // forall neighbor of A or B
        for (const auto& neighbor : adjList[new_point])
        {
            adjList[neighbor].erase(A);
            adjList[neighbor].erase(B);
            adjList[neighbor].insert(new_point);
        }
    }

    void split_edge(LinePoint* A, LinePoint* B, LinePoint* new_point)
    {
        add_edge(A, new_point);
        add_edge(B, new_point);
        erase_edge(A, B);
    }

    bool are_linked(LinePoint* A, LinePoint* B, int max_depth)
    {
        if (A == B) return true;

        std::queue<std::pair<LinePoint*, int>> q;
        std::unordered_set<LinePoint*> visited;

        q.emplace(A, max_depth);
        visited.insert(A);

        while (!q.empty())
        {
            auto [current, depth] = q.front();
            q.pop();

            for (const auto& neighbor : adjList.at(current))
            {
                if (neighbor == B) return true;

                if (depth > 1 && visited.find(neighbor) == visited.end())
                {
                    q.emplace(neighbor, depth - 1);
                    visited.insert(neighbor);
                }
            }
        }

        return false;
    }

    template <typename Stream>
    static void export_point(Stream &stream, Vec3 pos)
    {
        const Vec3 drawing_center = {D_WIDTH/2.f, D_HEIGHT/2.f, 0.f};
        const Vec3 drawing_scale  = {1.f/FACTOR_X, -1.f/FACTOR_Y, 4.f/FACTOR_A};
        const Vec3 p = (pos - drawing_center) % drawing_scale;

        stream << "[" << p[0] << "," << p[1]  << "," << p[2] << "]";
    }

    void export_color(QTextStream &stream, LinePoint *p1, LinePoint *p2)
    {
        D_TYPE sum_cos=0., sum_sin=0., sum=0.;
        LinePoint *ps[2] = {p1, p2};
        for (int ip=0; ip<2; ip++)
            for (int ih=0; ih<D_HEADING; ih++)
            {
                if (ps[ip]->headings[ih])
                {
                    D_TYPE angle = (D_TYPE(ih) + .5f)/D_HEADING * 2.*M_PI;
                    sum += 1;
                    sum_cos += cos(angle);
                    sum_sin += sin(angle);
                }
            }
        
        //if (p1->line_id == p2->line_id)
        if (sum > D_TYPE(0.))
        {
            float angle = atan2(sum_sin, sum_cos)/M_PI*180.f + 180.f;
            while (angle >= 360.f) angle -= 360.f;
            RGB rgb = HSVtoRGB({angle, 1., 1.});
            stream << "[" << rgb.r << "," << rgb.g << "," << rgb.b << "]";
        }
        else
            stream << "[" << 1 << "," << 1 << "," << 1 << "]";
    }

    void export_json(string path, string additional_export)
    {
        QFile file(path.c_str());
        file.open(QIODevice::WriteOnly);
        QTextStream stream(&file);

        bool first = true;

        stream << "[";
        for (const auto& node : adjList)
        {
            for (const auto& neighbor : node.second)
            {
                if (node.first > neighbor)
                    continue;
                
                if (!first)
                    stream << ",";
                
                stream << "{";
                stream <<     "\"type\": 1,"; // line
                stream <<     "\"color\":";
                export_color(stream, node.first, neighbor);
                stream << ",";
                stream <<     "\"points\": [";
                export_point(stream, node.first->xyz);
                stream << ",";
                export_point(stream, neighbor->xyz);
                stream <<     "]";
                stream << "}";

                first = false;
            }

            // 1 edge -> it is a tip
            if (node.second.size() == 1)
            {
                if (!first)
                    stream << ",";
                
                stream << "{";
                stream <<     "\"type\": 2,"; // sphere
                stream <<     "\"color\": [1, 1, 1],";
                stream <<     "\"radius\": 0.2,";
                stream <<     "\"opacity\": 1,";
                stream <<     "\"position\": ";
                export_point(stream, node.first->xyz);
                stream << "}";
                
                first = false;
            }
        }

        stream << additional_export.c_str();
        stream << "]";

        file.close();
    }

    std::unordered_map<LinePoint*, std::set<LinePoint*>> adjList;
};



//=============================================================================
//
//=============================================================================

struct Edge
{
    LinePoint* nodes[2];
    D_TYPE collpase_cost;

    Edge(LinePoint *n1, LinePoint *n2)
    {
        this->nodes[0] = n1 < n2 ? n1 : n2;
        this->nodes[1] = n1 < n2 ? n2 : n1;
        
        // compue_cost
        this->collpase_cost = arma::norm(n1->xyz - n2->xyz);
    }

    LinePoint* get(int i) const
    {
        return this->nodes[i];
    }
};

struct EdgeComparator
{
    const Graph* graph;

    EdgeComparator(const Graph* graph) : graph(graph) {}

    bool operator()(const Edge& e1, const Edge& e2) const
    {
        // Primary sorting by collapse_cost
        if (e1.collpase_cost != e2.collpase_cost)
            return e1.collpase_cost < e2.collpase_cost;

        // Secondary sorting by nodes[0]
        if (e1.nodes[0] != e2.nodes[0])
            return e1.nodes[0] < e2.nodes[0];

        // Tertiary sorting by nodes[1]
        return e1.nodes[1] < e2.nodes[1];
    }
};



//=============================================================================
// main
//=============================================================================


void linking(Data *data, vector<Particle> *_particles, vector<Particle> *_particles0)
{
#if 0
    {
        vector<Line> lines;
        Line line;
        line.points.resize(3);
        line.points[0].xyz = {0,6,6};
        line.points[1].xyz = {1,6,6};
        line.points[2].xyz = {2,6,6};
        line.points[0].neighbors.resize(2);
        Vec3 xyz = {1,2,3};
        line.points[0].neighbors[0] = new Point(xyz, 0, 1.);
        line.points[0].neighbors[1] = new Point(xyz, 0, 1.);
        lines.push_back(line);

        Graph graph(lines);
        cout << "graph.adjList.size() " << graph.adjList.size() << endl;
        graph.save("/tmp/graph.bin");
        
        cout << "LOAD " << endl;
        Graph graph2("/tmp/graph.bin");
        cout << "graph2.adjList.size() " << graph2.adjList.size() << endl;
        graph2.export_json("res/loaded.json", "");
    }
    return;
#endif


    assert(_particles0->size() == _particles->size());


    const size_t nb_particles = _particles->size();
    ValueField value_field(data);

    //=========================================================================
    // Scale points
    //=========================================================================

    vector<Point> points(nb_particles);

    for (size_t p=0; p<nb_particles; p++)
    {
        points[p] = Point(// xyz
                          {(*_particles)[p].x * FACTOR_X,
                           (*_particles)[p].y * FACTOR_Y,
                           (*_particles)[p].a * FACTOR_A},
                          // h
                          int((*_particles)[p].h),
                          // w
                          (*_particles)[p].w,
                          // xyz 0
                          {(*_particles0)[p].x * FACTOR_X,
                           (*_particles0)[p].y * FACTOR_Y,
                           (*_particles0)[p].a * FACTOR_A});
    }


    //=========================================================================
    // Construct NearestNeighbors
    //=========================================================================

    NearestNeighbors nn(points);

    //=========================================================================
    // Compute lines
    //=========================================================================

    vector<Line> lines;

    const float distance = 2.f;
    const float sigma = distance / 1.7f;
    const int step_size_multiplication = 3; // augmenter pour diminuer step_size
    const float step_size = .25f/float(step_size_multiplication) * distance;

    for (size_t p=0; p<nb_particles; p+=10)
    {
        if (points[p].line_id >= 0)
            continue;
        
        Line line;
        set<Point*> line_particles;
        Point line_pos_start = points[p];
        D_TYPE heading = ((D_TYPE)line_pos_start.h + 0.5)/D_HEADING * 2.*M_PI;
        const D_TYPE cos_h = cos(heading);
        const D_TYPE sin_h = sin(heading);

        for (int way=0; way<2; way++)
        {
            int nb_empty_steps = 0;
            Point line_pos = line_pos_start;

            for (int ip=0; true; ip++)
            {
                vector<Point*> neighbors = nn.get_neighbors(line_pos, distance, 0);

                if (neighbors.size() > 2)
                {
                    Vec3 sum_pos(arma::fill::zeros);
                    float sum_w = 0.f;
                    bool all_in = true;

                    for (auto it=neighbors.begin(); it!=neighbors.end(); it++)
                    {
                        if ((*it)->line_id < 0)
                        {
                            Vec3 diff = (*it)->xyz - line_pos.xyz;
                            const float w = (*it)->w * exp(-0.5*arma::dot(diff, diff)/(sigma*sigma));

                            sum_w   += w;
                            sum_pos += w * (*it)->xyz;

                            if (line_particles.find(*it) == line_particles.end())
                            {
                                line_particles.insert(*it);
                                all_in = false;
                            }
                        }
                    }

                    if (ip==0 && sum_w < MIN_SUM_W)
                        break;

                    if ((all_in && !(ip==0 && way==1)) || sum_w < MIN_SUM_W2)
                    {
                        if (++nb_empty_steps >= step_size_multiplication*2)
                            break;
                    }
                    else
                        nb_empty_steps = 0;

                    Vec3 tangent = value_field.get_tangent(line_pos);

                    // projection de line_pos sur mean perpendiculairement à tangent
                    if (sum_w >= MIN_SUM_W2)
                    {
                        Vec3 mean_pos = sum_pos / sum_w;
                        line_pos.xyz += mean_pos-line_pos.xyz - arma::dot(mean_pos-line_pos.xyz, tangent) * tangent;
                    }

                    // set sign of tangent according to heading
                    tangent = (cos_h*tangent[0] + sin_h*tangent[1] > 0.f ? 1.f : -1.f) * (float)(way*2-1) * tangent;

                    Point new_pos(line_pos.xyz + step_size * tangent, line_pos.h, sum_w);

                    // next pos
                    if (way==0)
                    {
                        if (ip==0) {
                            line_pos.w = sum_w;
                            line.points.push_back(LinePoint(line_pos, neighbors, tangent));
                        }
                        line.points.push_back(LinePoint(new_pos, neighbors, tangent));
                    }
                    else
                    {
                        line.points.insert(line.points.begin(), LinePoint(new_pos, neighbors, tangent));
                    }

                    line_pos = new_pos;
                }
                else
                    break;
            }

            if (nb_empty_steps > 0)
            {
                if (way == 0)
                    line.points.erase(line.points.end() - (nb_empty_steps-1), line.points.end());
                else
                    line.points.erase(line.points.begin(), line.points.begin() + (nb_empty_steps-1));
            }
        }
        
        if (line.size() > 1)
        {
            int line_id = (int)lines.size();

            for (auto itp=line.points.begin(); itp!=line.points.end(); itp++)
                itp->line_id = line_id;

            line.update_bb();
            lines.push_back(line);
            
            cout << line_id << " " << line.size() << endl;

            for (auto itp=line_particles.begin(); itp!=line_particles.end(); itp++)
                (*itp)->line_id = line_id;
        }
    }
    
    size_t nb_lines = lines.size();

    //=========================================================================
    // Smooth lines
    //=========================================================================

    const D_TYPE sigma2 = 1.25;
    const D_TYPE max_distance2 = 3. * sigma2;

#if 1 // gradient descent

    for (int it=0; it<100 / 1; it++)
    {
        cout << it << endl;
        vector<Line> lines_next = lines;

        for (size_t l1=0; l1<nb_lines; l1++)
        {
            const size_t l1_size = lines[l1].size();

#pragma omp parallel for
            for (size_t p1=0; p1<lines[l1].size(); p1++)
            {
                //=============================================================
                // attraction cost
                //=============================================================

                Vec3 attraction = {0., 0., 0.};
                const Vec3 line_tan = lines_next[l1].points[p1 < l1_size-1 ? p1+1 : p1].xyz -
                                      lines_next[l1].points[p1 > 0         ? p1-1 :  0].xyz;
                const D_TYPE tan_norm = arma::norm(line_tan);

                if (tan_norm > D_TYPE(0.))
                {
                    for (size_t l2=0; l2<nb_lines; l2++)
                    {
                        // l2 candidate for p1 projection ?
                        if (l1 == l2 ||
                            compute_distance_h(lines[l1].points[0].h, lines[l2].points[0].h) > 1 ||
                            lines[l2].distance_to_bb(lines[l1].points[p1]) >= max_distance2)
                            continue;

                        // project p1 -> l2
                        D_TYPE closest_distance = -1.;
                        Vec3 closest_proj = {-1., -1., -1.};
                        const size_t l2_size = lines[l2].size();

                        for (size_t p2=0; p2<l2_size-1; p2++)
                        {
                            D_TYPE distance, t;
    #if 0
                            Vec3 proj = project_point_segment_3d(lines[l1].points[p1  ].xyz,
                                                                 lines[l2].points[p2  ].xyz,
                                                                 lines[l2].points[p2+1].xyz,
                                                                 &t, &distance, true);
    #else
                            Vec3 proj = project_point_segment_3d_perpendi(lines[l1].points[p1  ].xyz,
                                                                          line_tan,
                                                                          lines[l2].points[p2  ].xyz,
                                                                          lines[l2].points[p2+1].xyz,
                                                                          &t, &distance);
    #endif
                            if (p2==0 || distance < closest_distance)
                            {
                                closest_distance = distance;
                                closest_proj = proj;
                            }
                        }
                    
                        if (closest_distance < max_distance2)
                        {
                            D_TYPE w = exp(-D_TYPE(.5)*closest_distance*closest_distance/(sigma2*sigma2));
                            attraction += w * (closest_proj - lines[l1].points[p1].xyz);
                        }
                    } // for l2
                }

                if (arma::norm(attraction) > D_TYPE(0))
                    lines_next[l1].points[p1].xyz += D_TYPE(.05) * arma::normalise(attraction);

                //=============================================================
                // field
                //=============================================================

                Vec3 field_gradient = value_field.get_gradient(lines[l1].points[p1]);
                D_TYPE fiel_value = value_field.get_value(lines[l1].points[p1]);

                lines_next[l1].points[p1].xyz += D_TYPE(.2) * field_gradient/fiel_value;

                //=============================================================
                // segment length
                //=============================================================

                // TODO: éviter retractation en mesurant longueur au départ
                // TODO: éviter retractation en mesurant longueur au départ
                // TODO: éviter retractation en mesurant longueur au départ

                if (0 < p1 && p1 < lines[l1].size()-1)
                {
                    Vec3 s0 = lines_next[l1].points[p1  ].xyz - lines_next[l1].points[p1-1].xyz;
                    Vec3 s1 = lines_next[l1].points[p1+1].xyz - lines_next[l1].points[p1  ].xyz;
                
                    lines_next[l1].points[p1].xyz += D_TYPE(.05) * (s1 - s0);
                }

            } // for p1

            lines_next[l1].update_bb();
        } // for l1

        lines = lines_next;
    }
#endif
    //=========================================================================
    // Create graph from links
    //=========================================================================

    Graph graph(lines);

    //=========================================================================
    // Add inter-line edges
    //=========================================================================
#if 1
    const D_TYPE max_distance3 = D_TYPE(.2)*max_distance2;
    map<pair<LinePoint*,LinePoint*>, set<pair<D_TYPE, LinePoint*>>> points_to_add;

    for (size_t l1=0; l1<nb_lines; l1++)
    {
        for (size_t p1=0; p1<lines[l1].size(); p1++)
        {
            for (size_t l2=0; l2<nb_lines; l2++)
            {
                //=========================================================
                // l2 candidate for p1 projection ?
                //=========================================================

                if (l1 == l2 ||
                    compute_distance_h(lines[l1].points[0].h, lines[l2].points[0].h) > 1 ||
                    lines[l2].distance_to_bb(lines[l1].points[p1]) >= max_distance3)
                    continue;

                //=========================================================
                // 
                //=========================================================

                D_TYPE closest_distance = 1e20;
                LinePoint *closest_point_A=nullptr, *closest_point_B=nullptr;
                D_TYPE closest_t=0;

                for (size_t p2=0; p2<lines[l2].size()-1; p2++)
                {
                    D_TYPE t, distance;
                    Vec3 projection = project_point_segment_3d(lines[l1].points[p1  ].xyz,
                                                               lines[l2].points[p2  ].xyz,
                                                               lines[l2].points[p2+1].xyz, &t, &distance, true);
                    
                    if (p2==0 || distance < closest_distance)
                    {
                        closest_distance = distance;
                        closest_point_A = &lines[l2].points[p2  ];
                        closest_point_B = &lines[l2].points[p2+1];
                        closest_t = t;
                    }
                }
                
                if (closest_distance < max_distance3)
                {
                    if (closest_t <= D_TYPE(EPSILON))
                        graph.add_edge(&lines[l1].points[p1], closest_point_A);
                    else if (closest_t >= D_TYPE(1. - EPSILON))
                        graph.add_edge(&lines[l1].points[p1], closest_point_B);
                    else // the new point will actually be added at the end
                        points_to_add[pair<LinePoint*,LinePoint*>(closest_point_A, closest_point_B)].insert(pair<D_TYPE,LinePoint*>(closest_t, &lines[l1].points[p1]));
                }
            }
        }
    }
    
    // add all the new points to the graph
    for (const auto& [key, value_set] : points_to_add)
    {
        LinePoint* p2A = key.first;
        LinePoint* p2B = key.second;
        LinePoint *prev = p2A;
        for (const auto& [t, p1] : value_set)
        {
            LinePoint *new_point = new LinePoint(Point((D_TYPE(1.)-t)*p2A->xyz + t*p2B->xyz,
                                                       p2A->h,
                                                       (D_TYPE(1.)-t)*p2A->w + t*p2B->w));
            graph.split_edge(prev, p2B, new_point);
            graph.add_edge(p1, new_point);
            prev = new_point;
        }
    }
#endif

    //=========================================================================
    // Edge collapse
    //=========================================================================

#if 1
    
    set<Edge, EdgeComparator> edgeQueue{EdgeComparator(&graph)};

    for (const auto& pair : graph.adjList)
    {
        for (const auto& neighbor : pair.second)
            edgeQueue.insert(Edge(pair.first, neighbor));
    }

    cout << "edge collapse " << edgeQueue.size() << endl;

    for (int it=0; edgeQueue.size() > 10000; it++)
    {
        auto it_edge = edgeQueue.begin();
        Edge edge_to_collapse = *it_edge;

        if (it % 100 == 0)
            cout << it << " " << edge_to_collapse.collpase_cost << " " << edgeQueue.size() << endl;

        assert(edge_to_collapse.get(0) != edge_to_collapse.get(1));

        // Retire l'edge et ses connections de la queue

        edgeQueue.erase(it_edge);
        for (const auto& neighbor : graph.adjList[edge_to_collapse.get(0)])
            edgeQueue.erase(Edge(edge_to_collapse.get(0), neighbor));
        for (const auto& neighbor : graph.adjList[edge_to_collapse.get(1)])
            edgeQueue.erase(Edge(edge_to_collapse.get(1), neighbor));

        LinePoint *new_point = new LinePoint(edge_to_collapse.get(0),
                                             edge_to_collapse.get(1));

        graph.collapse_edge(edge_to_collapse.get(0), edge_to_collapse.get(1), new_point);

        // Met les edges du nouveau point
        for (const auto& neighbor : graph.adjList[new_point])
            edgeQueue.insert(Edge(new_point, neighbor));
    }

#endif

    //=========================================================================
    // add missing connexions
    //=========================================================================

    set<LinePoint*> tips;
    for (const auto& node : graph.adjList)
        if (node.second.size() == 1) // 1 edge -> it is a tip
            tips.insert(node.first);
    cout << tips.size() << " tips" << endl;

    
#if 1
    set<LinePoint*> debug_clouds;
    vector<tuple<Vec3,Vec3,Vec3,Vec3>> best_debug_stitching;

    {
        D_TYPE cos_h[D_HEADING];
        D_TYPE sin_h[D_HEADING];

        for (int ih=0; ih<D_HEADING; ih++)
        {
            D_TYPE angle = (D_TYPE(ih) + .5f)/D_HEADING * 2.*M_PI;
            cos_h[ih] = cos(angle);
            sin_h[ih] = sin(angle);
        }

        // Compute neightbors BB of all the graph points
        map<LinePoint*, tuple<Vec3,Vec3>> neighbors_bb_xyz0;
        for (const auto& node : graph.adjList)
            neighbors_bb_xyz0[node.first] = compute_bb(node.first->neighbors, true);

        // forall tip point
        for (const auto& tip1 : tips)
        {
            D_TYPE best_nb_direct_neighbors = 0.;
            LinePoint* best_point = nullptr;
            tuple<Vec3,Vec3> tip1_bb_xyz0 = neighbors_bb_xyz0[tip1];

            // forall the other points
#if 1
            for (const auto& node2 : graph.adjList)
            {
                LinePoint* point2 = node2.first;
#else
            for (const auto& point2 : tips)
            {
#endif
                if (tip1 != point2)
                {
                    const int distance_bb = distance_bb_bb_l1_grid(tip1_bb_xyz0, neighbors_bb_xyz0[point2]);
                    if (distance_bb > 1 || graph.are_linked(tip1, point2, 20))
                        continue;

                    vector<tuple<Vec3,Vec3,Vec3,Vec3>> debug_stitching;
                    D_TYPE nb_direct_neighbors = 0.;
                    for (const auto& p1 : tip1->neighbors)
                    {
                        D_TYPE best_w = 0.;
                        Point* best_other = nullptr;

                        for (const auto& p2 : point2->neighbors)
                        {
                            if (compute_l1_distance_grid_coord(p1->xyz0, p2->xyz0) <= 1 && compute_distance_h(p1->h, p2->h) <= D_HEADING/4)
                            {
                                Vec3 v = p2->xyz0 - p1->xyz0;
                                D_TYPE Vxy_n_sq = v[0]*v[0] + v[1]*v[1];
                                D_TYPE dot_prod_v = Vxy_n_sq > EPSILON ? (v[0]*cos_h[p1->h] + v[1]*sin_h[p1->h]) / sqrt(Vxy_n_sq) : D_TYPE(1.);
                                D_TYPE w = dot_prod_v * p2->w;
                                if (w > best_w)
                                {
                                    best_w = w;
                                    best_other = p2;
                                }
                            }
                        }
                        if (best_w > D_TYPE(0.))
                        {
                            nb_direct_neighbors += best_w;
                            debug_stitching.push_back(make_tuple(tip1->xyz, p1->xyz0, best_other->xyz0, point2->xyz));
                        }
                    }
                    
                    if (nb_direct_neighbors > best_nb_direct_neighbors)
                    {
                        best_nb_direct_neighbors = nb_direct_neighbors;
                        best_point = point2;
                        best_debug_stitching.insert(best_debug_stitching.end(), debug_stitching.begin(), debug_stitching.end());
                    }
                }
            }

            cout << "best_nb_direct_neighbors: " << best_nb_direct_neighbors << endl;

            if (best_nb_direct_neighbors > D_TYPE(0.))
            {
                LinePoint* new_connexion_A = new LinePoint();
                LinePoint* new_connexion_B = new LinePoint();
                *new_connexion_A = *tip1;
                *new_connexion_B = *best_point;
                graph.add_edge(new_connexion_A, new_connexion_B);
                
                debug_clouds.insert(tip1);
                debug_clouds.insert(best_point);
            }
        }
    }
#endif

    std::ostringstream stream;
#if 0
    {
        //for (const auto& tip : tips)
        for (const auto& tip : debug_clouds)
        {
            bool first = true;
            stream << ",{";
            stream <<     "\"type\": 0,";
            stream <<     "\"colors\": [["<< (float)(rand()%1001)/1000.f << ", " <<
                                             (float)(rand()%1001)/1000.f << ", " <<
                                             (float)(rand()%1001)/1000.f << "]],";
            stream <<     "\"points\": [";
            
            for (const auto& pt : tip->neighbors)
            {
                if (!first) stream << ",";
                Graph::export_point(stream, pt->xyz0);
                first = false;
            }

            stream <<     "]";
            stream << "}";
        }

        for (const auto& line : best_debug_stitching)
        {
            stream << ",{";
            stream <<     "\"type\": 1,";
            stream <<     "\"color\": [1, 1, 1],";
            stream <<     "\"points\": [";
            Graph::export_point(stream, std::get<0>(line));
            stream << ",";
            Graph::export_point(stream, std::get<1>(line));
            stream << ",";
            Graph::export_point(stream, std::get<2>(line));
            stream << ",";
            Graph::export_point(stream, std::get<3>(line));
            stream <<     "]";
            stream << "}";
        }
    }
#endif

    //=========================================================================
    // Export
    //=========================================================================

    graph.export_json("../results/last.json", stream.str());
}