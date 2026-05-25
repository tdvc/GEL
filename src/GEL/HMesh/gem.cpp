/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * Made by Thor Christiansen
 * ----------------------------------------------------------------------- */

#include "gem.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace CGLA;
using namespace HMesh;
using namespace Eigen;

 /* ----------------------------------------------------------------------- *
  * Finds all the edges of the set of faces (fs)
  * ----------------------------------------------------------------------- */
HMesh::HalfEdgeSet all_edges(const HMesh::Manifold& m, const HMesh::FaceSet& fs) {
    HMesh::HalfEdgeSet hs;
    for(auto f: fs)
        circulate_face_ccw(m, f, [&](HMesh::HalfEdgeID h) {
            hs.insert(h);
        });
    return hs;
}

 /* ----------------------------------------------------------------------- *
  * For every face in the faceset (fs), find its halfedge and its opposite halfedge
  * ----------------------------------------------------------------------- */
HMesh::HalfEdgeSet extended_patch_edges(const HMesh::Manifold &m, HMesh::FaceSet& fs) {
    HalfEdgeSet extend_patch_edges;
    for (auto h : all_edges(m, fs)) {
        extend_patch_edges.insert(h);
        extend_patch_edges.insert(m.walker(h).opp().halfedge());
    }
    return extend_patch_edges;
}

/* ----------------------------------------------------------------------- *
 * Finds all the vertices of the set of faces (fs)
 * ----------------------------------------------------------------------- */
HMesh::VertexSet all_verts(const HMesh::Manifold& m, const HMesh::FaceSet& fs) {
    HMesh::VertexSet vs;
    for(auto f: fs)
        circulate_face_ccw(m, f, [&](HMesh::VertexID v) {
            if(v != HMesh::InvalidVertexID && m.in_use(v))
              vs.insert(v);
        });
    return vs;
}

/* ----------------------------------------------------------------------- *
 * Finds all the boundary vertices of a set of faces (fs)
 * ----------------------------------------------------------------------- */
HMesh::VertexSet boundary_verts(const HMesh::Manifold& m, const HMesh::FaceSet& fs) {
    HMesh::VertexSet vs = all_verts(m, fs);

    HMesh::VertexSet vsb;
    for(auto v : vs)
    {
        int cnt = 0;
        circulate_vertex_ccw(m,v,[&](HMesh::FaceID f){if(fs.count(f)) ++cnt;});
        if(valency(m, v) > cnt) {
            vsb.insert(v);
        }
    }
    return vsb;
}

/* ----------------------------------------------------------------------- *
 * Find the set of boundary halfedges from a set of faces (fs) with disk topology
 * ----------------------------------------------------------------------- */
HMesh::HalfEdgeSet boundary_hes(const HMesh::Manifold &m, HMesh::FaceSet& fs) {
  HMesh::HalfEdgeSet hs = all_edges(m, fs);

  HMesh::HalfEdgeSet hsb;
  for(auto h : hs)
  {
      int cnt = 0;
      if(!fs.count(m.walker(h).face()) || !fs.count(m.walker(h).opp().face()))
          hsb.insert(h);
  }
  return hsb;
}


/* ----------------------------------------------------------------------- *
 * Purpose of function: To compute how many vertices the polygonal face f consists of
 * ----------------------------------------------------------------------- */
 int face_valency(HMesh::Manifold &m, HMesh::FaceID f) {
    if (!m.in_use(f)) {
        return 0;
    }
    int counter = 0;
    circulate_face_ccw(m, f, [&] (HMesh::VertexID vn) {
        counter += 1;
    });
    return counter;
}

/* ----------------------------------------------------------------------- *
 * Purpose of function: To find the boundary vertices of a set of faces (patch_faces)
 * in a counter clock-wise manner starting from the boundary vertex (ref_v)
 * ----------------------------------------------------------------------- */
std::vector<HMesh::VertexID> find_boundary_vertices(const HMesh::Manifold &m, HMesh::FaceSet patch_faces, HMesh::VertexID ref_v) {

    auto unordered_bd_vertices = boundary_verts(m, patch_faces);
    if (unordered_bd_vertices.find(ref_v) == unordered_bd_vertices.end()) {
        std::cout << "Reference vertex is not on the boundary!" << std::endl;
        assert(false);
    }

    std::vector<HMesh::VertexID> bd_vertices;
    bd_vertices.push_back(ref_v);
    
    HMesh::HalfEdgeSet bd_edges = boundary_hes(m, patch_faces);

    HMesh::HalfEdgeID next_edge;

    circulate_vertex_ccw(m, ref_v, [&] (HMesh::HalfEdgeID h) {
      if(bd_edges.find(h) != bd_edges.end()) {
        next_edge = h;
      }
    });

    while (m.walker(next_edge).vertex() != ref_v) {
        auto next_vertex = m.walker(next_edge).vertex();
        bd_vertices.push_back(next_vertex);
        circulate_vertex_ccw(m, next_vertex, [&] (HMesh::HalfEdgeID h) {
        if(bd_edges.find(h) != bd_edges.end()) {
            next_edge = h;
        }
        });
    }
    return bd_vertices;
}

/* ----------------------------------------------------------------------- *
 * Purpose of function: To find the boundary edges of a set of faces (patch_faces)
 * in a counter clock-wise manner starting from the boundary edge originating from boundary vertex (ref_v)
 * ----------------------------------------------------------------------- */
std::vector<HMesh::HalfEdgeID> find_boundary_edges_from_ref_v(const HMesh::Manifold &m, HMesh::FaceSet patch_faces, HMesh::VertexID ref_v) {
    std::vector<HMesh::HalfEdgeID> bd_edges_in_order;
    
    HMesh::HalfEdgeSet bd_edges = boundary_hes(m, patch_faces);

    HMesh::HalfEdgeID next_edge;

    circulate_vertex_ccw(m, ref_v, [&] (HMesh::HalfEdgeID h) {
      if(bd_edges.find(h) != bd_edges.end()) {
        next_edge = h;
      }
    });
    bd_edges_in_order.push_back(next_edge);

    while (m.walker(next_edge).vertex() != ref_v) {
        auto next_vertex = m.walker(next_edge).vertex();

        circulate_vertex_ccw(m, next_vertex, [&] (HMesh::HalfEdgeID h) {
            if(bd_edges.find(h) != bd_edges.end()) {
                next_edge = h;
            }
        });
        bd_edges_in_order.push_back(next_edge);
    }
    return bd_edges_in_order;
}


HMesh::FaceSet find_interior_faces(HMesh::Manifold &m, HMesh::HalfEdgeID h) {

//find which side is leaf
//do the queue bfs thing
//return the set
   HalfEdgeAttributeVector<int> touched(m.no_halfedges(), 0);
   FaceLoop l = trace_face_loop(m, touched, h);

   FaceSet curr_faces, interior_faces;

   for(auto h : l.hvec) {
        curr_faces.insert(m.walker(h).face());
   }

    HalfEdgeID h_above;
    HalfEdgeID h_below;

    h_above = m.walker(h).next().opp().next().halfedge();
    h_below = m.walker(h).prev().opp().prev().halfedge();


    //find interior face set
    FaceAttributeVector<int> face_status(m.no_faces(),0);
    FaceAttributeVector<int> face_visited(m.no_faces(),0);


    for(auto f : curr_faces) {
        face_status[f] = 1;
    }

    queue<FaceID> fq;
    FaceID leaf_face = InvalidFaceID;
    if(check_leaf(m, h, 0)) {
        leaf_face = m.walker(h_above).face();
    }
    else if(check_leaf(m,h,1)) {
        leaf_face = m.walker(h_below).face();
    }

    fq.push(leaf_face);
    face_visited[leaf_face] = 1;

    int limit = 10000000;
    int count = 0;


    while(!fq.empty()) {
        auto f = fq.front();
        fq.pop();
        interior_faces.insert(f);
        FaceSet nb_faces;
        circulate_face_ccw(m, f, [&](FaceID fn){
            nb_faces.insert(fn);
        });
        for(auto fn : nb_faces) {
          if(face_visited[fn] == 0 && face_status[fn] != 1 && interior_faces.find(fn) == interior_faces.end()) {
            fq.push(fn);
            face_visited[fn] = 1;
          }
        }
    
        count++;

    }

    return interior_faces;
}

// TC: The purpose of this function is the following: Sometimes it might happen, that the faces needed for an 
// extrusion are actually disconnected. This is probably a consequence of the fact that we do not decompose 
// a block of faceloops at a time but a single face loop at a time and that we always trace to find the next face loop
// with the maximal cylindricity. This happens for the quadmeshes/bone.obj mesh. Therefore, we need to find out which faces
// are disconnected and group them.
std::vector<HMesh::FaceSet> find_groups_of_faces(HMesh::Manifold &m, HMesh::FaceSet fs) {
    std::vector<HMesh::FaceSet> face_groups;

    std::list<FaceID> face_list;
    for (auto f : fs) {
        face_list.push_front(f);
    }

    while (!face_list.empty()) {
        // A group of faces
        FaceSet new_face_group;

        auto start_face = face_list.front();
        face_list.pop_front();

        std::queue<HMesh::FaceID> faces_to_vist;
        faces_to_vist.push(start_face);
        while (!faces_to_vist.empty()) {
            auto f = faces_to_vist.front();
            faces_to_vist.pop();

            new_face_group.insert(f);

            // TC: Explore its neighbour faces
            circulate_face_ccw(m, f, [&](HMesh::FaceID fn){
                if (std::find(face_list.begin(), face_list.end(), fn) != face_list.end()) {
                    faces_to_vist.push(fn);
                    face_list.remove(fn);
                }
            });

        }
        face_groups.push_back(new_face_group);
    }
    return face_groups;
}


// TC: Simply just check that the faces in the face_set are connected.
// Find groups of faces and select the largest group
HMesh::FaceSet check_faceset_connectivity(HMesh::Manifold &m, HMesh::FaceSet face_set) {

    auto face_groups = find_groups_of_faces(m, face_set);

    FaceSet largest_face_group; 
    for (int ii = 0; ii < face_groups.size(); ii++) {
        auto group = face_groups[ii];
        if (group.size() > largest_face_group.size()) {
            largest_face_group = group;
        }
    }

    return largest_face_group;
}


double geodesic_curvature(HMesh::Manifold& m, HMesh::Walker& w) {
    CGLA::Vec3d v = m.pos(w.vertex()) - m.pos(w.opp().vertex());
    HMesh::Walker wn = w.next().next();
    CGLA::Vec3d vn = m.pos(wn.opp().vertex())-m.pos(wn.vertex());
    return acos(std::min(1.0,std::max(-1.0,dot(v,vn)/(length(v)*length(vn)))));
};

MatrixXd patch_laplacian_matrix(HMesh::Manifold &m, HMesh::VertexSet verts) {
   int num_verts = verts.size();
   ArrayXXd A = ArrayXXd::Zero(num_verts,num_verts);
   Eigen::Triplet<double> T;
   std::vector<Eigen::Triplet<double>> triplet_list;
   std::map<HMesh::VertexID, int> vert_index;
   int interior_vert_index = 0;
   for(auto v : verts) {
     vert_index.insert(std::make_pair(v, interior_vert_index));
     interior_vert_index++;
   }

   for (auto v: verts) {
     int v_id_1 = vert_index.find(v)->second;
     A (v_id_1,v_id_1) = 1;
     Eigen::Triplet<double> Tv (v_id_1,v_id_1, 1);
     triplet_list.push_back(Tv);
     circulate_vertex_ccw(m,v, [&](VertexID vn){
       if(verts.find(vn) != verts.end()) {
         Eigen::Triplet<double> Tvn (v_id_1,vert_index.find(vn)->second,-1);
         A (v_id_1 , vert_index.find(vn)->second) = double(-1/double(valency(m,v)));
         triplet_list.push_back(Tvn);
       }
     });
   }
   //A.setFromTriplets(triplet_list.begin(), triplet_list.end());
   return A.matrix();
}

MatrixXd patch_laplacian_coeff(HMesh::Manifold &m, HMesh::VertexSet verts) {

  int num_verts = verts.size();
  VertexSet bd_verts;

  for (auto v: verts) {
    //Eigen::Triplet<double> Tv (int(v.get_index()),int(v.get_index()), 1);
    circulate_vertex_ccw(m, v, [&](HMesh::VertexID vn){
      if(verts.find(vn) == verts.end()) {
        bd_verts.insert(vn);
      }
    });
  }

  Eigen::MatrixXd bd_coords(bd_verts.size(),3);

  std::map<HMesh::VertexID, int> bd_vert_ids;
  int bd_vert_index = 0;
  for(auto v : bd_verts) {
    bd_vert_ids.insert(std::make_pair(v, bd_vert_index));
    bd_coords (bd_vert_index, 0) = m.pos(v)[0];
    bd_coords (bd_vert_index, 1) = m.pos(v)[1];
    bd_coords (bd_vert_index, 2) = m.pos(v)[2];
//    Vector3d(m.pos(v)[0],m.pos(v)[1],m.pos(v)[2]);
    bd_vert_index++;
  }

  std::map<HMesh::VertexID, int> vert_ids;
  int interior_vert_index = 0;
  for(auto v : verts) {
    vert_ids.insert(std::make_pair(v, interior_vert_index));
    interior_vert_index++;
  }


  int num_bd_verts = bd_verts.size();

  ArrayXXd A = ArrayXXd::Zero(num_verts, num_bd_verts);

  //Eigen::SparseMatrix<double> A(num_verts,num_bd_verts);
  Eigen::Triplet<double> T;
  std::vector<Eigen::Triplet<double>> triplet_list;
  std::vector<HMesh::VertexID> bd_vertices;
  for (auto v: verts) {
    //Eigen::Triplet<double> Tv (int(v.get_index()),int(v.get_index()), 1);
    int v_id_1 = vert_ids.find(v)->second;
    //triplet_list.push_back(Tv);
    circulate_vertex_ccw(m,v, [&](VertexID vn){
      if(bd_verts.find(vn) != bd_verts.end()) {
        bd_vertices.push_back(vn);
        A(v_id_1, bd_vert_ids.find(vn)->second) = double(-1/double(valency(m,v)));
        Eigen::Triplet<double> Tvn (v_id_1,bd_vert_ids.find(vn)->second,-1);
        triplet_list.push_back(Tvn);
      }
    });
  }

  MatrixXd coeffs(num_bd_verts, 3);

  coeffs = A.matrix()*bd_coords;

  return coeffs;

}

void smooth_faceset_lap_solve(HMesh::Manifold &m, HMesh::FaceSet faces) {

  HMesh::VertexSet verts;
  HMesh::VertexSet interior_vertices, bd_vertices;
  bd_vertices = boundary_verts(m, faces);
  for (auto f : faces)
      circulate_face_ccw(m , f, [&](VertexID v) {
        if(bd_vertices.find(v) == bd_vertices.end())
          interior_vertices.insert(v);
      });

  verts = interior_vertices;
  const int num_pts = verts.size();

  Eigen::MatrixXd lap_matrix = patch_laplacian_matrix(m, verts);

  Eigen::MatrixXd lap_coeffs = patch_laplacian_coeff(m, verts);

  Eigen::MatrixXd A(num_pts, num_pts);
  Eigen::MatrixXd B(num_pts, 3);

  A << lap_matrix;

  B << -1*lap_coeffs;

  Eigen::VectorXd B_x = B.col(0);
  Eigen::VectorXd B_y = B.col(1);
  Eigen::VectorXd B_z = B.col(2);


//bdcSvd(ComputeThinU | ComputeThinV)
  Eigen::VectorXd x_solve = A.matrix().fullPivHouseholderQr().solve(B_x.matrix());
  Eigen::VectorXd y_solve = A.matrix().fullPivHouseholderQr().solve(B_y.matrix());
  Eigen::VectorXd z_solve = A.matrix().fullPivHouseholderQr().solve(B_z.matrix());

  VertexAttributeVector<CGLA::Vec3d> new_pos =  m.positions_attribute_vector();
  int i = 0;
  for (auto v : verts) {
    new_pos[v][0] = x_solve(i);
    new_pos[v][1] = y_solve(i);
    new_pos[v][2] = z_solve(i);
    i++;
  }

  m.positions_attribute_vector() = new_pos;
}


// --- Robust helpers (median + MAD) ---
static double median(std::vector<double> v) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    size_t n = v.size(), mid = n / 2;
    std::nth_element(v.begin(), v.begin() + mid, v.end());
    double m = v[mid];
    if (n % 2 == 0) {
        auto it = std::max_element(v.begin(), v.begin() + mid);
        m = (m + *it) * 0.5;
    }
    return m;
}

static double mad(const std::vector<double>& x, double med) {
    std::vector<double> d;
    d.reserve(x.size());
    for (double xi : x) d.push_back(std::abs(xi - med));
    return median(std::move(d));
}

// --- Simple detector ---
struct SimpleResult {
    std::vector<int> high;   // robustly high values (includes plateaus)
    std::vector<int> peaks;  // local maxima (shape-based)
    std::vector<int> marked; // union
};

// k_high: robust "high" threshold
// min_bump: require peak to exceed its best neighbor by at least this amount (set 0 to disable)
SimpleResult detectHighAndPeaksSimpleCyclic(const std::vector<double>& x,
                                            double k_high = 2.0,
                                            double min_bump = 0.0)
{
    SimpleResult r;
    const int n = (int)x.size();
    if (n == 0) return r;

    // --- robust "high" threshold (unchanged) ---
    double med = median(x);
    double m = mad(x, med);
    double scale = 1.4826 * m;
    if (!(scale > 0)) scale = 1e-12;

    double thr = med + k_high * scale;

    for (int i = 0; i < n; ++i) {
        if (x[i] > thr) r.high.push_back(i);
    }

    // --- cyclic peaks (wrap-around neighbors) ---
    // For n==1 there are no meaningful peaks.
    if (n >= 2) {
        for (int i = 0; i < n; ++i) {
            int il = (i - 1 + n) % n; // left neighbor (wrap)
            int ir = (i + 1) % n;     // right neighbor (wrap)

            if (x[i] > x[il] && x[i] > x[ir]) {
                double bump = x[i] - std::max(x[il], x[ir]);
                if (bump >= min_bump) r.peaks.push_back(i);
            }
        }
    }

    // --- union (unchanged) ---
    r.marked = r.high;
    r.marked.insert(r.marked.end(), r.peaks.begin(), r.peaks.end());
    std::sort(r.marked.begin(), r.marked.end());
    r.marked.erase(std::unique(r.marked.begin(), r.marked.end()), r.marked.end());

    return r;
}


// Update the extrusion boundary vertex to the vertex with the highest curvature
// Just use turning angle as a measure of curvature
std::tuple<double, HMesh::VertexID, std::vector<std::pair<HMesh::VertexID, double>>> compute_new_bd_v_with_high_curvature(HMesh::Manifold m, HMesh::FaceSet faces, std::map<HMesh::VertexID, CGLA::Vec2d> vertex_uv_map, HMesh::HalfEdgeSet possible_bd_edges, double cutoff) {

    // Handle fail-case in case we get an invalid input from the LLM
    if (faces.size() == 0) {
         std::vector<std::pair<HMesh::VertexID, double>> empty_vector;
        return {-1.0, HMesh::InvalidVertexID, empty_vector};
    }

    // Remove faces from the manifold
    for (auto f : m.faces()) {
        if (faces.find(f) == faces.end()) {
            m.remove_face(f);
        }
    }
    // Split faces into triangles
    for (auto f : faces) {
        if (m.in_use(f)) {
            m.split_face_by_vertex(f);
        }
    }
    faces.clear();
    for (auto f : m.faces()) {
        faces.insert(f);
    }

    faces = check_faceset_connectivity(m, faces);

    HMesh::VertexID temp_ref_v = *boundary_verts(m, faces).begin();

    HMesh::VertexSet allowable_bd_vs;
    //cout << "The allowable_bd_vs are: ";
    for (auto h : possible_bd_edges) {
        if (m.in_use(h)) {
            allowable_bd_vs.insert(m.walker(h).vertex());
            //cout << m.walker(h).vertex() << " ";
        }
    }
    //cout << endl;

    std::vector<double> curvature_measurements;
    std::vector<HMesh::VertexID> candidate_vertices;

    //auto test_fs = check_faceset_connectivity(m, faces);

    auto patch_bd_edges = find_boundary_edges_from_ref_v(m, faces, temp_ref_v);
    HMesh::VertexID bd_v = HMesh::InvalidVertexID;
    double max_curvature = -1.0;


    for (int ii = 0; ii < patch_bd_edges.size(); ii++) {
        auto h = patch_bd_edges[ii];
        auto next_h = patch_bd_edges[(ii+1)%patch_bd_edges.size()];
        auto opp_h = m.walker(h).opp().halfedge();

        auto v_bd = m.walker(h).vertex();
        auto v_next_id = m.walker(next_h).vertex();

        // Compute the angle defect to compute the 
        std::vector<HMesh::VertexID> vns; 
        circulate_vertex_ccw(m, v_bd,[&](VertexID vn){
            vns.push_back(vn);
        });

        auto it = std::find(vns.begin(), vns.end(), v_next_id);
        if (it != vns.end()) {
            std::rotate(vns.begin(), it, vns.end());
        }

        double angle_defect = 0.0;
        for (int ii = 0; ii < vns.size()-1; ii++) {
            auto v = vns[ii];
            auto v_next = vns[ii+1];

            auto edge1_vector = normalize(m.pos(v) -  m.pos(v_bd));
            auto edge2_vector = normalize(m.pos(v_next) - m.pos(v_bd));

            angle_defect += acosf(std::min(std::max(dot(edge1_vector, edge2_vector), -1.0), 1.0));
        }
        
        double kappa = M_PI - angle_defect;

        if (allowable_bd_vs.find(m.walker(h).vertex()) != allowable_bd_vs.end()) {
            curvature_measurements.push_back(kappa);
            candidate_vertices.push_back(m.walker(h).vertex());

            if (kappa > max_curvature ) {
                max_curvature = kappa;
                bd_v = m.walker(h).vertex();
            }
        }
        
    }

    std::vector<std::pair<HMesh::VertexID, double>> peak_bd_vertices;
    // Find the outlier
    //auto angle_outliers = madOutliers(angle_measurements, cutoff);
    auto res = detectHighAndPeaksSimpleCyclic(curvature_measurements);
    for (int idx : res.marked) {
        peak_bd_vertices.push_back( std::make_pair(candidate_vertices[idx], curvature_measurements[idx]) );
    }

    // If the bd_v is empty, then we do not insert the element
    if (bd_v != HMesh::InvalidVertexID) {
        peak_bd_vertices.push_back(std::make_pair(bd_v, max_curvature));
    }

    return {max_curvature, bd_v, peak_bd_vertices};
}


HMesh::FaceID find_patch_center(HMesh::Manifold &m, HMesh::FaceSet fs) {

    int curr_dist = 0, max_dist = -1;

    HMesh::FaceID patch_center;

    for (auto f : fs) {

        curr_dist = -1;

        HMesh::FaceID left_f, right_f, top_f, bottom_f;

        HMesh::HalfEdgeID right = m.walker(f).halfedge();
        HMesh::HalfEdgeID top = m.walker(right).next().halfedge();
        HMesh::HalfEdgeID left = m.walker(top).next().halfedge();
        HMesh::HalfEdgeID bottom = m.walker(left).next().halfedge();
        int left_flag = 0, right_flag = 0, top_flag = 0, bottom_flag = 0;

        int dist_limit = 1000;
        int count = 0;

        while (left_flag == 0 && right_flag == 0 && top_flag == 0 && bottom_flag == 0 ) {

            left_f = m.walker(left).opp().face();

            right_f = m.walker(right).opp().face();

            top_f = m.walker(top).opp().face();

            bottom_f = m.walker(bottom).opp().face();

            if(fs.find(left_f) != fs.end()) {
                left = m.walker(left).opp().next().next().halfedge();
            }
            else
                left_flag = 1;

            if(fs.find(right_f) != fs.end()) {
                right = m.walker(right).opp().next().next().halfedge();
            }
            else
                right_flag = 1;

            if(fs.find(top_f) != fs.end()) {
                top = m.walker(top).opp().next().next().halfedge();
            }
            else
                top_flag = 1;

            if(fs.find(bottom_f) != fs.end()) {
                bottom = m.walker(bottom).opp().next().next().halfedge();
            }
            else
                bottom_flag = 1;

            curr_dist++;
            count++;
            if(count > dist_limit) {
              break;
            }

        }

        if(curr_dist > max_dist) {
            max_dist = curr_dist;
            patch_center = f;
        }
    }

    return patch_center;

}
