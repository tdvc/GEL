/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * Made by Thor Christiansen
 * ----------------------------------------------------------------------- */

#include "gem.h"
#include <GEL/HMesh/face_loop.h>
#include <GEL/HMesh/Manifold.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/segment_segment_intersect.h>
#include <GEL/HMesh/Delaunay_triangulate.h>

using namespace CGLA;
using namespace HMesh;
using namespace Eigen;

 
/* ----------------------------------------------------------------------- *
 * Function, which just determines, whether a string is actually an integer
 * ----------------------------------------------------------------------- */
bool isInteger(const std::string& str) {
    if (str.empty() || str == " ") {
        return false; // Empty strings are not integers
    }

    size_t start = 0;

    // Check for optional '+' or '-' sign
    if (str[0] == '+' || str[0] == '-') {
        start = 1;
    }

    // Traverse the rest of the string
    for (size_t i = start; i < str.length(); ++i) {
        if (str[i] < '0' || str[i] > '9') { // Check if character is a digit
            return false;
        }
    }

    return true;
}

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
  * Finds all the edges of the set of faces (fs), but does not return both halfedges
  * ----------------------------------------------------------------------- */
HMesh::HalfEdgeSet all_edges_wo_duplicates(const HMesh::Manifold& m, const HMesh::FaceSet& fs) {
    HMesh::HalfEdgeSet hs;
    for(auto f: fs)
        circulate_face_ccw(m, f, [&](HalfEdgeID h) {
            if (hs.find(m.walker(h).opp().halfedge()) == hs.end()) {
                hs.insert(h);
            }
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
 * Finds all the interior vertices on a set of faces
 * ----------------------------------------------------------------------- */
HMesh::VertexSet find_interior_vertices(const HMesh::Manifold &m, HMesh::FaceSet patch_faces) {
    HMesh::VertexSet interior_verts;

    // Get all the vertices in the patch
    HMesh::VertexSet all_patch_verts = all_verts(m, patch_faces);

    // Get the boundary vertices
    HMesh::VertexSet bd_verts = boundary_verts(m, patch_faces);

    // The interior vertices are those that are not on the boundary
    for (auto v : all_patch_verts) {
        if (bd_verts.find(v) == bd_verts.end()) {
            interior_verts.insert(v);
        }
    }

    return interior_verts;
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
std::vector<HMesh::VertexID> ccw_ordered_bd_vertices(const HMesh::Manifold &m, HMesh::FaceSet patch_faces, HMesh::VertexID ref_v) {

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
std::vector<HMesh::HalfEdgeID> ccw_ordered_bd_edges(const HMesh::Manifold &m, HMesh::FaceSet patch_faces, HMesh::VertexID ref_v) {
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


// The purpose of this function is to track boundary loops given boundary edges
std::vector<std::vector<HMesh::HalfEdgeID>> find_boundary_edge_loops(HMesh::Manifold &m, HMesh::FaceSet face_set, HMesh::HalfEdgeSet bd_edges) {
    std::vector<std::vector<HMesh::HalfEdgeID>> bd_loops;

    std::vector<HMesh::HalfEdgeID> edge_loop; 


    do {
        auto edge_begin = *bd_edges.begin();

        edge_loop.clear();
        edge_loop.push_back(edge_begin);
        
        HMesh::HalfEdgeID next_edge;
        auto incident_vertex = m.walker(edge_begin).vertex();
        circulate_vertex_ccw(m, incident_vertex, [&] (HalfEdgeID h) {
            if (bd_edges.find(h) != bd_edges.end()) {
                next_edge = h;
            }
        }); 
        
        while (next_edge != edge_begin && !bd_edges.empty()) {
 
            edge_loop.push_back(next_edge);

            auto incident_vertex = m.walker(next_edge).vertex();
            circulate_vertex_ccw(m, incident_vertex, [&] (HalfEdgeID h) {
                if (bd_edges.find(h) != bd_edges.end()) {
                    next_edge = h;
                }
            });
        }
        for (auto h : edge_loop) {
            bd_edges.erase(h);
        }


        bd_loops.push_back(edge_loop);

        if (bd_edges.empty()) {
            break;
        }
    } while(true);

    return bd_loops;
}


// Ensure that the the face_set has disk topology, so there are no holes in the middle. If there is a hole in the middle, 
// we add the faces that make up the hole in the middle to the face_set
HMesh::FaceSet check_topology(HMesh::Manifold &m, HMesh::FaceSet face_set) {

    HMesh::FaceSet new_face_set;
    for (auto f : face_set) {
        new_face_set.insert(f);
    }

    // As a first step, just add all of those faces of m, which are not in face_set, but only has neighbouring faces in face_set. These should be added anyway
    for (auto f : m.faces()) {

        // Check that it is not already in the set
        if (face_set.find(f) == face_set.end()) {
            bool should_be_included = true;
            circulate_face_ccw(m, f, [&] (FaceID fn) {
                if (face_set.find(fn) == face_set.end()) {
                    // One neighbour face was not in face_set, so we do not included, because it is not totally sourrounded
                    should_be_included = false;
                }
            }); 
            if (should_be_included) {
                new_face_set.insert(f); // We add the face
            }
        }
    }

    // As a second step: There might be a triangle in the faceset, which has a vertex on the boundary, and if that is the case, then it should be included.
    // It also causes an error, when we try to find the boundary edges from a vertex, which causes the code to crash. This problem can be solved in two different ways: 
    // 1) Find the boundary edges of the face-set and detect, if any vertex appears twice
    // 2) Find any triangle face that has neighbours in the face-set and then ensure that the vertex is on the boundary. 
    // Based on this 1) seems like the best solution
    auto bd_edges = boundary_hes(m, new_face_set);
    std::map<HMesh::VertexID, std::vector<HMesh::HalfEdgeID>> v_to_h_map;
    for (auto h : bd_edges) {
        auto v = m.walker(h).vertex();
        if (v_to_h_map.find(v) == v_to_h_map.end()) {
            std::vector<HMesh::HalfEdgeID> tmp = {h};
            v_to_h_map.insert(std::make_pair(v, tmp));
        }
        else {
            v_to_h_map.find(v)->second.push_back(h);
        }
    }
    // Loop through all those entries that have 2 edges
    HMesh::FaceSet potential_fs;
    for (auto it : v_to_h_map) {
        // We found a vertex that appears twice
        if (it.second.size() == 2) {
            //cout << "Vertex: " << it.first << " has a face that is not included " << endl;
            //cout.flush();
            // Circulate all the faces of the vertex and include the face, which has two edges whose opposite edges are in bd_edges
            circulate_vertex_ccw(m, it.first, [&] (FaceID fn) {
                // Ensure that the face is not already in the face set
                if (new_face_set.find(fn) == new_face_set.end()) {
                    int counter = 0;
                    circulate_face_ccw(m, fn, [&] (HalfEdgeID h) {
                        auto opp_h = m.walker(h).opp().halfedge();
                        if (bd_edges.find(opp_h) != bd_edges.end()) {
                            counter += 1;
                        }
                    }); 
                    if (counter == 2) {
                        //cout << "Adding face fn " << fn << " because it was on the boundary " << endl;
                        //cout.flush();
                        new_face_set.insert(fn);
                    }
                }
            }); 
        }
    }

    bd_edges = boundary_hes(m, new_face_set);

    auto bd_loops = find_boundary_edge_loops(m, new_face_set, bd_edges);

    while (bd_loops.size() != 1) {

        double min_size = std::numeric_limits<double>::infinity();
        std::vector<HMesh::HalfEdgeID> smallest_edge_loop;
        for (auto edge_loop : bd_loops) {

            if (double(edge_loop.size()) < min_size) {
                min_size = double(edge_loop.size()); 
                smallest_edge_loop = edge_loop;
            }
        }

        for (auto h : smallest_edge_loop) {
            new_face_set.insert(m.walker(h).face());
            new_face_set.insert(m.walker(h).opp().face());
        }
        bd_edges = boundary_hes(m, new_face_set);

        bd_loops = find_boundary_edge_loops(m, new_face_set, bd_edges);

    }

    return new_face_set;
}

HMesh::FaceSet extract_face_set(std::vector<std::tuple<HMesh::FaceSet, bool, bool, std::string, int>>& faces_2_be_extruded) {
    
    HMesh::FaceSet fs;
    for (auto it : faces_2_be_extruded) {
        for (auto f : std::get<0>(it)) {
            fs.insert(f);
        }
    }
    return fs;
}

void close_holes_in_fs(HMesh::Manifold &m, 
                        std::vector<std::tuple<HMesh::FaceSet, bool, bool, std::string, int>> &faces_2_be_extruded,
                        std::map<int, std::pair<int, std::tuple<HMesh::Manifold, 
                                                                                                    HMesh::FaceSet, 
                                                                                                    HMesh::FaceSet, 
                                                                                                    std::map<HMesh::VertexID,CGLA::Vec2d>,
                                                                                                    HMesh::VertexID, 
                                                                                                    std::map<HMesh::FaceID, bool>,
                                                                                                    HMesh::Manifold, 
                                                                                                    std::vector<HMesh::HalfEdgeID>,
                                                                                                    std::map<HMesh::VertexID,CGLA::Vec2d>
                                                                                                    >>>& face_map) {



    HMesh::FaceSet fs = extract_face_set(faces_2_be_extruded);
    HMesh::FaceSet new_fs = check_topology(m, fs);

    //cout << "Size of fs: " << fs.size() << endl;
    //cout << "Size of new_fs: " << new_fs.size() << endl;
    HMesh::FaceSet added_faces;
    std::set_difference(
        new_fs.begin(), new_fs.end(),
        fs.begin(), fs.end(),
        std::inserter(added_faces, added_faces.begin())
    );

    // Figure out which faces_2_be_extruded the added faces belong to
    for (auto f : added_faces) {
        //cout << "This face has been added: " << f << endl;

        for (auto &it : faces_2_be_extruded) {
            auto tuple_value = face_map.find(std::get<4>(it))->second.second;
            
            auto face_loop_faces = std::get<1>(tuple_value);
            auto base_patch_faces = std::get<2>(tuple_value);

            if (base_patch_faces.find(f) != base_patch_faces.end() && std::get<1>(it)) {
                std::get<0>(it).insert(f);
            }
            else if (face_loop_faces.find(f) != face_loop_faces.end() && !std::get<1>(it)) {
                std::get<0>(it).insert(f);
            }
            else {
                //cout << "Could not add face: " << f << endl;
            }
        }
    }
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

/* ----------------------------------------------------------------------- *
 * Function which finds the edge, which connects the two vertices v1 and v2, if such one exists
 * ----------------------------------------------------------------------- */
HMesh::HalfEdgeID find_connecting_edge(HMesh::Manifold m, HMesh::VertexID v1, HMesh::VertexID v2) {

    HMesh::HalfEdgeID connecting_edge;
    circulate_vertex_ccw(m, v1, [&](HalfEdgeID h) {
        if (m.walker(h).vertex() == v2) {
            connecting_edge = h;
        }
    });
    if (connecting_edge == HMesh::InvalidHalfEdgeID) {
        std::cout << "vertices v1: " << v1 << " and v2: " << v2 << " are not connected" << std::endl;
        assert(false);
    }
    return connecting_edge;
}

/* ----------------------------------------------------------------------- *
 * Function which finds the face, which is shared by the two vertices v1 and v2
 * ----------------------------------------------------------------------- */
HMesh::FaceID find_shared_face(HMesh::Manifold &m, HMesh::VertexID v1, HMesh::VertexID v2) {

    if (!m.in_use(v1) || !m.in_use(v2)) {
        return HMesh::InvalidFaceID;
    }

    HMesh::FaceID shared_face = HMesh::InvalidFaceID;

    HMesh::FaceSet fn;
    circulate_vertex_ccw(m, v1, [&] (HMesh::FaceID f) {
        if (m.in_use(f)) {
            fn.insert(f);
        }
    });
    circulate_vertex_ccw(m, v2, [&] (HMesh::FaceID f) {
        if (fn.find(f) != fn.end() && m.in_use(f)) {
            shared_face = f;
        }
    });
    return shared_face;
}

/* ----------------------------------------------------------------------- *
 * Check if the two vertices can be connected
 * ----------------------------------------------------------------------- */
bool vertex_vertex_connection(HMesh::Manifold m, HMesh::VertexID v1, HMesh::VertexID v2) {
    if (are_vertices_connected(m, v1, v2)) {
        return true;
    }
    else {
        auto shared_face = find_shared_face(m, v1, v2);
        if (shared_face == HMesh::InvalidFaceID) {
            return false;
        }
        else {
            return true;
        }
    }
    return false;
}

/* ----------------------------------------------------------------------- *
 * Check if a vertex can be connected to a point on an edge h
 * ----------------------------------------------------------------------- */
bool vertex_edge_connection(HMesh::Manifold m, HMesh::VertexID v1, HMesh::HalfEdgeID h) {

    HMesh::FaceSet fs;
    circulate_vertex_ccw(m, v1, [&](FaceID fn) {
        fs.insert(fn);
    });
    auto edges = extended_patch_edges(m, fs);
    if (edges.find(h) != edges.end()) {
        return true;
    }
    return false;
}

/* ----------------------------------------------------------------------- *
 * Find the area of a 2D polygon, where the vertices are given by a Nx2 matrix
 * ----------------------------------------------------------------------- */
double area_of_polygon(Eigen::MatrixXd curr_loop_V) {
    if (curr_loop_V.rows() < 4) {
        return 0.0;
    }

    double area = 0.0;
    int n = curr_loop_V.rows();
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        area += curr_loop_V(i,0) * curr_loop_V(j,1);
        area -= curr_loop_V(j,0) * curr_loop_V(i,1);
    }
    area = fabs(area) / 2.0;
    return area;
}

std::vector<double> seg_circle_ts(CGLA::Vec2d a, CGLA::Vec2d b) {

    double eps = 1e-5;

    std::vector<double> ts;

    auto d = b - a;
    auto A = dot(d,d);
    auto B = 2.0 * dot(a, d);
    auto Cq = dot(a, a) - 1.0;
    if (A < eps && A > -eps) {
        return ts;
    }
    double disc = B * B - 4.0 * A * Cq;
    if (disc < eps) {
        return ts;
    }
    auto sqrt_disc = sqrt(std::max(0.0, disc));
    double t1 = (-B - sqrt_disc) / (2.0 * A);
    double t2 = (-B + sqrt_disc) / (2.0 * A);
    if (eps < t1 && t1 < 1.0 - eps) {
        ts.push_back(t1);
    }
    if (eps < t2 && t2 < 1.0 - eps) {
        ts.push_back(t2);
    }
    std::sort(ts.begin(), ts.end());
    return ts;
}

double sector_area(CGLA::Vec2d a, CGLA::Vec2d b) {
    double eps = 1e-5;
    CGLA::Vec2d uu, vv = CGLA::Vec2d(0.0, 0.0);
    if (length(a) > eps) {
        uu = a;
        normalize(uu);
    }
    if (length(b) > eps) {
        vv = b;
        normalize(vv);
    }
    double theta = atan2(uu[0]*vv[1] - uu[1]*vv[0], dot(uu, vv));
    return 0.5 * theta;
}

/* ----------------------------------------------------------------------- *
 * The purpose of this function is to find out, whether the curve curr_loop_V is segmenting faces on the face-loop or on the base-patch
 * ----------------------------------------------------------------------- */
double area_of_curve_inside_circle(Eigen::MatrixXd curr_loop_V) {
    if (curr_loop_V.rows() < 4) {
        return 0.0;
    }

    std::vector<CGLA::Vec2d> curve;
    for (int i = 0; i < curr_loop_V.rows(); i++) {
        curve.push_back(CGLA::Vec2d(curr_loop_V(i,0), curr_loop_V(i,1)));
    }

    double inside_area = 0.0;

    for (int ii = 0; ii < curve.size(); ii++) {
        auto a = curve[ii];
        auto b = curve[(ii+1)%curve.size()];
        
        std::vector<double> ts;
        ts.push_back(0.0);
        auto intersections = seg_circle_ts(a, b);
        for (auto t : intersections) {
            ts.push_back(t);
        }
        ts.push_back(1.0);
        for (int jj = 0; jj < ts.size() - 1; jj++) {
            auto t0 = ts[jj];
            auto t1 = ts[jj+1];
            auto p0 = a + t0 * (b - a);
            auto p1 = a + t1 * (b - a);
            auto mid = p0 + 0.5 * (p1 - p0);
            if (length(mid) <= 1.0 + 1e-12) {
                inside_area += 0.5 * (p0[0]*p1[1] - p0[1]*p1[0]);
            }
            else {
                inside_area += sector_area(p0, p1);
            }
        }
    }

    return fabs(inside_area);
}

/* ----------------------------------------------------------------------- *
 * Check if a 2D-point C is on a line between 2D-point A and 2D-point B
 * ----------------------------------------------------------------------- */
bool is_point_on_the_line(CGLA::Vec2d a, CGLA::Vec2d b, CGLA::Vec2d c) {
    double eps = 1e-5;
    
    // Determinant check
    double determinant = (a[0] - b[0])*(c[1] - b[1]) - (a[1] - b[1])*(c[0] - b[0]);
    
    // Bounding box check
    bool inside_x = (std::min(b[0],c[0]) <= a[0] && std::max(b[0],c[0]) >= a[0]);
    bool inside_y = (std::min(b[1],c[1]) <= a[1] && std::max(b[1],c[1]) >= a[1]);
    
    if (fabs(determinant) < eps && inside_x && inside_y) {
        return true;
    }
    return false;
}

/* ----------------------------------------------------------------------- *
 * Does an halfedge h in the patch m intersect a line between the two 2D points A and B.
 * 2D-coordinates for the vertices of the mesh m are in v_uv_map.
 * ----------------------------------------------------------------------- */
std::pair<bool, double> do_edge_intersect_curve_t_value(HMesh::Manifold &m, HMesh::HalfEdgeID h1, CGLA::Vec2d a, CGLA::Vec2d b, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map) {
    Matrix <double, 1, 3> p;
    Matrix <double, 1, 3> r;
    Matrix <double, 1, 3> q;
    Matrix <double, 1, 3> s;

    double eps = 1e-6;

    double t, u;

    auto p_pos = v_uv_map.find(m.walker(h1).opp().vertex())->second;
    auto r_pos = v_uv_map.find(m.walker(h1).vertex())->second;

    p << p_pos[0], p_pos[1], 0.0;
    r << r_pos[0], r_pos[1], 0.0;

    auto q_pos = a;
    auto s_pos = b;

    q << q_pos[0], q_pos[1], 0.0;
    s << s_pos[0], s_pos[1], 0.0;

    bool intersection = (igl::segment_segment_intersect(p, r - p, q, s - q, t, u) && t > eps&& t < 1.0-eps);

    //return std::make_pair(intersection, Vec2d(p(0,0) + t * (r(0,0) - p(0,0)) , p(0,1) + t * (r(0,1) - p(0,1))));
    return std::make_pair(intersection, t);
}

/* ----------------------------------------------------------------------- *
 * Check if a 2D-point C is on a line between 2D-point A and 2D-point B
 * ----------------------------------------------------------------------- */
bool do_edges_intersect(HMesh::Manifold &m, HMesh::HalfEdgeID h1, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map, HMesh::VertexID v1, HMesh::VertexID v2) {

    return do_edge_intersect_curve_t_value(m, h1, v_uv_map.find(v1)->second, v_uv_map.find(v1)->second, v_uv_map).first;
}

/* ----------------------------------------------------------------------- *
 * Checks if two vertices v1 and v2 are connected in the mesh m
 * ----------------------------------------------------------------------- */
bool are_vertices_connected(HMesh::Manifold &m, HMesh::VertexID v1, HMesh::VertexID v2) {
    
    if (!m.in_use(v1) || !m.in_use(v2) || v1 == HMesh::InvalidVertexID || v2 == HMesh::InvalidVertexID) {
        return false;
    }
    else if (v1 == v2) {
        return true;
    }

    bool connected = false;
    
    circulate_vertex_ccw(m, v1,[&](VertexID vn){
        if (vn == v2) {
            connected = true;
        }
    });
    return connected;
}

/* ----------------------------------------------------------------------- *
 * Removes an edge between two faces. This also removes the two incident faces, 
 * and a new face is then inserted
 * ----------------------------------------------------------------------- */
HMesh::FaceID remove_edge_add_face(HMesh::Manifold &m, HMesh::HalfEdgeID h) {

    if (!m.in_use(h) || boundary(m,h)) {
        h = m.walker(h).opp().halfedge();
    }

    // Assert that the vertices at the end of the edge h have valency bigger than 1, otherwise we cannot merge the two faces
    assert(valency(m, m.walker(h).vertex()) > 1 && valency(m, m.walker(h).opp().vertex()) > 1);
    //assert(m.in_use(h) && m.in_use(m.walker(h).opp().halfedge()));
    assert(m.in_use(h) && !boundary(m,h));

    std::vector<CGLA::Vec3d> pts;
    auto h_orig = h;
    while (m.walker(h).next().halfedge() != h_orig) {
        pts.push_back(m.pos(m.walker(h).vertex()));
        h = m.walker(h).next().halfedge();
    }
    h_orig = m.walker(h_orig).opp().halfedge();
    h = h_orig;
    while (m.walker(h).next().halfedge() != h_orig) {
        pts.push_back(m.pos(m.walker(h).vertex()));
        h = m.walker(h).next().halfedge();
    }
    m.remove_edge(h_orig);
    auto new_face = m.add_face(pts);
    stitch_mesh(m, 1e-5);
    return new_face;
}

/* ----------------------------------------------------------------------- *
 * Computes barycentric coordinates of a 2D point p given three 2D points V0, V1 and V2
 * ----------------------------------------------------------------------- */
std::vector<double> compute_barycentric_coordinates(Matrix<double, 1, 2> V0, Matrix<double, 1, 2> V1, Matrix<double, 1, 2> V2, Matrix<double, 1, 2> p) {

    std::vector<double> coordinates;

    double A0 = (V1[0] - p[0]) * (V2[1] - p[1]) - (V1[1] - p[1]) * (V2[0] - p[0]);
    double A1 = (V2[0] - p[0]) * (V0[1] - p[1]) - (V2[1] - p[1]) * (V0[0] - p[0]);
    double A2 = (V0[0] - p[0]) * (V1[1] - p[1]) - (V0[1] - p[1]) * (V1[0] - p[0]);

    double total = A0 + A1 + A2;

    coordinates.push_back(A0 / (total));
    coordinates.push_back(A1 / (total));
    coordinates.push_back(A2 / (total));
    return coordinates;
}

/* ----------------------------------------------------------------------- *
 * Finds the barycentric coordinates as well as the vertices of the triangle face, which contains the 2D-point "point"
 * ----------------------------------------------------------------------- */
std::tuple<FaceID, VertexID, VertexID, VertexID, double, double, double> locate_point_in_face(HMesh::Manifold &m, HMesh::FaceSet base_patch_faces, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map, CGLA::Vec2d point) {
    
    HMesh::FaceID inside_face = HMesh::InvalidFaceID;
    HMesh::VertexID v1, v2, v3;
    double alpha, beta, gamma;

    double eps = 1e-5;

    for (auto f : base_patch_faces) {


        if (face_valency(m, f) != 3) {
            continue;
        }

        // Every face is a triangle
        std::vector<CGLA::Vec2d> uv_points;
        std::vector<HMesh::VertexID> vertices;
        circulate_face_ccw(m, f, [&] (VertexID vn) {
            uv_points.push_back(v_uv_map.find(vn)->second);
            vertices.push_back(vn);

        });
        //auto pA = m.pos(m.walker(h).vertex());
        auto pA = uv_points[0];
        
        //auto pB = centre(m,f);
        auto pB = uv_points[1];

        //auto pC = m.pos(m.walker(h).opp().vertex());
        auto pC = uv_points[2];

        Matrix<double, 1, 2> p;
        p << point[0], point[1];

        Matrix<double, 1, 2> Va;
        Va << pA[0], pA[1];

        Matrix<double, 1, 2> Vb;
        Vb << pB[0], pB[1];

        Matrix<double, 1, 2> Vc;
        Vc << pC[0], pC[1];

        bool inside = true;
        
        auto coordinates = compute_barycentric_coordinates(Va, Vb, Vc, p);
        for (auto c : coordinates) {
            if(c > 1 + eps || c < -eps || std::isnan(c) || std::signbit(c)) {
                inside = false;
            }
        }

        if (inside) {
            inside_face = f;
            v1 = vertices[0];
            v2 = vertices[1];
            v3 = vertices[2];

            alpha = coordinates[0];
            beta = coordinates[1];
            gamma = coordinates[2];
        }

    }
    return std::make_tuple(inside_face, v1, v2, v3, alpha, beta, gamma);
}

/* ----------------------------------------------------------------------- *
 * Checks whether a 2D point is on an edge, and if so where on the edge it is
 * ----------------------------------------------------------------------- */
std::pair<HMesh::HalfEdgeID, double> locate_point_on_edge(HMesh::Manifold&m, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map, CGLA::Vec2d p) {

    HMesh::HalfEdgeID edge;
    double s = 0.0;

    double eps = 1e-5;

    for (auto h : m.halfedges()) {

        // Make sure that we do not run over any edge that is not possible
        if (v_uv_map.find(m.walker(h).opp().vertex()) != v_uv_map.end() && v_uv_map.find(m.walker(h).vertex()) != v_uv_map.end()) {

            auto A = v_uv_map.find(m.walker(h).opp().vertex())->second;
            auto B = v_uv_map.find(m.walker(h).vertex())->second;

            auto AB = B - A;
            auto AP = p - A;

            double t = dot(AP, AB) / powf(length(AB),2.0);

            t = std::max(0.0, std::min(1.0, t));

            auto Q = A + t * AB;

            double dist = length(Q - p);

            if (dist < eps) {
                edge = h;
                s = t;
            }
        }
    }
    return {edge, s};
}

bool vertex_face_connection(HMesh::Manifold m, HMesh::FaceSet base_patch_faces, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map, HMesh::VertexID v1, CGLA::Vec2d p) {
    
    auto [f_to_split, tmp1, tmp2, tmp3, alpha, beta, gamma] = locate_point_in_face(m, base_patch_faces, v_uv_map, p);

    HMesh::FaceSet fs;
    circulate_vertex_ccw(m, v1, [&](FaceID fn) {
        fs.insert(fn);
    });
    if (fs.find(f_to_split) != fs.end()) {
        return true;
    }
    return false;
}

CGLA::Vec3d compute_new_v_pos_3d(HMesh::Manifold &m_copy, HMesh::FaceSet& base_patch_faces_copy, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map_copy, CGLA::Vec2d new_v_uv) {

    auto [f_to_split, v1, v2, v3, alpha, beta, gamma] = locate_point_in_face(m_copy, base_patch_faces_copy, v_uv_map_copy, new_v_uv);

    // Compute the new position
    CGLA::Vec3d new_pos_3d = alpha * m_copy.pos(v1) + beta * m_copy.pos(v2) + gamma * m_copy.pos(v3);

    // Split the face
    auto new_v = m_copy.split_face_by_vertex(f_to_split); 

    v_uv_map_copy.insert(std::make_pair(new_v, new_v_uv));

    base_patch_faces_copy.erase(f_to_split);
    circulate_vertex_ccw(m_copy, new_v, [&] (HMesh::FaceID fn) {
        if (fn != HMesh::InvalidFaceID) {
            base_patch_faces_copy.insert(fn);
        }
    });

    return new_pos_3d;
}

void delaunay_triangulate_each_single_face2(HMesh::Manifold &m, HMesh::Manifold &m_copy, HMesh::FaceSet& base_patch_faces, HMesh::FaceSet& base_patch_faces_copy, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map_copy, std::vector<std::tuple<int, HMesh::FaceSet, HMesh::VertexID, bool, std::vector<HMesh::VertexID>>> &curve_enclosed_faces, std::map<HMesh::VertexID, std::tuple<CGLA::Vec2d, HMesh::VertexID, HMesh::VertexID, double>>& new_vertices) {

    std::map<int, HMesh::FaceSet> non_base_patch_faces;

    for (auto &it : curve_enclosed_faces) {
        auto extrusion_id = std::get<0>(it); //.first;
        auto& yellow_faces = std::get<1>(it); //.second;

        HMesh::FaceSet faces_to_triangulate;

        for (auto f : yellow_faces) {
            if (base_patch_faces.find(f) != base_patch_faces.end()) {
                
                // Check that it is a base-patch face and then that it has more than 3 vertices
                if (face_valency(m, f) > 3) {
                    faces_to_triangulate.insert(f);
                }
            }
        }

        // Split the faces
        for (auto f : faces_to_triangulate) {

            //cout << endl;
            //cout << "Triangulating face f: " << f << endl;

            int no_vertices = face_valency(m, f);

            // So we actually know that the face is in use
            if (no_vertices > 3) {

                HMesh::FaceSet newly_created_faces;

                //cout << "no_vertices: " << no_vertices << endl;
                //cout.flush();

                std::vector<CGLA::Vec3d> mesh_points;

                std::vector<int> edge_list;

                std::map<HMesh::VertexID, int> v_id_map;
                std::map<int, HMesh::VertexID> id_v_map;
                int counter = 0; 

                HMesh::VertexSet face_vertices;

                circulate_face_ccw(m, f, [&] (VertexID vn) {
                    face_vertices.insert(vn);
                    //cout << "Face vertex: " << vn << " with uv: " << v_uv_map.find(vn)->second << endl;
                    //cout.flush();

                    v_id_map.insert({vn, counter});
                    id_v_map.insert({counter, vn});

                    mesh_points.push_back(CGLA::Vec3d(v_uv_map.find(vn)->second[0], v_uv_map.find(vn)->second[1], 0.0));

                    edge_list.push_back(counter);
                    edge_list.push_back((counter+1)%no_vertices);
                    
                    counter += 1;
                });

                // Do the Delaunay triangulation of this particular face
                CGLA::Vec3d X_axis = CGLA::Vec3d(1.0, 0.0, 0.0);
                CGLA::Vec3d Y_axis = CGLA::Vec3d(0.0, 1.0, 0.0);

                auto result = constrained_Delaunay_triangulate(mesh_points, edge_list, X_axis, Y_axis);
                auto triangle_edges = std::get<1>(result);
                auto polygon_points = std::get<2>(result);

                // If triangle_edges, which is a map of all vertices and their outgoing edges, has more elements than mesh_points, 
                // then the constrained_Delaunay_triangulation code has added Steiner points. We don't know how to handle this yet, 
                // and therefore we just place a vertex in the middle of the face 
                HMesh::FaceID f_to_split;
                if (polygon_points.size() > mesh_points.size()) {
  
                    CGLA::Vec3d new_v_3d_pos;

                    // Loop over all the new inserted points
                    HMesh::FaceSet fs_2_split;
                    fs_2_split.insert(f);

                    for (int ii = mesh_points.size(); ii < polygon_points.size(); ii++) {
                        auto new_v_uv = polygon_points[ii];

                        //cout << "The new vertex uv coordinates are: " << new_v_uv << endl;

                        HMesh::VertexID new_v;

                        // The first time that we insert a point
                        if (ii == mesh_points.size()) {
                            //cout << "Inserting the first Steiner point in the middle of the face..." << f << endl;
                            f_to_split = f;
                            //cout << "The face to split is: " << f_to_split << endl;
                            new_v = m.split_face_by_vertex(f_to_split);
                        }
                        else {
                            // Detect which face to split
                            // Only look at those faces that have been created, because all Steiner points will of course be inside one of these faces
                            auto result = locate_point_in_face(m, fs_2_split, v_uv_map, new_v_uv);
                            f_to_split = std::get<0>(result);
                            new_v = m.split_face_by_vertex(f_to_split);
                        }

                        //cout << "Before computing the new vertex position, the face to split is: " << f_to_split << endl;
                        //cout.flush();

                        m.pos(new_v) = compute_new_v_pos_3d(m_copy, base_patch_faces_copy, v_uv_map_copy, new_v_uv);

                        //cout << "The new vertex 3D coordinates are: " << m.pos(new_v) << endl;

                        new_vertices.insert({new_v, std::make_tuple(new_v_uv, HMesh::InvalidVertexID, HMesh::InvalidVertexID, -1.0)});

                        v_uv_map.insert(std::make_pair(new_v, new_v_uv));

                        // Add new faces
                        circulate_vertex_ccw(m, new_v, [&] (FaceID fn) {
                            base_patch_faces.insert(fn);
                            yellow_faces.insert(fn);
                            fs_2_split.insert(fn);
                        });

                        // Delete the old face
                        base_patch_faces.erase(f_to_split);
                        yellow_faces.erase(f_to_split);
                        fs_2_split.erase(f_to_split);

                        //cout << "Inserted new vertex: " << new_v << " with uv: " << new_v_uv << " and 3D position: " << m.pos(new_v) << endl;
                        //cout.flush();
                    }     
                }
                else {
                    // No steiner points were added and we just triangualte it.

                    newly_created_faces.insert(f);

                    for (auto v : face_vertices) {
                        auto v_id = v_id_map.find(v)->second;
                        if (triangle_edges.find(v_id) != triangle_edges.end()) {
                            for (auto vn_id : triangle_edges.find(v_id)->second) {
                                if (id_v_map.find(vn_id) != id_v_map.end()) {
                                    auto vn = id_v_map.find(vn_id)->second;

                                    //cout << "Are vertices connected: " << v << " and " << vn << "? " << are_vertices_connected(m, v, vn) << endl;

                                    // Check that they are not already connected
                                    if (!are_vertices_connected(m, v, vn)) {
                                        //cout << "Connecting vertex: " << v << " and " << vn << endl;

                                        auto shared_face = find_shared_face(m, v, vn);
                                        if (shared_face != HMesh::InvalidFaceID && m.in_use(shared_face) && newly_created_faces.find(shared_face) != newly_created_faces.end()) {
                                        //if (shared_face != HMesh::InvalidFaceID && m.in_use(shared_face)) {

                                            auto new_face = m.split_face_by_edge(shared_face, v, vn);

                                            //cout << "Adding new face: " << new_face << " by splitting face: " << shared_face << " between vertices: " << v << " and " << vn << endl;
                                            //cout.flush();

                                            // Insert the newly created face
                                            newly_created_faces.insert(new_face);

                                            // Insert into the base patch faces
                                            base_patch_faces.insert(new_face);

                                            yellow_faces.insert(new_face);

                                        }
                                    }

                                } 
                            }

                        }
                    }
                }
            }
    
        }
    }


}


