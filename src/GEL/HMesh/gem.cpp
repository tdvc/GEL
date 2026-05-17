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
   map<VertexID, int> vert_index;
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

