/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * Made by Thor Christiansen 
 * ----------------------------------------------------------------------- */

#include "HarmonicMap.h"
#include "gem.h"

#include <igl/harmonic.h>
#include <igl/barycentric_coordinates.h>
#include <igl/segment_segment_intersect.h>
#include <igl/per_face_normals.h>
#include <igl/map_vertices_to_circle.h>

#include "assert.h"


/* ----------------------------------------------------------------------- *
 * Given a 2D point p with coordinates (uv), the 2D positions of the vertices in the parameterized patch V,
 * and a matrix specifying which triangle faces are made up of vertices F, the purpose of this function is to compute, which
 * triangle face contains the point p, and the vertices of this triangle face (vertex_ids), and the barycentric coordinates
 * of the point p (bary_coords)
 * ----------------------------------------------------------------------- */
int find_barycentric_coords(Eigen::MatrixXd p, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::VectorXd& bary_coords, Eigen::VectorXi& vertex_ids) {

  double eps = 1e-6;

  Eigen::MatrixXd curr_bary_coords;
  for(int i = 0; i < F.rows(); i++) {
    Eigen::MatrixXd Va = V.row(F(i,0));
    Eigen::MatrixXd Vb = V.row(F(i,1));
    Eigen::MatrixXd Vc = V.row(F(i,2));

    igl::barycentric_coordinates(p, Va, Vb, Vc, curr_bary_coords);

    bool inside = true;

    for(int j = 0; j < 3; j++) {
      double c = curr_bary_coords(j);
      if(c > 1 + eps || c < -eps || std::isnan(c))
        inside = false;
    }

    if(inside) {

      Eigen::VectorXi v_ids(3);
      vertex_ids.resize(3);
      bary_coords.resize(3);
      vertex_ids(0) = F(i,0);
      vertex_ids(1) = F(i,1);
      vertex_ids(2) = F(i,2);
      bary_coords(0) = curr_bary_coords(0,0);
      bary_coords(1) = curr_bary_coords(0,1);
      bary_coords(2) = curr_bary_coords(0,2);

      return i;
    }
  }
  return -1;
}

CGLA::Vec2d HarmonicMap::find_bd_intersection(CGLA::Vec2d uv) {

  for(int i = 0; i < bnd_uv.rows(); i++) {
    //double t;

    int curr_index = i;
    int next_index = (i+1)%bnd_uv.rows();

    Eigen::Matrix <double, 1, 3> p;
    Eigen::Matrix <double, 1, 3> r;
    Eigen::Matrix <double, 1, 3> q;
    Eigen::Matrix <double, 1, 3> s;

    p << uv[0], uv[1], 0.0;
    r << 0.0, 0.0, 0.0;
    //r << boundary_face_uv_centres(i,0), boundary_face_uv_centres(i,1), 0.0;
    q << bnd_uv(i, 0), bnd_uv(i, 1), 0.0;
    s << bnd_uv(next_index, 0), bnd_uv(next_index, 1), 0;

    double t, u;

    if(igl::segment_segment_intersect(p, r - p, q, s - q, t, u)) {


      return CGLA::Vec2d(p(0,0) + t*(-p(0,0)), p(0,1) + t*(-p(0,1)));
    }

  }
  return CGLA::Vec2d(0.0);
}

int HarmonicMap::uv_to_face(CGLA::Vec2d v_uv) {

  Eigen::Matrix<double, 1, 2> p;
  p << v_uv[0], v_uv[1];
  Eigen::VectorXd bary_coords;
  Eigen::VectorXi vertex_ids;
  int face_id = find_barycentric_coords(p, V_uv, F, bary_coords, vertex_ids);

  if(face_id == -1) {
    v_uv = find_bd_intersection(v_uv);
    Eigen::Matrix<double, 1, 2> bd_p;
    bd_p << v_uv[0], v_uv[1];

    face_id = find_barycentric_coords(bd_p, V_uv, F, bary_coords, vertex_ids);
  }
  return face_id;

}

/* ----------------------------------------------------------------------- *
 * Given two 2D points p1 and p2, find the intersection between the line segment between p1 and p2 
 * and the unit circle - Used in later functions to determine the 3D position of a point in the parameterized map
 * ----------------------------------------------------------------------- */
CGLA::Vec2d intersectRayWithUnitCircle(const CGLA::Vec2d& p1, const CGLA::Vec2d& p2) {
    CGLA::Vec2d v = p2 - p1;            // direction vector

    double dx = v[0];
    double dy = v[1];

    double x1 = p1[0];
    double y1 = p1[1];

    // Quadratic coefficients
    double A = dx*dx + dy*dy;
    double B = 2.0 * (x1*dx + y1*dy);
    double C = x1*x1 + y1*y1 - 1.0;

    double disc = B*B - 4.0*A*C;
    if (disc < 0.0 || A == 0.0) {
        return CGLA::Vec2d(0.0, 0.0);  // No real intersection
    }

    double sqrtDisc = std::sqrt(disc);

    // Two roots
    double t1 = (-B - sqrtDisc) / (2.0 * A);
    double t2 = (-B + sqrtDisc) / (2.0 * A);

    // We want the smallest t > 0
    double t = std::numeric_limits<double>::infinity();
    if (t1 > 0.0 && t1 < t) {
      t = t1;
    }
    if (t2 > 0.0 && t2 < t) {
      t = t2;
    }

    if (!std::isfinite(t)) {
        return CGLA::Vec2d(0.0, 0.0); // Intersection is behind p1
    }  

    return p1 + v * t;
}

/* ----------------------------------------------------------------------- *
 * The purpose of this function is to compute the 3D position of a given a 2D point with coordinates (u,v), 
 * based on the Harmonic Map parameterization of the mesh
 * ----------------------------------------------------------------------- */
CGLA::Vec3d HarmonicMap::uv_to_vertex(CGLA::Vec2d v_uv) {

  Eigen::Matrix<double, 1, 2> p;
  p << v_uv[0], v_uv[1];
  Eigen::VectorXd bary_coords;
  Eigen::VectorXi vertex_ids;

  int face_id = find_barycentric_coords(p, V_uv, F, bary_coords, vertex_ids);

  double x, y, z; 

  if(face_id != -1) {

    int va, vb, vc;
    double alpha, beta, gamma;

    va = vertex_ids(0);
    vb = vertex_ids(1);
    vc = vertex_ids(2);

    alpha = bary_coords(0);
    beta = bary_coords(1);
    gamma = bary_coords(2);


    x = alpha * V(va, 0) + beta * V(vb, 0) + gamma * V(vc,0);
    y = alpha * V(va, 1) + beta * V(vb, 1) + gamma * V(vc,1);
    z = alpha * V(va, 2) + beta * V(vb, 2) + gamma * V(vc,2);

    return CGLA::Vec3d(x, y, z);

  }
  else {

    // Find the Edge intersection
    auto edge_id_intersection = find_bd_edge_intersection(v_uv);

    // Find the vertices of the boundary face
    auto v1 = bd_F_vertices(edge_id_intersection, 0);
    auto v2 = bd_F_vertices(edge_id_intersection, 1);
    auto v3 = bd_F_vertices(edge_id_intersection, 2);

    // Find the length of the edge in radiance
    double alpha = atan2(V_uv(v1, 1), V_uv(v1, 0));
    alpha = std::fmod(alpha + 2.0 * M_PI, 2.0 * M_PI);

    double beta;
    if (edge_id_intersection == bd_F_vertices.rows() - 1) {
        beta = 2.0 * M_PI;
    }
    else {
        beta = atan2(V_uv(v2, 1), V_uv(v2, 0));
        beta = std::fmod(beta + 2.0 * M_PI, 2.0 * M_PI);
    }

    double edge_rad_length = beta - alpha;

    // Find the intersection between the third vertex and the v_uv point
    CGLA::Vec2d intersection_p = intersectRayWithUnitCircle(CGLA::Vec2d(V_uv(v3, 0), V_uv(v3,1)), v_uv);

    // The intersection point in radians
    double rad_length_intersection_p = atan2(intersection_p[1], intersection_p[0]);
    rad_length_intersection_p = std::fmod(rad_length_intersection_p + 2.0 * M_PI, 2.0 * M_PI);

    double t = (rad_length_intersection_p - alpha) / edge_rad_length;
    
    CGLA::Vec3d v_pos_1 = CGLA::Vec3d(V(v1,0),V(v1,1),V(v1,2));
    CGLA::Vec3d v_pos_2 = CGLA::Vec3d(V(v2,0),V(v2,1),V(v2,2));

    auto bd_v_ext_disp = (1.0 - t) * v_pos_1 + t * v_pos_2;

    auto v_pos_3 = CGLA::Vec3d(V(v3,0),V(v3,1),V(v3,2));

    // Compute the proportional distance from v3 to v_uv
    double s = length(v_uv - CGLA::Vec2d(V_uv(v3, 0), V_uv(v3, 1))) / length(intersection_p - CGLA::Vec2d(V_uv(v3, 0), V_uv(v3, 1)));

    CGLA::Vec3d v_pos = (1.0 - s) * v_pos_3 + s * bd_v_ext_disp;

    return v_pos;
  }
  return CGLA::Vec3d(x, y, z);
}

/* ----------------------------------------------------------------------- *
 * The purpose of this function is to compute a basis centered at the 3D location of a vertex on the boundary of the patch.
 * The z-axis is the normalized patch normal, the y-axis is given as the vector from the boundary point to the 
 * ----------------------------------------------------------------------- */
CGLA::Mat3x3d HarmonicMap::get_patch_frame(CGLA::Vec3d bd_pt) {

  CGLA::Vec3d y_vec = normalize(patch_centre - bd_pt);
  CGLA::Vec3d z_vec = normalize(patch_normal);
  y_vec -= dot(z_vec,y_vec)*z_vec;
  y_vec = normalize(y_vec);
  CGLA::Vec3d x_vec = normalize(cross(y_vec, z_vec));
  CGLA::Mat3x3d patch_frame = CGLA::Mat3x3d(x_vec, y_vec, z_vec);
  return patch_frame;
}

/* ----------------------------------------------------------------------- *
 * The purpose of this function is to find the index of the boundary edge of the mesh,
 * which intersects the line segment between the origin of the parameterized mesh
 * and the 2D point with coordinates (u,v)
 * ----------------------------------------------------------------------- */
int HarmonicMap::find_bd_edge_intersection(CGLA::Vec2d uv) {

  for(int i = 0; i < bnd_uv.rows(); i++) {

    int curr_index = i;
    int next_index = (i+1)%bnd_uv.rows();

    Eigen::Matrix <double, 1, 3> p;
    Eigen::Matrix <double, 1, 3> r;
    Eigen::Matrix <double, 1, 3> q;
    Eigen::Matrix <double, 1, 3> s;

    p << uv[0], uv[1], 0.0;
    r << 0.0, 0.0, 0.0;

    q << bnd_uv(i, 0), bnd_uv(i, 1), 0.0;
    s << bnd_uv(next_index, 0), bnd_uv(next_index, 1), 0;

    double t, u;

    if(igl::segment_segment_intersect(p, r - p, q, s - q, t, u)) {

      return i;
    }

  }
  return -1;
}

/* ----------------------------------------------------------------------- *
 * The purpose of this function is to compute the 3D displacement vector of a 2D point
 * with coordinates (u,v) in the paramterized patch
 * ----------------------------------------------------------------------- */
CGLA::Vec3d HarmonicMap::interp_disp_vec_extended(CGLA::Vec2d v_uv, CGLA::Vec3d source_pt) {
  // The vector which we return
  CGLA::Vec3d v_ext_disp;

  int va, vb, vc;
  double u, v, w;

  // Find the barycentric coordinates and the vertices of the triangle
  Eigen::Matrix<double, 1, 2> p;
  p << v_uv[0], v_uv[1];

  Eigen::VectorXd bary_coords;
  Eigen::VectorXi vertex_ids;
  int face_id = find_barycentric_coords(p, V_uv, F, bary_coords, vertex_ids);

  // The v_uv point is inside a triangle face on the source patch
  if (face_id != -1) {

    va = vertex_ids(0);
    vb = vertex_ids(1);
    vc = vertex_ids(2);

    u = bary_coords(0);
    v = bary_coords(1);
    w = bary_coords(2);

    CGLA::Vec3d v_ext_a = CGLA::Vec3d(V_ext(va,0),V_ext(va,1),V_ext(va,2))  - source_pt;
    CGLA::Vec3d v_ext_b = CGLA::Vec3d(V_ext(vb,0),V_ext(vb,1),V_ext(vb,2))  - source_pt;
    CGLA::Vec3d v_ext_c = CGLA::Vec3d(V_ext(vc,0),V_ext(vc,1),V_ext(vc,2))  - source_pt;

    v_ext_disp = u*v_ext_a + v*v_ext_b + w*v_ext_c;

  }
  // The point is not inside a triangle face on the source patch, so we need to find the intersection between the line segment from (0,0) to (u,v), 
  // and an edge on the boundary of the mesh. We then interpolate the displacement vectors.
  else {

    auto edge_id_intersection = find_bd_edge_intersection(v_uv);

    // Find the vertices of the boundary face
    auto v1 = bd_F_vertices(edge_id_intersection, 0);
    auto v2 = bd_F_vertices(edge_id_intersection, 1);
    auto v3 = bd_F_vertices(edge_id_intersection, 2);

    // Find the length of the edge in radiance
    double alpha = atan2(V_uv(v1, 1), V_uv(v1, 0));
    alpha = std::fmod(alpha + 2.0 * M_PI, 2.0 * M_PI);

    double beta;
    if (edge_id_intersection == bd_F_vertices.rows() - 1) {
        beta = 2.0 * M_PI;
    }
    else {
        beta = atan2(V_uv(v2, 1), V_uv(v2, 0));
        beta = std::fmod(beta + 2.0 * M_PI, 2.0 * M_PI);
    }

    double edge_rad_length = beta - alpha;

    // Find the intersection between the third vertex and the v_uv point
    CGLA::Vec2d intersection_p = intersectRayWithUnitCircle(CGLA::Vec2d(V_uv(v3, 0), V_uv(v3,1)), v_uv);

    // The intersection point in radians
    double rad_length_intersection_p = atan2(intersection_p[1], intersection_p[0]);
    rad_length_intersection_p = std::fmod(rad_length_intersection_p + 2.0 * M_PI, 2.0 * M_PI);

    double t = (rad_length_intersection_p - alpha) / edge_rad_length;
    
    CGLA::Vec3d v_ext_1 = CGLA::Vec3d(V_ext(v1,0),V_ext(v1,1),V_ext(v1,2))  - source_pt;
    CGLA::Vec3d v_ext_2 = CGLA::Vec3d(V_ext(v2,0),V_ext(v2,1),V_ext(v2,2))  - source_pt;

    auto bd_v_ext_disp = (1.0 - t) * v_ext_1 + t * v_ext_2;

    auto v3_ext_disp = CGLA::Vec3d(V_ext(v3,0),V_ext(v3,1),V_ext(v3,2))  - source_pt;

    // Compute the proportional distance from v3 to v_uv
    double s = length(v_uv - CGLA::Vec2d(V_uv(v3, 0), V_uv(v3, 1))) / length(intersection_p - CGLA::Vec2d(V_uv(v3, 0), V_uv(v3, 1)));

    v_ext_disp = (1.0 - s) * v3_ext_disp + s * bd_v_ext_disp;
  }
  return v_ext_disp;
}

/* ----------------------------------------------------------------------- *
 * The purpose of this function is to find the 2D coordinates (u,v) of a vertex with vertex-id (v)
 * in the paramterized patch. Since the vertex was there, when the paramterization was computed,
 * it should always have corresponding 2D (u,v) coordinates.
 * ----------------------------------------------------------------------- */
CGLA::Vec2d HarmonicMap::patch_vertex_uv(HMesh::VertexID v) {
  int igl_v = -1;
  if(vertex_igl_map.find(v) != vertex_igl_map.end()) {
    igl_v = vertex_igl_map.find(v)->second;
  }

  if(igl_v != -1) {
    return CGLA::Vec2d(V_uv(igl_v, 0), V_uv(igl_v, 1));
  }
  else {
    std::cout << "Issue. Not possible to find vertex v" << std::endl;
    return CGLA::Vec2d(100, 100);
  }
}

/* ----------------------------------------------------------------------- *
 * Compute the Harmonic Map of a mesh (m) given a set of faces (patch_faces) with disk topology,
 * and vertex with vertex-id (ref_v) for determining the orientation of the harmonic map
 * The string extrusion_name is just included for debugging purposes.
 * ----------------------------------------------------------------------- */
void HarmonicMap::init_HarmonicMap(HMesh::Manifold m, HMesh::FaceSet patch_faces, HMesh::VertexID ref_v, std::string extrusion_name) {

  bd_v = ref_v;

  // -------------------------
  // Prior statistics
  // -------------------------
  patch_area = 0.0;

  for (auto f : patch_faces) {
    patch_area += area(m, f);
  }

  bd_perim = 0.0;
  for(auto h : boundary_hes(m, patch_faces)) { 
    bd_perim += length(m.pos(m.walker(h).vertex()) - m.pos(m.walker(h).opp().vertex()));
  }

  for (auto f : patch_faces) {
    interior_faces.insert(f);
  }

  // -------------------------
  // Remove all faces
  // -------------------------
  HMesh::FaceSet non_patch_faces;
  for (auto f : m.faces()) {
    if (patch_faces.find(f) == patch_faces.end()) {
      non_patch_faces.insert(f);
    }
  }

  for (auto f : non_patch_faces) {
    m.remove_face(f);
  }

  // -------------------------
  // Split patch faces
  // -------------------------
  // Here we split all faces into triangles. This is necessary, as the libigl code only accepts triangles to compute the Harmonic Map. 
  // However, init_HarmonicMap is also used to do an extrusion, when we built the mesh, and here there might be triangles, and maybe 
  // we do not have to split them, because they are already triangles. Therefore, we could do a check to see, if a face was actually 
  // a triangle face or something else, and then do not split if it is a triangle face. Because every time we add edges in the Harmonic Map,
  // it becomes skewed.
  for (auto f : patch_faces) {

    // Added on 10.1.2026
    if (face_valency(m,f) != 3) {
      auto f_2_split = f;
      auto new_v = m.split_face_by_vertex(f);

      f_v_map.insert({f_2_split, new_v});
    }
  }

  patch_faces.clear();
  for (auto f : m.faces()) {
    patch_faces.insert(f);
  }

  // -------------------------
  // Compute the boundary vertices
  // -------------------------
  HMesh::VertexSet bd_vs = boundary_verts(m, patch_faces);

  assert(bd_vs.find(ref_v) != bd_vs.end());

  bd_vertices = ccw_ordered_bd_vertices(m, patch_faces, ref_v);

  std::vector<HMesh::HalfEdgeID> bd_edges = ccw_ordered_bd_edges(m,  patch_faces, ref_v);

  // -------------------------
  // Set the matrices
  // -------------------------
  HMesh::VertexSet patch_vertices = all_verts(m, patch_faces);

  V.resize(patch_vertices.size(), 3); // All the vertices + 1 vertex in the middle of each face.

  F.resize(patch_faces.size(), 3); // The number of triangle faces we can make in each polygon. This number is equal to the number of vertices of each patch face.

  V_uv.resize(patch_vertices.size(), 2);

  // Vertex counter
  int v_counter = 0;
  for(auto v : patch_vertices) {
    vertex_igl_map.insert(std::make_pair(v, v_counter));
    igl_vertex_map.insert(std::make_pair(v_counter, v));

    V(v_counter, 0) = m.pos(v)[0];
    V(v_counter, 1) = m.pos(v)[1];
    V(v_counter, 2) = m.pos(v)[2];

    v_counter++;
  }


  // Face counter
  int f_counter = 0;
  for(auto f : patch_faces) {

    int ii = 0;
    circulate_face_ccw(m, f, [&] (HMesh::VertexID vn) {
      F(f_counter, ii) = vertex_igl_map.find(vn)->second;
      ii++;
    });
    face_igl_map.insert(std::make_pair(f, f_counter));
    igl_face_map.insert(std::make_pair(f_counter, f));


    f_counter++;
  }

  // find boundary vertices
  Eigen::VectorXi bnd(bd_vertices.size());
  int i = 0;
  for (auto v : bd_vertices) {
    bnd(i) = vertex_igl_map.find(v)->second;
    i++;
  }

  // ------------------------------------------
  // Compute the Harmonic Map
  // ------------------------------------------
  bnd_uv.resize(bd_vertices.size(), 2);

  igl::map_vertices_to_circle(V,bnd,bnd_uv);

  igl::harmonic(V, F,bnd,bnd_uv,1,V_uv); // The original one
  bool hasNaN = V_uv.array().isNaN().any();
  if (hasNaN) {
    igl::harmonic(F,bnd,bnd_uv,1,V_uv);
  }

  std::map<HMesh::VertexID, CGLA::Vec2d> v_uv_map;

  for (auto v : patch_vertices) {
    auto uv = V_uv.row(vertex_igl_map.find(v)->second);

    v_uv_map.insert({v, CGLA::Vec2d(uv(0), uv(1))});
  }

  // ------------------------------------------
  // Find the boundary face vertices
  // ------------------------------------------
  bd_F_vertices.resize(bd_edges.size(), 3);

  int bd_F_counter = 0;
  for (auto h : bd_edges) {

    auto v1 = m.walker(h).opp().vertex();
    auto v2 = m.walker(h).vertex();
    auto v3 = m.walker(h).next().vertex();

    bd_F_vertices(bd_F_counter,0) = vertex_igl_map.find(v1)->second;
    bd_F_vertices(bd_F_counter,1) = vertex_igl_map.find(v2)->second;
    bd_F_vertices(bd_F_counter,2) = vertex_igl_map.find(v3)->second;

    bd_F_counter++;
  }

  // ------------------------------------------
  // Save the extrusion mesh
  // ------------------------------------------
  if (!extrusion_name.empty()) {

    for (auto f : m.faces()) {
      if (patch_faces.find(f) == patch_faces.end()) {
        m.remove_face(f);
      }
    }

    for (auto v : m.vertices()) {
      if (vertex_igl_map.find(v) == vertex_igl_map.end()) {
        continue;
      }
      else {
        auto igl_v = vertex_igl_map.find(v)->second; 
        m.pos(v) = CGLA::Vec3d(V_uv(igl_v, 0), 0.0, V_uv(igl_v, 1));
      }
    }
  }


  // TC: Change 16th of June 2025
  patch_centre = uv_to_vertex(CGLA::Vec2d(0.0));

  centre_face_id = uv_to_face(CGLA::Vec2d(0.0));

  igl::per_face_normals(V, F, N_faces);


  patch_normal = CGLA::Vec3d(0.0);

  for(int j = 0; j < N_faces.rows(); j++) {
    patch_normal[0] += N_faces(j,0);
    patch_normal[1] += N_faces(j,1);
    patch_normal[2] += N_faces(j,2);
  }

  patch_normal /= N_faces.rows();
  patch_normal = normalize(patch_normal);

}

/* ----------------------------------------------------------------------- *
 * An overloaded function which recomputes the Harmonic Map 
 * ----------------------------------------------------------------------- */
void HarmonicMap::recompute_HarmonicMap(HMesh::Manifold m, HMesh::FaceSet patch_faces, HMesh::VertexID ref_v, int extrusion_id) {

  std::string extrusion_name = "recomputing_harmonic_map_" + std::to_string(extrusion_id);
  
  init_HarmonicMap(m, patch_faces, ref_v, extrusion_name);
}

/* ----------------------------------------------------------------------- *
 * Also an overloaded function which recomputes the Harmonic Map 
 * ----------------------------------------------------------------------- */
HarmonicMap::HarmonicMap(HMesh::Manifold m, HMesh::FaceSet patch_faces, HMesh::VertexID ref_v, std::string extrusion_name) {  

  init_HarmonicMap(m, patch_faces, ref_v, extrusion_name);
}

/* ----------------------------------------------------------------------- *
 * Another function which takes other arguments but also computes the Harmonic Map
 * ----------------------------------------------------------------------- */
HarmonicMap::HarmonicMap(HMesh::Manifold m, HMesh::FaceSet patch_faces, HMesh::FaceID centre_face, HMesh::FaceID right_face) {  

  //find ref v
  HMesh::HalfEdgeID bd_h = HMesh::InvalidHalfEdgeID;
  circulate_face_ccw(m, centre_face, [&]  (HMesh::HalfEdgeID h) {
    if(m.walker(h).opp().face() == right_face)
      bd_h = h;
    if(patch_faces.find(m.walker(h).opp().face()) == patch_faces.end())
      bd_h = h;
  });

  assert(bd_h != InvalidHalfEdgeID);

  while(patch_faces.find(m.walker(bd_h).opp().face()) != patch_faces.end()) {
    bd_h = m.walker(bd_h).opp().next().next().halfedge();
  }
  HMesh::VertexID ref_v = m.walker(bd_h).vertex();


  bd_v = ref_v;

  init_HarmonicMap(m, patch_faces, ref_v, "");
}


CGLA::Vec2d find_circle_intersection(CGLA::Vec2d P, CGLA::Vec2d d, CGLA::Vec2d C, double r) {
  // P : Start point of ray
  // d : direction of ray
  // C : Center of circle
  // r = radius of circle
  auto a = dot(d, d);
  auto m = P - C;
  auto b = 2.0 * dot(m, d);
  auto c = dot(m, m) - r * r;
  auto disc = b * b - 4.0 * a * c;

  auto sqrt_disc = sqrt(disc);
  auto t1 = (-b - sqrt_disc) / (2.0 * a);
  auto t2 = (-b + sqrt_disc) / (2.0 * a);
  if (t1 > 0) {
    return P + t1 * d;
  }
  else {
    return P + t2 * d;
  }

}


/* ----------------------------------------------------------------------- *
 * Compute a Harmonic Map of the base-patch faces which have disk topology. Add the face-loop faces afterwards.
 * ----------------------------------------------------------------------- */
std::map<HMesh::VertexID, CGLA::Vec2d> HarmonicMap::compute_normal_harmonic_map_with_face_loop(HMesh::Manifold m, std::string extrusion_name, bool save_extrusion, HMesh::VertexID correct_bd_v) {

  std::map<HMesh::VertexID, CGLA::Vec2d> uv_map;

  HMesh::FaceSet patch_faces = interior_faces;

  // TC: FaceSet containing extruded patch as well as face loop
  HMesh::FaceSet extruded_patch_faces;
  for (auto f : patch_faces) {
    extruded_patch_faces.insert(f);
  }

  HMesh::VertexID ref_v;
  if (correct_bd_v != HMesh::InvalidVertexID) {
    ref_v = correct_bd_v;
  }
  else {
    ref_v = bd_v;
  }
  // TC: Find boundary vertices
  std::vector<HMesh::VertexID> bd_vertices = ccw_ordered_bd_vertices(m, patch_faces, ref_v);

  // Change made on 4.9.2025
  std::vector<HMesh::HalfEdgeID> bd_edges = ccw_ordered_bd_edges(m, patch_faces, ref_v);

  // TC: Find Face loop faces
  for (auto h : bd_edges) {
    if (patch_faces.find(m.walker(h).opp().face()) == patch_faces.end() && m.walker(h).opp().face() != HMesh::InvalidFaceID) {
      extruded_patch_faces.insert(m.walker(h).opp().face());
    }
  }

  // TC: Remove faces that are not in the full face set
  HMesh::FaceSet faces_to_remove;
  for (auto f : m.faces()) {
    if (extruded_patch_faces.find(f) == extruded_patch_faces.end()) {
      //m.remove_face(f);
      faces_to_remove.insert(f);
    }
  }

  // Set new bd_edges
  bd_edges = ccw_ordered_bd_edges(m, patch_faces, ref_v);

  // Set new bd_vertices
  bd_vertices = ccw_ordered_bd_vertices(m, patch_faces, ref_v);

  for (auto f : faces_to_remove) {
    m.remove_face(f);
  }

  if (save_extrusion) {
      //bool save_status = obj_save(extrusion_name + "_before_normal_harmonic.obj", m);
  }

  Eigen::MatrixXd V_harmonic, V_uv_harmonic;
  Eigen::MatrixXd F_centres_harmonic;
  Eigen::MatrixXd bnd_uv_harmonic;

  Eigen::MatrixXi F_harmonic;

  std::map<HMesh::VertexID, int> vertex_igl_map_harmonic;
  std::map<int, HMesh::VertexID> igl_vertex_map_harmonic;

  std::map<HMesh::FaceID, int> face_igl_map_harmonic;
  std::map<int, HMesh::FaceID> igl_face_map_harmonic;

  HMesh::VertexSet patch_vertices_harmonic;

  int v_counter = 0;
  for (auto f : patch_faces) {
    circulate_face_ccw(m, f, [&] (HMesh::VertexID v) {
      patch_vertices_harmonic.insert(v);
      v_counter += 1;
    });
  }

  V_harmonic.resize(patch_vertices_harmonic.size() + patch_faces.size(), 3);
  V_uv_harmonic.resize(patch_vertices_harmonic.size() + patch_faces.size(), 2);

  F_harmonic.resize(v_counter, 3); // The number of triangle faces we can make in each polygon. This number is equal to the number of vertices of each patch face.

  int i = 0;
  for(auto v : patch_vertices_harmonic) {
    vertex_igl_map_harmonic.insert(std::make_pair(v, i));
    igl_vertex_map_harmonic.insert(std::make_pair(i, v));

    V_harmonic(i, 0) = m.pos(v)[0];
    V_harmonic(i, 1) = m.pos(v)[1];
    V_harmonic(i, 2) = m.pos(v)[2];

    i++;
  }


  int j = 0;

  // triangulate face set

  for(auto f : patch_faces) {

    V_harmonic(i, 0) = centre(m, f)[0];
    V_harmonic(i, 1) = centre(m, f)[1];
    V_harmonic(i, 2) = centre(m, f)[2];

    circulate_face_ccw(m, f, [&] (HMesh::HalfEdgeID h) {
      F_harmonic(j, 0) = vertex_igl_map_harmonic.find(m.walker(h).opp().vertex())->second;
      F_harmonic(j, 1) = vertex_igl_map_harmonic.find(m.walker(h).vertex())->second;
      F_harmonic(j, 2) = i;
      face_igl_map_harmonic.insert(std::make_pair(f, j));
      igl_face_map_harmonic.insert(std::make_pair(j, f));
      j++;
    });

    i++;
  }
  // TC: Finding boundary vertices
  std::vector<HMesh::VertexID> bd_verts = ccw_ordered_bd_vertices(m, patch_faces, ref_v);
  //std::reverse(bd_verts.begin(), bd_verts.end());

  Eigen::VectorXi bnd_harmonic(bd_verts.size());
  i = 0;
  for (auto v : bd_verts) {
//    cout<<v<<endl;
    bnd_harmonic(i) = vertex_igl_map_harmonic.find(v)->second;
    i++;
  }

  bnd_uv_harmonic.resize(bd_verts.size(), 2);
  igl::map_vertices_to_circle(V_harmonic, bnd_harmonic, bnd_uv_harmonic);

  igl::harmonic(V_harmonic, F_harmonic, bnd_harmonic, bnd_uv_harmonic, 1, V_uv_harmonic); // original one
  bool hasNaN = V_uv_harmonic.array().isNaN().any();
  if (hasNaN) {
    igl::harmonic(F_harmonic, bnd_harmonic, bnd_uv_harmonic, 1, V_uv_harmonic);
  }
  
  //V_uv_harmonic = compute_intrinsic_harmonic_map(F_harmonic, V_harmonic, bnd_harmonic, bnd_uv_harmonic);

  for (auto v : patch_vertices_harmonic) {
    int index = vertex_igl_map_harmonic.find(v)->second;
    double u_coor = V_uv_harmonic(index, 0);
    double v_coor = V_uv_harmonic(index, 1);
    m.pos(v) = CGLA::Vec3d(u_coor, 0.0, v_coor);
  }

  if (patch_faces.size() != extruded_patch_faces.size()) {

    // Map the boundary faces to the circle with radius r = 1.1
    double radius = 1.1;
    bd_edges = ccw_ordered_bd_edges(m, patch_faces, ref_v);
    auto bd_edges_set = boundary_hes(m, patch_faces);
    for (int ii = 0; ii < bd_edges.size(); ii++) {

      auto h = bd_edges[ii];
      auto next_h = m.walker(h).next().halfedge();

      // Case 1: The incident face is not double sided 
      if (bd_edges_set.find(next_h) == bd_edges_set.end()) {

        next_h = bd_edges[(ii+1)%bd_edges.size()];
        auto vec1 = m.pos(m.walker(h).vertex()) - m.pos(m.walker(h).opp().vertex());
        auto vec2 = m.pos(m.walker(next_h).vertex()) - m.pos(m.walker(h).vertex());

        auto vec1_2D = CGLA::Vec2d(vec1(0), vec1(2));
        auto vec2_2D = CGLA::Vec2d(vec2(0), vec2(2));

        auto vec = normalize(vec1_2D + vec2_2D);
        auto normal_vec = CGLA::Vec2d(-vec(1), vec(0));

        auto P_pos = m.pos(m.walker(h).vertex());
        auto P = CGLA::Vec2d(P_pos(0), P_pos(2));

        normal_vec.normalize();
        normal_vec *= std::copysign(1.0, dot(P,normal_vec));


        auto new_v_pos = find_circle_intersection(P, normal_vec, CGLA::Vec2d(0.0, 0.0), radius);
        auto incident_edge = m.walker(h).opp().next().next().halfedge();
        m.pos(m.walker(incident_edge).vertex()) = CGLA::Vec3d(new_v_pos[0], 0.0, new_v_pos[1]);

        // In case the faces on the boundary are disconnected for some reason, we do this: 
        auto disconected_vertex = m.walker(next_h).opp().next().vertex();
        if (m.walker(incident_edge).vertex() != disconected_vertex) {
          m.pos(disconected_vertex) = CGLA::Vec3d(new_v_pos[0], 0.0, new_v_pos[1]);
        }

      }
      // Case 2: The incident face is double sided - e.g. two of its edges have a boundary face, which is not part of patch_faces
      // But they have two different faces
      else if (m.walker(h).opp().face() != m.walker(next_h).opp().face()) {

        auto vec1 = m.pos(m.walker(h).vertex()) - m.pos(m.walker(h).opp().vertex());
        auto vec2 = m.pos(m.walker(next_h).vertex()) - m.pos(m.walker(h).vertex());

        auto vec1_2D = CGLA::Vec2d(vec1(0), vec1(2));
        auto vec2_2D = CGLA::Vec2d(vec2(0),vec2(2));

        auto vec = normalize(vec1_2D + vec2_2D);
        auto normal_vec = CGLA::Vec2d(-vec(1), vec(0));

        auto P_pos = m.pos(m.walker(h).vertex());
        auto P = CGLA::Vec2d(P_pos(0), P_pos(2));

        // Change the sign of the normal vector, so it points outwards
        // Quick change of normal vector - 16.9.2025
        //auto opp_P_vertex = m.pos(m.walker(h).prev().opp().vertex());
        //normal_vec = normalize(P - CGLA::Vec2d(opp_P_vertex(0), opp_P_vertex(2)));
        normal_vec *= std::copysign(1.0, dot(P,normal_vec));

        auto new_v_pos = find_circle_intersection(P, normal_vec, CGLA::Vec2d(0.0, 0.0), radius);
        auto incident_edge = m.walker(h).opp().next().next().halfedge();
        m.pos(m.walker(incident_edge).vertex()) = CGLA::Vec3d(new_v_pos[0], 0.0, new_v_pos[1]);
        
        auto disconected_vertex = m.walker(h).next().opp().next().vertex();
        if (m.walker(incident_edge).vertex() != disconected_vertex) {
          m.pos(disconected_vertex) = CGLA::Vec3d(new_v_pos[0], 0.0, new_v_pos[1]);
        }

      }

    }
  }


  if (save_extrusion) {
     
      //cout << "Saving the normal harmonic map extrusion as: " << extrusion_name + "_normal_harmonic_map.obj" << endl;
      //cout.flush();
      //bool save_status = obj_save(extrusion_name + "_normal_harmonic_map.obj", m);
  }

  for (auto v : all_verts(m, extruded_patch_faces)) {
  //for (auto v : all_verts(m, patch_faces)) {
    int index = vertex_igl_map_harmonic.find(v)->second;
    //CGLA::Vec2d uv(V_uv_lscm(index, 0), V_uv_lscm(index, 1));
    double x_coor = m.pos(v)[0];
    double z_coor = m.pos(v)[2];
    CGLA::Vec2d uv(x_coor, z_coor);
    uv_map.insert(std::make_pair(v, uv));
  }

  return uv_map;
}
