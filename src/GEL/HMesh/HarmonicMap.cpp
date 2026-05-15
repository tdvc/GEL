/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * Made by Thor Christiansen 
 * ----------------------------------------------------------------------- */

#include <HarmonicMap.h>




/* ----------------------------------------------------------------------- *
 * Given a 2D point with coordinates (uv), the point might be outside the parameterized mesh
 * if the mesh is very coarse. E.g. if the mesh is a right angled square, then the boundary vertices
 * will be located at (0,0), (0,1), (-1,0) and (-1,-1) in the harmonic map, and the 2D point might then be outside
 * the boundary. Consequently, we create a line going from the centre of the map (0,0) to the 2D point, and find out which
 * of the boundary edges intersect this line. This function thus returns the intersection point, and 
 * this intersection point will be used in downstream applications like finding the 3D position of the 2D point.
 * ----------------------------------------------------------------------- */
Vec2d HarmonicMap::find_bd_intersection(Vec2d uv) {

    // Loop thorugh each of the edges on the boundary of the paramterized patch
  for(int i = 0; i < bnd_uv.rows(); i++) {
    int curr_index = i;
    int next_index = (i+1)%bnd_uv.rows();

    Matrix <double, 1, 3> p;
    Matrix <double, 1, 3> r;
    Matrix <double, 1, 3> q;
    Matrix <double, 1, 3> s;

    p << uv[0], uv[1], 0.0;
    r << 0.0, 0.0, 0.0;

    q << bnd_uv(i, 0), bnd_uv(i, 1), 0.0;
    s << bnd_uv(next_index, 0), bnd_uv(next_index, 1), 0;

    double t, u;

    // Find the intersection point
    if(igl::segment_segment_intersect(p, r - p, q, s - q, t, u)) {

      return Vec2d(p(0,0) + t*(-p(0,0)), p(0,1) + t*(-p(0,1)));
    }
  }
  return Vec2d(0.0);
}



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
  for(auto h : h_boundary_hes(m, patch_faces)) { 
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

  /*
  cout << "Saving the mesh inside init_HarmonicMap as: mesh_" << extrusion_name << ".obj" << endl;
  cout.flush();
  obj_save("mesh_" + extrusion_name + ".obj", m);
  */

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

  /*
  cout << "Saving the mesh inside init_HarmonicMap as: mesh_after_split_" << extrusion_name << ".obj" << endl;
  cout.flush();
  obj_save("mesh_after_split_" + extrusion_name + ".obj", m);

  cout << "After splitting the faces, how many are there: " << patch_faces.size(); 
  cout << "Are they all in use: " << endl;
  for (auto f : patch_faces) {
    cout << "f: " << f << " is in use: " << m.in_use(f) << endl;
  }
  cout.flush();
  */

  // -------------------------
  // Compute the boundary vertices
  // -------------------------
  VertexSet bd_vs = h_boundary_verts(m, patch_faces);

  assert(bd_vs.find(ref_v) != bd_vs.end());

  bd_vertices = find_boundary_vertices(m, patch_faces, ref_v);

  std::vector<HMesh::HalfEdgeID> bd_edges = find_boundary_edges_from_ref_v(m,  patch_faces, ref_v);

  // -------------------------
  // Set the matrices
  // -------------------------
  VertexSet patch_vertices = all_verts(m, patch_faces);

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
    circulate_face_ccw(m, f, [&] (VertexID vn) {
      F(f_counter, ii) = vertex_igl_map.find(vn)->second;
      ii++;
    });
    face_igl_map.insert(std::make_pair(f, f_counter));
    igl_face_map.insert(std::make_pair(f_counter, f));


    f_counter++;
  }

  // find boundary vertices
  //cout << "Printing the boundary vertices: " << endl;
  VectorXi bnd(bd_vertices.size());
  int i = 0;
  for (auto v : bd_vertices) {
    //cout << "v: " << v << " m.pos(v): " << m.pos(v) << " V.row(i): " << V.row(vertex_igl_map.find(v)->second) << endl;
    bnd(i) = vertex_igl_map.find(v)->second;
    i++;
  }

  // ------------------------------------------
  // Compute the Harmonic Map
  // ------------------------------------------
  bnd_uv.resize(bd_vertices.size(), 2);

  igl::map_vertices_to_circle(V,bnd,bnd_uv);

  //cout << "bnd: " << bnd << endl;
  //cout << "bnd_uv: " << bnd_uv << endl;

  igl::harmonic(V, F,bnd,bnd_uv,1,V_uv); // The original one
  bool hasNaN = V_uv.array().isNaN().any();
  if (hasNaN) {
    igl::harmonic(F,bnd,bnd_uv,1,V_uv);
  }
  //
  //V_uv = compute_intrinsic_harmonic_map(F, V, bnd, bnd_uv);

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

  /*
  cout << "V: " << endl;
  for (int ii = 0; ii < V.rows(); ii++) {
    cout << "[" << V(ii,0) << "," << V(ii,1) << "," << V(ii,2) << "]," << endl;
  }
  
  cout << "F: " << endl;
  for (int ii = 0; ii < F.rows(); ii++) {
    cout << "[" << F(ii,0) << "," << F(ii,1) << "," << F(ii,2) << "]," << endl;
  }

  cout << "V_uv: " << V_uv << endl;
  cout.flush();
  */

  // ------------------------------------------
  // Save the extrusion mesh
  // ------------------------------------------
  if (!extrusion_name.empty()) {

    for (auto f : m.faces()) {
      if (patch_faces.find(f) == patch_faces.end()) {
        m.remove_face(f);
      }
    }

    //obj_save(("pre_" + extrusion_name + ".obj"), m);

    for (auto v : m.vertices()) {
      if (vertex_igl_map.find(v) == vertex_igl_map.end()) {
        continue;
      }
      else {
        auto igl_v = vertex_igl_map.find(v)->second; 
        m.pos(v) = Vec3d(V_uv(igl_v, 0), 0.0, V_uv(igl_v, 1));
      }
    }
    //cout << "Saving the mesh as: " << extrusion_name << endl;
    //cout.flush();
    //obj_save(extrusion_name + ".obj", m);
  }


  // TC: Change 16th of June 2025
  //cout << "About to set the patch_centre " << endl;
  //cout.flush();
  patch_centre = uv_to_vertex(Vec2d(0.0));

  centre_face_id = uv_to_face(Vec2d(0.0));

  igl::per_face_normals(V, F, N_faces);


  patch_normal = Vec3d(0);

  for(int j = 0; j < N_faces.rows(); j++) {
    patch_normal[0] += N_faces(j,0);
    patch_normal[1] += N_faces(j,1);
    patch_normal[2] += N_faces(j,2);
  }


  patch_normal /= N_faces.rows();
  patch_normal = normalize(patch_normal);

}

void HarmonicMap::recompute_HarmonicMap(HMesh::Manifold m_old, HMesh::FaceSet patch_faces, HMesh::VertexID ref_v, int extrusion_id) {  //:  m(m), patch_faces(patch_faces) {

  //cout << "Recomputing harmonic_map for extrusion: " << to_string(extrusion_id) << endl;

  std::string extrusion_name = "recomputing_harmonic_map_" + to_string(extrusion_id);
  
  HMesh::Manifold m = m_old; 
  init_HarmonicMap(m, patch_faces, ref_v, extrusion_name);

}

HarmonicMap::HarmonicMap(HMesh::Manifold m_old, HMesh::FaceSet patch_faces, HMesh::VertexID ref_v, string extrusion_name) {  //:  m(m), patch_faces(patch_faces) {

  HMesh::Manifold m = m_old;
  init_HarmonicMap(m, patch_faces, ref_v, extrusion_name);
}

HarmonicMap::HarmonicMap(HMesh::Manifold m_old, HMesh::FaceSet patch_faces, HMesh::FaceID centre_face, HMesh::FaceID right_face) {  //:  m(m), patch_faces(patch_faces) {

  HMesh::Manifold m = m_old;

  //find ref v
  HalfEdgeID bd_h = InvalidHalfEdgeID;
  circulate_face_ccw(m, centre_face, [&]  (HalfEdgeID h) {
    if(m.walker(h).opp().face() == right_face)
      bd_h = h;
    if(patch_faces.find(m.walker(h).opp().face()) == patch_faces.end())
      bd_h = h;
  });

  assert(bd_h != InvalidHalfEdgeID);

  while(patch_faces.find(m.walker(bd_h).opp().face()) != patch_faces.end()) {
    bd_h = m.walker(bd_h).opp().next().next().halfedge();
  }
  VertexID ref_v = m.walker(bd_h).vertex();


  bd_v = ref_v;

  init_HarmonicMap(m, patch_faces, ref_v, "");
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
  if (correct_bd_v != InvalidVertexID) {
    ref_v = correct_bd_v;
  }
  else {
    ref_v = bd_v;
  }
  // TC: Find boundary vertices
  std::vector<VertexID> bd_vertices = find_boundary_vertices(m, patch_faces, ref_v);

  // Change made on 4.9.2025
  vector<HalfEdgeID> bd_edges = find_boundary_edges_from_ref_v(m, patch_faces, ref_v);

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
  bd_edges = find_boundary_edges_from_ref_v(m, patch_faces, ref_v);

  // Set new bd_vertices
  bd_vertices = find_boundary_vertices(m, patch_faces, ref_v);

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

  VertexSet patch_vertices_harmonic;

  int v_counter = 0;
  for (auto f : patch_faces) {
    circulate_face_ccw(m, f, [&] (VertexID v) {
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

    circulate_face_ccw(m, f, [&] (HalfEdgeID h) {
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
  vector<VertexID> bd_verts = find_boundary_vertices(m, patch_faces, ref_v);
  //std::reverse(bd_verts.begin(), bd_verts.end());

  VectorXi bnd_harmonic(bd_verts.size());
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
    bd_edges = find_boundary_edges_from_ref_v(m, patch_faces, ref_v);
    auto patch_edges = extended_patch_edges(m, patch_faces);
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


        auto new_v_pos = find_cirle_intersection(P, normal_vec, CGLA::Vec2d(0.0, 0.0), radius);
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

        auto new_v_pos = find_cirle_intersection(P, normal_vec, CGLA::Vec2d(0.0, 0.0), radius);
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
