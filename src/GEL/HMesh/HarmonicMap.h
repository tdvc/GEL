/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * Made by Thor Christiansen
 * ----------------------------------------------------------------------- */

/**
 * @file HarmonicMap.h
 * @brief Compute a Harmonic Map of a mesh with disk topology
 */

#ifndef __MeshEditE__HarmonicMap__
#define __MeshEditE__HarmonicMap__

#include <iostream>
#include <GEL/CGLA/Vec2d.h>
#include <GEL/CGLA/Vec2f.h>
#include <GEL/CGLA/Vec3f.h>
#include <GEL/HMesh/Manifold.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>

class HarmonicMap
{

    Eigen::MatrixXd V, V_uv, V_ext, bnd_uv, N_faces;
  //(patch_vertices.size(), 3);
    Eigen::MatrixXi F; //(patch_faces.size()*2, 3);
    Eigen::MatrixXd F_centres;
    Eigen::MatrixXi bd_F_vertices; // Matrix storing the vretices of the boundary faces

    // TC: Added by TC
    std::vector<HMesh::VertexID> bd_vertices;
    //std::vector<CGLA::Vec3d> Vertex_normals; // A vector to store the normals for each vertex on the boundary of the patch
    std::map<int, CGLA::Vec3d> vertex_normals; // A map to store the normals for each vertex on the boundary of the patch
    HMesh::FaceSet interior_faces; // All the faces that make up the patch

    CGLA::Vec3d patch_centre, patch_normal;

    int centre_face_id;
    double bd_perim, patch_area;
    HMesh::VertexID bd_v;

    HMesh::HalfEdgeID bd_h; // The halfedge pointing to bd_v.

    std::map<HMesh::VertexID, int> vertex_igl_map;
    std::map<int, HMesh::VertexID> igl_vertex_map;

    std::map<HMesh::FaceID, int> face_igl_map;
    std::map<int, HMesh::FaceID> igl_face_map;

    std::vector<CGLA::Vec2d> bd_uvs;

    // Face to vertex map
    std::map<HMesh::FaceID, HMesh::VertexID> f_v_map;

    void init_HarmonicMap(HMesh::Manifold m, HMesh::FaceSet patch_faces, HMesh::VertexID ref_v, std::string extrusion_name);

public:

    int get_centre_face_id() { return centre_face_id; }
    void set_centre_face_id(int new_centre_face_id) {centre_face_id = new_centre_face_id;}

    double get_bd_perim() { return bd_perim; }
    void set_bd_perim(double new_bd_perim) {bd_perim = new_bd_perim;}

    double get_patch_area() { return patch_area;}
    void set_patch_area(double new_patch_area) {patch_area = new_patch_area;}
    
    HMesh::VertexID get_bd_v() { return bd_v; }
    void set_bd_v(HMesh::VertexID new_bd_v) {bd_v = new_bd_v;}

    HMesh::HalfEdgeID get_bd_h() { return bd_h; }
    void set_bd_h(HMesh::HalfEdgeID new_bd_h) {bd_h = new_bd_h;}
    
    // TC: Added by TC
    Eigen::MatrixXd get_V() { return V; }
    void set_V(Eigen::MatrixXd new_V) { V = new_V; }
    
    Eigen::MatrixXd get_V_uv() { return V_uv;}
    void set_V_uv(Eigen::MatrixXd new_V_uv) {V_uv = new_V_uv;}

    Eigen::MatrixXd get_V_ext() { return V_ext;}
    void set_V_ext(Eigen::MatrixXd new_V_ext) {V_ext = new_V_ext;}
    
    Eigen::MatrixXd get_bnd_uv() { return bnd_uv;}
    void set_bnd_uv(Eigen::MatrixXd new_bnd_uv) {bnd_uv = new_bnd_uv;}

    Eigen::MatrixXd get_N_faces() {return N_faces;}
    void set_N_faces(Eigen::MatrixXd new_N_faces) {N_faces = new_N_faces;}
    
    Eigen::MatrixXi get_F() { return F; }
    void set_F(Eigen::MatrixXi new_F) { F = new_F; }
    
    Eigen::MatrixXd get_F_centres() { return F_centres;}
    void set_F_centres(Eigen::MatrixXd new_F_centres) {F_centres = new_F_centres;}

    Eigen::MatrixXi get_bd_F_vertices() { return bd_F_vertices;}
    void set_bd_F_vertices(Eigen::MatrixXi new_bd_F_vertices) {bd_F_vertices = new_bd_F_vertices;}

    const std::vector<HMesh::VertexID>& get_bd_vertices() const { return bd_vertices;}
    void set_bd_vertices(std::vector<HMesh::VertexID> new_bd_vertices) {bd_vertices = new_bd_vertices;}

    const HMesh::FaceSet& get_patch_faces() const {return interior_faces;}
    void set_patch_faces(HMesh::FaceSet new_patch_faces) {interior_faces = new_patch_faces;}
    
    CGLA::Vec3d get_patch_centre() { return patch_centre;}
    void set_patch_centre(CGLA::Vec3d new_patch_centre) {patch_centre = new_patch_centre;}
    
    CGLA::Vec3d get_patch_normal() { return patch_normal;}
    void set_patch_normal(CGLA::Vec3d new_patch_normal) {patch_normal = new_patch_normal;}
    
    const std::map<HMesh::VertexID, int>& get_vertex_igl_map() const { return vertex_igl_map;}
    void set_vertex_igl_map(std::map<HMesh::VertexID, int> new_vertex_igl_map) {vertex_igl_map = new_vertex_igl_map;}
    
    const std::map<int, HMesh::VertexID>& get_igl_vertex_map() const { return igl_vertex_map;}
    void set_igl_vertex_map(std::map<int, HMesh::VertexID> new_igl_vertex_map) {igl_vertex_map = new_igl_vertex_map;}
    
    const std::map<HMesh::FaceID, int>& get_face_igl_map() const { return face_igl_map;}
    void set_face_igl_map(std::map<HMesh::FaceID, int> new_face_igl_map) {face_igl_map = new_face_igl_map;}
    
    const std::map<int, HMesh::FaceID>& get_igl_face_map() const { return igl_face_map;}
    void set_igl_face_map(std::map<int, HMesh::FaceID> new_igl_face_map) {igl_face_map = new_igl_face_map;}
    
    const std::vector<CGLA::Vec2d>& get_bd_uvs() const {return bd_uvs;}
    void set_bd_uvs(std::vector<CGLA::Vec2d> new_bd_uvs) {bd_uvs = new_bd_uvs;}

    const std::map<HMesh::FaceID, HMesh::VertexID>& get_f_v_map() const {return f_v_map;}
    void set_f_v_map(std::map<HMesh::FaceID, HMesh::VertexID> new_f_v_map) {f_v_map = new_f_v_map;}


    void input_ext_vertices(std::vector<HMesh::VertexID> v_ids, std::vector<CGLA::Vec3d> ext_pos, std::map<HMesh::FaceID, CGLA::Vec3d> f_centre_map) {
      //std::cout << "Inside input_ext_vertices: " << std::endl;
      // Set the extruded positions for the ordinary vertices in the patch
      V_ext.resize(V.rows(), 3);
      for(int i = 0; i < v_ids.size(); i++) {
        //std::cout << "For vertex: " << v_ids[i] << " ext_pos is: " << ext_pos[i] <<std::endl;
        int igl_id = vertex_igl_map.find(v_ids[i])->second;
        V_ext(igl_id, 0) = ext_pos[i][0];
        V_ext(igl_id, 1) = ext_pos[i][1];
        V_ext(igl_id, 2) = ext_pos[i][2];
      }

      // Set the extruded positions for the vertices in the middle of the faces
      for (auto it : f_centre_map) {
        auto f_id = it.first;
        auto ext_pos = it.second;

        auto v_id = f_v_map.find(f_id)->second;
        int igl_id = vertex_igl_map.find(v_id)->second;
        V_ext(igl_id, 0) = ext_pos[0];
        V_ext(igl_id, 1) = ext_pos[1];
        V_ext(igl_id, 2) = ext_pos[2];
      }
    }


    void save_loop(Eigen::MatrixXd loop_V);


    CGLA::Vec2d patch_vertex_uv(HMesh::VertexID v);

    std::map<HMesh::VertexID, CGLA::Vec2d> vertex_uv_map(HMesh::Manifold m,  std::string extrusion_name, bool save_extrusion);

    std::map<HMesh::VertexID, CGLA::Vec2d> compute_normal_harmonic_map_with_face_loop(HMesh::Manifold m, std::string extrusion_name, bool save_extrusion, HMesh::VertexID new_bd_v);

    std::map<HMesh::VertexID, CGLA::Vec2d> compute_normal_harmonic_map_with_face_loop_from_triangle_faces(HMesh::Manifold m, std::string extrusion_name, bool save_extrusion, HMesh::VertexID new_bd_v);

    void rotate_uvs(Eigen::Vector3d src_v, Eigen::Vector3d tgt_v);

    Eigen::MatrixXd uv_loop_from_faceset(const HMesh::Manifold& m, HMesh::FaceSet fs);

    CGLA::Vec2d find_bd_intersection(CGLA::Vec2d uv);

    int find_bd_edge_intersection(CGLA::Vec2d uv);

    CGLA::Vec3d find_disp_vec(CGLA::Vec2d bd_uv, CGLA::Vec2d v_uv);

    CGLA::Vec3d interp_disp_vec(CGLA::Vec2d v_uv, CGLA::Vec3d bd_pt, int va, int vb, int vc, double u, double v, double w);

    CGLA::Vec3d interp_disp_vec_extended(CGLA::Vec2d v_uv, CGLA::Vec3d bd_pt);

    CGLA::Vec3d interp_disp_vec_extended_edge(CGLA::Vec2d v_uv, HMesh::VertexID vA, HMesh::VertexID vB, double s, CGLA::Vec3d source_pt);

    void uv_to_bary_igl(CGLA::Vec2d v_uv, int& va, int& vb, int& vc, double& u, double& v, double& w);

    void vertex_to_bary(CGLA::Vec3d v_pos, HMesh::VertexID& va, HMesh::VertexID& vb, HMesh::VertexID& vc, double& u, double& v, double& w);

    void uv_to_bary(CGLA::Vec2d v_uv, HMesh::VertexID& va, HMesh::VertexID& vb, HMesh::VertexID& vc, double& u, double& v, double& w);

    int uv_to_face(CGLA::Vec2d v);

    CGLA::Mat3x3d get_patch_frame(CGLA::Vec3d bd_pt);

    CGLA::Mat3x3d source_patch_frame(int va, int vb, int vc, double u, double v, double w);
    CGLA::Mat3x3d target_patch_frame(HMesh::VertexID v, CGLA::Vec3d v_pos);
    std::vector<CGLA::Vec3d> find_normal_and_edge_vector_uv_coor(CGLA::Vec2d uv_coor);

    CGLA::Vec2d vertex_to_uv(CGLA::Vec3d v_pos);

    CGLA::Vec3d uv_to_vertex(CGLA::Vec2d v_uv);

    std::map<HMesh::VertexID, CGLA::Vec2d> get_v_uv_map(HMesh::Manifold m, std::string extrusion_name, bool save_extrusion);

    void recompute_HarmonicMap(HMesh::Manifold m, HMesh::FaceSet patch_faces, HMesh::VertexID ref_v, int extrusion_id);

    HarmonicMap(HMesh::Manifold m, HMesh::FaceSet patch_faces, HMesh::VertexID ref_v, std::string extrusion_name = "");

    HarmonicMap(HMesh::Manifold m, HMesh::FaceSet patch_faces, HMesh::FaceID centre_face, HMesh::FaceID right_face);


    

    HarmonicMap() {};



};


#endif /* defined(__MeshEditE__RadialMap__) */
