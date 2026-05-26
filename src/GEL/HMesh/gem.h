/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * Made by Thor Christiansen
 * ----------------------------------------------------------------------- */

#ifndef gem_hpp
#define gem_hpp


#include <GEL/HMesh/HMesh.h>
#include <GEL/HMesh/Manifold.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>

bool isInteger(const std::string& str);

HMesh::VertexSet all_verts(const HMesh::Manifold& m, const HMesh::FaceSet& fs);

HMesh::HalfEdgeSet all_edges(const HMesh::Manifold& m, const HMesh::FaceSet& fs);

HMesh::HalfEdgeSet all_edges_wo_duplicates(const HMesh::Manifold& m, const HMesh::FaceSet& fs);

int face_valency(HMesh::Manifold &m, HMesh::FaceID f);

HMesh::VertexSet boundary_verts(const HMesh::Manifold& m, const HMesh::FaceSet& fs);

HMesh::HalfEdgeSet boundary_hes(const HMesh::Manifold &m, HMesh::FaceSet& fs);

HMesh::HalfEdgeSet extended_patch_edges(const HMesh::Manifold &m, HMesh::FaceSet& fs);

std::vector<HMesh::VertexID> ccw_ordered_bd_vertices(const HMesh::Manifold &m, HMesh::FaceSet patch_faces, HMesh::VertexID ref_v);

std::vector<HMesh::HalfEdgeID> find_boundary_edges_from_ref_v(const HMesh::Manifold &m, HMesh::FaceSet patch_faces, HMesh::VertexID ref_v);

std::vector<HMesh::HalfEdgeID> ccw_ordered_bd_edges(const HMesh::Manifold &m, HMesh::FaceSet patch_faces, HMesh::VertexID ref_v);

double geodesic_curvature(HMesh::Manifold& m, HMesh::Walker& w);

void smooth_faceset_lap_solve(HMesh::Manifold &m, HMesh::FaceSet faces);

HMesh::FaceID find_patch_center(HMesh::Manifold &m, HMesh::FaceSet fs);

std::vector<HMesh::FaceSet> find_groups_of_faces(HMesh::Manifold &m, HMesh::FaceSet fs);

HMesh::FaceSet find_interior_faces(HMesh::Manifold &m, HMesh::HalfEdgeID h);

std::tuple<double, HMesh::VertexID, std::vector<std::pair<HMesh::VertexID, double>>> compute_new_bd_v_with_high_curvature(HMesh::Manifold m, HMesh::FaceSet faces, std::map<HMesh::VertexID, CGLA::Vec2d> vertex_uv_map, HMesh::HalfEdgeSet possible_bd_edges, double cutoff);

HMesh::HalfEdgeID find_connecting_edge(HMesh::Manifold m, HMesh::VertexID v1, HMesh::VertexID v2);

HMesh::FaceID find_shared_face(HMesh::Manifold &m, HMesh::VertexID v1, HMesh::VertexID v2);

double area_of_polygon(Eigen::MatrixXd curr_loop_V);

double area_of_curve_inside_circle(Eigen::MatrixXd curr_loop_V);

bool is_point_on_the_line(CGLA::Vec2d a, CGLA::Vec2d b, CGLA::Vec2d c);

std::pair<bool, double> do_edge_intersect_curve_t_value(HMesh::Manifold &m, HMesh::HalfEdgeID h1, CGLA::Vec2d a, CGLA::Vec2d b, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map);

bool do_edges_intersect(HMesh::Manifold &m, HMesh::HalfEdgeID h1, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map, HMesh::VertexID v1, HMesh::VertexID v2);

bool are_vertices_connected(HMesh::Manifold &m, HMesh::VertexID v1, HMesh::VertexID v2);

HMesh::FaceID remove_edge_add_face(HMesh::Manifold &m, HMesh::HalfEdgeID h);

std::tuple<HMesh::FaceID, HMesh::VertexID, HMesh::VertexID, HMesh::VertexID, double, double, double> locate_point_in_face(HMesh::Manifold &m, HMesh::FaceSet base_patch_faces, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map, CGLA::Vec2d point);

std::pair<HMesh::HalfEdgeID, double> locate_point_on_edge(HMesh::Manifold&m, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map, CGLA::Vec2d p);

void delaunay_triangulate_each_single_face2(HMesh::Manifold &m, HMesh::Manifold &m_copy, HMesh::FaceSet& base_patch_faces, HMesh::FaceSet& base_patch_faces_copy, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map_copy, std::vector<std::tuple<int, HMesh::FaceSet, HMesh::VertexID, bool, std::vector<HMesh::VertexID>>> &curve_enclosed_faces, std::map<HMesh::VertexID, std::tuple<CGLA::Vec2d, HMesh::VertexID, HMesh::VertexID, double>>& new_vertices);

#endif /* gem_hpp */