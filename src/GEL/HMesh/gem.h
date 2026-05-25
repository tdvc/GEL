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


HMesh::VertexSet all_verts(const HMesh::Manifold& m, const HMesh::FaceSet& fs);

int face_valency(HMesh::Manifold &m, HMesh::FaceID f);

HMesh::VertexSet boundary_verts(const HMesh::Manifold& m, const HMesh::FaceSet& fs);

HMesh::HalfEdgeSet boundary_hes(const HMesh::Manifold &m, HMesh::FaceSet& fs);

HMesh::HalfEdgeSet extended_patch_edges(const HMesh::Manifold &m, HMesh::FaceSet& fs);

std::vector<HMesh::VertexID> find_boundary_vertices(const HMesh::Manifold &m, HMesh::FaceSet patch_faces, HMesh::VertexID ref_v);

std::vector<HMesh::HalfEdgeID> find_boundary_edges_from_ref_v(const HMesh::Manifold &m, HMesh::FaceSet patch_faces, HMesh::VertexID ref_v);

double geodesic_curvature(HMesh::Manifold& m, HMesh::Walker& w);

void smooth_faceset_lap_solve(HMesh::Manifold &m, HMesh::FaceSet faces);

HMesh::FaceID find_patch_center(HMesh::Manifold &m, HMesh::FaceSet fs);

std::vector<HMesh::FaceSet> find_groups_of_faces(HMesh::Manifold &m, HMesh::FaceSet fs);

HMesh::FaceSet find_interior_faces(HMesh::Manifold &m, HMesh::HalfEdgeID h);

std::tuple<double, HMesh::VertexID, std::vector<std::pair<HMesh::VertexID, double>>> compute_new_bd_v_with_high_curvature(HMesh::Manifold m, HMesh::FaceSet faces, std::map<HMesh::VertexID, CGLA::Vec2d> vertex_uv_map, HMesh::HalfEdgeSet possible_bd_edges, double cutoff);


#endif /* gem_hpp */