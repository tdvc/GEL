/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * Made by Thor Christiansen
 * ----------------------------------------------------------------------- */

/**
 * @file feq.h
 * @brief Compute the decomposition of a FEQ mesh
 */

#pragma once


#ifndef GEL_FEQ_H
#define GEL_FEQ_H

#include <iostream>
#include <fstream>


struct Extrusion {

    HMesh::FaceID origin_face;
    HMesh::FaceID right_face;

    // TC: The responsible_faces is the face_set which an extrusion contributed to for the finished mesh
    // stack_faces in this case are actually the face loop faces
    HMesh::FaceSet base_face_set, stack_faces, responsible_faces;

    HMesh::VertexID bd_v;
    
    HMesh::HalfEdgeID bd_h;

    CGLA::Vec3d bd_v_pos;

    CGLA::Vec3d right_vec;

    CGLA::Vec2d bd_vs;

    // Information for the bd_v responsible faceset
    int bd_vs_responsible_extrusion;
    HMesh::FaceSet ref_v_face_set;

    int id = -1;

    // TC: Variable made by TC
    double scale = 1.0;

    int stack_size = 0;

    HMesh::HalfEdgeID curr_h;

    // TC: The purpose of this boolean variable is basically to make a flag saying that an extrusion k + 1
    // is using the base_face_set of extrusion k. In that way we avoid having to rely on the yellow_vertices and face-loops, 
    // which might not work immediately, if hte base-patch might be kind a weird.
    std::map<int, bool> contributing_extrusions_entire_base_patch;

    // TC: Variable for storing which extrusions contribute with faces
    std::vector<int> contributing_extrusions;

    // A vector that stores which face ID belongs to which extrusion
    HMesh::FaceAttributeVector<int> extrusion_id;

    std::vector<int> child_nodes; // TC: The child nodes of the extrusion

    std::stack<HarmonicMap> hmap_stack;

    // This variable is used to store the boundary curves with (uv-coordinates) for each faceset that 
    // we need to select. For convenience, we store each boundary curve as a Nx2 Eigen::Matrix with doubles.
    std::map<int, std::vector<Eigen::MatrixXd>> next_gen_extrusion_loops;

    // This variable is used to store the boundary curves with (uv-coordinates) for each faceset that 
    // we need to select. For convenience, we store each boundary curve as a Nx2 Eigen::Matrix with doubles.
    std::map<int, std::vector<HMesh::FaceSet>> next_gen_extrusion_face_sets;

    std::map<int, CGLA::Vec2d> next_gen_extrusion_bd_vs;

    // The state of the mesh after killing the extrusion
    HMesh::Manifold m_state;

    // The state of the mesh before killing the extrusion
    HMesh::Manifold m_state_before;

};


#endif