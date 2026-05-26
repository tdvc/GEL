/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * Made by Thor Christiansen
 * ----------------------------------------------------------------------- */

/**
 * @file extrusion.h
 * @brief All functions related to an extrusion
 */

#pragma once

#ifndef GEL_EXTRUSION_H
#define GEL_EXTRUSION_H

#include <GEL/CGLA/CGLA.h>
#include <GEL/HMesh/HMesh.h>
#include <GEL/Geometry/Graph.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/Geometry/KDTree.h>
#include <stack>
#include <GEL/HMesh/HarmonicMap.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <GEL/Geometry/KDTree.h>
#include <GEL/HMesh/load.h>

 /* ----------------------------------------------------------------------- *
  * Struct to store an extrusion element
  * ----------------------------------------------------------------------- */
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


/* ----------------------------------------------------------------------- *
 * The Generic Extrusion which is used to discretize the 2D continuous coordinates
 * of the Harmonic Map paramterized patches
 * ----------------------------------------------------------------------- */
class Generic_Extrusion{
    
    // The Manifold
    HMesh::Manifold m;
    
    // Set the orientation
    HMesh::FaceID origin_face;
    HMesh::FaceID right_face;
    HMesh::FaceID bd_f;
    HMesh::VertexID ref_v;
    HMesh::HalfEdgeID bd_h;
    
    // The faces in the middle of the shape
    HMesh::FaceSet base_face_set;
    
    // Boundary vertices of the base fase set
    std::vector<HMesh::VertexID> boundary_vertices;

    // Boundary vertices on the entire extrusion (so the rim of the generic extrusion)
    HMesh::VertexSet rim_vertices;
    
    // All the faces on the generic extrusion
    HMesh::FaceSet curr_ext_faces;
    
    // The positions of the mesh
    HMesh::VertexAttributeVector<CGLA::Vec3d> pos;
    
    // For the uv-map
    std::map<HMesh::VertexID, std::pair<int, HMesh::VertexID>> base_loop_v_map;
    
    // For the vertex to uv_map
    std::map<HMesh::VertexID, CGLA::Vec2d> vertex_uv_map;
 
    // GEL kDtree
    Geometry::KDTree<CGLA::Vec2d, int> uv_tree;
    
    
public:
    
    // Initialize the class function
    Generic_Extrusion();


    // TC: The purpose of this function is
    std::vector<HMesh::VertexID> get_neighbours(Eigen::MatrixXd curr_loop_V) {
            
        // Store the vertices of the found boundary loop
        std::vector<HMesh::VertexID> closest_vertices;

        // Loop through all the vertices in the loop
        for (int ii = 0; ii < curr_loop_V.rows(); ii++) {
            CGLA::Vec2d uv, uv_closest_point;
            int index;
            double inf = std::numeric_limits<double>::infinity();
            
            uv[0] = curr_loop_V(ii, 0);
            uv[1] = curr_loop_V(ii, 1);
            uv_tree.closest_point(uv, inf, uv_closest_point, index);
            
            // Insert the vertex ID into the list
            closest_vertices.push_back(HMesh::VertexID(index));
        }
        return closest_vertices;
    }
    
    HMesh::VertexID get_closest_vertex(CGLA::Vec2d uv_coordinate) {
        CGLA::Vec2d uv_closest_point;
        int index;
        double inf = std::numeric_limits<double>::infinity();
        uv_tree.closest_point(uv_coordinate, inf, uv_closest_point, index);
        
        return HMesh::VertexID(index);
    }

    HMesh::VertexID get_closest_vertex_on_the_rim(CGLA::Vec2d uv_coordinate) {
        CGLA::Vec2d uv_closest_point;
        int index;
        double inf = std::numeric_limits<double>::infinity();
        uv_tree.closest_point(uv_coordinate, inf, uv_closest_point, index);
        auto v = HMesh::VertexID(index);
        HMesh::VertexID v_on_rim;

        circulate_vertex_ccw(m, v,[&](HMesh::VertexID vn){
            if (rim_vertices.find(vn) != rim_vertices.end()) {
              v_on_rim = vn; 
            }
        });
        
        return v_on_rim;
    }
    
    CGLA::Vec2d vertex_to_uv(HMesh::VertexID v_index) {
        return vertex_uv_map.find(v_index)->second;
    }
    
    // Get functions
    std::map<HMesh::VertexID, CGLA::Vec2d> get_vertex_uv_map() {return vertex_uv_map;}
    
    HMesh::FaceSet get_curr_ext_faces() {return curr_ext_faces;}
     
    HMesh::VertexAttributeVector<CGLA::Vec3d> get_pos() {return pos;}
    
    HMesh::Manifold get_m() {return m;}
    
    std::vector<HMesh::VertexID> get_boundary_vertices() {return boundary_vertices;}

    HMesh::VertexSet get_rim_vertices() {return rim_vertices;}
    
    HMesh::FaceSet get_base_face_set() {return base_face_set;}
};



#endif