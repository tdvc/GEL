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

#include <GEL/CGLA/CGLA.h>
#include <GEL/HMesh/HMesh.h>
#include <GEL/Geometry/Graph.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/Geometry/KDTree.h>
#include <stack>
#include <GEL/HMesh/HarmonicMap.h>
#include <GEL/HMesh/extrusion.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <GEL/Geometry/KDTree.h>
#include <GEL/HMesh/load.h>


/* ----------------------------------------------------------------------- *
 * A simple data entity, which stores information about the previous and next extrusions
 * ----------------------------------------------------------------------- */
struct DAG_node {
  int id;
  std::deque<std::pair<int, std::string>> contributing_extrusions;
  std::string extrusion_text;
  std::string contributing_extrusion_text;
  std::set<int> next_extrusions;
  HMesh::FaceSet loop_faces;

  HMesh::Manifold base_base_mesh;
  std::vector<HMesh::HalfEdgeID> bd_edges_base_base_patch;
  std::map<HMesh::VertexID,CGLA::Vec2d> v_uv_map_base_base_patch;
};

/* ----------------------------------------------------------------------- *
 * A Simple graph data structure, which stores nodes
 * ----------------------------------------------------------------------- */
class Extrusion_DAG{
  private: 
      std::vector<DAG_node> nodes;

  public: 
      Extrusion_DAG() = default;

      // Method to add an item
    void push_back(const DAG_node& node) {
      nodes.push_back(node);
    }

    // Overload operator[] for access
    DAG_node& operator[](size_t index) {
        return nodes[index];
    }

    const DAG_node& operator[](size_t index) const {
        return nodes[index];
    }

    // Size of the collection
    size_t size() const {
        return nodes.size();
    }

    // Optional: begin/end to support range-based for loops
    auto begin() { return nodes.begin(); }
    auto end() { return nodes.end(); }

    auto begin() const { return nodes.begin(); }
    auto end() const { return nodes.end(); }

    // Optional: clear, erase, etc.
    void clear() {
        nodes.clear();
    }

    void pop_back() {
        nodes.pop_back();
    }

    bool empty() const {
        return nodes.empty();
    }
};



std::pair< std::map<int,Extrusion>, std::map<HMesh::FaceID, std::tuple<int, std::vector<HMesh::HalfEdgeID>, bool, HMesh::FaceSet>> > kill_individual_extrusions_hmap(HMesh::Manifold &m, 
                                                                                                                                                                    HMesh::HalfEdgeID h, 
                                                                                                                                                                    int pos_flag, 
                                                                                                                                                                    Generic_Extrusion &gen_ext, 
                                                                                                                                                                    HMesh::VertexID start_vertex);


#endif