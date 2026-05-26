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

HMesh::VertexSet find_interior_vertices(const HMesh::Manifold &m, HMesh::FaceSet patch_faces);

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

bool vertex_vertex_connection(HMesh::Manifold m, HMesh::VertexID v1, HMesh::VertexID v2);

bool vertex_edge_connection(HMesh::Manifold m, HMesh::VertexID v1, HMesh::HalfEdgeID h);

double area_of_polygon(Eigen::MatrixXd curr_loop_V);

double area_of_curve_inside_circle(Eigen::MatrixXd curr_loop_V);

bool is_point_on_the_line(CGLA::Vec2d a, CGLA::Vec2d b, CGLA::Vec2d c);

std::pair<bool, double> do_edge_intersect_curve_t_value(HMesh::Manifold &m, HMesh::HalfEdgeID h1, CGLA::Vec2d a, CGLA::Vec2d b, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map);

bool do_edges_intersect(HMesh::Manifold &m, HMesh::HalfEdgeID h1, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map, HMesh::VertexID v1, HMesh::VertexID v2);

bool are_vertices_connected(HMesh::Manifold &m, HMesh::VertexID v1, HMesh::VertexID v2);

HMesh::FaceID remove_edge_add_face(HMesh::Manifold &m, HMesh::HalfEdgeID h);

std::tuple<HMesh::FaceID, HMesh::VertexID, HMesh::VertexID, HMesh::VertexID, double, double, double> locate_point_in_face(HMesh::Manifold &m, HMesh::FaceSet base_patch_faces, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map, CGLA::Vec2d point);

std::pair<HMesh::HalfEdgeID, double> locate_point_on_edge(HMesh::Manifold&m, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map, CGLA::Vec2d p);

bool vertex_face_connection(HMesh::Manifold m, HMesh::FaceSet base_patch_faces, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map, HMesh::VertexID v1, CGLA::Vec2d p);

void delaunay_triangulate_each_single_face2(HMesh::Manifold &m, HMesh::Manifold &m_copy, HMesh::FaceSet& base_patch_faces, HMesh::FaceSet& base_patch_faces_copy, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map, std::map<HMesh::VertexID, CGLA::Vec2d>& v_uv_map_copy, std::vector<std::tuple<int, HMesh::FaceSet, HMesh::VertexID, bool, std::vector<HMesh::VertexID>>> &curve_enclosed_faces, std::map<HMesh::VertexID, std::tuple<CGLA::Vec2d, HMesh::VertexID, HMesh::VertexID, double>>& new_vertices);


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
                                                                                                    >>>& face_map);

HMesh::FaceSet check_faceset_connectivity(HMesh::Manifold &m, HMesh::FaceSet face_set);

HMesh::FaceSet extract_face_set(std::vector<std::tuple<HMesh::FaceSet, bool, bool, std::string, int>>& faces_2_be_extruded);


// MST for finding faces
class MST{
  // The Curve
  Eigen::MatrixXd curve;

  // The Manifold
  HMesh::Manifold m;

  // The Heap
  std::priority_queue<
        std::pair<double,HMesh::HalfEdgeID>,
        std::vector<std::pair<double,HMesh::HalfEdgeID>>,
        std::greater<std::pair<double,HMesh::HalfEdgeID>>
    > heap;


  // Allowable edges
  // The edges we can use - either they are part of the base-patch or they are part of the face-loop
  HMesh::HalfEdgeSet allowable_edges;

  // MST edges
  HMesh::HalfEdgeSet mst_edges;

  // non-MST edges
  HMesh::HalfEdgeSet non_mst_edges;

  // non-useable edges
  HMesh::HalfEdgeSet unusable_edges;

  double edge_threshold;

  // Average edge length
  double avg_edge_length;

  // components
  std::map<HMesh::VertexID, HMesh::VertexID> components;

  // GEL kDtree
  Geometry::KDTree<CGLA::Vec2d, int> curve_tree;

  // Gradient vectors
  std::map<int, CGLA::Vec2d> gradient_vectors;

  // No of points pr. edge
  int no_edge_points = 10;

  double eps = 1e-5;

  // The positions of the vertex_uv_map
  std::map<HMesh::VertexID, CGLA::Vec2d> vertex_uv_map;

  HMesh::FaceSet curr_ext_faces;

  // Get the Ground truth curve
  std::vector<std::vector<double>> gt_curve;

  std::vector<std::vector<std::vector<double>>> gt_curve_segments;

  private:
    DTWComputer comp_;
    
  public:

    // Initialize the MST
    // Fix the issue in that m is actually the whole mesh and not the 2D harmonic map of the extrusion
    MST(HMesh::Manifold m, HMesh::FaceSet& curr_ext_faces, std::map<HMesh::VertexID, CGLA::Vec2d>& vertex_uv_map, Eigen::MatrixXd curr_loop_V);

    HMesh::HalfEdgeSet all_edges(HMesh::Manifold& m, HMesh::FaceSet fs) {
        HMesh::HalfEdgeSet hs;
        for(auto f: fs) {
            circulate_face_ccw(m, f, [&](HMesh::HalfEdgeID h) {
                hs.insert(h);
                hs.insert(m.walker(h).opp().halfedge());
            });
          }
        return hs;
    }

    std::vector<std::vector<double>> subsample_curve(std::vector<std::vector<double>> curve) {
      
      std::vector<std::vector<double>> subsampled_curve;
      int n = curve.size();

      for (int ii = 0; ii < n; ii++) {
        auto p1 = curve[ii];
        auto p2 = curve[(ii + 1) % n]; // wrap around to the first point
      
        for (int t = 0; t < no_edge_points; t++) {
            auto u = p1[0] + (p2[0] - p1[0]) * (static_cast<double>(t) / static_cast<double>(no_edge_points));
            auto v = p1[1] + (p2[1] - p1[1]) * (static_cast<double>(t) / static_cast<double>(no_edge_points));

            std::vector<double> element;
            element.push_back(u);
            element.push_back(v);
            subsampled_curve.push_back(element);
        }
      }
      return subsampled_curve;
    }

    std::vector<std::vector<double>> subsample_segment(std::vector<std::vector<double>> segment) {
      std::vector<std::vector<double>> subsampled_segment;
      int n = curve.size();

      auto p1 = segment[0];
      auto p2 = segment[1]; // wrap around to the first point
    
      for (int t = 0; t < no_edge_points; t++) {
          auto u = p1[0] + (p2[0] - p1[0]) * (static_cast<double>(t) / static_cast<double>(no_edge_points));
          auto v = p1[1] + (p2[1] - p1[1]) * (static_cast<double>(t) / static_cast<double>(no_edge_points));

          std::vector<double> element;
          element.push_back(u);
          element.push_back(v);
          subsampled_segment.push_back(element);
      }
      return subsampled_segment;
      
    }

    static std::pair<std::vector<double>, size_t>
    flatten2D_interleaved(const std::vector<std::vector<double>>& S) {
      if (S.empty()) throw std::invalid_argument("Series is empty");

      // Case A: N×2  (rows are samples)
      const bool maybe_Nx2 = (S[0].size() == 2);
      // Case B: 2×N  (rows/cols are dimensions)
      const bool maybe_2xN = (S.size() == 2);

      if (!maybe_Nx2 && !maybe_2xN)
        throw std::invalid_argument("Expected N×2 or 2×N");

      std::vector<double> flat;
      size_t N = 0;

      if (maybe_Nx2) {
        // Validate all rows have 2 elements
        for (const auto& row : S) {
          if (row.size() != 2) throw std::invalid_argument("N×2: every row must have size 2");
        }
        N = S.size();
        flat.reserve(N * 2);
        for (const auto& row : S) {
          flat.push_back(row[0]);
          flat.push_back(row[1]);
        }
      } else { // 2×N
        if (S.size() != 2) throw std::invalid_argument("2×N: outer size must be 2");
        if (S[0].size() != S[1].size())
          throw std::invalid_argument("2×N: both dimensions must have same length");
        N = S[0].size();
        flat.reserve(N * 2);
        for (size_t i = 0; i < N; ++i) {
          flat.push_back(S[0][i]); // x_i
          flat.push_back(S[1][i]); // y_i
        }
      }

      return {std::move(flat), N};
    }

    double dtw2d_dependent(const std::vector<std::vector<double>>& A,
                       const std::vector<std::vector<double>>& B) {
        auto [Af, n] = flatten2D_interleaved(A);
        auto [Bf, m] = flatten2D_interleaved(B);

        DTWSettings s = dtw_settings_default();
        s.inner_dist= 1; // Use Euclidean distance
        //s.psi = 0;
        // Optional tuning, e.g.:
        // s.use_pruning = 1;
        // s.window = 0;      // 0 = no window (Sakoe–Chiba band)
        // s.psi = 0;         // free begin/end
        // s.max_dist = -1.0; // pruning threshold (-1 disables)
        //dtw_distances_ndim_ptrs()


        return dtw_distance_ndim(Af.data(), n, Bf.data(), m, 2, &s);
    }


    // Compute distance between two vertices
    double compute_dist(HMesh::VertexID u, HMesh::VertexID v) {

      double alpha = 50.0;

      auto pos_u = vertex_uv_map.find(u)->second; //vertex_uv_map.find(u);
      auto pos_v = vertex_uv_map.find(v)->second; //vertex_uv_map.find(v);

      auto edge_vector = pos_v - pos_u;
      auto length_edge_vector = length(edge_vector);
      edge_vector.normalize();
      auto edge_normal = CGLA::Vec2d(-edge_vector[1], edge_vector[0]);
      edge_normal.normalize();
      //std::cout << "Edge normal: " << edge_normal << std::endl;
      
      double max_dist = 0.0;
      double min_dist = std::numeric_limits<double>::infinity();

      double sum = 0.0;
      for (int t = 0; t < no_edge_points; t++) {
        CGLA::Vec2d edge_point = pos_u + (pos_v - pos_u) * (static_cast<double>(t) / static_cast<double>(no_edge_points - 1));
        
        CGLA::Vec2d uv_closest_point;
        int index;
        double dist = std::numeric_limits<double>::infinity();
        curve_tree.closest_point(edge_point, dist, uv_closest_point, index);

        if (dist > max_dist) {
          max_dist = dist;
        }
        if (dist < min_dist) {
          min_dist = dist;
        }

        auto vp = gradient_vectors.find(index)->second;
       
        //double value = dot(edge_normal, gradient_vectors.find(index)->second);
        // Alternative formulation, where we find the parrallel vectors
        double value = dot(edge_normal, vp);
        //sum -= value;
        sum += (-dot(edge_normal, vp) + alpha * powf(dist,2.0));
      }
     // return (length_edge_vector/no_edge_points) * sum; // + powf(max_dist - min_dist, 2.0);
      //return sum + powf(max_dist - min_dist, 2.0);
      return sum; 
    }

    void get_nodes(HMesh::HalfEdgeID h, HMesh::VertexID& u, HMesh::VertexID& v) {
        u = m.walker(h).opp().vertex();
        v = m.walker(h).vertex();
    }

    void sort_edges() {

      //HMesh::HalfEdgeSet patch_edges = all_edges(m, curr_ext_faces);

      /*
      std::cout << "The allowable_edges are" << std::endl;
      for (auto h : allowable_edges) {
        std::cout << h << " ";
      }
      std::cout << std::endl;
      */

      //std::cout << "The dist for each edge is: " << std::endl;
      for (auto h : allowable_edges) {
        HMesh::VertexID u, v;
        get_nodes(h, u, v);

        double dist = compute_dist(u, v);
        //std::cout << "Edge h " << h << " between: " << v << vertex_uv_map.find(v)->second << " and: " << u << vertex_uv_map.find(u)->second << " has dist: " << dist << std::endl;
        heap.push({dist, h});
      }
    }

    bool is_connected(HMesh::HalfEdgeID h) {
      HMesh::VertexID u, v;
      get_nodes(h, u, v);

      return components.find(u) != components.end() && components.find(v) != components.end() && components.find(u)->second == components.find(v)->second;
    }

    void connect_nodes(HMesh::HalfEdgeID h) {
      HMesh::VertexID u, v;
      get_nodes(h, u, v);
      auto u_ID = components.find(u)->second;
      auto v_ID = components.find(v)->second;
      if (u_ID != v_ID) {
        for (auto& pair : components) {
          if (pair.second == u_ID) {
            pair.second = v_ID;
          }
        }
      }
    }

    void insert(HMesh::HalfEdgeID h) {
      connect_nodes(h);
      mst_edges.insert(h);
    }

    void remove_edge(HMesh::HalfEdgeID h) {
      if (mst_edges.find(h) != mst_edges.end()) {
        mst_edges.erase(h);
      }
      non_mst_edges.insert(h);
    }

    double compute_winding_number(CGLA::Vec2d p) {
      
      double winding_number = 0.0;
      for (int ii = 0; ii < curve.rows(); ii++) {
        auto p1 = curve.row(ii);
        auto p2 = curve.row((ii + 1) % curve.rows()); // wrap around to the first point

        CGLA::Vec2d v1 = CGLA::Vec2d(p1[0], p1[1]) - p;
        CGLA::Vec2d v2 = CGLA::Vec2d(p2[0], p2[1]) - p;

        v1.normalize();
        v2.normalize();

        double angle = acos(dot(v1, v2));
        winding_number += angle;



        // Accumulate the angle
        //winding_number += angle;
      }
      return winding_number / (2.0 * M_PI);

    }

    bool is_point_inside_curve(CGLA::Vec2d p) {
      double winding_number = compute_winding_number(p);
      return fabs(winding_number) > 0.90;
    }

    // The purpose of this function is to check, whether the edge is inside the curve by a magnitude of the 
    // average edge-length. If that is the case, then we should not use this edge in the tree.
    bool is_edge_useable(HMesh::HalfEdgeID h) {

      HMesh::VertexID u, v;
      get_nodes(h, u, v);

      auto pos_u = vertex_uv_map.find(u)->second; //vertex_uv_map.find(u);
      auto pos_v = vertex_uv_map.find(v)->second; //vertex_uv_map.find(v);

      CGLA::Vec2d uv_closest_point;
      int index;
      double dist_u = std::numeric_limits<double>::infinity();
      curve_tree.closest_point(pos_u, dist_u, uv_closest_point, index);

      double dist_v = std::numeric_limits<double>::infinity();
      curve_tree.closest_point(pos_v, dist_v, uv_closest_point, index);

      //double threshold = median_edge_length(m)*0.75;
      double threshold = edge_threshold*0.75;

      if (dist_u > threshold && dist_v > threshold) {
        // Check if the edge is inside or outside the curve by computing the winding number
        /*
        if (is_point_inside_curve(pos_u) || is_point_inside_curve(pos_v)) {
          return false;
        }
        */
        return false;
      }

      return true;
    }

    void Kruskal() {
      //std::cout << "Running Kruskal's algorithm" << std::endl;
      sort_edges();
      //std::cout << "Adding the following edges " << std::endl;
      while (!heap.empty()) {
        auto entry = heap.top();
        heap.pop();
        HMesh::HalfEdgeID h = entry.second;
        double dist = entry.first;
        
        // Insert the edge into the MST
        if (!is_connected(h)) {
            //std::cout << h << " ";
            insert(h);
        }
        else {
          auto opp_h = m.walker(h).opp().halfedge();

          if (non_mst_edges.find(h) == non_mst_edges.end() && non_mst_edges.find(opp_h) == non_mst_edges.end()) {
            non_mst_edges.insert(h);
          }

        } 
      }
      //std::cout << std::endl;
    }

    std::vector<int> findUndirectedCycle(HMesh::HalfEdgeID h) {
        // TC: This function finds one cycle in an undirected graph represented by halfedges
        // The function returns a vector of vertex IDs that form the cycle

        // TC: Using DFS to find a cycle
        // TC: The input is a halfedge ID, and we will find the cycle starting from this halfedge

      std::vector<std::pair<int,int>> edges;

      // Insert MST edges
      for (auto h : mst_edges) {
          HMesh::VertexID u, v;
          get_nodes(h, u, v);
          edges.push_back(std::make_pair(u.index, v.index));
      }
      // Insert the new edge
      HMesh::VertexID u, v;
      get_nodes(h, u, v);
      edges.push_back(std::make_pair(u.index, v.index));
        
        std::map<int, std::vector<int>> adj;

        for (auto [u, v] : edges) {
            adj[u].push_back(v);
            adj[v].push_back(u);
        }

        std::map<int, int> parent;
        std::set<int> visited;
        std::vector<int> cycle;

        std::function<bool(int,int)> dfs = [&](int u, int p)->bool {
            visited.insert(u);
            parent[u] = p;
            for (int v : adj.find(u)->second) {
                if (v == p) continue;
                if (!visited.count(v)) {
                    if (dfs(v, u)) return true;
                } 
                else {
                  // reconstruct cycle
                  std::vector<int> path;
                  path.push_back(v);
                  int cur = u;
                  while (cur != v && cur != -1) {
                      path.push_back(cur);
                      cur = parent.find(cur)->second;
                  }
                  path.push_back(v);
                  cycle = std::move(path);
                  return true;
                }
            }
            return false;
        };

        for (auto &kv : adj) {
            if (!visited.count(kv.first)) {
                parent[kv.first] = -1;
                if (dfs(kv.first, -1)) return cycle;
            }
        }
        return {};
    }

  void rotate_vector(std::vector<int>& v, int index) {
      if (index < 0 || index > static_cast<int>(v.size())) return; // safety check
      std::rotate(v.begin(), v.begin() + index, v.end());
  }

  std::vector<std::vector<double>> cycle_to_curve(std::vector<int> cycle) {
    std::vector<std::vector<double>> candidate_curve;
    for (int ii = 0; ii < cycle.size(); ii++) {
        CGLA::Vec2d uv = vertex_uv_map.find(HMesh::VertexID(cycle[ii]))->second;
        std::vector<double> element;
        element.push_back(uv[0]);
        element.push_back(uv[1]);
        candidate_curve.push_back(element);
    }
    return candidate_curve;
  }


  std::map<int, HMesh::HalfEdgeID> point_index_to_edge_map(std::vector<HMesh::VertexID> cycle) {
    std::map<int, HMesh::HalfEdgeID> pt_edge_map;
    
    int counter = 0;
    for (int ii = 0; ii < cycle.size(); ii++) {
      HMesh::VertexID u = cycle[ii];
      HMesh::VertexID v = cycle[(ii + 1) % cycle.size()];

      HMesh::HalfEdgeID edge; // Halfedge from u to v
      circulate_vertex_cw(m, u, [&](HMesh::HalfEdgeID h) {
        if (m.walker(h).vertex() == v) {
          edge = h;
        }
      });

      for (int t = 0; t < no_edge_points; t++) {
        if (t > 0) {
          pt_edge_map.insert({counter, edge});
        }
        counter++;
      }
    }
    return pt_edge_map;
  }

  std::pair<double, std::vector<int>> compute_DTW(HMesh::HalfEdgeID h) {
    std::vector<int> best_cycle;
    auto cycle = findUndirectedCycle(h);
    if (cycle.empty()) { 
      return {std::numeric_limits<double>::infinity(), best_cycle};
    }

    cycle.pop_back();

    double dtw_distance = std::numeric_limits<double>::infinity();

    auto original_cycle = cycle; 

    for (int ii = 0; ii < original_cycle.size(); ii++) {
      cycle = original_cycle;
      std::rotate(cycle.begin(), cycle.begin() + ii, cycle.end());

      // Subsample curve
      auto candidate_curve = cycle_to_curve(cycle);
      candidate_curve = subsample_curve(candidate_curve);

      // TC: Compute the DTW distance
      double value1 = dtw2d_dependent(gt_curve, candidate_curve);
      std::reverse(candidate_curve.begin() + 1, candidate_curve.end());
      double value2 = dtw2d_dependent(gt_curve, candidate_curve);

      double dist = std::min(value1, value2);
      if (dist < dtw_distance) {
        dtw_distance = dist;
        if (value2 < value1) {
          std::reverse(cycle.begin() + 1, cycle.end());
        }
        best_cycle.clear();
        for (auto v : cycle) {
          best_cycle.push_back(v);
        }
      }

    }
    return {dtw_distance, best_cycle};


  }

  std::tuple<std::vector<HMesh::VertexID>, double, HMesh::HalfEdgeID> find_curve() {
    Kruskal();
    //std::cout << "Computed Kruskal" << std::endl;
    double min_dist = std::numeric_limits<double>::infinity();
    std::vector<HMesh::VertexID> min_cycle;
    std::vector<int> cycle;
    HMesh::HalfEdgeID min_edge;

    // Loop over all the edges taht are not in the MST
    for (auto h : non_mst_edges) {
      //std::cout << "Processing edge: " << h;
      auto [dist, cycle] = compute_DTW(h);
      //std::cout << " Dist is: " << dist << std::endl;
      if (dist < min_dist) {
        min_dist = dist;
        min_cycle.clear();
        min_edge = h;
        for (auto v : cycle) {
          min_cycle.push_back(HMesh::VertexID(v));
        }
      }
    }
    //std::cout << "min_edge is: " << min_edge << " with distance: " << min_dist << std::endl;
    return std::make_tuple(min_cycle, min_dist, min_edge);
  }

  // Function that checks, whether a more suitable curve exists, by removing and adding edges
  std::vector<HMesh::VertexID> refine_curve() {
    auto [min_cycle, min_dist, min_edge] = find_curve();
    std::vector<int> cycle;
    for (auto v : min_cycle) {
      cycle.push_back(v.index);
    }
    auto candidate_curve = cycle_to_curve(cycle);
    candidate_curve = subsample_curve(candidate_curve);
    double current_dist = dtw2d_dependent(gt_curve, candidate_curve);

    auto p_to_h_map = point_index_to_edge_map(min_cycle);

    auto [Af, n] = flatten2D_interleaved(gt_curve);
    auto [Bf, m] = flatten2D_interleaved(candidate_curve);

    DTWSettings s = dtw_settings_default();
    s.inner_dist= 1; // Use Euclidean distance
    
    idx_t length_i;
    idx_t to_i[n + m];
    idx_t from_i[n + m];

    auto result = dtw_warping_path_ndim(Af.data(), n, Bf.data(), m, from_i, to_i, &length_i, 2, &s); 
    
    std::map<HMesh::HalfEdgeID, double> edge_cost;

    for (int ii = 0; ii < length_i; ii++) {
      if (p_to_h_map.find(to_i[ii]) != p_to_h_map.end()) {
        HMesh::HalfEdgeID h = p_to_h_map.find(to_i[ii])->second;

        CGLA::Vec2d gt_pt = CGLA::Vec2d(gt_curve[from_i[ii]][0], gt_curve[from_i[ii]][1]);
        CGLA::Vec2d cand_pt = CGLA::Vec2d(candidate_curve[to_i[ii]][0], candidate_curve[to_i[ii]][1]);

        if (edge_cost.find(h) == edge_cost.end()) {
          edge_cost.insert({h, 0.0});
        }
        edge_cost.find(p_to_h_map.find(to_i[ii])->second)->second += length(gt_pt - cand_pt);
      }
    }

    double max_edge_cost = -1.0;
    HMesh::HalfEdgeID max_edge;
    for (auto [h, cost] : edge_cost) {
      if (cost > max_edge_cost) {
        max_edge_cost = cost;
        max_edge = h;
      }
    }

    remove_edge(max_edge);
    insert(min_edge);

    return min_cycle;
}
  

  // --------------------
  // Get functions
  // --------------------

  HMesh::HalfEdgeSet get_mst_edges() {return mst_edges;}
  HMesh::HalfEdgeSet get_unusable_edges() {return unusable_edges;}
};


#endif /* gem_hpp */