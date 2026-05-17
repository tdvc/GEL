
#include "feq.h"
#include "gem.h"
#include "face_loop.h"
#include "HarmonicMap.h"

using namespace std;
using namespace CGLA;
using namespace Geometry;
using namespace HMesh;
using namespace Eigen;


// A global variable
std::map<int,Extrusion> extrusion_tree;

// TC: Basically, this function checks, whether the face-loop of h1 is alligned to the face-loop of h2. The two face-loops are alligned, 
// if every face in loop 1 has a neighbour in loop2 and vice versa. 
// If we start with h_begin = 1461 on a body of the DFAUST data set, it will be clear, that this function will evaluate to false in the for-loop below, 
// because not every face of h2 face-loop will be neighbour to a face of h1 face-loop
bool aligned_face_loop(Manifold &m, HalfEdgeID h1, HalfEdgeID h2) {

//fix aligned face loop.

    // TC: In order for a two faces-loops to be aligned, the face-loop of h1 should not contain any faces of the faces from the face-loop of h2
    // Otherwise h2's face-loop intersects h1's face-loop, and this is what Karran checks here.
    if(intersect_face_loop(m, h1, h2)) {
      return false;
    }

    // TC: For this function, Karran traces out the faces in the face-loop of an edge - e.g. edge h1. 
    // For each face f in the face-loop Karran checks, whether the face above f_above it or the face below it f_below are also part of the face-lop.
    // If f_above or f_below are also in the face-loop, then the face-loop is self-adjacent and it is therefore definitely not alligned with the face-loop of h2.
    if(self_adjacent(m, h1) || self_adjacent(m, h2)) {
      return false;
    }

    // ----------------------------------------------------
    // TC: In this part below Karran also checks, whether the face-loop of h2 (loop2) is alligned to the face-loop of h1 (loop1). 
    // The faces of loop1 are called curr_faces. 
    // Then Karran simply checks, whether every face in loop2 has a neighbour face that is part loop1 (named curr_faces)
    // Because he does not know which direction to check, he checks for faces above loop2 and below loop2. If the face-loops do not allign,
    // then there will be a case, where both a face above loop 2 and a face below loop 2 are not part of loop 1. They need both to be correct,
    // because if only one face either above or below is correct, then it might just be the wrong direction
    // ----------------------------------------------------
    HalfEdgeAttributeVector<int> touched(m.no_halfedges(), 0);

    FaceLoop loop1 = trace_face_loop(m, touched, h1);

    FaceLoop loop2 = trace_face_loop(m, touched, h2);

    FaceSet above_faces, curr_faces;

    int above_flag = 0;
    int below_flag = 0;

    for(auto h : loop1.hvec) {
        curr_faces.insert(m.walker(h).face());
    }
    for(auto h : loop2.hvec) {
        HMesh::FaceID f = m.walker(h).face();
        HMesh::FaceID above_face = m.walker(h).next().opp().face();
        HMesh::FaceID below_face = m.walker(h).prev().opp().face();

        if((curr_faces.find(above_face) == curr_faces.end()) && (curr_faces.find(below_face) == curr_faces.end()))
             return false;
    }

    // If the two face-loops are aligned, we return true.
    return true;
}

/* ----------------------------------------------------------------------- *
 * TC: This function checks, whether a face-loop is self-adjacent. It basically means that for a face-loop on a FEQ-mesh, 
 * if the neighbour face above a face in the face-loop, or a neighbour face below a face in the face-loop is also part of the face-loop
 * then the face-loop is self-adjacent. If the face-loop is self-adjacent the set of face in the face-loop does not create a nice ring around
 * the feature of the FEQ-mesh.
 * ----------------------------------------------------------------------- */
bool self_adjacent(HMesh::Manifold& m, HMesh::HalfEdgeID h) {
  HalfEdgeAttributeVector<int> touched(m.no_halfedges(), 0);
  FaceLoop curr_faceloop = trace_face_loop(m, touched, h);
  FaceSet loop_faces;

  for(auto h_loop : curr_faceloop.hvec) {
    loop_faces.insert(m.walker(h_loop).face());
  }

  for(auto h_loop : curr_faceloop.hvec) {
    HMesh::FaceID check_top = m.walker(h_loop).next().opp().face();
    HMesh::FaceID check_bottom = m.walker(h_loop).prev().opp().face();

    if(loop_faces.find(check_top) != loop_faces.end()) {
      return true;
    }

    if(loop_faces.find(check_bottom) != loop_faces.end()) {
      return true;
    }

  }
  return false;
}

/* ----------------------------------------------------------------------- *
 * If the function returns true, then the face-loop from h2 intersect the face-loop from h1.
 * ----------------------------------------------------------------------- */
bool intersect_face_loop(HMesh::Manifold &m, HMesh::HalfEdgeID h1, HMesh::HalfEdgeID h2) {
   HalfEdgeAttributeVector<int> touched(m.no_halfedges(), 0);

    FaceLoop loop1 = trace_face_loop(m, touched, h1);

    FaceLoop loop2 = trace_face_loop(m, touched, h2);

    HMesh::FaceSet above_faces, curr_faces;

    int above_flag = 0;
    int below_flag = 0;

    for(auto h : loop1.hvec) {
        curr_faces.insert(m.walker(h).face());
    }

    for(auto h : loop2.hvec) {
        HMesh::FaceID f = m.walker(h).face();
        if(curr_faces.find(f) != curr_faces.end())
            return true;
    }
    return false;
}

/* ----------------------------------------------------------------------- *
 * TC: The purpose of this function is the following. We check, whether the face-loop of edge h is the last face-loop we can decompose, because the next edge
 * of h (m.walker(h).next().opp().next().halfedge()) points to a face, which is part of face set that becomes a base-patch of the face-loop (also called leaf). 
 * The way that we check, whether the next edge points to a face, which is a base-patch and not a regular face-loop, is the following: 
 * 1) We go through all edges in the face-loop of edge h. 
 * 2) For each edge, we find the next edge, and this next edge is a part of a face-loop. The faces in this face-loop should intersect the face-loop of edge h. 
 * If there is no self-intersection with face-loop of edge h, we need to check, whether the new-face loop is self-adjacent. If the face-loop is not self-adjacent, 
 * then the faces cannot be a base-patch. Consequently, every face in the base-patch needs to either be part of a face-loop that intersects that last face-loop or be part of face-loop that
 * is self-adjacent.
 * ----------------------------------------------------------------------- */
bool check_leaf(HMesh::Manifold& m, HMesh::HalfEdgeID h, int pos_flag) {

    HMesh::HalfEdgeID curr_h = h;
    HMesh::HalfEdgeID h_next, curr_h_next;

    int leaf_flag = 0;

    HMesh::HalfEdgeAttributeVector<int> touched(m.no_halfedges(), 0);

    FaceLoop curr_faceloop = trace_face_loop(m, touched, curr_h);

    for (auto h_iter : curr_faceloop.hvec) {
        if(pos_flag == 0)
            curr_h_next = m.walker(h_iter).next().opp().next().halfedge();
        else
            curr_h_next = m.walker(h_iter).prev().opp().prev().halfedge();

        // It is a requirement that all faces in the next loop intersect the face-loop of edge h (also called curr_h), or that they are part of face-loop that intersects.
        if(!intersect_face_loop(m, h_iter, curr_h_next)) {
            leaf_flag = 1;
        }

        if(self_adjacent(m, curr_h_next)) {
          leaf_flag = 0;
        }
    }

    if(leaf_flag == 0) {
        return true;
    }
    else {
        return false;
    }
}


/* ----------------------------------------------------------------------- *
 * This function is the same as the check_leaf_stack function except for the fact that it returns
 * a std::vector of edges instead of a boolean value. But it just keep going from one-face loop to another, as long
 * as the two face-loops are aligned.
 * ----------------------------------------------------------------------- */
std::vector<HMesh::HalfEdgeID> find_face_loop_stack(HMesh::Manifold &m, HMesh::HalfEdgeID curr_h, int pos_flag) {

    HMesh::HalfEdgeID h_next;
    std::vector<HMesh::HalfEdgeID> feature_edges;
    feature_edges.push_back(curr_h);

    if(pos_flag == 0)
        h_next = m.walker(curr_h).next().opp().next().halfedge();
    else
        h_next = m.walker(curr_h).prev().opp().prev().halfedge();

    while(aligned_face_loop(m, curr_h, h_next)) {
        curr_h = h_next;
        feature_edges.push_back(curr_h);
        if(pos_flag == 0)
            h_next = m.walker(curr_h).next().opp().next().halfedge();
        else
            h_next = m.walker(curr_h).prev().opp().prev().halfedge();
    }

    return feature_edges;
}

/* ----------------------------------------------------------------------- *
 * Given a HalfEdgeID h, this function traces the face-loop, which h is part of
 * and finds all faces that make up a patch above the face-loop. The direction of the 
 * HalfEdge indicates, which way is 'up', and therefore which set of faces is above
 * the face-loop
 * ----------------------------------------------------------------------- */
HMesh::FaceSet find_interior_faces(HMesh::Manifold &m, HMesh::HalfEdgeID h) {


   HalfEdgeAttributeVector<int> touched(m.no_halfedges(), 0);
   FaceLoop l = trace_face_loop(m, touched, h);

   HMesh::FaceSet curr_faces, interior_faces;

   for(auto h : l.hvec) {
        curr_faces.insert(m.walker(h).face());
   }

    HMesh::HalfEdgeID h_above, h_below;

    h_above = m.walker(h).next().opp().next().halfedge();
    h_below = m.walker(h).prev().opp().prev().halfedge();

    //find interior face set
    FaceAttributeVector<int> face_status(m.no_faces(),0);
    FaceAttributeVector<int> face_visited(m.no_faces(),0);


    for(auto f : curr_faces) {
        face_status[f] = 1;
    }

    std::queue<HMesh::FaceID> fq;
    HMesh::FaceID leaf_face = HMesh::InvalidFaceID;
    if(check_leaf(m, h, 0)) {
        leaf_face = m.walker(h_above).face();
    }
    else if(check_leaf(m,h,1)) {
        leaf_face = m.walker(h_below).face();
    }

    fq.push(leaf_face);
    face_visited[leaf_face] = 1;

    int count = 0;

    while(!fq.empty()) {
        auto f = fq.front();
        fq.pop();
        interior_faces.insert(f);
        FaceSet nb_faces;
        circulate_face_ccw(m, f, [&](FaceID fn){
            nb_faces.insert(fn);
        });
        for(auto fn : nb_faces) {
          if(face_visited[fn] == 0 && face_status[fn] != 1 && interior_faces.find(fn) == interior_faces.end()) {
            fq.push(fn);
            face_visited[fn] = 1;
          }
        }
    
        count++;
    }
    return interior_faces;
}



/* ----------------------------------------------------------------------- *
 * Function which stores the next extrusion
 * ----------------------------------------------------------------------- */
HMesh::FaceSet store_next_gen_extrusions_hmap(HMesh::Manifold &m, Extrusion& ext) {

    VertexAttributeVector<Vec3d> old_pos =  m.positions_attribute_vector();

    HMesh::FaceSet base_face_set = ext.base_face_set; //get_faceset_from_patch_coords(m, ext.origin_face, ext.right_face, ext.patch_coords, true);

    HMesh::FaceSet full_face_set;

    for(auto f : base_face_set)  {
        full_face_set.insert(f);
    }

    for(auto f : ext.stack_faces) {
        full_face_set.insert(f);
    }


    // Find boundary of interior face

    smooth_faceset_lap_solve(m, base_face_set);

    // TC: This function basically computes the (u,v) coordinates for the base-face set
    // as well as finding the boundary edge and the reference vertex
    // Moreover, it also finds the patch normal, the 3D position of the patch centre based on the
    // (u,v) map, and the id of the face in the centre of the patch according to the (0,0) point
    // in the (u,v) coordinates.
    
    HarmonicMap hmap(m, base_face_set, ext.origin_face, ext.right_face);


    map<int, HMesh::FaceSet> ext_face_sets;

    auto extrusion_id = ext.extrusion_id;

    // TC: This function only changes something, if this code has been run before, or an extrusion has been carried out before. So it is only for the following extrusions that extrusion_id
    // will be different from -1.
    // Note: extrusion_id is actually extrusion_id_full from the kill_extrusion_hmap function.
    // For extrusion N Karran basically checks, if any of the faces in all the face loops that make up the extrusion also appear in the base-patch of another
    // extrusion. If so, he creates a set of the Face-ids, which appear in the other extrusion and give this set an id.
    for(auto f : full_face_set) {
    if (extrusion_id[f] != -1) {
        int next_ext_id = extrusion_id[f];
        if(ext_face_sets.find(next_ext_id) != ext_face_sets.end()) {
        ext_face_sets.find(next_ext_id)->second.insert(f);
        }
        else {
        HMesh::FaceSet next_fs;
        next_fs.insert(f);
        ext_face_sets.insert(std::make_pair(next_ext_id, next_fs));
        }
    }
    }
        
    // TC: Thiss function seems to find the (u,v) coordinates for all vertices in the base-face set.
    // The (u,v) coordinates are 2D parameter space coordinates in a circle it seems?
    map<HMesh::VertexID, CGLA::Vec2d> vertex_uv_map;

    //vertex_uv_map = hmap.compute_harmonic_map_with_face_loop(m, extrusion_name, save_mesh, HMesh::InvalidVertexID);
    vertex_uv_map = hmap.compute_normal_harmonic_map_with_face_loop(m, "", false, HMesh::InvalidVertexID);

    std::cout << "Inside store_next_gen_extrusions_hmap for extrusion " << ext.id << std::endl;

        // TC: The following function stores the extrusion loops, but this is only necessary, if we are killing an extrusion, which is not like a
    // finger on a hand or something similar to this.
      for(auto it = ext_face_sets.begin(); it != ext_face_sets.end(); it++) {
        int next_ext_id = it->first;
        FaceSet next_fs = it->second;

        Extrusion next_ext = extrusion_tree.find(next_ext_id)->second;

        // TC: If the next extrusion is sitting right on top of the this extrusion meaning that it uses all the faces, then we 
        // don't really need the face-loop or yellow vertices. Then next_fs should be equal to base_face_set.
        // next_fs is the face set that the next extrusion needs, and base_fase_set is the faces on top of the current extrusion.
       
        bool use_entire_base_patch = false;
        if (next_fs == base_face_set) {
            use_entire_base_patch = true;
        }
        extrusion_tree[next_ext_id].contributing_extrusions_entire_base_patch.insert(std::make_pair(ext.id, use_entire_base_patch));

        // TC: We might risk that the faces in the set next_fs are disconnected, so we need to find the groups of faces
        // This happens with the wrench.obj file
        auto groups_of_face_sets = find_groups_of_faces(m, next_fs);
        auto [temp1, temp2, face_groups] = split_groups_of_faces(groups_of_face_sets, base_face_set);

        HMesh::VertexID current_ref_v;

        // 
        std::vector<Eigen::MatrixXd> fs_loop_group;
        std::vector<HMesh::FaceSet> fs_group;

        HMesh::FaceSet current_ref_v_face_set;

        for (int kk = 0; kk < face_groups.size(); kk++) {
            auto faceset = face_groups[kk].first;

            // TC: It seems like Karran basically finds all boundary vertices with a specific valency and insert that into the boundary_vertex. I don't know what the * operation actually does. It seems like it gives him the first elemetn.
            HMesh::VertexID start_v = *boundary_verts(m, faceset).begin();


            // TC: It seems like this just finds the boundary vertices of the curve.
            std::vector<HMesh::VertexID> bd_verts = find_boundary_vertices(m, faceset, start_v);

            Eigen::MatrixXd fs_loop;

            fs_loop.resize(bd_verts.size(), 2);

            for(int i = 0; i < bd_verts.size(); i++) {

                if(bd_verts[i] == next_ext.bd_v || length(old_pos[bd_verts[i]] - next_ext.bd_v_pos) < 1e-5) {
                    current_ref_v = bd_verts[i];
                    current_ref_v_face_set = faceset;
                }

                if(vertex_uv_map.find(bd_verts[i]) != vertex_uv_map.end()) {

                    Vec2d curr_bd_uv = vertex_uv_map.find(bd_verts[i])->second;

                    fs_loop(i, 0) = curr_bd_uv[0];
                    fs_loop(i, 1) = curr_bd_uv[1];
                }

            }

            fs_loop_group.push_back(fs_loop);
            fs_group.push_back(faceset);
        }


        // Insert the face loop group into the next_gen_extrusion_loops map, and overwrite the previous one, if it exists.
        auto it1 = ext.next_gen_extrusion_loops.find(next_ext_id);
        if (it1 != ext.next_gen_extrusion_loops.end()) {
            it1->second = fs_loop_group;
        }
        else {
            ext.next_gen_extrusion_loops.insert(std::make_pair(next_ext_id, fs_loop_group));
        }
        // Insert the face loop group into the next_gen_extrusion_face_sets map, and overwrite the previous one, if it exists.
        auto it2 = ext.next_gen_extrusion_face_sets.find(next_ext_id);
        if (it2 != ext.next_gen_extrusion_face_sets.end()) {
            it2->second = fs_group;
        }
        else {
            ext.next_gen_extrusion_face_sets.insert(std::make_pair(next_ext_id, fs_group));
        }
        
        // TC: Insert the ids of the extrusions that are responsible for the faces for next_ext_id.
        extrusion_tree[next_ext_id].contributing_extrusions.push_back(ext.id);

        if(current_ref_v != InvalidVertexID) {


            auto element = ext.next_gen_extrusion_bd_vs.find(next_ext_id);
            if (element != ext.next_gen_extrusion_bd_vs.end()) {
                element->second = vertex_uv_map.find(current_ref_v)->second;
            }
            else {
                ext.next_gen_extrusion_bd_vs.insert(std::make_pair(next_ext_id, vertex_uv_map.find(current_ref_v)->second));
            }

            // TC: comment
            extrusion_tree[next_ext_id].bd_vs = vertex_uv_map.find(current_ref_v)->second;
            extrusion_tree[next_ext_id].bd_vs_responsible_extrusion = ext.id;
            extrusion_tree[next_ext_id].ref_v_face_set = current_ref_v_face_set;
        }
      }

      m.positions_attribute_vector() = old_pos;


      HMesh::FaceSet fs;
      return fs;
}


/* ----------------------------------------------------------------------- *
 * Update the next extrusion only based on one previous extrusion
 * ----------------------------------------------------------------------- */
void compute_new_bd_v_from_prev_ext(Extrusion &contrib_ext, Extrusion &next_ext) {
    // Find upward pointing edge
    HalfEdgeSet extend_patch_edges;
    HalfEdgeSet patch_edges = all_edges(contrib_ext.m_state, contrib_ext.base_face_set);
    for (auto h : patch_edges) {
        extend_patch_edges.insert(h);
        extend_patch_edges.insert(contrib_ext.m_state.walker(h).opp().halfedge());
    }
    HMesh::HalfEdgeID upward_edge;
    circulate_vertex_ccw(contrib_ext.m_state, contrib_ext.bd_v, [&](HMesh::HalfEdgeID h) {
        if (extend_patch_edges.find(h) == extend_patch_edges.end()) {
            upward_edge = contrib_ext.m_state.walker(h).opp().halfedge();
        }
    });

    stack<HarmonicMap> transform_stack = extrusion_tree[next_ext.id].hmap_stack;
    auto new_HarmonicMap = transform_stack.top();
    transform_stack.pop();
    HMesh::VertexID new_ref_v;

    auto bd_vertices = new_HarmonicMap.get_bd_vertices();
    
    auto incident_vertex = extrusion_tree[next_ext.id].m_state.walker(upward_edge).vertex();
    std::queue<HMesh::HalfEdgeID> edge_queue;
    circulate_vertex_ccw(extrusion_tree[next_ext.id].m_state, incident_vertex, [&](HMesh::HalfEdgeID h) {
        edge_queue.push(h);
    });


    while (true) {
        upward_edge = edge_queue.front();
        edge_queue.pop();
        auto incident_vertex = extrusion_tree[next_ext.id].m_state.walker(upward_edge).vertex();
        auto it = std::find(bd_vertices.begin(), bd_vertices.end(), incident_vertex);
        if (it != bd_vertices.end()) {
                new_ref_v = incident_vertex;
                break;
        } else {
            upward_edge = extrusion_tree[next_ext.id].m_state.walker(upward_edge).next().opp().next().halfedge();
            edge_queue.push(upward_edge);
        }
    }
    
    // --------------------------
    // Find the correct boundary vertex, because the vertex bd_v in contrib_ext is the one, after the extrusion contrib_ext has been killed. 
    // But when we need to compute the (u,v) coordinates for new_ref_v, which is the bd_v for the extrusion next_ext, then we need to find correct_bd_v for the extrusion contrib_ext
    // --------------------------
    HMesh::VertexSet temp_bd_verts = boundary_verts(contrib_ext.m_state_before, contrib_ext.base_face_set);
    HMesh::VertexID correct_bd_v;
    for (auto v : temp_bd_verts) {
        circulate_vertex_ccw(contrib_ext.m_state_before, v, [&] (HMesh::HalfEdgeID h) {
            if (contrib_ext.m_state_before.walker(h).vertex() == contrib_ext.bd_v) {
                correct_bd_v = v;
            }
        });
    }

    auto hmap = contrib_ext.hmap_stack.top();
    bool save_mesh = false;
    std::string extrusion_name = "ext_" + std::to_string(next_ext.id);
    if (next_ext.id >= 100 && save_mesh_patches) {
        save_mesh = true;
    }
    //std::map<HMesh::VertexID, CGLA::Vec2d> vertex_uv_map = hmap.compute_harmonic_map_with_face_loop(contrib_ext.m_state_before, extrusion_name, save_mesh, correct_bd_v);
    HMesh::Manifold m = contrib_ext.m_state_before;
    smooth_faceset_lap_solve(m, contrib_ext.base_face_set);
    std::map<HMesh::VertexID, CGLA::Vec2d> vertex_uv_map = hmap.compute_normal_harmonic_map_with_face_loop(m, extrusion_name, save_mesh, correct_bd_v);

    // --------------------------
    // Update the extrusion element
    // --------------------------
    HMesh::HalfEdgeID bd_h;
    HMesh::HalfEdgeSet next_ext_patch_edges = extended_patch_edges(extrusion_tree[next_ext.id].m_state, extrusion_tree[next_ext.id].base_face_set);
    circulate_vertex_ccw(extrusion_tree[next_ext.id].m_state, new_ref_v, [&](HMesh::HalfEdgeID h) {
        auto opp_h = extrusion_tree[next_ext.id].m_state.walker(h).opp().halfedge();
        if (extrusion_tree[next_ext.id].m_state.walker(opp_h).vertex() == new_ref_v && next_ext_patch_edges.find(opp_h) == next_ext_patch_edges.end()) {
            bd_h = opp_h;
        }
    });
    update_extrusion_element_with_bd_v(extrusion_tree[next_ext.id].m_state, next_ext, bd_h, vertex_uv_map.find(new_ref_v)->second);

}

/* ----------------------------------------------------------------------- *
 * The purpose of this function is to look at an extrusion and then change all the boundary vertices
 * It makes most sense to run this function after we have killed all the extrusions on our mesh, as we updated from the last extrusion
 *  and upwards.
 * ----------------------------------------------------------------------- */
void update_extrusions(HMesh::Manifold &m, HMesh::VertexID start_vertex) {

    std::queue<int> ext_queue;
    ext_queue.push(extrusion_tree.size() - 1);

    std::set<int> ext_visited;

    while (!ext_queue.empty()) {
        auto ext_id = ext_queue.front();
        ext_visited.insert(ext_id);
        ext_queue.pop();
        auto ext = extrusion_tree[ext_id];

        // Special case: If the extrusion is the last one in the tree
        // Then just compute the harmonic map and find the boundary vertex with the most curvature
        if (ext.contributing_extrusions.empty()) {
            HarmonicMap hmap(m, ext.base_face_set, ext.origin_face, ext.right_face);

            std::map<HMesh::VertexID, CGLA::Vec2d> vertex_uv_map = hmap.compute_normal_harmonic_map_with_face_loop(m, "start_patch", true, HMesh::InvalidVertexID);
   
            double angle; HMesh::VertexID new_bd_v; std::vector<std::pair<HMesh::VertexID, double>> peak_bd_vertices;
            if (start_vertex == HMesh::InvalidVertexID) {
                //cout << "Not using the provided start_vertex" << endl;
                double cutoff = 1.0;
                
                auto result = compute_new_bd_v_with_high_curvature(extrusion_tree[ext_id].m_state, 
                                                                    extrusion_tree[ext_id].base_face_set, 
                                                                    vertex_uv_map, boundary_hes(extrusion_tree[ext_id].m_state, 
                                                                        extrusion_tree[ext_id].base_face_set), cutoff);
                angle = std::get<0>(result);
                new_bd_v = std::get<1>(result);
                peak_bd_vertices = std::get<2>(result);
            }
            else {
                //cout << "Using the provided start_vertex: " << start_vertex << endl;
                new_bd_v = start_vertex;
            }
            
            HMesh::HalfEdgeID bd_h;
            for (auto h : boundary_hes(m, ext.base_face_set)) {
                if (m.walker(h).vertex() == new_bd_v) {
                    bd_h = h;
                }
            }
            update_extrusion_element_with_bd_v(extrusion_tree[ext_id].m_state, ext, bd_h, vertex_uv_map.find(new_bd_v)->second);
        }
        else {

            auto contrib_ext = extrusion_tree[extrusion_tree[ext_id].bd_vs_responsible_extrusion];

            bool use_entire_base_patch = extrusion_tree[ext_id].contributing_extrusions_entire_base_patch.find(contrib_ext.id)->second;
 
            // If the extrusion depend on the entire base patch of its parent extrusion
            // Then just move the boundary vertex X edge loops up
            if (extrusion_tree[ext_id].contributing_extrusions.size() == 1 && use_entire_base_patch && extrusion_tree[ext_id].bd_vs_responsible_extrusion == contrib_ext.id) {
                compute_new_bd_v_from_prev_ext(extrusion_tree[contrib_ext.id], extrusion_tree[ext_id]);
            }
            // If the extrusion does not depend on the entire base patch of its parent extrusion
            // Then we need to find the hmap of the parent extrusion, find the face set, then find the boundary vertices and select the one with biggest curvature
            else {
                compute_new_bd_v_from_partial_prev_ext(ext);
            }
        }
        // Add the next gen extrusions to the queue
        // But only add them, if we have already visited the extrusion, which is responsible for its bd_vs
        for (auto it = extrusion_tree[ext_id].next_gen_extrusion_loops.begin(); it != extrusion_tree[ext_id].next_gen_extrusion_loops.end(); ++it) {
            auto next_ext_id = it->first;

            // Check that all contributng extrusions have been visited
            bool all_contrib_visited = true;
            for (auto contrib_ext_id : extrusion_tree[next_ext_id].contributing_extrusions) {
                if (ext_visited.find(contrib_ext_id) == ext_visited.end()) {
                    all_contrib_visited = false;
                }
            }
            if (all_contrib_visited) {
                ext_queue.push(next_ext_id);
                continue;
            }
        }
         
    }

}


/* ----------------------------------------------------------------------- *
 * The Purpose of this function is to kill one extrusion element at a time
 * ----------------------------------------------------------------------- */
std::pair< std::map<int,Extrusion>, std::map<HMesh::FaceID, std::tuple<int, std::vector<HMesh::HalfEdgeID>, bool, HMesh::FaceSet>> > kill_individual_extrusions_hmap(HMesh::Manifold &m, 
                                                                                                                                                                    HMesh::HalfEdgeID h, 
                                                                                                                                                                    int pos_flag, 
                                                                                                                                                                    Generic_Extrusion &gen_ext, 
                                                                                                                                                                    HMesh::VertexID start_vertex) {
    
    //std::cout << "Inside kill_individual_extrusions_hmap" << std::endl;


    // 1) The face-ID X which becomes a face-loop face, after it has been a base-patch face
    // 2) The ID of the base-base extrusion
    // 3) Vector of halfedges on the boundary with the last edge being the one that is on the boundary
    // 4) bool - Whether we can keep tracking the edges, on the boundary.
    // 5) FaceSet - The face-set of the base-patch, where the face-ID X is a face-loop face.
    // E.g. for an extrusion X is a face-loop face and we then store the base-patch of this extrusion as the face set.
    std::map<HMesh::FaceID, std::tuple<int, std::vector<HMesh::HalfEdgeID>, bool, HMesh::FaceSet>> face_loop_prev_edge;
        
    //basically we need to go feature by feature. for each feature we first find patch center and x dir, and associated coords.
    
    //FaceAttributeVector<int> extrusion_id(-1);
    
    //TC:  These just appear to be vectors
    extrusion_id_full = FaceAttributeVector<int>(m.no_faces(), -1);
    
    // TC:This appears to be some sort of dictionary, where you can insert elements in
    // a tuple <key, value>, where the elements are sorted by the key-value.
    // Here the value is of type extrusion, so it seems like this is some sort of variable
    // that Karran defined.
    extrusion_tree.clear();

    int global_ext_id = 0;
    
    // TC: Some half-edge id values
    HMesh::HalfEdgeID curr_h, h_next, global_h;
    global_h = h;
    
    // TC: curr_patch_faces seems to be a C++ set, so some sort of list.
    HMesh::FaceSet curr_patch_faces;
    int curr_stack_size, break_flag = 0;
    // TC: Just a variable to be used later on for reducing the feature edges.
    int num_edges;
    // TC: An instance of the struct Extrusion, which contains some information
    // about the extrusion listed in the paper.
    Extrusion ext;
    
    int ext_count = 0;
    
    HMesh::HalfEdgeSet boundary_edges;
    double source_boundary_perim = 0.0;
    
    // TC: Vector to store the feature edges
    std::vector<HMesh::HalfEdgeID> global_feature_edges;
    std::vector<HMesh::HalfEdgeID> feature_edges, feature_edges_above, feature_edges_below;
    
    int counter = 0;

    bool save_intermediate_meshes = false;

    // Variable needed to store, which faces that were base-patch faces, changed into face-loop faces
    HMesh::FaceSet change_of_base_patch_faces_to_face_loop_faces;

    bool local_flag = false;
    
    //TC: Now we begin a while loop here.
    do {
        // TC: This seems to be a struct that contains an extrusion node, which is part of a list of multiple extrusion nodes
        // necessary for creating the extrusion graph
        // This is outside of the loop
        Extrusion curr_ext;
        
        // TC: It appears that the purpose of this function is to go through the face-loops
        // and see if the last face loop is a leaf face loop.
        if(global_feature_edges.empty() && check_leaf_stack(m, h, pos_flag)) {
            
            local_flag = true;
            
            // TC: The purpose of this function is to find a list of half-edges, which is on
            // the path of alligned face loops to the leaf face loop - the end of the feature.
            global_feature_edges = find_face_loop_stack(m, h, pos_flag);
            check_extrusion_on_leaf(m, global_feature_edges);     

    
        }
        
        // TC: It seems like the purpose of this part of this code is to find other features - e.g. if we set pos_flag to
        // another value, so we were not removing a simple feature but many different features.
        // So we go into this else-statement, if we the feature we are currently looking at did not have a leaf-face loop.
        else if (global_feature_edges.empty()) {    

            local_flag = false;
            
            // TC: It seems like the purpose of this function is to make a list of other features, if multiple features exist
            // and we cannot find a leaf face-loop
            curr_h = find_next_feature(m, h, pos_flag);
            
            auto feature_edges_above = find_face_loop_stack(m, curr_h , 0);
            auto feature_edges_below = find_face_loop_stack(m, curr_h , 1);
            
            //cout<<"found"<<endl;
            
            // TC: So either feature_edges_above or feature_edges_below will be 1, because curr_h comes from a face_loop, where there is self-intersection,
            // so in order to find the correct stack_size, we just add them together and subtract 1, because then we will be left with the correct
            // stack_size.
            curr_stack_size = feature_edges_above.size() + feature_edges_below.size() - 1;
            
            if(check_leaf(m, feature_edges_above.back(), 0)) {
                curr_h = feature_edges_above.back();
                // TC: If feature_edges_above contains the leaf edge, then we should put all edges into the global_feature_edges
                for (auto h : feature_edges_above) {
                    global_feature_edges.push_back(h);
                }
            }
            else if (check_leaf(m, feature_edges_below.back(), 1)) {
                curr_h = feature_edges_below.back();
                // TC: Contrary if feature_edges_below contains the leaf edge, then we should put those eges into hte global_feature_edges.
                for (auto h : feature_edges_below) {
                    global_feature_edges.push_back(h);
                }
            }
            //cout << "curr_h is: " << curr_h << endl;
            check_extrusion_on_leaf(m, global_feature_edges);
        }


        feature_edges_above.clear();
        feature_edges_below.clear();
        feature_edges_above.push_back(global_feature_edges.back());
        feature_edges_below.push_back(global_feature_edges.back());

        curr_h = global_feature_edges.back();

        global_feature_edges.pop_back();
        //global_feature_edges.clear(); // TC: Added by TC 27.3.2025 - 17:50
        
        // TC: So either feature_edges_above or feature_edges_below will be 1, because curr_h comes from a face_loop, where there is self-intersection,
        // so in order to find the correct stack_size, we just add them together and subtract 1, because then we will be left with the correct
        // stack_size.
        curr_stack_size = feature_edges_above.size() + feature_edges_below.size() - 1;
        
        HMesh::FaceSet loop_faces;
        
        // Just insert all the face loops into the vector loop_faces.
        for (auto h_loop : feature_edges_above) {
            HalfEdgeAttributeVector<int> touched(m.no_halfedges(), 0);
            FaceLoop curr_loop = trace_face_loop(m, touched, h_loop);
            for(auto h_loop_iter : curr_loop.hvec)
                loop_faces.insert(m.walker(h_loop_iter).face());
        }
        
        for (auto h_loop : feature_edges_below) {
            HalfEdgeAttributeVector<int> touched(m.no_halfedges(), 0);
            FaceLoop curr_loop = trace_face_loop(m, touched, h_loop);
            for(auto h_loop_iter : curr_loop.hvec)
                loop_faces.insert(m.walker(h_loop_iter).face());
        }
        
        
        curr_ext.stack_faces = loop_faces;
        
        curr_ext.stack_size = curr_stack_size;

        
        // TC: This function basically finds the interior faces of a patch
        curr_patch_faces = find_interior_faces(m, curr_h);
        
        // TC: Just assigning the interior face set to the extrusion.
        curr_ext.base_face_set = curr_patch_faces;

        // TC: This just appears to be a global variable, which keeps track of the extrusion number, so it is 0 to begin with,
        // and then it is increased with every extrusion. NOTE: First the value of global_ext_id is assigned to the extrusion id,
        // and then it is increased, so the first extrusion will have id 0, the second id 1 etc.
        curr_ext.id = global_ext_id++;
        

        // ----------------------------------------------------
        // Check when a base-patch face becomes a face-loop face, so we can mark it in the variable "face_patch_flag", which is a dictionary, where the key is the face id and the value is a boolean variable, which is true if the face is a base-patch face and false if it is a face-loop face. We need to do this, because when we kill the feature, we need to know which faces that were base-patch faces, but now are face-loop faces, because they are part of the stack of face loops that we are killing.
        HalfEdgeAttributeVector<int> touched2(m.no_halfedges(), 0);
        FaceLoop new_face_loop = trace_face_loop(m, touched2, curr_h);
        
        for (auto f : curr_ext.stack_faces) {
            if (face_patch_flag.find(f) != face_patch_flag.end() && face_patch_flag[f] == true) {
                change_of_base_patch_faces_to_face_loop_faces.insert(f); // Doesn't do anything. Just for debugging purposes

                for (auto h : new_face_loop.hvec) {
                    if (m.walker(h).face() == f) {
                        std::vector<HalfEdgeID> prev_edges = {m.walker(h).prev().opp().halfedge()};
                        face_loop_prev_edge.insert(std::make_pair(f, std::make_tuple(curr_ext.id, prev_edges, true, curr_patch_faces)));
                    }
                }
            }
        }

        // If we need to decompose a face-loop and the edges in the face-loop are in face_loop_prev_edge, 
        // then we need to update the previous edge of the face loop, because the face loop will be decomposed and therefore the previous edge will change. So we need to update it to make sure that we can still find the correct previous edge of the face loop after the decomposition.
        for (auto &pair : face_loop_prev_edge) {

            auto edge = std::get<1>(pair.second).back();
            auto should_update = std::get<2>(pair.second);

            auto base_patch_faces = std::get<3>(pair.second);

            if (curr_ext.base_face_set == base_patch_faces && curr_ext.stack_faces.find(m.walker(edge).face()) != curr_ext.stack_faces.end()) {
                auto new_edge = m.walker(edge).prev().prev().opp().halfedge();
                std::get<1>(pair.second).push_back(new_edge);
                std::get<2>(pair.second) = true;
                std::get<0>(pair.second) = curr_ext.id;
            }

        }
        

        // TC: Assign a boolean variable to each face in the patch - If it is a face-loop face or a base-patch faces. 
        // When we come all the way to the base-base patch, all faces should have been marked
        for (auto f : curr_ext.stack_faces) {
            if (face_patch_flag.find(f) == face_patch_flag.end()) {
                face_patch_flag.insert(std::make_pair(f, false)); // Mark as face-loop face
            }
            else {
                face_patch_flag[f] = false; // Mark as face-loop face
            }
        }
        for (auto f : curr_ext.base_face_set) {
            if (face_patch_flag.find(f) == face_patch_flag.end()) {
                face_patch_flag.insert(std::make_pair(f, true)); // Mark as base-patch face
            }
            else {
                face_patch_flag[f] = true; // Mark as base-patch face
            }
        }

        // -----------------------------------------------------

        
        // The quad face in the middle of the interior face set.
        // TC: Karran finds this origin_face (the center of the patch) by basically looping
        // through all the faces in the patch, then counting the number of steps you have to
        // take in the top-direction, left-direction, bottom-direction and right-direction
        // of the quad face to ensure that you reach the boundary. The face, where you have
        // to take the maximum number of steps will be the origin face (patch center).
        curr_ext.origin_face = find_patch_center(m, curr_patch_faces);
        //cout << "The patch center is: " << curr_ext.origin_face << endl;
        
        // TC: What does this right face do?
        curr_ext.right_face = m.walker(m.walker(curr_ext.origin_face).halfedge()).opp().face();
        

        
        curr_ext.m_state_before = m;

        curr_ext.extrusion_id = extrusion_id_full;
        
        // TC:
        //store_next_gen_extrusions_hmap(m, curr_ext, extrusion_id_full, extrusion_tree);
        HMesh::FaceSet returned_faces = store_next_gen_extrusions_hmap(m, curr_ext);
        
        //first need to make sure that all these functions are working as expected one by one
        //find interior faces, find patch center, find patch coords
        //print extrusion stuff and see if it makes sense
        
        // TC: So these are the interior faces. It seems like we just give them an id as to say, which extrusion
        // operation the interior faces belong to?
        for(auto ext_f : curr_patch_faces) {
            extrusion_id_full[ext_f] = curr_ext.id;
        }

        
        //cout<<"extrusion id stored, going to kill"<<endl<<h<<endl;
        
        stack<HarmonicMap> face_transform_stack;

        // TC: It seems like the purpose of this function is to remove all the face-loops and
        // return all the positions of the face-vertices in the face-loop?
        // local_flag is True
        // Edge h is the selected edge
        // curr_ext is the extrusion
        // If we cannot kill the edge h, which we initially selected, we need to kill a feature on that mesh, and therefore
        // we need to kill the stack from curr_h.
        if(local_flag) {
            //face_transform_stack = kill_feature_local_hmap(m, h, pos_flag, curr_ext);
            face_transform_stack = kill_feature_local_hmap(m, curr_h, pos_flag, curr_ext);
        }
        else {
            face_transform_stack = kill_individual_feature_hmap(m, curr_h, pos_flag, curr_ext);
        }

        // The current state of the mesh
        curr_ext.m_state = m;

        // TC: Compute the scale of the extrusion -- Flytning til før kill_feature_local_hmap
        source_boundary_perim = 0.0;
        boundary_edges = boundary_hes(m, curr_patch_faces);
        for (auto h : boundary_edges) {
            source_boundary_perim += length(m, h);
        }
        curr_ext.scale = source_boundary_perim;

        ext_count++;
        //curr_ext.per_face_transform_stack = per_coord_transform_stack(m, face_transform_stack, curr_ext);
        curr_ext.hmap_stack = face_transform_stack;
        
        //FaceSet ext_patch_faces = get_faceset_from_patch_coords(m, curr_ext.origin_face, curr_ext.right_face, curr_ext.patch_coords, true);
        
        // TC: Here Karran basically initializes the harmonic map. The map is simply initialized using the faces in the base patch "base_face_set",
        // and the origin face and the right face. In the function he can then initialize the other variable such as the boundary vertex, which he
        // finds as the boundary edge on the right face.
        //HarmonicMap hmap(m, curr_ext.base_face_set, curr_ext.origin_face, curr_ext.right_face, HMesh::Manifold(m));
        

        HarmonicMap hmap(m, curr_ext.base_face_set, curr_ext.origin_face, curr_ext.right_face);

         
        curr_ext.bd_v = hmap.get_bd_v();

        curr_ext.bd_v_pos = m.pos(curr_ext.bd_v);
        
        // TC: Insert responsible faces
        for(auto f : curr_ext.base_face_set) {
            curr_ext.responsible_faces.insert(f);
        }
        for (auto f : curr_ext.stack_faces) {
            curr_ext.responsible_faces.insert(f);
        }
        
        
        // TC: Insert the extrusion into the tree of extrusions necessary to kill the feature.
        // We insert a list of two items, where the first item is the id of the killed extrusion and the
        // other item is the extrusion itself.
        extrusion_tree.insert(std::make_pair(curr_ext.id, curr_ext));
        
        // TC: For debug purposes
        if (save_intermediate_meshes) {
            obj_save("test_mesh_" + to_string(counter) + ".obj", m);
        }
    
        ext = curr_ext;
        
        //return ext;
        counter+=1;
        
        // Find vertices on the generic extrusion
        //auto yellow_vertices = find_yellow_vertices(ext, gen_ext);
        //find_yellow_loops(ext, gen_ext, yellow_vertices);

        // TC: Så vi breaker kun, når local_flag er True.
        if(local_flag && !m.in_use(global_h)) {
        //if(local_flag) {
            break;
        }
    } while(true);
    std::cout<< "returning" <<endl;
    std::cout.flush();


    // TC: Update child nodes
    for (int ii = extrusion_tree.size() - 1; ii >= 0; ii--) {
        for (auto contrib_ext : extrusion_tree[ii].contributing_extrusions) {
            extrusion_tree[contrib_ext].child_nodes.push_back(extrusion_tree[ii].id);
        }
    }

    // TC: Updating the responsible faces for each extrusion
    for (int ii = extrusion_tree.size() - 1; ii >= 0; ii--) {
        for (auto contrib_ext : extrusion_tree[ii].contributing_extrusions) {
            HMesh::FaceSet new_responsible_faces;
            for (auto f : extrusion_tree[contrib_ext].responsible_faces) {
                if (extrusion_tree[ii].responsible_faces.find(f) == extrusion_tree[ii].responsible_faces.end()) {
                    new_responsible_faces.insert(f);
                }
            }
            extrusion_tree[contrib_ext].responsible_faces = new_responsible_faces;
        }
    }

    // Updating the extrusions with a new boundary vertex as refernece point
    update_extrusions(m, start_vertex);

    return {extrusion_tree, face_loop_prev_edge};
}