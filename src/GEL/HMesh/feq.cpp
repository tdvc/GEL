
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
// Variable to store, whether a face is part of the base-base patch, or is on a face-loop
std::map<HMesh::FaceID, bool> face_patch_flag;

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


// TC: Basically, this function checks, whether the face-loop of h1 is alligned to the face-loop of h2. The two face-loops are alligned, 
// if every face in loop 1 has a neighbour in loop2 and vice versa. 
// If we start with h_begin = 1461 on a body of the DFAUST data set, it will be clear, that this function will evaluate to false in the for-loop below, 
// because not every face of h2 face-loop will be neighbour to a face of h1 face-loop
bool aligned_face_loop(HMesh::Manifold &m, HMesh::HalfEdgeID h1, HMesh::HalfEdgeID h2) {

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

bool check_contains(HMesh::Manifold& m, HMesh::HalfEdgeID h, HMesh::HalfEdgeID invalid_edge) {

    HalfEdgeAttributeVector<int> touched(m.no_halfedges(), 0);
    FaceLoop curr_faceloop = trace_face_loop(m, touched, h);

    for(auto h_iter: curr_faceloop.hvec) {
        //cout<<"H iter, opp and invalid are: "<<h_iter<<" "<<m.walker(h_iter).opp().halfedge()<<" "<<invalid_edge<<endl;
        if(h_iter == invalid_edge || m.walker(h_iter).opp().halfedge() == invalid_edge) {
            //cout<<"False"<<endl;
            return false;
        }
    }
    return true;
}

// TC: The purpose of this function is to filter out those edges, which would create dangling faces
std::vector<HMesh::HalfEdgeID> filter_the_leaf_edges(HMesh::Manifold &m, std::vector<HMesh::HalfEdgeID> leaf_edges, int pos_flag) {
    
    std::vector<HMesh::HalfEdgeID> filtered_leaf_edges;

    // Edge -> (base_patch_faces, face_loop_faces)
    std::map<HMesh::HalfEdgeID, std::pair<HMesh::FaceSet, HMesh::FaceSet>> h_to_fs_set;

    HalfEdgeSet edges_2_keep; // To ensure that there are no duplicates 
    for (auto h : leaf_edges) {

        bool keep_edge = true;
        for (auto hh : edges_2_keep) {
            if (!check_contains(m, h, hh)) {
                keep_edge = false;
            }
        }
        if (keep_edge) {
            edges_2_keep.insert(h);
        }
    }

    // Do the check for dangling faces
    for (auto h : edges_2_keep) {
        auto feature_edges = find_face_loop_stack(m, h, pos_flag);
    
        auto base_patch_faces = find_interior_faces(m, feature_edges.back());

        HMesh::FaceSet face_loop_faces;
        for (auto hh : feature_edges) {
            HalfEdgeAttributeVector<int> touched(m.no_halfedges(), 0);
            FaceLoop curr_faceloop = trace_face_loop(m, touched, hh);

            face_loop_faces.insert(curr_faceloop.face_loop_faces.begin(), curr_faceloop.face_loop_faces.end());
        }

        h_to_fs_set.insert(std::make_pair(h, std::make_pair(base_patch_faces, face_loop_faces)));
    }

    for (auto h : edges_2_keep) {
        bool are_all_base_patch_faces_covered = false;
        auto base_patch_faces = h_to_fs_set.find(h)->second.first;

        for (auto hh : edges_2_keep) {
            if (h == hh) {
                continue;
            }
            else {
                auto face_loop_faces = h_to_fs_set.find(hh)->second.second;

                if (std::includes(face_loop_faces.begin(), face_loop_faces.end(), base_patch_faces.begin(), base_patch_faces.end())) {
                    are_all_base_patch_faces_covered = true;
                }
            }
        }

        
        if (are_all_base_patch_faces_covered) {
            filtered_leaf_edges.push_back(h);
        }
        
    }

    if (filtered_leaf_edges.size() == 0) {
        for (auto h : edges_2_keep) {
            filtered_leaf_edges.push_back(h);
        }
    }


    return filtered_leaf_edges;
}


/* ----------------------------------------------------------------------- *
 * In this function Karran basically checks, whether h is at the bottom of a single feature, which can be decomposed. 
 * ----------------------------------------------------------------------- */
bool check_leaf_stack(HMesh::Manifold& m, HMesh::HalfEdgeID h, int pos_flag) {

    HalfEdgeID curr_h = h;
    HalfEdgeID h_next, curr_h_next;

    int leaf_flag = 0;

    if(pos_flag == 0) {
        h_next = m.walker(curr_h).next().opp().next().halfedge();
    }
    else {  
        h_next = m.walker(curr_h).prev().opp().prev().halfedge();
    }

    while(aligned_face_loop(m, curr_h, h_next)) {
        curr_h = h_next;
        if(pos_flag == 0)
            h_next = m.walker(curr_h).next().opp().next().halfedge();
        else
            h_next = m.walker(curr_h).prev().opp().prev().halfedge();
    }

    // In order to be a leaf, then every face in the leaf/base-patch should be part of a face-loop, which intersects the face-loop of curr_h or be a self-adjacent face-loop
    return check_leaf(m, curr_h, pos_flag);

}


// TC: Sometimes we might encounter the situation that we have a flat plane and on one side of the plane there might be multiple extrusions on that plane
// but Karran's framework does not detect that, because the check_leaf function does not handle this case, since some edges are not next to the extrusion.
// An example of this behavoiur is experieneced with the feline_demo, where two extrusions are treated as a base-fase-set, which creates some unfortunate 
// situations. Try decomposing the feline_demo with pos_flag = 0 and h = 6651
void check_extrusion_on_leaf(HMesh::Manifold &m, std::vector<HMesh::HalfEdgeID> &edges) {

    HMesh::HalfEdgeID curr_h, curr_h_next, leaf_edge;
    curr_h = edges.back();

    HalfEdgeAttributeVector<int> touched(m.no_halfedges(), 0);
    FaceLoop curr_faceloop = trace_face_loop(m, touched, curr_h);

    for (auto h_iter : curr_faceloop.hvec) {
        curr_h_next = m.walker(h_iter).next().opp().next().halfedge();

        if (!intersect_face_loop(m, curr_h, curr_h_next) && !self_adjacent(m, curr_h_next)) {
            leaf_edge = curr_h_next;
            break;
        }

    }
    if (leaf_edge != InvalidHalfEdgeID) {
        edges.push_back(leaf_edge);
        curr_h_next = m.walker(leaf_edge).next().opp().next().halfedge();
        while (!intersect_face_loop(m, curr_h, curr_h_next) && !self_adjacent(m, curr_h_next)) {
            edges.push_back(curr_h_next);
            curr_h_next = m.walker(curr_h_next).next().opp().next().halfedge();
        }

    }
}


// TC: The purpose of this function is simply to find face-loops above the face-loop of edge h
// which we can trace out in order to find new potential face-loops which we can decompose.
// So edge h will be part of the last face-loop which are aligned.
std::vector<HMesh::HalfEdgeID> find_unique_jns(Manifold&m, HalfEdgeID h, int pos_flag) {

    std::vector<HalfEdgeID> unique_jns;

    HMesh::HalfEdgeID h_next, curr_h, curr_h_next;

    curr_h = h;
    int visited;

    FaceLoop curr_faceloop;
    HalfEdgeAttributeVector<int> touched(m.no_halfedges(), 0);

    curr_faceloop = trace_face_loop(m, touched, curr_h);

    
    // TC: The purpose of this function is to find edges in the next face-loop, which goes in other
    // directions than the face-loop of edge h
    for (auto h_iter : curr_faceloop.hvec) {
        HMesh::HalfEdgeSet h_next_set;
        if(pos_flag == 0) {
          //FaceID loop_f =
          HMesh::HalfEdgeID top_h = m.walker(h_iter).next().opp().halfedge();
          HMesh::VertexID center_v = m.walker(top_h).vertex();

          HMesh::HalfEdgeID end_h = m.walker(h_iter).opp().prev().halfedge();

          assert(center_v == m.walker(end_h).vertex());

          HMesh::HalfEdgeID next_h = m.walker(top_h).next().opp().halfedge();

          h_next_set.insert(m.walker(next_h).opp().halfedge());

          while(m.walker(next_h).next().opp().halfedge() != end_h) {
            next_h = m.walker(next_h).next().opp().halfedge();
            assert(center_v == m.walker(next_h).vertex());
            h_next_set.insert(m.walker(next_h).opp().halfedge());
          }
        }
            //curr_h_next = m.walker(h_iter).next().opp().next().halfedge();
        else {
          HMesh::HalfEdgeID bottom_h = m.walker(h_iter).prev().halfedge();
          HMesh::VertexID center_v = m.walker(bottom_h).vertex();

          HMesh::HalfEdgeID end_h = m.walker(h_iter).opp().next().opp().halfedge();

          assert(center_v == m.walker(end_h).vertex());

          HMesh::HalfEdgeID next_h = m.walker(bottom_h).opp().prev().halfedge();

          h_next_set.insert(m.walker(next_h).halfedge());

          while(m.walker(next_h).opp().prev().halfedge() != end_h) {
            next_h = m.walker(next_h).opp().prev().halfedge();
            assert(center_v == m.walker(next_h).vertex());
            h_next_set.insert(m.walker(next_h).halfedge());
          }

        }



        // TC: The purpose of this function is the following: We check, whether the next edge curr_h_next is part of face-loop,
        // that another edge in unique_jns is also part of. If curr_h_next is part of this face-loop, which an edge jn_h is also part of
        // then we don't push curr_h_next to the vector unique_jns. However, if the edge curr_h_next is part of a face-loop, which
        // has not been seen before, then we include it in the vector unique_jns.
        for(auto curr_h_next : h_next_set) {

          //cout<<curr_h_next<<endl;

          // TC: Check that curr_h_next does not self-intersect the face-loop of edge h_iter, which is an edge in the face-loop
          // of edge h = curr_h
          if(intersect_face_loop(m, h_iter, curr_h_next)) {
            continue;
            //goto failed;
          }

          visited = 0;
          for (auto jn_h : unique_jns) {
            if(!check_contains(m, curr_h_next, jn_h)) {
                visited = 1;
            }
          }

          if(visited == 1) {
            continue;//goto failed;
          }
          else {
            unique_jns.push_back(curr_h_next);
          }
        }

    }

    return unique_jns;
}

HMesh::HalfEdgeID find_next_feature(Manifold& m, HalfEdgeID h, int pos_flag) {

//    HalfEdgeSet feature_set;
//    feature_set.insert(h);
    std::queue<HMesh::HalfEdgeID> feature_queue;
    feature_queue.push(h);

    std::vector<HMesh::HalfEdgeID> feature_edges, unique_jns, leaf_edges;

    HMesh::HalfEdgeID curr_h;

    while(!feature_queue.empty()) {

        curr_h = feature_queue.front();
        feature_edges = find_face_loop_stack(m, curr_h, pos_flag);
 

        // ----------------------------------------
        // TC: If feature_edges does not contain any edges, because the mesh is empty, 
        // or if that last inserted edge into feature_edges (the one which could be a leaf extrusion) is actually not a leaf
        // extrusion, then we need to run the find_unique_jns function, with the argument being the last edge, where the face-loops
        // where alligned.
        if(feature_edges.size() == 0) {
//            feature_set.erase(curr_h_it);
            feature_queue.pop();
            continue;
        }

        if(check_leaf(m, feature_edges.back(), pos_flag)) {
            leaf_edges.push_back(curr_h);
//            feature_set.erase(curr_h_it);

            feature_queue.pop();
            continue;
        }

        unique_jns = find_unique_jns(m, feature_edges.back(), pos_flag);

        //cout<<"found unique jns"<<endl;

        int contains_flag = 0;
        for (auto h_jn : unique_jns) {
            //cout<<h_jn<<endl;
 /**           contains_flag = 0;
            for(auto h_set : feature_set)
               if(check_contains(m, h_set, h_jn))
                   contains_flag = 1;
            if(contains_flag == 0) **/
            feature_queue.push(h_jn);
        }

//        feature_set.erase(curr_h_it);
        feature_queue.pop();
    }

    leaf_edges = filter_the_leaf_edges(m, leaf_edges, pos_flag);
    

    // ---------------------------------------
    // TC: In the following we go through all the leaf edges, which we have detected in the code above. Each leaf-edge is on a feature, 
    // and for each leaf-edge leaf_h we trace all the aligned face-loops above it and below it. With these face-loops we then compute a cylindiricty 
    // measure, and then we dide that cylindiricty measure with the total number of face-loops. The leaf-edge with the biggest cylindiricty measure is then
    // used and returned as the best_leaf_edge.
    // Karran computes the curr_cyl in a strange way: He devides it by (1 + 10*avg_radius + 10*curr_face_loop.area), but I am not quite sure, why he does that?
    // ---------------------------------------
    //cout<<"Queue traversal done, checking leaves"<<endl;
    FaceLoop curr_face_loop;
    HalfEdgeAttributeVector<int> touched(m.no_halfedges(), 0);
    float curr_cyl;
    float eps = 1e-6;
    float max_cyl = 0;
    std::vector<HMesh::HalfEdgeID> above_stack, below_stack;
    HMesh::HalfEdgeID best_leaf_edge;

    for(auto leaf_h : leaf_edges) {//

        //if(valency(m, m.walker(leaf_h).vertex()) != 4 || valency(m, m.walker(leaf_h).opp().vertex()) != 4)
        //  continue;
        curr_cyl = 0.0;
        //cout<<leaf_h<<endl<<m.walker(leaf_h).opp().halfedge()<<endl;
        if(leaf_h == InvalidHalfEdgeID) {
            continue;
        }
        //curr_cyl = trace_face_loop(m, touched, curr_h).cylindricity;
        above_stack = find_face_loop_stack(m, leaf_h, 0);//;

        below_stack = find_face_loop_stack(m, leaf_h, 1);
 
        for(auto stack_h : above_stack) {
            //cout<<stack_h<<endl;
            HalfEdgeAttributeVector<int> h_touched(m.no_halfedges(), 0);
            curr_face_loop = trace_face_loop(m, h_touched, stack_h);
            double avg_radius = 0.0;
            for(auto loop_h : curr_face_loop.hvec)
              avg_radius += length(curr_face_loop.center - m.pos(m.walker(loop_h).vertex()));
            avg_radius /= curr_face_loop.hvec.size();

            curr_cyl += curr_face_loop.cylindricity/(1 + 10*avg_radius + 10*curr_face_loop.area);
            //curr_cyl += curr_face_loop.perimeter; //1.0/curr_face_loop.cylindricity; //curr_face_loop.perimeter; //curr_face_loop.radius / (1 + curr_face_loop.cylindricity); //curr_face_loop.cylindricity; ///(1 + 10*avg_radius + 10*curr_face_loop.area);

        }
        for(auto stack_h : below_stack) {
        //cout<<stack_h<<endl;
            HalfEdgeAttributeVector<int> h_touched(m.no_halfedges(), 0);
            curr_face_loop = trace_face_loop(m, h_touched, stack_h);
            double avg_radius = 0.0;
            for(auto loop_h : curr_face_loop.hvec)
              avg_radius += length(curr_face_loop.center - m.pos(m.walker(loop_h).vertex()));
            avg_radius /= curr_face_loop.hvec.size();

            curr_cyl += curr_face_loop.cylindricity/(1 + 10*avg_radius + 10*curr_face_loop.area);
            //curr_cyl += curr_face_loop.perimeter; //1.0/curr_face_loop.cylindricity; //curr_face_loop.perimeter;//curr_face_loop.radius / (1 + curr_face_loop.cylindricity); //curr_face_loop.cylindricity; ///(1 + 10*avg_radius + 10*curr_face_loop.area);

        }
        if(abs(curr_cyl) > 1e-6) {
            curr_cyl /= (above_stack.size() + below_stack.size());
        }
        //cout<<curr_cyl<<endl;

        if(curr_cyl > max_cyl) {
            max_cyl = curr_cyl;
            best_leaf_edge = leaf_h;
        }
    }
    return best_leaf_edge;
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
 * Split groups of faces into faces that belong to the base patch and faces that belong to the boundary loops
 * In this way we can handle faces that need to be found on the previous extrusion
 * ----------------------------------------------------------------------- */
std::tuple<std::vector<HMesh::FaceSet>, std::vector<HMesh::FaceSet>, std::vector<std::pair<HMesh::FaceSet, bool>>> split_groups_of_faces(std::vector<HMesh::FaceSet> face_groups, HMesh::FaceSet base_face_set) {
    std::vector<std::pair<HMesh::FaceSet, bool>> split_face_groups;
    std::vector<FaceSet> face_loop_groups;
    std::vector<FaceSet> base_patch_groups;

    // True if it is base-patch faces, and false if it is face_loop faces

    for (auto fg : face_groups) {
        FaceSet base_patch_faces;
        FaceSet face_loop_faces;
        for (auto f : fg) {
            if (base_face_set.find(f) != base_face_set.end()) {
                base_patch_faces.insert(f);
            }
            else {
                face_loop_faces.insert(f);
            }
        }
        if (!base_patch_faces.empty()) {
            split_face_groups.push_back(std::make_pair(base_patch_faces, true));
            base_patch_groups.push_back(base_patch_faces);
        }
        if (!face_loop_faces.empty()) {
            split_face_groups.push_back(std::make_pair(face_loop_faces, false));
            face_loop_groups.push_back(face_loop_faces);
        }
    }

    return {face_loop_groups, base_patch_groups, split_face_groups};
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

void update_extrusion_element_with_bd_v(HMesh::Manifold &m, Extrusion & ext, HMesh::HalfEdgeID bd_h, CGLA::Vec2d bd_vs) {
    extrusion_tree[ext.id].bd_h = bd_h;
    extrusion_tree[ext.id].bd_v = m.walker(bd_h).vertex();
    extrusion_tree[ext.id].bd_v_pos = m.pos(m.walker(bd_h).vertex());
    extrusion_tree[ext.id].bd_vs = bd_vs;

    std::stack<HarmonicMap> transform_stack = ext.hmap_stack;
    auto new_HarmonicMap = transform_stack.top();
    transform_stack.pop();
    new_HarmonicMap.set_bd_h(bd_h);
    new_HarmonicMap.set_bd_v(m.walker(bd_h).vertex());
    new_HarmonicMap.recompute_HarmonicMap(ext.m_state, 
                                                ext.base_face_set, 
                                                m.walker(bd_h).vertex(), ext.id);
    transform_stack.push(new_HarmonicMap);
    extrusion_tree[ext.id].hmap_stack = transform_stack;
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

void compute_new_bd_v_from_partial_prev_ext(Extrusion & ext) {

    // A map to store which extrusions have which faces
    std::map<int, HMesh::FaceSet> contrib_ext_face_sets;

    HMesh::VertexSet boundary_vertices = boundary_verts(ext.m_state, ext.base_face_set);
    HMesh::HalfEdgeSet boundary_edges = boundary_hes(ext.m_state, ext.base_face_set);
    HalfEdgeSet patch_edges = extended_patch_edges(ext.m_state, ext.base_face_set);


    double max_angle = 0.0;
    HMesh::VertexID new_bd_v = HMesh::InvalidVertexID;
    CGLA::Vec2d new_bd_vs;
    int bd_vs_responsible_extrusion = -1;
    std::vector<HMesh::FaceSet> bd_v_face_groups;
    std::map<HMesh::VertexID, CGLA::Vec2d> bd_v_uv_map;

    // Variables needed to ensure that the boundary vertex is from a face-set that comes from a base-patch
    bool is_ref_v_already_set = false;
    double cutoff = 1.5;

    // Edge and face to store the boundary vertex
    HMesh::HalfEdgeID bd_h;
    HMesh::FaceID bd_f;
    HMesh::FaceSet current_ref_v_face_set;
    
    // Find faces from all contributing extrusions
    for (auto contrib_ext_id : ext.contributing_extrusions) {

        //cout << "Investigating contributing extrusion with id: " << contrib_ext_id << endl;

        auto contrib_ext = extrusion_tree[contrib_ext_id];

        HMesh::Manifold m = contrib_ext.m_state_before;
        HMesh::Manifold m_unsmoothed = contrib_ext.m_state_before;

        // Find all the faces that belong to the contributing extrusion
        HMesh::FaceSet base_face_set = contrib_ext.base_face_set;
        HMesh::FaceSet full_face_set;
        for(auto f : base_face_set)  {
            full_face_set.insert(f);
        }

        for(auto f : contrib_ext.stack_faces) {
            full_face_set.insert(f);
        }

        smooth_faceset_lap_solve(m, base_face_set);

        //cout << "Fetching extrusion with id: " << contrib_ext.id << endl;
        stack<HarmonicMap> transform_stack = contrib_ext.hmap_stack;
        auto hmap = transform_stack.top();

        std::map<int, HMesh::FaceSet> ext_face_sets;
        auto extrusion_id = contrib_ext.extrusion_id;
        for(auto f : full_face_set) {
            if (extrusion_id[f] != -1) {
                int curr_ext_id = extrusion_id[f];
                if(ext_face_sets.find(curr_ext_id) != ext_face_sets.end()) {
                    ext_face_sets.find(curr_ext_id)->second.insert(f);
                }
                else {
                    FaceSet curr_fs;
                    curr_fs.insert(f);
                    ext_face_sets.insert(std::make_pair(curr_ext_id, curr_fs));
                }
            }
        }

        // --------------------------
        // Find the correct boundary vertex, because we have selected the one after the extrusion is killed
        // We cannot use contrib_ext.bd_v, because bd_v was found on contrib_ext.m_state and not on contrib_ext.m_state_before
        // and thus bd_v might be different on contrib_ext.m_state_before than on contrib_ext.m_state. Consequently, we need to 
        // use the boundary edge instead.
        // --------------------------
        HMesh::VertexSet temp_bd_verts = boundary_verts(m, base_face_set);
        HMesh::VertexID correct_bd_v;
        for (auto v : temp_bd_verts) {
            circulate_vertex_ccw(m, v, [&] (HMesh::HalfEdgeID h) {
                if (m.walker(h).vertex() == m.walker(contrib_ext.bd_h).vertex()) {
                    correct_bd_v = v;
                }
            });
        }

        bool save_mesh = true;
        std::string extrusion_name = "ext_" + std::to_string(contrib_ext.id);
        //if (contrib_ext.id >= 0 && save_mesh_patches) {
        //    save_mesh = true;
        //}
  
        //std::map<HMesh::VertexID, CGLA::Vec2d> vertex_uv_map = hmap.compute_harmonic_map_with_face_loop(m, extrusion_name, save_mesh, correct_bd_v);
        std::map<HMesh::VertexID, CGLA::Vec2d> vertex_uv_map = hmap.compute_normal_harmonic_map_with_face_loop(m, extrusion_name, save_mesh, correct_bd_v);

        for(auto it = ext_face_sets.begin(); it != ext_face_sets.end(); it++) {
            int next_ext_id = it->first;
            FaceSet next_fs = it->second;

            Extrusion next_ext = extrusion_tree.find(next_ext_id)->second;

            // Only include the faces of the extrusion that matches ext 
            if (next_ext.id == ext.id) {

                auto groups_of_face_sets = find_groups_of_faces(m, next_fs);
                auto [face_loop_groups, base_patch_groups, face_groups] = split_groups_of_faces(groups_of_face_sets, base_face_set);
                
                // If some of the face sets in face_groups come from the base-patch, where they are the result of previous face-loop faces that have been extruded, we need to detect that,
                // because then they should be split up
                // However, we should only split it them up, if we don't need the entire base-patch
                // ------------------------------------------------------------------------------------------
                bool use_entire_base_patch = ext.contributing_extrusions_entire_base_patch.find(contrib_ext.id)->second;
                
                /*
                if (!use_entire_base_patch) {
                    auto face_groups_copy = face_groups; 
                    std::vector<HMesh::FaceSet> new_face_sets;
                    std::set<int> face_groups_to_remove;
                    for (int kk = 0; kk < face_groups.size(); kk++) {
                        auto grouped_faces = group_f_loop_faces(face_groups[kk].first, f_loop_map, f_loop_counter, contrib_ext.id);

                        cout << "Size of grouped_faces is: " << grouped_faces.size() << endl;
                        for (auto new_fs : grouped_faces) {
                            cout << "New face set: " << endl;
                            for (auto f : new_fs.first) {
                                cout << f << " ";
                            }
                            cout << endl;
                        }

                        f_loop_counter += 1;

                        if (grouped_faces.size() > 1) {
                            face_groups_to_remove.insert(kk);
                            for (auto new_fs : grouped_faces) {
                                new_face_sets.push_back(new_fs.first);
                            }
                        }
                    }
                    for (auto kk : face_groups_to_remove) {
                        face_groups_copy.erase(face_groups_copy.begin() + kk);
                    }
                    for (auto new_fs : new_face_sets) {
                        face_groups_copy.push_back(std::make_pair(new_fs, true));
                    }
                    face_groups.clear(); face_groups = face_groups_copy;


                    cout << "face_groups after change: " << endl;
                    for (auto faceset : face_groups) {
                        for (auto f : faceset.first) {
                            cout << f << " ";
                        }
                        cout << endl;
                    }
                    cout << endl;
                }   
                */
                // ------------------------------------------------------------------------------------------

                int index_of_bd_v_responsible_fs = -1;                


                // Only use the base-patch groups to find the new boundary vertex, because since we cut the mesh, then this might be more accurate 
                // than placing on the face-loop. However, if this group is empty, then we need to use the face-loop groups instead.
                //for (auto faceset : ref_v_responsible_face_groups) {
                for (int kk = 0; kk < face_groups.size(); kk++) {

                    auto faceset = face_groups[kk].first;
                    auto inside_base_patch = face_groups[kk].second;

                    if ((is_ref_v_already_set || base_patch_groups.size() > 0) && !inside_base_patch) {
                        continue;
                    }
                    // Change the cutoff
                    if (inside_base_patch) {
                        cutoff = 0.7;
                        is_ref_v_already_set = false;
                    }
                    else {
                        cutoff = 1.5;
                    }
                    // Change the max-angle
                    if (base_patch_groups.size() > 0 && is_ref_v_already_set == false) {
                        max_angle = 0.0;
                    }
                    /*
                    cout << "Before we begin new_bd_v is: " << new_bd_v << endl;
                    cout << "From " << contrib_ext_id << " we find the new ref_v using faceset " << endl;
                    for (auto f : faceset) {
                        cout << f << " ";
                    }
                    cout << endl;
                    */
                    /*
                    cout << "And the boundary_edges are: " << endl;
                    for (auto h : boundary_edges) {
                        cout << h << " ";
                    }
                    cout << endl;
                    */

                    //obj_save("mesh_used_to_compute_new_bd_v_for_ext_" + to_string(ext.id) + "_using_contrib_ext_id_" + to_string(contrib_ext_id) + ".obj", m);

                    // Find the new boundary vertex with the highest curvature
                    auto [tmp1, tmp2, peak_bd_vertices] = compute_new_bd_v_with_high_curvature(m_unsmoothed, 
                            faceset, vertex_uv_map, boundary_edges, cutoff);

                    // Go through all peak_bd_vertices
                    for (auto peak : peak_bd_vertices) {
                        auto bd_v = peak.first;
                        auto angle = peak.second;
   

                        // Do a preliminary check to see if the edge that points to the boundary vertex actually exists
                        HMesh::HalfEdgeID local_bd_h, candidate_bd_h;
                        for (auto h : boundary_hes(m, faceset)) {
                            if (m.walker(h).vertex() == bd_v) {
                                local_bd_h = h;
                            }
                        }
                        
                        circulate_vertex_ccw(ext.m_state, ext.m_state.walker(local_bd_h).vertex(), [&](HMesh::HalfEdgeID h) {
                            auto opp_h = ext.m_state.walker(h).opp().halfedge();
                            if (m.walker(opp_h).vertex() == bd_v && patch_edges.find(opp_h) == patch_edges.end()) {
                                candidate_bd_h = opp_h;
                            }
                        });

                        /*
                        cout << "bd_v != InvalidVertexID is: " << (bd_v != InvalidVertexID) << endl;
                        cout << "candidate_bd_h != HMesh::InvalidHalfEdgeID is: " << (candidate_bd_h != HMesh::InvalidHalfEdgeID) << endl;
                        cout << "angle > max_angle is: " << (angle > max_angle) << endl;
                        cout << "!is_ref_v_already_set is: " << (!is_ref_v_already_set) << endl;
                        */
                        if (bd_v != InvalidVertexID && candidate_bd_h != HMesh::InvalidHalfEdgeID && angle > max_angle && !is_ref_v_already_set) {
                            max_angle = angle;
                            new_bd_v = bd_v;
                            new_bd_vs = vertex_uv_map.find(bd_v)->second;
                            bd_vs_responsible_extrusion = contrib_ext.id;
                            HMesh::HalfEdgeID local_bd_h;
                            index_of_bd_v_responsible_fs = kk;
                            current_ref_v_face_set = faceset;

                            // Select the bd_h as the edge on the face-loop which points towards the bd_v
                            for (auto h : boundary_hes(m, faceset)) {
                                if (m.walker(h).vertex() == new_bd_v) {
                                    local_bd_h = h;
                                }
                            }
                            /*
                            cout << "The maximum angle is now: " << angle << endl;
                            cout << "The new_bd_v is: " << new_bd_v << endl;
                            cout << "local_bd_h is: " << local_bd_h << endl;
                            cout.flush();
                            */

                            circulate_vertex_ccw(ext.m_state, ext.m_state.walker(local_bd_h).vertex(), [&](HMesh::HalfEdgeID h) {
                                auto opp_h = ext.m_state.walker(h).opp().halfedge();
                                if (m.walker(opp_h).vertex() == bd_v && patch_edges.find(opp_h) == patch_edges.end()) {
                                    bd_h = opp_h;
                                }
                            });
                            
                            // The very special case with quadmeshes/bunny.obj, where the next extrusion (extrusion 37) only depends on one previous extrusion
                            // namely extrusion 38, but the extrusion 37 only depends on a subset of the faces of the base-patch of extrusion 38.
                            if (bd_h == HMesh::InvalidHalfEdgeID) {
                                cout << "Problem in 'computing_new_bd_v_from_partial_prev_ext' in decomposing the mesh - The mesh seems to not have a proper face-loop structure" << endl;
                                assert(false);
                            }
                        }
                    }
                    if (inside_base_patch) {
                        is_ref_v_already_set = true;
                    } 
                }

      
 
                std::vector<MatrixXd> fs_loop_group;
                std::vector<HMesh::FaceSet> fs_group;

                for (int kk = 0; kk < face_groups.size(); kk++) {
                    auto faceset = face_groups[kk].first;

                    // --------------------------
                    // Change the curves that make up the face sets
                    // If the face set contains the new_bd_v, then we should change the order of the vertices in the curve
                    // --------------------------
                    vector<VertexID> bd_verts;

                    auto unordered_bd_vertices = boundary_verts(m, faceset);

                    if (unordered_bd_vertices.find(new_bd_v) != unordered_bd_vertices.end()) {
                        bd_verts = find_boundary_vertices(m, faceset, new_bd_v);

                    }
                    else {
                        bd_verts = find_boundary_vertices(m, faceset, *unordered_bd_vertices.begin());
                    }
                    MatrixXd fs_loop;

                    fs_loop.resize(bd_verts.size(), 2);

                    for(int i = 0; i < bd_verts.size(); i++) {

                        if(vertex_uv_map.find(bd_verts[i]) != vertex_uv_map.end()) {

                            Vec2d curr_bd_uv = vertex_uv_map.find(bd_verts[i])->second;

                            fs_loop(i, 0) = curr_bd_uv[0];
                            fs_loop(i, 1) = curr_bd_uv[1];
                        }

                    }
                    //cout << "The fs_loop is: " << fs_loop << endl;
                    fs_loop_group.push_back(fs_loop);
                    fs_group.push_back(faceset);
                }
                // Rearrange the elements in the vector, so index_of_bd_v_responsible_fs becomes the last element in the vector
                if (index_of_bd_v_responsible_fs >= 0) {
                    std::rotate(fs_loop_group.begin() + index_of_bd_v_responsible_fs, fs_loop_group.begin() + index_of_bd_v_responsible_fs + 1, fs_loop_group.end());
                    std::rotate(fs_group.begin() + index_of_bd_v_responsible_fs, fs_group.begin() + index_of_bd_v_responsible_fs + 1, fs_group.end());
                }
                extrusion_tree[contrib_ext.id].next_gen_extrusion_loops[ext.id] = fs_loop_group;
                extrusion_tree[contrib_ext.id].next_gen_extrusion_face_sets[ext.id] = fs_group;

            }

        }
    }

    // --------------------------
    // Update the extrusion element - including changing the bd_vs_responsible_extrusion
    // --------------------------
    // If the extrusion depends on more than one extrusion, then new_bd_v might actually occur in multiple contributing extrusions.
    // This can happen in the following scenario: Extrusion 139 depends on Extrusion 140 and Extrusion 141, but extrusion 140 depends entirely
    // on extrusion 141. If the new_bd_v is then in extrusion 141, it might also appear in Extrusion 140, because when we kill the face-loop in extrusion 140
    // the vertex ID might be moved down, thus making it appear twice in both extrusions. We thus need to store the bd_h instead, as this will be unique. 
    //cout << "Inside the function the extrusion is updated with bd_h: " << bd_h << " and new_bd_vs: " << new_bd_vs << endl;
    update_extrusion_element_with_bd_v(extrusion_tree[ext.id].m_state, ext, bd_h, new_bd_vs);
    extrusion_tree[ext.id].bd_vs_responsible_extrusion = bd_vs_responsible_extrusion;

    extrusion_tree[ext.id].ref_v_face_set = current_ref_v_face_set;
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


HarmonicMap  kill_selected_loop_hmap(Manifold& m, HalfEdgeID h, Extrusion ext) {

    //save_face_loop_collapse_animation(m, h);

    HMesh::FaceSet interior_faces = ext.base_face_set;

    VertexSet patch_vertices;

    FaceLoop above_loop;
    FaceLoop below_loop;
    FaceLoop l;

    HalfEdgeAttributeVector<int> touched(m.no_halfedges(), 0);

    FaceAttributeVector<int> visited(m.no_faces(), 0);

    l = trace_face_loop(m, touched, h);

    HMesh::HalfEdgeID h_above;
    HMesh::HalfEdgeID h_below;

    if(!check_leaf(m, h , 0) && !check_leaf(m, h, 1)) {
        HarmonicMap hmap;
        return hmap;
    }
    //find interior face set

    HMesh::Manifold m_old = HMesh::Manifold(m);

    std::vector<CGLA::Vec3d> old_vertex_positions;

    std::vector<CGLA::Vec3d> old_face_centres;

    std::map<HMesh::FaceID, CGLA::Vec3d> f_centre_map;

    for (auto f : interior_faces) {
        l.interior.insert(f);
        circulate_face_ccw(m, f, [&] (VertexID v) {
          patch_vertices.insert(v);
          old_vertex_positions.push_back(m.pos(v));
        });
        old_face_centres.push_back(centre(m,f));

        f_centre_map.insert({f, centre(m,f)});
    }

    std::vector<CGLA::Vec3d> ext_pos;

    for(auto v: patch_vertices) {
      ext_pos.push_back(m.pos(v));
    }

    if(check_leaf(m, h, 0)) {
        FaceLoop l_rev = trace_face_loop(m, touched, m.walker(l.hvec[0]).next().next().halfedge());
        remove_face_loop(m, l_rev);
    }
    else {
        remove_face_loop(m, l);
    }

    smooth_faceset_lap_solve(m, interior_faces);

    std::vector<CGLA::Vec3d> new_vertex_positions;
    std::vector<HMesh::VertexID> new_vertex_ids;

    for (auto f : interior_faces) {
        l.interior.insert(f);
        circulate_face_ccw(m, f, [&] (VertexID v) {
          new_vertex_positions.push_back(m.pos(v));
          new_vertex_ids.push_back(v);
        });
    }

    // TC: Commented out 5.6.2025

    HarmonicMap hmap(m, interior_faces, ext.origin_face, ext.right_face);

    hmap.input_ext_vertices(new_vertex_ids, old_vertex_positions, f_centre_map);

    return hmap;
}


// TC: Kill the individual features
std::stack<HarmonicMap> kill_individual_feature_hmap(Manifold& m, HalfEdgeID h, int pos_flag, Extrusion ext) {

    std::stack<HarmonicMap> extrusion_stack;


    if(!check_leaf_stack(m, h, 0) && !check_leaf_stack(m, h , 1))
        h = find_next_feature(m, h, pos_flag);

    //cout<<h<<endl;
    auto feature_edges_above = find_face_loop_stack(m, h , 0);
    auto feature_edges_below = find_face_loop_stack(m, h , 1);
    
    // TC: Functions made by TC to only kill individual face loops
    if (feature_edges_above.size() > 1) {
        feature_edges_above.erase(feature_edges_below.begin(), feature_edges_below.end() - 1);
    }
    if (feature_edges_below.size() > 1) {
        feature_edges_below.erase(feature_edges_below.begin()+1, feature_edges_below.end());
    }
    

    std::vector<HMesh::HalfEdgeID> feature_edges;

    for (auto h_iter : feature_edges_above) {
        feature_edges.push_back(h_iter);
    }
    for (auto h_iter : feature_edges_below) {
        feature_edges.push_back(h_iter);
    }
    HMesh::HalfEdgeID h_above, h_below;
    bool work_done;

    work_done = false;
    do {
        work_done = false;
        for(auto h: feature_edges) {
            if(!m.in_use(h))
                continue;

            if(!check_leaf(m,h, 0) && !check_leaf(m, h, 1))
                continue;

            //auto vecs = kill_selected_loop_param_coord(m,h, ext);

            auto vecs = kill_selected_loop_hmap(m, h, ext);


            //auto per_face_transform = face_to_coord_disp_vecs(m, vecs, ext);
            work_done = true;
            extrusion_stack.push(vecs);
            //cout<<vecs.first.first<<endl<<vecs.first.second<<endl;
            break;
        }
    } while (work_done);

    return extrusion_stack;
}


std::stack<HarmonicMap> kill_feature_local_hmap(Manifold& m, HalfEdgeID h, int pos_flag, Extrusion ext) {

    std::stack<HarmonicMap> extrusion_stack;

    if(!check_leaf_stack(m, h, pos_flag))
      return extrusion_stack;

    std::vector<HMesh::HalfEdgeID> feature_edges = find_face_loop_stack(m, h , pos_flag);

    reverse(feature_edges.begin(), feature_edges.end());

    bool work_done;

    work_done = false;
    do {
        work_done = false;
        for(auto h: feature_edges) {
            if(!m.in_use(h)) {
                continue;
            }

            auto vecs = kill_selected_loop_hmap(m, h, ext);

            work_done = true;
            extrusion_stack.push(vecs);
            break;
        }
    } while (work_done);

    return extrusion_stack;
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
    auto extrusion_id_full = FaceAttributeVector<int>(m.no_faces(), -1);
    
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