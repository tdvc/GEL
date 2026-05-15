#include "feq.h"
#include "gem.h"
#include "HarmonicMap.h"




/* ----------------------------------------------------------------------- *
 * Given a HalfEdgeID h, this function traces the face-loop, which h is part of
 * and finds all faces that make up a patch above the face-loop. The direction of the 
 * HalfEdge indicates, which way is 'up', and therefore which set of faces is above
 * the face-loop
 * ----------------------------------------------------------------------- */
HMesh::FaceSet find_interior_faces(Manifold &m, HMesh::HalfEdgeID h) {


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

    queue<HMesh::FaceID> fq;
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
        auto yellow_vertices = find_yellow_vertices(ext, gen_ext);
        find_yellow_loops(ext, gen_ext, yellow_vertices);

        // TC: Så vi breaker kun, når local_flag er True.
        if(local_flag && !m.in_use(global_h)) {
        //if(local_flag) {
            break;
        }
    } while(true);
    cout<<"returning"<<endl;
    cout.flush();


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