

#include <GEL/HMesh/extrusion.h>
#include <GEL/HMesh/gem.h>

using namespace Eigen;
using namespace std;
using namespace HMesh;

// TC: Purpose: To read a matrix of integers from the extrusion file
Eigen::MatrixXi read_MatrixXi(ifstream& extrusion_file) {
    string line;
    getline(extrusion_file, line);
    istringstream stream(line);
    int rows, cols;
    stream >> rows >> cols;
    Eigen::MatrixXi my_matrix(rows, cols);
    for (int i = 0; i < rows; i++) {
        getline(extrusion_file, line);
        istringstream stream(line);
        std::vector<int> numbers;
        int num;
        while (stream >> num) {
            numbers.push_back(num);
        }
        //stream >> my_matrix(i,0) >> my_matrix(i,1) >> my_matrix(i,2);
        for (size_t j = 0; j < numbers.size(); ++j) {
            my_matrix(i, j) = numbers[j];
        }
    }
    if (rows == 0) {
        getline(extrusion_file, line);
    }
    return my_matrix;
}

// TC: Purpose: To read a matrix of doubles from the extrusion file
Eigen::MatrixXd read_MatrixXd(ifstream& extrusion_file) {
    string line;
    getline(extrusion_file, line);
    istringstream stream(line);
    int rows, cols;
    stream >> rows >> cols;
    Eigen::MatrixXd my_matrix(rows, cols);
    for (int i = 0; i < rows; i++) {
        getline(extrusion_file, line);
        istringstream stream(line);
        std::vector<double> numbers;
        double num;
        while (stream >> num) {
            numbers.push_back(num);
        }
        for (size_t j = 0; j < numbers.size(); ++j) {
            my_matrix(i, j) = numbers[j];
        }
    }
    if (rows == 0) {
        getline(extrusion_file, line);
    }
    return my_matrix;
}

/* ----------------------------------------------------------------------- *
 * Function which loads an extrusion
 * ----------------------------------------------------------------------- */
Extrusion load_extrusion(string file_name) {
    // The loaded extrusion:
    Extrusion ext;
    
    ifstream extrusion_file(file_name);
    if (!extrusion_file.is_open()) {
        cerr << "Failed to open the file!" << endl;
        return ext;
    }
    
    string line;
    // --------------------------
    // origin face
    // --------------------------
    getline(extrusion_file, line);
    getline(extrusion_file, line);
    size_t origin_face = stoi(line);
    auto origin_face_id = HMesh::FaceID(origin_face);
    ext.origin_face = origin_face_id;
    // --------------------------
    // right face
    // --------------------------
    getline(extrusion_file, line);
    getline(extrusion_file, line);
    size_t right_face = stoi(line);
    auto right_face_id = HMesh::FaceID(right_face);
    ext.right_face = right_face_id;
    // --------------------------
    // base_face_set
    // --------------------------
    getline(extrusion_file, line);
    getline(extrusion_file, line);
    stringstream str_base_face(line);
    string number;
    HMesh::FaceSet base_face_set;
    while (std::getline(str_base_face, number, ',')) {
        size_t face = stoul(number);
        auto f_id = HMesh::FaceID(face);
        base_face_set.insert(f_id);
    }
    ext.base_face_set = base_face_set;
    // --------------------------
    // stack faces
    // --------------------------
    getline(extrusion_file, line);
    getline(extrusion_file, line);
    stringstream str_stack_faces(line);
    HMesh::FaceSet stack_faces;
    while (std::getline(str_stack_faces, number, ',')) {
        size_t face = stoul(number);
        auto f_id = HMesh::FaceID(face);
        stack_faces.insert(f_id);
    }
    ext.stack_faces = stack_faces;
    // --------------------------
    // responsible_faces
    // --------------------------
    getline(extrusion_file, line);
    getline(extrusion_file, line);
    stringstream str_responsible_faces(line);
    HMesh::FaceSet responsible_faces;
    while (std::getline(str_responsible_faces, number, ',')) {
        size_t face = stoul(number);
        auto f_id = HMesh::FaceID(face);
        responsible_faces.insert(f_id);
    }
    ext.responsible_faces = responsible_faces;
    // --------------------------
    // bd_v
    // --------------------------
    getline(extrusion_file, line);
    getline(extrusion_file, line);
    size_t vertex = stoul(line);
    auto bd_v = HMesh::VertexID(vertex);
    ext.bd_v = bd_v;
    // --------------------------
    // bd_h
    // --------------------------
    getline(extrusion_file, line);
    getline(extrusion_file, line);
    size_t edge = stoul(line);
    auto bd_h = HMesh::HalfEdgeID(edge);
    ext.bd_h = bd_h;
    // --------------------------
    // bd_v_pos
    // --------------------------
    getline(extrusion_file, line);
    getline(extrusion_file, line);
    double x, y, z;
    char discard;
    istringstream bd_v_pos(line);
    bd_v_pos >> discard >> x >> y >> z >> discard;
    ext.bd_v_pos = CGLA::Vec3d(x, y, z);
    // --------------------------
    // right_vec
    // --------------------------
    getline(extrusion_file, line);
    getline(extrusion_file, line);
    istringstream right_vec(line);
    right_vec >> discard >> x >> y >> z >> discard;
    ext.right_vec = CGLA::Vec3d(x, y, z);
    // --------------------------
    // bd_vs
    // --------------------------
    getline(extrusion_file, line);
    getline(extrusion_file, line);
    istringstream bd_vs(line);
    bd_vs >> discard >> x >> y >> discard;
    ext.bd_vs = CGLA::Vec2d(x, y);
    // --------------------------
    // bd_vs_responsible_extrusion
    // --------------------------
    getline(extrusion_file, line);
    getline(extrusion_file, line);
    ext.bd_vs_responsible_extrusion = stoi(line);
    // --------------------------
    // id
    // --------------------------
    getline(extrusion_file, line);
    getline(extrusion_file, line);
    ext.id = stoi(line);
    // --------------------------
    // scale
    // --------------------------
    getline(extrusion_file, line);
    getline(extrusion_file, line);
    ext.scale = stod(line);
    // --------------------------
    // stack_size
    // --------------------------
    getline(extrusion_file, line);
    getline(extrusion_file, line);
    ext.stack_size = stoi(line);
    // --------------------------
    // hmap_stack / Harmonic Map
    // --------------------------
    for (int kk = 0; kk < ext.stack_size; kk++) {
        HarmonicMap map;
        getline(extrusion_file, line);
        map.set_V(read_MatrixXd(extrusion_file)); // Matrix V

        getline(extrusion_file, line);
        map.set_V_uv(read_MatrixXd(extrusion_file)); // Matrix V_uv

        getline(extrusion_file, line);
        map.set_V_ext(read_MatrixXd(extrusion_file)); // Matrix V_ext

        getline(extrusion_file, line);
        map.set_bnd_uv(read_MatrixXd(extrusion_file)); // Matrix bnd_uv

        getline(extrusion_file, line);
        map.set_N_faces(read_MatrixXd(extrusion_file)); // Matrix N_faces

        getline(extrusion_file, line);
        map.set_F(read_MatrixXi(extrusion_file)); // Matrix F

        getline(extrusion_file, line);
        map.set_F_centres(read_MatrixXd(extrusion_file)); // Matrix F_centres

        getline(extrusion_file, line);
        map.set_bd_F_vertices(read_MatrixXi(extrusion_file)); // Matrix bd_F_vertices

        // patch centre
        getline(extrusion_file, line);
        getline(extrusion_file, line);
        istringstream patch_centre(line);
        patch_centre >> discard >> x >> y >> z >> discard;
        map.set_patch_centre(CGLA::Vec3d(x, y, z));

        // Patch normal
        getline(extrusion_file, line);
        getline(extrusion_file, line);
        istringstream patch_normal(line);
        patch_normal >> discard >> x >> y >> z >> discard;
        map.set_patch_normal(CGLA::Vec3d(x, y, z));

        // Centre face id
        getline(extrusion_file, line);
        getline(extrusion_file, line);
        map.set_centre_face_id(stoi(line));

        // bd_perim
        getline(extrusion_file, line);
        getline(extrusion_file, line);
        map.set_bd_perim(stod(line));

        // patch_area
        getline(extrusion_file, line);
        getline(extrusion_file, line);
        map.set_patch_area(stod(line));

        // bd_v
        getline(extrusion_file, line);
        getline(extrusion_file, line);
        vertex = stoi(line);
        bd_v = HMesh::VertexID(vertex);
        map.set_bd_v(bd_v);

        // bd_h
        getline(extrusion_file, line);
        getline(extrusion_file, line);
        size_t edge = stoi(line);
        auto bd_h = HMesh::HalfEdgeID(edge);
        map.set_bd_h(bd_h);

        // vertex_igl_map
        getline(extrusion_file, line);
        getline(extrusion_file, line);
        int num_elements = stoi(line);
        std::map<HMesh::VertexID, int> vertex_igl_map;
        for (int ii = 0; ii < num_elements; ii++) {
            getline(extrusion_file, line);
            istringstream stream(line);
            size_t v_id;
            int value;
            stream >> v_id >> value;
            vertex_igl_map[HMesh::VertexID(v_id)] = value;
        }
        map.set_vertex_igl_map(vertex_igl_map);

        // igl_vertex_map
        getline(extrusion_file, line);
        getline(extrusion_file, line);
        num_elements = stoi(line);
        std::map<int, HMesh::VertexID> igl_vertex_map;
        for (int ii = 0; ii < num_elements; ii++) {
            getline(extrusion_file, line);
            istringstream stream(line);
            size_t v_id;
            int value;
            stream >> value >> v_id;
            igl_vertex_map[value] = HMesh::VertexID(v_id);
        }
        map.set_igl_vertex_map(igl_vertex_map);

        // face_igl_map
        getline(extrusion_file, line);
        getline(extrusion_file, line);
        num_elements = stoi(line);
        std::map<HMesh::FaceID, int> face_igl_map;
        for (int ii = 0; ii < num_elements; ii++) {
            getline(extrusion_file, line);
            istringstream stream(line);
            size_t face_id;
            int value;
            stream >> face_id >> value;
            face_igl_map[HMesh::FaceID(face_id)] = value;
        }
        map.set_face_igl_map(face_igl_map);

        // igl_face_map
        getline(extrusion_file, line);
        getline(extrusion_file, line);
        num_elements = stoi(line);
        std::map<int, HMesh::FaceID> igl_face_map;
        for (int ii = 0; ii < num_elements; ii++) {
            getline(extrusion_file, line);
            istringstream stream(line);
            size_t face_id;
            int value;
            stream >> value >> face_id;
            igl_face_map[value] = HMesh::FaceID(face_id);
        }
        map.set_igl_face_map(igl_face_map);

        // bd_uvs
        getline(extrusion_file, line);
        getline(extrusion_file, line);
        num_elements = stoi(line);
        std::vector<CGLA::Vec2d> bd_uvs;
        for (int ii = 0; ii < num_elements; ii++) {
            getline(extrusion_file, line);
            istringstream vector2D(line);
            vector2D >> discard >> x >> y >> discard;
            bd_uvs.push_back(CGLA::Vec2d(x, y));
        }
        map.set_bd_uvs(bd_uvs);

        // f_v_map
        getline(extrusion_file, line);
        getline(extrusion_file, line);
        num_elements = stoi(line);
        std::map<HMesh::FaceID, HMesh::VertexID> f_v_map;
        for (int ii = 0; ii < num_elements; ii++) {
            getline(extrusion_file, line);
            istringstream stream(line);
            size_t face_id;
            size_t vertex_id;
            stream >> face_id >> vertex_id;
            f_v_map[HMesh::FaceID(face_id)] = HMesh::VertexID(vertex_id);
        }
        map.set_f_v_map(f_v_map);

        // Insert the extrusion into the stack
        ext.hmap_stack.push(map);
    }
    extrusion_file.close();
    //cout << "Done loading extrusion from file: " << file_name << endl;
    return ext;
}


/* ----------------------------------------------------------------------- *
 * Function used for initializing a Generic Extrusion
 * ----------------------------------------------------------------------- */
Generic_Extrusion::Generic_Extrusion() {
    
    if (!load(("gen_ext.obj"),m)) {
        cout << "Could not load the generic extrusion" << endl;
    };
    
    // Get the positions
    pos = m.positions_attribute_vector();
    
    
    for (auto f : m.faces()) {
        curr_ext_faces.insert(f);
    }
    
    for(auto v : all_verts(m, curr_ext_faces)) {
        if(vertex_uv_map.find(v) == vertex_uv_map.end()) {
            CGLA::Vec2d uv = CGLA::Vec2d(m.pos(v)[0], m.pos(v)[1]);
            vertex_uv_map.insert(std::make_pair(v, uv));
      }
    }

    // Initialize the kDtree
    for (auto it = vertex_uv_map.begin(); it != vertex_uv_map.end(); it++) {
        uv_tree.insert(it->second, it->first.index);
    }
    uv_tree.build();

    // Find rim vertices (boundary vertices on the entire extrusion)
    // Note: There are no particular order, because we don't have a 'start vertex'
    rim_vertices = boundary_verts(m, curr_ext_faces);
}