//
//  Delaunay_triangulate.hpp
//  GEL
//
//  Created by Jakob Andreas Bærentzen on 13/09/2018.
//  Copyright © 2018 J. Andreas Bærentzen. All rights reserved.
//

#ifndef Delaunay_triangulate_hpp
#define Delaunay_triangulate_hpp

#include <vector>
#include <GEL/HMesh/Manifold.h>

namespace HMesh {
    HMesh::Manifold Delaunay_triangulate(const std::vector<CGLA::Vec3d>& pts3d, const CGLA::Vec3d& X_axis = CGLA::Vec3d(1,0,0), const CGLA::Vec3d& Y_axis = CGLA::Vec3d(0,1,0));

    std::tuple<HMesh::Manifold, std::map<int, std::set<int>>, std::vector<CGLA::Vec2d>> constrained_Delaunay_triangulate(const std::vector<CGLA::Vec3d>& pts3d, std::vector<int>& edge_list, const CGLA::Vec3d& X_axis = CGLA::Vec3d(1,0,0), const CGLA::Vec3d& Y_axis = CGLA::Vec3d(0,1,0));
}
#endif /* Delaunay_triangulate_hpp */
