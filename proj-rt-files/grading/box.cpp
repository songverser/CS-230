#include <limits>
#include "box.h"
#include "mesh.h"

// Return whether the ray intersects this box.
bool Box::Intersection(const Ray& ray) const
{
    TODO;
    // using triangle;
    //    3   2
    //  7   6
    //    0   1
    //  4   5
    // Mesh mesh;
    // mesh.vertices.push_back(lo);
    // mesh.vertices.push_back({lo[0],hi[1],lo[2]});
    // mesh.vertices.push_back({lo[0],hi[1],hi[2]});
    // mesh.vertices.push_back({lo[0],lo[1],hi[2]});
    // mesh.vertices.push_back({hi[0],lo[1],lo[2]});
    // mesh.vertices.push_back({hi[0],hi[1],lo[2]});
    // mesh.vertices.push_back(hi);
    // mesh.vertices.push_back({hi[0],lo[1],hi[2]});
    // mesh.triangles.push_back({0,1,3});
    // mesh.triangles.push_back({1,2,3});
    // mesh.triangles.push_back({0,1,4});
    // mesh.triangles.push_back({1,4,5});
    // mesh.triangles.push_back({0,3,4});
    // mesh.triangles.push_back({3,4,7});
    // mesh.triangles.push_back({4,5,7});
    // mesh.triangles.push_back({5,6,7});
    // mesh.triangles.push_back({1,5,6});
    // mesh.triangles.push_back({1,2,6});
    // mesh.triangles.push_back({2,3,6});
    // mesh.triangles.push_back({3,6,7});
    // double dist;
    // for(int i=0; i<12; i++){
    //     if(mesh.Intersect_Triangle(ray,i,dist)){
    //         return true;
    //     }
    // }

	// there are three tmin and tmax along 3 axis. at least some part should interset each other 
	double tmin = -100000000;
	double tmax = 100000000;
	
	for (int i = 0; i < 3; i++)
	{
		// if along one direction dir=0
		if (std::abs(ray.direction[i]) < small_t)
		{
			if (ray.endpoint[i] < lo[i] || ray.endpoint[i] > hi[i])
				return false;
		}
		else
		{
			const float ood = 1.0f / ray.direction[i];
			// t1 to be lo, t2 to be hi
			float t1 = (lo[i] - ray.endpoint[i]) * ood;
			float t2 = (hi[i] - ray.endpoint[i]) * ood;
			if (t1 > t2) { float tmp = t1; t1 = t2; t2 = tmp; }
			
			if (t1 > tmin) tmin = t1;
			if (t2 < tmax) tmax = t2;
 
			// not intersect
			if (tmin > tmax) return false;
		}
	}
	
	return true;
}

// Compute the smallest box that contains both *this and bb.
Box Box::Union(const Box& bb) const
{
    Box box;
    TODO;
    box.lo = lo;
    box.hi = hi;
    if(vequal(bb.lo, bb.hi)){
        return box;
    } else if(vequal(box.lo, box.hi)){
        return bb;
    }
    for(int j=0;j<3;j++){
        if(box.lo[j]>bb.lo[j]){
            box.lo[j] = bb.lo[j];
        }
        if(box.hi[j] < bb.hi[j]){
            box.hi[j] = bb.hi[j];
        }
    }

    return box;
}

// Enlarge this box (if necessary) so that pt also lies inside it.
void Box::Include_Point(const vec3& pt)
{
    TODO;
    for(int j=0;j<3;j++){
        if(pt[j]<lo[j]){
            lo[j] = pt[j];
        }
        if(pt[j] > hi[j]){
            hi[j] = pt[j];
        }
    }
}

// Create a box to which points can be correctly added using Include_Point.
void Box::Make_Empty()
{
    lo.fill(std::numeric_limits<double>::infinity());
    hi=-lo;
}
