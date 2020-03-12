#include "plane.h"
#include "ray.h"
#include <cfloat>
#include <limits>

// Intersect with the half space defined by the plane.  The plane's normal
// points outside.  If the ray starts on the "inside" side of the plane, be sure
// to record a hit with t=0 as the first entry in hits.
Hit Plane::Intersection(const Ray& ray, int part) const
{
    TODO;
    Hit hit;
    hit.dist=0;
    double length;
    length = dot((x1-ray.endpoint),normal)/dot(ray.direction,normal);
    if(length > small_t){
        hit.dist=length;
        hit.part = part;
    }
    return hit;
}

vec3 Plane::Normal(const vec3& point, int part) const
{
    return normal;
}

// There is not a good answer for the bounding box of an infinite object.
// The safe thing to do is to return a box that contains everything.
Box Plane::Bounding_Box(int part) const
{
    Box b;
    b.hi.fill(std::numeric_limits<double>::infinity());
    b.lo=-b.hi;
    // b.hi.fill(0);
    // b.lo.fill(0);
    
    return b;
}
