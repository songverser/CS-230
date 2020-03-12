#include "sphere.h"
#include "ray.h"
#include "simplemath.h"

// Determine if the ray intersects with the sphere
Hit Sphere::Intersection(const Ray& ray, int part) const
{
    TODO;
    // (x-o)^2=r^2   a+kv=x
    // (kv+a-o)^2=r^2
    // ray-> vec3 endpoint direction
    // sphere-> center radius
    vec3 d = ray.direction;
    vec3 e = ray.endpoint;
    double a = dot(d, d);
    double b = 2*dot(d, (e-center));
    double c = dot(e-center, e-center)-radius*radius;
    // std::cout<<"direction: "<<d<<" endpoint: "<<e<<" center: "<<center<<" radius: "<<radius<<std::endl;
    sol solu;
    Hit ans;
    ans.dist=0;
    solu = QuadEq(a,b,c);
    if(solu.solvable == true){
        // ans.object = self;
        ans.dist = minpositive(solu.solution[0],solu.solution[1]); // if both <= small_t set to 0
        ans.part = part;
    }
    return ans;
}

vec3 Sphere::Normal(const vec3& point, int part) const
{
    vec3 normal;
    TODO; // compute the normal direction
    normal = (point-center).normalized();
    return normal;
}

Box Sphere::Bounding_Box(int part) const
{
    Box box;
    box.lo = {center[0]-radius,center[1]-radius,center[2]-radius};
    box.hi = {center[0]+radius,center[1]+radius,center[2]+radius};
    TODO; // calculate bounding box
    return box;
}
