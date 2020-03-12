#include "reflective_shader.h"
#include "ray.h"
#include "render_world.h"

vec3 Reflective_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& normal,int recursion_depth) const
{
    vec3 color;
    //phong*(1-reflectivity) + reflecthit*reflectivity

    color += shader->Shade_Surface(ray,intersection_point,normal,recursion_depth-1)*(1-reflectivity);//from this surface
    if(recursion_depth>1){
        Ray reflectray;
        reflectray.direction = (2*dot(normal,-ray.direction)*normal+ray.direction).normalized();
        reflectray.endpoint = intersection_point;
        Hit hit = world.Closest_Intersection(reflectray);
        if (hit.object != NULL) // has obstacles
        {
            vec3 refintepoint = hit.dist*reflectray.direction+reflectray.endpoint;
            vec3 hitnormal = hit.object->Normal(refintepoint,1);//part ?
            color += hit.object->material_shader->Shade_Surface(reflectray, refintepoint, hitnormal, recursion_depth-1)*reflectivity;
        }
    }


    TODO; // determine the color
    return color;
}
