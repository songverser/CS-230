#include "light.h"
#include "phong_shader.h"
#include "ray.h"
#include "render_world.h"
#include "object.h"

vec3 Phong_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& normal,int recursion_depth) const
{
    vec3 color;
    vec3 colorambi = {0,0,0};
    vec3 colordiffu = {0,0,0};
    vec3 colorspec = {0,0,0};
    // color = ambient + diffuse + specular
    //       = RL + RL*max(n*l,0) + RL*max(v*r,0)^a
    for(int i=0;i<world.lights.size();i++){
        // ray is view direction
        Ray raylight;
        vec3 reflection;// reflection to raylight
        raylight.endpoint = intersection_point;
        vec3 vector_to_light =  (world.lights[i]->position - raylight.endpoint);
        raylight.direction =vector_to_light.normalized();
        // std::cout<<vector_to_light.magnitude()<< " :lightdist intersectdist: "<<world.Closest_Intersection(raylight).dist<<std::endl;
        colorambi = color_ambient*world.ambient_color*world.ambient_intensity;
        if (world.Closest_Intersection(raylight).dist>vector_to_light.magnitude() || !world.enable_shadows) // no obstacles
        {
            vec3 emitlight = world.lights[i]->Emitted_Light(vector_to_light);
            colordiffu += color_diffuse * emitlight*std::max(dot(normal,raylight.direction),double(0)); // n*l

            // wrong?
            reflection = 2*dot(normal,raylight.direction)*normal-raylight.direction; //2(n*l)*n-l
            colorspec += color_specular*emitlight*std::pow(std::max(dot(reflection,-ray.direction),double(0)),specular_power); // (n*l)^p
        }
    }
    
    // color=color_ambient+color_diffuse+color_specular;
    color=colorambi+colordiffu+colorspec;

    TODO; //determine the color
    return color;
}
