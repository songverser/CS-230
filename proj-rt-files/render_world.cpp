#include "render_world.h"
#include "flat_shader.h"
#include "object.h"
#include "light.h"
#include "ray.h"

extern bool disable_hierarchy;

Render_World::Render_World()
    :background_shader(0),ambient_intensity(0),enable_shadows(true),
    recursion_depth_limit(3)
{}

Render_World::~Render_World()
{
    delete background_shader;
    for(size_t i=0;i<objects.size();i++) delete objects[i];
    for(size_t i=0;i<lights.size();i++) delete lights[i];
}

// Find and return the Hit structure for the closest intersection.  Be careful
// to ensure that hit.dist>=small_t.
Hit Render_World::Closest_Intersection(const Ray& ray)
{
    TODO;
    Hit ans;
    ans.object = NULL;
    ans.dist = std::numeric_limits<double>::max();
    for(int i = 0; i<objects.size();i++){
        Hit tans;
        tans = objects[i]->Intersection(ray,-1); // part < 0 test for all part ?
        tans.object = objects[i];
        if(tans.dist == 0){
            continue;
        }
        else if(tans.dist < ans.dist){
            ans = tans;
        }
    }
    return ans;
}

// set up the initial view ray and call
void Render_World::Render_Pixel(const ivec2& pixel_index)
{
    TODO; // set up the initial view ray here
    Ray ray; // turn image position to world position
    vec2 uv;
    uv[0] = (pixel_index-camera.number_pixels/2)[0]+0.5;
    uv[1] = (pixel_index-camera.number_pixels/2)[1]+0.5;
    uv *= camera.pixel_size;
    vec3 indexposition = camera.film_position+double(uv[0])*camera.horizontal_vector+double(uv[1])*camera.vertical_vector;
    ray.endpoint = camera.position;
    ray.direction = (indexposition - camera.position).normalized();
    // std::cout<<ray.endpoint << ray.direction <<std::endl; //test
    vec3 color=Cast_Ray(ray,recursion_depth_limit); // failed
    camera.Set_Pixel(pixel_index,Pixel_Color(color)); 
}

void Render_World::Render()
{
    if(!disable_hierarchy)
        Initialize_Hierarchy();

    for(int j=0;j<camera.number_pixels[1];j++)
        for(int i=0;i<camera.number_pixels[0];i++)
            Render_Pixel(ivec2(i,j));
}

// cast ray and return the color of the closest intersected surface point,
// or the background color if there is no object intersection
vec3 Render_World::Cast_Ray(const Ray& ray,int recursion_depth)
{
    vec3 color;
    TODO; // determine the color here
    Hit h = Closest_Intersection(ray);
    if(h.object == NULL){
        color = {0,0,0};
        return color;
    }
    vec3 intepoint = ray.endpoint + h.dist*ray.direction;
    vec3 normal = h.object->Normal(intepoint,1); //part unknown
    color = h.object->material_shader->Shade_Surface(ray, intepoint, normal, recursion_depth);
    return color;
}

void Render_World::Initialize_Hierarchy()
{
    TODO; // Fill in hierarchy.entries; there should be one entry for
    // each part of each object.
    // hierarchy.entries = ;
    return;
    hierarchy.Reorder_Entries();
    hierarchy.Build_Tree();
}
