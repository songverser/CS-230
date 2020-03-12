#include "render_world.h"
#include "flat_shader.h"
#include "object.h"
#include "light.h"
#include "ray.h"
#include <queue>
#include <stack>

#include "simplemath.h"
#include <chrono>

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
    // do box intersection with hirarchy tree here
    //part1 bfs
    int flag = 0;
    if(hierarchy.tree.size()){
        std::queue<int> que;
        que.push(0);
        while(que.size()>0){ 
            int cur = que.front();
                            
            que.pop(); // true
            if(vequal(hierarchy.tree[cur].hi, hierarchy.tree[cur].lo)){
                continue;
            }
            if(hierarchy.tree[cur].Intersection(ray)){

                if(cur*2+1 >= hierarchy.tree.size()){ // at leaf
                    int level = gettreelevel(cur);
                    int ind = cur-(pow(2,level-1)-1);
                    Hit tans;
                    tans = hierarchy.entries[ind].obj->Intersection(ray, hierarchy.entries[ind].part);
                    tans.object = hierarchy.entries[ind].obj;
                    if(tans.dist == 0){ // no intersection, no hit will <small_t
                        continue;
                    } else if(tans.dist < ans.dist){
                        ans = tans;
                    }
                } else{
                    que.push(cur*2+1);
                    que.push(cur*2+2);
                }
            }
        }
    }
    //part2
    for(int i = 0; i<objouthirarchy.size();i++){
        Hit tans;
        tans = objouthirarchy[i]->Intersection(ray,1); // part later
        tans.object = objouthirarchy[i];
        if(tans.dist == 0){ // because of floating point error, using 0 to represent not hit, no hit will <small_t
            continue;
        } else if(tans.dist < ans.dist){
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
    vec3 color=Cast_Ray(ray,recursion_depth_limit); // failed
    camera.Set_Pixel(pixel_index,Pixel_Color(color)); 
}

void Render_World::Render()
{
    auto start = std::chrono::high_resolution_clock::now();

    if(!disable_hierarchy)
        Initialize_Hierarchy();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( end - start ).count();
    std::cout <<"Duration init: "<< duration<<std::endl;

    start = std::chrono::high_resolution_clock::now();

    for(int j=0;j<camera.number_pixels[1];j++)
        for(int i=0;i<camera.number_pixels[0];i++)
            Render_Pixel(ivec2(i,j));
    
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>( end - start ).count();
    std::cout <<"Duration render: "<< duration<<std::endl;
}

// cast ray and return the color of the closest intersected surface point,
// or the background color if there is no object intersection
vec3 Render_World::Cast_Ray(const Ray& ray,int recursion_depth)
{
    vec3 color;
    TODO; // determine the color here


    // auto start = std::chrono::high_resolution_clock::now();

    Hit h = Closest_Intersection(ray);

    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>( end - start ).count();
    // std::cout <<"Duration one inte: "<< duration<<std::endl;

    // start = std::chrono::high_resolution_clock::now();

    if(h.object == NULL){
        //using background shader
        if(background_shader){
            color = background_shader->Shade_Surface(ray, {0,0,0}, {0,0,0}, recursion_depth);
            return color;
        }
        color = {0,0,0};
        return color;
    }
    vec3 intepoint = ray.endpoint + h.dist*ray.direction;
    vec3 normal = h.object->Normal(intepoint,h.part); //part unknown
    // changing normal direction
    if(dot(normal,ray.direction)>0){
        std::cout<<"normal in wrong direction"<<std::endl;
        normal = -normal;
    }
    normal = normal.normalized();//maybe not normalized
    color = h.object->material_shader->Shade_Surface(ray, intepoint, normal, recursion_depth);
    
    // end = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::microseconds>( end - start ).count();
    // std::cout <<"Duration one rest: "<< duration<<std::endl;
    return color;

}

void Render_World::Initialize_Hierarchy()
{
    TODO; // Fill in hierarchy.entries; there should be one entry for
    // each part of each object.
    // return; // test
    for(int i=0;i<objects.size();i++){
        for(int j=0;j<objects[i]->number_parts;j++){
            Entry tentry;
            tentry.obj = objects[i];
            tentry.part = j;
            tentry.box = objects[i]->Bounding_Box(j);
            // std::cout<<"magnitude box: "<< tentry.box.lo<<" "<< tentry.box.hi<<" magnitude "<<((tentry.box.hi-tentry.box.lo).magnitude()<1000000)<<std::endl;
            if((tentry.box.hi-tentry.box.lo).magnitude()<1000000 && (tentry.box.hi-tentry.box.lo).magnitude()>0){ 
                // inf*inf ?
                hierarchy.entries.push_back(tentry);
            } 
            else{
                objouthirarchy.push_back(objects[i]);
            }
        }
    }
    // finish inserting entries
    
    // std::cout<< hierarchy.entries.size()<<" sizesss "<<objouthirarchy.size() << std::endl;
    hierarchy.Reorder_Entries();
    hierarchy.Build_Tree();
    // std::cout<<"hierarchy finish"<<std::endl;
}
