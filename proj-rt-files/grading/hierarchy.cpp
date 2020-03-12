#include <algorithm>
#include "hierarchy.h"

#include <chrono>


bool cmp(Entry aa, Entry bb){
    Box a;
    Box b;
    a = aa.box;
    b = bb.box;
    // if((a.lo[0]+a.lo[1]+a.lo[2])-(b.lo[0]+b.lo[1]+b.lo[2])<0){ //longest one
    // vec3 cut = a.lo-b.lo;
    // cut = cut*cut*cut;
    // if(cut[0]+cut[1]+cut[2]<0){
    if((a.lo[1]-b.lo[1]<0)){ //fastest
        return true;
    }else{
        return false;
    }
}

// Reorder the entries vector so that adjacent entries tend to be nearby.
// You may want to implement box.cpp first.
void Hierarchy::Reorder_Entries()
{
    if(!entries.size()) return;
    // later lo 1st hi 2nd

    auto start = std::chrono::high_resolution_clock::now();

    std::sort(entries.begin(),entries.end(),cmp);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( end - start ).count();
    std::cout <<"Duration: "<< duration<<std::endl;

    TODO;
}

// Populate tree from entries.
void Hierarchy::Build_Tree()
{
    if(!entries.size()) return;
    int si = entries.size();
    int i=0;
    while(pow(2,i)<si){
        i++;
    }
    int nsqr = pow(2,i);
    for(i=0;i<si;i++){
        tree.push_back(entries[i].box); // insert entry to leaf of tree
    }
    for(i=si;i<nsqr;i++){
        Box tbox;
        tbox.lo = {0,0,0};
        tbox.hi = {0,0,0};
        tree.push_back(tbox);
    }
    nsqr=nsqr/2;
    while(nsqr>=1){
        std::vector<Box> tboxes;
        for(i=0;i<nsqr;i++){
            tboxes.push_back(tree[i*2].Union(tree[i*2+1]));
        }
        std::vector<Box> AB = tboxes;
        AB.insert(AB.end(), tree.begin(), tree.end());
        tree = AB;
        nsqr = nsqr/2;
    }
    // for(i=0;i<tree.size();i++){
    //     std::cout<<"    hi: "<<tree[i].hi<<" lo: "<<tree[i].lo;
    // }

    TODO;
}

// Return a list of candidates (indices into the entries list) whose
// bounding boxes intersect the ray.
void Hierarchy::Intersection_Candidates(const Ray& ray, std::vector<int>& candidates) const
{
    TODO;
}
