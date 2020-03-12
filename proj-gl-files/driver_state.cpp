#include "driver_state.h"
#include <cstring>
#include <math.h>
#include <algorithm>
#include <numeric>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=0;
    state.image_depth=0;
    std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
    state.image_color = new pixel[state.image_width*state.image_height];
    state.image_depth = new float[state.image_width*state.image_height];
    for(int i=0; i<state.image_width*state.image_height; i++){
        state.image_depth[i]=1000;
        state.image_color[i]=make_pixel(0*255,0*255,0*255);
    }

}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    std::cout<<"TODO: implement rendering."<<std::endl;
    // state type
    if(type == render_type::triangle){
        // std::cout<<state.num_vertices<<std::endl;
        for(int k=0;k<state.num_vertices/3;k++){
            data_geometry* in[3];
            data_vertex* dv = new data_vertex[3];
            for(int i=0; i<3; i++){
                dv[i].data = new float[state.floats_per_vertex];
                for(int j=0; j<state.floats_per_vertex; j++){
                    dv[i].data[j] = state.vertex_data[(i+k*3)*state.floats_per_vertex+j];
                    // std::cout<<"dv dataj: "<<dv[i].data[j]<<std::endl;
                }
                in[i]=new data_geometry;
                in[i]->data = dv[i].data;
                state.vertex_shader(dv[i], *in[i], state.uniform_data);
                // std::cout<<in[i]->data[0]<<in[i]->data[1]<<std::endl;
            }
            for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                    in[i]->data[j] = in[i]->gl_Position[j]/in[i]->gl_Position[3];
                }
            }
            const data_geometry* in2[3] = {in[0],in[1],in[2]};
            rasterize_triangle(state, in2);
            // clip_triangle(state,in2,0);
            delete[] dv[0].data; delete dv[1].data; delete dv[2].data;
            delete[]* in;
        }
    }
    if(type == render_type::fan){
        data_geometry* in[3];
        data_vertex* dv = new data_vertex[3];
        // origin of fan
        dv[0].data = new float[state.floats_per_vertex];
        for(int j=0; j<state.floats_per_vertex; j++){
            dv[0].data[j] = state.vertex_data[j];
        }
        in[0] = new data_geometry;
        in[0]->data = dv[0].data;
        state.vertex_shader(dv[0], *in[0], state.uniform_data);
        for(int j=0;j<3;j++){
            in[0]->data[j] = in[0]->gl_Position[j]/in[0]->gl_Position[3];
        }
        // fans
        for(int k=0;k<state.num_vertices-2; k++){
            // std::cout<<"fan count"<<k<<std::endl;
            for(int i=1; i<3; i++){
                dv[i].data = new float[state.floats_per_vertex];
                for(int j=0; j<state.floats_per_vertex; j++){
                    dv[i].data[j] = state.vertex_data[(k+i)*state.floats_per_vertex+j];
                    // std::cout<<"dv dataj: "<<dv[i].data[j]<<std::endl;
                }
                in[i]=new data_geometry;
                in[i]->data = dv[i].data;
                state.vertex_shader(dv[i], *in[i], state.uniform_data);
                // std::cout<<in[i]->data[0]<<in[i]->data[1]<<std::endl;
            }
            for(int i=1;i<3;i++){
                for(int j=0;j<3;j++){
                    in[i]->data[j] = in[i]->gl_Position[j]/in[i]->gl_Position[3];
                }
            }
            const data_geometry* in2[3] = {in[0],in[1],in[2]};
            // clip_triangle(state,in2,0);
            rasterize_triangle(state, in2);
        }        
    }
    if(type == render_type::strip){
        data_geometry* in[3];
        data_vertex* dv = new data_vertex[3];
        // strip
        for(int k=0;k<state.num_vertices-2; k++){
            for(int i=0; i<3; i++){
                dv[i].data = new float[state.floats_per_vertex];
                for(int j=0; j<state.floats_per_vertex; j++){
                    dv[i].data[j] = state.vertex_data[(k+i)*state.floats_per_vertex+j];
                    // std::cout<<"dv dataj: "<<dv[i].data[j]<<std::endl;
                }
                in[i]=new data_geometry;
                in[i]->data = dv[i].data;
                state.vertex_shader(dv[i], *in[i], state.uniform_data);
                // std::cout<<in[i]->data[0]<<in[i]->data[1]<<std::endl;
            }
            for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                    in[i]->data[j] = in[i]->gl_Position[j]/in[i]->gl_Position[3];
                }
            }
            const data_geometry* in2[3] = {in[0],in[1],in[2]};
            // clip_triangle(state,in2,0);
            rasterize_triangle(state, in2);
            // delete[] dv[0].data,dv[1].data,dv[2].data,dv;
            // delete[]* in;
        }
    }
    if(type == render_type::indexed){
        data_geometry* in[3];
        data_vertex* dv = new data_vertex[3];
        // strip
        for(int k=0;k<state.num_triangles; k++){
            // std::cout<<"tris: "<<k<<std::endl;
            for(int i=0; i<3; i++){
                int ind = state.index_data[k*3+i];
                dv[i].data = new float[state.floats_per_vertex];
                for(int j=0; j<state.floats_per_vertex; j++){
                    dv[i].data[j] = state.vertex_data[ind*state.floats_per_vertex+j];
                    // std::cout<<"dv dataj: "<<dv[i].data[j]<<std::endl;
                }
                in[i]=new data_geometry;
                in[i]->data = dv[i].data;
                state.vertex_shader(dv[i], *in[i], state.uniform_data);
                // std::cout<<in[i]->data[0]<<in[i]->data[1]<<std::endl;
            }
            for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                    in[i]->data[j] = in[i]->gl_Position[j]/in[i]->gl_Position[3];
                }
            }
            const data_geometry* in2[3] = {in[0],in[1],in[2]};
            // clip_triangle(state,in2,0);
            rasterize_triangle(state, in2);
            // delete[] dv[0].data,dv[1].data,dv[2].data,dv;
            // delete[]* in;
        }
    }
}

float rayplaneintec(vec3 e, vec3 d, vec3 p0, vec3 n){
    float k = dot((p0-e),n)/dot(d,n);
    return fabs(k);
}
float area(double x1, double y1, double x2, double y2, double x3, double y3){
    return (float)fabs((x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2);
}

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    // left up right down fordword backword
    // left
    data_geometry* ina = new data_geometry;
    data_geometry* inb = new data_geometry;
    int inoout[3],sum=0;
    vec3 e,dir,d,p0,n;
    // find outlier
    if(face==0){
        for(int i=0;i<3;i++){
            if(in[i]->data[0]<-1){
                inoout[i] = 0;
            } else {
                inoout[i] = 1;
            }
        }
        sum = std::accumulate(inoout,inoout+3,0);
        p0 = vec3(-1,0,0);
        n = vec3(1,0,0);
    } else if (face==1){
        for(int i=0;i<3;i++){
            if(in[i]->data[1]>1){
                inoout[i] = 0;
            } else {
                inoout[i] = 1;
            }
        }
        sum = std::accumulate(inoout,inoout+3,0);
        p0 = vec3(0,1,0);
        n = vec3(0,-1,0);
    }
    else if(face==2){
        for(int i=0;i<3;i++){
            if(in[i]->data[0]>1){
                inoout[i] = 0;
            } else {
                inoout[i] = 1;
            }
        }
        sum = std::accumulate(inoout,inoout+3,0);
        p0 = vec3(1,0,0);
        n = vec3(-1,0,0);
    } else if (face==3){
        for(int i=0;i<3;i++){
            if(in[i]->data[1]<-1){
                inoout[i] = 0;
            } else {
                inoout[i] = 1;
            }
        }
        sum = std::accumulate(inoout,inoout+3,0);
        p0 = vec3(0,-1,0);
        n = vec3(0,1,0);
    } else if(face == 4){
        for(int i=0;i<3;i++){
            if(in[i]->data[2]<-1){
                inoout[i] = 0;
            } else {
                inoout[i] = 1;
            }
        }
        sum = std::accumulate(inoout,inoout+3,0);
        p0 = vec3(0,0,-1);
        n = vec3(0,0,1);
    }else if (face == 5) {
        for(int i=0;i<3;i++){
            if(in[i]->data[2]>1){
                inoout[i] = 0;
            } else {
                inoout[i] = 1;
            }
        }
        sum = std::accumulate(inoout,inoout+3,0);
        p0 = vec3(0,0,1);
        n = vec3(0,0,-1);
    }

    // do the cut
    if(sum == 2){
        for(int i=0;i<3;i++){
            if(inoout[i]==0){
                float* data0 = in[(i)%3]->data;
                float* data1 = in[(i+1)%3]->data;
                float* data2 = in[(i+2)%3]->data;

                vec3 colorintec,color0,color1,color2;
                color0 = vec3(data0[3],data0[4],data0[5]);
                color1 = vec3(data1[3],data1[4],data1[5]);
                color2 = vec3(data2[3],data2[4],data2[5]);

                e = vec3(data0[0],data0[1],data0[2]);
                dir = vec3(data1[0],data1[1],data1[2]);
                d = dir-e;
                d = d.normalized();
                float length = rayplaneintec(e, d, p0, n);
                vec3 intec = e+d*length;
                vec3 colora = (1-length/(dir-e).magnitude())*color0+(length/(dir-e).magnitude())*color1;
                if(state.interp_rules[3]==interp_type::smooth){
                    float alpha = 1-length/(dir-e).magnitude(), beta = length/(dir-e).magnitude();
                    float wa = in[i]->gl_Position[3], wb = in[(i+1)%3]->gl_Position[3];
                    float alphanew = (alpha/wa)/(alpha/wa+beta/wb), betanew = (beta/wb)/(alpha/wa+beta/wb);
                    colora = alphanew*color0+betanew*color1;
                }
                ina->data = new float[6];
                ina->data[0]=intec[0]; ina->data[1]=intec[1]; ina->data[2]=intec[2];ina->data[3]=colora[0];ina->data[4]=colora[1];ina->data[5]=colora[2];
                ina->gl_Position[3] = (1-length/(dir-e).magnitude())*in[i]->gl_Position[3]+length/(dir-e).magnitude()*in[(i+1)%3]->gl_Position[3];

                e = vec3(data0[0],data0[1],data0[2]);
                dir = vec3(data2[0],data2[1],data2[2]);
                d = dir-e;
                d = d.normalized();
                length = rayplaneintec(e, d, p0, n);
                intec = e+d*length;
                vec3 colorb = (1-length/(dir-e).magnitude())*color0+length/(dir-e).magnitude()*color2;
                if(state.interp_rules[3]==interp_type::smooth){
                    float alpha = 1-length/(dir-e).magnitude(), beta = length/(dir-e).magnitude();
                    float wa = in[i]->gl_Position[3], wb = in[(i+2)%3]->gl_Position[3];
                    float alphanew = (alpha/wa)/(alpha/wa+beta/wb), betanew = (beta/wb)/(alpha/wa+beta/wb);
                    colorb = alphanew*color0+betanew*color2;
                }
                inb->data = new float[6];
                inb->data[0]=intec[0]; inb->data[1]=intec[1]; inb->data[2]=intec[2];inb->data[3]=colorb[0];inb->data[4]=colorb[1];inb->data[5]=colorb[2];
                inb->gl_Position[3] = (1-length/(dir-e).magnitude())*in[i]->gl_Position[3]+length/(dir-e).magnitude()*in[(i+2)%3]->gl_Position[3];
                

                // std::cout<<inb->gl_Position[3]<<" "<<in[i]->gl_Position[3]<<" "<<in[(i+2)%3]->gl_Position[3]<<std::endl;
                // std::cout<<std::endl<<ina->data[0] <<" " << ina->data[1]<<" "<<ina->data[2]<<" "<<ina->data[3]<<" "<<ina->data[4]<<" "<<ina->data[5]<<" "<<std::endl;
                // std::cout<<inb->data[0] <<" " << inb->data[1]<<" "<<inb->data[2]<<" "<<inb->data[3]<<" "<<inb->data[4]<<" "<<inb->data[5]<<" face "<<face<<std::endl;

                const data_geometry* inans1[3] = {ina,in[(i+2)%3],in[(i+1)%3]};
                const data_geometry* inans2[3] = {inb,ina,in[(i+2)%3]};
                clip_triangle(state,inans1,face+1);
                clip_triangle(state,inans2,face+1);
                delete[] ina->data;
                delete[] inb->data;
            }
        }
    }else if (sum == 1){
        for(int i=0;i<3;i++){
            if(inoout[i]==1){
                float* data0 = in[(i)%3]->data;
                float* data1 = in[(i+1)%3]->data;
                float* data2 = in[(i+2)%3]->data;

                vec3 colorintec,color0,color1,color2;
                color0 = vec3(data0[3],data0[4],data0[5]);
                color1 = vec3(data1[3],data1[4],data1[5]);
                color2 = vec3(data2[3],data2[4],data2[5]);

                e = vec3(data0[0],data0[1],data0[2]);
                dir = vec3(data1[0],data1[1],data1[2]);
                d = dir-e;
                d = d.normalized();
                float length = rayplaneintec(e, d, p0, n);
                vec3 intec = e+d*length;
                vec3 colora = (1-length/(dir-e).magnitude())*color0+length/(dir-e).magnitude()*color1;
                if(state.interp_rules[3]==interp_type::smooth){
                    float alpha = 1-length/(dir-e).magnitude(), beta = length/(dir-e).magnitude();
                    float wa = in[i]->gl_Position[3], wb = in[(i+1)%3]->gl_Position[3];
                    float alphanew = (alpha/wa)/(alpha/wa+beta/wb), betanew = (beta/wb)/(alpha/wa+beta/wb);
                    colora = alphanew*color0+betanew*color1;
                }
                ina->data = new float[6];
                ina->data[0]=intec[0]; ina->data[1]=intec[1]; ina->data[2]=intec[2];ina->data[3]=colora[0];ina->data[4]=colora[1];ina->data[5]=colora[2];
                ina->gl_Position[3] = (1-length/(dir-e).magnitude())*in[i]->gl_Position[3]+length/(dir-e).magnitude()*in[(i+1)%3]->gl_Position[3];

                e = vec3(data0[0],data0[1],data0[2]);
                dir = vec3(data2[0],data2[1],data2[2]);
                d = dir-e;
                d = d.normalized();
                length = rayplaneintec(e, d, p0, n);
                intec = e+d*length;
                vec3 colorb = (1-length/(dir-e).magnitude())*color0+length/(dir-e).magnitude()*color2;
                if(state.interp_rules[3]==interp_type::smooth){
                    float alpha = 1-length/(dir-e).magnitude(), beta = length/(dir-e).magnitude();
                    float wa = in[i]->gl_Position[3], wb = in[(i+2)%3]->gl_Position[3];
                    float alphanew = (alpha/wa)/(alpha/wa+beta/wb), betanew = (beta/wb)/(alpha/wa+beta/wb);
                    colorb = alphanew*color0+betanew*color2;
                }
                inb->data = new float[6];
                inb->data[0]=intec[0]; inb->data[1]=intec[1]; inb->data[2]=intec[2];inb->data[3]=colorb[0];inb->data[4]=colorb[1];inb->data[5]=colorb[2];
                inb->gl_Position[3] = (1-length/(dir-e).magnitude())*in[i]->gl_Position[3]+length/(dir-e).magnitude()*in[(i+2)%3]->gl_Position[3];

                // std::cout<<std::endl<<ina->data[0] <<" " << ina->data[1]<<" "<<ina->data[2]<<" "<<ina->data[3]<<" "<<ina->data[4]<<" "<<ina->data[5]<<" "<<std::endl;
                // std::cout<<inb->data[0] <<" " << inb->data[1]<<" "<<inb->data[2]<<" "<<inb->data[3]<<" "<<inb->data[4]<<" "<<inb->data[5]<<" face "<<face<<std::endl;

                const data_geometry* inans1[3] = {ina,inb,in[i]};
                clip_triangle(state,inans1,face+1);
                
                delete[] ina->data;
                delete[] inb->data;
            }
        }
    }else if(sum == 0){
        return;
    }else if(sum == 3){
        clip_triangle(state,in,face+1);
    }
    delete[] ina;
    delete[] inb;

    // std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    // clip_triangle(state,in,face+1);
}


float min3(float a, float b, float c){
    if(a<=b && a<=c) return a;
    else if(b<=a&& b<=c) return b;
    else return c;
}
float max3(float a, float b, float c){
    if(a>=b && a>=c) return a;
    else if(b>=a&& b>=c) return b;
    else return c;
}


// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{

    // std::cout<<in[0]->data[0]<<" "<<in[0]->data[1]<<" "<<in[0]->data[2]<<std::endl;
    // std::cout<<in[1]->data[0]<<" "<<in[1]->data[1]<<" "<<in[1]->data[2]<<std::endl;
    // std::cout<<in[2]->data[0]<<" "<<in[2]->data[1]<<" "<<in[2]->data[2]<<std::endl;
    // std::cout<<"TODO: implement rasterization"<<std::endl;
    data_fragment df;// = new data_fragment;
    data_output* dataout = new data_output[state.image_width*state.image_height];

    
    float xmin = (min3(in[0]->data[0],in[1]->data[0],in[2]->data[0])+1)/2*(float)state.image_width;
    float ymin = (min3(in[0]->data[1],in[1]->data[1],in[2]->data[1])+1)/2*(float)state.image_height;
    float xmax = (max3(in[0]->data[0],in[1]->data[0],in[2]->data[0])+1)/2*(float)state.image_width;
    float ymax = (max3(in[0]->data[1],in[1]->data[1],in[2]->data[1])+1)/2*(float)state.image_height;
    // std::cout<<xmin<<" "<<xmax<<" "<<ymin<<" "<<ymax<<std::endl;
    // one by one fragment shader
    for(float i=std::max(0,(int)floor(ymin)); i<std::min(state.image_height,(int)ceil(ymax)); i++){
        for(float j=std::max(0,(int)floor(xmin)); j<std::min(state.image_width,(int)ceil(xmax)); j++){
            int position = i*state.image_width+j;
            // transform
            float h,w;
            float kh = 2/(float)state.image_height;
            float kw = 2/(float)state.image_width;
            h = kh*i-1;
            w = kw*j-1;

            float areab=area(in[0]->data[0],in[0]->data[1],in[1]->data[0],in[1]->data[1],in[2]->data[0],in[2]->data[1]);
            float area1=area(w,h,in[1]->data[0],in[1]->data[1],in[2]->data[0],in[2]->data[1]);
            float area2=area(w,h,in[0]->data[0],in[0]->data[1],in[2]->data[0],in[2]->data[1]);
            float area3=area(w,h,in[0]->data[0],in[0]->data[1],in[1]->data[0],in[1]->data[1]);
            float alpha = area1/areab; float beta = area2/areab; float gamma = area3/areab;
            float depth = alpha*in[0]->data[2]+beta*in[1]->data[2]+gamma*in[2]->data[2];
            // df.data = new float[state.floats_per_vertex];
            // std::cout<<areab<<" "<<area1<<" "<<area2<<" "<<area3<<" "<<depth<<std::endl;
            if(round((area1+area2+area3)*10000) == round(areab*10000)){
                if(state.floats_per_vertex==3 || state.interp_rules[3]==interp_type::flat){
                    df.data[0]=(in[0]->data[3]);df.data[1]=(in[0]->data[4]);df.data[2]=(in[0]->data[5]);
                    state.fragment_shader(df, dataout[position], state.uniform_data);
                    // std::cout<<dataout[position].output_color[0]<<dataout[position].output_color[1]<<dataout[position].output_color[2]<<std::endl;
                } else if(state.interp_rules[3]==interp_type::noperspective){
                    vec3 colora = vec3(in[0]->data[3],in[0]->data[4],in[0]->data[5]);
                    vec3 colorb = vec3(in[1]->data[3],in[1]->data[4],in[1]->data[5]);
                    vec3 colorc = vec3(in[2]->data[3],in[2]->data[4],in[2]->data[5]);
                    
                    df.data = alpha*colora+beta*colorb+gamma*colorc;
                    state.fragment_shader(df, dataout[position], state.uniform_data);
                    // std::cout<<dataout[position].output_color[0]<<dataout[position].output_color[1]<<dataout[position].output_color[2]<<std::endl;
                } else if(state.interp_rules[3]==interp_type::smooth){
                    vec3 colora = vec3(in[0]->data[3],in[0]->data[4],in[0]->data[5]);
                    vec3 colorb = vec3(in[1]->data[3],in[1]->data[4],in[1]->data[5]);
                    vec3 colorc = vec3(in[2]->data[3],in[2]->data[4],in[2]->data[5]);
                    float wa = in[0]->gl_Position[3], wb=in[1]->gl_Position[3], wc=in[2]->gl_Position[3];
                    float alphanew = (alpha/wa)/(alpha/wa+beta/wb+gamma/wc);
                    float betanew = (beta/wb)/(alpha/wa+beta/wb+gamma/wc);
                    float gammanew = (gamma/wc)/(alpha/wa+beta/wb+gamma/wc);
                    
                    df.data = alphanew*colora+betanew*colorb+gammanew*colorc;
                    state.fragment_shader(df, dataout[position], state.uniform_data);
                    // std::cout<<dataout[position].output_color[0]<<dataout[position].output_color[1]<<dataout[position].output_color[2]<<std::endl;
                }
                if(depth<state.image_depth[position] && depth>=-1 && depth <=1){
                    state.image_color[position] = make_pixel(dataout[position].output_color[0]*255,dataout[position].output_color[1]*255,dataout[position].output_color[2]*255);
                    state.image_depth[position]=depth;
                    // std::cout<<in[0]->data[0]<<" "<<in[0]->data[1]<<" "<<in[0]->data[2]<<std::endl;
                }
            }
        }
    }
    delete[] dataout;
}

