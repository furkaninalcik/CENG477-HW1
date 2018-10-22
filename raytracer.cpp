#include <iostream>
#include "parser.h"
#include "ppm.h"


#include <cmath> // MAKE SURE THAT IT IS ALLOWED



typedef unsigned char RGB[3];

//using namespace parser;

////////////////we may want to use struct instead of class

struct Vec3f // Is ": parser::Vec3f" necesssary?
{

    float x, y, z;

    Vec3f(){
        //printf("\n empty constructor \n");
    }

    Vec3f(parser::Vec3f vector) : x(vector.x), y(vector.y), z(vector.z) {

    }
    Vec3f(float x, float y, float z) : x(x), y(y), z(z) {}
    
    Vec3f operator * (float d) const { 
        //printf("Distance: %lf\n", d );
        //printf("MULTIPLICATION\n");
        return Vec3f(x*d, y*d, z*d); 
    }

    Vec3f operator + (Vec3f v) const { 
        return Vec3f(x+v.x, y+v.y, z+v.z); 
    }

    Vec3f operator - (Vec3f v) const { 
        return Vec3f(x-v.x, y-v.y, z-v.z); 
    }

    Vec3f operator = (parser::Vec3f vector) const { 
        printf("Assignment! \n" );

        return Vec3f(vector); 
    }

    Vec3f operator-() const {
        Vec3f v;
        v.x = -x;
        v.y = -y;
        v.z = -z;
        return v;
   }

    Vec3f normalize() const {
        double norm = sqrt(x*x + y*y + z*z);
        return Vec3f(x/norm,y/norm,z/norm);
    }
    
};

class Ray{


    private:
        

    public:
        
        Vec3f e;
        Vec3f d;

        Ray(){
            //printf("empty ray constructor\n");
        }
        Ray(Vec3f origin, Vec3f direction ){

            //printf("ray constructor\n");


            e = origin;
            d = direction;
        }
        Vec3f RayVect(float t){

            Vec3f v2 = Vec3f(d.x*t, d.y*t, d.z*t);

            Vec3f result = Vec3f(e.x + v2.x, e.y + v2.y, e.z + v2.z );

            return result;
        }
       
        //~Ray();
    
};

Vec3f crossProduct(Vec3f u , Vec3f v ){

    Vec3f result = Vec3f( (u.y*v.z - u.z*v.y) , (u.z*v.x - u.x*v.z) , (u.x*v.y - u.y*v.x) );

    return result;

}

float dotProduct(Vec3f u , Vec3f v ){

    return (u.x*v.x + u.y*v.y + u.z*v.z);

}

bool intersection(Ray ray, parser::Sphere sphere, Vec3f center , float& t){

    Vec3f e = ray.e;
    Vec3f d = ray.d;

    float r = sphere.radius; // radius of the sphere

    float a = dotProduct(d,d);           // a is A in the equation -> At^2 + Bt + C = 0 // 
    float b = 2*dotProduct(d,e-center);       // b is B in the equation -> At^2 + Bt + C = 0 // 
    float c = dotProduct(e-center,e-center) - r*r; // c is C in the equation -> At^2 + Bt + C = 0 // 

    float discriminant = b*b - 4*a*c;

    if (discriminant < 0.005) // 
    {
        return false;
    }
    else{
        float x0 = -b - sqrt(discriminant); // one of the real roots of the equation
        float x1 = -b + sqrt(discriminant); // one of the real roots of the equation
        t = (x0 < x1) ? x0 : x1;
        return true;        
    }

    //Vec3f c = sphere.vertex_data[scene.center_vertex_id]; // center of the sphere


}



int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    Ray gazeRay = Ray(scene.cameras[0].position , scene.cameras[0].gaze); // the eye ray which is perpendicular to the image plane

    Vec3f e = scene.cameras[0].position; // camera position, the origin of the rays we trace

    Vec3f w = scene.cameras[0].gaze; // camera gaze vector in xyz coordinates
    Vec3f v = scene.cameras[0].up; // camera gaze vector in xyz coordinates
    Vec3f u = crossProduct(v,w); 

    printf("u vector: %lf , %lf , %lf\n" , u.x , u.y , u.z );

    Vec3f s;
    
    float s_u,s_v;

    int n_x = scene.cameras[0].image_width;
    int n_y = scene.cameras[0].image_height;

    float distance = scene.cameras[0].near_distance; 

    float l = scene.cameras[0].near_plane.x;
    float r = scene.cameras[0].near_plane.y;
    float b = scene.cameras[0].near_plane.z;
    float t = scene.cameras[0].near_plane.w;

    printf("width: %d \n"  , n_x);
    printf("height: %d \n" , n_y);
    printf("l: %lf , r: %lf , b: %lf , t: %lf  \n", l, r, b, t  );


    // slide -> http://saksagan.ceng.metu.edu.tr/courses/ceng477/files/pdf/week_02.pdf ------- page 13/49

    //find the coordanates of the point "q" (the point at the top-left of image plane )


    Vec3f m = e + (-w) * distance ;  // m is the intersection point of the gazeRay and the image plane

    Vec3f q = m + u*l + v*t; //  

    

    /////////// VECTOR OPERATIONS TEST ///////////
        
        Vec3f testVector  = Vec3f(1.2,-2.4,11.1);
        Vec3f testVector2 = Vec3f(0.2,0.4,1.1);

        
        int testDistance = 2;

        Vec3f result1 = (-testVector) * testDistance;
        Vec3f result = (testVector) -  testVector2;

        printf("x = %lf y = %lf z = %lf \n", result1.x, result1.y, result1.z  );
        printf("x = %lf y = %lf z = %lf \n", result.x, result.y, result.z  );

        
        
    /////////// VECTOR OPERATIONS TEST ///////////





    //find the coordanates of the point "s" (the point we look through in ray tracing)


    Ray eyeRay ;

    printf("test\n");


    // IMAGE CREATION FOR TESTING

    int width = scene.cameras[0].image_width;
    int height = scene.cameras[0].image_height;
    int columnWidth = width / 8;

    unsigned char* image = new unsigned char [width * height * 3];
    /*
    int i = 0;
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int colIdx = x / columnWidth;
            image[i++] = BAR_COLOR[colIdx][0];
            image[i++] = BAR_COLOR[colIdx][1];
            image[i++] = BAR_COLOR[colIdx][2];
        }
    }
    */

    // IMAGE CREATION FOR TESTING

    int index = 0;

    for (int i = 0; i < n_x; ++i)
    {
        for (int j = 0; j < n_y; ++j)
        {
            s_u = (r - l)*(i + 0.5)/n_x;
            s_v = (t - b)*(j + 0.5)/n_y;

            //printf(" s_u : %lf\n", s_u  );
            //printf(" s_v : %lf\n", s_v  );


            //printf("s_u: %lf \n" , (r - l)*(i + 0.5)/n_x );
            //printf("s_v: %lf \n" , (t - b)*(j + 0.5)/n_y );

            s = q + (u * s_u) - (v * s_v);


            eyeRay = Ray(e, s-e);


            std::vector<parser::Mesh>     meshes    = scene.meshes;
            std::vector<parser::Triangle> triangles = scene.triangles;
            std::vector<parser::Sphere>   spheres   = scene.spheres;


            float t;

            if (intersection(eyeRay, spheres[0], scene.vertex_data[spheres[0].center_vertex_id-1] ,t )){
                image[index++] = 255;
                image[index++] = 255;
                image[index++] = 255;
            }
            else{
                image[index++] = 20;
                image[index++] = 20;
                image[index++] = 20;
            } 

        }
    }

    printf("TOP LEFT  :\n"    );
    printf("s_u: %lf \n" , (r - l)*(0 + 0.5)/n_x );
    printf("s_v: %lf \n" , (t - b)*(0 + 0.5)/n_y );

    printf("BOTTOM RIGHT  :\n"    );
    printf("s_u: %lf \n" , (r - l)*(n_x-1 + 0.5)/n_x );
    printf("s_v: %lf \n" , (t - b)*(n_y-1 + 0.5)/n_y );

    printf("NUMBERS: \n");
    printf("l: %lf \n"  , l);
    printf("r: %lf \n"  , r);
    printf("b: %lf \n"  , b);
    printf("t: %lf \n"  , t);
    printf("n_x: %d \n"  , n_x);
    printf("n_y: %d \n"  , n_y);


    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    /*
    const RGB BAR_COLOR[8] =
    {
        { 255, 255, 255 },  // 100% White
        { 255, 255,   0 },  // Yellow
        {   0, 255, 255 },  // Cyan
        {   0, 255,   0 },  // Green
        { 255,   0, 255 },  // Magenta
        { 255,   0,   0 },  // Red
        {   0,   0, 255 },  // Blue
        {   0,   0,   0 },  // Black
    };

    int width = 640, height = 480;
    int columnWidth = width / 8;

    unsigned char* image = new unsigned char [width * height * 3];

    int i = 0;
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int colIdx = x / columnWidth;
            image[i++] = BAR_COLOR[colIdx][0];
            image[i++] = BAR_COLOR[colIdx][1];
            image[i++] = BAR_COLOR[colIdx][2];
        }
    }


    */


    for (int i = 0; i < scene.cameras.size(); ++i)
    {
        std::cout << scene.cameras[i].image_name << std::endl;

        const char* filename =  scene.cameras[i].image_name.c_str();

        //const char* filepath = std::strcat("TestOutputs/", filename);

        write_ppm(filename, image, width, height);

        
    }
    


}
