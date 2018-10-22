#include <iostream>
#include "parser.h"
#include "ppm.h"

typedef unsigned char RGB[3];

//using namespace parser;

////////////////we may want to use struct instead of class

struct Vec3f // Is ": parser::Vec3f" necesssary?
{

    float x, y, z;

    Vec3f(){
        printf("\n empty constructor \n");
    }

    Vec3f(parser::Vec3f vector) : x(vector.x), y(vector.y), z(vector.z) {

    }
    Vec3f(float x, float y, float z) : x(x), y(y), z(z) {}
    
    Vec3f operator * (float d) const { 
        printf("Distance: %lf\n", d );
        printf("MULTIPLICATION\n");
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
    
};

class Ray{


    private:
        parser::Vec3f e;
        parser::Vec3f d;

    public:
        
        Ray(){
            //printf("empty ray constructor\n");
        }
        Ray(parser::Vec3f origin, parser::Vec3f direction ){

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



void intersection(Ray r){

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

    Vec3f s;
    
    float s_u,s_v;

    int n_x = scene.cameras[0].image_width;
    int n_y = scene.cameras[0].image_height;

    float distance = scene.cameras[0].near_distance; 

    int l = scene.cameras[0].near_plane.x;
    int r = scene.cameras[0].near_plane.y;
    int b = scene.cameras[0].near_plane.z;
    int t = scene.cameras[0].near_plane.w;

    printf("width: %d \n"  , n_x);
    printf("height: %d \n" , n_y);


    // slide -> http://saksagan.ceng.metu.edu.tr/courses/ceng477/files/pdf/week_02.pdf ------- page 13/49

    //find the coordanates of the point "q" (the point at the top-left of image plane )


    Vec3f m =  (-w) * distance ;  // m is the intersection point of the gazeRay and the image plane

    Vec3f q = m + u*((float)l) + v*((float)t); // BE CAREFULL about type casting -> "(float)l" 

    

    /////////// VECTOR OPERATIONS TEST ///////////
        /*
        Vec3f testVector = Vec3f(1.2,-2.4,11.1);

        int testDistance = 2;

        Vec3f result = (-testVector) * testDistance;

        printf("x = %lf y = %lf y = %lf \n", result.x, result.y, result.z  );

        */
        
    /////////// VECTOR OPERATIONS TEST ///////////





    //find the coordanates of the point "s" (the point we look through in ray tracing)



    printf("test\n");

    for (int i = 0; i < n_x; ++i)
    {
        for (int j = 0; j < n_y; ++j)
        {
            s_u = (r - l)*(i + 0.5)/n_x;
            s_v = (t - b)*(j + 0.5)/n_y;


            //printf("s_u: %lf \n" , (r - l)*(i + 0.5)/n_x );
            //printf("s_v: %lf \n" , (t - b)*(j + 0.5)/n_y );

            s = q + u * s_u - v * s_v;

            //Ray eyeRay = Ray();

            intersection(Ray());


        }
    }

    printf("TOP LEFT  :\n"    );
    printf("s_u: %lf \n" , (r - l)*(0 + 0.5)/n_x );
    printf("s_v: %lf \n" , (t - b)*(0 + 0.5)/n_y );

    printf("BOTTOM RIGHT  :\n"    );
    printf("s_u: %lf \n" , (r - l)*(n_x-1 + 0.5)/n_x );
    printf("s_v: %lf \n" , (t - b)*(n_y-1 + 0.5)/n_y );

    printf("NUMBERS: \n");
    printf("l: %d \n"  , l);
    printf("r: %d \n"  , r);
    printf("b: %d \n"  , b);
    printf("t: %d \n"  , t);
    printf("n_x: %d \n"  , n_x);
    printf("n_y: %d \n"  , n_y);


    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.

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





    for (int i = 0; i < scene.cameras.size(); ++i)
    {
        std::cout << scene.cameras[i].image_name << std::endl;

        const char* filename =  scene.cameras[i].image_name.c_str();

        //const char* filepath = std::strcat("TestOutputs/", filename);

        write_ppm(filename, image, width, height);

        
    }
    


}
