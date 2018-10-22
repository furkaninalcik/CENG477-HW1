#include <iostream>
#include "parser.h"
#include "ppm.h"

typedef unsigned char RGB[3];

using namespace parser;

////////////////we may want to use struct instead of class

class Ray{


    private:
        parser::Vec3f e;
        parser::Vec3f d;

    public:
        
        Ray(){
            printf("empty ray constructor\n");
        }
        Ray(parser::Vec3f origin, parser::Vec3f direction ){

            printf("ray constructor\n");


            e = origin;
            d = direction;
        }
        //~Ray();
    
};


int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    Ray midEyeRay = Ray(scene.cameras[0].position , scene.cameras[0].gaze);

    int n_x = scene.cameras[0].image_width;
    int n_y = scene.cameras[0].image_height;

    int l = scene.cameras[0].near_plane.x;
    int r = scene.cameras[0].near_plane.y;
    int b = scene.cameras[0].near_plane.z;
    int t = scene.cameras[0].near_plane.w;

    printf("width: %d \n"  , n_x);
    printf("height: %d \n" , n_y);


    //Ray RayArray[800*800] ;

    printf("test\n");

    for (int i = 0; i < n_x; ++i)
    {
        for (int j = 0; j < n_y; ++j)
        {
            // s_u = (r – l)(i + 0.5)/n_x
            // s_v = (t – b)(j + 0.5)/n_y


            printf("s_u: %lf \n" , (r - l)*(i + 0.5)/n_x );
            printf("s_v: %lf \n" , (t - b)*(j + 0.5)/n_y );
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
