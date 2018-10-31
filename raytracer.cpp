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
        float norm = sqrt(x*x + y*y + z*z);
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

Vec3f clamp(Vec3f vector) {
  Vec3f v ;
  v.x = (vector.x > 255) ? 255 : (vector.x < 0) ? 0 : vector.x;
  v.y = (vector.y > 255) ? 255 : (vector.y < 0) ? 0 : vector.y;
  v.z = (vector.z > 255) ? 255 : (vector.z < 0) ? 0 : vector.z;
  return v;
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
        float x0 = (-b - sqrt(discriminant))/(2*a); // one of the real roots of the equation
        float x1 = (-b + sqrt(discriminant))/(2*a); // one of the real roots of the equation
        t = (x0 < x1) ? x0 : x1;
        //printf("t1 %lf \n", x0 );
        //printf("t2 %lf \n", x1 );
        return true;        
    }

    //Vec3f c = sphere.vertex_data[scene.center_vertex_id]; // center of the sphere
}



bool intersection(Ray ray, parser::Face face, parser::Scene scene,  float& t, Vec3f& surfaceNormal){

    Vec3f e = ray.e; // origin 
    Vec3f d = ray.d; // direction

    Vec3f p ; // the ray-plane intersection point (may or may not be inside the triangle) 

    float gama, beta; // variables for barycentric coordinates


    Vec3f v1 = scene.vertex_data[face.v0_id - 1];
    Vec3f v2 = scene.vertex_data[face.v1_id - 1];
    Vec3f v3 = scene.vertex_data[face.v2_id - 1];
    /*
    printf("VERTEX 1 : %lf , %lf  , %lf \n" , v1.x, v1.y , v1.z );
    printf("VERTEX 2 : %lf , %lf  , %lf \n" , v2.x, v2.y , v2.z );
    printf("VERTEX 3 : %lf , %lf  , %lf \n" , v3.x, v3.y , v3.z );
    */




    // calculating plane normal


    Vec3f normalVector = crossProduct( v3-v2 , v2-v1);  // BE CAREFULL ABOUT THE ORDER OF THE VERTICES
    surfaceNormal = normalVector; // TO BE USED BY SHADING PART OF THE CODE

    if (dotProduct(normalVector , d)  < 0.000001) // if plane and ray are parallel 
    {
        return false;
    }

    t = (dotProduct((v1 - e),normalVector))/(dotProduct(d,normalVector)); // calculating t to find the ray-plane intersection point "p"


    //printf("t : %lf \n" , t);

    p = e + d * t;


    //printf("TEST1\n");

    /*
    if (t <= 0.000001) // t_min
    {
        return false;
    }
    */

    //printf("TEST2\n");

    /////////////////////////////////////////////

    //calculating the barycentric coordanates
    

    /*

    https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates


    // Compute barycentric coordinates (u, v, w) for
    // point p with respect to triangle (a, b, c)
    void Barycentric(Point p, Point a, Point b, Point c, float &u, float &v, float &w)
    {
        Vector v0 = b - a, v1 = c - a, v2 = p - a;
        float d00 = Dot(v0, v0);
        float d01 = Dot(v0, v1);
        float d11 = Dot(v1, v1);
        float d20 = Dot(v2, v0);
        float d21 = Dot(v2, v1);
        float denom = d00 * d11 - d01 * d01;
        v = (d11 * d20 - d01 * d21) / denom;
        w = (d00 * d21 - d01 * d20) / denom;
        u = 1.0f - v - w;
    }

    */


    //a = v1 
    //b = v2 
    //c = v3 
    //v0 = v_21 
    //v1 = v_31 
    //v2 = v_p1 

    Vec3f v_21 = v2-v1;
    Vec3f v_31 = v3-v1;
    Vec3f v_p1 = p-v1;

    float p1 = dotProduct(v_21, v_21);
    float p2 = dotProduct(v_21, v_31);
    float p3 = dotProduct(v_31, v_31);
    float p4 = dotProduct(v_p1, v_21);
    float p5 = dotProduct(v_p1, v_31);


    float den = p1*p3 - p2*p2; // denominator

    gama = (p3*p4 - p2*p5) / den; // GAMA OR BETA ???

    //printf("GAMA : %lf \n", gama);

    if (gama < 0 || gama > 1 )
    {
        return false;
    }

    //printf("TEST3\n");


    beta = (p1*p5 - p2*p4) / den; // BETA OR GAMA ???

    if (beta < 0 || beta > 1-gama)
    {
        return false;
    }

    //printf("TEST4\n");



    return true;
}





int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    Ray gazeRay = Ray(scene.cameras[0].position , scene.cameras[0].gaze); // the eye ray which is perpendicular to the image plane

    Vec3f e = scene.cameras[0].position; // camera position, the origin of the rays we trace

    Vec3f w = scene.cameras[0].gaze; // camera gaze vector in xyz coordinates
    Vec3f v = scene.cameras[0].up; // camera up vector in xyz coordinates
    Vec3f u = crossProduct(v,-w); 

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


    Vec3f m = e + (w) * distance ;  // m is the intersection point of the gazeRay and the image plane

    Vec3f q = m + u*l + v*t; //  

    

    /////////// VECTOR OPERATIONS TEST ///////////
      
    /*
        
        Vec3f testVector  = Vec3f(1.2,-2.4,11.1);
        Vec3f testVector2 = Vec3f(0.2,0.4,1.1);

        
        int testDistance = 2;

        Vec3f result1 = (-testVector) * testDistance;
        Vec3f result = (testVector) -  testVector2;

        printf("x = %lf y = %lf z = %lf \n", result1.x, result1.y, result1.z  );
        printf("x = %lf y = %lf z = %lf \n", result.x, result.y, result.z  );

    */
        
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

    Vec3f surfaceNormal; // "intersection" function will assign this variable 


    for (int i = 0; i < n_x; ++i)
    {
        for (int j = 0; j < n_y; ++j)
        {
            s_u = (r - l)*(j + 0.5)/n_x;
            s_v = (t - b)*(i + 0.5)/n_y;

            //printf(" s_u : %lf\n", s_u  );
            //printf(" s_v : %lf\n", s_v  );


            //printf("s_u: %lf \n" , (r - l)*(i + 0.5)/n_x );
            //printf("s_v: %lf \n" , (t - b)*(j + 0.5)/n_y );

            s = q + (u * s_u) - (v * s_v);


            eyeRay = Ray(e, (s-e).normalize());


            std::vector<parser::Mesh>     meshes    = scene.meshes;
            std::vector<parser::Triangle> triangles = scene.triangles;
            std::vector<parser::Sphere>   spheres   = scene.spheres;


            float t;

            bool sphereIntersection = false;
            bool triangleIntersection = false;
            bool faceIntersection = false;


            Vec3f lightPosition  = scene.point_lights[0].position; // for testing 
            Vec3f lightIntensity = scene.point_lights[0].intensity; // for testing 

            for (int i = 0; i < spheres.size(); ++i)
            {
                Vec3f center = scene.vertex_data[spheres[i].center_vertex_id-1]; // center of the sphere 
                if (intersection(eyeRay, spheres[i], center ,t )){


                    Vec3f pointOnTheSphere  = eyeRay.e + eyeRay.d*t; 
                    Vec3f sphereSurfaceNormal = (pointOnTheSphere - center) * (1.0 / spheres[i].radius);





                    Vec3f vectorToLight = lightPosition - pointOnTheSphere ; 

                    float lightDistance = sqrt(dotProduct(vectorToLight,vectorToLight));

                    float cosTheta = dotProduct(vectorToLight.normalize(), sphereSurfaceNormal.normalize());

                    cosTheta = (cosTheta < 0) ? 0 : cosTheta;
                    

                    Vec3f diffuseShadingParams = scene.materials[spheres[i].material_id-1].diffuse; // for RGB values -> between 0 and 1

                    Vec3f irradiance = lightIntensity * (1.0/(lightDistance*lightDistance));


                    float diffuseShadingRed   = diffuseShadingParams.x * cosTheta * irradiance.x; 
                    float diffuseShadingGreen = diffuseShadingParams.y * cosTheta * irradiance.y; 
                    float diffuseShadingBlue  = diffuseShadingParams.z * cosTheta * irradiance.z; 

                    Vec3f diffuseShading = Vec3f(diffuseShadingRed,diffuseShadingGreen,diffuseShadingBlue);


                    Vec3f halfWayVector = ((-eyeRay.d).normalize() + vectorToLight.normalize()).normalize();

                    float cosAlpha = dotProduct(halfWayVector.normalize(), sphereSurfaceNormal.normalize()); // for specular shading

                    cosAlpha = (cosAlpha < 0) ? 0 : cosAlpha;


                    Vec3f specularShadingParams = scene.materials[spheres[i].material_id-1].specular; // for RGB values -> between 0 and 1
                    float phong_exponent = scene.materials[spheres[i].material_id-1].phong_exponent; // for RGB values -> between 0 and 1
                    float cosAlphaWithPhong = pow(cosAlpha,phong_exponent); 
                    //printf("Specular : %lf %lf %lf  \n", specularShadingParams.x, specularShadingParams.y, specularShadingParams.z   );


                    float specularShadingRed   = specularShadingParams.x * cosAlphaWithPhong * irradiance.x; 
                    float specularShadingGreen = specularShadingParams.y * cosAlphaWithPhong * irradiance.y; 
                    float specularShadingBlue  = specularShadingParams.z * cosAlphaWithPhong * irradiance.z; 

                    Vec3f specularShading = Vec3f(specularShadingRed,specularShadingGreen,specularShadingBlue);



                    Vec3f diffuseAndSpecular = clamp(diffuseShading+specularShading);

                    image[index++] = diffuseAndSpecular.x;
                    image[index++] = diffuseAndSpecular.y;
                    image[index++] = diffuseAndSpecular.z;
                    sphereIntersection = true;

                    break;
                }
                
            }
            for (int i = 0; i < triangles.size(); ++i)
            {
                if(!sphereIntersection && intersection(eyeRay, triangles[i].indices, scene ,t , surfaceNormal)){

                    //printf("Triangle is hit!\n");

                    Vec3f pointOnTheTriangle    = eyeRay.e + eyeRay.d*t; 

                    Vec3f vectorToLight = lightPosition - pointOnTheTriangle ; 

                    // CONTINUE WITH THE TRIANGLE SHADING

                    image[index++] = 15;
                    image[index++] = 115;
                    image[index++] = 70;
                    triangleIntersection = true;
                    break;
                }
            }

            bool breakLoop = false;

            for (int i = 0; i < scene.meshes.size(); ++i)
            {
                if(breakLoop == true){
                    break;
                }
                for (int j = 0; j < scene.meshes[i].faces.size(); ++j)
                {
                    if (!sphereIntersection && !triangleIntersection && intersection(eyeRay, scene.meshes[i].faces[j], scene ,t , surfaceNormal))
                    {




                        Vec3f pointOnTheMesh    = eyeRay.e + eyeRay.d*t; 

                        Vec3f vectorToLight = -(lightPosition - pointOnTheMesh) ;


                        float lightDistance = sqrt(dotProduct(vectorToLight,vectorToLight));

                        float cosTheta = dotProduct(vectorToLight.normalize(), surfaceNormal.normalize());

                        //printf("COSTHETA: %lf \n", cosTheta );


                        cosTheta = (cosTheta < 0) ? 0 : cosTheta;
                        

                        Vec3f diffuseShadingParams = scene.materials[meshes[i].material_id-1].diffuse; // for RGB values -> between 0 and 1


                        //printf("Diffuse parameters: %lf , %lf , %lf \n", diffuseShadingParams.x, diffuseShadingParams.y, diffuseShadingParams.z );

                        Vec3f irradiance = lightIntensity * (1.0/(lightDistance*lightDistance));


                        float diffuseShadingRed   = diffuseShadingParams.x * cosTheta * irradiance.x; 
                        float diffuseShadingGreen = diffuseShadingParams.y * cosTheta * irradiance.y; 
                        float diffuseShadingBlue  = diffuseShadingParams.z * cosTheta * irradiance.z; 

                        Vec3f diffuseShading = Vec3f(diffuseShadingRed,diffuseShadingGreen,diffuseShadingBlue);


                        //printf("Diffuse: %lf , %lf , %lf \n", diffuseShading.x, diffuseShading.y, diffuseShading.z );


                        Vec3f halfWayVector = ((-eyeRay.d).normalize() + vectorToLight.normalize()).normalize();

                        float cosAlpha = dotProduct(halfWayVector.normalize(), surfaceNormal.normalize()); // for specular shading

                        cosAlpha = (cosAlpha < 0) ? 0 : cosAlpha;


                        Vec3f specularShadingParams = scene.materials[meshes[i].material_id-1].specular; // for RGB values -> between 0 and 1
                        float phong_exponent = scene.materials[meshes[i].material_id-1].phong_exponent; // for RGB values -> between 0 and 1
                        float cosAlphaWithPhong = pow(cosAlpha,phong_exponent); 
                        //printf("Specular : %lf %lf %lf  \n", specularShadingParams.x, specularShadingParams.y, specularShadingParams.z   );


                        float specularShadingRed   = specularShadingParams.x * cosAlphaWithPhong * irradiance.x; 
                        float specularShadingGreen = specularShadingParams.y * cosAlphaWithPhong * irradiance.y; 
                        float specularShadingBlue  = specularShadingParams.z * cosAlphaWithPhong * irradiance.z; 

                        Vec3f specularShading = Vec3f(specularShadingRed,specularShadingGreen,specularShadingBlue);

                        //printf("Specular: %lf , %lf , %lf \n", specularShading.x, specularShading.y, specularShading.z );


                        Vec3f diffuseAndSpecular = clamp(diffuseShading+specularShading);



                        image[index++] = diffuseAndSpecular.x;
                        image[index++] = diffuseAndSpecular.y;
                        image[index++] = diffuseAndSpecular.z;

                        printf("FACE INTERSECTION! \n");

                        faceIntersection = true;
                        breakLoop = true;
                        break;

                        
                    }
                     
                }
            }


            if (!sphereIntersection && !triangleIntersection && !faceIntersection)
            {
                image[index++] = scene.background_color.x;
                image[index++] = scene.background_color.y;
                image[index++] = scene.background_color.z;
            } 

            //printf("INDEX : %d \n " , index);




        }
    }
/*
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

*/
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
