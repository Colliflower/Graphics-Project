/*
  CSC418 - RayTracer code - Winter 2017 - Assignment 3&4

  Written Dec. 9 2010 - Jan 20, 2011 by F. J. Estrada
  Freely distributable for adacemic purposes only.

  Uses Tom F. El-Maraghi's code for computing inverse
  matrices. You will need to compile together with
  svdDynamic.c

  You need to understand the code provided in
  this file, the corresponding header file, and the
  utils.c and utils.h files. Do not worry about
  svdDynamic.c, we need it only to compute
  inverse matrices.

  You only need to modify or add code in sections
  clearly marked "TO DO"
*/

#include "utils.h"

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
struct image *env_img;
int MAX_DEPTH;

void buildScene(void)
{
 // Sets up all objects in the scene. This involves creating each object,
 // defining the transformations needed to shape and position it as
 // desired, specifying the reflectance properties (albedos and colours)
 // and setting up textures where needed.
 // Light sources must be defined, positioned, and their colour defined.
 // All objects must be inserted in the object_list. All light sources
 // must be inserted in the light_list.
 //
 // To create hierarchical objects:
 //   Copy the transform matrix from the parent node to the child, and
 //   apply any required transformations afterwards.
 //
 // NOTE: After setting up the transformations for each object, don't
 //       forget to set up the inverse transform matrix!

 struct object3D *o;
 struct pointLS *l;
 struct point3D p;

env_img = readPPMimage("Env_Map/logcabin.ppm");

    ////////////SNOWGLOBE/////////////////////
    o=newSphere(.2,1,.6,0.15,1,1,1,.1,1.4,22);
    Scale(o,7,7,7);
    Translate(o,0,-2,15);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);
    
    o=newSphere(0,1,.6,0.15,1,1,1,0,1,22);;
    Scale(o,5.5,5.5,5.5);
    Translate(o,0,-2,15);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);

    o=newCylinder(.2,1,.6,.2,.2,.2,.8,1,1,40);
    Scale(o,6.8,6.8,-0.5);
    RotateX(o,PI/2); 
    Translate(o,0,-4.7,15);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);
 
    o=newSphere(.2,1,.6,.5,1,1,1,1,1,5);
    Scale(o,6.8,6.8,0.9);
    RotateX(o,PI/2); 
    Translate(o,0,-4.8,15);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);

    ////////////////TREE///////////////////
    double treex = -2;
    double treez = 15;
    //Trunk
    o=newCylinder(.4,1,0,0.4,.2,.1,0,1,1,1);
    Scale(o,.65,.65,-1);
    RotateX(o,PI/2); 
    Translate(o,treex,-4.3,treez);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);
    //Lower part of tree
    o=newCone(.4,1,0,0.4,0,.5,.2,1,1,1);
    Translate(o,0,0,-1);
    Scale(o,2,2,3);
    RotateX(o,PI/2); 
    RotateY(o,PI);
    Translate(o,treex,-3.2,treez);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    loadTexture(o,"tree_tex2.ppm");
    insertObject(o,&object_list);
    //Upper part of tree
    o=newCone(.4,1,0,0.4,0,.5,.2,1,1,1);
    Translate(o,0,0,-1);
    Scale(o,2,2,3);
    Scale(o,0.8,0.8,0.8);
    RotateX(o,PI/2);
    RotateY(o,PI); 
    Translate(o,treex,-1.5,treez);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    loadTexture(o,"tree_tex.ppm");
    insertObject(o,&object_list);

    ////////////////SNOWMAN///////////////////
    double snowmanx = 2;
    double snowmanz = 14;
    /////////////HAT////////////////////////
    o=newCylinder(.4,1,0,0.4,0,0,0,1,1,1);
    Scale(o,.65,.65,-.1);
    RotateX(o,PI/2);
    RotateZ(o,-PI/12); 
    Translate(o,0.6*sin(PI/12),-1.8+0.5*cos(PI/12),+0.5*sin(PI/12));
    RotateY(o,PI/7);
    Translate(o,snowmanx,0,snowmanz);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);

    o=newCylinder(.4,1,0,0.4,0,0,0,1,1,1);
    Scale(o,.4,.4,-.6);
    RotateX(o,PI/2);
    RotateZ(o,-PI/12); 
    Translate(o,0.6*sin(PI/12),-1.8+0.5*cos(PI/12),+0.5*sin(PI/12));
    RotateY(o,PI/7);
    Translate(o,snowmanx,0,snowmanz);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);
   
    ///////////HEAD/////////////
    //EYE
    o=newSphere(.4,1,0,0.4,0,0,0,1,1,1);
    Scale(o,0.1,0.1,0.01);
    Translate(o,-0.6*sin(PI/8),-1.8+0.6*sin(PI/8),-0.6*cos(PI/8));
    RotateY(o,PI/7);
    Translate(o,snowmanx,0,snowmanz);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);
    //EYE
    o=newSphere(.4,1,0,0.4,0,0,0,1,1,1);
    Scale(o,0.1,0.1,0.01);
    Translate(o,0.6*sin(PI/8),-1.8+0.6*sin(PI/8),-0.6*cos(PI/8));
    RotateY(o,PI/7);
    Translate(o,snowmanx,0,snowmanz);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);
    //NOSE
    o=newCone(.4,1,0,0.4,1,0.6,0,1,1,1);
    Translate(o,0,0,-1);
    Scale(o,0.15,0.15,.5);
    Translate(o,0,-1.8-0.6*sin(PI/36),-0.6*cos(PI/36));
    RotateY(o,PI/7);
    Translate(o,snowmanx,0,snowmanz);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);
    //HEAD
    o=newSphere(.5,1,0,0.4,1,1,1,1,1,1);
    Scale(o,0.6,0.6,0.6);
    Translate(o,0,-1.8,0);
    RotateY(o,PI/7);
    Translate(o,snowmanx,0,snowmanz);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);

    ////////////////TORSO
    //BUTTON
    o=newSphere(.4,1,0,0.4,0,0,0,1,1,1);
    Scale(o,0.1,0.1,0.1);
    Translate(o,0,-2.8+0.8*sin(PI/6),-0.8*cos(PI/6));
    RotateY(o,PI/7);
    Translate(o,snowmanx,0,snowmanz);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);
    //BUTTON
    o=newSphere(.4,1,0,0.4,0,0,0,1,1,1);
    Scale(o,0.1,0.1,0.1);
    Translate(o,0,-2.8-0.8*sin(PI/12),-0.8*cos(PI/12));
    RotateY(o,PI/7);
    Translate(o,snowmanx,0,snowmanz);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);
    //TORSO
    o=newSphere(.5,1,0,0.4,1,1,1,1,1,1);
    Scale(o,0.8,0.8,0.8);
    Translate(o,0,-2.8,0);
    RotateY(o,PI/7);
    Translate(o,snowmanx,0,snowmanz);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);

    /////////////////BASE
    //BUTTON
    o=newSphere(.4,1,0,0.4,0,0,0,1,1,1);
    Scale(o,0.1,0.1,0.1);
    Translate(o,0,-3.8+sin(PI/12),-cos(PI/12));
    RotateY(o,PI/7);
    Translate(o,snowmanx,0,snowmanz);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);
    //BASE
    o=newSphere(.5,1,0,0.4,1,1,1,1,1,1);
    Scale(o,1,1,1);
    Translate(o,0,-3.8,0);
    RotateY(o,PI/7);
    Translate(o,snowmanx,0,snowmanz);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);
    

    ///////////////////////WOODFLOOR////////////////////
    o=newPlane(.2,1,0,.5,1,1,.9,1,1,19);
    RotateX(o,PI/2);
    Scale(o,20,20,20);
    Translate(o,-.1,-4.7,20);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    loadTexture(o,"wood_tex.ppm");
    insertObject(o,&object_list);

/*
 env_img = readPPMimage("building.ppm");
    o=newSphere(.2,1,0,0.4,1,1,1,.1,1.23,17);
    //RotateY(o,PI);
    Scale(o,6,6,6);
    Translate(o,0,-2,15);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);
    
    o=newSphere(.2,1,0,0.2,1,1,1,.1,1.23,17);
    //RotateY(o,PI);
    Scale(o,3,3,3);
    Translate(o,0,-2,11);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);
    
    o=newPlane(.2,1,0,1,1,1,.9,1,1,19);
    RotateX(o,PI/2);
    Scale(o,20,20,20);
    Translate(o,-.1,-4.2,20);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    loadTexture(o,"wood_tex.ppm");
    insertObject(o,&object_list); */

 addAreaLight(4,20,0,-1,0,0,18,4,3,20,.99,.99,1,&object_list,&light_list,0);
 //addAreaLight(17,1,0,-1,0,0,18,21.5,10,4,.99,.99,1,&object_list,&light_list,0);
 
 //addAreaLight(1,15,.01,-1,0,0,16,0,9,9,.99,.99,1,&object_list,&light_list,0);
    
 /* SPACE SCENE
 env_img = readPPMimage("Env_Map/milky.ppm");
    
 o=newSphere(.05,1,0,0.01,1,1,.25,1,1,40);
 Scale(o,100,100,100);
 RotateX(o,PI/2);
 //RotateY(o,PI);
 Translate(o,0,-200,100);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 loadTexture(o,"Textures/earth.ppm");
 insertObject(o,&object_list);
    
 o=newSphere(.05,1,0,0.01,1,1,.25,1,1,40);
 Scale(o,4,4,4);
 RotateX(o,0.1);
 RotateY(o,PI);
 Translate(o,10,2,20);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 loadTexture(o,"Textures/moon.ppm");
 insertObject(o,&object_list);
    
 o=newCylinder(.05,1,0,1,.5,.5,.5,1,1,50);
 Scale(o,1/sqrt(2.0),1/sqrt(2.0),3);
 //RotateY(o,PI);
 RotateX(o,PI/3.0);
 Translate(o,10,0,0);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
    
 o=newCone(.05,1,0,1,.5,.5,.5,1,1,50);
 Scale(o,1,1,1.5);
 //RotateY(o,PI);
 Translate(o,0,0,-3/2.0);
 RotateX(o,PI/3.0);
 Translate(o,10,0,0);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 //Insert a single point light source.
 p.px=0;
 p.py=100;
 p.pz=10.0;
 p.pw=1;
 l=newPLS(&p,.95,.95,.95);
 insertPLS(l,&light_list);
  
 addAreaLight(20,20,0,1,0,0,100,10,9,9,1,1,0.95,&object_list,&light_list);*/

 // End of simple scene for Assignment 3
 // Keep in mind that you can define new types of objects such as cylinders and parametric surfaces,
 // or, you can create code to handle arbitrary triangles and then define objects as surface meshes.
 //
 // Remember: A lot of the quality of your scene will depend on how much care you have put into defining
 //           the relflectance properties of your objects, and the number and type of light sources
 //           in the scene.
}

void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double a, double b, struct colourRGB *col)
{
 // This function implements the shading model as described in lecture. It takes
 // - A pointer to the first object intersected by the ray (to get the colour properties)
 // - The coordinates of the intersection point (in world coordinates)
 // - The normal at the point
 // - The ray (needed to determine the reflection direction to use for the global component, as well as for
 //   the Phong specular component)
 // - The current racursion depth
 // - The (a,b) texture coordinates (meaningless unless texture is enabled)
 //
 // Returns:
 // - The colour for this ray (using the col pointer)
 //

 struct colourRGB tmp_col;	// Accumulator for colour components
 double R,G,B;			// Colour for the object in R G and B

 // This will hold the colour as we process all the components of
 // the Phong illumination model
 tmp_col.R=0;
 tmp_col.G=0;
 tmp_col.B=0;

 if (obj->texImg==NULL)		// Not textured, use object colour
 {
  R = obj->col.R;
  G = obj->col.G;
  B = obj->col.B;
 }
 else
 {
  // Get object colour from the texture given the texture coordinates (a,b), and the texturing function
  // for the object. Note that we will use textures also for Photon Mapping.
  obj->textureMap(obj->texImg,a,b,&R,&G,&B);
 }
 struct pointLS *currLight = light_list;
 struct point3D s,r;
 double d,rb;
 double scalar;
 double I;
 double alpha = 2;
 double norm;
 int count = 0;
 double lambda;
 double da,db;
 if (obj->isLightSource) {
  tmp_col.R = R;
  tmp_col.G = G;
  tmp_col.B = B;
 }
 else {
  while (currLight) {
   memcpy(&s, &currLight->p0, sizeof(point3D));
   subVectors(p,&s);
   s.pw = 0;
   struct ray3D *shadow_ray = newRay(p,&s);
   struct object3D *o;
   struct point3D inter_p;
   struct point3D inter_n;
   findFirstHit(shadow_ray,&lambda,obj,&o,&inter_p,&inter_n,&da,&db,1);
   if(lambda > 0 && lambda < 1) {
    tmp_col.R += currLight->col.R*R*obj->alb.ra;
    tmp_col.G += currLight->col.G*G*obj->alb.ra;
    tmp_col.B += currLight->col.B*B*obj->alb.ra;
   }
   else {
    normalize(&s);
    d = dot(&s,n);
    if (obj->frontAndBack)
     d = fabs(d);
    else if (d < 0)
     d = 0;
    r.px = -s.px + 2*d*n->px;
    r.py = -s.py + 2*d*n->py;
    r.pz = -s.pz + 2*d*n->pz;
    normalize(&r);
    rb = -dot(&ray->d,&r);
    if (obj->frontAndBack)
     rb = fabs(rb);
    else if (rb < 0)
     rb = 0;
    tmp_col.R += currLight->col.R*(R*(obj->alb.ra + obj->alb.rd*d) + obj->alb.rs*pow(rb,obj->shinyness));
    tmp_col.G += currLight->col.G*(G*(obj->alb.ra + obj->alb.rd*d) + obj->alb.rs*pow(rb,obj->shinyness));
    tmp_col.B += currLight->col.B*(B*(obj->alb.ra + obj->alb.rd*d) + obj->alb.rs*pow(rb,obj->shinyness));
   }
   currLight = currLight->next;
   count++;
   free(shadow_ray);
  }
  tmp_col.R = tmp_col.R/count;
  tmp_col.G = tmp_col.G/count;
  tmp_col.B = tmp_col.B/count;
 }
 if (tmp_col.R > 1)
  tmp_col.R = 1;
 if (tmp_col.G > 1)
  tmp_col.G = 1;
 if (tmp_col.B > 1)
  tmp_col.B = 1;
 col->R = tmp_col.R;
 col->G = tmp_col.G;
 col->B = tmp_col.B;
 //////////////////////////////////////////////////////////////
 // TO DO: Implement this function. Refer to the notes for
 // details about the shading model.
 //////////////////////////////////////////////////////////////
 // Be sure to update 'col' with the final colour computed here!
 return;

}

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b, char shadowFlag)
{
 // Find the closest intersection between the ray and any objects in the scene.
 // It returns:
 //   - The lambda at the intersection (or < 0 if no intersection)
 //   - The pointer to the object at the intersection (so we can evaluate the colour in the shading function)
 //   - The location of the intersection point (in p)
 //   - The normal at the intersection point (in n)
 //
 // Os is the 'source' object for the ray we are processing, can be NULL, and is used to ensure we don't 
 // return a self-intersection due to numerical errors for recursive raytrace calls.
 //

 /////////////////////////////////////////////////////////////
 // TO DO: Implement this function. See the notes for
 // reference of what to do in here
 /////////////////////////////////////////////////////////////
 struct object3D *currObj = object_list;
 int flag = 0;
 double currLambda = DBL_MAX;
 double nextLambda;
 struct point3D intersect;
 struct point3D normal;
 double currA,currB;
 while (currObj) {
  currObj->intersect(currObj, ray, &nextLambda,  &intersect, &normal, &currA, &currB);
  if (nextLambda > 0 && currObj != Os && nextLambda < currLambda && (!currObj->isLightSource || !shadowFlag)) {
   currLambda = nextLambda;
   *obj = currObj;
   *p = intersect;
   *n = normal;
   flag = 1;
   *a = currA;
   *b = currB;
  }
  currObj = currObj->next;
 }
 if (!flag)
  *lambda = -1;
 else
  *lambda = currLambda;
}

void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os)
{
 // Ray-Tracing function. It finds the closest intersection between
 // the ray and any scene objects, calls the shading function to
 // determine the colour at this intersection, and returns the
 // colour.
 //
 // Os is needed for recursive calls to ensure that findFirstHit will
 // not simply return a self-intersection due to numerical
 // errors. For the top level call, Os should be NULL. And thereafter
 // it will correspond to the object from which the recursive
 // ray originates.
 //

 double lambda;		// Lambda at intersection
 double a,b;		// Texture coordinates
 struct object3D *obj;	// Pointer to object at intersection
 struct point3D p;	// Intersection point
 struct point3D n;	// Normal at intersection
 struct colourRGB I;	// Colour returned by shading function
 int num_refls = 1;

 if (depth>MAX_DEPTH)	// Max recursion depth reached. Return invalid colour.
 {
  col->R=-1;
  col->G=-1;
  col->B=-1;
  return;
 }

 ///////////////////////////////////////////////////////
 // TO DO: Complete this function. Refer to the notes
 // if you are unsure what to do here.
 ///////////////////////////////////////////////////////
 
 findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b, 0);
 if (lambda > 0) {
  rtShade(obj, &p, &n, ray, depth, a, b, col);
  struct colourRGB refcol;
  struct colourRGB collector;
  memset(&collector, 0, sizeof(struct colourRGB));
  struct point3D r,*u,*v;
  struct ray3D *reflray;
  struct ray3D *refrray;
  double d  = dot(&ray->d, &n);
  if (obj->alpha <1 && obj->r_index != 0) {
   double A = 1/obj->r_index;
   double dot_prod  = dot(&ray->d, &n);
   if(dot_prod>0){ //We're inside the object, swap the normal for both reflection and refraction
    A = 1/A;
    n.px = -n.px;
    n.py = -n.py;
    n.pz = -n.pz;
    dot_prod  = dot(&ray->d, &n);
   }
   double radicand = 1 - A*A*(1-dot_prod*dot_prod);
   if(radicand >= 0){
    double B = A*dot_prod + sqrt(radicand);
    r.px = A*ray->d.px - B*n.px;
    r.py = A*ray->d.py - B*n.py;
    r.pz = A*ray->d.pz - B*n.pz;
    r.pw = 0;
    normalize(&r);
    refrray = newRay(&p, &r);
    rayTrace(refrray, depth, &refcol, obj);
    if (refcol.R != -1) {
     col->R = (1-obj->alpha)*refcol.R + obj->alpha*col->R;
     col->G = (1-obj->alpha)*refcol.G + obj->alpha*col->G;
     col->B = (1-obj->alpha)*refcol.B + obj->alpha*col->B;
    }
    free(refrray);
   }
  }
  if (obj->alb.rg != 0) {
   for (int i = 0; i < num_refls; i++) {
    r.px = ray->d.px - 2*d*n.px;
    r.py = ray->d.py - 2*d*n.py;
    r.pz = ray->d.pz - 2*d*n.pz;
    normalize(&r);
    u = cross(&r,&n);
    v = cross(&r,u);
    double theta, phi, roughness, x, y, z;
    roughness = (1-tanh(obj->shinyness/5-3))/2;
    theta = M_PI*roughness*drand48();
    phi = 2*M_PI*roughness*drand48();
    x = sin(theta)*cos(phi);
    y = sin(theta)*sin(phi);
    z = cos(theta);
    r.px = u->px*x + v->px*y + r.px*z;
    r.py = u->py*x + v->py*y + r.py*z;
    r.pz = u->pz*x + v->pz*y + r.pz*z;
    r.pw = 0;
    free(u);
    free(v);
    normalize(&r);
    reflray = newRay(&p, &r);
    rayTrace(reflray, depth+1, &refcol, obj);
    if (refcol.R != -1) {
     collector.R += refcol.R;
     collector.G += refcol.G;
     collector.B += refcol.B;
    }
    free(reflray);
   }
  }
  col->R += obj->alb.rg*collector.R/(double)num_refls;
  col->G += obj->alb.rg*collector.G/(double)num_refls;
  col->B += obj->alb.rg*collector.B/(double)num_refls;
  if (col->R > 1)
   col->R = 1;
  if (col->G > 1)
   col->G = 1;
  if (col->B > 1)
   col->B = 1;
 }
 else {
  if(env_img!=NULL){
   double R,G,B;
   b = acos(ray->d.py)/PI;
   a = atan2(ray->d.px,ray->d.pz)/(2*PI) + 0.5;
   texMap(env_img,a,b,&R,&G,&B);
   col->R=R;
   col->G=G;
   col->B=B;
  } else {
   col->R=-1;
   col->G=-1;
   col->B=-1;
  }
 }
}

int main(int argc, char *argv[])
{
 // Main function for the raytracer. Parses input parameters,
 // sets up the initial blank image, and calls the functions
 // that set up the scene and do the raytracing.
 struct image *im;	// Will hold the raytraced image
 struct view *cam;	// Camera and view for this scene
 int sx,sy;		// Size of the raytraced image
 int antialiasing;	// Flag to determine whether antialiaing is enabled or disabled
 char output_name[1024];	// Name of the output file for the raytraced .ppm image
 struct point3D e;		// Camera view parameters 'e', 'g', and 'up'
 struct point3D g;
 struct point3D up;
 double du, dv;			// Increase along u and v directions for pixel coordinates
 struct colourRGB background;   // Background colour
 int i,j,k,l;			// Counters for pixel coordinates
 unsigned char *rgbIm;
 char VR_flag = 0;
 char DOF = 1;
 double F = 14; //Focal Length
 double R = .25; //Apature Radius
 struct point3D C; //Focal Point
 double theta,r; //Used to sample from circle

 if (argc<5)
 {
  fprintf(stderr,"RayTracer: Can not parse input parameters\n");
  fprintf(stderr,"USAGE: RayTracer size rec_depth antialias output_name\n");
  fprintf(stderr,"   size = Image size (both along x and y)\n");
  fprintf(stderr,"   rec_depth = Recursion depth\n");
  fprintf(stderr,"   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
  fprintf(stderr,"   output_name = Name of the output file, e.g. MyRender.ppm\n");
  exit(0);
 }
 if (argc>5) {
  fprintf(stderr,"VR MODE ENGAGED\n");
  VR_flag = (atoi(argv[5]) > 0);
 }
 else {
  fprintf(stderr,"If you want 360 images add a 1 as an additional argument\n");
 }
 sy=atoi(argv[1]);
 sx=(1+VR_flag)*sy;
 MAX_DEPTH=atoi(argv[2]);
 if (atoi(argv[3])==0) antialiasing=0; else antialiasing=1;
 strcpy(&output_name[0],argv[4]);

 fprintf(stderr,"Rendering image at %d x %d\n",sx,sy);
 fprintf(stderr,"Recursion depth = %d\n",MAX_DEPTH);
 if (!antialiasing) fprintf(stderr,"Antialising is off\n");
 else fprintf(stderr,"Antialising is on\n");
 fprintf(stderr,"Output file name: %s\n",output_name);

 object_list=NULL;
 light_list=NULL;

 // Allocate memory for the new image
 im=newImage(sx, sy);
 if (!im)
 {
  fprintf(stderr,"Unable to allocate memory for raytraced image\n");
  exit(0);
 }
 else rgbIm=(unsigned char *)im->rgbdata;

 ///////////////////////////////////////////////////
 // TO DO: You will need to implement several of the
 //        functions below. For Assignment 3, you can use
 //        the simple scene already provided. But
 //        for Assignment 4 you need to create your own
 //        *interesting* scene.
 ///////////////////////////////////////////////////
 buildScene();		// Create a scene. This defines all the
			// objects in the world of the raytracer

 //////////////////////////////////////////
 // TO DO: For Assignment 3 you can use the setup
 //        already provided here. For Assignment 4
 //        you may want to move the camera
 //        and change the view parameters
 //        to suit your scene.
 //////////////////////////////////////////

 // Mind the homogeneous coordinate w of all vectors below. DO NOT
 // forget to set it to 1, or you'll get junk out of the
 // geometric transformations later on.

 // Camera center is at (0,0,-1)
 e.px=0;
 e.py=0;
 e.pz=0;
 e.pw=1;

 // To define the gaze vector, we choose a point 'pc' in the scene that
 // the camera is looking at, and do the vector subtraction pc-e.
 // Here we set up the camera to be looking at the origin, so g=(0,0,0)-(0,0,-1)
 g.px=0;
 g.py=0;
 g.pz=1;
 g.pw=1;

 // Define the 'up' vector to be the Y axis
 up.px=0;
 up.py=1;
 up.pz=0;
 up.pw=1;

 // Set up view with given the above vectors, a 4x4 window,
 // and a focal length of -1 (why? where is the image plane?)
 // Note that the top-left corner of the window is at (-2, 2)
 // in camera coordinates.
 cam=setupView(&e, &g, &up, -3, -2, 2, 4);

 if (cam==NULL)
 {
  fprintf(stderr,"Unable to set up the view and camera parameters. Our of memory!\n");
  cleanup(object_list,light_list);
  deleteImage(im);
  exit(0);
 }

 // Set up background colour here
 background.R=0;
 background.G=0;
 background.B=0;

 // Do the raytracing
 //////////////////////////////////////////////////////
 // TO DO: You will need code here to do the raytracing
 //        for each pixel in the image. Refer to the
 //        lecture notes, in particular, to the
 //        raytracing pseudocode, for details on what
 //        to do here. Make sure you undersand the
 //        overall procedure of raytracing for a single
 //        pixel.
 //////////////////////////////////////////////////////
 du=cam->wsize/(sx-1);		// du and dv. In the notes in terms of wl and wr, wt and wb,
 dv=-cam->wsize/(sx-1);		// here we use wl, wt, and wsize. du=dv since the image is
				// and dv is negative since y increases downward in pixel
				// coordinates and upward in camera coordinates.

 fprintf(stderr,"View parameters:\n");
 fprintf(stderr,"Left=%f, Top=%f, Width=%f, f=%f\n",cam->wl,cam->wt,cam->wsize,cam->f);
 fprintf(stderr,"Camera to world conversion matrix (make sure it makes sense!):\n");
 printmatrix(cam->C2W);
 fprintf(stderr,"World to camera conversion matrix\n");
 printmatrix(cam->W2C);
 fprintf(stderr,"\n");

 
 fprintf(stderr,"Rendering row: ");
 double aa_res = atof(argv[3]);
 if (!antialiasing)
  aa_res = 1;
 #pragma omp parallel for schedule(dynamic) num_threads(4)
 for (int j=0;j<sy;j++)		// For each of the pixels in the image
 {
  for (int i=0;i<sx;i++)
  {
   if(i==0){
    fprintf(stderr,"%d/%d, ",1+j,sy);
   }
   struct colourRGB pixelcol;
   memset(&pixelcol,0,sizeof(struct colourRGB));
   for (int k=0;k<aa_res;k++) {
    for (int l=0;l<aa_res;l++) {
     struct point3D p0;
     double jitter_x, jitter_y;
     struct point3D pc,d;		// Point structures to keep the coordinates of a pixel and
                           // the direction or a ray
     struct ray3D *ray;		// Structure to keep the ray from e to a pixel
     struct colourRGB col;		// Return colour for raytraced pixels
     double theta,phi;
     //struct object3D Os;
     ///////////////////////////////////////////////////////////////////
     // TO DO - complete the code that should be in this loop to do the
     //         raytracing!
     ////////////////a///////////////////////////////////////////////////
     if (antialiasing) {
      jitter_x = ((double)k + drand48())/aa_res;
      jitter_y = ((double)l + drand48())/aa_res;
     }
     else {
      jitter_x = 0;
      jitter_y = 0;
     }
     if (VR_flag) {
      theta = 2*PI*(-((double)sx)/2 + i + jitter_x + 0.5)/(double)sx;
      phi = PI*(-((double)sy)/2 + j + jitter_y + 0.5)/(double)sy;
      d.px = cos(theta)*cos(phi);
      d.py = sin(phi);
      d.pz = sin(theta)*cos(phi);
      d.pw = 0;
      p0.px = 0;
      p0.py = 0;
      p0.pz = 0;
      p0.pw = 1;
     }
     else {
      if(DOF){
       theta = drand48()*2*PI;
       r = R*sqrt(drand48());
       double rx = r*cos(theta);
       double ry = r*sin(theta);
       d.px = 1.5*(-((double)sx)/2 + i + jitter_x + 0.5)/(double)sx;
       d.py = 1.5*(-((double)sy)/2 + j + jitter_y + 0.5)/(double)sy;
       d.pz = -1;
       d.pw = 0;
       matVecMult(cam->C2W, &d);
       normalize(&d);
       C.px = F*d.px;
       C.py = F*d.py;
       C.pz = F*d.pz;
       matVecMult(cam->W2C, &C);
       p0.px = rx;
       p0.py = ry;
       p0.pz = 0;
       p0.pw = 1;
       d.px = C.px - p0.px;
       d.py = C.py - p0.py;
       d.pz = C.pz - p0.pz;
       d.pw = 0;
      } else {
       d.px = 1.5*(-((double)sx)/2 + i + jitter_x + 0.5)/(double)sx;
       d.py = 1.5*(-((double)sy)/2 + j + jitter_y + 0.5)/(double)sy;
       d.pz = -1;
       d.pw = 0;
       p0.px = 0;
       p0.py = 0;
       p0.pz = 0;
       p0.pw = 1;
      }
     }
     matVecMult(cam->C2W, &p0);
     matVecMult(cam->C2W, &d);
     //printf("%f, %f, %f\n",d.px,d.py,d.pz);
     normalize(&d);
     ray = newRay(&p0,&d);
     rayTrace(ray, 0, &col, NULL);
     /*col.R = (char)ceil(cos(20*PI*(i+1)/(double)sx+PI/2)) ^ (char)ceil(cos(20*PI*(j+1)/(double)sx+PI/2));
     col.G = col.R;
     col.B = col.R;*/
     free(ray);
     if (col.R >= 0) {
      pixelcol.R += col.R;
      pixelcol.G += col.G;
      pixelcol.B += col.B;
     }
    }
   }
   ((unsigned char *)im->rgbdata)[((sy-j-1)*sx + i)*3] = (unsigned char)(255*pixelcol.R/pow(aa_res,2));
   ((unsigned char *)im->rgbdata)[((sy-j-1)*sx + i)*3 + 1] = (unsigned char)(255*pixelcol.G/pow(aa_res,2));
   ((unsigned char *)im->rgbdata)[((sy-j-1)*sx + i)*3 + 2] = (unsigned char)(255*pixelcol.B/pow(aa_res,2));
  } // end for i
 } // end for j

 
 fprintf(stderr,"\nDone!\n");

 // Output rendered image
 imageOutput(im,output_name);

 // Exit section. Clean up and return.
 cleanup(object_list,light_list);		// Object and light lists
 deleteImage(im);				// Rendered image
 free(cam);					// camera view
 exit(0);
}
