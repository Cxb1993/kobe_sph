/***************************************/
/* Include files                       */
/***************************************/
#include "colors.inc"
#include "shapes.inc"

/***************************************/
/* Parameter                           */
/***************************************/
#declare Arg1=0;
#declare Arg2=0;
#declare Arg3=0;
#declare Fnum = str(clock, -5, 0);
#declare Filename = concat("data", Fnum)
#debug concat("*** Filename = ", Filename, " ***\n")
#fopen MyFile Filename read

/***************************************/

/* Camera                              */
/***************************************/
// Declare camera object
camera {
  //location <100.0, 100.0, -100.0>
  //look_at  <5.0, 5.0, 5.0> 
  location <50.0, 50.0, 50.0>
  look_at  <0.0, 0.0, 0.0> 
  angle 30
  sky <0.0, 0.0, 1.0>
}

/***************************************/
/* Light source                        */
/***************************************/
// Declare light source
light_source { <150, 0, 0> color 1.5*White }
light_source { <0, 150, 0> color 1.5*White }
light_source { <0, 0, 150> color 1.5*White }


/***************************************/
/* Mapping Objects                     */
/***************************************/
object{
  Disk_X
  pigment{color Red}
  scale <5 ,0.5, 0.5>
  translate <5 ,0, 0>
}
object{
  Cone_X
  pigment{color Red}
  scale <2 ,1.5, 1.5>
  translate <12 ,0, 0>
}

object{
  Disk_X
  pigment{color Blue}
  scale <5 ,0.5, 0.5>
  translate <5 ,0, 0>
  rotate 90*z
}
object{
  Cone_X
  pigment{color Blue}
  scale <2 ,1.5, 1.5>
  translate <12 ,0, 0>
  rotate 90*z
}

object{
  Disk_X
  pigment{color Yellow}
  scale <5 ,0.5, 0.5>
  translate <5 ,0, 0>
  rotate -90*y
}
object{
  Cone_X
  pigment{color Yellow}
  scale <2 ,1.5, 1.5>
  translate <12 ,0, 0>
  rotate -90*y
}

#while (defined(MyFile))
  #read (MyFile,Arg1,Arg2,Arg3)
  object { 
    Sphere        
    pigment {
      color Green
    }
    //scale <0.5,0.5,0.5>        
    scale <0.4,0.4,0.4>        
    translate <Arg1 ,Arg2, Arg3>
  }
#end

