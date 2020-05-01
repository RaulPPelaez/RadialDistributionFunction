/*Raul P. Pelaez 2016. vector types algebra*/
#ifndef VECTOR_OVERLOADS_H
#define VECTOR_OVERLOADS_H
#include"defines.h"
#ifdef GPUMODE
#include <cuda_runtime.h>
#define VECATTR inline __host__ __device__
#else
struct double2{double x,y;};
struct double3{double x,y,z;};
struct double4{double x,y,z,w;};

struct float2{float x,y;};
struct float3{float x,y,z;};
struct float4{float x,y,z,w;};

struct int2{int x,y;};
struct int3{int x,y,z;};
struct int4{int x,y,z,w;};

#define VECATTR inline
VECATTR int2 make_int2(int x, int y){return {x,y};}
VECATTR int3 make_int3(int x, int y, int z){return {x,y,z};}
VECATTR int4 make_int4(int x, int y, int z, int w){return {x,y,z,w};}

VECATTR float2 make_float2(float x, float y){return {x,y};}
VECATTR float3 make_float3(float x, float y, float z){return {x,y,z};}
VECATTR float4 make_float4(float x, float y, float z, float w){return {x,y,z,w};}

VECATTR double2 make_double2(double x, double y){return {x,y};}
VECATTR double3 make_double3(double x, double y, double z){return {x,y,z};}
VECATTR double4 make_double4(double x, double y, double z, double w){return {x,y,z,w};}


#endif
#include <cmath>

#ifndef SINGLE_PRECISION
typedef  double  real ;
typedef  double2 real2;
typedef  double3 real3;
typedef  double4 real4;
#else
typedef  float  real ;
typedef  float2 real2;
typedef  float3 real3;
typedef  float4 real4;
#endif

typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long long int ullint;





  /////////////////////FLOAT2///////////////////////////////

  VECATTR int2 make_int2(float2 a){return make_int2((int)a.x, (int)a.y);}
  VECATTR float2 make_float2(float a){return make_float2(a, a);}

  VECATTR float2 make_float2(int2 a){return make_float2(a.x, a.y);}
  VECATTR float2 make_float2(float2 a){return make_float2(a.x, a.y);}

  VECATTR  float2 operator +(const float2 &a, const float2 &b){return make_float2(
										  a.x + b.x,
										  a.y + b.y);
  }
  VECATTR  void operator +=(float2 &a, const float2 &b){
    a.x += b.x;
    a.y += b.y;
  }
  VECATTR  float2 operator +(const float2 &a, const float &b){return make_float2(
										 a.x + b,
										 a.y + b);
  }
  VECATTR  float2 operator +(const float &b, const float2 &a){return a+b;}
  VECATTR  void operator +=(float2 &a, const float &b){
    a.x += b;
    a.y += b;
  }

  VECATTR  float2 operator -(const float2 &a, const float2 &b){return make_float2(
										  a.x - b.x,
										  a.y - b.y);
  }

  VECATTR  void operator -=(float2 &a, const float2 &b){
    a.x -= b.x;
    a.y -= b.y;
  }

  VECATTR  float2 operator -(const float2 &a, const float &b){return make_float2(
										 a.x - b,
										 a.y - b);
  }

  VECATTR  float2 operator -(const float &b, const float2 &a){return make_float2(
										 b-a.x,
										 b-a.y);
  }
  VECATTR  void operator -=(float2 &a, const float &b){
    a.x -= b;
    a.y -= b;
  }
  VECATTR  float2 operator *(const float2 &a, const float2 &b){
    return make_float2(a.x * b.x,
		       a.y * b.y);
  }
  VECATTR  void operator *=(float2 &a, const float2 &b){
    a.x *= b.x;
    a.y *= b.y;
  }
  VECATTR  float2 operator *(const float2 &a, const float &b){
    return make_float2(a.x * b,
		       a.y * b);
  }
  VECATTR  float2 operator *(const float &b, const float2 &a){
    return make_float2(a.x * b,
		       a.y * b);
  }
  VECATTR  void operator *=(float2 &a, const float &b){
    a.x *= b;
    a.y *= b;
  }
  VECATTR  float2 operator /(const float2 &a, const float2 &b){
    return make_float2(a.x / b.x,
		       a.y / b.y);
  }
  VECATTR  void operator /=(float2 &a, const float2 &b){
    a.x /= b.x;
    a.y /= b.y;
  }
  VECATTR  float2 operator /(const float2 &a, const float &b){
    return (1.0f/b)*a;
  }
  VECATTR  float2 operator /(const float &b, const float2 &a){
    return make_float2(b / a.x,
		       b / a.y);
  }
  VECATTR  void operator /=(float2 &a, const float &b){
    a *= 1.0f/b;
  }

  VECATTR  float2 floorf(const float2 &a){
    return make_float2(floorf(a.x), floorf(a.y));
  }


  /////////////////////FLOAT3///////////////////////////////

  VECATTR int3 make_int3(float3 a){return make_int3((int)a.x, (int)a.y, (int)a.z);}
  VECATTR float3 make_float3(float a){return make_float3(a, a, a);}

  VECATTR float3 make_float3(int3 a){return make_float3(a.x, a.y, a.z);}
  VECATTR float3 make_float3(float3 a){return make_float3(a.x, a.y, a.z);}

  VECATTR  float3 operator +(const float3 &a, const float3 &b){return make_float3(
										  a.x + b.x,
										  a.y + b.y,
										  a.z + b.z);
  }
  VECATTR  void operator +=(float3 &a, const float3 &b){
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
  }
  VECATTR  float3 operator +(const float3 &a, const float &b){return make_float3(
										 a.x + b,
										 a.y + b,
										 a.z + b);
  }
  VECATTR  float3 operator +(const float &b, const float3 &a){return a+b;}
  VECATTR  void operator +=(float3 &a, const float &b){
    a.x += b;
    a.y += b;
    a.z += b;
  }

  VECATTR  float3 operator -(const float3 &a, const float3 &b){return make_float3(
										  a.x - b.x,
										  a.y - b.y,
										  a.z - b.z);
  }

  VECATTR  void operator -=(float3 &a, const float3 &b){
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
  }

  VECATTR  float3 operator -(const float3 &a, const float &b){return make_float3(
										 a.x - b,
										 a.y - b,
										 a.z - b);
  }

  VECATTR  float3 operator -(const float &b, const float3 &a){return make_float3(
										 b-a.x,
										 b-a.y,
										 b-a.z);
  }
  VECATTR  void operator -=(float3 &a, const float &b){
    a.x -= b;
    a.y -= b;
    a.z -= b;
  }
  VECATTR  float3 operator *(const float3 &a, const float3 &b){
    return make_float3(
		       a.x * b.x,
		       a.y * b.y,
		       a.z * b.z
		       );
  }
  VECATTR  void operator *=(float3 &a, const float3 &b){
    a.x *= b.x;
    a.y *= b.y;
    a.z *= b.z;
  }
  VECATTR  float3 operator *(const float3 &a, const float &b){
    return make_float3(
		       a.x * b,
		       a.y * b,
		       a.z * b
		       );
  }
  VECATTR  float3 operator *(const float &b, const float3 &a){
    return make_float3(
		       a.x * b,
		       a.y * b,
		       a.z * b
		       );
  }
  VECATTR  void operator *=(float3 &a, const float &b){
    a.x *= b;
    a.y *= b;
    a.z *= b;
  }
  VECATTR  float3 operator /(const float3 &a, const float3 &b){
    return make_float3(
		       a.x / b.x,
		       a.y / b.y,
		       a.z / b.z
		       );
  }
  VECATTR  void operator /=(float3 &a, const float3 &b){
    a.x /= b.x;
    a.y /= b.y;
    a.z /= b.z;
  }
  VECATTR  float3 operator /(const float3 &a, const float &b){
    return (1.0f/b)*a;
  }
  VECATTR  float3 operator /(const float &b, const float3 &a){
    return make_float3(
		       b / a.x,
		       b / a.y,
		       b / a.z
		       );
  }
  VECATTR  void operator /=(float3 &a, const float &b){
    a *= 1.0f/b;
  }

  VECATTR  float3 floorf(const float3 &a){return make_float3(floorf(a.x), floorf(a.y), floorf(a.z));}

  /////////////////////FLOAT4///////////////////////////////


  VECATTR float4 make_float4(float a){return make_float4(a,a,a,a);}

  VECATTR float4 make_float4(float3 a){return make_float4(a.x, a.y, a.z, 0);}
  VECATTR float4 make_float4(float4 a){return make_float4(a.x, a.y, a.z, a.w);}

  VECATTR  float4 operator +(const float4 &a, const float4 &b){return make_float4(
										  a.x + b.x,
										  a.y + b.y,
										  a.z + b.z,
										  a.w + b.w);
  }
  VECATTR  void operator +=(float4 &a, const float4 &b){
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    a.w += b.w;
  }
  VECATTR  float4 operator +(const float4 &a, const float &b){return make_float4(
										 a.x + b,
										 a.y + b,
										 a.z + b,
										 a.w + b);
  }
  VECATTR  float4 operator +(const float &b, const float4 &a){return a+b;}
  VECATTR  void operator +=(float4 &a, const float &b){
    a.x += b;
    a.y += b;
    a.z += b;
    a.w += b;
  }

  VECATTR  float4 operator -(const float4 &a, const float4 &b){return make_float4(
										  a.x - b.x,
										  a.y - b.y,
										  a.z - b.z,
										  a.w - b.w);
  }

  VECATTR  void operator -=(float4 &a, const float4 &b){
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    a.w -= b.w;
  }

  VECATTR  float4 operator -(const float4 &a, const float &b){return make_float4(
										 a.x - b,
										 a.y - b,
										 a.z - b,
										 a.w - b);
  }

  VECATTR  float4 operator -(const float &b, const float4 &a){return make_float4(
										 b-a.x,
										 b-a.y,
										 b-a.z,
										 b-a.w);
  }
  VECATTR  void operator -=(float4 &a, const float &b){
    a.x -= b;
    a.y -= b;
    a.z -= b;
    a.w -= b;
  }
  VECATTR  float4 operator *(const float4 &a, const float4 &b){
    return make_float4(a.x * b.x,
		       a.y * b.y,
		       a.z * b.z,
		       a.w * b.w);
  }
  VECATTR  void operator *=(float4 &a, const float4 &b){
    a.x *= b.x;
    a.y *= b.y;
    a.z *= b.z;
    a.w *= b.w;
  }
  VECATTR  float4 operator *(const float4 &a, const float &b){
    return make_float4(a.x * b,
		       a.y * b,
		       a.z * b,
		       a.w * b);
  }
  VECATTR  float4 operator *(const float &b, const float4 &a){
    return make_float4(a.x * b,
		       a.y * b,
		       a.z * b,
		       a.w * b);
  }
  VECATTR  void operator *=(float4 &a, const float &b){
    a.x *= b;
    a.y *= b;
    a.z *= b;
    a.w *= b;
  }
  VECATTR  float4 operator /(const float4 &a, const float4 &b){
    return make_float4(a.x / b.x,
		       a.y / b.y,
		       a.z / b.z,
		       a.w / b.w);
  }
  VECATTR  void operator /=(float4 &a, const float4 &b){
    a.x /= b.x;
    a.y /= b.y;
    a.z /= b.z;
    a.w /= b.w;
  }
  VECATTR  float4 operator /(const float4 &a, const float &b){
    return (1.0f/b)*a;
  }
  VECATTR  float4 operator /(const float &b, const float4 &a){
    return make_float4(b / a.x,
		       b / a.y,
		       b / a.z,
		       b / a.w);
  }
  VECATTR  void operator /=(float4 &a, const float &b){
    a *= 1.0f/b;
  }


  VECATTR  float4 floorf(const float4 &a){
    return make_float4(floorf(a.x), floorf(a.y), floorf(a.z), floorf(a.w));
  }

  VECATTR float dot(float4 a, float4 b){return a.x*b.x + a.y*b.y + a.z*b.z + a.w*b.w;}



  /////////////////REAL4////////////////////////////////
  VECATTR real4 make_real4(real x, real y, real z, real w){
#ifdef SINGLE_PRECISION
    return make_float4(x,y,z,w);
#else
    return make_double4(x,y,z,w);
#endif
  }

  VECATTR real4 make_real4(real s){return make_real4(s, s, s, s);}
  VECATTR real4 make_real4(real3 a){ return make_real4(a.x, a.y, a.z, 0.0f);}
  VECATTR real4 make_real4(real3 a, real w){ return make_real4(a.x, a.y, a.z, w);}

#ifdef SINGLE_PRECISION
  VECATTR real4 make_real4(double3 a, real w){return make_real4(a.x, a.y, a.z, w);}
#else
  VECATTR real4 make_real4(float3 a, real w){ return make_real4(a.x, a.y, a.z, w);}
#endif

  VECATTR real4 make_real4(int4 a){ return make_real4(real(a.x), real(a.y), real(a.z), real(a.w));}


  //////////////////REAL3///////////////////////////


  VECATTR real3 make_real3(real x, real y, real z){
#ifdef SINGLE_PRECISION
    return make_float3(x,y,z);
#else
    return make_double3(x,y,z);
#endif
  }

  VECATTR real3 make_real3(real s){ return make_real3(s, s, s);}
  VECATTR real3 make_real3(real3 a){return make_real3(a.x, a.y, a.z);}

#ifdef SINGLE_PRECISION
  VECATTR real3 make_real3(double3 a){return make_real3(a.x, a.y, a.z);}
  VECATTR real3 make_real3(double4 a){return make_real3(a.x, a.y, a.z);}
#else
  VECATTR real3 make_real3(float3 a){return make_real3(a.x, a.y, a.z);}
  VECATTR real3 make_real3(float4 a){return make_real3(a.x, a.y, a.z);}
#endif
  VECATTR real3 make_real3(real4 a){ return make_real3(a.x, a.y, a.z);}

  VECATTR real3 make_real3(real2 a, real z){return make_real3(a.x, a.y, z);}
  VECATTR real3 make_real3(int3 a){ return make_real3(real(a.x), real(a.y), real(a.z));}


  //////////////////REAL2///////////////////////////


  VECATTR real2 make_real2(real x, real y){
#ifdef SINGLE_PRECISION
    return make_float2(x,y);
#else
    return make_double2(x,y);
#endif
  }

  VECATTR real2 make_real2(real s){ return make_real2(s, s);}
  VECATTR real2 make_real2(real2 a){return make_real2(a.x, a.y);}
  VECATTR real2 make_real2(real4 a){return make_real2(a.x, a.y);}
  VECATTR real2 make_real2(int3 a){ return make_real2(real(a.x), real(a.y));}



  ////////////////DOUBLE PRECISION//////////////////////
#ifndef SINGLE_PRECISION
  VECATTR double3 make_double3(real3 a){return make_double3(a.x, a.y, a.z);}
#endif
  VECATTR float4 make_float4(double4 a){return make_float4(float(a.x), float(a.y), float(a.z), float(a.w));}

  VECATTR double4 make_double4(double s){ return make_double4(s, s, s, s);}
  VECATTR double4 make_double4(double3 a){return make_double4(a.x, a.y, a.z, 0.0f);}
  VECATTR double4 make_double4(double3 a, double w){return make_double4(a.x, a.y, a.z, w);}
  VECATTR double4 make_double4(int4 a){return make_double4(double(a.x), double(a.y), double(a.z), double(a.w));}

  VECATTR double4 make_double4(float4 a){return make_double4(double(a.x), double(a.y), double(a.z), double(a.w));}

  //////DOUBLE4///////////////
  VECATTR  double4 operator +(const double4 &a, const double4 &b){
    return make_double4(a.x + b.x,
			a.y + b.y,
			a.z + b.z,
			a.w + b.w
			);
  }
  VECATTR  void operator +=(double4 &a, const double4 &b){
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    a.w += b.w;
  }
  VECATTR  double4 operator +(const double4 &a, const double &b){
    return make_double4(
			a.x + b,
			a.y + b,
			a.z + b,
			a.w + b
			);
  }
  VECATTR  double4 operator +(const double &b, const double4 &a){
    return a+b;
  }
  VECATTR  void operator +=(double4 &a, const double &b){
    a.x += b;
    a.y += b;
    a.z += b;
    a.w += b;
  }

  VECATTR  double4 operator -(const double4 &a, const double4 &b){
    return make_double4(
			a.x - b.x,
			a.y - b.y,
			a.z - b.z,
			a.w - b.w
			);
  }
  VECATTR  void operator -=(double4 &a, const double4 &b){
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    a.w -= b.w;
  }
  VECATTR  double4 operator -(const double4 &a, const double &b){
    return make_double4(
			a.x - b,
			a.y - b,
			a.z - b,
			a.w - b
			);
  }
  VECATTR  double4 operator -(const double &b, const double4 &a){
    return make_double4(
			b - a.x,
			b - a.y,
			b - a.z,
			b - a.w
			);
  }
  VECATTR  void operator -=(double4 &a, const double &b){
    a.x -= b;
    a.y -= b;
    a.z -= b;
    a.w -= b;
  }
  VECATTR  double4 operator *(const double4 &a, const double4 &b){
    return make_double4(
			a.x * b.x,
			a.y * b.y,
			a.z * b.z,
			a.w * b.w
			);
  }
  VECATTR  void operator *=(double4 &a, const double4 &b){
    a.x *= b.x;
    a.y *= b.y;
    a.z *= b.z;
    a.w *= b.w;
  }
  VECATTR  double4 operator *(const double4 &a, const double &b){
    return make_double4(
			a.x * b,
			a.y * b,
			a.z * b,
			a.w * b
			);
  }
  VECATTR  double4 operator *(const double &b, const double4 &a){
    return a*b;
  }
  VECATTR  void operator *=(double4 &a, const double &b){
    a.x *= b;
    a.y *= b;
    a.z *= b;
    a.w *= b;
  }
  VECATTR  double4 operator /(const double4 &a, const double4 &b){
    return make_double4(
			a.x / b.x,
			a.y / b.y,
			a.z / b.z,
			a.w / b.w
			);
  }
  VECATTR  void operator /=(double4 &a, const double4 &b){
    a.x /= b.x;
    a.y /= b.y;
    a.z /= b.z;
    a.w /= b.w;
  }
  VECATTR  double4 operator /(const double4 &a, const double &b){return (1.0/b)*a;}
  VECATTR  double4 operator /(const double &b, const double4 &a){
    return make_double4(
			b / a.x,
			b / a.y,
			b / a.z,
			b / a.w
			);
  }
  VECATTR  void operator /=(double4 &a, const double &b){
    a *= 1.0/b;
  }

  VECATTR double dot(double4 a, double4 b)
  {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
  }
  VECATTR double length(double4 v)
  {
    return sqrt(dot(v, v));
  }
  VECATTR double4 normalize(double4 v)
  {
    double invLen = 1.0/sqrt(dot(v, v));
    return v * invLen;
  }
  VECATTR double4 floorf(double4 v)
  {
    return make_double4(floor(v.x), floor(v.y), floor(v.z), floor(v.w));
  }

  /////////////////////DOUBLE3///////////////////////////////

  VECATTR int3 make_int3(double3 a){
    return make_int3((int)a.x, (int)a.y, (int)a.z);
  }
  VECATTR double3 make_double3(double a){
    return make_double3(a, a, a);
  }

  VECATTR double3 make_double3(int3 a){
    return make_double3(a.x, a.y, a.z);
  }
  VECATTR double3 make_double3(float3 a){
    return make_double3(a.x, a.y, a.z);
  }

  VECATTR  double3 operator +(const double3 &a, const double3 &b){
    return make_double3(
			a.x + b.x,
			a.y + b.y,
			a.z + b.z
			);
  }
  VECATTR  void operator +=(double3 &a, const double3 &b){
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
  }
  VECATTR  double3 operator +(const double3 &a, const double &b){
    return make_double3(
			a.x + b,
			a.y + b,
			a.z + b
			);
  }
  VECATTR  double3 operator +(const double &b, const double3 &a){
    return a+b;
  }
  VECATTR  void operator +=(double3 &a, const double &b){
    a.x += b;
    a.y += b;
    a.z += b;
  }

  VECATTR  double3 operator -(const double3 &a, const double3 &b){
    return make_double3(
			a.x - b.x,
			a.y - b.y,
			a.z - b.z
			);
  }
  VECATTR  void operator -=(double3 &a, const double3 &b){
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
  }
  VECATTR  double3 operator -(const double3 &a, const double &b){
    return make_double3(
			a.x - b,
			a.y - b,
			a.z - b
			);
  }
  VECATTR  double3 operator -(const double &b, const double3 &a){
    return make_double3(
			b-a.x,
			b-a.y,
			b-a.z
			);
  }
  VECATTR  void operator -=(double3 &a, const double &b){
    a.x -= b;
    a.y -= b;
    a.z -= b;
  }
  VECATTR  double3 operator *(const double3 &a, const double3 &b){
    return make_double3(
			a.x * b.x,
			a.y * b.y,
			a.z * b.z
			);
  }
  VECATTR  void operator *=(double3 &a, const double3 &b){
    a.x *= b.x;
    a.y *= b.y;
    a.z *= b.z;
  }
  VECATTR  double3 operator *(const double3 &a, const double &b){
    return make_double3(
			a.x * b,
			a.y * b,
			a.z * b
			);
  }
  VECATTR  double3 operator *(const double &b, const double3 &a){
    return a*b;
  }
  VECATTR  void operator *=(double3 &a, const double &b){
    a.x *= b;
    a.y *= b;
    a.z *= b;
  }
  VECATTR  double3 operator /(const double3 &a, const double3 &b){
    return make_double3(
			a.x / b.x,
			a.y / b.y,
			a.z / b.z
			);
  }
  VECATTR  void operator /=(double3 &a, const double3 &b){
    a.x /= b.x;
    a.y /= b.y;
    a.z /= b.z;
  }
  VECATTR  double3 operator /(const double3 &a, const double &b){return (1.0/b)*a;}

  VECATTR  double3 operator /(const double &b, const double3 &a){
    return make_double3(
			b / a.x,
			b / a.y,
			b / a.z
			);
  }
  VECATTR  void operator /=(double3 &a, const double &b){

    a *= 1.0/b;

  }

  //DOUBLE2


  VECATTR  double2 operator -(const double2 &a, const double2 &b){
    return make_double2(
			a.x - b.x,
			a.y - b.y
			);
  }
  VECATTR  void operator -=(double2 &a, const double2 &b){
    a.x -= b.x;
    a.y -= b.y;
  }
  VECATTR  double2 operator -(const double2 &a, const double &b){
    return make_double2(
			a.x - b,
			a.y - b
			);
  }
  VECATTR  double2 operator -(const double &b, const double2 &a){
    return make_double2(
			b - a.x,
			b - a.y
			);
  }
  VECATTR  void operator -=(double2 &a, const double &b){a.x -= b; a.y -= b;}

  VECATTR  double2 operator *(const double2 &a, const double2 &b){
    return make_double2(a.x * b.x, a.y * b.y);
  }
  VECATTR  void operator *=(double2 &a, const double2 &b){
    a.x *= b.x;
    a.y *= b.y;
  }
  VECATTR  double2 operator *(const double2 &a, const double &b){
    return make_double2(a.x * b, a.y * b);
  }
  VECATTR  double2 operator *(const double &b, const double2 &a){
    return a*b;
  }
  VECATTR  void operator *=(double2 &a, const double &b){
    a.x *= b;
    a.y *= b;
  }

  VECATTR  double2 floorf(const double2 &a){
    return make_double2(floor(a.x), floor(a.y));
  }

  ////////////////////////////

  VECATTR double3 floorf(double3 v){return make_double3(floor(v.x), floor(v.y), floor(v.z));}



  VECATTR double dot(const double3 &a, const double3 &b){return a.x * b.x + a.y * b.y + a.z * b.z;}
  VECATTR float dot(const float3 &a, const float3 &b){return a.x * b.x + a.y * b.y + a.z * b.z;}
  VECATTR int dot(const int3 &a, const int3 &b){return a.x * b.x + a.y * b.y + a.z * b.z;}


  VECATTR double length(double3 v){return sqrt(dot(v, v));}
  VECATTR double3 normalize(double3 v)
  {
    double invLen = 1.0/sqrt(dot(v, v));
    return v * invLen;
  }

  VECATTR double3 cross(double3 a, double3 b){
    return make_double3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
  }


  //////////////////////////////////////////////////////////


  /****************************************************************************************/


  ///////////INT3/////////////////

  VECATTR int3 make_int3(int a){return make_int3(a,a,a);}

  VECATTR int3 operator /(int3 a, int3 b){
    return make_int3( a.x/b.x, a.y/b.y, a.z/b.z);
  }

  VECATTR  int3 operator +(const int3 &a, const int3 &b){
    return make_int3(a.x + b.x, a.y + b.y, a.z + b.z);
  }

  VECATTR  int3 operator +(const int3 &a, const int &b){
    return make_int3(a.x + b, a.y + b, a.z + b);
  }

  VECATTR  int3 operator *(const int3 &a, const int &b){return make_int3(a.x*b, a.y*b, a.z*b);}
  VECATTR  int3 operator *(const int &b, const int3 &a){return a*b;}

#endif
