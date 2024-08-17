#pragma once
#include <vector>
#include <memory>
#include <cmath>
#include <execution>
#include <algorithm>

#define MY_BLACK       float3(0.,0.,0.)   
#define MY_RED         float3(0.93, 0   , 0.01)
#define MY_LIGHT_RED   float3(1,  0.5 , 0.5)
#define MY_ORANGE      float3(1   , 0.45, 0.21)
#define MY_YELLOW      float3(0.9 , 0.9 , 0   )
#define MY_LIGHT_YELLOW float3(1 , 1 , 0.8 )  
#define MY_LIGHT_ROSE  float3(1 , 0.8 , 0.8 )  
#define MY_DARK_YELLOW float3(0.45 , 0.45, 0   )
#define MY_GREEN       float3(0.01, 0.98, 0.01)
#define MY_BLUE        float3(0   , 0.1 , 1   )
#define MY_LIGHT_BLUE  float3(0.5 , 0.75 , 1   )
#define MY_INDIGO      float3(0.3 , 0   , 0.53)
#define MY_VIOLET      float3(0.58, 0   , 0.83)
#define MY_DARK_RED    float3(0.5, 0   , 0.01)
#define MY_ROSE        float3(1., 0.  , 1.)
#define MY_PURPLE      float3(0.6, 0  , 0.6)
#define MY_DARK_GREEN  float3(0.01, 0.5, 0.01)
#define MY_DARK_BLUE   float3(0   , 0.1 , 0.5   )
#define MY_GOLD        float3(0.864 , 0.7 , 0.325 )
#define MY_DARK_GOLD   float3(0.4,    0.3 , 0.16 )
#define MY_SILVER      float3(0.75,0.75,0.75)
#define MY_PEATCH      float3(1, 0.89, 0.705)
#define MY_BRONZE      float3(0.8, 0.5, 0.2)
#define MY_GRAY        float3(0.5,0.5,0.5) 
#define MY_LIGHT_GRAY  float3(0.8,0.8,0.8) 
#define MY_DARK_GRAY   float3(0.2,0.2,0.2) 
#define MY_WHITE       float3(1.,1.,1.)
#define MY_DARK_CHERRY float3(0.57,0.11,0.26)
#define MY_LIGHT_CHERRY float3(0.87,0.2,0.4)

static const float gravity = 1;
static const float M_PI = 3.1415926535897932384;