#pragma once
#include "common.h"
#include "math_helpers.h"
#include <algorithm>
#include <cmath>
#include <iostream>
namespace muni
{
  struct Dielectric
  {
    Vec3f eval() const {}
    float F(Vec3f wi, Vec3f h)
    {
      float c = std::abs(dot(wi, h));
        }
  };
}