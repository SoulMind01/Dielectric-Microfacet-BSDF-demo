#pragma once
#include "common.h"
#include "math_helpers.h"
#include "sampler.h"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace muni
{
  /**
   * Check if a value is falsy and print a message if it is.
   */
  bool check(Vec3f x, std::string name)
  {
    if (x.x == 0 || x.y == 0 || x.z == 0)
    {
      printf("%s is falsy value: %f %f %f\n", name.c_str(), x.x, x.y, x.z);
      return true;
    }
    else if (std::isnan(x.x) || std::isnan(x.y) || std::isnan(x.z))
    {
      printf("%s is NaN: %f %f %f\n", name.c_str(), x.x, x.y, x.z);
      return true;
    }
    return false;
  }
  bool check(float x, std::string name)
  {
    if (x == 0)
    {
      printf("%s is falsy value: %f\n", name.c_str(), x);
      return true;
    }
    else if (std::isnan(x))
    {
      printf("%s is NaN: %f\n", name.c_str(), x);
      return true;
    }
    return false;
  }
  float tot_G = 0, tot_D = 0, cnt = 0, cnt1 = 0;
  float tot_F = 0;
  float tot_h = 0;
  struct Lambertian
  {
    Vec3f albedo;

    /** Evaluates the BRDF for the Lambertian material.
      \return The BRDF (fr) value.
  */
    Vec3f eval() const
    {
      // =============================================================================================
      // TODO: Implement this function
      return albedo / M_PI;
      // =============================================================================================
    }

    /** Samples the BRDF for the Lambertian material.
      \param[in] normal The normal of the surface.
      \param[in] u A random number in (0,1)^2.
      \return A tuple containing the sampled direction in world space and the PDF.
    */
    std::tuple<Vec3f, float> sample(Vec3f normal, Vec2f u) const
    {
      // =============================================================================================
      // TODO: Implement this function
      // =============================================================================================
      float u1 = u.x, u2 = u.y;
      float z = u1;
      float phi = 2 * M_PI * u2;
      float r = sqrt(1 - z * z);
      float x = r * cos(phi);
      float y = r * sin(phi);
      Vec3f wi = Vec3f{x, y, z};
      wi = from_local(wi, normal);
      return std::make_tuple(wi, 1.0f / (2 * M_PI));
    }
    /** Computes the PDF for the Lambertian material.
      \param[in] wo The outgoing direction in world space.
      \param[in] wi The light incident direction in world space.
      \return The PDF value.
    */
    float pdf(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const
    {
      // =============================================================================================
      // TODO: Implement this function
      return 1.0f / 2 * M_PI;
      // =============================================================================================
    }
  };

  struct Microfacet
  {
    const float roughness; // roughness parameter
    // refraction indices for RGB channels
    float n1;
    float n2;

    /** Computes the Fresnel term for the microfacet material.
      \param[in] wi The light incident direction in local space.
      \return The Fresnel term.
    */
    float F(Vec3f wi) const
    {
      // =============================================================================================
      // TODO: Implement this function
      float R0 = (n1 - n2) * (n1 - n2) / ((n1 + n2) * (n1 + n2));
      float cos_theta = wi.z;
      float ans_F = R0 + (1.0f - R0) * (float)pow(1.0f - cos_theta, 5.0f);
      return ans_F;
      // =============================================================================================
    }
    /** Computes the Beckmann normal distribution function for the microfacet material.
      \param[in] h The half vector in local space.
      \return The normal distribution function.
    */
    float D(Vec3f h) const
    {
      // =============================================================================================
      // TODO: Implement this function
      tot_h += h.z;
      cnt1++;
      float theta_h = acos(h.z);
      float ans_D = exp(-tan(theta_h) * tan(theta_h) / (roughness * roughness)) / (M_PI * roughness * roughness * pow(cos(theta_h), 4.0f));
      return ans_D;
      // =============================================================================================
    }

    /** Computes the shadowing-masking function for the microfacet material.
      \param[in] wo The outgoing direction in local space.
      \param[in] wi The light incident direction in local space.
      \return The shadowing-masking value.
    */
    float G(Vec3f wo, Vec3f wi) const
    {
      // =============================================================================================
      // TODO: Implement this function
      float G_wo = 1.0f / (1.0f + lambda(wo));
      float G_wi = 1.0f / (1.0f + lambda(wi));
      return G_wo * G_wi;
      // =============================================================================================
    }

    float lambda(Vec3f w) const
    {
      float a = 1.0f / (roughness * abs(w.x / w.z)); //?
      if (a < 1.6f)
      {
        return (1 - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a);
      }
      else
      {
        return 0;
      }
    }

    /** Evaluates the BRDF for the microfacet material.
      \param[in] wo_world The outgoing direction in world space.
      \param[in] wi_world The light incident direction in world space.
      \param[in] normal The normal of the surface.
      \return The BRDF (fr) value.
    */
    Vec3f eval(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const
    {
      Vec3f wo = to_local(wo_world, normal);
      Vec3f wi = to_local(wi_world, normal);
      // =============================================================================================
      // TODO: Implement this function
      Vec3f h = normalize(wo + wi);
      float ans_F = F(wi);
      float ans_G = G(wo, wi);
      float ans_D = D(h);
      // debug
      tot_F = tot_F + ans_F;
      tot_G += ans_G;
      tot_D += ans_D;
      cnt++;
      float ans = ans_F * ans_G * ans_D / (4.0f * dot(wo, Vec3f{0.0f, 0.0f, 1.0f}) * dot(wi, Vec3f{0.0f, 0.0f, 1.0f}));
      return ans * Vec3f{1.0f, 1.0f, 1.0f};
      // =============================================================================================
    }
    /** Computes the PDF for the microfacet material.
      \param[in] wo The outgoing direction in world space.
      \param[in] wi The light incident direction in world space.
      \param[in] normal The normal of the surface.
      \return The PDF value.
    */
    float pdf(Vec3f wo_world, Vec3f wi_world, Vec3f normal, bool debug = false) const
    {
      // =============================================================================================
      // TODO: Implement this function
      Vec3f wh = normalize(wo_world + wi_world);
      float theta_h = acos(wh.z);
      float phi_h = atan2(wh.y, wh.x);
      float p_phi = 1.0f / (2.0f * M_PI);
      float p_theta = 2.0f * sin(theta_h) * exp(-tan(theta_h) * tan(theta_h) / (roughness * roughness)) / (roughness * roughness * pow(cos(theta_h), 3.0f));
      float p_wh = p_phi * p_theta / sin(theta_h);
      float p_wi = p_wh / (4.0f * dot(wo_world, wh));
      return p_wi;
      // =============================================================================================
    }

    float pdf_theta(float theta_h) const
    {
      return 2 * sin(theta_h) * exp(-tan(theta_h) * tan(theta_h) / (roughness * roughness)) / (roughness * roughness * pow(cos(theta_h), 3));
    }

    /** Samples the BRDF for the microfacet material.
      \param[in] wo_world The outgoing direction in world space.
      \param[in] normal The normal of the surface.
      \param[in] u A random number in (0,1)^2.
      \return A tuple containing the sampled direction in world space and the PDF.
    */
    std::tuple<Vec3f, float> sample(Vec3f wo_world, Vec3f normal, Vec2f u, bool debug = false) const
    {
      // =============================================================================================
      // TODO: Implement this function
      float u1 = u.x, u2 = u.y;
      float phi_h = 2.0f * M_PI * u1;
      float p_phi_h = 1.0f / (2.0f * M_PI);
      float theta_h = atan(sqrt(-roughness * roughness * log(1.0f - u2)));
      float tan_theta_h = sin(theta_h) / cos(theta_h);
      float p_theta_h = 2.0 * sin(theta_h) * exp(-tan_theta_h * tan_theta_h / (roughness * roughness)) / (roughness * roughness * pow(cos(theta_h), 3));
      Vec3f wh = Vec3f{sin(theta_h) * cos(phi_h), sin(theta_h) * sin(phi_h), cos(theta_h)};
      Vec3f wh_world = from_local(wh, normal);
      Vec3f wi_world = mirror_reflect(wo_world, wh_world);
      if (dot(wi_world, normal) < 0.0f)
        wi_world = -wi_world;
      float p_omega_h = p_theta_h * p_phi_h / sin(theta_h);
      float p_omega_i = p_omega_h / (4.0f * dot(wo_world, wh_world));
      return {wi_world, p_omega_i};
      // =============================================================================================
    }
  };
  struct Dielectric_BSDF
  {
    const float roughness;
    const float n1;
    const float n2;
    bool check_vec_invalid(const Vec3f &v, std::string name) const
    {
      if (std::abs(v.x) > 1.01f || std::abs(v.y) > 1.01f || std::abs(v.z) > 1.01f)
      {
        printf("%s is invalid(abs > 1): %f %f %f\n", name.c_str(), v.x, v.y, v.z);
        return true;
      }
      if ((float)length_squared(v) > 1.01f)
      {
        printf("length_squared: %f\n", length_squared(v));
        printf("%s is invalid(length > 1): %f %f %f\n", name.c_str(), v.x, v.y, v.z);
        return true;
      }
      if (std::isnan(v.x) || std::isnan(v.y) || std::isnan(v.z))
      {
        printf("%s is NaN: %f %f %f\n", name.c_str(), v.x, v.y, v.z);
        return true;
      }
      return false;
    }
    float F(Vec3f wi, bool debug = false) const
    {
      float cos_theta = std::abs(wi.z);
      if (cos_theta == 1.0f)
      {
        return 0;
      }
      float eta_i = wi.z > 0.0f ? n1 : n2;
      float eta_t = wi.z > 0.0f ? n2 : n1;
      float sin_theta = sqrt(std::max(0.0f, 1.0f - cos_theta * cos_theta));
      float sin_theta_t = eta_i / eta_t * sin_theta;
      // Total internal reflection
      if (sin_theta_t >= 1.0f)
      {
        return 1.0f;
      }
      float cos_theta_t = sqrt(std::max(0.0f, 1.0f - sin_theta_t * sin_theta_t));
      float R_parl = (eta_t * cos_theta - eta_i * cos_theta_t) / (eta_t * cos_theta + eta_i * cos_theta_t);
      float R_perp = (eta_i * cos_theta - eta_t * cos_theta_t) / (eta_i * cos_theta + eta_t * cos_theta_t);
      // Debug
      if (debug)
      {
        printf("eta_i: %f\n", eta_i);
        printf("eta_t: %f\n", eta_t);
        printf("cos_theta: %f\n", cos_theta);
        printf("cos_theta_t: %f\n", cos_theta_t);
        printf("R_parl: %f\n", R_parl);
        printf("R_perp: %f\n", R_perp);
      }
      return 0.5f * (R_parl * R_parl + R_perp * R_perp);
    }
    float D(Vec3f h, bool debug = false) const
    {
      float theta_h = acos(h.z);
      float ans_D = exp(-tan(theta_h) * tan(theta_h) / (roughness * roughness)) / (M_PI * roughness * roughness * pow(cos(theta_h), 4.0f));
      return ans_D;
    }
    float lambda(Vec3f w, bool debug) const
    {
      float absTanTheta = std::abs(tan(acos(w.z)));
      if (debug)
      {
        check(absTanTheta, "absTanTheta");
        printf("w.z: %f\n", w.z);
        printf("acos(w.z): %f\n", acos(w.z));
        printf("absTanTheta: %f\n", absTanTheta);
      }
      if (std::isinf(absTanTheta))
      {
        return 0.0f;
      }
      // Compute alpha for direction w
      float a = 1.0f / (roughness * absTanTheta);
      if (a >= 1.6f)
      {
        return 0.0f;
      }
      return (1.0f - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a);
    }
    float G1(Vec3f w, bool debug) const
    {
      return 1.0f / (1.0f + lambda(w, debug));
    }
    float G(Vec3f wo, Vec3f wi, bool debug = false) const
    {
      return G1(wo, debug) * G1(wi, debug);
    }
    bool is_valid_refraction(Vec3f wo, Vec3f wi, Vec3f wh, float etai, float etao) const
    {
      if (dot(wo, wh) * dot(wi, wh) > 0.0f)
      {
        printf("wo, wi are in the same hemisphere\n");
        return false;
      }
      float cos_theta_o = dot(wo, wh);
      float cos_theta_i = dot(wi, wh);
      float sin_theta_o = sqrt(std::max(0.0f, 1.0f - cos_theta_o * cos_theta_o));
      float sin_theta_i = sqrt(std::max(0.0f, 1.0f - cos_theta_i * cos_theta_i));
      float diff = etai * sin_theta_o - etao * sin_theta_i;
      if (diff > 0.01f)
      {
        printf("diff: %f\n", diff);
        printf("does not satisfy Snell's law\n");
        return false;
      }
      return true;
    }
    Vec3f computeHalfVector(Vec3f wo, Vec3f wi, bool debug = false) const
    {
      Vec3f wh;
      if (wi.z * wo.z >= 0.0f)
      {
        if (debug)
          printf("compute a half vector for a reflected ray\n");
        // Reflected
        wh = wo + wi;
        if (wh.x != 0 || wh.y != 0 || wh.z != 0)
        {
          if (debug)
            printf("wh is not 0\n");
          wh = normalize(wh);
        }
        else
        {
          float phi_o = atan2(wo.y, wo.x);
          float phi_i = atan2(wi.y, wi.x);
          float theta_o = acos(wo.z);
          float theta_i = acos(wi.z);
          if (theta_i + theta_o - M_PI < 0.001f)
          {
            if (debug)
              printf("wo, wi are parallel\n");
            wh = Vec3f{0.0f, 0.0f, 1.0f};
            return wh;
          }
          else
          {
            float phi_h = (phi_o + phi_i) / 2.0f;
            float theta_h = theta_o + theta_i == 0 ? 0.5f * M_PI : (theta_o + theta_i) / 2.0f;
            if (theta_h == 0.5f * M_PI)
            {
              wh = Vec3f{0.0f, 0.0f, 1.0f};
            }
            else
            {
              wh = Vec3f{sin(theta_h) * cos(phi_h), sin(theta_h) * sin(phi_h), cos(theta_h)};
              wh = normalize(wh);
            }
          }
        }
      }
      else
      { // Refracted
        if (debug)
          printf("compute a half vector for a refracted ray\n");
        float eta_i = wo.z > 0.0f ? n1 : n2;
        float eta_t = wo.z > 0.0f ? n2 : n1;
        wh = -normalize(wo * eta_t + eta_i * wi);
        if (std::isnan(wh.x) || std::isnan(wh.y) || std::isnan(wh.z))
        {
          printf("wo: %f %f %f\n", wo.x, wo.y, wo.z);
          printf("wi: %f %f %f\n", wi.x, wi.y, wi.z);
          printf("eta_i: %f\n", eta_i);
          printf("eta_t: %f\n", eta_t);
          printf("wh: %f %f %f\n", wh.x, wh.y, wh.z);
          throw std::runtime_error("wh is NaN");
        }
        if (debug)
        {
          printf("(compute)wo: %f %f %f\n", wo.x, wo.y, wo.z);
          printf("(compute)wi: %f %f %f\n", wi.x, wi.y, wi.z);
          printf("(compute)wh: %f %f %f\n", wh.x, wh.y, wh.z);
          // if (!is_valid_refraction(wo, wi, wh, eta_i, eta_t))
          // {
          //   throw std::runtime_error("wh is invalid");
          // }
        }
      }
      return wh;
    }

    Vec3f eval(Vec3f wo_world, Vec3f wi_world, Vec3f normal, bool debug = 0) const
    {
      // New method
      Vec3f wo = to_local(wo_world, normal);
      Vec3f wi = to_local(wi_world, normal);
      if (wi.x > 1.01f || wi.y > 1.01f || wi.z > 1.01f)
      {
        printf("wi_world: %f %f %f\n", wi_world.x, wi_world.y, wi_world.z);
        printf("wi: %f %f %f\n", wi.x, wi.y, wi.z);
        throw std::runtime_error("wi is invalid");
      }
      Vec3f wh = computeHalfVector(wo, wi);
      float F = this->F(to_local(wo, wh));
      float D = this->D(wh);
      float G = this->G(wo, wi);
      float fr = F * D * G / (4.0f * std::abs(wi.z * wo.z));
      float eta_i = wo.z > 0.0f ? n1 : n2;
      float eta_t = wo.z > 0.0f ? n2 : n1;
      float eta = eta_i / eta_t;
      // Modify: change dot(wi, normal) to dot(wi_world, normal) since wi is in local space
      float term1 = dot(wi, wh) * dot(wo, wh) / (dot(wi_world, normal) * dot(wo_world, normal));
      term1 = std::abs(term1);
      // Modify: term2 is used to normalize the result but wo and wi are already normalized, so term2 is not needed
      float term2 = (float)pow(dot(wh, wo) + eta * dot(wh, wi), 2);
      float ft = term1 * eta * eta * (1.0f - F) * D * G / term2;
      ft = 1.0f;
      float res = ft;
      if (debug && fr < 0.01f)
      {
        printf("F, D, G: %f %f %f\n", F, D, G);
        printf("dot(wh, normal): %f\n", wh.z);
        printf("fr: %f\n", fr);
        // inspect each term
        // printf("normal: %f %f %f\n", normal.x, normal.y, normal.z);
        // printf("wo: %f %f %f\n", wo.x, wo.y, wo.z);
        // printf("wi: %f %f %f\n", wi.x, wi.y, wi.z);
        // printf("wh: %f %f %f\n", wh.x, wh.y, wh.z);
        // printf("term1: %f\n", term1);
        // printf("dot(wi, wh): %f\n", dot(wi, wh));
        // printf("dot(wo, wh): %f\n", dot(wo, wh));
        // printf("dot(wi_world, normal): %f\n", dot(wi_world, normal));
        // printf("dot(wo_world, normal): %f\n", dot(wo_world, normal));
        // printf("eta: %f\n", eta);
        // printf("wo + eta * wi: %f %f %f\n", wo.x + eta * wi.x, wo.y + eta * wi.y, wo.z + eta * wi.z);
        // printf("term2: %f\n", term2);
        // printf("(1 - F) * D * G: %f\n", (1.0f - F) * D * G);
        // printf("F: %f\n", F);
        // check(F, "F");
        // check(D, "D");
        // if (check(G, "G"))
        //   this->G(wo, wi, true);
        // printf("res: %f\n", res);
        // throw std::runtime_error("res is invalid");
      }
      // if (debug)
      // {
      //   printf("fr: %f\n", fr);
      //   printf("ft: %f\n", ft);
      // }
      return (fr + ft) * Vec3f{1.0f, 1.0f, 1.0f};

      // Vec3f wo = to_local(wo_world, normal);
      // Vec3f wi = to_local(wi_world, normal);
      // // Judge if the ray is reflected or refracted
      // // If the ray is reflected, they should be in the same hemisphere
      // // If the ray is refracted, they should be in different hemisphere
      // bool is_reflected = wo.z * wi.z > 0.0f;
      // Vec3f wh = wo + wi;
      // if (!is_reflected)
      // {
      //   float eta_i = wo.z > 0.0f ? n1 : n2;
      //   float eta_t = wo.z > 0.0f ? n2 : n1;
      //   wh = -normalize(wo + eta_i / eta_t * wi);
      // }
      // // Handle the degenerate case
      // if (wh.x == 0 && wh.y == 0 && wh.z == 0)
      // {
      //   // In this case, wo and wi are parallel, so assign wh to the perpendicular vector
      //   if (wo.x != 0 || wo.y != 0)
      //   {
      //     wh = Vec3f{-wo.y, wo.x, 0};
      //   }
      //   else
      //   {
      //     wh = Vec3f{0, -wo.z, wo.y};
      //   }
      //   wh = normalize(wh);
      // }
      // // Debug
      // if (debug)
      // {
      //   // check wo, wi, wh
      //   check(wo, "wo");
      //   check(wi, "wi");
      //   if (check(wh, "wh"))
      //     printf("wo, wi, wh: %f %f %f, %f %f %f, %f %f %f\n", wo.x, wo.y, wo.z, wi.x, wi.y, wi.z, wh.x, wh.y, wh.z);
      // }
      // float F = this->F(wo);
      // float D = this->D(wh, debug);
      // float G = this->G(wo, wi);
      // // float cos_theta_o = std::abs(wo.z);
      // // float cos_theta_i = std::abs(wi.z);
      // float cos_theta_o = wo.z;
      // float cos_theta_i = wi.z;
      // float fr = F * D * G / std::abs(4.0f * cos_theta_o * cos_theta_i);
      // bool entering = wo.z > 0.0f;
      // float eta_i = entering ? n1 : n2;
      // float eta_t = entering ? n2 : n1;
      // float eta = eta_i / eta_t;
      // wh = computeHalfVector(wo, wi, debug);
      // // Debug
      // // float ft = eta_t * eta_t * (1.0f - F) * D * G / (float)pow(eta_i * dot(wo, wh) + eta_t * dot(wi, wh), 2.0f);
      // float ft = this->D(wh) * G * (1.0f - F) / (float)pow(dot(wo, wh) + eta * dot(wi, wh), 2.0f) * dot(wi, wh) * dot(wo, wh) / (cos_theta_o * cos_theta_i);
      // ft = std::abs(ft);

      // if (debug)
      // {
      //   printf("fr, ft: %f %f\n", fr, ft);
      //   // Check ft recursively
      //   if (ft == 0)
      //   {
      //     printf("ft is 0\n");
      //     // Check D, G, F
      //     check(this->D(wh), "D");
      //     check(G, "G");
      //     check(1.0 - F, "1 - F");
      //     float term1 = pow(dot(wo, wh) + eta * dot(wi, wh), 2.0f);
      //     float term2 = std::abs(dot(wi, wh)) * dot(wo, wh);
      //     float term3 = cos_theta_o * cos_theta_i;
      //     check(term1, "term1");
      //     check(term2, "term2");
      //     check(term3, "term3");
      //     printf("ft = D * G * (1 - F) / term1 * term2 / term3\n");
      //     printf("term1, term2, term3: %f %f %f\n", term1, term2, term3);
      //     printf("D, G, 1 - F: %f %f %f\n", this->D(wh), G, 1.0 - F);
      //   }
      // }
      // return (fr + ft) * Vec3f{1.0f, 1.0f, 1.0f};
    }
    float pdf(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const
    {
      // Debug: sample a uniform whole sphere
      return 1.0f / (4.0f * M_PI);
      // Vec3f wo = to_local(wo_world, normal);
      // Vec3f wi = to_local(wi_world, normal);
      // float eta_i = n1, eta_o = n2;
      // if (dot(wo, normal) < 0.0f)
      // {
      //   std::swap(eta_i, eta_o);
      // }
      // if (dot(wo, normal) * dot(wi, normal) > 0.0f)
      // {
      //   // Reflect
      //   Vec3f h = normalize(wo + wi);
      //   float weight = std::abs(dot(wo, h) * G(wo, wi) / dot(wo, normal) / dot(h, normal));
      //   float p_wi = weight != 0 ? eval(wo_world, wi_world, normal).x / weight
      //                            : 100;
      //   return p_wi;
      // }
      // else
      // {
      //   // Refract
      //   float eta_i = n1, eta_o = n2;
      //   if (dot(wo, normal) < 0.0f)
      //   {
      //     std::swap(eta_i, eta_o);
      //   }
      //   Vec3f h = -normalize(eta_i * wi + eta_o * wo);
      //   float weight = std::abs(dot(wo, h) * G(wo, wi) / dot(wo, normal) / dot(h, normal));
      //   float p_wi = weight != 0 ? eval(wo_world, wi_world, normal).x / weight
      //                            : 100;
      //   return p_wi;
      // }
    }

    /**
     * @brief Computes the refracted direction given an incident direction, surface normal, and the ratio of indices of refraction.
     *
     * @param wi The incident direction.
     * @param n The surface normal, in the same hemisphere as wi.
     * @param eta The ratio of indices of refraction in the incident and transmitted media, respectively.
     * @param wt The refracted direction, returned by this function.
     *
     * @return True if a valid refracted ray was returned in wt, false in the case of total internal reflection.
     */
    bool Refract(const Vec3f &wi, const Vec3f &n, float eta, Vec3f &wt, bool debug = false) const
    {
      if (check(wi, "wi"))
      {
        throw std::runtime_error("in refract function: wi is invalid");
      }
      float cosThetaI = wi.z;
      float sin2ThetaI = std::max(0.0f, 1.0f - cosThetaI * cosThetaI);
      float sin2ThetaT = eta * eta * sin2ThetaI;
      if (sin2ThetaT >= 1)
      {
        printf("cosThetaI: %f\n", cosThetaI);
        printf("sin2ThetaI: %f\n", sin2ThetaI);
        printf("eta: %f\n", eta);
        printf("sin2ThetaT: %f\n", sin2ThetaT);
        return false;
      }
      float cosThetaT = std::sqrt(1 - sin2ThetaT);
      wt = eta * -wi + (eta * cosThetaI - cosThetaT) * n;
      if (check(wt, "wt"))
      {
        printf("wi: %f %f %f\n", wi.x, wi.y, wi.z);
        printf("eta: %f\n", eta);
        printf("cosThetaI: %f\n", cosThetaI);
        printf("cosThetaT: %f\n", cosThetaT);
        printf("n: %f %f %f\n", n.x, n.y, n.z);
        throw std::runtime_error("in refract function: wt is invalid");
      }
      wt = normalize(wt);
      return true;
    }

    Vec3f Faceforward(const Vec3f &n, const Vec3f &v) const
    {
      return (dot(n, v) < 0.0f) ? -n : n;
    }

    std::tuple<Vec3f, float> sample(Vec3f wo_world, Vec3f normal, Vec2f u, bool debug = false) const
    {
      // Debug: sample a uniform whole sphere
      // float u1 = u.x, u2 = u.y;
      // float phi = 2 * M_PI * u1;
      // float theta = acos(1 - 2 * u2);
      // Vec3f wi = Vec3f{sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
      // wi = from_local(wi, normal);
      // return std::make_tuple(wi, 1.0f / (2 * M_PI));
      // Debug: sample a microfacet normal
      // float u1 = u.x, u2 = u.y;
      // float phi_h = 2.0f * M_PI * u1;
      // float theta_h = atan(sqrt(-roughness * roughness * log(1.0f - u2)));
      // Vec3f wh = Vec3f{sin(theta_h) * cos(phi_h), sin(theta_h) * sin(phi_h), cos(theta_h)};
      float logSample = std::log(1 - u.x);
      if (std::isinf(logSample))
        logSample = 0;
      float tan2Theta = -roughness * roughness * logSample;
      float phi = u.y * 2.0f * M_PI;
      float cos_theta_h = 1.0f / sqrt(1.0f + tan2Theta);
      float sin_theta_h = (float)sqrt(std::max(0.0f, 1.0f - cos_theta_h * cos_theta_h));
      float p_theta_h = 2.0f * sin_theta_h * exp(-tan2Theta / (roughness * roughness)) / (roughness * roughness * pow(cos_theta_h, 3));
      Vec3f wo = to_local(wo_world, normal);
      Vec3f wh = Vec3f{sin_theta_h * cos(phi), sin_theta_h * sin(phi), cos_theta_h};
      if (wh.z * wo.z < 0.0f)
      {
        wh = -wh;
      }

      wo = to_local(wo, wh);
      float p_wi = p_theta_h / (4.0f * std::abs(dot(wo, wh)));

      float Fresnel = F(wo);
      float tmp = UniformSampler::next1d();
      if (tmp <= Fresnel)
      {
        // Reflect
        Vec3f wi = mirror_reflect(-wo, wh);
        wi = from_local(wi, wh);
        wi = from_local(wi, normal);
        if (check_vec_invalid(wi, "wi"))
          throw std::runtime_error("(reflect) wi is invalid");
        return {wi, p_wi};
      }
      else
      {
        // Refract
        bool entering = wo.z > 0.0f;
        float eta_i = entering ? n1 : n2;
        float eta_t = entering ? n2 : n1;
        Vec3f wi;
        check_vec_invalid(wi, "wi");
        if (!Refract(wo, Faceforward(wh, wo), eta_i / eta_t, wi))
        {
          printf("tmp: %f\n", tmp);
          printf("Fresnel: %f\n", Fresnel);
          printf("wo: %f %f %f\n", wo.x, wo.y, wo.z);
          printf("normal: %f %f %f\n", normal.x, normal.y, normal.z);
          printf("wo * normal: %f\n", dot(wo, normal));
          throw std::runtime_error("Total internal reflection");
          return {Vec3f{0.0f, 0.0f, 0.0f}, 0.0f};
        }
        // Debug
        if (debug)
        {
          printf("wi: %f %f %f\n", wi.x, wi.y, wi.z);
          printf("wo: %f %f %f\n", wo.x, wo.y, wo.z);
          printf("normal: %f %f %f\n", normal.x, normal.y, normal.z);
        }
        check_vec_invalid(wi, "wi");
        wi = from_local(wi, wh);
        check_vec_invalid(wi, "wi");
        wi = from_local(wi, normal);
        if (check_vec_invalid(wi, "wi"))
          throw std::runtime_error("(refract) wi is invalid");
        return {wi, p_wi};
      }
    }
  };
}; // namespace muni
