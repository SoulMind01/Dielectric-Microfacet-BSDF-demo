#pragma once
#include "common.h"
#include "math_helpers.h"
#include "sampler.h"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace muni
{

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
  Microfacet tmp = Microfacet{0.1f, 1.5f, 1.0f};
  struct Dielectric_BSDF
  {
    const float roughness;
    const float n1;
    const float n2;
    float F(Vec3f wi) const
    {
      float cos_theta = acos(wi.z);
      bool entering = cos_theta > 0.0f;
      float eta_i = n1, eta_t = n2;
      if (!entering)
      {
        std::swap(eta_i, eta_t);
      }
      float sin_theta = sqrt(std::max(0.0f, 1.0f - cos_theta * cos_theta));
      float sin_theta_t = eta_i / eta_t * sin_theta;
      float cos_theta_t = sqrt(std::max(0.0f, 1.0f - sin_theta_t * sin_theta_t));
      float R_parl = (eta_t * cos_theta - eta_i * cos_theta_t) / (eta_t * cos_theta + eta_i * cos_theta_t);
      float R_perp = (eta_i * cos_theta - eta_t * cos_theta_t) / (eta_i * cos_theta + eta_t * cos_theta_t);
      return 0.5f * (R_parl * R_parl + R_perp * R_perp);
    }
    float D(Vec3f h) const
    {
      float theta_h = acos(h.z);
      float tan_theta_h = tan(theta_h);
      float cos_theta = cos(theta_h);
      float res = exp(-tan_theta_h * tan_theta_h / (roughness * roughness)) / (M_PI * roughness * roughness * pow(cos_theta, 4.0f));
      return res;
    }
    float lambda(Vec3f w) const
    {
      float absTanTheta = std::abs(tan(acos(w.z)));
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
    float G1(Vec3f w) const
    {
      return 1.0f / (1.0f + lambda(w));
    }
    float G(Vec3f wo, Vec3f wi) const
    {
      return G1(wo) * G1(wi);
    }
    Vec3f eval(Vec3f wo_world, Vec3f wi_world, Vec3f normal, bool debug = 0) const
    {
      Vec3f wo = to_local(wo_world, normal);
      Vec3f wi = to_local(wi_world, normal);
      Vec3f wh = normalize(wo + wi);
      float F = this->F(wi);
      float D = this->D(wh);
      float G = this->G(wo, wi);
      float cos_theta_o = std::abs(wo.z);
      float cos_theta_i = std::abs(wi.z);
      float fr = F * D * G / (4.0f * cos_theta_o * cos_theta_i);
      bool entering = wo.z > 0.0f;
      float eta_i = entering ? n1 : n2;
      float eta_t = entering ? n2 : n1;
      float ft = eta_t * eta_t * (1.0f - F) * D * G / pow(eta_i * dot(wo, wh) + eta_t * dot(wi, wh), 2.0f);
      // Debug
      if (debug)
        ft = 0.5f;
      return (fr + ft) * Vec3f{1.0f, 1.0f, 1.0f};
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
    bool Refract(const Vec3f &wi, const Vec3f &n, float eta, Vec3f &wt) const
    {
      float cosThetaI = dot(n, wi);
      float sin2ThetaI = std::max(0.0f, 1.0f - cosThetaI * cosThetaI);
      float sin2ThetaT = eta * eta * sin2ThetaI;
      if (sin2ThetaT >= 1)
      {
        return false;
      }
      float cosThetaT = std::sqrt(1 - sin2ThetaT);
      wt = eta * -wi + (eta * cosThetaI - cosThetaT) * n;
      return true;
    }

    Vec3f Faceforward(const Vec3f &n, const Vec3f &v) const
    {
      return (dot(n, v) < 0.0f) ? -n : n;
    }

    std::tuple<Vec3f, float> sample(Vec3f wo_world, Vec3f normal, Vec2f u) const
    {
      // Debug: sample a uniform whole sphere
      float u1 = u.x, u2 = u.y;
      float phi = 2 * M_PI * u1;
      float theta = acos(1 - 2 * u2);
      Vec3f wi = Vec3f{sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
      wi = from_local(wi, normal);
      return std::make_tuple(wi, 1.0f / (2 * M_PI));
      // Vec3f wo = to_local(wo_world, normal);
      // float Fresnel = F(wo);
      // if (UniformSampler::next1d() < Fresnel)
      // {
      //   // Reflect
      //   Vec3f wi = Vec3f{-wo.x, -wo.y, wo.z};
      //   float p_wi = Fresnel;
      //   return {from_local(wi, normal), p_wi};
      // }
      // else
      // {
      //   // Refract
      //   bool entering = wo.z > 0.0f;
      //   float eta_i = entering ? n1 : n2;
      //   float eta_t = entering ? n2 : n1;
      //   Vec3f wi;
      //   if (!Refract(wo, Faceforward(Vec3f{0.0f, 0.0f, 1.0f}, wo), eta_i / eta_t, wi))
      //   {
      //     return {wi, 0.01f};
      //   }
      //   float p_wi = 1.0f - Fresnel;
      //   return {from_local(wi, normal), p_wi};
      // }
    }
  };
}; // namespace muni
