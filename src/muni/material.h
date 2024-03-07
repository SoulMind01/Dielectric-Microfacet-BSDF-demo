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
    float F(Vec3f wi, Vec3f h) const
    {
      // pbr
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
      // float c = std::abs(dot(wi, h));
      // float eta_i = n1, eta_t = n2;
      // if (dot(wi, h) < 0.0f)
      // {
      //   std::swap(eta_i, eta_t);
      // }
      // float g = eta_t * eta_t / (eta_i * eta_i) - 1.0f + c * c;
      // if (g < 0.0f)
      // {
      //   return 1.0f;
      // }
      // g = sqrt(g);
      // float res = 0.5f * ((g - c) / (g + c)) * ((g - c) / (g + c)) * (1.0f + ((c * (g + c) - 1.0f) / (c * (g - c) + 1.0f)) * ((c * (g + c) - 1.0f) / (c * (g - c) + 1.0f)));
      // return res;
    }
    float D(Vec3f h) const
    {
      float x = h.z > 0.0f ? 1.0f : 0.0f;
      float theta_h = acos(h.z);
      float res = x * exp(-tan(theta_h) * tan(theta_h) / (roughness * roughness)) / (M_PI * roughness * roughness * pow(cos(theta_h), 4.0f));
      return res;
    }
    float G1(Vec3f w, Vec3f h) const
    {
      float x = dot(w, h) / w.z > 0.0f ? 1.0f : 0.0f;
      float theta_h = acos(h.z);
      float a = 1.0f / (roughness * tan(theta_h));
      if (a < 1.6f)
      {
        return x * (3.535f * a + 2.181f * a * a) / (1.0f + 2.276f * a + 2.577f * a * a);
      }
      return x;
    }
    float G(Vec3f wo, Vec3f wi, Vec3f h) const
    {
      return G1(wo, h) * G1(wi, h);
    }
    Vec3f eval(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const
    {
      Vec3f wo = to_local(wo_world, normal);
      Vec3f wi = to_local(wi_world, normal);
      Vec3f hr = normalize(wo + wi);
      float eta_i = n1, eta_o = n2;
      if (dot(wo, normal) < 0.0f)
      {
        std::swap(eta_i, eta_o);
      }
      Vec3f ht = -normalize(eta_i * wo + eta_o * wi);
      float fr = F(wo, hr) * G(wo, wi, hr) * D(hr) / (4.0f * std::abs(dot(wi, normal)) * std::abs(dot(wo, normal)));
      float ft = std::abs(dot(wo, ht) * dot(wi, ht) / dot(wo, normal) / dot(wi, normal)) * eta_o * eta_o * (1.0f - F(wi, ht)) * G(wo, wi, ht) * D(ht) / (float)pow(eta_i * dot(wo, ht) + eta_o * dot(wi, ht), 2.0f);
      return (fr + ft) * Vec3f{1.0f, 1.0f, 1.0f};
    }
    float pdf(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const
    {
      Vec3f wo = to_local(wo_world, normal);
      Vec3f wi = to_local(wi_world, normal);
      float eta_i = n1, eta_o = n2;
      if (dot(wo, normal) < 0.0f)
      {
        std::swap(eta_i, eta_o);
      }
      // return 1.0f / 2.0f / M_PI;
      if (dot(wo, normal) * dot(wi, normal) > 0.0f)
      {
        // Reflect
        Vec3f h = normalize(wo + wi);
        float weight = std::abs(dot(wo, h) * G(wo, wi, h) / dot(wo, normal) / dot(h, normal));
        float p_wi = weight != 0 ? eval(wo_world, wi_world, normal).x / weight
                                 : 100;
        return p_wi;
      }
      else
      {
        // Refract
        float eta_i = n1, eta_o = n2;
        if (dot(wo, normal) < 0.0f)
        {
          std::swap(eta_i, eta_o);
        }
        Vec3f h = -normalize(eta_i * wi + eta_o * wo);
        float weight = std::abs(dot(wo, h) * G(wo, wi, h) / dot(wo, normal) / dot(h, normal));
        float p_wi = weight != 0 ? eval(wo_world, wi_world, normal).x / weight
                                 : 100;
        return p_wi;
      }
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
      // // Debug: return a random direction on a sphere
      // float u1 = u.x, u2 = u.y;
      // float theta = acos(1 - 2 * u1);
      // float phi = 2 * M_PI * u2;
      // Vec3f wi = Vec3f{sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
      // wi = from_local(wi, normal);
      // return std::make_tuple(wi, 1.0f / (2 * M_PI));

      // Vec3f wo = to_local(wo_world, normal);
      // float theta_h = atan(sqrt(-roughness * roughness * log(1.0f - u.x)));
      // float phi_h = 2.0f * M_PI * u.y;
      // Vec3f h = Vec3f{sin(theta_h) * cos(phi_h), sin(theta_h) * sin(phi_h), cos(theta_h)};
      // // Evaluate the Fresnel term to determine reflection or refraction
      // float Fresnel = F(wo, h);
      // float eta_i = n1, eta_o = n2;
      // if (dot(wo, normal) < 0.0f)
      // {
      //   std::swap(eta_i, eta_o);
      // }
      // float c = dot(wo, h);
      // float temp = 1.0f + eta_i / eta_o * (c * c - 1.0f);
      // if (UniformSampler::next1d() < Fresnel || temp < 0)
      // {
      //   // Reflect
      //   Vec3f wi = mirror_reflect(-wo, h);
      //   float weight = std::abs(dot(wo, h) * G(wo, wi, h) / dot(wo, normal) / dot(h, normal));
      //   Vec3f wi_world = from_local(wi, normal);
      //   float p_wi = eval(wo_world, wi_world, normal).x / weight;
      //   return {wi_world, p_wi};
      // }
      // else
      // {
      //   // Refract
      //   float sign = dot(wo, normal) < 0.0f ? 1.0f : -1.0f;
      //   Vec3f wi = h * (eta_i / eta_o * c - sign * (float)sqrt(temp)) - eta_i / eta_o * wo;
      //   float weight = std::abs(dot(wo, h) * G(wo, wi, h) / dot(wo, normal) / dot(h, normal));
      //   Vec3f wi_world = from_local(wi, normal);
      //   float p_wi = eval(wo_world, wi_world, normal).x / weight;
      //   return {wi_world, p_wi};
      // }
      Vec3f wo = to_local(wo_world, normal);
      float Fresnel = F(wo, Vec3f{0.0f, 0.0f, 1.0f});
      if (UniformSampler::next1d() < Fresnel)
      {
        // Reflect
        Vec3f wi = Vec3f{-wo.x, -wo.y, wo.z};

        float p_wi = Fresnel;
        return {from_local(wi, normal), p_wi};
      }
      else
      {
        // Refract
        bool entering = wo.z > 0.0f;
        float eta_i = entering ? n1 : n2;
        float eta_t = entering ? n2 : n1;
        Vec3f wi;
        if (!Refract(wo, Faceforward(Vec3f{0.0f, 0.0f, 1.0f}, wo), eta_i / eta_t, wi))
        {
          return {wi, 0.01f};
        }
        float p_wi = 1.0f - Fresnel;
        return {from_local(wi, normal), p_wi};
      }
    }
  };
}; // namespace muni
