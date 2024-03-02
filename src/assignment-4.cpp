#include "material.h"
#include "muni/camera.h"
#include "muni/common.h"
#include "muni/image.h"
#include "muni/material.h"
#include "muni/math_helpers.h"
#include "muni/obj_loader.h"
#include "muni/ray_tracer.h"
#include "muni/sampler.h"
#include "muni/scenes/box.h"
#include "muni/triangle.h"
#include "ray_tracer.h"
#include "spdlog/spdlog.h"
#include "triangle.h"
#include <cmath>
#include <iostream>

using namespace muni;

RayTracer::Octree octree{};

/** Offset the ray origin to avoid self-intersection.
    \param[in] ray_pos The original ray origin.
    \param[in] normal The normal of the surface at the hit point.
    \return The offset ray origin.
*/
Vec3f offset_ray_origin(Vec3f ray_pos, Vec3f normal)
{
    return ray_pos + EPS * normal;
}

/** Check if the triangle is an emitter.
    \param[in] tri The triangle to check
    \return True if the triangle is an emitter, false otherwise.
*/
bool is_emitter(const Triangle &tri) { return tri.emission != Vec3f{0.0f}; }

/** Evaluate the radiance of the area light. We **do not** check whether the hit
 point is on the light source, so make sure
 *  the hit point is on the light source before calling this function.
    \param[in] light_dir The **outgoing** direction from the light source to the
 scene. \return The radiance of the light source.
*/
Vec3f eval_area_light(const Vec3f light_dir)
{
    if (dot(light_dir, BoxScene::light_normal) > 0.0f)
        return BoxScene::light_color;
    return Vec3f{0.0f};
}

/** Sample a point on the area light with a uniform distribution.
    \param[in] samples A 2D uniform random sample.
    \return A tuple containing the sampled position, the normal of the light
 source, and the L_BRDF value.
*/
std::tuple<Vec3f, Vec3f, float> sample_area_light(Vec2f samples)
{
    // =============================================================================================
    // TODO: Implement this function
    // =============================================================================================
    // Sample a point on the area light with a uniform distribution
    float x = muni::BoxScene::light_x, y = muni::BoxScene::light_y,
          z = muni::BoxScene::light_z;
    Vec3f p = Vec3f{x + muni::BoxScene::light_len_x * samples.x,
                    y + muni::BoxScene::light_len_y * samples.y, z};
    Vec3f n = muni::BoxScene::light_normal;
    float L_BRDF =
        1.0f / (muni::BoxScene::light_len_x * muni::BoxScene::light_len_y);
    return std::make_tuple(p, n, L_BRDF);
}

Vec3f shade_with_light_sampling(Triangle tri, Vec3f p, Vec3f wo)
{
    // =============================================================================================
    // TODO: Implement this function
    // Please refer to lecture 9, page 20&21 for the details of the
    // implementation. spdlog::warn(
    //     "Part 2 is not implemented yet.");  // Remove this line after you
    // finish the implementation.

    /*
        shade(p, wo)
        // Contribution from the light source.
        Uniformly sample the light at x’ (L_light(x’) = 1 / A)
        L_dir = L_i * f_r * cos θ / L_light(ω from p to x’)
        // Contribution from other reflectors.
        L_indir = 0.0
        Test Russian Roulette with probability P_RR
        Randomly sample the hemisphere toward ωi (pdf_brdf)
        Trace a ray r(p, ωi)
        If ray r hit a non-emitting object at q
        L_indir = shade(q, -ωi) * f_r * cos θ / pdf_brdf / P_RR
        Return L_dir + L_indir
    */

    bool flag = std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id]);
    // Contribution from the light source
    Vec3f L_dir{0.0f};
    // Uniformly sample the light at x
    auto [x, n, L_light] = sample_area_light(UniformSampler::next2d());
    L_light = length_squared(x - p) * L_light / dot(n, normalize(p - x));
    // Avoid self-intersection
    p = offset_ray_origin(p, tri.face_normal);
    // Shoot a ray from p to x
    Vec3f wi = normalize(x - p);
    // If the ray is not blocked in the middle
    float max_t = length(x - p);
    bool hit = RayTracer::any_hit(p, wi, max_t, octree, BoxScene::triangles);
    if (!hit)
    {
        Vec3f f_r;
        if (flag)
        {
            f_r = std::get<Lambertian>(BoxScene::materials[tri.material_id]).eval();
        }
        else
        {
            f_r = std::get<Microfacet>(BoxScene::materials[tri.material_id]).eval(wo, wi, tri.face_normal);
        }
        float cosine = dot(wi, tri.face_normal);
        if (cosine < 0.0f)
            cosine = 0.0f;
        L_dir = eval_area_light(-wi) * f_r * cosine / L_light;
        // printf("L_dir: %f %f %f\n", L_dir.x, L_dir.y, L_dir.z);
    }
    // Contribution from other reflectors
    Vec3f L_indir{0.0f};
    // Test Russian Roulette with probability p_rr = 0.8f
    const float p_rr = 0.8f;
    float ksi = UniformSampler::next1d();
    if (ksi > p_rr)
        return L_dir + L_indir;
    // Randomly choose one direction wi with uniform hemisphere sampling
    /*
        auto lambertian = std::get<Lambertian>(BoxScene::materials[tri.material_id]);
        auto [wj, pdf_brdf] = lambertian.sample(tri.face_normal, UniformSampler::next2d());
    */

    Vec3f wj;
    float pdf_brdf;
    if (flag)
    {
        auto lambertian = std::get<Lambertian>(BoxScene::materials[tri.material_id]);
        std::tie(wj, pdf_brdf) = lambertian.sample(tri.face_normal, UniformSampler::next2d());
    }
    else
    {
        auto microfacet = std::get<Microfacet>(BoxScene::materials[tri.material_id]);
        std::tie(wj, pdf_brdf) = microfacet.sample(wo, tri.face_normal, UniformSampler::next2d());
    }
    // Trace the new ray
    auto [hit2, t_min2, nearest_tri2] =
        RayTracer::closest_hit(p, wj, octree, BoxScene::triangles);
    // If the ray hit a non-emitting object at q
    if (hit2 && !is_emitter(nearest_tri2))
    {
        Vec3f q = p + t_min2 * wj;
        // Vec3f f_r =
        //     (std::get<Lambertian>(BoxScene::materials[tri.material_id])).eval();
        Vec3f f_r;
        if (flag)
        {
            f_r = std::get<Lambertian>(BoxScene::materials[tri.material_id]).eval();
        }
        else
        {
            f_r = std::get<Microfacet>(BoxScene::materials[tri.material_id]).eval(wo, wj, tri.face_normal);
        }
        float cosine = dot(wj, tri.face_normal);
        if (cosine < 0.0f)
            cosine = 0.0f;
        L_indir = shade_with_light_sampling(nearest_tri2, q, -wj) * f_r * cosine /
                  pdf_brdf / p_rr;
    }
    return L_dir + L_indir;
    // =============================================================================================
}

Vec3f path_tracing_with_light_sampling(Vec3f ray_pos, Vec3f ray_dir)
{
    const auto [is_ray_hit, t_min, nearest_tri] =
        RayTracer::closest_hit(ray_pos, ray_dir, octree, BoxScene::triangles);
    if (!is_ray_hit)
        return Vec3f{0.0f};
    const Vec3f hit_position = ray_pos + t_min * ray_dir;
    if (is_emitter(nearest_tri))
        return eval_area_light(-ray_dir);
    return shade_with_light_sampling(nearest_tri, hit_position, -ray_dir);
}

std::tuple<Vec3f, float> sample_uniform_hemisphere(Vec3f normal, Vec2f u)
{
    // =============================================================================================
    // TODO: Implement this function
    // =============================================================================================
    // Sample a direction on the hemisphere with a uniform distribution
    // u1 and u2 are two uniform random numbers in the range (0, 1)
    // z = u1, phi = 2 * PI * u2, r = sqrt(1 - z * z)
    // x = r * cos(phi), y = r * sin(phi)
    // The sampled direction is (x, y, z)
    // The L_BRDF value is 1 / (2 * PI)
    // the original normal is (0, 0, 1), so we need to transform the sampled
    // direction to the world space
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

Vec3f shade_with_MIS(Triangle tri, Vec3f p, Vec3f wo)
{
    Vec3f L_dir{0.0f};
    Vec3f L_dir_light{0.0f};
    
    // Uniformly sample the light at x
    auto [x, light_normal, L_light] = sample_area_light(UniformSampler::next2d());
    // Shoot a ray from p to x
    Vec3f ray_dir = normalize(x - p);

    float pdf_brdf = std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id])
                         ? std::get<Lambertian>(BoxScene::materials[tri.material_id]).pdf(wo, ray_dir, tri.face_normal)
                         : std::get<Microfacet>(BoxScene::materials[tri.material_id]).pdf(wo, ray_dir, tri.face_normal);

    float dis_p_x = distance(p, x);
    float pdf_light = L_light / dot(light_normal, -ray_dir) * dis_p_x * dis_p_x;
    bool is_hit_L = RayTracer::any_hit(p, ray_dir, dis_p_x, octree, BoxScene::triangles);
    // If the ray is not blocked in the middle

    if (!is_hit_L)
    {
        Vec3f f_r;        
        if (std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id]))
        {
            f_r = std::get<Lambertian>(BoxScene::materials[tri.material_id]).eval();
        }
        if (std::holds_alternative<Microfacet>(BoxScene::materials[tri.material_id]))
        {
            f_r = std::get<Microfacet>(BoxScene::materials[tri.material_id]).eval(wo, ray_dir, tri.face_normal);
        }
        L_dir_light = eval_area_light(-ray_dir) * f_r *
                      std::abs(dot(tri.face_normal, ray_dir)) / pdf_light;
    }

    // Contribution from other reflectors
    Vec3f L_dir_brdf{0.0f};
    Vec3f L_indir{0.0f};
    const float p_rr = 0.8f;

    // Trace the new ray
    Vec3f wi{0.0f};
    float pdf = 0.0;
    if (std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id]))
    {
        std::tie(wi, pdf) = std::get<Lambertian>(BoxScene::materials[tri.material_id]).sample(tri.face_normal, UniformSampler::next2d());
    }
    else
    {
        std::tie(wi, pdf) = std::get<Microfacet>(BoxScene::materials[tri.material_id]).sample(wo, tri.face_normal, UniformSampler::next2d());
    }
    float pdf_light_prime = 0.0f;
    auto [is_hit, t, hit_tri] = RayTracer::closest_hit(p, wi, octree, BoxScene::triangles);

    // If the ray hit a non-emitting object at q
    if (is_hit)
    {
        Vec3f q_hit_point = p + t * wi;
        Vec3f f_r;
        if (std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id]))
        {
            f_r = std::get<Lambertian>(BoxScene::materials[tri.material_id]).eval();
        }
        else
        {
            f_r = std::get<Microfacet>(BoxScene::materials[tri.material_id]).eval(wo, wi, tri.face_normal);
        }
        if (is_emitter(hit_tri))
        {
            pdf_light_prime = BoxScene::inv_light_area * distance(q_hit_point, p) * distance(q_hit_point, p) /
                              std::abs(dot(hit_tri.face_normal, -wi));
            L_dir_brdf = eval_area_light(-wi) * f_r *
                         std::abs(dot(tri.face_normal, wi)) / pdf;
        }
        else if (UniformSampler::next1d() < p_rr)
        {
            L_indir = shade_with_MIS(hit_tri, q_hit_point, -wi) * f_r *
                      std::abs(dot(tri.face_normal, wi)) / pdf / p_rr;
        }
    }

    float weight_light = pdf_light / (pdf_light + pdf_brdf);
    float weight_brdf = pdf / (pdf + pdf_light_prime);

    L_dir = L_dir_light * weight_light + L_dir_brdf * weight_brdf;

    return L_dir + L_indir;
    // =============================================================================================
}

Vec3f path_tracing_with_MIS(Vec3f ray_pos, Vec3f ray_dir)
{
    const auto [is_ray_hit, t_min, nearest_tri] =
        RayTracer::closest_hit(ray_pos, ray_dir, octree, BoxScene::triangles);
    if (!is_ray_hit)
        return Vec3f{0.0f};
    const Vec3f hit_position = ray_pos + t_min * ray_dir;
    if (is_emitter(nearest_tri))
        return eval_area_light(-ray_dir);

    return shade_with_MIS(nearest_tri, hit_position, -ray_dir);
}

int main(int argc, char **argv)
{
    const unsigned int max_spps[] = {4, 32, 128, 256};

    spdlog::info("\n"
                 "----------------------------------------------\n"
                 "Welcome to CS 190I Assignment 4: Microfacet Materials\n"
                 "----------------------------------------------");
    const unsigned int max_spp = 1; // original: 32
    const unsigned int image_width = 512;
    const unsigned int image_height = 512;
    // Some prepereations
    Image image{.width = image_width,
                .height = image_height,
                .pixels = std::vector<Vec3f>(image_width * image_height)};
    Camera camera{.vertical_field_of_view = 38.6f,
                  .aspect = static_cast<float>(image_width) / image_height,
                  .focal_distance = 0.8f,
                  .position = Vec3f{0.278f, 0.8f, 0.2744f},
                  .view_direction = Vec3f{0.0f, -1.0f, 0.0f},
                  .up_direction = Vec3f{0.0f, 0.0f, 1.0f},
                  .right_direction = Vec3f{-1.0f, 0.0f, 0.0f}};
    camera.init();
    UniformSampler::init(190);

    // =============================================================================================
    // Change the material ID after you have implemented the Microfacet BRDF
    // Diffuse
    const int bunny_material_id = 6;
    // Iron
    // const int bunny_material_id = 5;
    // Gold
    // const int bunny_material_id = 6;

    // Load the scene
    // If program can't find the bunny.obj file, use xmake run -w . or move the bunny.obj file to the
    // same directory as the executable file.

    const std::string obj_path = "./bunny.obj";
    std::vector<Triangle> obj_triangles = load_obj(obj_path, bunny_material_id);
    BoxScene::triangles.insert(BoxScene::triangles.end(),
                               std::make_move_iterator(obj_triangles.begin()),
                               std::make_move_iterator(obj_triangles.end()));
    octree.build_octree(BoxScene::triangles);

    for (unsigned int max_spp : max_spps)
    {
        // =============================================================================================
        // Path Tracing with light sampling
        spdlog::info("Path Tracing with light sampling: rendering started!");
        // for (int y = 0; y < image.height; y++)
        // {
        //     if (y % 50 == 0)
        //     {
        //         spdlog::info("Rendering row {} / {} \r", y, image.height);
        //     }
        //     for (int x = 0; x < image.width; x++)
        //     {
        //         image(x, y) = Vec3f{0.0f};
        //         for (int sample = 0; sample < max_spp; sample++)
        //         {
        //             const float u = (x + UniformSampler::next1d()) / image.width;
        //             const float v = (y + UniformSampler::next1d()) / image.height;
        //             Vec3f ray_direction = camera.generate_ray(u, (1.0f - v));
        //             image(x, y) += clamp(path_tracing_with_light_sampling(
        //                                      camera.position, ray_direction),
        //                                  Vec3f(0.0f), Vec3f(50.0f));
        //         }
        //         image(x, y) /= (float)max_spp;
        //     }
        // }
        // // printf("%3f", MY/count);
        // spdlog::info("Path Tracing with light sampling: Rendering finished!");
        // // iron
        // if (bunny_material_id == 5)
        //     image.save_with_tonemapping("./path_tracing_with_light_sampling_iron_" + std::to_string(max_spp) + "spp.png");
        // // gold
        // else
        //     image.save_with_tonemapping("./path_tracing_with_light_sampling_gold_" + std::to_string(max_spp) + "spp.png");

        // =============================================================================================
        // Path Tracing with MIS
        spdlog::info("Path Tracing with MIS: Rendering started!");
        for (int y = 0; y < image.height; y++)
        {
            if (y % 50 == 0)
            {
                spdlog::info("Rendering row {} / {} \r", y, image.height);
            }
            for (int x = 0; x < image.width; x++)
            {
                image(x, y) = Vec3f{0.0f};
                for (int sample = 0; sample < max_spp; sample++)
                {
                    const float u = (x + UniformSampler::next1d()) / image.width;
                    const float v = (y + UniformSampler::next1d()) / image.height;
                    Vec3f ray_direction = camera.generate_ray(u, (1.0f - v));
                    image(x, y) +=
                        clamp(path_tracing_with_MIS(camera.position, ray_direction),
                              Vec3f(0.0f), Vec3f(50.0f));
                }
                image(x, y) /= (float)max_spp;
            }
        }
        spdlog::info("Path Tracing with MIS: Rendering finished!");
        // iron
        if (bunny_material_id == 5)
            image.save_with_tonemapping("./path_tracing_with_MIS_iron_" + std::to_string(max_spp) + "spp.png");
        // gold
        else
            image.save_with_tonemapping("./path_tracing_with_MIS_gold_" + std::to_string(max_spp) + "spp.png");
    }

    // =============================================================================================
    return 0;
}
