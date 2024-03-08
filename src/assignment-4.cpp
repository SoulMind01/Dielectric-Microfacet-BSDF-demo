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

Vec3f shade(Triangle tri, Vec3f p, Vec3f wo)
{
    const float p_rr = 0.8f;
    float ksi = UniformSampler::next1d();
    if (ksi > p_rr)
        return Vec3f{0.0f};

    auto [wi, pdf] =
        sample_uniform_hemisphere(tri.face_normal, UniformSampler::next2d());

    // Avoid self-intersection
    p = offset_ray_origin(p, tri.face_normal);

    // Trace the new ray
    auto [hit, t_min, nearest_tri] =
        RayTracer::closest_hit(p, wi, octree, BoxScene::triangles);

    Vec3f L_i = eval_area_light(-wi);
    float cosine = dot(wi, tri.face_normal);
    if (cosine < 0.0f)
        cosine = 0.0f;
    bool flag = std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id]);
    if (flag)
    {
        Vec3f f_r = std::get<Lambertian>(BoxScene::materials[tri.material_id]).eval();
        if (hit)
        {
            if (is_emitter(nearest_tri))
            {
                return L_i * f_r * cosine / pdf / p_rr;
            }
            else
            {
                Vec3f q = p + t_min * wi;
                return shade(nearest_tri, q, -wi) * f_r * cosine / pdf / p_rr;
            }
        }
    }
    else
    {
        Vec3f f_r = std::get<Dielectric_BSDF>(BoxScene::materials[tri.material_id]).eval(wo, wi, tri.face_normal);
        if (hit)
        {
            if (is_emitter(nearest_tri))
            {
                return L_i * f_r * cosine / pdf / p_rr;
            }
            else
            {
                Vec3f q = p + t_min * wi;
                return shade(nearest_tri, q, -wi) * f_r * cosine / pdf / p_rr;
            }
        }
    }
    return Vec3f{0.0f};
    printf("material_id: %d\n", tri.material_id);
    printf("holds_alternative: %d\n", std::holds_alternative<Microfacet>(BoxScene::materials[tri.material_id]));
    printf("holds_alternative: %d\n", std::holds_alternative<Dielectric_BSDF>(BoxScene::materials[tri.material_id]));
    printf("holds_alternative: %d\n", std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id]));
    printf("flag: %d\n", flag);
    throw std::runtime_error("Unknown material type");
}

Vec3f path_tracing(Vec3f ray_pos, Vec3f ray_dir)
{
    const auto [is_ray_hit, t_min, nearest_tri] =
        RayTracer::closest_hit(ray_pos, ray_dir, octree, BoxScene::triangles);
    if (!is_ray_hit)
        return Vec3f{0.0f};
    const Vec3f hit_position = ray_pos + t_min * ray_dir;
    if (is_emitter(nearest_tri))
        return eval_area_light(-ray_dir);
    return shade(nearest_tri, hit_position, -ray_dir);
}

int main(int argc, char **argv)
{
    const unsigned int max_spps[] = {64};
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
    const int bunny_material_id = 7;
    // Iron: 5, Gold: 6, Glass: 7

    const std::string obj_path = "./sphere.obj";
    const Vec3f bias = Vec3f{0.275f, -0.35f, 0.25f};
    const float scale = 0.2f;
    std::vector<Triangle> obj_triangles = load_obj(obj_path, bunny_material_id);
    for (auto &i : obj_triangles)
    {
        i.v0 = i.v0 * scale + bias;
        i.v1 = i.v1 * scale + bias;
        i.v2 = i.v2 * scale + bias;
    }
    BoxScene::triangles.insert(BoxScene::triangles.end(),
                               std::make_move_iterator(obj_triangles.begin()),
                               std::make_move_iterator(obj_triangles.end()));
    octree.build_octree(BoxScene::triangles);

    for (unsigned int max_spp : max_spps)
    {
        // Path Tracing with light sampling
        spdlog::info("Path Tracing with light sampling: rendering started!");
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
                    image(x, y) += clamp(path_tracing(
                                             camera.position, ray_direction),
                                         Vec3f(0.0f), Vec3f(50.0f));
                }
                image(x, y) /= (float)max_spp;
            }
        }
        spdlog::info("Path Tracing with light sampling: Rendering finished!");
        image.save_with_tonemapping("./path_tracing" + std::to_string(max_spp) + ".png");
        return 0;
    }
}
