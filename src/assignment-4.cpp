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

bool same_triangle(const Triangle &a, const Triangle &b)
{
    return a.v0 == b.v0 && a.v1 == b.v1 && a.v2 == b.v2;
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

std::tuple<Vec3f, float> sample_whole_sphere()
{
    float u1 = UniformSampler::next1d();
    float u2 = UniformSampler::next1d();
    float theta = acos(1 - 2 * u1);
    float phi = 2 * M_PI * u2;
    Vec3f wi = Vec3f{sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
    return std::make_tuple(wi, 1.0f / (4 * M_PI));
}
Vec3f shade(Triangle tri, Vec3f p, Vec3f wo, bool debug)
{
    const float p_rr = 0.8f;
    float ksi = UniformSampler::next1d();
    if (ksi > p_rr)
        return Vec3f{0.0f};
    Vec3f wi;
    float pdf;
    if (std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id]))
        std::tie(wi, pdf) = std::get<Lambertian>(BoxScene::materials[tri.material_id]).sample(tri.face_normal, UniformSampler::next2d());
    else
        std::tie(wi, pdf) = std::get<Dielectric_BSDF>(BoxScene::materials[tri.material_id]).sample(wo, tri.face_normal, UniformSampler::next2d());

    // Avoid self-intersection
    if (dot(wi, tri.face_normal) < 0.0f)
        p = offset_ray_origin(p, -1.0f * tri.face_normal);
    else
        p = offset_ray_origin(p, tri.face_normal);

    // Trace the new ray
    auto [hit, t_min, nearest_tri] =
        RayTracer::closest_hit(p, wi, octree, BoxScene::triangles);
    if (hit && same_triangle(nearest_tri, tri))
    {
        throw std::runtime_error("Self-intersection!");
        return Vec3f{0.0f};
    }

    Vec3f L_i = eval_area_light(-wi);
    float cosine = std::abs(dot(wi, tri.face_normal));
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
                return shade(nearest_tri, q, -wi, debug) * f_r * cosine / pdf / p_rr;
            }
        }
    }
    else
    {
        // Debug
        bool param = debug;
        Vec3f f_r = std::get<Dielectric_BSDF>(BoxScene::materials[tri.material_id]).eval(wo, wi, tri.face_normal, param);
        if (hit)
        {
            if (is_emitter(nearest_tri))
            {
                return L_i * f_r * cosine / pdf / p_rr;
            }
            else
            {
                Vec3f q = p + t_min * wi;
                // Debug
                Vec3f tmp = shade(nearest_tri, q, -wi, debug);
                // Judge green
                // if (debug && tmp.y > tmp.x && tmp.y > tmp.z)
                // {
                //     printf("tmp: %f, %f, %f\n", tmp.x, tmp.y, tmp.z);
                //     printf("f_r: %f, %f, %f\n", f_r.x, f_r.y, f_r.z);
                //     printf("cosine: %f\n", cosine);
                //     printf("pdf: %f\n", pdf);
                //     printf("p_rr: %f\n", p_rr);
                // }
                return tmp * f_r * cosine / pdf / p_rr;
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

Vec3f path_tracing(Vec3f ray_pos, Vec3f ray_dir, bool debug = false)
{
    const auto [is_ray_hit, t_min, nearest_tri] =
        RayTracer::closest_hit(ray_pos, ray_dir, octree, BoxScene::triangles);
    if (!is_ray_hit)
        return Vec3f{0.0f};
    const Vec3f hit_position = ray_pos + t_min * ray_dir;
    if (is_emitter(nearest_tri))
        return eval_area_light(-ray_dir);
    return shade(nearest_tri, hit_position, -ray_dir, debug);
}
void test_bsdf_sample(Camera);
void test_compute_half_vector();
void test_geometry();
void test_sphere(Camera);
int main(int argc, char **argv)
{
    const unsigned int max_spps[] = {32};
    const unsigned int image_width = 256;
    const unsigned int image_height = 256;
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

    const std::string obj_path = "./sphere_high.obj";
    Vec3f bias, scale;
    if (obj_path == "./sphere_high.obj")
    {
        bias = Vec3f{0.275f, -0.35f, 0.25f};
        scale = Vec3f{0.2f};
    }
    else if (obj_path == "./bunny.obj")
    {
        bias = Vec3f{0.0f, -0.0f, 0.0f};
        scale = Vec3f{1.0f};
    }
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
    bool debug = 0;
    if (debug)
    {
        // test_compute_half_vector();
        // test_bsdf_sample(camera);
        // test_geometry();
        test_sphere(camera);
        // Shoot a ray to the center of the image and print the hit triangle position
        // const float u = 0.5f;
        // const float v = 0.5f;
        // Vec3f ray_direction = camera.generate_ray(u, (1.0f - v));
        // printf("ray_direction: %f, %f, %f\n", ray_direction.x, ray_direction.y, ray_direction.z);
        // const auto [is_ray_hit, t_min, nearest_tri] =
        //     RayTracer::closest_hit(camera.position, ray_direction, octree, BoxScene::triangles);
        // if (is_ray_hit)
        // {
        //     Vec3f hit_position = camera.position + t_min * ray_direction;
        //     spdlog::info("The hit triangle position: ({}, {}, {})", hit_position.x, hit_position.y, hit_position.z);
        //     // shoot another ray based on the hit position
        //     printf("1hit_position: %f, %f, %f\n", hit_position.x, hit_position.y, hit_position.z);
        //     // hit_position = offset_ray_origin(hit_position, nearest_tri.face_normal);
        //     printf("2hit_position: %f, %f, %f\n", hit_position.x, hit_position.y, hit_position.z);
        //     Vec3f ray_direction = hit_position - camera.position;
        //     printf("ray_direction: %f, %f, %f\n", ray_direction.x, ray_direction.y, ray_direction.z);
        //     printf("1: tri_normal: %f, %f, %f\n", nearest_tri.face_normal.x, nearest_tri.face_normal.y, nearest_tri.face_normal.z);
        //     const auto [is_ray_hit, t_min, nearest_tri] =
        //         RayTracer::closest_hit(hit_position, ray_direction, octree, BoxScene::triangles);
        //     if (is_ray_hit)
        //     {
        //         printf("2: tri_normal: %f, %f, %f\n", nearest_tri.face_normal.x, nearest_tri.face_normal.y, nearest_tri.face_normal.z);

        //         const Vec3f hit_position1 = hit_position + t_min * ray_direction;
        //         spdlog::info("The hit triangle position: ({}, {}, {})", hit_position1.x, hit_position1.y, hit_position1.z);
        //     }
        //     else
        //     {
        //         spdlog::info("No triangle is hit!");
        //     }
        // }
        // else
        // {
        //     spdlog::info("No triangle is hit!");
        // }

        // Dielectric_BSDF dielectric_bsdf{
        //     .roughness = 0.0f,
        //     .n1 = 1.0f,
        //     .n2 = 1.5f};
        // // Test the Fresnel function
        // Vec3f w_perp = Vec3f{0, 0, 1};
        // Vec3f w_parl = Vec3f{1, 0, 0};
        // float cos_theta = cos(8.0f / 9.0f * M_PI / 2.0f);
        // Vec3f w_common = normalize(Vec3f{sqrt(1 - cos_theta * cos_theta), 0, cos_theta});
        // float test_F_perp = dielectric_bsdf.F(w_perp);
        // float test_F_parl = dielectric_bsdf.F(w_parl);
        // float test_F_neg_perp = dielectric_bsdf.F(-w_perp);
        // float test_F_neg_parl = dielectric_bsdf.F(-w_parl);
        // float test_F_common = dielectric_bsdf.F(w_common);
        // float test_F_neg_common = dielectric_bsdf.F(-w_common, debug = true);
        // spdlog::info("Fresnel-perp: {}", test_F_perp);
        // spdlog::info("Fresnel-parl: {}", test_F_parl);
        // spdlog::info("Fresnel-neg_perp: {}", test_F_neg_perp);
        // spdlog::info("Fresnel-neg_parl: {}", test_F_neg_parl);
        // printf("w_common: %f, %f, %f\n", w_common.x, w_common.y, w_common.z);
        // spdlog::info("Fresnel-common: {}", test_F_common);
        // spdlog::info("Fresnel-neg_common: {}", test_F_neg_common);
    }
    if (!debug)
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
                        bool debug = false;
                        if (x >= 160 && x <= 192 && y >= 128 && y <= 160)
                        {
                            debug = true;
                        }
                        const float u = (x + UniformSampler::next1d()) / image.width;
                        const float v = (y + UniformSampler::next1d()) / image.height;
                        Vec3f ray_direction = camera.generate_ray(u, (1.0f - v));
                        image(x, y) += clamp(path_tracing(
                                                 camera.position, ray_direction, debug),
                                             Vec3f(0.0f), Vec3f(50.0f));
                    }
                    image(x, y) /= (float)max_spp;
                }
            }
            image(180, 144) = Vec3f(10000.0f, 0.0f, 0.0f);
            spdlog::info("Path Tracing with light sampling: Rendering finished!");
            image.save_with_tonemapping("./path_tracing" + std::to_string(max_spp) + ".png");
            return 0;
        }
}

void test_bsdf_sample(Camera camera)
{
    spdlog::info("\nStart testing the sample function of the Dielectric_BSDF class!");
    // Shoot a ray to the center of the image, get the hit triangle (dielectric) and test its sample function, count how many transmitted rays are there
    const float u = 0.5f;
    const float v = 0.5f;
    Vec3f ray_direction = camera.generate_ray(u, (1.0f - v));
    const auto [is_ray_hit, t_min, nearest_tri] =
        RayTracer::closest_hit(camera.position, ray_direction, octree, BoxScene::triangles);
    if (!(is_ray_hit && std::holds_alternative<Dielectric_BSDF>(BoxScene::materials[nearest_tri.material_id])))
    {
        throw std::runtime_error("The hit triangle is not a dielectric!");
    }
    const Vec3f hit_position = camera.position + t_min * ray_direction;
    const Vec3f wo = -ray_direction;
    const Vec3f normal = nearest_tri.face_normal;
    const Vec3f p = offset_ray_origin(hit_position, normal);
    int transmitted_rays = 0;
    for (int i = 0; i < 1000; i++)
    {
        Vec3f wi;
        float pdf;
        std::tie(wi, pdf) = std::get<Dielectric_BSDF>(BoxScene::materials[nearest_tri.material_id]).sample(wo, normal, UniformSampler::next2d());
        // if wi and wo are not in the same hemisphere, it's a transmitted ray
        if (dot(wi, normal) * dot(wo, normal) < 0.0f)
        {
            transmitted_rays++;
        }
    }
    spdlog::info("First shoot: the number of transmitted rays: {}", transmitted_rays);

    // Shoot a ray from the hit position to the front, get the hit triangle (dielectric) and test its sample function, count how many transmitted rays are there
    Vec3f ray_direction1 = hit_position - camera.position;
    const auto [is_ray_hit1, t_min1, nearest_tri1] =
        RayTracer::closest_hit(hit_position, ray_direction1, octree, BoxScene::triangles);
    if (!(is_ray_hit1 && std::holds_alternative<Dielectric_BSDF>(BoxScene::materials[nearest_tri1.material_id])))
    {
        throw std::runtime_error("Second hit: the hit triangle is not a dielectric!");
    }
    const Vec3f hit_position1 = hit_position + t_min1 * ray_direction1;
    const Vec3f wo1 = -ray_direction1;
    const Vec3f normal1 = nearest_tri1.face_normal;
    const Vec3f p1 = offset_ray_origin(hit_position1, normal1);
    int transmitted_rays1 = 0;
    for (int i = 0; i < 1000; i++)
    {
        Vec3f wi;
        float pdf;
        std::tie(wi, pdf) = std::get<Dielectric_BSDF>(BoxScene::materials[nearest_tri1.material_id]).sample(wo1, normal1, UniformSampler::next2d());
        // if wi and wo are not in the same hemisphere, it's a transmitted ray
        if (dot(wi, normal1) * dot(wo1, normal1) < 0.0f)
        {
            transmitted_rays1++;
        }
    }
    spdlog::info("Second shoot: the number of transmitted rays: {}", transmitted_rays1);
    // The second ray should hit another side of the sphere, the material id should be 7
    if (nearest_tri1.material_id != 7)
    {
        // Print the material of the back wall of the box, print the position of the triangle
        printf("material_id: %d\n", nearest_tri1.material_id);
        printf("position: %f, %f, %f\n", hit_position1.x, hit_position1.y, hit_position1.z);
        throw std::runtime_error("The second hit: the hit triangle is not the sphere!");
    }
    spdlog::info("The second hit: the hit triangle is the sphere!");
    // Shoot a ray from the second hit position to the back, get the hit triangle (back wall) and test its sample function, count how many transmitted rays are there

    Vec3f sphere_front = hit_position;
    Vec3f sphere_back = hit_position1;
    Triangle tri_front = nearest_tri;
    Triangle tri_back = nearest_tri1;

    Vec3f ray_direction2 = sphere_back - sphere_front;
    const auto [is_ray_hit2, t_min2, nearest_tri2] =
        RayTracer::closest_hit(hit_position1, ray_direction2, octree, BoxScene::triangles);
    if (!(is_ray_hit2 && nearest_tri2.material_id == 0))
    {
        throw std::runtime_error("The third hit: the hit triangle is not the back wall!");
    }
    spdlog::info("The third hit: the hit triangle is the back wall!");

    Dielectric_BSDF glass = std::get<Dielectric_BSDF>(BoxScene::materials[tri_front.material_id]);
    // Test the fr value at the front position
    Vec3f front_direction = normalize(sphere_back - sphere_front);
    Vec3f back_direction = normalize(sphere_front - sphere_back);
    float fs = glass.eval(back_direction, front_direction, tri_front.face_normal, true).x;
    if (!fs)
    {
        spdlog::info("Get small fs value at the front position: {}", fs);
        glass.eval(back_direction, front_direction, tri_front.face_normal, true);
    }
    printf("fs: %f\n", fs);
}
void test_compute_half_vector()
{
    spdlog::info("\nStart testing the computeHalfVector function!");
    float critical_angle = asin(BoxScene::Glass.n1 / BoxScene::Glass.n2) - 0.05;
    float perp_angle = 0;
    float grazing_angle = 0.5f * M_PI;
    float common_angle = (critical_angle + perp_angle) / 2;

    Vec3f perp = Vec3f{0, 0, 1};
    Vec3f grazing = Vec3f{1, 0, 0};
    Vec3f common = Vec3f{sqrt(1 - cos(common_angle) * cos(common_angle)), 0, cos(common_angle)};

    Vec3f wo, wi;
    // All test cases expect to see the half vector is 0, 0, 1
    // Test case: wo at top, wi at bottom
    wo = perp;
    wi = -perp;
    Vec3f half_vector = BoxScene::Glass.computeHalfVector(wo, wi);
    if (half_vector != Vec3f{0, 0, 1})
    {
        printf("half_vector: %f, %f, %f\n", half_vector.x, half_vector.y, half_vector.z);
        throw std::runtime_error("Test case perp failed!");
    }
    spdlog::info("Perp test passed!");
    // Test case: wo at right, wi at left
    wo = grazing;
    wi = -grazing;
    half_vector = BoxScene::Glass.computeHalfVector(wo, wi, true);
    if (half_vector != Vec3f{0, 0, 1})
    {
        printf("half_vector: %f, %f, %f\n", half_vector.x, half_vector.y, half_vector.z);
        throw std::runtime_error("Test case graz failed!");
    }
    spdlog::info("Grazing test passed!");
    // Test case: wo at common, wi at common
    wo = common;
    wi = -common;
    wi.z = wo.z;
    half_vector = BoxScene::Glass.computeHalfVector(wo, wi);
    if (half_vector != Vec3f{0, 0, 1})
    {
        printf("half_vector: %f, %f, %f\n", half_vector.x, half_vector.y, half_vector.z);
        throw std::runtime_error("Test case common failed!");
    }
    spdlog::info("Common test passed!");
    // Test case: refraction
    // simulate the refraction, where wo is near the critical angle, wi is near {0, 0, -1}
    wo = Vec3f{std::sqrt(1 - cos(critical_angle) * cos(critical_angle)), 0, cos(critical_angle)};
    wi = normalize(Vec3f{0, 0.05, -0.85});
    half_vector = BoxScene::Glass.computeHalfVector(wo, wi);
    printf("half_vector: %f, %f, %f\n", half_vector.x, half_vector.y, half_vector.z);
    spdlog::info("Refraction test passed!");
}

void test_geometry()

{
    spdlog::info("\nStart testing the geometry function!");
    Dielectric_BSDF glass = BoxScene::Glass;
    Vec3f w = Vec3f{0, 0, 1};
    Vec3f x = Vec3f{1, 0, 0};
    int n = 8;
    // for (int i = 0; i <= n; i++)
    // {
    //     Vec3f wx = normalize(x * i + w * (n - i));
    //     printf("wx: %f, %f, %f\n", wx.x, wx.y, wx.z);
    //     printf("F(wx): %f\n", glass.F(wx));
    //     printf("D(wx): %f\n", glass.D(wx));
    // }
}
void test_sphere(Camera camera)
{
    // Shoot a ray to (180, 144) and test its eval function
    spdlog::info("\nStart testing the eval function of the Dielectric_BSDF class!");
    int x = 180, y = 144;
    const float u = (x + UniformSampler::next1d()) / 256;
    const float v = (y + UniformSampler::next1d()) / 256;
    Vec3f ray_direction = camera.generate_ray(u, (1.0f - v));
    const auto [is_ray_hit, t_min, nearest_tri] =
        RayTracer::closest_hit(camera.position, ray_direction, octree, BoxScene::triangles);
    if (!(is_ray_hit && std::holds_alternative<Dielectric_BSDF>(BoxScene::materials[nearest_tri.material_id])))
    {
        throw std::runtime_error("The hit triangle is not a dielectric!");
    }
    const Vec3f hit_position = camera.position + t_min * ray_direction;
    const Vec3f wo = -ray_direction;
    Dielectric_BSDF glass = std::get<Dielectric_BSDF>(BoxScene::materials[nearest_tri.material_id]);
    // Sample the wi and cound how many reflected rays are there
    int reflected_rays = 0;
    for (int i = 0; i < 1000; i++)
    {
        Vec3f wi;
        float pdf;
        std::tie(wi, pdf) = glass.sample(wo, nearest_tri.face_normal, UniformSampler::next2d());
        // if wi and wo are in the same hemisphere, it's a reflected ray
        if (dot(wi, nearest_tri.face_normal) * dot(wo, nearest_tri.face_normal) > 0.0f)
        {
            reflected_rays++;
        }
    }
    spdlog::info("The number of reflected rays: {}", reflected_rays);
}