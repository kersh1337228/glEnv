#ifndef GL_ENV_LIGHTING_HPP
#define GL_ENV_LIGHTING_HPP

#include "linalg.hpp"

struct LightSource {
    float intensity;
    vec4 color;
    //
protected:
    explicit LightSource(
        const vec4& color = {1., 1., 1., 1.},
        const float intensity = 1
    ): color(color), intensity(intensity) {}
};

struct PointLight : virtual public LightSource {
    vec3 pos;
    //
    explicit PointLight(
        const vec3& pos,
        const vec4& color = {1., 1., 1., 1.},
        const float intensity = 1
    ): pos(pos), LightSource(color, intensity) {}
};

struct DirectionalLight : virtual public LightSource {
    vec3 dir;
    //
    explicit DirectionalLight(
        const vec3& dir,
        const vec4& color = {1., 1., 1., 1.},
        const float intensity = 1
    ): dir(dir.normalize()), LightSource(color, intensity) {}
};

struct SpotLight : public PointLight, public DirectionalLight {
    float inner, outer;
    //
    explicit SpotLight(
        const vec3& pos,
        const vec3& dir,
        const float inner,
        const float outer,
        const vec4& color = {1., 1., 1., 1.},
        const float intensity = 1
    ): inner(inner), outer(outer), PointLight(pos), DirectionalLight(dir), LightSource(color, intensity) {}
};

#endif
