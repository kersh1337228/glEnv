#ifndef GL_TEST_PRIMITIVES_HPP

#define GL_TEST_PRIMITIVES_HPP

#include <glad/glad.h>
#include <stbi/stb_image.h>
#include "linalg.hpp"

namespace primitives {
    class Material {
    private:
        unsigned int diffuse, specular;
        float shininess;
    public:
        explicit Material(
            const std::string& material_path,
            float shininess = 32.0f
        ): shininess(shininess) {
            if (material_path.length()) {
                Material::load_texture(GL_TEXTURE0, diffuse, material_path + "diffuse.png");
                Material::load_texture(GL_TEXTURE1, specular, material_path + "specular.png");
            }
        }
        static void load_texture(
            GLenum glTexture,
            unsigned int& texture,
            const std::string& texture_path
        ) {
            glGenTextures(1, &texture);
            glActiveTexture(glTexture);
            glBindTexture(GL_TEXTURE_2D, texture);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
            glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY, 16.0f);
            int width, height, channels;
            unsigned char* data = stbi_load(texture_path.c_str(), &width, &height, &channels, 0);
            if (data) {
                glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
                stbi_image_free(data);
            } else
                std::cerr << "Failed to load texture (" + texture_path + ")" << std::endl;
            glGenerateMipmap(GL_TEXTURE_2D);
        }
        void bind(const Shader& shader) const {
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, diffuse);
            glActiveTexture(GL_TEXTURE1);
            glBindTexture(GL_TEXTURE_2D, specular);
            shader.setUniform("material.diffuse", 0);
            shader.setUniform("material.specular", 1);
            shader.setUniform("material.shininess", shininess);
        }
    };

    class Primitive {
    protected:
        vec3 translation;
        vec3 scaling;
        quat rotation;
        Material material;
    public:
        explicit Primitive(
            const vec3& translation = {0., 0., 0.},
            const vec3& scaling = {1., 1., 1.},
            const quat& rotation = {1., 0., 0., 0.},
            const std::string& material_path = ""
        ): translation(translation), scaling(scaling), rotation(rotation), material(Material(material_path)) {};
        void translate(const float dx, const float dy, const float dz) {
            translation.x += dx;
            translation.y += dy;
            translation.z += dz;
        }
        void translate(const vec3& dl) {
            translation += dl;
        }
        void scale(const float sx, const float sy, const float sz) {
            scaling.x *= sx;
            scaling.y *= sy;
            scaling.z *= sz;
        }
        void scale(const vec3& ds) {
            scaling *= ds;
        }
        void rotate(const float angle, const vec3& axis) {
            rotation = rotation.quatmul(quat::rotate(angle, axis));
        }
        void rotate(const float pitch, const float yaw, const float roll) {
            rotation = rotation.quatmul(
                quat::rotate(pitch, {1., 0., 0.}).quatmul(
                    quat::rotate(yaw, {0., 1., 0.}).quatmul(
                        quat::rotate(roll, {0., 0., 1.})
                    )
                )
            );
        }
        [[nodiscard]] mat4 model() const {
            return mat4::translate(translation).matmul(
                mat4::scale(scaling).matmul(
                    rotation.to_mat4()
                )
            );
        }
        [[nodiscard]] virtual const float* get_vertices() const {
            return {};
        };
        [[nodiscard]] virtual const unsigned int* get_indices() const {
            return {};
        };
        virtual void draw(const Shader& shader) const {
            material.bind(shader);
        };
    };

    class Cube : public Primitive {
    public:
        static const float vertices[288];
        static const unsigned int indices[36];
        explicit Cube(
            const vec3& translate = {0., 0., 0.},
            const vec3& scale = {1., 1., 1.},
            const quat& rotation = {1., 0., 0., 0.},
            const std::string& material_path = ""
        ): Primitive(translate, scale, rotation, material_path) {}
        [[nodiscard]] const float* get_vertices() const override {
            return Cube::vertices;
        }
        [[nodiscard]] const unsigned int* get_indices() const override {
            return Cube::indices;
        }
        void draw(const Shader& shader) const override {
            this->Primitive::draw(shader);
            glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, nullptr);
        }
    };

    const float Cube::vertices[288] = {
        -0.5f, -0.5f, -0.5f, 0.0f,  0.0f, -1.0f,  0.0f,  0.0f,
        0.5f, -0.5f, -0.5f, 0.0f,  0.0f, -1.0f,  1.0f,  0.0f,
        0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  1.0f,  1.0f,
        0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  1.0f,  1.0f,
        -0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  0.0f,  1.0f,
        -0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  0.0f,  0.0f,

        -0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  0.0f,  0.0f,
        0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  1.0f,  0.0f,
        0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  1.0f,  1.0f,
        0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  1.0f,  1.0f,
        -0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  0.0f,  1.0f,
        -0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  0.0f,  0.0f,

        -0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f,  1.0f,  0.0f,
        -0.5f,  0.5f, -0.5f, -1.0f,  0.0f,  0.0f,  1.0f,  1.0f,
        -0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,  0.0f,  1.0f,
        -0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,  0.0f,  1.0f,
        -0.5f, -0.5f,  0.5f, -1.0f,  0.0f,  0.0f,  0.0f,  0.0f,
        -0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f,  1.0f,  0.0f,

        0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f,  1.0f,  0.0f,
        0.5f,  0.5f, -0.5f,  1.0f,  0.0f,  0.0f,  1.0f,  1.0f,
        0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,  0.0f,  1.0f,
        0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,  0.0f,  1.0f,
        0.5f, -0.5f,  0.5f,  1.0f,  0.0f,  0.0f,  0.0f,  0.0f,
        0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f,  1.0f,  0.0f,

        -0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,  0.0f,  1.0f,
        0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,  1.0f,  1.0f,
        0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,  1.0f,  0.0f,
        0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,  1.0f,  0.0f,
        -0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,  0.0f,  0.0f,
        -0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,  0.0f,  1.0f,

        -0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,  0.0f,  1.0f,
        0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,  1.0f,  1.0f,
        0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,  1.0f,  0.0f,
        0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,  1.0f,  0.0f,
        -0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,  0.0f,  0.0f,
        -0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,  0.0f,  1.0f
    };

    const unsigned int Cube::indices[36] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
        24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35
    };
}

#endif
