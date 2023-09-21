#ifndef GL_TEST_GEOMETRY_HPP
#define GL_TEST_GEOMETRY_HPP

#include "../libs/glad-2.0/include/glad/glad.h"
#include "../libs/assimp-master/include/assimp/Importer.hpp"
#include "../libs/assimp-master/include/assimp/scene.h"
#include "../libs/assimp-master/include/assimp/postprocess.h"
#include "../libs/stbi/include/stbi/stb_image.h"
//
#include "../../../../Software/IDE/CLion/bin/mingw/lib/gcc/x86_64-w64-mingw32/13.1.0/include/c++/iostream"
#include "../../../../Software/IDE/CLion/bin/mingw/lib/gcc/x86_64-w64-mingw32/13.1.0/include/c++/utility"
#include "../../../../Software/IDE/CLion/bin/mingw/lib/gcc/x86_64-w64-mingw32/13.1.0/include/c++/vector"
#include "../../../../Software/IDE/CLion/bin/mingw/lib/gcc/x86_64-w64-mingw32/13.1.0/include/c++/unordered_map"
//
#include "linalg.hpp"

struct Vertex {
    vec3 pos, norm, tg, bitg;
    vec2 texCoord;
};

struct Texture {
    unsigned int id;
    aiTextureType type;
};

class Material {
private:
    std::vector<Texture> textures;
    static std::unordered_map<std::string, Texture> textures_loaded;
    static constexpr aiTextureType tts[4] = {
        aiTextureType_DIFFUSE,
        aiTextureType_SPECULAR,
        aiTextureType_HEIGHT,
        aiTextureType_AMBIENT
    };
public:
    explicit Material(const aiMaterial* const material, const std::string& dir) {
        aiString str; std::string path;
        std::unordered_map<std::string, Texture>::const_iterator used;
        Texture texture{};
        for (const aiTextureType& tt : tts)
            for(unsigned int i = 0; i < material->GetTextureCount(tt); i++) {
                material->GetTexture(tt, i, &str);
                path = std::string(str.C_Str());
                used = textures_loaded.find(path);
                if (used != textures_loaded.end()) {
                    textures.emplace_back(used->second);
                } else {
                    texture = Texture(Material::load_texture(dir + "/" + str.C_Str()), tt);
                    textures.emplace_back(texture);
                    textures_loaded.emplace(path, texture);
                }
            }
    }

    static unsigned int load_texture(const std::string& path) {
        unsigned int texture;
        glGenTextures(1, &texture);
        int width, height, channels;
        unsigned char* data = stbi_load(path.c_str(), &width, &height, &channels, 0);
        if (data) {
            const GLint format = channels == 4 ? GL_RGBA : channels == 3 ? GL_RGB : GL_RED;
            glBindTexture(GL_TEXTURE_2D, texture);
            glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
            stbi_image_free(data);
            glGenerateMipmap(GL_TEXTURE_2D);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
            glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY, 16.0f);
        } else {
            std::cerr << "Failed to load texture (" + path + ")" << std::endl;
            stbi_image_free(data);
        }
        return texture;
    }

    void bind(const Shader& shader) const noexcept {
        unsigned int diffuse = 1, specular = 1, normal = 1, height = 1;
        for (unsigned int i = 0; i < textures.size(); ++i) {
            glActiveTexture(GL_TEXTURE0 + i);
            switch (textures[i].type) {
                case aiTextureType_DIFFUSE:
                    shader.setUniform(std::string("material.").append("diffuse").append(
                        std::to_string(diffuse++)
                    ), i);
                    break;
                case aiTextureType_SPECULAR:
                    shader.setUniform(std::string("material.").append("specular").append(
                        std::to_string(specular++)
                    ), i);
                    break;
                case aiTextureType_HEIGHT:
                    shader.setUniform(std::string("material.").append("normal").append(
                        std::to_string(normal++)
                    ), i);
                    break;
                case aiTextureType_AMBIENT:
                    shader.setUniform(std::string("material.").append("height").append(
                        std::to_string(height++)
                    ), i);
                    break;
                default:
                    break;
            }
            glBindTexture(GL_TEXTURE_2D, textures[i].id);
        }
    }
};
std::unordered_map<std::string, Texture> Material::textures_loaded = {};

class Mesh {
private:
    GLsizei count;
    Material material;
    unsigned int vao = 0, vebo[2] = {0, 0};
public:
    Mesh(
        const std::vector<Vertex>& vertices,
        const std::vector<unsigned int>& indices,
        const aiMaterial* const material,
        const std::string& dir
    ): count(static_cast<GLsizei>(indices.size())), material(Material(material, dir)) {
        glGenVertexArrays(1, &vao);
        glGenBuffers(2, vebo);
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vebo[0]);
        glBufferData(
            GL_ARRAY_BUFFER,
            static_cast<GLsizeiptr>(vertices.size() * sizeof(Vertex)),
            &vertices[0], GL_STATIC_DRAW
        );
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vebo[1]);
        glBufferData(
            GL_ELEMENT_ARRAY_BUFFER,
            static_cast<GLsizeiptr>(count * sizeof(unsigned int)),
            &indices[0], GL_STATIC_DRAW
        );
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(
            0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
            reinterpret_cast<void*>(0)
        );
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(
            1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
            reinterpret_cast<void*>(offsetof(Vertex, norm))
        );
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(
            2, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
            reinterpret_cast<void*>(offsetof(Vertex, tg))
        );
        glEnableVertexAttribArray(3);
        glVertexAttribPointer(
            3, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
            reinterpret_cast<void*>(offsetof(Vertex, bitg))
        );
        glEnableVertexAttribArray(4);
        glVertexAttribPointer(
            4, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex),
            reinterpret_cast<void*>(offsetof(Vertex, texCoord))
        );
    }

    void draw(const Shader& shader) const noexcept {
        material.bind(shader);
        glBindVertexArray(vao);
        glDrawElements(GL_TRIANGLES, count, GL_UNSIGNED_INT, nullptr);
    }

    inline void clear() const noexcept {
        glDeleteVertexArrays(1, &vao);
        glDeleteBuffers(2, vebo);
    }
};

struct Transform {
    vec3 translation = vec3::null();
    vec3 scaling = vec3::identity();
    quat rotation = quat::identity();
    //
    void translate(const float dx, const float dy, const float dz) {
        translation.x() += dx;
        translation.y() += dy;
        translation.z() += dz;
    }

    void translate(const vec3& dl) {
        translation += dl;
    }

    void scale(const float sx, const float sy, const float sz) {
        scaling.x() *= sx;
        scaling.y() *= sy;
        scaling.z() *= sz;
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

    [[nodiscard]] inline mat4 model() const noexcept {
        return {
            {
                scaling.x() * (1 - 2 * (rotation.v.y() * rotation.v.y() + rotation.v.z() * rotation.v.z())),
                scaling.x() * 2 * (rotation.v.x() * rotation.v.y() - rotation.v.z() * rotation.w),
                scaling.x() * 2 * (rotation.v.x() * rotation.v.z() + rotation.v.y() * rotation.w),
                translation.x()
            }, {
                scaling.y() * 2 * (rotation.v.x() * rotation.v.y() + rotation.v.z() * rotation.w),
                scaling.y() * (1 - 2 * (rotation.v.x() * rotation.v.x() + rotation.v.z() * rotation.v.z())),
                scaling.y() * 2 * (rotation.v.y() * rotation.v.z() - rotation.v.x() * rotation.w),
                translation.y()
            }, {
                scaling.z() * 2 * (rotation.v.x() * rotation.v.z() - rotation.v.y() * rotation.w),
                scaling.z() * 2 * (rotation.v.y() * rotation.v.z() + rotation.v.x() * rotation.w),
                scaling.z() * (1 - 2 * (rotation.v.x() * rotation.v.x() + rotation.v.y() * rotation.v.y())),
                translation.z()
            }, {0.0f, 0.0f, 0.0f, 1.0f}
        };
    }
};

class Model {
private:
    std::vector<Mesh> meshes;
    Transform transform;
public:
    explicit Model(const std::string& path) {
        Assimp::Importer importer;
        const aiScene* scene = importer.ReadFile(
            path, aiProcess_Triangulate
            | aiProcess_GenSmoothNormals
            | aiProcess_CalcTangentSpace
        );
        if(!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
            throw std::ios::failure(std::string("ERROR::ASSIMP::").append(importer.GetErrorString()));
        const std::string dir = path.substr(0, path.find_last_of('/'));
        for (aiMesh **mesh = scene->mMeshes, **eom = scene->mMeshes + scene->mNumMeshes; mesh != eom; ++mesh) {
            std::vector<Vertex> vertices;
            bool normals = (*mesh)->HasNormals(), tnbt = (*mesh)->HasTangentsAndBitangents(), tc = (*mesh)->mTextureCoords[0];
            for (unsigned int i = 0; i < (*mesh)->mNumVertices; ++i)
                vertices.emplace_back(
                    vec3((*mesh)->mVertices[i].x, (*mesh)->mVertices[i].y, (*mesh)->mVertices[i].z),
                    normals ? vec3((*mesh)->mNormals[i].x, (*mesh)->mNormals[i].y, (*mesh)->mNormals[i].z) : vec3(),
                    tnbt ? vec3((*mesh)->mTangents[i].x, (*mesh)->mTangents[i].y, (*mesh)->mTangents[i].z) : vec3(),
                    tnbt ? vec3((*mesh)->mBitangents[i].x, (*mesh)->mBitangents[i].y, (*mesh)->mBitangents[i].z) : vec3(),
                    tc ? vec2((*mesh)->mTextureCoords[0][i].x, (*mesh)->mTextureCoords[0][i].y) : vec2()
                );
            std::vector<unsigned int> indices;
            for (aiFace* face = (*mesh)->mFaces, *eof = (*mesh)->mFaces + (*mesh)->mNumFaces; face != eof; ++face)
                indices.insert(indices.end(), face->mIndices, face->mIndices + face->mNumIndices);
            meshes.emplace_back(vertices, indices, scene->mMaterials[(*mesh)->mMaterialIndex], dir);
        }
    }

    void draw(const Shader& shader) const noexcept {
        shader.setUniform("model", transform.model());
        for (const Mesh& mesh : meshes)
            mesh.draw(shader);
    }

    ~Model() {
        for (const Mesh& mesh : meshes)
            mesh.clear();
    };
};

#endif
