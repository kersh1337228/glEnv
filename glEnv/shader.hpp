#ifndef GL_ENV_SHADER_HPP
#define GL_ENV_SHADER_HPP

#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include "linalg.hpp"

class Shader {
public:
    Shader(const std::string& vPath, const std::string& fPath) {
        this->vShader = Shader::compile(GL_VERTEX_SHADER, vPath);
        this->fShader = Shader::compile(GL_FRAGMENT_SHADER, fPath);
        this->program = Shader::link(this->vShader, this->fShader);
    }
    void setShader(GLenum type, const std::string& sPath) {
        const unsigned int shader = Shader::compile(type, sPath);
        switch (type) {
            case GL_VERTEX_SHADER:
                this->vShader = shader;
                break;
            case GL_GEOMETRY_SHADER:
                break;
            case GL_FRAGMENT_SHADER:
                this->fShader = shader;
                break;
            default:
                throw std::invalid_argument("Shader_ type passed does not exists");
        }
        this->program = Shader::link(this->vShader, this->fShader);
    }
    void use() const {
        glUseProgram(this->program);
    }
    ~Shader() {
        glDeleteProgram(this->program);
    }
    [[maybe_unused]] void setUniform(const std::string& name, const vec4& v) const {
        glUniform4f(
            glGetUniformLocation(this->program, name.c_str()),
            v.x(), v.y(), v.z(), v.w()
        );
    }
    [[maybe_unused]] void setUniform(const std::string& name, const vec3& v) const {
        glUniform3f(
            glGetUniformLocation(this->program, name.c_str()),
            v.x(), v.y(), v.z()
        );
    }
    [[maybe_unused]] void setUniform(const std::string& name, const float v) const {
        glUniform1f(glGetUniformLocation(this->program, name.c_str()), v);
    }
    [[maybe_unused]] void setUniform(const std::string& name, const int v) const {
        glUniform1i(glGetUniformLocation(this->program, name.c_str()), v);
    }
    [[maybe_unused]] void setUniform(const std::string& name, const unsigned int v) const {
        glUniform1ui(glGetUniformLocation(this->program, name.c_str()), v);
    }
    [[maybe_unused]] void setUniform(const std::string& name, const mat4& v) const {
        glUniformMatrix4fv(
            glGetUniformLocation(this->program, name.c_str()),
            1, GL_TRUE, v.raw()
        );
    }
private:
    unsigned int vShader, fShader, program;
    static std::map<GLenum, std::string> typeNames;

    static unsigned int compile(GLenum type, const std::string& path) {
        unsigned int shader = glCreateShader(type);
        std::stringstream ss;
        std::ifstream file(path, std::ios::in);
        if (file.is_open()) {
            ss << file.rdbuf();
            std::string raw = ss.str();
            const char* code = raw.c_str();
            glShaderSource(shader, 1, &code, nullptr);
            file.close();
        } else
            throw std::ios::failure("ERROR::SHADER::" + Shader::typeNames[type] + "::OPENING_FAILED");
        glCompileShader(shader);
        int success;
        glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
        if(!success) {
            int bufSize = 512;
            char infoLog[bufSize];
            glGetShaderInfoLog(shader, bufSize, nullptr, infoLog);
            std::cerr << "ERROR::SHADER::" << Shader::typeNames[type] << "::COMPILATION_FAILED\n" << infoLog << std::endl;
        }
        return shader;
    }
    static unsigned int link(unsigned int vShader, unsigned int fShader) {
        const unsigned int program = glCreateProgram();
        glAttachShader(program, vShader);
        glAttachShader(program, fShader);
        glLinkProgram(program);
        int success;
        glGetProgramiv(program, GL_LINK_STATUS, &success);
        if(!success) {
            int bufSize = 512;
            char infoLog[bufSize];
            glGetProgramInfoLog(program, bufSize, nullptr, infoLog);
            std::cerr << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
        }
        glDeleteShader(vShader);
        glDeleteShader(fShader);
        return program;
    }
};

std::map<GLenum, std::string> Shader::typeNames = {
    {GL_VERTEX_SHADER, "VERTEX"},
    {GL_GEOMETRY_SHADER, "GEOMETRY"},
    {GL_FRAGMENT_SHADER, "FRAGMENT"},
};

#endif
