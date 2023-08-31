#include <glad/glad.h>
#include <GLFW/glfw3.h>
//
#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
//
#include "glEnv/shader.hpp"
#include "glEnv/camera.hpp"
#include "glEnv/geometry.hpp"
#include "glEnv/lighting.hpp"

extern "C" {
    __declspec(dllexport) unsigned long NvOptimusEnablement = 1;
}

// Constants
static const unsigned int width = 800.;
static const unsigned int height = 600.;
static Camera camera(
    {0., 0., 3.},
    {0., 0., -1.},
    {0., 1., 0.}
);
static double dt = 0, t = 0;
static double xp = width / 2., yp = height / 2.;

void keyboard_handler(GLFWwindow* window) {
    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
    if(glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.translate(0., 0., -1.);
    if(glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.translate(-1., 0., 0.);
    if(glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.translate(0., 0., 1.);
    if(glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.translate(1., 0., 0.);
}

void mouse_handler(GLFWwindow* window, double x, double y) {
    camera.rotate(yp - y, xp - x);
    xp = x; yp = y;
}

int main() {
    // GLFW Init
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 4);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    // Window
    GLFWwindow* window = glfwCreateWindow(width, height, "gl_test", nullptr, nullptr);
    if (window == nullptr) {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    // GLAD Init
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }
    // Viewport automatic resize
    glfwSetFramebufferSizeCallback(window, [](GLFWwindow* window, int width, int height) -> void {
        glViewport(0, 0, width, height);
    });
    glfwSetCursorPosCallback(window, mouse_handler);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glEnable(GL_STENCIL_TEST);
    glStencilFunc(GL_NOTEQUAL, 1, 0xFF);
    glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
    // Scene
    Shader shader("../assets/shaders/phong.vert", "../assets/shaders/phong.frag");
    Shader outline("../assets/shaders/outline.vert", "../assets/shaders/outline.frag");
    Model m("../assets/models/backpack/backpack.obj");
    // Setting up shader program
    glClearColor(0.01568627450980392, 0.10196078431372549, 0.25098039215686274, 1.);
    DirectionalLight dl({0.2, 1., 0.3}, {1., 1., 1., 1.});
    PointLight pl({0., 0., 0.}, {1., 1., 0.8, 1.});
    SpotLight sl (
        {1.5, 2.5, 1.}, {1.5, 3.5, 0.2},
        std::cos(radians(12.5)), std::cos(radians(14.)),
        {1., 1., 0.5, 1.}
    );
    // Runtime loop
    const mat4 model = mat4::identity();
    const mat4 projection = mat4::perspective(
        camera.FoV,
        static_cast<float>(width) / static_cast<float>(height),
        0.1, 100.
    );
    while (!glfwWindowShouldClose(window)) {
        keyboard_handler(window);
        // Frame time
        const double t0 = glfwGetTime();
        dt = t0 - t; t = t0;
        camera.speed = 7 * dt;
        camera.sensitivity = 0.5 * dt;
        pl.pos = {static_cast<float>(3 * cos(t0)), 0.3, static_cast<float>(3 * sin(t0))};
        //
        shader.use();
        shader.setUniform("pl.color", pl.color);
        shader.setUniform("dl.color", dl.color);
        shader.setUniform("sl.inner", sl.inner);
        shader.setUniform("sl.outer", sl.outer);
        shader.setUniform("sl.color", sl.color);
        shader.setUniform("projection", projection);
        // Matrices
        const mat4 view = camera.view();
        shader.setUniform("view", view);
        shader.setUniform("norm", (view.matmul(model)).transpose().inverse());
        // Lighting
        shader.setUniform("pl.pos", vec3(view.matmul(vec4(pl.pos, 1.))));
        shader.setUniform("dl.dir", vec3(view.matmul(vec4(dl.dir, 0.))));
        shader.setUniform("sl.pos", vec3(view.matmul(vec4(sl.pos, 1.))));
        shader.setUniform("sl.dir", vec3(view.matmul(vec4(sl.dir, 0.))));
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
        // Draw call
        glStencilFunc(GL_ALWAYS, 1, 0xFF);
        glStencilMask(0xFF);
        m.draw(shader);
        // Outline
        glStencilFunc(GL_EQUAL, 0, 0xFF);
        glStencilMask(0x00);
        glDisable(GL_DEPTH_TEST);
        outline.use();
        outline.setUniform("view", view);
        outline.setUniform("projection", projection);
        outline.setUniform("outlineWidth", 0.004f);
        outline.setUniform("outlineColor", {1., 1., 0., 1.});
        m.draw(outline);
        glStencilMask(0xFF);
        glStencilFunc(GL_ALWAYS, 0, 0xFF);
        glEnable(GL_DEPTH_TEST);
        // Swap
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glfwTerminate();
}
