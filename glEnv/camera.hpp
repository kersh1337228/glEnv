#ifndef GL_TEST_CAMERA_HPP

#define GL_TEST_CAMERA_HPP

#include "linalg.hpp"

class Camera {
protected:
    vec3 pos, x, y, z, up;
public:
    float FoV, speed, sensitivity;
    //
    Camera(
        const vec3& pos, const vec3& dir, const vec3& up,
        const float FoV = PI / 4., const float speed = 0.05f, const float sensitivity = 0.001f
    ): pos(pos), z(-dir.normalize()), up(up), FoV(FoV), speed(speed), sensitivity(sensitivity) {
        x = up.cross(z).normalize();
        y = z.cross(x);
    }
    void rotate(float pitch = 0, float yaw = 0) {
        if (std::abs(z.dot(up)) > cos(PI / 180.))
            pitch = 0;
        const quat rot = quat::rotate(pitch * sensitivity, x).quatmul(
                quat::rotate(yaw * sensitivity, y)
        );
        z = rot.quatmul({0., z.x(), z.y(), z.z()}).quatmul(rot.conj()).v;
        x = up.cross(z).normalize();
        y = z.cross(x);
    }
    void rotate(float pitch, float yaw, float roll) {
        const quat rot = quat::rotate(pitch * sensitivity, x).quatmul(
            quat::rotate(yaw * sensitivity, y).quatmul(
                quat::rotate(roll * sensitivity, z)
            )
        );
        x = rot.quatmul({0., x.x(), x.y(), x.z()}).quatmul(rot.conj()).v;
        y = rot.quatmul({0., y.x(), y.y(), y.z()}).quatmul(rot.conj()).v;
        z = rot.quatmul({0., z.x(), z.y(), z.z()}).quatmul(rot.conj()).v;
    }
    void rotate(float angle, const vec3& axis) {
        const quat rot = quat::rotate(angle, axis);
        x = rot.quatmul({0., x.x(), x.y(), x.z()}).quatmul(rot.conj()).v;
        y = rot.quatmul({0., y.x(), y.y(), y.z()}).quatmul(rot.conj()).v;
        z = rot.quatmul({0., z.x(), z.y(), z.z()}).quatmul(rot.conj()).v;
    }
    virtual void translate(float dx, float dy, float dz) {
        pos += (x * dx + y * dy + z * dz) * speed;
    }
    [[nodiscard]] mat4 view() const {
        return {
            {x.x(), x.y(), x.z(), -pos.x() * x.x() - pos.y() * x.y() - pos.z() * x.z()},
            {y.x(), y.y(), y.z(), -pos.x() * y.x() - pos.y() * y.y() - pos.z() * y.z()},
            {z.x(), z.y(), z.z(), -pos.x() * z.x() - pos.y() * z.y() - pos.z() * z.z()},
            {0.0f, 0.0f, 0.0f, 1.0f}
        };
    }
};

class FPCamera : public Camera {
public:
    FPCamera(
            const vec3& pos, const vec3& dir, const vec3& up,
            float speed = 0.05f, float sensitivity = 0.001f
    ): Camera(pos, dir, up, speed, sensitivity) {}
    void translate(float dx, float dy, float dz) override {
        pos += (x * dx + y * dy + z * dz) * speed;
        pos.y() = 0;
    }
};

#endif


