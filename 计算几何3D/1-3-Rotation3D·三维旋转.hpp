#pragma once
#include "1-1-Geo3D·三维浮点基础.hpp"

/*
 * 三维旋转与欧拉角变换
 *
 * Overview:
 *     提供三维空间中绕任意单位轴旋转的 Rodrigues 公式（点变换与旋转矩阵），以及 ZYX
 *     内旋顺序（yaw-pitch-roll）欧拉角与 3x3 旋转矩阵的互转。
 *
 * API:
 *     rodrigues_rotate(p, axis, angle)        — 绕单位轴 axis 旋转 angle 弧度，
 *                                                O(1)。
 *     rotation_matrix_about_axis(axis, angle) — 绕单位轴旋转的 3x3 旋转矩阵，O(1)。
 *     matrix_apply(R, p)                      — 3x3 矩阵作用于三维点/向量，O(1)。
 *     euler_zyx_to_matrix(yaw, pitch, roll)   — ZYX 欧拉角转 3x3 旋转矩阵，O(1)。
 *     matrix_to_euler_zyx(R)                  — 3x3 旋转矩阵反解 ZYX 欧拉角，O(1)。
 *
 * Notes:
 *     1. axis 必须为单位向量（内部含 AST 检查，非单位向量请先调用 axis.unit()）。
 *     2. 旋转方向遵循右手定则：拇指指向 axis 正方向，四指弯曲方向为 angle > 0。
 *     3. ZYX 欧拉角对应 R = Rz(yaw) * Ry(pitch) * Rx(roll)；万向节死锁时约定 yaw =
 *        0。
 *
 * Related:
 *     1-1-Geo3D·三维浮点基础.hpp — 提供三维点/向量及 safe_asin。
 */

namespace Geo3D {

// Rodrigues 公式：点 p 绕单位向量 axis 按右手定则旋转 angle 弧度
template<typename T = ld>
Point<T> rodrigues_rotate(Point<T> p, Point<T> axis, T angle) {
    AST(is_zero(axis.len2() - 1));
    T c = std::cos(angle), s = std::sin(angle);
    return p * c + axis.cross(p) * s + axis * (axis.dot(p) * (1 - c));
}

// 构造绕单位轴 axis 旋转 angle 弧度的 3x3 正交旋转矩阵
template<typename T = ld>
std::array<std::array<T, 3>, 3> rotation_matrix_about_axis(Point<T> axis, T angle) {
    AST(is_zero(axis.len2() - 1));
    T c = std::cos(angle), s = std::sin(angle), omc = 1 - c;
    T x = axis.x, y = axis.y, z = axis.z;
    return {{
        {c + omc * x * x, omc * x * y - s * z, omc * x * z + s * y},
        {omc * y * x + s * z, c + omc * y * y, omc * y * z - s * x},
        {omc * z * x - s * y, omc * z * y + s * x, c + omc * z * z}
    }};
}

template<typename T = ld>
Point<T> matrix_apply(const std::array<std::array<T, 3>, 3>& R, Point<T> p) {
    return {
        R[0][0] * p.x + R[0][1] * p.y + R[0][2] * p.z,
        R[1][0] * p.x + R[1][1] * p.y + R[1][2] * p.z,
        R[2][0] * p.x + R[2][1] * p.y + R[2][2] * p.z
    };
}

// ZYX 内旋欧拉角转旋转矩阵 R = Rz(yaw) * Ry(pitch) * Rx(roll)
template<typename T = ld>
std::array<std::array<T, 3>, 3> euler_zyx_to_matrix(T yaw, T pitch, T roll) {
    T cy = std::cos(yaw), sy = std::sin(yaw);
    T cp = std::cos(pitch), sp = std::sin(pitch);
    T cr = std::cos(roll), sr = std::sin(roll);
    return {{
        {cy * cp, cy * sp * sr - sy * cr, cy * sp * cr + sy * sr},
        {sy * cp, sy * sp * sr + cy * cr, sy * sp * cr - cy * sr},
        {-sp, cp * sr, cp * cr}
    }};
}

// 3x3 旋转矩阵反解 ZYX 欧拉角 (yaw, pitch, roll)，pitch = +-pi/2 万向节死锁时置 yaw
// = 0
template<typename T = ld>
std::tuple<T, T, T> matrix_to_euler_zyx(const std::array<std::array<T, 3>, 3>& R) {
    T sp = -R[2][0];
    T pitch = safe_asin(sp);
    if (std::abs(sp - 1) < (T)1e-12) {
        return {0, (T)PI / 2, std::atan2(R[0][1], R[0][2])};
    }
    if (std::abs(sp + 1) < (T)1e-12) {
        return {0, -(T)PI / 2, std::atan2(-R[0][1], -R[0][2])};
    }
    return {std::atan2(R[1][0], R[0][0]), pitch, std::atan2(R[2][1], R[2][2])};
}

} // namespace Geo3D
