#pragma once
#include "../1-base/Point·三维点.hpp"
#include "../1-base/PointFP·浮点三维点.hpp"

/*
 * 三维旋转
 *
 * Overview:
 * 	Rodrigues 旋转公式:绕任意单位轴 axis 旋转 angle 弧度;
 * 	提供点和矩阵两种返回值的形式;附欧拉角 ZYX 顺序互转矩阵。
 *
 * API:
 * 	rodrigues_rotate(p, axis, angle) -> Point<ld>: 绕单位轴 axis 旋转 angle 后 p 的新坐标。浮点。O(1)。
 * 	rotation_matrix_about_axis(axis, angle) -> array<array<ld, 3>, 3>: 绕单位轴旋转的 3×3 矩阵。浮点。O(1)。
 * 	euler_zyx_to_matrix(yaw, pitch, roll) -> array<array<ld, 3>, 3>: ZYX 顺序 (先 yaw 后 pitch 后 roll) 的合成旋转矩阵。浮点。O(1)。
 * 	matrix_to_euler_zyx(m) -> tuple<ld, ld, ld> -> (yaw, pitch, roll): 反解;奇异时 pitch = ±π/2,任选 yaw 0。浮点。O(1)。
 * 	matrix_apply(m, p) -> Point<ld>: 用 3×3 矩阵对向量 p 作线性变换。浮点。O(1)。
 *
 * Notes:
 * 	axis 必须为单位向量;非单位向量应先调用 PointFP::normalize(注释强调,不在内部自动归一)。
 * 	rotate 角度约定:右手定则,从 axis 正方向看,逆时针旋转 angle > 0。
 * 	欧拉角采用 ZYX 内旋约定 (与 Tait-Bryan 顺序对应);定义见 rotation_zyx.pdf 文档。
 *
 * Related:
 * 	PointFP·浮点三维点.hpp::normalize: axis 预归一化。
 */
namespace Geo3D {

// Rodrigues 旋转公式:p' = p cos + (axis × p) sin + axis (axis · p)(1 - cos)
template<typename T>
requires std::is_floating_point_v<T>
Point<T> rodrigues_rotate(Point<T> p, Point<T> axis, T angle) {
	AST(is_zero(axis.len2() - (T)1));                            // 假设 axis 已归一:AST 确保不归一就报错
	T c = std::cos(angle);
	T s = std::sin(angle);
	T dot = axis.dot(p);
	Point<T> cross = axis.cross(p);
	return p * c + cross * s + axis * (dot * ((T)1 - c));
}

// 旋转矩阵 R = I cos + (1 - cos) axis.axisᵀ + sin [axis]×
// [axis]× 是叉积的反对称矩阵:[[0, -az, ay], [az, 0, -ax], [-ay, ax, 0]]
template<typename T>
requires std::is_floating_point_v<T>
std::array<std::array<T, 3>, 3> rotation_matrix_about_axis(Point<T> axis, T angle) {
	AST(is_zero(axis.len2() - (T)1));
	T c = std::cos(angle);
	T s = std::sin(angle);
	T omc = (T)1 - c;
	T ax = axis.x, ay = axis.y, az = axis.z;
	std::array<std::array<T, 3>, 3> R = {{
		{ c + omc * ax * ax,        omc * ax * ay - s * az,  omc * ax * az + s * ay },
		{ omc * ay * ax + s * az,  c + omc * ay * ay,        omc * ay * az - s * ax },
		{ omc * az * ax - s * ay,  omc * az * ay + s * ax,  c + omc * az * az       }
	}};
	return R;
}

// 矩阵作用向量
template<typename T>
requires std::is_floating_point_v<T>
Point<T> matrix_apply(const std::array<std::array<T, 3>, 3>& R, Point<T> p) {
	return Point<T>(
		R[0][0] * p.x + R[0][1] * p.y + R[0][2] * p.z,
		R[1][0] * p.x + R[1][1] * p.y + R[1][2] * p.z,
		R[2][0] * p.x + R[2][1] * p.y + R[2][2] * p.z
	);
}

// ZYX 顺序(yaw-pitch-roll)合成旋转:R = Rz(yaw) · Ry(pitch) · Rx(roll)
template<typename T>
requires std::is_floating_point_v<T>
std::array<std::array<T, 3>, 3> euler_zyx_to_matrix(T yaw, T pitch, T roll) {
	T cy = std::cos(yaw),   sy = std::sin(yaw);
	T cp = std::cos(pitch), sp = std::sin(pitch);
	T cr = std::cos(roll),  sr = std::sin(roll);
	std::array<std::array<T, 3>, 3> R = {{
		{ cy * cp,                cy * sp * sr - sy * cr,   cy * sp * cr + sy * sr },
		{ sy * cp,                sy * sp * sr + cy * cr,   sy * sp * cr - cy * sr },
		{ -sp,                    cp * sr,                  cp * cr               }
	}};
	return R;
}

// 反解:从矩阵得 ZYX 欧拉角;pitch 在 ±π/2 处奇,约定 pitch = π/2 时 yaw = roll = 0
template<typename T>
requires std::is_floating_point_v<T>
std::tuple<T, T, T> matrix_to_euler_zyx(const std::array<std::array<T, 3>, 3>& R) {
	T sp = -R[2][0];
	T pitch = std::asin(clamp_unit(sp));
	if (std::abs(sp - (T)1) < (T)1e-12) {
		// gimbal lock:pitch = π/2,yaw = 0 仅 roll 由矩阵前两行定
		T roll = std::atan2(-R[0][1], R[0][2]);
		return {(T)0, (T)PI / 2, roll};
	}
	if (std::abs(sp + (T)1) < (T)1e-12) {
		// pitch = -π/2,yaw = 0 仅 roll 由矩阵前两行定
		T roll = std::atan2(R[0][1], -R[0][2]);
		return {(T)0, -(T)PI / 2, roll};
	}
	T yaw  = std::atan2(R[1][0], R[0][0]);
	T roll = std::atan2(R[2][1], R[2][2]);
	return {yaw, pitch, roll};
}

} // namespace Geo3D
