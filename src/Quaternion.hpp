//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Motocyte
//
//        Copyright (C) 2019 Satya N.V. Arjunan
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// Motocyte is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// Motocyte is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with Motocyte -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Satya N. V. Arjunan <satya.arjunan@gmail.com>
// based on the Java cell motility simulation code by Mark N. Read
//


#ifndef __Quaternion_hpp
#define __Quaternion_hpp

#include <Vector.hpp>

class Quaternion {
public:
  Quaternion(const double w, const Vector<double>& v):
    w_(w),
    v_(v) {}

  ~Quaternion() {};

  static Quaternion represent_rotation(const double angle,
                                       Vector<double> v) {
    double len2(v.x*v.x + v.y*v.y + v.z*v.z);
    if (len2 != 1) {
      v = v/sqrt(len2);
    }
    double a(angle/2);
    return Quaternion(cos(a), v*sin(a));
  }

  Quaternion inverse() {
    double len2(w_*w_ + v_.x*v_.x + v_.y*v_.y + v_.z*v_.z);
    return Quaternion(w_/len2, v_/(-len2));
  }

  Quaternion normalise() {
    double len2(w_*w_ + v_.x*v_.x + v_.y*v_.y + v_.z*v_.z);
    w_ = w_/len2;
    v_ = v_/len2;
    return *this;
  }

  Quaternion multiply(const Quaternion& q) {
    const double new_w(w_*q.w_ - v_.x*q.v_.x - v_.y*q.v_.y - v_.z*q.v_.z);
    const Vector<double> new_v(
                         w_*q.v_.x + v_.x*q.w_ + v_.y*q.v_.z - v_.z*q.v_.y,
                         w_*q.v_.y - v_.x*q.v_.z + v_.y*q.w_ + v_.z*q.v_.x,
                         w_*q.v_.z + v_.x*q.v_.y - v_.y*q.v_.x + v_.z*q.w_);
    return Quaternion(new_w, new_v);
  }

  Quaternion passive_rotate(const double angle, const Vector<double>& vec) {
    // y_ = q(yq') passive rotation
    Quaternion q(represent_rotation(angle, vec));
    Quaternion q_(q.inverse());
    Quaternion yq_(multiply(q_));
    Quaternion y_(q.multiply(yq_).normalise());
    return y_;
  }

	static Quaternion face_vector(const Vector<double> dir) {
    const Vector<double> xaxis(1, 0, 0);  
    // the normal is perpendicular to the plane on which bounce and facing
    // vectors lie. 
    Vector<double> normal(xaxis.cross_product(dir));
    // happens when dir is parallel to xaxis. 
    if (normal.magnitude_squared() == 0) {
      // select y axis as default, but any vector perpendicular to x will do.
      normal = Vector<double>(0, 1, 0);
    }
    const double angle(xaxis.angle_between_vectors(dir)); // in radians.
    return Quaternion::represent_rotation(angle, normal);
  }

  const Vector<double> transform(const Vector<double>& vec) {
    const double x(w_*w_*vec.x + 2*v_.y*w_*vec.z - 2*v_.z*w_*vec.y +
                   v_.x*v_.x*vec.x + 2*v_.y*v_.x*vec.y + 2*v_.z*v_.x*vec.z -
                   v_.z*v_.z*vec.x - v_.y*v_.y*vec.x);
    const double y(2*v_.x*v_.y*vec.x +   v_.y*v_.y*vec.y + 2*v_.z*v_.y*vec.z +
                   2*w_*v_.z*vec.x - v_.z*v_.z*vec.y +   w_*w_*vec.y -
                   2*v_.x*w_*vec.z - v_.x*v_.x*vec.y);
    const double z(2*v_.x*v_.z*vec.x + 2*v_.y*v_.z*vec.y +   v_.z*v_.z*vec.z -
                   2*w_*v_.y*vec.x - v_.y*v_.y*vec.z + 2*w_*v_.x*vec.y -
                   v_.x*v_.x*vec.z + w_*w_*vec.z);
    return Vector<double>(x, y, z);
	}

  double w_;
  Vector<double> v_;
};

#endif /* __Quaternion_hpp */


