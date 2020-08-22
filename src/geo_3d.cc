// GeoLib - 2D/3D Geometry
// Copyright (C) 2020 Christopher Sieh (stelzo@steado.de)

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include "geo.h"
#include <math.h>

#include <sstream>

using namespace geo;

Vec3f::Vec3f(const Vec3f &from, const Vec3f &to) : _x(to._x - from._x), _y(to._y - from._y), _z(to._z - from._z) {}

std::string Vec3f::to_string() const
{
    std::stringstream ss;
    ss << "[Vec3] "
       << "x: " << _x << " -- "
       << "y: " << _y << " -- "
       << "z: " << _z;

    return ss.str();
}

Vec3f::~Vec3f() {}

Vec3f::operator std::string() const
{
    return to_string();
}

Vec3f::Vec3f(float x, float y, float z) : _x(x), _y(y), _z(z) {}

Vec3f::Vec3f() : _x(0), _y(0), _z(0) {}

Vec3f::Vec3f(const Vec3f &v) : _x(v._x), _y(v._y), _z(v._z) {}