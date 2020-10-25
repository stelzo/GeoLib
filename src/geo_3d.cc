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

Shape3::Shape3(const std::vector<std::vector<Vec3f>> &layers) : vertices(layers) {}

/*
bool Shape3::contains(const Vec3f &p)
{
    size_t nvert = vertices.size();
    size_t i, j;
    i = j = 0;
    bool c = false;

    for (i = 0, j = nvert - 1; i < nvert; j = i++)
    {
        if (((vertices[i].y() > p.y()) != (vertices[j].y() > p.y())) &&
            (p.x() < (vertices[j].x() - vertices[i].x()) * (p.y() - vertices[i].y()) / (vertices[j].y() - vertices[i].y()) + vertices[i].x()))
        {
            c = !c;
        }
    }

    return c;
}
*/

/*
void Shape3::smooth(float alpha, float tension, unsigned int distance)
{
    size_t size = vertices.size();
    std::vector<Vec2f> smoothed;

    for (unsigned int i = 0; i < size - 3; i++)
    {
        float t01 = pow(Vec2f(vertices[i], vertices[i + 1]).length(), alpha);
        float t12 = pow(Vec2f(vertices[i + 1], vertices[i + 2]).length(), alpha);
        float t23 = pow(Vec2f(vertices[i + 2], vertices[i + 3]).length(), alpha);

        float x1, y1, x2, y2;

        x1 = (1.0 - tension) * (vertices[i + 2].x() - vertices[i + 1].x() + t12 * ((vertices[i + 1].x() - vertices[i].x()) / t01 - (vertices[i + 2].x() - vertices[i].x()) / (t01 + t12)));
        y1 = (1.0 - tension) * (vertices[i + 2].y() - vertices[i + 1].y() + t12 * ((vertices[i + 1].y() - vertices[i].y()) / t01 - (vertices[i + 2].y() - vertices[i].y()) / (t01 + t12)));

        x2 = (1.0 - tension) * (vertices[i + 2].x() - vertices[i + 1].x() + t12 * ((vertices[i + 3].x() - vertices[i + 2].x()) / t23 - (vertices[i + 3].x() - vertices[i + 1].x()) / (t12 + t23)));
        y2 = (1.0 - tension) * (vertices[i + 2].y() - vertices[i + 1].y() + t12 * ((vertices[i + 3].y() - vertices[i + 2].y()) / t23 - (vertices[i + 3].y() - vertices[i + 1].y()) / (t12 + t23)));

        float _a_x, _a_y, _b_x, _b_y, _c_x, _c_y, _d_x, _d_y;
        _a_x = 2.0 * (vertices[i + 1].x() - vertices[i + 2].x()) + x1 + x2;
        _a_y = 2.0 * (vertices[i + 1].y() - vertices[i + 2].y()) + y1 + y2;
        _b_x = -3.0 * (vertices[i + 1].x() - vertices[i + 2].x()) - x1 - x1 - x2;
        _b_y = -3.0 * (vertices[i + 1].y() - vertices[i + 2].y()) - y1 - y1 - y2;
        _c_x = x1;
        _c_y = y1;
        _d_x = vertices[i + 1].x();
        _d_y = vertices[i + 1].y();

        for (int j = 0; j < 100; j = j + distance)
        {
            float t = (float)j / 100;
            float x_s = _a_x * t * t * t + _b_x * t * t + _c_x * t + _d_x;
            float y_s = _a_y * t * t * t + _b_y * t * t + _c_y * t + _d_y;

            smoothed.emplace_back(x_s, y_s);
        }
    }

    vertices = smoothed;
}
*/

void Shape3::move(const Vec3f& movement)
{
    for (auto& layer : vertices) {
        for (auto& vertex : layer) {
            vertex += movement;
        }
    }
}

void Shape3::rotate(float rad, const Vec3f& pivot)
{
    for (auto& layer : vertices) {
        for (auto &vertex : layer) {
            vertex = vertex.rotate(rad, pivot);
        }
    }
}

Shape3::Shape3(const Shape3 &v) : vertices(v.vertices) {}

Shape3 Shape3::scaled (float dist) {
    std::vector<std::vector<Vec3f>> sized;

    for (int j = 0; j < (int) vertices.size(); ++j) {
        Vec3f con(vertices[j].back(), vertices[j][1]);
        Vec3f rn(con.y(), -con.x(), con.z());
        sized[j].push_back(rn.normalize() * dist + vertices[j][0]);

        for (int i = 1; i < (int) vertices[j].size() - 1; ++i) {
            Vec3f connection(vertices[j][i - 1], vertices[j][i + 1]);
            Vec3f rnorm(connection.y(), -connection.x(), connection.z());
            sized[j].push_back(rnorm.normalize() * dist + vertices[j][i]);
        }


        Vec3f c(vertices[j][vertices[j].size() - 2], vertices[j][0]);
        Vec3f r(c.y(), -c.x(), c.z());
        sized[j].push_back(r.normalize() * dist + vertices[j][vertices[j].size() - 1]);
    }
    return Shape3(sized);
}

std::string Shape3::to_string () const {
    std::stringstream ss;
    ss << "[Shape3]\n";
    for (auto layer : vertices) {
        ss << "[Shape3 new layer]\n";
        for (auto vertex : layer) {
            ss << vertex.to_string() << "\n";
        }
    }

    return ss.str();
}

/*
double Shape3::volume () {
    double phalf = 0, nhalf = 0;

    for (int i = 0; i < (int) vertices.size(); ++i) {
        phalf += vertices[i].x() * vertices[(i + 1) % vertices.size()].y();
        nhalf += vertices[(i + 1) % vertices.size()].x() * vertices[i].y();
    }

    return std::abs(phalf - nhalf) / 2;
}
*/