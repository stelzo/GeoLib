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

Vec2f::Vec2f(const Vec2f &from, const Vec2f &to) : _x(to._x - from._x), _y(to._y - from._y) {}

std::string Vec2f::to_string() const
{
    std::stringstream ss;
    ss << "[Vec2] "
       << "x: " << _x << " -- "
       << "y: " << _y;

    return ss.str();
}

Vec2f::~Vec2f() {}

Vec2f::operator std::string() const
{
    return to_string();
}

Vec2f::Vec2f(float x, float y) : _x(x), _y(y) {}

Vec2f::Vec2f() : _x(0), _y(0) {}

Vec2f::Vec2f(const Vec2f &v) : _x(v._x), _y(v._y) {}

Polygon2::Polygon2(const std::vector<std::vector<Vec2f>> &sides)
{
    for (auto side : sides)
    {
        vertices.insert(vertices.end(), side.begin(), side.end());
    }
}

Polygon2::Polygon2(const std::vector<Vec2f> &vertices) : vertices(vertices) {}

bool Polygon2::contains(const Vec2f &p)
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

void Polygon2::smooth(float alpha, float tension, unsigned int distance)
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

Polygon2::Polygon2(const Polygon2 &v) : vertices(v.vertices) {}