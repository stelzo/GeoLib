//
//  geo_2d.hpp
//  libgeo
//
//  Created by Christopher Sieh on 21.06.20.
//  Copyright © 2020 steado. All rights reserved.
//

#ifndef geo_2d_hpp
#define geo_2d_hpp

#include <geo_base.h>
#include <string>

namespace geo
{

    class Vec2f
    {

    private:
        float _x;
        float _y;

    public:
        Vec2f(float x, float y);
        Vec2f();
        Vec2f(const Vec2f &v);
        Vec2f(const Vec2f &from, const Vec2f &to);

        float x() const;

        float y() const;

        Vec2f operator+(const Vec2f &o) const;

        Vec2f &operator+=(const Vec2f &o);

        Vec2f &operator=(const Vec2f &o);

        Vec2f &operator-=(const Vec2f &o);

        Vec2f operator-(const Vec2f &o) const;

        Vec2f operator-();

        Vec2f operator-(float s) const;

        Vec2f operator/(float s) const;

        double dot(const Vec2f &v) const;

        bool empty() const;

        double lengthSquared() const;

        Vec2f projected_point(const Vec2f &v);

        Vec2f closest_vec2_to(const Vec2f &v);

        Vec2f rotate(float rad, const Vec2f &pivot = Vec2f(0, 0));

        Vec2f reflect(const Vec2f &v, const Vec2f &n);

        double length() const;

        Vec2f normalize() const;

        Vec2f operator*(float s) const;

        double angle(const Vec2f &v, bool sign_by_x_coord = true);

        bool equals(const Vec2f &v) const;

        bool operator==(const Vec2f &b) const;

        bool operator!=(const Vec2f &b) const;

        operator std::string() const;

        std::string to_string() const;
    };

} // namespace geo

#endif /* geo_2d_hpp */