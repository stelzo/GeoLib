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

    template <class _T>
    class Vec2
    {

    private:
        _T _x;
        _T _y;

    public:
        Vec2(_T x, _T y);
        Vec2();
        Vec2(const Vec2<_T> &v);
        Vec2(const Vec2<_T> &from, const Vec2<_T> &to);

        _T x() const;

        _T y() const;

        Vec2<_T> operator+(const Vec2<_T> &o) const;

        Vec2<_T> &operator+=(const Vec2<_T> &o);

        Vec2<_T> &operator=(const Vec2<_T> &o);

        Vec2<_T> &operator-=(const Vec2<_T> &o);

        Vec2<_T> operator-(const Vec2<_T> &o) const;

        Vec2<_T> operator-();

        Vec2<_T> operator-(_T s) const;

        Vec2<_T> operator/(_T s) const;

        double dot(const Vec2<_T> &v) const;

        bool empty() const;

        double lengthSquared() const;

        Vec2<_T> projected_point(const Vec2<_T> &v);

        Vec2<_T> closest_vec2_to(const Vec2<_T> &v);

        Vec2<_T> rotate(float rad, const Vec2<_T> &pivot = Vec2<_T>(0, 0));

        Vec2<_T> reflect(const Vec2<_T> &v, const Vec2<_T> &n);

        double length() const;

        Vec2<_T> normalize() const;

        Vec2<_T> operator*(_T s) const;

        double angle(const Vec2<_T> &v, bool sign_by_x_coord = true);

        bool equals(const Vec2<_T> &v) const;

        bool operator==(const Vec2<_T> &b) const;

        bool operator!=(const Vec2<_T> &b) const;

        operator std::string() const;

        std::string to_string() const;
    };

} // namespace geo

// template -> need to include the cc directly after header for linking reasons
#include "geo_2d.cc"

#endif /* geo_2d_hpp */