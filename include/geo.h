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

#ifndef geo_h
#define geo_h

#include <string>

namespace geo
{
    // Converts degrees to radians.
    //
    // @param degrees
    // @return radians
    double deg2rad(double deg);

    // Converts radians to degrees.
    //
    // @param radians
    // @return degrees
    double rad2deg(double rad);

    // Vec2f represents a 2D vector in space which also can be a point.
    // It implements basic arithmetic and helper functions, to create and manipulate
    // points and vectors.
    //
    // For position information precision float was chosen because of its small memory footprint
    // and therefore better caching in algorithms and use cases where milimeter precision
    // is enough when 1 meter equals 1.0f.
    //
    // This class should never have any kind of setter methods. Algorithms and conversions need
    // to be implemented inside the library.
    class Vec2f
    {

    private:
        float _x;
        float _y;

        // Squared length for utilization in multiple algorithms
        //
        // @return length squared
        double _lengthSquared() const;

    public:
        // Constructor with initialization of both positions.
        //
        // @param x: position in x dimension
        // @param y: position in y dimension
        Vec2f(float x, float y);

        // Constructor initializes with x = y = 0.
        //
        Vec2f();

        // Copy constructor initializes with x and y from v.
        //
        // @param vector to include from
        Vec2f(const Vec2f &v);

        // Constructs a vector from a point to another point.
        //
        // @param from: origin
        // @param to: destination
        Vec2f(const Vec2f &from, const Vec2f &to);

        // Getter for x position.
        //
        // @return x position
        float x() const;

        // Getter for y position.
        //
        // @return y position
        float y() const;

        // Addition of both x and y positions.
        //
        // @param other vector
        // @return new vector with the values from the addition
        Vec2f operator+(const Vec2f &o) const;

        // Addition with another vector but the calling class gets the result.
        //
        // @param other vector
        // @return *this for concatenation
        Vec2f &operator+=(const Vec2f &o);

        // Assigns the calling class to the values of the given vector.
        //
        // @param other vector
        // @return *this for concatenation
        Vec2f &operator=(const Vec2f &o);

        // Substraction with another vector but the calling class gets the result.
        //
        // @param other vector
        // @return *this for concatenation
        Vec2f &operator-=(const Vec2f &o);

        // Substraction with another vector.
        //
        // @param other vector
        // @return new vector with the result of the substraction.
        Vec2f operator-(const Vec2f &o) const;

        // Unary inversion. All x and y change the sign.
        //
        // @return new vector with inverted direction
        Vec2f operator-() const;

        // Substraction with a scalar. All elements get substracted by s.
        //
        // @param scalar
        Vec2f operator-(float s) const;

        // Division by a scalar. All elements get divided by s.
        //
        // @param scalar
        Vec2f operator/(float s) const;

        // Dot product with another vector. Mind the order. a.dot(b) == a . b.
        //
        // @param vector
        double dot(const Vec2f &v) const;

        // Checks whether a vector has values near zero (std::limits::epsilon).
        //
        // @return true if values are near zero, else false
        bool zero() const;

        // Calculates the projected point from calling vector on vector v.
        // This is also the closest point possible from the destination the
        // calling vector is pointing to to the line that v describes.
        //
        // @param vector to project on
        // @return point on v
        Vec2f projected_point(const Vec2f &v) const;

        // Calculates the direct vector to any point on v.
        // The resulting vector is has the smallest length possible.
        //
        // @param vector that describes the line
        // @return vector from destination of calling vector to v
        Vec2f closest_vec2_to(const Vec2f &v) const;

        // Rotates a vector around a pivot (origin is default).
        // Positive radians mean clockwise,
        // negative radians mean anticlockwise rotation.
        //
        // @param radians to rotate
        // @param pivot to rotate around
        // @return rotated vector
        Vec2f rotate(float rad, const Vec2f &pivot = Vec2f(0, 0)) const;

        // Reflects a vector on another vector and returns the resulting bounce vector
        // originated from v in the normal direction of v.
        //
        // @param vector v
        // @param normal n
        // @return vector reflected on v
        Vec2f reflect(const Vec2f &v, const Vec2f &n) const;

        // Length of the vector in units.
        //
        // @return length
        double length() const;

        // Normalizes vector of the calling vector and returns it.
        // Does not change the calling vector.
        //
        // @return normalized vector, length 1
        Vec2f normalize() const;

        // Scales the vector with a scalar s.
        //
        // @param scalar
        // @return scaled vector
        Vec2f operator*(float s) const;

        // Calculates the angle between the calling vector and v.
        // Positive results mean clockwise,
        // negative results mean anticlockwise rotation.
        //
        // Signing mechanism written by Patrick Hoffmann.
        //
        // @param v: vector/point to find angle to.
        // @param sign_coord: whether to check the x or y position for the
        // resulting sign.
        // @return angle in radians
        double angle_signed(const Vec2f &v, bool sign_coord_x = true) const;

        // Calculates the angle between the calling vector and v.
        //
        // @param vector to find angle to.
        // @return angle in radians
        double angle(const Vec2f &v) const;

        // Check equality of two vectors. Same functionality as operator==
        // but mainly used for testing.
        //
        // @param vector v to check equality against.
        // @return true if same values, else false
        bool equals(const Vec2f &v) const;

        // Check equality of two vectors.
        //
        // @param vector v to check equality against.
        // @return true if same values, else false
        bool operator==(const Vec2f &b) const;

        // Check inequality of two vectors.
        //
        // @param vector v to check equality against.
        // @return true if different values, else false
        bool operator!=(const Vec2f &b) const;

        // Implicit string conversion. Minimal string representation
        // of a vector. Same functionality as to_string().
        //
        // @return string representation of the vector.
        operator std::string() const;

        // Minimal string representation of a vector.
        //
        // @return string representation of the vector.
        std::string to_string() const;
    };

} // namespace geo

#endif // geo_h