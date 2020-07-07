#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "geo.h"

#include <string>
#include <sstream>

using namespace geo;

TEST_CASE("rad degree conversion", "[base]")
{
    CHECK(deg2rad(1) == Approx(0.0174533));
    CHECK(deg2rad(400) == Approx(6.98132));

    CHECK(rad2deg(3) == Approx(171.887));
    CHECK(deg2rad(-1) == Approx(-0.0174533));
}

TEST_CASE("construct", "[vec2]")
{
    Vec2f a(2, 2);
    Vec2f c(3, 3);

    Vec2f b;
    Vec2f e(a);
    Vec2f d(a, c);

    CHECK(a.x() == Approx(2));
    CHECK(a.y() == Approx(2));

    CHECK(c.x() == Approx(3));
    CHECK(c.y() == Approx(3));

    CHECK(b.x() == Approx(0));
    CHECK(b.y() == Approx(0));

    CHECK(e.x() == Approx(2));
    CHECK(e.y() == Approx(2));

    CHECK(d.x() == Approx(1));
    CHECK(d.y() == Approx(1));

    CHECK_FALSE(a.zero());
    CHECK_FALSE(c.zero());

    CHECK(b.zero());
    CHECK_FALSE(e.zero());

    CHECK_FALSE(d.zero());

    auto y = Vec2f(1,1);
    auto z = Vec2f(1,1);

    CHECK(y == z);
    CHECK_FALSE(y == a);
    CHECK(y != a);

    std::string str_res = "[Vec2] x: 2 -- y: 2";

    CHECK(a.to_string() == str_res);

    std::stringstream ss;
    ss << (std::string)a;
    CHECK(ss.str() == str_res);
}

TEST_CASE("basic arithm", "[vec2]")
{

    Vec2f a(2, 2);
    Vec2f b(4, 4);

    // make sure the basics for the following tests are ok
    CHECK_FALSE(a.zero());
    CHECK_FALSE(b.zero());


    // tests depend on working equality operator
    Vec2f t(2, 2);
    CHECK(a.equals(t));

    SECTION("plus")
    {
        Vec2f res(6, 6);
        auto c = a + b;
        CHECK(c.equals(res));
    }

    SECTION("minus")
    {
        Vec2f res(2, 2);
        Vec2f c = b - a;
        CHECK(c.equals(res));
    }

    SECTION("scale (multiplication)")
    {
        Vec2f res(8, 8);
        Vec2f c = (b * 2);
        CHECK(c.equals(res));
    }

    SECTION("inline plus")
    {
        Vec2f res(6, 6);
        a += b;
        CHECK(a.equals(res));
    }

    SECTION("inline minus")
    {
        Vec2f res(2, 2);
        b -= a;
        CHECK(b.equals(res));
    }

    SECTION("unary minus -- invert")
    {
        Vec2f res(-2, -2);
        a = -a;
        CHECK(a.equals(res));
    }

    SECTION("minus scalar")
    {
        Vec2f res(0, 0);
        a = a - 2;
        CHECK(a.equals(res));
    }

    SECTION("divide scalar")
    {
        Vec2f res(1, 1);
        a = a / 2;
        CHECK(a.equals(res));
    }
}

TEST_CASE("functions", "[vec2]")
{

    Vec2f a(2, 2);
    Vec2f b(4, 4);

    SECTION("dot product")
    {
        CHECK(a.dot(b) == Approx(16));
    }

    SECTION("rotate clockwise around point")
    {
        auto c = a.rotate(M_PI/2, b);
        Vec2f res(2, 6);
        CHECK(c.equals(res));
    }

    SECTION("rotate counter-clockwise around origin")
    {
        auto c = a.rotate(-M_PI/2);
        Vec2f res(-2, 2);
        CHECK(c.equals(res));
    }

    SECTION("length")
    {
        CHECK(a.length() == Approx(2.82843));
    }

    SECTION("normalize")
    {
        CHECK(a.normalize().length() == Approx(1));
    }

    SECTION("angle zero")
    {
        CHECK(a.angle(b) == Approx(0));
    }

    SECTION("angle vector")
    {
        CHECK(a.angle(Vec2f(2, 4)) == Approx(0.321751));
    }

    SECTION("reflection")
    {
        auto c = -Vec2f(-1, 1).reflect(Vec2f(0, 1));
        CHECK(c.equals(Vec2f(1,1)));
    }

    SECTION("projected point")
    {
        CHECK(a.projected_point(Vec2f(2, 4)).equals(Vec2f(1.2, 2.4)));
    }

}
