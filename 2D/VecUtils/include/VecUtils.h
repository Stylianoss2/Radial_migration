#ifndef VecUtils_h
#define VecUtils_h

#include <cmath>

struct Vec2 
{
    double x, y;
};

struct Cyl 
{
    double R, phi;
};

inline Vec2 operator+(const Vec2 &v1, const Vec2 &v2) 
{
    return {v1.x + v2.x, v1.y + v2.y};
}
inline Vec2 operator-(const Vec2 &v1, const Vec2 &v2) 
{
    return {v1.x - v2.x, v1.y - v2.y};
}
inline Vec2 operator*(Vec2 const& v, double scalar)
{
    return {v.x * scalar, v.y * scalar};
}
inline Vec2 operator/(Vec2 const& v, double scalar) 
{
    return {v.x / scalar, v.y / scalar}; 
}

inline Cyl operator+(const Cyl &c1, const Cyl &c2)
{
    return {c1.R + c2.R, c1.phi + c2.phi};
}
inline Cyl operator-(const Cyl &c1, const Cyl &c2) 
{
    return {c1.R - c2.R, c1.phi - c2.phi};
}
inline Cyl operator*(const Cyl &c, double scalar)
{
    return {c.R * scalar, c.phi * scalar};
}
inline Cyl operator/(const Cyl &c, double scalar) 
{
    return { c.R / scalar, c.phi / scalar};
}

//////Conveersion between coordinate functions////////////
inline Vec2 cyl_to_cart(const Cyl &cyl) 
{
    double x = cyl.R * std::cos(cyl.phi);
    double y = cyl.R * std::sin(cyl.phi);
    return {x, y};
}

inline Cyl cart_to_cyl(const Vec2 &pos) 
{
    double R = std::hypot(pos.x, pos.y);  
    double phi = std::atan2(pos.y, pos.x);

    if (phi < 0.0) 
        phi += 2.0 * M_PI;

    return {R, phi};
}

inline Cyl cart_to_cyl_vel(const Vec2& pos, const Vec2& vel)
{
    const double R = std::hypot(pos.x, pos.y);
    if (R < 1e-12) return {0.0, 0.0};

    const double vR   = (pos.x * vel.x + pos.y * vel.y) / R;
    const double vphi = (pos.x * vel.y - pos.y * vel.x) / R;
    return { vR, vphi };
}

inline Vec2 cart_to_cyl_inertial(const Vec2& pos, const Vec2& vel)
{
    double R = std::hypot(pos.x, pos.y);
    if (R < 1e-12) return {0.0, 0.0};

    double c = pos.x / R;
    double s = pos.y / R;

    double vR   =  c * vel.x + s * vel.y;
    double vphi = -s * vel.x + c * vel.y;

    return {vR, vphi};
}

inline Vec2 cyl_to_cart_vel(double vR, double vphi, const Vec2& pos)
{
    double phi = std::atan2(pos.y, pos.x);
    double vx = vR * cos(phi) - vphi * sin(phi);
    double vy = vR * sin(phi) + vphi * cos(phi);

    return { vx, vy }; 
}

#endif