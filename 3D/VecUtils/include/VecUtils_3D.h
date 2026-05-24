#pragma once
#include <cmath>

struct Vec3 
{
    double x, y, z;
};

inline Vec3 operator+(const Vec3 &v1, const Vec3 &v2) 
{
    return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
}
inline Vec3 operator-(const Vec3 &v1, const Vec3 &v2) 
{
    return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
}
inline Vec3 operator*(Vec3 const& v, double scalar)
{
    return { v.x * scalar, v.y * scalar, v.z * scalar };
}
inline Vec3 operator/(Vec3 const& v, double scalar) 
{
    return { v.x / scalar, v.y / scalar, v.z / scalar };
}

struct Cyl 
{
    double R, phi, z;
};

inline Cyl operator+(const Cyl &c1, const Cyl &c2)
{
    return { c1.R + c2.R, c1.phi + c2.phi, c1.z + c2.z };
}
inline Cyl operator-(const Cyl &c1, const Cyl &c2) 
{
    return { c1.R - c2.R, c1.phi - c2.phi, c1.z - c2.z };
}
inline Cyl operator*(const Cyl &c, double scalar)
{
    return { c.R * scalar, c.phi * scalar, c.z * scalar };
}
inline Cyl operator/(const Cyl &c, double scalar) 
{
    return { c.R / scalar, c.phi / scalar, c.z / scalar };
}

//////Conversion between coordinate functions////////////

inline Vec3 cyl_to_cart_vel(const Cyl &pos, const Cyl &vel)
{
    double cos_phi = std::cos(pos.phi);
    double sin_phi = std::sin(pos.phi);

    double vx = vel.R * cos_phi - vel.phi * sin_phi;
    double vy = vel.R * sin_phi + vel.phi * cos_phi;
    double vz = vel.z;

    return { vx, vy, vz };
}

inline Vec3 cyl_to_cart(const Cyl &cyl) 
{
    double x = cyl.R * std::cos(cyl.phi);
    double y = cyl.R * std::sin(cyl.phi);
    double z = cyl.z;

    return { x, y, z };
}

inline Cyl cart_to_cyl(const Vec3 &pos) 
{
    double R = std::sqrt(pos.x * pos.x + pos.y * pos.y);       
    double phi = std::atan2(pos.y, pos.x);                   
    if (phi < 0) phi += 2.0 * M_PI;  // φ = 0-2π
    double z = pos.z;                                        

    return { R, phi, z };
}

inline Vec3 cross_product(const Vec3& a, const Vec3& b)
{
    return 
    {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x 
    };
}

inline Cyl cart_to_cyl_vel(const Vec3& pos, const Vec3& vel)
{
    double R = std::sqrt(pos.x * pos.x + pos.y * pos.y);

    if (R < 1e-12)
    {
        return {0.0, 0.0, vel.z};
    }

    double v_R = (pos.x * vel.x + pos.y * vel.y) / R;
    double v_phi = (pos.x * vel.y - pos.y * vel.x ) / R;

    return {v_R, v_phi, vel.z}; 
}
