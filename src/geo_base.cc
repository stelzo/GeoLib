//
//  geo_base.cpp
//  libgeo
//
//  Created by Christopher Sieh on 21.06.20.
//  Copyright © 2020 steado. All rights reserved.
//

#include "geo_base.h"
#include <math.h>

namespace geo
{
    double deg2rad(float deg)
    {
        return deg * M_PI / 180;
    }

    double rad2deg(float rad)
    {
        return rad * 180 / M_PI;
    }

} // namespace geo
