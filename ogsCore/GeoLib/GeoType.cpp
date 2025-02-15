/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file GeoType.cpp
 *
 * Created on 2010-12-01 by Thomas Fischer
 */

#include "GeoType.h"

namespace GeoLib {

GEOTYPE convertGeoType (const std::string& geo_type_str)
{
    if (geo_type_str.compare ("POINT") == 0) return POINT;
    if (geo_type_str.compare ("POLYLINE") == 0) return POLYLINE;
    if (geo_type_str.compare ("POLYGON") == 0) return POLYGON;
    if (geo_type_str.compare ("SURFACE") == 0) return SURFACE;
    if (geo_type_str.compare ("VOLUME") == 0) return VOLUME;
    if (geo_type_str.compare ("DOMAIN") == 0) return GEODOMAIN;
    return INVALID;
}

std::string convertGeoTypeToString (GEOTYPE geo_type)
{
    if (geo_type == POINT) return "POINT";
    if (geo_type == POLYLINE) return "POLYLINE";
    if (geo_type == POLYGON) return "POLYGON";
    if (geo_type == SURFACE) return "SURFACE";
    if (geo_type == GEODOMAIN) return "DOMAIN";
    return "INVALID";
}

} // end namespace GEOLIB
