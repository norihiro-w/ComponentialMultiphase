/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PolynomialOrder.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

namespace FemLib
{

/// Polynomial order
struct PolynomialOrder
{
    enum type {
        Linear = 1,
        Quadratic = 2,
        Cubic = 3,
        Quartic = 4,
        INVALID = -1
    };
};

} //end
