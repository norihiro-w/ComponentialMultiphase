/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXWrapped3DVectorFunction.h
 *
 * Created on 2012-09-19 by Norihiro Watanabe
 */

#pragma once

#include "MeshLib/Core/CoordinateSystem.h"

#include "TXPosition.h"
#include "ITXFunction.h"


namespace NumLib
{

template <class F_VECTOR>
class TXWrapped3DVectorFunction : public ITXFunction
{
public:
    explicit TXWrapped3DVectorFunction(F_VECTOR* f_vec, const MeshLib::CoordinateSystem &coord_type)
    : _f_vec(f_vec), _coord_type(coord_type)
    {
        ITXFunction::isConst(f_vec->isConst());
        ITXFunction::isTemporallyConst(f_vec->isTemporallyConst());
        ITXFunction::isSpatiallyConst(f_vec->isSpatiallyConst());
    };

    virtual ~TXWrapped3DVectorFunction() {};

    void resetVectorFunction(F_VECTOR* f_vec) {_f_vec = f_vec; };

    virtual void eval(const TXPosition x, DataType &v) const OGS_DECL_OVERRIDE
    {
        DataType tmp_v;
        _f_vec->eval(x, tmp_v);
        if (tmp_v.array().size()>0) {
            v = DataType::Zero(3, 1);
            if (_coord_type.hasX())
                v(0,0) = tmp_v(_coord_type.getIndexOfX());
            if (_coord_type.hasY())
                v(1,0) = tmp_v(_coord_type.getIndexOfY());
            if (_coord_type.hasZ())
                v(2,0) = tmp_v(_coord_type.getIndexOfZ());
        }
    }

    virtual TXWrapped3DVectorFunction<F_VECTOR>* clone() const OGS_DECL_OVERRIDE
    {
        return new TXWrapped3DVectorFunction<F_VECTOR>(_f_vec, _coord_type);
    }

private:
    F_VECTOR* _f_vec;
    MeshLib::CoordinateSystem _coord_type;
};

} //end


