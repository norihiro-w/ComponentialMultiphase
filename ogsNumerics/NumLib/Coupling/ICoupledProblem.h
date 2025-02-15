/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ICoupledProblem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "NumLib/IOSystem/IIOSystem.h"
#include "NumLib/IOSystem/INamedIO.h"
#include "NumLib/Coupling/Algorithm/IConvergenceCheck.h"

namespace NumLib
{

/**
 * \brief Interface class of coupling problems
 */
class ICoupledSystem : public IIOSystem, public INamedIO
{
public:
    ICoupledSystem() : _id(0) {};
    /// 
    virtual ~ICoupledSystem() {};

    /// solve
    virtual int solve() = 0;

    /// check consistency
    virtual bool check() const= 0;

    /// return a convergence check for this system
    virtual IConvergenceCheck* getConvergenceChecker() = 0;

    size_t getID() const {return _id;};
    void setID(size_t id) {_id=id;};
private:
    size_t _id;
};

}
