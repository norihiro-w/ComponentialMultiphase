/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file THMCSimulator.h
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include "GeoProcessBuilder.h"

// forward declaration
namespace NumLib
{
class ITransientCoupledSystem;
}

namespace ogs6
{

// forward declaration
class SimulationInfo;

/**
 * \brief initialize OGS
 * @param argc
 * @param argv
 */
void ogsInit(int argc, char* argv[]);

/**
 * \brief finalize OGS
 */
void ogsExit();

/**
 * \brief OGS simulator class
 */
class THMCSimulator
{
public:
    typedef GeoProcessBuilder ProcessBuilder;

    /**
     *
     * @param argc
     * @param argv
     */
    THMCSimulator(int argc, char* argv[]);

    ///
    ~THMCSimulator();

    /**
     *
     * @return error code
     */
    int execute();

private:
    bool checkInputFiles(const std::string& proj_path);

private:
    SimulationInfo* _sim_info;
    NumLib::ITransientCoupledSystem* _cpl_system;
};

} //end ogs6
