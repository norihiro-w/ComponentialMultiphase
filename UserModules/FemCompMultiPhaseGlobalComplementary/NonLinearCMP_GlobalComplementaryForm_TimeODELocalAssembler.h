/**
* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file NonLinearCMP_2P2C_TimeODELocalAssembler.h
*
* Created on 2015-06-19 by Yonghui HUANG
*/

#ifndef NON_LINEAR_CMP_GLOBALCOMPLEMENTARYFORM_TIME_ODE_LOCAL_ASSEMBLER_H
#define NON_LINEAR_CMP_GLOBALCOMPLEMENTARYFORM_TIME_ODE_LOCAL_ASSEMBLER_H

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "MaterialLib/PorousMedia.h"
#include "Ogs6FemData.h"
#include "EOS_GlobalComplementaryForm.h"




/**
* \brief Local assembly of time ODE components for linear transport in porous media
* This class is same as the MassTransportTimeODELocalAssembler class
* onlye difference is that no compound information is provided and no molecular
* diffusion is included in the assembly.
*/
template <class T, class T_NODAL_FUNCTION_SCALAR, class T_FUNCTION_DATA>
class NonLinearCMP_GlobalComplementaryForm_TimeODELocalAssembler : public T
{
public:
	typedef MathLib::LocalVector LocalVectorType;
	typedef MathLib::LocalMatrix LocalMatrixType;

	NonLinearCMP_GlobalComplementaryForm_TimeODELocalAssembler(FemLib::LagrangeFeObjectContainer* feObjects, T_FUNCTION_DATA* func_data, const MeshLib::CoordinateSystem &problem_coordinates)
		: _feObjects(*feObjects), _function_data(func_data), _problem_coordinates(problem_coordinates)
	{
		_EOS = new EOS_GlobalComplementaryForm();
	};
	
	virtual ~NonLinearCMP_GlobalComplementaryForm_TimeODELocalAssembler() {
		BaseLib::releaseObject(_EOS);
	};
	T_FUNCTION_DATA* get_function_data(void)
	{
		return _function_data;
	}

protected:
	
	virtual void assembleODE(const NumLib::TimeStep & /*time*/, const MeshLib::IElement &e, const LocalVectorType & u1, const LocalVectorType & u0, LocalMatrixType & localM, LocalMatrixType & localK, LocalVectorType & localF)
	{
		// -------------------------------------------------
		// current element
		FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
		// integration method 
		FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
		// now get the sizes
		const std::size_t n_dim = e.getDimension();
		const std::size_t mat_id = e.getGroupID();
		const std::size_t ele_id = e.getID();
		const std::size_t n_nodes = e.getNumberOfNodes(); // number of connecting nodes
		const std::size_t n_gsp = q->getNumberOfSamplingPoints(); // number of Gauss points
		std::size_t node_id(0);  // index of the node

		// get the instance of porous media class
		MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
		// get the fluid property for gas phase
		MaterialLib::Fluid* fluid1 = Ogs6FemData::getInstance()->list_fluid[0];
		// get the fluid property for liquid phase
		MaterialLib::Fluid* fluid2 = Ogs6FemData::getInstance()->list_fluid[1];
		// get the first (light) component - hydrogen
		MaterialLib::Compound* component1 = Ogs6FemData::getInstance()->list_compound[0];
		// get the second (heavy) component - water
		MaterialLib::Compound* component2 = Ogs6FemData::getInstance()->list_compound[1];

		const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e.getID());
		double geo_area = 1.0;
		pm->geo_area->eval(e_pos, geo_area);
		const bool hasGravityEffect =  _problem_coordinates.hasZ();//detect the gravity term based on the mesh structure

		if (hasGravityEffect) {
			vec_g = LocalVectorType::Zero(_problem_coordinates.getDimension());
			vec_g[_problem_coordinates.getIndexOfY()] = 9.81;
		}
		/*
		*SECONDARY VARIABLES initialize
		*/
		PC = LocalVectorType::Zero(n_nodes);
		rho_L_h = LocalVectorType::Zero(n_nodes);
		rho_G_h = LocalVectorType::Zero(n_nodes);
		dPC_dSg = LocalVectorType::Zero(n_nodes);

		rho_L_h_gp = LocalMatrixType::Zero(1, 1);
		rho_G_h_gp = LocalMatrixType::Zero(1, 1);
		PC_gp = LocalMatrixType::Zero(1, 1);
		dPC_dSg_gp = LocalMatrixType::Zero(1, 1);
		// Loop over each node on this element e
		// calculate the secondary variables' values on each node 
		// calculate the derivative of secondary variables based on P and X on each node
		// store the value in the vector respectively

		M = LocalMatrixType::Zero(2, 3);//method 2
		D = LocalMatrixType::Zero(2, 3);//method 2
		H = LocalVectorType::Zero(2);//for gravity term
		tmp = LocalMatrixType::Zero(1, 1);
		Input = LocalVectorType::Zero(3);
		localMass_tmp = LocalMatrixType::Zero(n_nodes, n_nodes); // for method2
		localDispersion_tmp = LocalMatrixType::Zero(n_nodes, n_nodes); //for method2
		localDispersion = LocalMatrixType::Zero(n_nodes, n_nodes);
		LocalJacb = LocalMatrixType::Zero(3*n_nodes, 3 * n_nodes);
		localGravity_tmp = LocalVectorType::Zero(n_nodes);//tmp vect for gravity 
		localRes = LocalVectorType::Zero(n_nodes);//tmp vect for gravity 
		/*
		* Loop for each Gauss Point
		*/
		Kr_L_gp = 0.0;
		Kr_G_gp = 0.0; lambda_L = 0.0; lambda_G = 0.0;
		Lambda_h = 0.0;
		//-------------------Begin assembly on each gauss point---------------------------
		for (j = 0; j < n_gsp; j++)
		{
			// calc on each gauss point
			// get the sampling point
			q->getSamplingPoint(j, gp_x);
			// compute basis functions
			fe->computeBasisFunctions(gp_x);
			// get the coordinates of current location 
			fe->getRealCoordinates(real_x);
			// transfer to Gauss point position
			NumLib::TXPosition gp_pos(NumLib::TXPosition::IntegrationPoint, e.getID(), j, real_x);
			double fac = geo_area*fe->getDetJ()*q->getWeight(j);
			// get the porosity
			pm->porosity->eval(gp_pos, poro);
			// get the intrinsic permeability
			pm->permeability->eval(gp_pos, K_perm);
			// evaluate viscosity of each fluid phase
			fluid1->dynamic_viscosity->eval(gp_pos, mu_G);
			fluid2->dynamic_viscosity->eval(gp_pos, mu_L);

			//fluid1->molar_mass->eval(gp_pos, M_G);
			//fluid2->molar_mass->eval(gp_pos, M_L);

			fluid1->density->eval(gp_pos, rho_G_std);
			fluid2->density->eval(gp_pos, rho_L_std);

			// evaluate componential diffusion coefficients
			component1->molecular_diffusion->eval(gp_pos, D_G);
			component2->molecular_diffusion->eval(gp_pos, D_L);

			// evaluation of the shape function
			LocalMatrixType &Np = *fe->getBasisFunction();
			LocalMatrixType &dNp = *fe->getGradBasisFunction();


			P_gp = Np*u1.head(n_nodes);
			X_gp = Np*u1.block(n_nodes, 0, n_nodes, 1);
			S_gp = Np*u1.tail(n_nodes);
			Input(0) = P_gp(0, 0);
			Input(1) = X_gp(0, 0);
			Input(2) = S_gp(0, 0);

			_EOS->set_env_condition(Input);//apply the primary variables as the input variables for Eos
			rho_L_h_gp(0, 0) = _EOS->get_RHO_L_H(S_gp(0, 0));
			rho_G_h_gp(0,0) = _EOS->get_RHO_G_H(S_gp(0, 0));
			PC_gp(0,0) = pm->getPc_bySat(S_gp(0, 0));
			dPC_dSg_gp(0,0) = pm->Deriv_dPCdS(S_gp(0, 0));
			isinf = std::isfinite(dPC_dSg_gp(0, 0));
			if (isinf == 0)
			{
				dPC_dSg_gp(0, 0) = 0.0;
			}
			else
			{
				dPC_dSg_gp(0, 0) = dPC_dSg_gp(0, 0);
			}

			/*if (C_h*(P_gp(0, 0) + PC_gp(0, 0)) <= X_gp(0,0)*M_G*rho_l_std / M_L){//then rho_L^h=C_h*PG
				//Calc each entry of the mass matrix
				M(0, 0) = poro*(1 - S_gp(0, 0))*C_h + poro*S_gp(0, 0)*C_v;
				M(0, 1) = 0.0;
				M(0, 2) = poro*(rho_G_h_gp(0, 0) - rho_L_h_gp(0, 0)) + poro*((1 - S_gp(0, 0))*C_h*dPC_dSg_gp(0, 0)) + poro*S_gp(0, 0)*C_v*dPC_dSg_gp(0, 0);
			}
			else if (C_h*(P_gp(0, 0) + PC_gp(0, 0)) > X_gp(0,0)*M_G*rho_l_std / M_L){//then rho_L^h=M^h*rho_L^std*X_L^h/M^w
				M(0, 0) = poro*S_gp(0, 0)*C_v;
				M(0, 1) = poro*(1 - S_gp(0, 0))*M_G*rho_l_std / M_L;
				M(0, 2) = poro*(rho_G_h_gp(0, 0) - rho_L_h_gp(0, 0)) + poro*S_gp(0, 0)*C_v*dPC_dSg_gp(0, 0);
			}*/
			M(0, 0) = poro*S_gp(0, 0)*C_v;
			M(0, 1) = poro*(1 - S_gp(0, 0))*M_G*rho_l_std / M_L;
			M(0, 2) = poro*(rho_G_h_gp(0, 0) - rho_L_h_gp(0, 0)) + poro*S_gp(0, 0)*C_v*dPC_dSg_gp(0, 0);
			M(1, 0) = 0.0;
			M(1, 1) = 0.0;
			M(1, 2) = -poro*rho_l_std;//
			
			//-------------debugging------------------------
			//std::cout << "M=" << std::endl;
			//std::cout << M << std::endl;
			//--------------end debugging-------------------

			//assembly the mass matrix
			for (ii = 0; ii < 2; ii++){
				for (jj = 0; jj < 3; jj++){
					tmp(0, 0) = M(ii, jj);
					localMass_tmp.setZero();
					fe->integrateWxN(j, tmp, localMass_tmp);
					localM.block(n_nodes*ii, n_nodes*jj, n_nodes, n_nodes) += localMass_tmp;
				}
			}
			//-------------debugging------------------------
			//std::cout << "localM=" << std::endl;
			//std::cout << localM << std::endl;
			//-------------End debugging------------------------
			//calculate the relative permeability 
			Kr_L_gp = pm->getKr_L_bySg(S_gp(0, 0)); // change the Relative permeability calculation into porous media file			
			Kr_G_gp = pm->getKr_g_bySg(S_gp(0, 0));
			lambda_L = K_perm*Kr_L_gp / mu_L;
			lambda_G = K_perm*Kr_G_gp / mu_G;
			/*if (C_h*(P_gp(0, 0) + PC_gp(0, 0)) <= X_gp(0, 0)*M_G*rho_l_std / M_L){//then rho_L^h=C_h*PG
				//Calc each entry of the mass matrix
				D(0, 0) = rho_L_h_gp(0)*lambda_L + rho_G_h_gp(0)*lambda_G + poro * (1 - S_gp(0))*D_L*C_h;
				D(0, 1) = 0.0;
				D(0, 2) = poro * (1 - S_gp(0))*D_L*C_h*dPC_dSg_gp(0,0);
				D(1, 0) = lambda_L*rho_l_std - poro*pow(1 - S_gp(0), 1)*D_L*C_h;
				D(1, 1) = 0.0;
				D(1, 2) = -poro * (1 - S_gp(0))*D_L*C_h*dPC_dSg_gp(0,0);
			}
			else if (C_h*(P_gp(0, 0) + PC_gp(0, 0)) > X_gp(0, 0)*M_G*rho_l_std / M_L){//then rho_L^h=M^h*rho_L^std*X_L^h/M^w
				D(0, 0) = rho_L_h_gp(0)*lambda_L + rho_G_h_gp(0)*lambda_G;
				D(0, 1) = poro*M_G*pow(1 - S_gp(0), 1)*rho_L_std*D_L / (M_L);
				D(0, 2) = rho_G_h_gp(0)*lambda_G*dPC_dSg_gp(0,0);
				D(1, 0) = lambda_L*rho_l_std;
				D(1, 1) = -poro*M_G*pow(1 - S_gp(0), 1)*rho_L_std*D_L / (M_L);
				D(1, 2) = 0.0;
			}*/
			
			//Calc each entry of the Laplace Matrix
			D(0, 0) = rho_L_h_gp(0, 0)*lambda_L + rho_G_h_gp(0, 0)*lambda_G;
			D(0, 1) = poro*M_G*pow(1 - S_gp(0, 0), 1)*rho_L_std*D_L *(1 / 1 + X_gp(0,0)) / (M_L);
			D(0, 2) = rho_G_h_gp(0, 0)*lambda_G*dPC_dSg_gp(0,0);
			D(1, 0) = lambda_L*rho_l_std;
			D(1, 1) = -poro*M_G*pow(1 - S_gp(0, 0), 1)*rho_L_std*D_L *(1 / 1 + X_gp(0, 0)) / (M_L);
			D(1, 2) = 0.0;
			//-------------debugging------------------------
			//std::cout << "D=" << std::endl;
			//std::cout << D << std::endl;
			//--------------end debugging-------------------
			//
			for (ii = 0; ii < 2; ii++){
				for (jj = 0; jj < 3; jj++){
					tmp(0, 0) = D(ii, jj);
					localDispersion_tmp.setZero();
					fe->integrateDWxDN(j, tmp, localDispersion_tmp);
					localK.block(n_nodes*ii, n_nodes*jj, n_nodes, n_nodes) += localDispersion_tmp;
				}
			}
			//-------------debugging------------------------
			//std::cout << "localK=" << std::endl;
			//std::cout << localK << std::endl;
			//--------------end debugging-------------------
			//----------------assembly the gravity term--------------------------------
			H(0) = -rho_L_h_gp(0, 0)*lambda_L*(rho_L_h_gp(0, 0) + rho_L_std) - rho_G_h_gp(0, 0)*rho_G_h_gp(0,0)*lambda_G; //-pow(RHO_G, 2)*lambda_G - pow(RHO_L, 2)*lambda_L;
			H(1) = -rho_l_std*lambda_L*(rho_L_h_gp(0, 0) + rho_L_std);
			//-------------debugging------------------------
			//std::cout << "H=" << std::endl;
			//std::cout << H << std::endl;
			//--------------end debugging-------------------
			if (hasGravityEffect) {
				// F += dNp^T * H* gz
				
				for (int idx = 0; idx < 2; idx++){
					tmp(0, 0) = H(idx);
					localGravity_tmp.setZero();
					//fe->integrateDWxvec_g(j, tmp, localGravity_tmp, vec_g);
					localGravity_tmp = fac * dNp.transpose()*tmp(0,0)*vec_g;
					localF.block(n_nodes*idx, 0, n_nodes, 1) += localGravity_tmp;

				}
			}
			//----------------end assembly--------------------------------
			//std::cout << "H=" << std::endl;
			//std::cout << localF << std::endl;

		}
		//----------------end assemby gaus---------------------------------
		//----------------calc jacobian for complementary contrains--------------
		/*
		for (i = 0; i < n_nodes; i++)
		{
			node_id = e.getNodeID(i);
			localF[i + 2 * n_nodes] = -std::min(u1[i+2*n_nodes],C_h*(u1[i] + PC[i]) - u1[i + n_nodes] * M_G*rho_L_std / M_L);
			if (u1[i + 2 * n_nodes] <=C_h*(u1[i] + PC[i]) - u1[i + n_nodes] * M_G*rho_l_std / M_L){//then rho_L^h=C_h*PG
			//Calc each entry of the mass matrix
				LocalJacb(2 * n_nodes+i, 2 * n_nodes + i) = 1;
			}
			else if (u1[i + 2 * n_nodes] > C_h*(u1[i] + PC[i]) - u1[i + n_nodes] * M_G*rho_l_std / M_L){//then rho_L^h=M^h*rho_L^std*X_L^h/M^w
				LocalJacb(i + 2 * n_nodes, i) = C_h;
				LocalJacb(i + 2 * n_nodes, n_nodes + i) = -M_G*rho_L_std / M_L;
				LocalJacb(i + 2 * n_nodes, 2 * n_nodes + i) = C_h*dPC_dSg[i];
			}
		}
		*/
		//std::cout << "H=" << std::endl;
		//std::cout << LocalJacb<< std::endl;
		_function_data->get_elem_M_matrix()[ele_id] = localM;
		_function_data->get_elem_K_matrix()[ele_id] = localK;
		
		//_function_data->get_elem_J_matrix()[ele_id] = LocalJacb;
	}

private:
	FemLib::LagrangeFeObjectContainer _feObjects;
	MeshLib::CoordinateSystem _problem_coordinates;
	/**
	  * pointer to the function data class
	  */
	T_FUNCTION_DATA* _function_data;
	std::size_t i, j, ii, jj;

	double real_x[3], gp_x[3];

	LocalVectorType vec_g; // gravity term 
	LocalVectorType H;//this is for local gravity vector 
	LocalVectorType localGravity_tmp;
	LocalVectorType localRes;

	LocalVectorType PC;
	LocalVectorType rho_L_h;
	LocalVectorType rho_G_h;
	LocalVectorType dPC_dSg;

	LocalVectorType Input;

	
	LocalMatrixType rho_L_h_gp;
	LocalMatrixType rho_G_h_gp;
	LocalMatrixType PC_gp;
	LocalMatrixType dPC_dSg_gp;
	//LocalMatrixType PG_gp;
	//LocalMatrixType PL_gp;

	//primary variable on gauss point
	LocalMatrixType P_gp;
	LocalMatrixType X_gp;
	LocalMatrixType S_gp;


	LocalMatrixType M;//Local Mass Matrix
	LocalMatrixType D;//Local Laplace Matrix
	LocalMatrixType tmp;//method 2


	LocalMatrixType localMass_tmp; //for method2
	LocalMatrixType localDispersion_tmp; //for method2
	LocalMatrixType localDispersion;
	LocalMatrixType matN;
	LocalMatrixType LocalJacb;


	double Kr_L_gp, Kr_G_gp;
	double lambda_L, lambda_G;
	double poro, K_perm, mu_L, mu_G, D_G, D_L;
	double rho_G_std;
	double rho_L_std;
	double Lambda_h;
	double Var_a;
	const double M_L = 0.01;
	const double M_G = 0.002;

	double RHO_L;// RHO_L=rho_L_std+rho_L_h
	double RHO_G;// RHO_G=rho_G_h+rho_G_w

	const double C_h = 1.53e-8; 
	const double C_v = 7.938638137741564e-07;
	const double C_w = 3.969319068870782e-06;//M_L/RT
	const double rho_l_std = 1000.0;
	EOS_GlobalComplementaryForm * _EOS;
	std::size_t isinf;
};



#endif  // end of ifndef