/* Copyright (C) 2002-2012 The SSI Project. All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
   3. Neither the name of the project nor the names of its contributors 
      may be used to endorse or promote products derived from this software 
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE SCALABLE SOFTWARE INFRASTRUCTURE PROJECT
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE SCALABLE SOFTWARE INFRASTRUCTURE
   PROJECT BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
   OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_CONFIG_H
	#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN32_H
	#include "lis_config_win32.h"
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
        #include <malloc.h>
#endif
#include <math.h>
#include <string.h>
#include <stdarg.h>
#ifdef USE_SSE2
	#include <emmintrin.h>
#endif
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

#define NWORK 3
#undef __FUNC__
#define __FUNC__ "lis_erqi_check_params"
LIS_INT lis_erqi_check_params(LIS_ESOLVER esolver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_erqi_malloc_work"
LIS_INT lis_erqi_malloc_work(LIS_ESOLVER esolver)
{
	LIS_VECTOR	*work;
	LIS_INT			i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_erqi_malloc_work::work" );
	if( work==NULL )
	{
		LIS_SETERR_MEM(worklen*sizeof(LIS_VECTOR));
		return LIS_ERR_OUT_OF_MEMORY;
	}
	if( esolver->eprecision==LIS_PRECISION_DEFAULT )
	{
		for(i=0;i<worklen;i++)
		{
			err = lis_vector_duplicate(esolver->A,&work[i]);
			if( err ) break;
		}
	}
	else
	{
		for(i=0;i<worklen;i++)
		{
			err = lis_vector_duplicateex(LIS_PRECISION_QUAD,esolver->A,&work[i]);
			if( err ) break;
		}
	}
	if( i<worklen )
	{
		for(j=0;j<i;j++) lis_vector_destroy(work[j]);
		lis_free(work);
		return err;
	}
	esolver->worklen = worklen;
	esolver->work    = work;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_erqi"
LIS_INT lis_erqi(LIS_ESOLVER esolver)
{
  LIS_MATRIX        A;
  LIS_VECTOR        x, Ax;
  LIS_SCALAR        evalue, ievalue;
  LIS_SCALAR        xAx, xx, mu, lshift;
  LIS_INT               emaxiter;
  LIS_REAL          tol;
  LIS_INT               iter,iter2,i,n,gn,is,ie,output;
  LIS_INT               nprocs,my_rank;
  LIS_REAL          nrm2,resid;
  LIS_VECTOR        z,q;
  LIS_SOLVER        solver;
  double	    times,itimes,ptimes,p_c_times,p_i_times;

  LIS_INT		    err;
  LIS_PRECON        precon;
  LIS_INT           nsol, precon_type;
  char              solvername[128], preconname[128];

  LIS_DEBUG_FUNC_IN;

  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  lshift = esolver->lshift;
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];

  A = esolver->A;
  x = esolver->x;
  if (esolver->options[LIS_EOPTIONS_INITGUESS_ONES] ) 
    {
      lis_vector_set_all(1.0,x);
    }
  evalue = 1.0;
  z = esolver->work[0];
  q = esolver->work[1];
  Ax = esolver->work[2];

  iter=0;
  ievalue = 1/(evalue);
  if( A->my_rank==0 ) printf("local shift = %e\n", lshift);
  if (lshift != 0) lis_matrix_shift_diagonal(A, lshift);
  lis_solver_create(&solver);
  lis_solver_set_option("-p ilu -maxiter 10",solver);
  lis_solver_set_optionC(solver);
  lis_solver_get_solver(solver, &nsol);
  lis_solver_get_precon(solver, &precon_type);
  lis_get_solvername(nsol, solvername);
  lis_get_preconname(precon_type, preconname);
  printf("solver     : %s %d\n", solvername, nsol);
  printf("precon     : %s %d\n", preconname, precon_type);

  /* create preconditioner */
  solver->A = A;
  err = lis_precon_create(solver, &precon);
  if( err )
    {
      lis_solver_work_destroy(solver);
      solver->retcode = err;
      return err;
    }

  lis_vector_nrm2(x, &nrm2);
  lis_vector_scale(1/nrm2, x);
  lis_matvec(A, x, Ax);
  lis_vector_dot(x, Ax, &xAx);
  lis_vector_dot(x, x, &xx);
  mu = xAx / xx;

  while (iter<emaxiter)
    {
      iter = iter+1;
      lis_vector_nrm2(x, &nrm2);
      lis_vector_scale(1/nrm2, x);
      lis_matrix_shift_diagonal(A, -mu);
      lis_solve_kernel(A, x, z, solver, precon);
      lis_matrix_shift_diagonal(A, mu);
      lis_solver_get_iters(solver,&iter2);
      lis_vector_dot(x, z, &ievalue); 
      mu = mu + 1/ievalue;
      lis_vector_axpyz(-ievalue,x,z,q); 
      lis_vector_nrm2(q, &resid); 
      resid = fabs(resid/ievalue);
      if( output )
	{
	  if( output & LIS_EPRINT_MEM ) esolver->residual[iter] = resid;
	  if( output & LIS_EPRINT_OUT && A->my_rank==0 ) printf("iter: %5d  residual = %e\n", iter, resid);
	}
      lis_vector_copy(z,x);
      lis_solver_get_timeex(solver,&times,&itimes,&ptimes,&p_c_times,&p_i_times);
      esolver->ptimes += solver->ptimes;
      esolver->itimes += solver->itimes;
      esolver->p_c_times += solver->p_c_times;
      esolver->p_i_times += solver->p_i_times;
      if( tol >= resid ) 
	{
	  esolver->retcode    = LIS_SUCCESS;
	  esolver->iter       = iter;
	  esolver->resid      = resid;
	  esolver->evalue[0] = mu;
	  if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
	  lis_solver_destroy(solver); 
	  LIS_DEBUG_FUNC_OUT;
	  return LIS_SUCCESS;
	}
    }

  lis_precon_destroy(precon);

  esolver->retcode    = LIS_MAXITER;
  esolver->iter       = iter;
  esolver->resid      = resid;
  esolver->evalue[0] = mu;
  if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
  lis_solver_destroy(solver); 
  LIS_DEBUG_FUNC_OUT;
  return LIS_MAXITER;
}

#ifdef USE_QUAD_PRECISION
#undef __FUNC__
#define __FUNC__ "lis_erqi_quad"
LIS_INT lis_erqi_quad(LIS_ESOLVER esolver)
{
  LIS_MATRIX        A;
  LIS_VECTOR        x, Ax;
  LIS_SCALAR        evalue, ievalue;
  LIS_SCALAR        xAx, xx, mu, lshift;
  LIS_INT               emaxiter;
  LIS_REAL          tol;
  LIS_INT               iter,iter2,i,n,gn,is,ie,output;
  LIS_INT               nprocs,my_rank;
  LIS_REAL          nrm2,resid;
  LIS_QUAD_PTR      qdot_xz;
  LIS_VECTOR        z,q;
  LIS_SOLVER        solver;
  double	    times,itimes,ptimes,p_c_times,p_i_times;

  LIS_INT		    err;
  LIS_PRECON        precon;

  LIS_DEBUG_FUNC_IN;

  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  lshift = esolver->lshift;
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];

  A = esolver->A;
  x = esolver->x;
  if (esolver->options[LIS_EOPTIONS_INITGUESS_ONES] ) 
    {
      lis_vector_set_all(1.0,x);
    }
  evalue = 1.0;
  z = esolver->work[0];
  q = esolver->work[1];
  Ax = esolver->work[2];

  LIS_QUAD_SCALAR_MALLOC(qdot_xz,0,1);

  iter=0;
  ievalue = 1/(evalue);
  if( A->my_rank==0 ) printf("local shift = %e\n", lshift);
  if (lshift != 0) lis_matrix_shift_diagonal(A, lshift);
  lis_solver_create(&solver);
  lis_solver_set_option("-p ilu -precision quad -maxiter 10",solver);
  lis_solver_set_optionC(solver);

  /* create preconditioner */
  solver->A = A;
  err = lis_precon_create(solver, &precon);
  if( err )
    {
      lis_solver_work_destroy(solver);
      solver->retcode = err;
      return err;
    }

  lis_vector_nrm2(x, &nrm2);
  lis_vector_scale(1/nrm2, x);
  lis_matvec(A, x, Ax);
  lis_vector_dot(x, Ax, &xAx);
  lis_vector_dot(x, x, &xx);
  mu = xAx / xx;

  while (iter<emaxiter)
    {
      iter = iter+1;
      lis_vector_nrm2(x, &nrm2);
      lis_vector_scale(1/nrm2, x);
      lis_matrix_shift_diagonal(A, -mu);
      lis_solve_kernel(A, x, z, solver, precon);
      lis_matrix_shift_diagonal(A, mu);
      lis_solver_get_iters(solver,&iter2);
      lis_vector_dotex_mmm(x, z, &qdot_xz);
      lis_quad_minus((LIS_QUAD *)qdot_xz.hi);
      lis_vector_axpyzex_mmmm(qdot_xz,x,z,q);
      lis_quad_minus((LIS_QUAD *)qdot_xz.hi);
      ievalue = qdot_xz.hi[0];
      mu = mu + 1/ievalue;
      lis_vector_axpyz(-ievalue,x,z,q); 
      lis_vector_nrm2(q, &resid); 
      resid = fabs(resid/ievalue);

      if( output )
	{
	  if( output & LIS_EPRINT_MEM ) esolver->residual[iter] = resid;
	  if( output & LIS_EPRINT_OUT && A->my_rank==0 ) printf("iter: %5d  residual = %e\n", iter, resid);
	}
      lis_vector_copy(z,x);
      lis_solver_get_timeex(solver,&times,&itimes,&ptimes,&p_c_times,&p_i_times);
      esolver->ptimes += solver->ptimes;
      esolver->itimes += solver->itimes;
      esolver->p_c_times += solver->p_c_times;
      esolver->p_i_times += solver->p_i_times;
      if( tol >= resid ) 
	{
	  esolver->retcode    = LIS_SUCCESS;
	  esolver->iter       = iter;
	  esolver->resid      = resid;
	  esolver->evalue[0] = mu;
	  if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
	  lis_solver_destroy(solver); 
	  LIS_DEBUG_FUNC_OUT;
	  return LIS_SUCCESS;
	}
    }

  lis_precon_destroy(precon);

  esolver->retcode    = LIS_MAXITER;
  esolver->iter       = iter;
  esolver->resid      = resid;
  esolver->evalue[0] = mu;
  if (lshift != 0) 
    {
      lis_matrix_shift_diagonal(A, -lshift);
    }
  lis_solver_destroy(solver); 
  LIS_DEBUG_FUNC_OUT;
  return LIS_MAXITER;
}
#endif
