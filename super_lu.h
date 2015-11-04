/*
 * http://crd-legacy.lbl.gov/~xiaoye/SuperLU
 */

#ifndef __H__UG__SUPER_LU_SOLVER__
#define __H__UG__SUPER_LU_SOLVER__

#include "lib_algebra/operator/linear_solver/external_solvers/external_solvers.h"
#include <vector>


namespace ug{

struct SuperLUConfiguration
{
	bool bPrintStat;
	bool equil;
	enum
	{
	CPT_NATURAL, CPT_MMD_ATA, CPT_MMD_AT_PLUS_A, CPT_COLAMD
	} colPerm;
};


IExternalSolverImplementation *CreateSuperLUImplementation(SuperLUConfiguration &config);


template<typename TAlgebra>
class SuperLUSolver : public IExternalSolver<TAlgebra>
{
	SuperLUConfiguration config;
	IExternalSolverImplementation *impl;
public:
	SuperLUSolver()
	{
		config.bPrintStat = false;
		config.equil = true;
		config.colPerm = SuperLUConfiguration::CPT_COLAMD;

		impl = CreateSuperLUImplementation(config);
	}

	void print_stat(bool b)
	{
		config.bPrintStat = b;
	}

	void equil(bool b)
	{
		config.equil = b;
	}

	void col_perm_natural()
	{
		config.colPerm = SuperLUConfiguration::CPT_NATURAL;
	}

	void col_perm_mdo_ATA()
	{
		config.colPerm = SuperLUConfiguration::CPT_MMD_ATA;
	}

	void col_perm_mdo_AT_plus_A()
	{
		config.colPerm = SuperLUConfiguration::CPT_MMD_AT_PLUS_A;
	}

	void col_perm_approximate()
	{
		config.colPerm = SuperLUConfiguration::CPT_COLAMD;
	}


	virtual bool double_apply(CPUAlgebra::vector_type &c, const CPUAlgebra::vector_type &d)
	{
		return impl->apply(c, d);
	}

	virtual const char *double_name() const
	{
		return impl->name();
	}

	virtual void double_init(const CPUAlgebra::matrix_type &mat)
	{
		impl->init(mat);
	}

};

} // end namespace ug

#endif

