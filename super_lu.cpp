/*
 * super_lu.cpp
 *
 *  Created on: 16.01.2014
 *      Author: mrupp
 */

#include "super_lu.h"

// fix warning
#undef TRUE
#undef FALSE
#include "slu_ddefs.h"

namespace ug{


class SuperLUImplementation : public IExternalSolverImplementation
{
	bool m_bInited;
	std::vector<double> rhs, nzval;
	/* row permutations from partial pivoting */
	/* column permutation vector */
	std::vector<int> perm_r, perm_c, colind, rowptr;

	SuperMatrix SuperLU_A, SuperLU_L, SuperLU_U, SuperLU_B;

	SuperLUConfiguration &config;

public:
	SuperLUImplementation(SuperLUConfiguration &_config) : config(_config), trans(NOTRANS)
	{
		m_bInited=false;
	}

	~SuperLUImplementation()
	{
		destroy();
	}



	void destroy()
	{
		if(m_bInited)
		{
			// DO NOT destroy SuperLU_A and AA since we supplied all pointers!!!
			Destroy_SuperMatrix_Store(&SuperLU_B);
			memset(&SuperLU_B, 0, sizeof(SuperMatrix));
			Destroy_SuperNode_Matrix(&SuperLU_L);
			memset(&SuperLU_L, 0, sizeof(SuperMatrix));
			Destroy_CompCol_Matrix(&SuperLU_U);
			memset(&SuperLU_U, 0, sizeof(SuperMatrix));

		    SUPERLU_FREE (etree);
		    Destroy_CompCol_Permuted(&AC);

			m_bInited = false;

		}
	}

	void get_options(superlu_options_t opt)
	{
		//opt.PrintStat = config.bPrintStat ? YES : NO;
		opt.Equil = config.equil ? YES : NO;
		switch(config.colPerm)
		{
		case SuperLUConfiguration::CPT_NATURAL: opt.ColPerm = NATURAL; break;
		case SuperLUConfiguration::CPT_MMD_ATA: opt.ColPerm = MMD_ATA; break;
		case SuperLUConfiguration::CPT_MMD_AT_PLUS_A: opt.ColPerm = MMD_AT_PLUS_A; break;
		case SuperLUConfiguration::CPT_COLAMD: opt.ColPerm = COLAMD; break;
		}
	}

	superlu_options_t options;



	int *etree;
	SuperMatrix *AA;/* A in SLU_NC format used by the factorization routine.*/
	SuperMatrix AC; /* Matrix postmultiplied by Pc */
	trans_t  trans;

	void
	dgssvA(superlu_options_t *options, SuperMatrix *A, int *perm_c, int *perm_r,
	      SuperMatrix *L, SuperMatrix *U, SuperMatrix *B,
	      SuperLUStat_t *stat, int *info )
	{

	    //DNformat *Bstore;

	    int      lwork = 0;

	    /* Set default values for some parameters */
	    int      panel_size;     /* panel size */
	    int      relax;          /* no of columns in a relaxed snodes */
	    int      permc_spec;
	     trans = NOTRANS;


	     /* Convert A to SLU_NC format when necessary. */
	    if ( A->Stype == SLU_NR ) {
		NRformat *Astore = (NRformat*)A->Store;
		AA = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
		dCreate_CompCol_Matrix(AA, A->ncol, A->nrow, Astore->nnz,
				(double*)Astore->nzval, Astore->colind, Astore->rowptr,
				       SLU_NC, A->Dtype, A->Mtype);
		trans = TRANS;
	    } else {
	        if ( A->Stype == SLU_NC ) AA = A;
	    }

	    /*
	     * Get column permutation vector perm_c[], according to permc_spec:
	     *   permc_spec = NATURAL:  natural ordering
	     *   permc_spec = MMD_AT_PLUS_A: minimum degree on structure of A'+A
	     *   permc_spec = MMD_ATA:  minimum degree on structure of A'*A
	     *   permc_spec = COLAMD:   approximate minimum degree column ordering
	     *   permc_spec = MY_PERMC: the ordering already supplied in perm_c[]
	     */
	    permc_spec = options->ColPerm;
	    if ( permc_spec != MY_PERMC && options->Fact == DOFACT )
	      get_perm_c(permc_spec, AA, perm_c);

	    etree = intMalloc(A->ncol);

	    sp_preorder(options, AA, perm_c, etree, &AC);

	    panel_size = sp_ienv(1);
	    relax = sp_ienv(2);

	    /* Compute the LU factorization of A. */
	    dgstrf(options, &AC, relax, panel_size, etree,
	            NULL, lwork, perm_c, perm_r, L, U, stat, info);
	}
	void
	dgssvB(superlu_options_t *options, SuperMatrix *A, int *perm_c, int *perm_r,
	      SuperMatrix *L, SuperMatrix *U, SuperMatrix *B,
	      SuperLUStat_t *stat, int *info )
	{
	    dgstrs (trans, L, U, perm_c, perm_r, B, stat, info);
	}



	virtual bool init(const CPUAlgebra::matrix_type &A)
	{

		destroy();
		PROFILE_BEGIN_GROUP(SuperLU_Preprocess, "algebra SuperLU");
		typedef CPUAlgebra::matrix_type::const_row_iterator row_it;
		typedef CPUAlgebra::matrix_type::value_type value_type;

		if( A.num_rows() == 0 || A.num_cols() == 0) return true;

		size_t numRows, numCols;
		A.copy_crs(numRows, numCols, nzval, rowptr, colind);
		THROW_IF_NOT_EQUAL(numRows, numCols);
		size_t N = numRows;
		size_t nnz = nzval.size();
		if(N > 40000 && nnz > 400000) { UG_LOG("SuperLU preprocess, N = " << N << ", nnz = " << nnz << "... "); }



		SuperLUStat_t stat;

		dCreate_CompRow_Matrix(&SuperLU_A, N, N, nnz, &nzval[0], &colind[0], &rowptr[0], SLU_NR, SLU_D, SLU_GE);

		//rhs = doubleMalloc(N);
		rhs.resize(N);
		for (size_t i = 0; i < N; ++i)
			rhs[i] = 0;
		dCreate_Dense_Matrix(&SuperLU_B, N, 1, &rhs[0], N, SLU_DN, SLU_D, SLU_GE);

		perm_r.resize(N+1);
		perm_c.resize(N+1);
		/* Set the default input options. */
		set_default_options(&options);

		get_options(options);

		StatInit(&stat);
		int info;

		dgssvA(&options, &SuperLU_A, &perm_c[0], &perm_r[0], &SuperLU_L, &SuperLU_U, &SuperLU_B, &stat, &info);

		/*dgssv_check_info(info, N);
		if(config.bPrintStat)
			StatPrint(&stat);
		StatFree(&stat);*/
		if(N > 40000 && nnz > 400000) { UG_LOG("done.\n"); }
		m_bInited = true;
		return true;
	}

	void dgssv_check_info(int info, size_t N)
	{
		if(info > 0)
		{
			if(info < (int)N)
			{
				UG_THROW("ERROR in SuperLU: U(i,i) with i=" << info << "is exactly zero. The factorization has\
							been completed, but the factor U is exactly singular,\
							so the solution could not be computed.");
			}
			else
			{ UG_THROW("ERROR in SuperLU: memory allocation failure");}
		}
		else if(info < 0)
		{
			UG_THROW("ERROR in SuperLU: info < 0 ???");
		}
	}

	virtual bool apply(CPUAlgebra::vector_type &c, const CPUAlgebra::vector_type &d)
	{
		PROFILE_BEGIN_GROUP(SuperLU_Apply, "algebra SuperLU");
		size_t N = c.size();
		if(N == 0) return true;
		double *b = (double*) ((DNformat*) SuperLU_B.Store)->nzval;

		for (size_t i = 0; i < N; ++i)
			b[i] = d[i];

		superlu_options_t options;
		set_default_options(&options);
		options.Fact = FACTORED;
		int info;
		SuperLUStat_t stat;
		StatInit(&stat);
		dgssvB(&options, &SuperLU_A, &perm_c[0], &perm_r[0], &SuperLU_L, &SuperLU_U, &SuperLU_B, &stat, &info);
		StatFree(&stat);
		dgssv_check_info(info, N);

		for (size_t i = 0; i < N; ++i)
			c[i] = b[i];

		return true;
	}

	virtual const char* name() const { return "SuperLU"; }
};


IExternalSolverImplementation *CreateSuperLUImplementation(SuperLUConfiguration &config)
{
	return new SuperLUImplementation(config);
}

}

