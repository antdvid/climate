/*******************************************************************
 * 		VCARTESIAN.c
 *******************************************************************/

#include <iFluid.h>
#include "solver.h"
#include "climate.h"
static void setSamplePoints(double*, double*, int);
static int find_state_at_crossing(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);
//----------------------------------------------------------------
//              RECTANGLE
//----------------------------------------------------------------

//RECTANGLE::RECTANGLE()
RECTANGLE::RECTANGLE(): index(-1), comp(-1)
{
}

void RECTANGLE::setCoords(
        double *crds,
        int dim)
{
        int i;
        for (i = 0; i < dim; ++i)
            coords[i] = crds[i];
}

// 		VCARTESIAN
//--------------------------------------------------------------------------

VCARTESIAN::~VCARTESIAN()
{
}

//---------------------------------------------------------------
//	initMesh
// include the following parts
// 1) setup edges:		
// 2) setup cell_center
//---------------------------------------------------------------
void VCARTESIAN::initMesh(void)
{
	int i,j,k,index;
	double crds[MAXD];
	int icoords[MAXD];
	int num_cells;
	int cell_index;
	RECTANGLE       rectangle;

	// init vertices,edges & cell_center
	if (debugging("trace")) printf("Entering initMesh()\n");
	FT_MakeGridIntfc(front);
	setDomain();

	num_cells = 1;
	for (i = 0; i < dim; ++i)
	{
	    num_cells *= (top_gmax[i] + 1);
	}
	cell_center.insert(cell_center.end(),num_cells,rectangle);
	
	// setup vertices
	// left to right, down to up
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	crds[0] = top_L[0] + top_h[0]*i;
		index = d_index1d(i,top_gmax);
	    	cell_center[index].setCoords(crds,dim);
	    	cell_center[index].icoords[0] = i;
	    }
	case 2:
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	crds[0] = top_L[0] + top_h[0]*i;
	    	crds[1] = top_L[1] + top_h[1]*j;
		index = d_index2d(i,j,top_gmax);
	    	cell_center[index].setCoords(crds,dim);
	    	cell_center[index].icoords[0] = i;
	    	cell_center[index].icoords[1] = j;
	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	crds[0] = top_L[0] + top_h[0]*i;
	    	crds[1] = top_L[1] + top_h[1]*j;
	    	crds[2] = top_L[2] + top_h[2]*k;
		index = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].setCoords(crds,dim);
	    	cell_center[index].icoords[0] = i;
	    	cell_center[index].icoords[1] = j;
	    	cell_center[index].icoords[2] = k;
	    }
	}

	setComponent();
	FT_FreeGridIntfc(front);
	if (debugging("trace")) printf("Leaving initMesh()\n");
}

void VCARTESIAN::setComponent(void)
{
	int i;

        static STATE *state = NULL;
        double *coords;
        int *icoords;
        int size = (int)cell_center.size();
        HYPER_SURF_ELEMENT *hse;
        HYPER_SURF *hs;
        double t[MAXD],point[MAXD];
        int n;

        if (state == NULL)
            FT_ScalarMemoryAlloc((POINTER*)&state, sizeof(STATE));

        for (i = 0; i < size; i++)
        {
            coords = cell_center[i].coords;
            if (cell_center[i].comp != -1 &&
                cell_center[i].comp != top_comp[i])
            {
                if (FT_FindNearestIntfcPointInRange(front,top_comp[i],coords,
                                INCLUDE_BOUNDARIES,point,t,&hse,&hs,2))
                {
                    if (!FrontNearestIntfcState(front,coords,top_comp[i],
                                (POINTER)state))
                    {
                        (void) printf("In setComponent()\n");
                        (void) printf("FrontNearestIntfcState() failed\n");
                        (void) printf("old_comp = %d new_comp = %d\n",
                                        cell_center[i].comp,top_comp[i]);
                        clean_up(ERROR);
                    }
                    field->vapor[i] = state->vapor;
                }
		else
                {
                    double temp_nb = 0.0;
                    int ii,jj,ic[MAXD],index;
                    icoords = cell_center[i].icoords;
                    n = 0;
                    for (ii = 0; ii < dim; ++ii)
                    {
                        for (jj = 0; jj < dim; ++jj) ic[jj] = icoords[jj];
                        ic[ii] = (icoords[ii] == 0) ? 0 : icoords[ii] - 1;
                        index = d_index(ic,top_gmax,dim);
                        if (cell_center[index].comp != -1 &&
                            cell_center[index].comp == top_comp[index])
                        {
                            temp_nb += field->vapor[index];
                            n++;
                        }
                        ic[ii] = (icoords[ii] == top_gmax[ii]) ? top_gmax[ii]
                                        : icoords[ii] + 1;
                        index = d_index(ic,top_gmax,dim);
                        if (cell_center[index].comp != -1 &&
                            cell_center[index].comp == top_comp[index])
                        {
                            temp_nb += field->vapor[index];
                            n++;
                        }
                    }
                    field->vapor[i] = temp_nb/n;
                }
	    }
	    cell_center[i].comp = top_comp[i];
        }
}	/* end setComponent */

static double computeVolumeMean(double* var, Front* front)
{
	int i;
	double mean = 0.0,dev;
	Deviation(var,front,mean,dev);
	return mean;
}

static double computeSatuVapor(double temp, double pres)
{
	double sat_vap_pre, sat_vap_rat;

	sat_vap_pre = 611.2*exp(17.67*(temp-273.15)/(temp-29.65));
        sat_vap_rat = 621.97 * sat_vap_pre/(pres-sat_vap_pre);
	return sat_vap_rat;
}

void VCARTESIAN::setInitialCondition(void)
{
	int i;
	double coords[MAXD], T1;
	INTERFACE *intfc = front->interf;
	POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	STATE *sl,*sr;
	PARAMS* eqn_params = (PARAMS*)front->extra2;
	IF_PARAMS* iFparams = (IF_PARAMS*)front->extra1;
	int c;
	
	short unsigned int seed[3] = {2,72,7172};
        GAUSS_PARAMS gauss_params;
        gauss_params.mu = eqn_params->qe;
        gauss_params.sigma = 0.2;

	FT_MakeGridIntfc(front);
        setDomain();
	eqn_params->qs = computeSatuVapor(eqn_params->T0,iFparams->ref_pres);

	/* Initialize states at the interface */
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            sr->vapor = eqn_params->qe;
	    sl->vapor = eqn_params->qs * 1.02;
        }

	// cell_center
	for (i = 0; i < cell_center.size(); i++)
	{
	    field->vapor[i] = 0;
	    field->supersat[i] = 0;
	    field->temperature[i] = 0;
	}
	for (i = 0; i < cell_center.size(); i++)
	{
	    getRectangleCenter(i,coords);
	    c = top_comp[i];
	    if (c == LIQUID_COMP2)
	    {
		if (eqn_params->init_vapor_state == RAND_STATE)
	    	    field->vapor[i] = gauss_center_limit(
				      (POINTER)&gauss_params,seed);
		else if (eqn_params->init_vapor_state == CONST_STATE)
		    field->vapor[i] = eqn_params->qe;
		else if (getInitialVaporState != NULL)
		    getInitialVaporState(&c,coords,field,i,dim,eqn_params);
		else
		{
		    printf("Unknown initial state %d for vapor\n",
			   eqn_params->init_vapor_state);
		    clean_up(ERROR);
		}
		if (getInitialTempState != NULL)
                    getInitialTempState(&c,coords,field,i,dim,eqn_params);
                else
                    field->temperature[i] = eqn_params->T0;
	    }
	}
	FT_ParallelExchGridArrayBuffer(field->vapor,front,NULL);
	FT_ParallelExchGridArrayBuffer(field->temperature,front,NULL);
	eqn_params->qv0 = computeVolumeMean(field->vapor,front);
	eqn_params->T0 = computeVolumeMean(field->temperature,front);
        printf("T0 = %f\n", eqn_params->T0);
	printf("qs = %f, qv0 = %20.14f\n",eqn_params->qs,eqn_params->qv0);
}	/* end setInitialCondition */


void VCARTESIAN::setParallelVapor(void)
{
        FILE *infile;
        int i,j,id,k,l,index,G_index;
        char fname[100];
        COMPONENT comp;
        double coords[MAXD];
        int size = (int)cell_center.size();
        int myid = pp_mynode();
        int numprocs = pp_numnodes();

        int G_icoords[MAXD],pp_icoords[MAXD],icoords[MAXD];
        int local_gmax[MAXD], global_gmax[MAXD];
        int G_size, L_size;
        PP_GRID *pp_grid = front->pp_grid;
        double *local_L = pp_grid->Zoom_grid.L;
        double *local_U = pp_grid->Zoom_grid.U;
        double *GV_buff, *V_buff;
	
	if(debugging("trace"))
	    printf("Entering setParallelVapor()\n");
        for (i = 0; i < dim; i++)
        {
            global_gmax[i] = pp_grid->Global_grid.gmax[i]-1;
            local_gmax[i] = pp_grid->Zoom_grid.gmax[i]-1;
        }
        FT_MakeGridIntfc(front);
        setDomain();
        G_size = 1;
        L_size = 1;
        for (i = 0; i < dim; i++)
        {
            G_size = G_size * (global_gmax[i]+1);
            L_size = L_size * (top_gmax[i]+1);
        }
        uni_array(&V_buff,L_size,sizeof(double));
        if (myid == 0)
        {
            uni_array(&GV_buff,G_size,sizeof(double));
	    /*setInitialVapor(front,LIQUID_COMP,GV_buff);*/
	    /*Please set GV_buff before sending out*/
            for (id = 0; id < numprocs; id++)
            {
                find_Cartesian_coordinates(id,pp_grid,pp_icoords);
		switch (dim)
		{
		    case 2:
                	for (j = jmin; j <= jmax; ++j)
                	for (i = imin; i <= imax; ++i)
                	{
                    	    icoords[0] = i;
                    	    icoords[1] = j;
                    	    G_icoords[0] = pp_icoords[0]*(local_gmax[0]+1)+icoords[0]-imin;
                    	    G_icoords[1] = pp_icoords[1]*(local_gmax[1]+1)+icoords[1]-jmin;
                    	    G_index = d_index(G_icoords,global_gmax,dim);
                    	    index = d_index(icoords,top_gmax,dim);
                    	    V_buff[index] = GV_buff[G_index];
                	}
			break;
		    case 3:
			for (k = kmin; k <= kmax; ++k)
                	for (j = jmin; j <= jmax; ++j)
                	for (i = imin; i <= imax; ++i)
                	{
                    	    icoords[0] = i;
                    	    icoords[1] = j;
                    	    icoords[2] = k;
                    	    G_icoords[0] = pp_icoords[0]*(local_gmax[0]+1)+icoords[0]-imin;
                    	    G_icoords[1] = pp_icoords[1]*(local_gmax[1]+1)+icoords[1]-jmin;
                    	    G_icoords[2] = pp_icoords[2]*(local_gmax[2]+1)+icoords[2]-kmin;
                    	    G_index = d_index(G_icoords,global_gmax,dim);
                    	    index = d_index(icoords,top_gmax,dim);
                    	    V_buff[index] = GV_buff[G_index];
                	}
			break;
            	    Default:
                  	printf("Unknown dim = %d\n",dim);
                	clean_up(ERROR);
		}
                if (id == 0)
                {
                    for (i = 0; i < L_size; i++)
                        field->vapor[i] = V_buff[i];
                }
                else
                {
                    pp_send(1,(POINTER)(V_buff),sizeof(double)*L_size,id);
                }
            }
            FT_FreeThese(1,GV_buff);
        }
        else
        {
            pp_recv(1,0,(POINTER)(V_buff),sizeof(double)*L_size);
            for (i = 0; i < L_size; i++)
            {
                field->vapor[i] = V_buff[i];
            }
        }

        FT_FreeThese(1,V_buff);
        setAdvectionDt();
	if(debugging("trace"))
	    printf("Leaving setParallelVapor()\n");
}

void VCARTESIAN::setIndexMap(COMPONENT sub_comp)
{
	static boolean first = YES;
	int i,j,k,ic,index;
	int llbuf[MAXD],uubuf[MAXD];
	int count;

	if (debugging("trace")) printf("Entering setIndexMap()\n");
	if (first)
	{
	    first = NO;
	    switch (dim)
	    {
	    case 1:
		count = (imax - imin + 1);
	    	FT_VectorMemoryAlloc((POINTER*)&i_to_I,top_gmax[0]+1,INT);
	    	FT_MatrixMemoryAlloc((POINTER*)&I_to_i,count,1,INT);
	    	break;
	    case 2:
		count = (imax - imin + 1)*(jmax - jmin + 1);
	    	FT_MatrixMemoryAlloc((POINTER*)&ij_to_I,top_gmax[0]+1,
					top_gmax[1]+1,INT);
	    	FT_MatrixMemoryAlloc((POINTER*)&I_to_ij,count,2,INT);
	    	break;
	    case 3:
		count = (imax - imin + 1)*(jmax - jmin + 1)*(kmax - kmin + 1);
	    	FT_TriArrayMemoryAlloc((POINTER*)&ijk_to_I,top_gmax[0]+1,
					top_gmax[1]+1,top_gmax[2]+1,INT);
	    	FT_MatrixMemoryAlloc((POINTER*)&I_to_ijk,count,3,INT);
	    	break;
	    }
	}

	index = 0;
	for (i = 0; i < dim; ++i)
	{
	    llbuf[i] = lbuf[i] != 0 ? lbuf[i] : 1;
	    uubuf[i] = ubuf[i] != 0 ? ubuf[i] : 1;
	}
	switch (dim)
	{
	case 1:
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index1d(i,top_gmax);
		if (sub_comp == NO_COMP || cell_center[ic].comp == sub_comp)
		{
	    	    i_to_I[i] = index + ilower;
	    	    index++;
		}
		else
		    i_to_I[i] = -1;
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)i_to_I);
	    break;
	case 2:
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (sub_comp == NO_COMP || cell_center[ic].comp == sub_comp)
		{
	    	    ij_to_I[i][j] = index + ilower;
		    I_to_ij[index][0] = i;
                    I_to_ij[index][1] = j;
	    	    index++;
		}
		else
		    ij_to_I[i][j] = -1;
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ij_to_I);
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (sub_comp == NO_COMP || cell_center[ic].comp == sub_comp)
		{
	    	    ijk_to_I[i][j][k] = index + ilower;
		    I_to_ijk[index][0] = i;
                    I_to_ijk[index][1] = j;
                    I_to_ijk[index][2] = k;
	    	    index++;
		}
		else
		    ijk_to_I[i][j][k] = -1;
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ijk_to_I);
	    break;
	}
	if (debugging("trace")) printf("Leaving setIndexMap()\n");

}	/* end setIndexMap */

void VCARTESIAN::computeAdvection()
{
	int i;
	COMPONENT sub_comp[2];

	sub_comp[0] = SOLID_COMP;
	sub_comp[1] = LIQUID_COMP2;

	for (i = 0; i < 2; ++i)
	{
		if(sub_comp[i] == SOLID_COMP)
		    continue;
	        setGlobalIndex(sub_comp[i]);
	        if (eqn_params->num_scheme == CRANK_NICOLSON)
		{
		    if(eqn_params->prob_type == PARTICLE_TRACKING)
	    		computeVaporSource();
		    else
	    		computeVolumeForce();

	    	    computeAdvectionCN(sub_comp[i],field->vapor,eqn_params->D);
		    
		    computeTemperatureSource();
	    	    computeAdvectionCN(sub_comp[i],field->temperature,eqn_params->D);
		}
	        else if (eqn_params->num_scheme == WENO_CRANK_NICOLSON)
		{
		    if(eqn_params->prob_type == PARTICLE_TRACKING)
	    		computeVaporSource();
		    else
	    		computeVolumeForce();
	    	    computeAdvectionWENO(sub_comp[i],field->vapor,eqn_params->D);

		    computeTemperatureSource();
	    	    computeAdvectionWENO(sub_comp[i],field->temperature,eqn_params->D);
		}
	}
}

void VCARTESIAN::computeAdvectionWENO(COMPONENT sub_comp,double* Temp,const double D)
{
	int i,j,k,l,m,ic,icn,I,I_nb,icoords[MAXD];
	int gmin[MAXD],ipn[MAXD];
	double crx_coords[MAXD];
	double T0,T_nb,lambda,coeff,coeff_nb,rhs;
	COMPONENT comp;
	PETSc solver;
	double *x;
        int num_iter = 0;
        double rel_residual = 0;
        boolean fr_crx_grid_seg;
        const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	double v[MAXD];

	start_clock("computeAdvectionWENO");
	if (debugging("trace")) printf("Entering computeAdvectionWENO()\n");
	//testParallel("adv",field->adv);
	//testParallel("vap",field->vapor);

	/*compute advection term with WENO5*/
	static WENO_FLUX *weno_flux = new WENO_FLUX(*front);
	weno_flux->adv = field->adv;
	weno_flux->U = field->vel;
	weno_flux->T = field->vapor;
	weno_flux->findStateAtCrossing = find_state_at_crossing;
	weno_flux->computeAdvection();

	/*extrapolate q at t = t_n+0.5*/
	/*to achieve 2-order*/
	double p0,p1;
	static double old_dt = 0.0;
	p0 = p1 = 0.0;
	if (old_dt != 0)
	{
	    p0 = 1.0 + m_dt/(2.0*old_dt);
	    p1 = -0.5*m_dt/old_dt;
	}
	old_dt = m_dt;
	/*end extrapolation*/
	
	for (i = 0; i < dim; ++i) gmin[i] = 0;

	setIndexMap(sub_comp);

	/*end computing advection term with WENO5*/
	start_clock("set_coefficients");
	switch(dim)
	{
	case 1:
	    solver.Create(ilower, iupper-1, 3, 3);
	    for (i = imin; i <= imax; ++i)
	    {
	    	icoords[0] = i;
	    	ic = d_index1d(i,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp != sub_comp) 
	    	    continue;
		I = i_to_I[i];
                T0 = Temp[ic];
		rhs = T0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*m_dt/sqr(top_h[l]);
		    coeff += lambda;
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index1d(ipn[0],top_gmax);
			I_nb = i_to_I[ipn[0]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,getStateVapor,
                                &T_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
			    T_nb = Temp[icn];
			    rhs -= coeff_nb*(T_nb - T0);
                        }
			else
			{
			    rhs -= coeff_nb*(2.0*T_nb - T0);
			}
                    }
		}
		solver.Add_A(I,I,coeff);
		solver.Add_b(I,rhs);
	    }
	    break;
	case 2:
	    solver.Create(ilower, iupper-1, 5, 5);
	    for (j = jmin; j <= jmax; ++j)
	    for (i = imin; i <= imax; ++i)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = top_comp[ic];
		I = ij_to_I[i][j];
	    	if (comp != sub_comp) 
	    	    continue;
                T0 = Temp[ic];
		rhs = T0;
		if (source != NULL)
		    rhs += m_dt*source[ic];
		/*advection term computed by WENO*/
		rhs += m_dt*(p0*field->adv[ic]
		    + p1*field->adv_old[ic]); 	
		coeff = 1.0;
	 	for (l = 0; l < dim; ++l) v[l] = 0.0;
                if (field->vel != NULL)
                {
                    for (l = 0; l < dim; ++l)
                        v[l] = field->vel[l][ic];
                }
                for (l = 0; l < dim; ++l)
                {
                    lambda = D*m_dt/sqr(top_h[l]);
		    coeff += lambda;
		    
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index2d(ipn[0],ipn[1],top_gmax);
                        I_nb = ij_to_I[ipn[0]][ipn[1]];
                        coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,getStateVapor,
                                &T_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                            solver.Add_A(I,I_nb,coeff_nb);
                            T_nb = Temp[icn];
                            rhs -= -0.5*lambda*(T_nb - T0);
                        }
                        else
                        {
                            rhs -= -0.5*lambda*(2.0*T_nb - T0);
                        }
                    }
                }
                solver.Add_A(I,I,coeff);
                solver.Add_b(I,rhs);
	    }
	    break;
	case 3:
	    solver.Create(ilower, iupper-1, 7, 7);
	    for (k = kmin; k <= kmax; ++k)
	    for (j = jmin; j <= jmax; ++j)
	    for (i = imin; i <= imax; ++i)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	icoords[2] = k;
	    	ic = d_index3d(i,j,k,top_gmax);
	    	comp = top_comp[ic];
		I = ijk_to_I[i][j][k];
	    	if (comp != sub_comp) 
	    	    continue;
                T0 = Temp[ic];
		rhs = T0;
                if (source != NULL)
                    rhs += m_dt*source[ic];
		/*advection term computed by WENO*/
		rhs += m_dt*(p0*field->adv[ic]
		    + p1*field->adv_old[ic]); 
		coeff = 1.0;
                for (l = 0; l < dim; ++l) v[l] = 0.0;
                if (field->vel != NULL)
                {
                    for (l = 0; l < dim; ++l)
                        v[l] = field->vel[l][ic];
                }
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*m_dt/sqr(top_h[l]);
		    coeff += lambda;
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
			I_nb = ijk_to_I[ipn[0]][ipn[1]][ipn[2]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,getStateVapor,
                                &T_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
                            T_nb = Temp[icn];
                            rhs -= -0.5*lambda*(T_nb - T0);
                        }
			else
			{
			    rhs -= coeff_nb*(2.0*T_nb - T0);
			}
                    }
		}
		solver.Add_A(I,I,coeff);
		solver.Add_b(I,rhs);
	    }
	    break;
	}
	stop_clock("set_coefficients");
	start_clock("petsc_solve");
	solver.SetMaxIter(10);   
	solver.SetTol(1e-20);   
	solver.Solve();
	solver.GetNumIterations(&num_iter);
        solver.GetFinalRelativeResidualNorm(&rel_residual);

	if (debugging("PETSc"))
	{
	    (void) printf("VCARTESIAN::computeAdvectionWENO: "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	}

	FT_VectorMemoryAlloc((POINTER*)&x,cell_center.size(),sizeof(double));
	solver.Get_x(x);
	stop_clock("petsc_solve");

	start_clock("scatter_data");
	switch (dim)
	{
        case 1:
            for (i = imin; i <= imax; i++)
            {
	    	I = i_to_I[i];
	    	ic = d_index1d(i,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
	    	    array[ic] = x[I-ilower];
	    	else
	    	    array[ic] = 0.0;
	    }
	    break;
        case 2:
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
	    	I = ij_to_I[i][j];
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
	    	    array[ic] = x[I-ilower];
	    	else
	    	    array[ic] = 0.0;
	    }
	    break;
        case 3:
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
	    	I = ijk_to_I[i][j][k];
	    	ic = d_index3d(i,j,k,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
	    	    array[ic] = x[I-ilower];
	    	else
	    	    array[ic] = 0.0;
	    }
	    break;
	}
	scatMeshArray();
	switch (dim)
	{
        case 1:
            for (i = 0; i <= top_gmax[0]; ++i)
            {
		ic = d_index1d(i,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
		    Temp[ic] = array[ic];
	    }
	    break;
        case 2:
            for (j = 0; j <= top_gmax[1]; ++j)
            for (i = 0; i <= top_gmax[0]; ++i)
            {
		ic = d_index2d(i,j,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
		    Temp[ic] = array[ic];
		field->adv_old[ic] = field->adv[ic];
	    }
	    break;
        case 3:
            for (k = 0; k <= top_gmax[2]; ++k)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (i = 0; i <= top_gmax[0]; ++i)
            {
		ic = d_index3d(i,j,k,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
		    Temp[ic] = array[ic];
		field->adv_old[ic] = field->adv[ic];
	    }
	    break;
	}
	stop_clock("scatter_data");
	FT_FreeThese(1,x);

	if (debugging("trace")) printf("Leaving computeAdvectionWENO()\n");
	stop_clock("computeAdvectionWENO");

}	/* end computeAdvectionCN */

void VCARTESIAN::testParallel(const char* ptag, double* Temp)
{

	FILE* outfile;
        char filename[200];
	int i,j,k,ic;
        sprintf(filename,"%s/%s-%d-nd%d",
        OutName(front),ptag,front->step,pp_mynode());
        outfile = fopen(filename,"w");
	if (pp_numnodes() > 1)
	{
        for (j = 0; j <= top_gmax[1]; ++j)
        for (i = 0; i <= top_gmax[0]; ++i)
        {
                ic = d_index2d(i,j,top_gmax);
                fprintf(outfile,"%20.19f\n",Temp[ic]);
        }
        fclose(outfile);
	}
	else
	{
        sprintf(filename,"%s/%s-%d-nd%d",
        OutName(front),ptag,front->step,0);
        outfile = fopen(filename,"w");
	for (j = 0; j <= top_gmax[1]; ++j)
        for (i = 0; i <= (imin+imax)/2+4; ++i)
        {
                ic = d_index2d(i,j,top_gmax);
                fprintf(outfile,"%20.19f\n",Temp[ic]);
        }
        fclose(outfile);

        sprintf(filename,"%s/%s-%d-nd%d",
        OutName(front),ptag,front->step,1);
        outfile = fopen(filename,"w");
	for (j = 0; j <= top_gmax[1]; ++j)
        for (i = (imin+imax)/2 - 3; i <= top_gmax[0]; ++i)
        {
                ic = d_index2d(i,j,top_gmax);
                fprintf(outfile,"%20.19f\n",Temp[ic]);
        }
        fclose(outfile);
	}
}

void VCARTESIAN::computeAdvectionCN(COMPONENT sub_comp,double* Temp,const double D)
{
	int i,j,k,l,m,ic,icn,I,I_nb,icoords[MAXD];
	int gmin[MAXD],ipn[MAXD];
	double crx_coords[MAXD];
	double T0,T_nb,lambda,coeff,coeff_nb,rhs;
	COMPONENT comp;
	PETSc solver;
	double *x;
        int num_iter = 0;
        double rel_residual = 0;
        boolean fr_crx_grid_seg;
        const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	double v[MAXD];
	double eta;

	start_clock("computeAdvectionCN");
	if (debugging("trace")) printf("Entering computeAdvectionCN()\n");

	for (i = 0; i < dim; ++i) gmin[i] = 0;

	setIndexMap(sub_comp);

	start_clock("set_coefficients");
	switch(dim)
	{
	case 1:
	    solver.Create(ilower, iupper-1, 3, 3);
	    for (i = imin; i <= imax; ++i)
	    {
	    	icoords[0] = i;
	    	ic = d_index1d(i,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp != sub_comp) 
	    	    continue;
		I = i_to_I[i];
                T0 = Temp[ic];
		rhs = T0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*m_dt/sqr(top_h[l]);
		    coeff += lambda;
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index1d(ipn[0],top_gmax);
			I_nb = i_to_I[ipn[0]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,getStateVapor,
                                &T_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
			    T_nb = Temp[icn];
			    rhs -= coeff_nb*(T_nb - T0);
                        }
			else
			{
			    rhs -= coeff_nb*(2.0*T_nb - T0);
			}
                    }
		}
		solver.Add_A(I,I,coeff);
		solver.Add_b(I,rhs);
	    }
	    break;
	case 2:
	    solver.Create(ilower, iupper-1, 5, 5);
	    for (j = jmin; j <= jmax; ++j)
	    for (i = imin; i <= imax; ++i)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = top_comp[ic];
		I = ij_to_I[i][j];
	    	if (comp != sub_comp) 
	    	    continue;
                T0 = Temp[ic];
		rhs = T0;
		if (source != NULL)
		    rhs += m_dt*source[ic];
		coeff = 1.0;
	 	for (l = 0; l < dim; ++l) v[l] = 0.0;
                if (field->vel != NULL)
                {
                    for (l = 0; l < dim; ++l)
                        v[l] = field->vel[l][ic];
                }
                for (l = 0; l < dim; ++l)
                {
                    lambda = D*m_dt/sqr(top_h[l]);
                    eta = v[l]*m_dt/top_h[l];
		    coeff += lambda;
		    
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index2d(ipn[0],ipn[1],top_gmax);
                        I_nb = ij_to_I[ipn[0]][ipn[1]];
                        coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,getStateVapor,
                                &T_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                            solver.Add_A(I,I_nb,coeff_nb);
                            T_nb = Temp[icn];
                            rhs -= -0.5*lambda*(T_nb - T0);
			    if(v[l] > 0 && m == 0)
				rhs += eta * (T_nb - T0);
			    if(v[l] < 0 && m == 1)
				rhs += eta * (T0 - T_nb);
                        }
                        else
                        {
                            rhs -= -0.5*lambda*(2.0*T_nb - T0);
                        }
                    }
                }
                solver.Add_A(I,I,coeff);
                solver.Add_b(I,rhs);
	    }
	    break;
	case 3:
	    solver.Create(ilower, iupper-1, 7, 7);
	    for (k = kmin; k <= kmax; ++k)
	    for (j = jmin; j <= jmax; ++j)
	    for (i = imin; i <= imax; ++i)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	icoords[2] = k;
	    	ic = d_index3d(i,j,k,top_gmax);
	    	comp = top_comp[ic];
		I = ijk_to_I[i][j][k];
	    	if (comp != sub_comp) 
	    	    continue;
                T0 = Temp[ic];
		rhs = T0;
                if (source != NULL)
                    rhs += m_dt*source[ic];
		coeff = 1.0;
                for (l = 0; l < dim; ++l) v[l] = 0.0;
                if (field->vel != NULL)
                {
                    for (l = 0; l < dim; ++l)
                        v[l] = field->vel[l][ic];
                }
		for (l = 0; l < dim; ++l)
		{
                    eta = v[l]*m_dt/top_h[l];
		    lambda = D*m_dt/sqr(top_h[l]);
		    coeff += lambda;
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
			I_nb = ijk_to_I[ipn[0]][ipn[1]][ipn[2]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,getStateVapor,
                                &T_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
                            T_nb = Temp[icn];
                            rhs -= -0.5*lambda*(T_nb - T0);
                            if(v[l] > 0 && m == 0)
                                rhs += eta * (T_nb - T0);
                            if(v[l] < 0 && m == 1)
                                rhs += eta * (T0 - T_nb);
                        }
			else
			{
			    rhs -= coeff_nb*(2.0*T_nb - T0);
			}
                    }
		}
		solver.Add_A(I,I,coeff);
		solver.Add_b(I,rhs);
	    }
	    break;
	}
	stop_clock("set_coefficients");
	start_clock("petsc_solve");
	solver.SetMaxIter(500);   
	solver.SetTol(1e-8);   
	solver.Solve();

	if (debugging("PETSc"))
	{
	    (void) printf("VCARTESIAN::computeAdvectionCN: "
	       		"num_iter = %d, rel_residual = %g \n", 
			num_iter, rel_residual);
	}

	FT_VectorMemoryAlloc((POINTER*)&x,cell_center.size(),sizeof(double));
	solver.Get_x(x);
	stop_clock("petsc_solve");

	start_clock("scatter_data");
	switch (dim)
	{
        case 1:
            for (i = imin; i <= imax; i++)
            {
	    	I = i_to_I[i];
	    	ic = d_index1d(i,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
	    	    array[ic] = x[I-ilower];
	    	else
	    	    array[ic] = 0.0;
	    }
	    break;
        case 2:
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
	    	I = ij_to_I[i][j];
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
	    	    array[ic] = x[I-ilower];
	    	else
	    	    array[ic] = 0.0;
	    }
	    break;
        case 3:
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
	    	I = ijk_to_I[i][j][k];
	    	ic = d_index3d(i,j,k,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
	    	    array[ic] = x[I-ilower];
	    	else
	    	    array[ic] = 0.0;
	    }
	    break;
	}
	scatMeshArray();
	switch (dim)
	{
        case 1:
            for (i = 0; i <= top_gmax[0]; ++i)
            {
		ic = d_index1d(i,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
		    Temp[ic] = array[ic];
	    }
	    break;
        case 2:
            for (j = 0; j <= top_gmax[1]; ++j)
            for (i = 0; i <= top_gmax[0]; ++i)
            {
		ic = d_index2d(i,j,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
		    Temp[ic] = array[ic];
	    }
	    break;
        case 3:
            for (k = 0; k <= top_gmax[2]; ++k)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (i = 0; i <= top_gmax[0]; ++i)
            {
		ic = d_index3d(i,j,k,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
		    Temp[ic] = array[ic];
	    }
	    break;
	}
	stop_clock("scatter_data");
	FT_FreeThese(1,x);

	if (debugging("trace")) printf("Leaving computeAdvectionCN()\n");
	stop_clock("computeAdvectionCN");
}	/* end computeAdvectionCN */

// for initial condition: 
// 		setInitialCondition();	
// this function should be called before solve()
// for the source term of the momentum equation: 	
// 		computeSourceTerm();
void VCARTESIAN::solve(double dt)
{
	if (debugging("trace")) printf("Entering solve()\n");
	m_dt = dt;
	start_clock("solve");

	setDomain();
        if (debugging("trace")) printf("Passing setDomain()\n");

	setComponent();
	if (debugging("trace")) printf("Passing setComponent()\n");
	
	computeSource();
	if (debugging("trace")) printf("Passing computeSource()\n");
	
	computeAdvection();
	if (debugging("trace")) printf("Passing computeAdvection()\n");
	
	computeSupersat();
	if (debugging("trace")) printf("Passing computeSupersat()\n");

	setAdvectionDt();
	if (debugging("trace")) printf("Passing setAdvectionDt()\n");
	stop_clock("solve");

	if (debugging("trace")) printf("Leaving solve()\n");
}

void VCARTESIAN::computeSupersat()
{	
	int i,j,k,size,index;
	double *temp  = eqn_params->field->temperature;
	double *vapor  = eqn_params->field->vapor;
	double sat_vap_pre, sat_vap_rat;
	double *super = eqn_params->field->supersat;

	/*compute coeffecient for condensation*/
        double Lh, Rv, Rd, rhoL, Kc, es, D, Cp, ksi;
        double T,p;
        double A;
	/*see instructions in climate.h*/
        Lh = eqn_params->Lh;
        Rv = eqn_params->Rv;
        Rd = eqn_params->Rd;
        rhoL = eqn_params->rho_l;
        Kc = eqn_params->Kc;
        D = eqn_params->D;
        ksi = Rd/Rv;

        T = 0;
        p = 0;
	for (index = 0; index < comp_size; index++)
        {
            T += field->temperature[index];
            p += field->pres[index];
        }
	size = comp_size;
#if defined(__MPI)
	pp_gsync();
        pp_global_sum(&T,1);
        pp_global_sum(&p,1);
        pp_global_isum(&size,1);
#endif
        T /= size;
        p /= size;

        es = 611.2*exp(17.67*(T-273.15)/(T-29.65));
        eqn_params->K = 1/((Lh/(Rv*T)-1)*Lh*rhoL/(Kc*T)+rhoL*Rv*T/(D*es));
	printf("Condensation coeffecient = %e\n",eqn_params->K);

	switch(dim)
	{
	    case 2:
		for (i = 0; i <= top_gmax[0]; ++i)
            	for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index = d_index2d(i,j,top_gmax);
		    sat_vap_pre = 611.2*exp(17.67*(temp[index]-273.15)
					           /(temp[index]-29.65));
		    sat_vap_rat = 621.97 * sat_vap_pre
					/(eqn_params->field->pres[index]
					  -sat_vap_pre);
                    if(top_comp[index] == SOLID_COMP)
                        super[index] = 0;
                    else
                        super[index] = vapor[index]/sat_vap_rat - 1;
		    if(eqn_params->field->pres[index] == 0)
			super[index] = 0;
		}	
		break;
	    case 3:
                for (i = 0; i <= top_gmax[0]; ++i)
                for (j = 0; j <= top_gmax[1]; ++j)
		for (k = 0; k <= top_gmax[2]; ++k)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    sat_vap_pre = 611.2*exp(17.67*(temp[index]-273.15)
						/(temp[index]-29.65));
                    sat_vap_rat = 621.97 * sat_vap_pre
					/(eqn_params->field->pres[index]
					  -sat_vap_pre);
		    if(top_comp[index] == SOLID_COMP)
			super[index] = 0;
		    else
                        super[index] = vapor[index]/sat_vap_rat - 1;
		    if(eqn_params->field->pres[index] == 0)
			super[index] = 0;
                }
		break;
	}
}


void VCARTESIAN::setAdvectionDt()
{
	double D,Dl,Ds;
	static double m_dt_expl,m_dt_impl;  // explicit and implicit time step
	static boolean first = YES;

	if (first)
	{
	    first = NO;
	    eqn_params = (PARAMS*)front->extra2;
	    D = eqn_params->D;
	    m_dt_expl = 0.5*sqr(hmin)/D/(double)dim;
	    m_dt_impl = 0.5*hmin/D/(double)dim;
	    min_dt = 0.1*sqr(hmin)/D/(double)dim;
	}

	m_dt = m_dt_expl;
	if (debugging("trace"))
	{
	    printf("In setAdvectionDt: m_dt = %24.18g min_dt = %f\n",
				m_dt,min_dt);
	}
}	/* end setAdvectionDt */

void VCARTESIAN::getVelocity(double *p, double *U)
{
	// locate the point
	int icoords[MAXD];
	int i,j,k;
	double c1[MAXD], c2[MAXD], denominator;

	if (!rect_in_which(p,icoords,top_grid))
	{
	    for (i=0; i<2; i++)
	    {
	    	U[i] = 0.0;
	    }
	    return;
	}

	switch (dim)
	{
	case 2:
	    break;
	}
}

void VCARTESIAN::getRectangleIndex(int index, int &i, int &j)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
}

void VCARTESIAN::getRectangleIndex(int index, int &i, int &j, int &k)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
	k = cell_center[index].icoords[2];
}


int VCARTESIAN::getRectangleComponent(int index)
{	
	return getComponent(cell_center[index].icoords);
}

void VCARTESIAN::getRectangleCenter(
	int index, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	    coords[i] = cell_center[index].coords[i];
}

void VCARTESIAN::getRectangleCenter(
	int index0, 
	int index1, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = 0.5*(cell_center[index0].coords[i] +
	    		     cell_center[index1].coords[i]);
	}
}

int VCARTESIAN::getComponent(
	double *p)
{
	return component(p,front->interf);
}

int VCARTESIAN::getComponent(
	int *icoords)
{
	int index;
	switch (dim)
	{
	case 1:
	    index = d_index1d(icoords[0],top_gmax);
	    return top_comp[index];
	case 2:
	    index = d_index2d(icoords[0],icoords[1],top_gmax);
	    return top_comp[index];
	case 3:
	    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
	    return top_comp[index];
	}
}

void VCARTESIAN::save(char *filename)
{
	
	RECT_GRID *rect_grid = front->rect_grid;
	INTERFACE *intfc    = front->interf;
		
	int i, j;
	int xmax = rect_grid->gmax[0];
	int ymax = rect_grid->gmax[1];
	double x, y;
	
	FILE *hfile = fopen(filename, "w");
	if(hfile==NULL)
	{
		printf("\n can't open %s in "
		       "SaveAsTecplot_rect_grid_and_interface().", filename);
		exit(0);
	}
	
	// secondly print out the interface
		
	if(exists_interface(intfc))
	{
	    CURVE		**cur;
	    CURVE		*curve;
	    BOND		*bond;
			
	    for(cur=intfc->curves; cur && *cur; cur++)
	    {
		curve = *cur;
		fprintf(hfile, "ZONE I=%d J=%d F=POINT \n", 
				curve->num_points, 1);
		bond=curve->first;
		fprintf(hfile, "%.4f %.4f \n",bond->start->_coords[0], 
				bond->start->_coords[1]);
		for(bond=curve->first; bond!=NULL; bond=bond->next)
		    fprintf(hfile, "%.4f %.4f \n",bond->end->_coords[0], 
		    		bond->end->_coords[1]);
		}					
	}		
	fclose(hfile);
}

VCARTESIAN::VCARTESIAN(Front &front):front(&front),field(NULL)
{
}

void VCARTESIAN::makeGridIntfc()
{
	static boolean first = YES;
	INTERFACE *grid_intfc;
	Table *T;
	static double *vapor;
	int i;

	FT_MakeGridIntfc(front);
	if (debugging("trace")) printf("Passed FT_MakeGridIntfc()\n");

	grid_intfc = front->grid_intfc;
	top_grid = &topological_grid(grid_intfc);
	lbuf = front->rect_grid->lbuf;
	ubuf = front->rect_grid->ubuf;
	top_gmax = top_grid->gmax;
	top_L = top_grid->L;
	top_U = top_grid->U;
	top_h = top_grid->h;
	dim = grid_intfc->dim;
	T = table_of_interface(grid_intfc);
	top_comp = T->components;
	eqn_params = (PARAMS*)front->extra2;
	hmin = top_h[0];
	for (i = 1; i < dim; ++i)
	    if (hmin > top_h[i]) hmin = top_h[i];

	switch (dim)
	{
	case 1:
            if (first)
            {
                FT_VectorMemoryAlloc((POINTER*)&array,top_gmax[0]+1,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&vapor,top_gmax[0]+1,FLOAT);
                first = NO;
            }	
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    eqn_params->field->vapor = vapor;
	    eqn_params = (PARAMS*)front->extra2;
	    break;
	case 2:
	    if (first)
	    {
	    	FT_VectorMemoryAlloc((POINTER*)&array,(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&vapor,(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    eqn_params->field->vapor = vapor;
	    break;
	case 3:
	    if (first)
	    {
	    	FT_VectorMemoryAlloc((POINTER*)&array,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&vapor,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
	    eqn_params->field->vapor = vapor;
	    break;
	}
}	/* end makeGridIntfc */

void VCARTESIAN::deleteGridIntfc()
{
	FT_FreeGridIntfc(front);
}

void VCARTESIAN::scatMeshArray()
{
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
}

void VCARTESIAN::setGlobalIndex(COMPONENT sub_comp)
{
	int i,j,k,ic;
	int num_nodes = pp_numnodes();
	int myid = pp_mynode();
	static boolean first = YES;

	if (first)
	{
	    first = NO;
	    FT_VectorMemoryAlloc((POINTER*)&n_dist,num_nodes,sizeof(int));
	}
	NLblocks = 0;
	switch (dim)
	{
	case 1:
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index1d(i,top_gmax);
		if (sub_comp != NO_COMP && cell_center[ic].comp != sub_comp) 
		    continue;
		NLblocks++;
	    }
	    break;
	case 2:
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (sub_comp != NO_COMP && cell_center[ic].comp != sub_comp) 
		    continue;
		NLblocks++;
	    }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (sub_comp != NO_COMP && cell_center[ic].comp != sub_comp) 
		    continue;
		NLblocks++;
	    }
	    break;
	}

	for (i = 0; i < num_nodes; ++i) n_dist[i] = 0;
	n_dist[myid] = NLblocks;
	pp_global_imax(n_dist,num_nodes);
	ilower = 0;
        iupper = n_dist[0];

        for (i = 1; i <= myid; i++)
        {
            ilower += n_dist[i-1];
            iupper += n_dist[i];
        }	
}


void VCARTESIAN::printFrontInteriorState(char *out_name)
{
        int i,j,k,l,index;
        char filename[100];
        FILE *outfile;
        INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        double *Temp = field->vapor;
        double *Supersat = field->supersat;

        sprintf(filename,"%s/state.ts%s-vapor",out_name,
                        right_flush(front->step,7));
#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */

        outfile = fopen(filename,"w");

        /* Initialize states at the interface */
        fprintf(outfile,"Interface states:\n");
        int count = 0;
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            count++;
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fprintf(outfile,"%24.18g %24.18g\n",sl->vapor,
                                sr->vapor);
            fprintf(outfile,"%24.18g %24.18g\n",sl->supersat,
                                sr->supersat);
        }
        fprintf(outfile,"\nInterior states:\n");
        switch (dim)
        {
        case 1:
            for (i = 0; i <= top_gmax[0]; ++i)
            {
                index = d_index1d(i,top_gmax);
                fprintf(outfile,"%24.18g\n",Temp[index]);
                fprintf(outfile,"%24.18g\n",Supersat[index]);
            }
            break;
        case 2:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            {
                index = d_index2d(i,j,top_gmax);
                fprintf(outfile,"%24.18g\n",Temp[index]);
                fprintf(outfile,"%24.18g\n",Supersat[index]);
            }
            break;
        case 3:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (k = 0; k <= top_gmax[2]; ++k)
            {
                index = d_index3d(i,j,k,top_gmax);
                fprintf(outfile,"%24.18g\n",Temp[index]);
                fprintf(outfile,"%24.18g\n",Supersat[index]);
            }
            break;
        }
        fclose(outfile);
}

void VCARTESIAN::readFrontInteriorState(char *restart_name)
{
        FILE *infile;
        int i,j,k,l,index;
        INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        double x;
        double *Temp = field->vapor;
        double *Supersat = field->supersat;

        char fname[100];
        sprintf(fname,"%s-vapor",restart_name);
        infile = fopen(fname,"r");

        /* Initialize states at the interface */
        next_output_line_containing_string(infile,"Interface states:");
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fscanf(infile,"%lf",&x);
            sl->vapor = x;
            fscanf(infile,"%lf",&x);
            sr->vapor = x;
            fscanf(infile,"%lf",&x);
            sl->supersat = x;
            fscanf(infile,"%lf",&x);
            sr->supersat = x;
        }

        FT_MakeGridIntfc(front);
        setDomain();
        /* Initialize states in the interior regions */

        next_output_line_containing_string(infile,"Interior states:");

        switch (dim)
        {
        case 1:
            for (i = 0; i <= top_gmax[0]; ++i)
            {
                index = d_index1d(i,top_gmax);
                fscanf(infile,"%lf",&Temp[index]);
                fscanf(infile,"%lf",&Supersat[index]);
            }
            break;
        case 2:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            {
                index = d_index2d(i,j,top_gmax);
                fscanf(infile,"%lf",&Temp[index]);
                fscanf(infile,"%lf",&Supersat[index]);
            }
            break;
        case 3:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (k = 0; k <= top_gmax[2]; ++k)
            {
                index = d_index3d(i,j,k,top_gmax);
                fscanf(infile,"%lf",&Temp[index]);
                fscanf(infile,"%lf",&Supersat[index]);
            }
            break;
        }
        fclose(infile);
}

void VCARTESIAN::setBoundary()
{
	int i,j,k,index0,index1;
	INTERFACE *intfc = front->interf;
	double *Temp = field->vapor;

	switch (dim)
	{
	case 1:
	    if (rect_boundary_type(intfc,0,0) == NEUMANN_BOUNDARY)
	    {
		index0 = d_index1d(0,top_gmax);
		index1 = d_index1d(1,top_gmax);
	    	Temp[index0]  = Temp[index1];
	    }
	    if (rect_boundary_type(intfc,0,1) == NEUMANN_BOUNDARY)
	    {
		index0 = d_index1d(top_gmax[0],top_gmax);
		index1 = d_index1d(top_gmax[0]-1,top_gmax);
	    	Temp[index0]  = Temp[index1];
	    }
	    break;
	case 2:
	    if (rect_boundary_type(intfc,0,0) == NEUMANN_BOUNDARY)
	    {
		for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index0 = d_index2d(0,j,top_gmax);
		    index1 = d_index2d(1,j,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    if (rect_boundary_type(intfc,0,1) == NEUMANN_BOUNDARY)
	    {
		for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index0 = d_index2d(top_gmax[0],j,top_gmax);
		    index1 = d_index2d(top_gmax[0]-1,j,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    if (rect_boundary_type(intfc,1,0) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		{
		    index0 = d_index2d(i,0,top_gmax);
		    index1 = d_index2d(i,1,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    if (rect_boundary_type(intfc,1,1) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		{
		    index0 = d_index2d(i,top_gmax[1],top_gmax);
		    index1 = d_index2d(i,top_gmax[1]-1,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    break;
	case 3:
	    if (rect_boundary_type(intfc,0,0) == NEUMANN_BOUNDARY)
	    {
		for (j = 0; j <= top_gmax[1]; ++j)
		for (k = 0; k <= top_gmax[2]; ++k)
		{
		    index0 = d_index3d(0,j,k,top_gmax);
		    index1 = d_index3d(1,j,k,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    if (rect_boundary_type(intfc,0,1) == NEUMANN_BOUNDARY)
	    {
		for (j = 0; j <= top_gmax[1]; ++j)
		for (k = 0; k <= top_gmax[2]; ++k)
		{
		    index0 = d_index3d(top_gmax[0],j,k,top_gmax);
		    index1 = d_index3d(top_gmax[0]-1,j,k,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    if (rect_boundary_type(intfc,1,0) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		for (k = 0; k <= top_gmax[2]; ++k)
		{
		    index0 = d_index3d(i,0,k,top_gmax);
		    index1 = d_index3d(i,1,k,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    if (rect_boundary_type(intfc,1,1) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		for (k = 0; k <= top_gmax[2]; ++k)
		{
		    index0 = d_index3d(i,top_gmax[1],k,top_gmax);
		    index1 = d_index3d(i,top_gmax[1]-1,k,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    if (rect_boundary_type(intfc,2,0) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index0 = d_index3d(i,j,0,top_gmax);
		    index1 = d_index3d(i,j,1,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    if (rect_boundary_type(intfc,2,1) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index0 = d_index3d(i,j,top_gmax[2],top_gmax);
		    index1 = d_index3d(i,j,top_gmax[2]-1,top_gmax);
	    	    Temp[index0]  = Temp[index1];
		}
	    }
	    break;
	}
}	/* end setBoundary */

void VCARTESIAN::vtk_plot3d(
        const char *varname,double *var)
{
	std::vector<int> ph_index;
	char *outname = front->out_name;
        int i,j,k,index;
        char dirname[256],filename[256];
        FILE *outfile;
        double coord_x,coord_y,coord_z,xmin,ymin,zmin;
        COMPONENT comp;
        int pointsx,pointsy,pointsz,num_points,num_cells,num_cell_list;
        int icoords[3],p_gmax[3];

        int ii,jj,kk;
        double ih,jh,kh;

        sprintf(dirname,"%s/vtk",OutName(front));
        if (pp_mynode() == 0)
        {
            if (!create_directory(dirname,YES))
            {
                screen("Cannot create directory %s\n",dirname);
                clean_up(ERROR);
            }
        }
        pp_gsync();
        sprintf(dirname,"%s/vtk.ts%s",dirname,right_flush(front->step,7));
        if (pp_numnodes() > 1)
            sprintf(dirname,"%s-nd%s",dirname,right_flush(pp_mynode(),4));
        if (!create_directory(dirname,YES))
        {
            screen("Cannot create directory %s\n",dirname);
            clean_up(ERROR);
        }

        //cell-based liquid phase
        ph_index.clear();

        sprintf(filename,"%s/%s.vtk",dirname,varname);
        outfile = fopen(filename,"w");
        fprintf(outfile,"# vtk DataFile Version 3.0\n");
        fprintf(outfile,"%s\n",varname);
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            if (cell_center[index].comp == LIQUID_COMP2)
                ph_index.push_back(index);
        }

        pointsx = top_gmax[0] + 2;
        pointsy = top_gmax[1] + 2;
        pointsz = top_gmax[2] + 2;
        num_points = pointsx*pointsy*pointsz;

        num_cells = (int)ph_index.size();
        num_cell_list = 9*num_cells;

        fprintf(outfile,"POINTS %d double\n", num_points);

        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

	for (i = 0; i <= top_gmax[0]; i++)
        for (j = 0; j <= top_gmax[1]; j++)
        {
            k = top_gmax[2];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] + top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (i = 0; i <= top_gmax[0]; i++)
        for (k = 0; k <= top_gmax[2]; k++)
        {
            j = top_gmax[1];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] + top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (j = 0; j <= top_gmax[1]; j++)
        for (k = 0; k <= top_gmax[2]; k++)
        {
            i = top_gmax[0];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] + top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

	for (i = 0; i <= top_gmax[0]; i++)
        {
            j = top_gmax[1];
            k = top_gmax[2];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] + top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] + top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (j = 0; j <= top_gmax[1]; j++)
        {
            i = top_gmax[0];
            k = top_gmax[2];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] + top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] + top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (k = 0; k <= top_gmax[2]; k++)
        {
            i = top_gmax[0];
            j = top_gmax[1];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].coords[0] + top_h[0]/2.0;
            coord_y = cell_center[index].coords[1] + top_h[1]/2.0;
            coord_z = cell_center[index].coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

	i = top_gmax[0];
        j = top_gmax[1];
        k = top_gmax[2];
        index = d_index3d(i,j,k,top_gmax);
        coord_x = cell_center[index].coords[0] + top_h[0]/2.0;
        coord_y = cell_center[index].coords[1] + top_h[1]/2.0;
        coord_z = cell_center[index].coords[2] + top_h[2]/2.0;
        fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);


        fprintf(outfile,"CELLS %i %i\n", num_cells,num_cell_list);
        for (i = 0; i < num_cells; i++)
        {
            int index0,index1,index2,index3,index4,index5,index6,index7;
            index = ph_index[i];
            icoords[0] = cell_center[index].icoords[0];
            icoords[1] = cell_center[index].icoords[1];
            icoords[2] = cell_center[index].icoords[2];
            index0 = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
            index1 =
                d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
            index2 =
                d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
            index3 =
                d_index3d(icoords[0]+1,icoords[1]+1,icoords[2],top_gmax);
            index4 =
                d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);
            index5 =
                d_index3d(icoords[0]+1,icoords[1],icoords[2]+1,top_gmax);
            index6 =
                d_index3d(icoords[0],icoords[1]+1,icoords[2]+1,top_gmax);
            index7 =
                d_index3d(icoords[0]+1,icoords[1]+1,icoords[2]+1,top_gmax);

            fprintf(outfile,"8 %i %i %i %i %i %i %i %i\n",
                index0,index1,index2,index3,index4,index5,index6,index7);
        }

	fprintf(outfile, "CELL_TYPES %i\n", num_cells);
        for (i = 0; i < num_cells; i++)
            fprintf(outfile,"11\n");

        fprintf(outfile, "CELL_DATA %i\n", num_cells);
        fprintf(outfile, "SCALARS %s double\n",varname);
        fprintf(outfile, "LOOKUP_TABLE default\n");
	if(var != NULL)
	{
	    for (i = 0; i < num_cells; i++)
            {
                index = ph_index[i];
                fprintf(outfile,"%f\n",var[index]);
            }
	}
	else
	{
	 	printf("The var is NULL!\n");
		clean_up(ERROR);
	}
        fclose(outfile);
}       /* end vtk_plot_vapor3d */

void VCARTESIAN::recordParticles()
{
	char fname[256];
	if (debugging("trace"))
            printf("Entering record particles\n");
        sprintf(fname,"%s/record_particles",OutName(front));
        if (!create_directory(fname,NO))
        {
            printf("Cannot create directory %s\n",fname);
            clean_up(ERROR);
        }
        sprintf(fname,"%s/particles-%4.2f",fname,front->time);
	recordParticles(fname,eqn_params->particle_array,eqn_params->num_drops);
}

#include<string>
static double (*getStateVel[MAXD])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};
void VCARTESIAN::recordParticles(char* filename, PARTICLE* particle_array, int num_particles)
{
	/*this function records all the information associated with particles*/
	/*using MPIIO to write data to one single file*/
	const int num_cols = 12;
	const int num_bytes = num_cols*sizeof(double); /*number of bytes per particle*/  
	int num_nodes = pp_numnodes();
	int my_node = pp_mynode();
	int total_particles = num_particles;
	long *offset = new long[num_nodes]();

	pp_global_isum(&total_particles,1);
	/*calculate offsets for each processor*/
	for (int i = my_node+1; i < num_nodes; ++i)
	    offset[i] = num_particles*num_bytes;
	pp_global_lsum(offset,num_nodes);
	
	/*Open MPI files*/
	MPI_Offset my_offset;
        MPI_File fh;
        MPI_Status status;
	bool file_open_error = true;
	file_open_error = MPI_File_open(MPI_COMM_WORLD, filename,
                     		MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	if (file_open_error != MPI_SUCCESS)
	{
	    printf("node %d cannot open file %s\n",my_node,filename);
	    clean_up(ERROR);
	}

	if (debugging("MPIIO"))	
	{
		printf("In recordParticles()\n");
		printf("%d particles\n",num_particles);
		printf("my_id = %d, offset = %ld bytes = %ld rows\n",
			my_node,offset[my_node],offset[my_node]/(num_cols*sizeof(double)));
		printf("file name = %s\n",filename);
	}
	/*writing header*/
	const int header_offset = 256*sizeof(char);
	char cstr[200];
	std::string str;
	if (my_node == 0)
	{
            MPI_File_seek(fh, 0, MPI_SEEK_SET);
	    str = "#Data type: double, offset: 256 chars\n";
	    MPI_File_write(fh,(void*)str.c_str(),str.size(),MPI_CHAR,&status); 
	    sprintf(cstr,"#%d rows: particles\n",total_particles);
	    str = cstr;
	    MPI_File_write(fh,(void*)str.c_str(),str.size(),MPI_CHAR,&status);
	    sprintf(cstr,"#%d columns: radius(1), coords(3), vel(3), supersaturation(1), temperature(1), vel(3)\n",
		    num_cols);
            str = cstr;
            MPI_File_write(fh,(void*)str.c_str(),str.size(),MPI_CHAR,&status);
	}
	/*writing data*/
	my_offset = header_offset+offset[my_node];
        MPI_File_seek(fh, my_offset, MPI_SEEK_SET);
	for (int i = 0; i < num_particles; ++i)
	{
	    double buffer[num_cols] = {0};
	    buffer[0] = particle_array[i].radius;
	    for (int j = 0; j < FT_Dimension(); ++j)
	    {
		buffer[j+1] = particle_array[i].center[j];
		buffer[j+4] = particle_array[i].vel[j];
		FT_IntrpStateVarAtCoords(front,LIQUID_COMP2,particle_array[i].center,
					field->vel[j],getStateVel[j],buffer+9+j,NULL);
	    }
	    FT_IntrpStateVarAtCoords(front,LIQUID_COMP2,particle_array[i].center,
					field->supersat,NULL,buffer+7,NULL);
	    FT_IntrpStateVarAtCoords(front,LIQUID_COMP2,particle_array[i].center,
					field->temperature,NULL,buffer+8,NULL);
	    MPI_File_write(fh, buffer, num_cols, MPI_DOUBLE, &status);
	}
	MPI_File_close(&fh);
}

void VCARTESIAN::recordGroupParticles()
{
	int i, j, index;
        int ic[MAXD];
        PARAMS* eqn_params = (PARAMS*)front->extra2;
        PARTICLE* particle_array = eqn_params->particle_array;
        double *coords,block_size[MAXD];
        static boolean first = YES;
        FILE *file;
        char fname[100];

        int group_id, num_drops = eqn_params->num_drops;
        int group_crds[MAXD], group_gmax[MAXD];
        static double *group_mean[2]; /*mean value for different groups*/
        static double *group_init[2];
        static int NGroup = 1;
        /*group_mean[0]:radius^3*/
        /*group_mean[1]:total number of droplets in group*/

        if (debugging("trace"))
            printf("Entering record clip \n");
        if (first)
        {
            NGroup  = 1;
            for (i = 0; i < dim; i++)
                NGroup *= eqn_params->Nclip[i];
            for (i = 0; i < 2; i++)
            {
                group_mean[i] = new double[NGroup];
                group_init[i] = new double[NGroup];
            }
        }
        for (i = 0; i < dim; i++)
        {
            block_size[i] = (front->pp_grid->Global_grid.U[i]-
                             front->pp_grid->Global_grid.L[i])/
                             eqn_params->Nclip[i];
            group_gmax[i] = eqn_params->Nclip[i] - 1;
        }

        for (j = 0; j < 2; j++)
        for (i = 0; i < NGroup; i++)
            group_mean[j][i] = 0.0;

        for (i = 0; i < num_drops; i++)
        {
            coords = particle_array[i].center;
	    /*exclude particles out of sample box*/
            for (j = 0; j < dim; j++)
	    ic[j] = floor((coords[j] - top_L[j] + 0.5*top_h[j])/top_h[j]);
            index = d_index(ic,top_gmax,dim);
	    if (field->sample[index] != 1)
		continue;
	    /*end of excluding*/
            for (j = 0; j < dim; j++)
                group_crds[j] = floor(coords[j]/block_size[j]);
            group_id = d_index(group_crds,group_gmax,dim);
            if (particle_array[i].radius > 0.0)
            {
                group_mean[0][group_id] += pow(particle_array[i].radius,3);
                group_mean[1][group_id] += 1.0;
            }
        }
        for (j = 0; j < 2; j++)
            pp_global_sum(group_mean[j],NGroup);

        for (i = 0; i < NGroup; i++)
            group_mean[0][i] /= group_mean[1][i];

        if (first)
        {
            first = NO;
            for (j = 0; j < 2; j++)
            for (i = 0; i < NGroup; i++)
                 group_init[j][i] = group_mean[j][i];
        }

        if (pp_mynode() == 0)
        {
            sprintf(fname,"%s/record-group",OutName(front));
            if (!create_directory(fname,NO))
            {
                printf("Cannot create directory %s\n",fname);
                clean_up(ERROR);
            }
            sprintf(fname,"%s/group-%f",fname,front->time);
            file = fopen(fname,"w");
            for (i = 0; i < NGroup; i++)
                fprintf(file,"%f %f %f\n",
                        front->time,
                        group_mean[0][i]/group_init[0][i],
                        group_mean[1][i]/group_init[1][i]);
            fclose(file);
        }
}

void VCARTESIAN::recordClusteringIndex()
{
	double Mean, Var, CL, sigma, CSL;
	static boolean first = YES;
	FILE *file;
	char fname[100];
	char *out_name = front->out_name;
	double *array = field->drops;
	int size = comp_size;
	sprintf(fname,"%s/cluster",out_name);
	if (first)
	{
	    if(pp_mynode() == 0)
	    {
	        file = fopen(fname,"w");
	        fprintf(file,"%%Clustering Index\n");
	        fprintf(file,"%%CL	sigma	CSL = CL/sigma\n");
	        fclose(file);
	    }
	    first = NO;
	}
	Deviation(array,size,Mean,Var);
	CL = Var/Mean - 1.0;
	sigma = sqrt(2.0/(size));
	CSL = CL/sigma;
	if(pp_mynode() == 0)
	{
	    file = fopen(fname,"a");
	    fprintf(file,"%f    %e    %f\n",CL,sigma,CSL);
	    fclose(file);
	}
	return;
}

double VCARTESIAN::computeReactTimeScale(double R0, double S0, int N)
{
	PARAMS *eqn_params = (PARAMS*)front->extra2;
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	
	int i,j,k,size,index;
	double R2, S, t, dt = 0.001;
	double Lh, Rv, Rd, rhoL, K, es, D, Cp, ksi;
	double T, p;  
	double A, B;
	Lh = eqn_params->Lh;
	Rv = eqn_params->Rv;
	Rd = eqn_params->Rd;
	rhoL = eqn_params->rho_l;
	K = eqn_params->Kc;
	D = eqn_params->D;
	Cp = eqn_params->Cp;
	ksi = Rd/Rv;

	T = 0;
	p = 0;
	size = 0;
        for (k = kmin; k < kmax; k++)
        for (j = jmin; j < jmax; j++)
        for (i = imin; i < imax; i++)
        {
	    size ++;
	    index = d_index3d(i,j,k,top_gmax);
	    T += field->temperature[index];
	    p += field->pres[index];
	}
#if defined(__MPI)
	pp_gsync();
	pp_global_sum(&T,1);
	pp_global_sum(&p,1);
	pp_global_isum(&size,1);
#endif
	T /= size;
	p /= size;

	es = 611.2*exp(17.67*(T-273.15)/(T-29.65));
	A = 1/((Lh/(Rv*T)-1)*Lh*rhoL/(K*T)+rhoL*Rv*T/(D*es));
	B = 4*PI*N*rhoL*(Rd*T/(ksi*es)+ksi*sqr(Lh)/(p*T*Cp))
	    /((Lh/(Rv*T)-1)*Lh*rhoL/(K*T)+rhoL*Rv*T/(D*es));
	t = 0;
	R2 = sqr(R0);
	S = S0;
	printf("A = %e, B = %e\n",A,B);
	while(R2 > 0 && S < -0.005)
	{
	    R2 += 2*A*S*dt;
	    S  /= 1+B*dt*sqrt(R2);
	    t  += dt;
	}
	return t;
}

double VCARTESIAN::computeDspRateLinear()
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        double Sij[MAXD][MAXD];
        double **vel = field->vel;
        double DspRat = 0.0;
        double ReduceBuff[2];
        double nu = iFparams->mu2/iFparams->rho2;
        int i,j,k,l,m,size; 
	int I_nb[MAXD][2],ipn[MAXD],gmin[MAXD],icoords[MAXD];
	int I0;
        const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
        for (i = 0; i < dim; i++)
            gmin[i] = 0;

        size = 1;
        switch(dim)
        {
	    case 2:
		size = (imax-imin+1)*(jmax-jmin+1);
		for (i = imin; i <= imax; i++)
		for (j = jmin; j <= jmax; j++)
		{
		    icoords[0] = i;
            	    icoords[1] = j;
		    I0 = d_index2d(i,j,top_gmax);
            	    for (m = 0; m < dim; ++m)
            	    for (l = 0; l < 2; ++l)
             	    {
                	if(next_ip_in_dir(icoords,dir[m][l],ipn,gmin,top_gmax))
                        I_nb[m][l] = d_index2d(ipn[0],ipn[1],top_gmax);
                    	else
                	{
                    	    printf("In computeDspRate(), cannot find next ip\n");
                    	    printf("gmin=[%d %d]\n",gmin[0],gmin[1]);
                    	    printf("gmax=[%d %d]\n",top_gmax[0],top_gmax[1]);
                    	    printf("icoords=[%d %d]\n",icoords[0],icoords[1]);
                    	    clean_up(ERROR);
                	}
            	    }
		    for (l = 0; l < dim; l++)
		    for (m = 0; m < dim; m++)
		    {
			DspRat += vel[l][I0]
				*(vel[l][I_nb[m][1]]+vel[l][I_nb[m][0]]
				- 2*vel[l][I0])
				/(sqr(top_h[m]));
		    }
		}
		break;
	    case 3:
		size = (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1);
		for (i = imin; i <= imax; i++)
		for (j = jmin; j <= jmax; j++)
		for (k = kmin; k <= kmax; k++)
		{
		    icoords[0] = i;
            	    icoords[1] = j;
		    icoords[2] = k;
		    I0 = d_index3d(i,j,k,top_gmax);
            	    for (m = 0; m < dim; ++m)
            	    for (l = 0; l < 2; ++l)
             	    {
                	if(next_ip_in_dir(icoords,dir[m][l],ipn,gmin,top_gmax))
                            I_nb[m][l] = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
                    	else
                	{
                    	    printf("In computeDspRate(), cannot find next ip\n");
                    	    printf("gmin=[%d %d]\n",gmin[0],gmin[1]);
                    	    printf("gmax=[%d %d]\n",top_gmax[0],top_gmax[1]);
                    	    printf("icoords=[%d %d]\n",icoords[0],icoords[1]);
                    	    clean_up(ERROR);
                	}
            	    }
		    for (l = 0; l < dim; l++)
		    for (m = 0; m < dim; m++)
		    {
			DspRat += vel[l][I0]
				*( vel[l][I_nb[m][1]]
				 + vel[l][I_nb[m][0]]
				 -2*vel[l][I0])
				/(sqr(top_h[m]));
		    }
		}

		break;
	    default:
		printf("In computeDspRateLinear(), unknown dim = %d\n",dim);
	        clean_up(ERROR);
	}
	pp_gsync();
	pp_global_isum(&size,1);
	pp_global_sum(&DspRat,1);
	DspRat = -nu*DspRat/size;
	return DspRat;
}

double VCARTESIAN::computeDspRate()
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	double Sij[MAXD][MAXD];
	double **vel = field->vel;
	double DspRat = 0.0;
	double ReduceBuff[2];
	double mu = iFparams->mu2/iFparams->rho2;
	int i,j,k,l,m,I_nb[MAXD][2],ipn[MAXD],gmin[MAXD],icoords[MAXD],size;
	const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	for (i = 0; i < dim; i++)
	{
            gmin[i] = 0;
	}

	size = 0;
	switch(dim)
	{
	
	case 2:
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    size ++;
	    icoords[0] = i;
            icoords[1] = j;
            for (m = 0; m < dim; ++m)
	    for (l = 0; l < 2; ++l)
            {
		if(next_ip_in_dir(icoords,dir[m][l],ipn,gmin,top_gmax))
                    I_nb[m][l] = d_index2d(ipn[0],ipn[1],top_gmax);
		else
		{
		    printf("In computeDspRate(), cannot find next ip\n");
	            printf("gmin=[%d %d]\n",gmin[0],gmin[1]);
		    printf("gmax=[%d %d]\n",top_gmax[0],top_gmax[1]);
		    printf("icoords=[%d %d]\n",icoords[0],icoords[1]);
		    clean_up(ERROR);
		}
            }
	    /*compute strain rate tensor Sij*/
	    for(l = 0; l < dim; ++l)
	    for(m = 0; m < dim; ++m)
	    {
		Sij[l][m]  = 0.5*(vel[l][I_nb[m][1]]-vel[l][I_nb[m][0]])/(2.0*top_h[m]);
		Sij[l][m] += 0.5*(vel[m][I_nb[l][1]]-vel[m][I_nb[l][0]])/(2.0*top_h[l]);
	    }
	    /*compute dissipation rate using: 2.0*mu*Sij^2*/ 
            for(l = 0; l < dim; ++l)
            for(m = 0; m < dim; ++m)
	    {
		DspRat += 2.0*mu*sqr(Sij[l][m]);
	    }
	}
	break;
	case 3:
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    size ++;
	    icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            for (m = 0; m < dim; ++m)
	    for (l = 0; l < 2; ++l)
            {
		if(next_ip_in_dir(icoords,dir[m][l],ipn,gmin,top_gmax))
                    I_nb[m][l] = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
		else
		{
		    printf("In computeDspRate(), cannot find next ip\n");
	            printf("gmin=[%d %d %d]\n",gmin[0],gmin[1],gmin[2]);
		    printf("gmax=[%d %d %d]\n",top_gmax[0],top_gmax[1],top_gmax[2]);
		    printf("icoords=[%d %d %d]\n",icoords[0],icoords[1],icoords[2]);
		    clean_up(ERROR);
		}
            }
	    /*compute strain rate tensor Sij*/
	    for(l = 0; l < dim; ++l)
	    for(m = 0; m < dim; ++m)
	    {
		Sij[l][m]  = 0.5*(vel[l][I_nb[m][1]]-vel[l][I_nb[m][0]])/(2.0*top_h[m]);
		Sij[l][m] += 0.5*(vel[m][I_nb[l][1]]-vel[m][I_nb[l][0]])/(2.0*top_h[l]);
	    }
	    /*compute dissipation rate using: 2.0*mu*Sij^2*/ 
            for(l = 0; l < dim; ++l)
            for(m = 0; m < dim; ++m)
	    {
		DspRat += 2.0*mu*sqr(Sij[l][m]);
	    }
	}

	break;
	}
	ReduceBuff[0] = DspRat;
	ReduceBuff[1] = size;
	pp_gsync();
        pp_global_sum(ReduceBuff,2);
	DspRat = ReduceBuff[0];
	size = ReduceBuff[1];
	DspRat = DspRat/size; 
	return DspRat;
}

void VCARTESIAN::recordMixingLine()
{
	PARAMS *eqn_params = (PARAMS*)front->extra2;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	PARTICLE *particle_array = eqn_params->particle_array;
	static boolean first = YES;
	double mu = iFparams->mu2, eta;
	double **vel = field->vel;
	double t_evap, t_mix, t_react;
	double Dev; /*Deviation of radius distribution*/
	double alpha; /**/
	double rv,r0,rm; /*dissipation rate*/
	double NL; /*transition scale number*/
	double env_supersat;  /*supersaturation in clear air*/
	int index, i, j, k, l, m, m1, m2, ic, I;
	int nzeros, nzeros_in_samplebox;
	int gmin[MAXD], ipn[MAXD], icoords[MAXD];
	const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	FILE *file;
	char fname[100];
	char *out_name = front->out_name;
	int  myid = pp_mynode();
	double L = front->pp_grid->Global_grid.U[0] 
		 - front->pp_grid->Global_grid.L[0];
	double ReduceBuff[4];

	static int max_array_size = 0;
	static double* radius_array;

	if (debugging("trace"))
	    printf("Entering record mixing line\n");
        if (eqn_params->num_drops > max_array_size)
        {
            max_array_size = eqn_params->num_drops;
            free_these(1,radius_array);
            FT_VectorMemoryAlloc((POINTER*)&radius_array,max_array_size,FLOAT);
        }

	double a3 = 1.0;
	for (i = 0; i < dim; i++) 
	{
	    gmin[i] = 0;
	    a3 *= top_h[i]; 
	}
	/*compute mean energy dissipation rate*/
	eqn_params->disp_rate = computeDspRate();
	
	/*compute Rv, t_evap*/
	rv = 0;
	r0 = 0;
	nzeros = 0;
	nzeros_in_samplebox = 0;
	for (i = 0; i < eqn_params->num_drops; i++)
	{   
	    radius_array[i] = particle_array[i].radius;
	    if (particle_array[i].radius > 0)
		nzeros ++;
	    for (j = 0; j < dim; j++)
                icoords[j] = floor((particle_array[i].center[j] 
			- top_L[j] + 0.5*top_h[j])/top_h[j]);
            index = d_index(icoords,top_gmax,dim);
	    if (particle_array[i].radius > 0 &&
		field->sample[index] == 1)
		nzeros_in_samplebox ++;
	    rv += pow(particle_array[i].radius,3.0);
	    r0 += pow(particle_array[i].R0,3.0);
	}

#if defined(__MPI__)
	ReduceBuff[0] = rv;
	ReduceBuff[1] = r0;
	ReduceBuff[2] = nzeros;
	ReduceBuff[3] = nzeros_in_samplebox;
	pp_gsync();
	pp_global_sum(ReduceBuff,4);
	rv = ReduceBuff[0];
	r0 = ReduceBuff[1];
	nzeros    = (int)ReduceBuff[2];
	nzeros_in_samplebox = (int)ReduceBuff[3];
#endif 
	rv /= nzeros;
	r0 /= nzeros;
	/*compute t_mix*/
	t_mix = pow(L*L/eqn_params->disp_rate,1.0/3.0);

	/*compute mean value and deviation of radius distribution*/
	Deviation(radius_array,eqn_params->num_drops,rm,Dev);

	env_supersat = eqn_params->qe/eqn_params->qs - 1.0;
	if (env_supersat < 0)
	{
	    t_evap  = sqr(rm)/(-eqn_params->K*env_supersat);
	    t_react = computeReactTimeScale(rm,env_supersat,nzeros); 
	}
	else
	{
	    t_evap = HUGE;
	    t_react = 0;
	}
	if (myid != 0) /*Post analysis is only done by the master processor*/
	    return; /*return if processor id is slaver processor*/

	if (first)
	{
	    sprintf(fname,"%s/mixing",out_name);
	    file = fopen(fname,"w");
	    fprintf(file,"%%t_mix    t_evap  t_react   Damkoehler\n");
	    fclose(file);

            sprintf(fname,"%s/RN",out_name);
            file = fopen(fname,"w");
            fprintf(file,"%%time    Rv    Rm    stdDev   N_total    N_in_samplebox\n");
            fclose(file);

            sprintf(fname,"%s/transition",out_name);
            file = fopen(fname,"w");
            fprintf(file,"%%time    epsilon    eta    t_react    NL\n");
            fclose(file);
	    first = NO;
	}
	sprintf(fname,"%s/mixing",out_name);
	file = fopen(fname,"a");
	fprintf(file,"%15.14f  %15.14f  %15.14f  %15.14f\n",
		t_mix,t_evap,t_react,t_mix/t_evap);
	fclose(file);
	/*plot mixing line: N/N0 VS. (Rv/Rv0)^3*/
        sprintf(fname,"%s/RN",out_name);
        file = fopen(fname,"a");
	fprintf(file,"%f  %15.14f  %15.14f  %15.14f  %d  %d\n", 
		front->time,pow(rv,1.0/3.0),rm,sqrt(Dev),nzeros,nzeros_in_samplebox);
	fclose(file);
        sprintf(fname,"%s/transition",out_name);
        file = fopen(fname,"a");
	eta = pow(pow(mu,3.0)/eqn_params->disp_rate,0.25);
	NL = pow(eqn_params->disp_rate,0.5)*pow(t_react,1.5)/eta;
	fprintf(file,"%f  %15.14f  %15.14f  %15.14f  %15.14f\n",
		front->time,eqn_params->disp_rate,eta,t_react,NL);
	fclose(file);
	if (debugging("trace"))
	    printf("Leaving record mixing line\n");
}

void VCARTESIAN::checkField()
{
	int i,j,k,index,count;
	printf("Enter checkField()\n");
	double prev = 0;

	count = 0;
	for(j = 0; j < comp_size; j++)
        {
                if(field->temperature[j] != prev)
		{
                    printf("%s[%d] = %20.14f\n",
                        "temperature",j,
                        field->temperature[j]);
		    prev = field->temperature[j];
		    count ++;
		}
		if(count > 20)
		{
		    printf("......\n");
		    break;
		}
        }

	count = 0;
        for(j = 0; j < comp_size; j++)
        {
		if(field->vapor[j] != prev)
                { 
		    printf("%s[%d] = %20.14f\n",
                        "vapor",j,
                        field->vapor[j]);
		    prev = field->vapor[j];
		    count ++;
		}
		if(count > 20)
                {
                    printf("......\n");
                    break;
                }
        }

	count = 0;
        for(j = 0; j < comp_size; j++)
        {
		if(field->supersat[j] != prev)
                {
		    printf("%s[%d] = %20.14f\n",
                        "supersat",j,
                        field->supersat[j]);
		    prev = field->supersat[j];
		    count++;
		}
		if(count > 20)
                {
                    printf("......\n");
                    break;
                }
        }
	
	count = 0;
        for(j = 0; j < comp_size; j++)
        {
		if(field->drops[j] != prev)
                {
		    printf("%s[%d] = %20.14f\n",
                        "drops",j,
                        field->drops[j]);
		    prev = field->drops[j];
		    count++;
		}
		if(count > 20)
                {
                    printf("......\n");
                    break;
                }
        }

	int s_count = 0, l_count = 0, o_count = 0;
	if (dim == 2)
	{
	    for(j = jmin; j <= jmax; j++)
            {
                int i = (int)(imax/2);
		index = d_index2d(i,j,top_gmax);
                    printf("%s[%d][%d] = (%f,%f)\n",
                        "velo",i,j,
                        field->vel[0][index],field->vel[1][index]);
            }
            for(i = imin; i <= imax; i++)
            {
		j = (int)(jmax/2);
		index = d_index2d(i,j,top_gmax);
                    printf("%s[%d][%d] = %f\n",
                        "pres",i,j,
                        field->pres[index]);
            }
	    index = d_index2d(imin,jmin,top_gmax);
	    printf("pres[%d][%d] = %f\n",imin,jmin,field->pres[index]);

            for (j = jmin; j <= jmax; ++j)
            for (i = imin; i <= imax; ++i)
            {
                index = d_index2d(i,j,top_gmax);
	        if(top_comp[index] == SOLID_COMP)
	    	    s_count ++; 
	        if(top_comp[index] == LIQUID_COMP2)
		    l_count ++;
	        if(field->vapor[index] == 0)
		    o_count ++;
	    }
	}
	else if (dim == 3)
	{
	    for(k = kmin; k <= kmax; ++k)
            {
                i = (int)(imax/2);
		j = (int)(jmax/2); 
		index = d_index3d(i,j,k,top_gmax);
                    printf("%s[%d][%d][%d] = (%f,%f,%f)\n",
                        "velo",i,j,k,
                        field->vel[0][index],field->vel[1][index],field->vel[2][index]);
            }
            for(i = imin; i <= imax; i++)
            {
		j = (int)(jmax/2);
		k = (int)(kmax/2);
		index = d_index3d(i,j,k,top_gmax);
                    printf("%s[%d][%d][%d] = %f\n",
                        "pres",i,j,k,
                        field->pres[index]);
            }
	    index = d_index3d(imin,jmin,kmin,top_gmax);
	    printf("pres[%d][%d][%d] = %f\n",imin,jmin,kmin,field->pres[index]);


	    for (k = kmin; k <= kmax; ++k)
            for (j = jmin; j <= jmax; ++j)
            for (i = imin; i <= imax; ++i)
            {
                index = d_index3d(i,j,k,top_gmax);
                if(top_comp[index] == SOLID_COMP)
                    s_count ++;
                if(top_comp[index] == LIQUID_COMP2)
                    l_count ++;
                if(field->vapor[index] == 0)
                    o_count ++;
            }
	}
}


void VCARTESIAN::recordField(char *outname, const char *varname)
{
	int i, j, k, index;
	FILE* outfile;
	char filename[256];
	double **vel = field->vel;
	
        sprintf(filename, "%s/record-%s",outname,varname);
        if (!create_directory(filename,NO))
        {
            printf("Cannot create directory %s\n",filename);
            clean_up(ERROR);
        }
	pp_gsync();

	if (pp_numnodes() > 1)
	{
	    sprintf(filename, "%s/record-%s-nd%03d",filename,varname,pp_mynode());
	    create_directory(filename,YES);
	}
        sprintf(filename,"%s/%s-%4.2f",filename,varname,front->time);
        outfile = fopen(filename,"w");
        switch (dim)
        {
            case 2:
                for (j = jmin; j <= jmax; j++)
                for (i = imin; i <= imax; i++)
                {
                    index = d_index2d(i,j,top_gmax);
		    if(strcmp(varname,"vapor") == 0)
                        fprintf(outfile,"%f\n",field->vapor[index]);
		    else if(strcmp(varname,"supersat") == 0)
			fprintf(outfile,"%f\n",field->supersat[index]);
		    else if(strcmp(varname,"temperature") == 0)
			fprintf(outfile,"%f\n",field->temperature[index]);
		    else if(strcmp(varname,"mrad") == 0)
			fprintf(outfile,"%15.14f\n",field->mrad[index]);
		    else if(strcmp(varname,"velocity") == 0)
			fprintf(outfile,"%f  %f\n",vel[0][index],vel[1][index]);
		    else
			printf("WARNING: Unknown field: %s\n",varname);
			
                }
                break;

                case 3:
                for (k = kmin; k <= kmax; k++)
                for (j = jmin; j <= jmax; j++)
                for (i = imin; i <= imax; i++)
                {
		    index = d_index3d(i,j,k,top_gmax);
                    if(strcmp(varname,"vapor") == 0)
                        fprintf(outfile,"%f\n",field->vapor[index]);
                    else if(strcmp(varname,"supersat") == 0)
                        fprintf(outfile,"%f\n",field->supersat[index]);
                    else if(strcmp(varname,"temperature") == 0)
                        fprintf(outfile,"%f\n",field->temperature[index]);
                    else if(strcmp(varname,"mrad") == 0)
                        fprintf(outfile,"%15.14f\n",field->mrad[index]);
                    else if(strcmp(varname,"velocity") == 0)
                        fprintf(outfile,"%f  %f  %f\n",
                            vel[0][index],vel[1][index],vel[2][index]);
                    else
                        printf("WARNING: Unknown field: %s\n",varname);
                }
                break;
        }
	fclose(outfile);
	return;
}

void VCARTESIAN::recordCondensationRate(char* outname){
	int i,index;
	double Cd_mean = 0.0;
	double *coords;
	PARTICLE* particle_array = eqn_params->particle_array;
	PARAMS* eqn_params = (PARAMS*)front->extra2;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	double rho_0 = iFparams->rho2;
	int num_drops = eqn_params->num_drops;
	double R, s; /*supersaturation at single droplet location*/
	double* L = front->pp_grid->Global_grid.L;
	double* U = front->pp_grid->Global_grid.U;
	double domain_volume = 1.0;
	int ic[MAXD];
	/*output file*/
	char filename[256];
	FILE* ofile;
	static boolean first = YES;

	for (i = 0; i < dim; i++)
	    domain_volume *= (U[i]-L[i]);

	for (i = 0; i < num_drops; i++)
	{
		coords = particle_array[i].center;
		R = particle_array[i].radius;
		rect_in_which(coords,ic,top_grid);
		index = d_index(ic,top_gmax,dim);
		s = field->supersat[index];
		FT_IntrpStateVarAtCoords(front,LIQUID_COMP,coords,
                                field->supersat,getStateSuper,&s,&s);
		Cd_mean += R * s;
	}
	pp_gsync();
	pp_global_sum(&Cd_mean,1);	
	Cd_mean *= 4.0*PI*eqn_params->rho_l*eqn_params->K
		   /(rho_0*domain_volume);

	/*output Cd_mean*/
	if (pp_mynode() == 0)
	{
	    sprintf(filename, "%s/Cd",outname);
	    if (first)
	    {
	        ofile = fopen(filename,"w");
	        fprintf(ofile,"%%time\tCondensationRate\n");
	        fclose(ofile);
	        first = NO;
	    }
	    ofile = fopen(filename,"a");
	    fprintf(ofile,"%f  %13.14e\n",front->time,Cd_mean);
	    fclose(ofile);
	}
}

void VCARTESIAN::recordPDF(char *outname, const char *varname)
{
	int i,j,k,index;
	FILE *outfile;
	FILE *superfile;
	char supername[256];
	char filename[256];
	double *PDF;
	double bin_size,mid_bin;
	double var_min,var_max;
	double mean_super,var_super,mean_temp,var_temp;
	double **vel = field->vel;
	int num_bins = 200;

	if(strcmp(varname,"all") == 0)
	{
	    recordPDF(outname,"vapor");
	    recordPDF(outname,"supersat");
	    recordPDF(outname,"temperature");
	    recordPDF(outname,"velocity");
	    recordPDF(outname,"numdensity");
	    recordPDF(outname,"cloud");
	    recordLagrangSupersat(outname);
	    return;
	}

	if (strcmp(varname,"velocity") == 0)
	{
	    recordPDF(outname,"xvel");
	    recordPDF(outname,"yvel");
	    if (dim == 3)
		recordPDF(outname,"zvel");	
	    return;
	}
	if (pp_mynode() == 0)
	{
            sprintf(filename, "%s/PDF-%s",outname,varname);

            if (!create_directory(filename,NO))
            {
                printf("Cannot create directory %s\n",filename);
                clean_up(ERROR);
            }
            sprintf(filename,"%s/%s-%4.2f",filename,varname,front->time);
            outfile = fopen(filename,"w");
	}
        if(strcmp(varname,"temperature") == 0)
	{
	    PDF = ComputePDF(field->temperature,comp_size,bin_size,num_bins,var_min,var_max);
	    Deviation(field->temperature,front,mean_temp,var_temp);
	    if (pp_mynode() == 0)
	    {
	        sprintf(supername, "%s/temperature",outname);
	        superfile = fopen(supername,"a");
		fprintf(superfile,"%f  %15.14f  %15.14f\n",front->time,mean_temp,var_temp);
		fclose(superfile);
	    }

	}
        else if(strcmp(varname,"vapor") == 0)
	{
	    PDF = ComputePDF(field->vapor,comp_size,bin_size,num_bins,var_min,var_max);
	    Deviation(field->vapor,front,mean_temp,var_temp);
            if (pp_mynode() == 0)
            {
                sprintf(supername, "%s/vapor",outname);
                superfile = fopen(supername,"a");
                fprintf(superfile,"%f  %15.14f  %15.14f\n",front->time,mean_temp,var_temp);
                fclose(superfile);
            }
	}
        else if(strcmp(varname,"supersat") == 0)
	{
	    PDF = ComputePDF(field->supersat,comp_size,bin_size,num_bins,var_min,var_max);
	    Deviation(field->supersat,front,mean_super,var_super);
	    if (pp_mynode() == 0)
	    {
	        sprintf(supername, "%s/supersat",outname);
	        superfile = fopen(supername,"a");
		fprintf(superfile,"%f  %15.14f  %15.14f\n",front->time,mean_super,var_super);
		fclose(superfile);
	    }
	}
	else if (strcmp(varname,"numdensity") == 0)
	{
	    recordClusteringIndex();
	    PDF = ComputePDF(field->drops,comp_size,bin_size,num_bins,var_min,var_max);
	}
	 else if(strcmp(varname,"cloud") == 0)
	{
	    PDF = ComputePDF(field->cloud,comp_size,bin_size,num_bins,var_min,var_max);
	}
        else if(strcmp(varname,"xvel") == 0)
	{
	    PDF = ComputePDF(field->vel[0],comp_size,bin_size,num_bins,var_min,var_max);
	}
        else if(strcmp(varname,"yvel") == 0)
	{
	    PDF = ComputePDF(field->vel[1],comp_size,bin_size,num_bins,var_min,var_max);
	}
        else if(strcmp(varname,"zvel") == 0)
	{
	    PDF = ComputePDF(field->vel[2],comp_size,bin_size,num_bins,var_min,var_max);
	}
	else
	{
	   printf("WARNING: Unknown field: %s\n",varname);
	   return; 
	}
	if (pp_mynode() == 0)
	{
	    for (i = 0; i < num_bins;i++)
	    {
	        mid_bin = var_min + (0.5+i)*bin_size;
	        fprintf(outfile,"%f  %f\n",mid_bin,PDF[i]);
	    }
	    fclose(outfile);
	}
}

void VCARTESIAN::recordTKE()
{
        double **vel = field->vel;
        int i,j,k,count,index,l;
        static boolean first = YES;
	double E, ReduceBuff[2];
        char fname[256];
        FILE *file;

        if (pp_mynode() == 0)
	{
            sprintf(fname,"%s/TKE",front->out_name);
	    if (first)
	    {
            	file = fopen(fname,"w");
                fprintf(file,"%%turbulence kinetic energy vs. time\n");
                fprintf(file,"%%time  TKE\n");
	    }
	    else
                file = fopen(fname,"a");
	}
	
	E = 0.0;
	count = 0;
	switch(dim)
	{
	    case 2:
                for (j = jmin; j <= jmax; ++j)
                for (i = imin; i <= imax; ++i)
                {
                    index = d_index2d(i,j,top_gmax);
		    for (l = 0.0; l < dim; l++)
		      E += 0.5*vel[l][index]*vel[l][index];
		    count++;	
		}
		break;
	    case 3:
                for (k = kmin; k <= kmax; ++k)
                for (j = jmin; j <= jmax; ++j)
		for (i = imin; i <= imax; ++i)
                {
                    index = d_index3d(i,j,k,top_gmax);
		    E += 0.5*vel[0][index]*vel[0][index];
                    E += 0.5*vel[1][index]*vel[1][index];
		    E += 0.5*vel[2][index]*vel[2][index];
		    count++;
		}
		break;
	    default:
		printf("Unknown dim = %d\n",dim);
                clean_up(ERROR);
	}
	/*MPI communications*/
	pp_gsync();
	ReduceBuff[0] = E;
	ReduceBuff[1] = count; 
	pp_global_sum(ReduceBuff,2);
	E = ReduceBuff[0];
	count = ReduceBuff[1];
	E = E/count;
        if (pp_mynode() == 0)
	{
	    fprintf(file,"%f  %20.19f\n",front->time,E);
	    fclose(file);
	}
	first = NO;
}


void VCARTESIAN::recordWaterBalance()
{
        double *vapor = field->vapor;
	double vapor_mass,liquid_mass,alpha;
	double rho,r,a3 = 1.0;
        int i,j,k,index,nzeros;
        static boolean first = YES;
	static double liquid_mass0 = 0.0;
	double water_balance[4];
	IF_PARAMS* iFparams = (IF_PARAMS*)front->extra1;
        PARAMS* eqn_params = (PARAMS*)front->extra2;
        PARTICLE* particle_array = eqn_params->particle_array;
        int num_drops = eqn_params->num_drops;
	double temp;
        char fname[256];
        FILE *file;
	double ReduceBuff[7];

	pp_gsync();
	if (pp_mynode() == 0)
	{
            sprintf(fname,"%s/water_balance",front->out_name);
            file = fopen(fname,"a");
	}

	for(i = 0; i < dim; i++)
            a3 *= top_h[i];

        if (first == YES && pp_mynode() == 0)
        {
	    fprintf(file,"%%water balance vs. time\n");
	    fprintf(file,"%%time  vapor  liquid  total  alpha\n");
        }
	/*compute water mass in vapor*/
        if (dim == 2)
        {
	    vapor_mass = 0.0;
	    for(j = jmin; j <= jmax; j++)
	    for(i = imin; i <= imax; i++)
	    {
                index = d_index2d(i,j,top_gmax);
		rho = iFparams->rho2;
		vapor_mass += rho * a3 *vapor[index]; 
	    }
        }
        else if (dim == 3)
        {
	    vapor_mass = 0.0;
	    int count = 0;
	    for(k = kmin; k <= kmax; k++)
            for(j = jmin; j <= jmax; j++)
            for(i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
		rho = iFparams->rho2;
                vapor_mass += rho * a3 *vapor[index];
            }
        }
	water_balance[0] = front->time;
        water_balance[1] = 0.001 * vapor_mass;
	/*compute water mass in droplets*/
	liquid_mass = 0.0;
	nzeros = 0;
	if (eqn_params->prob_type == PARTICLE_TRACKING)
	{
            for (i = 0; i < num_drops; i++)
            {
	        r = particle_array[i].radius;
		if (0 != r)
		    nzeros ++;
	        liquid_mass += 4.0/3.0*PI*r*r*r*particle_array[i].rho;
	    }
	    if (first) 
	    {
		liquid_mass0 = liquid_mass;
		first = NO; /*IMPORTANT: preserve initial LWC*/
		if (pp_mynode() == 0)
		    fclose(file);
		return;
	    }
	} 
	water_balance[2] = liquid_mass;
	water_balance[3] = liquid_mass + 0.001 * vapor_mass;

#if defined(__MPI__)
	for (i = 0; i < 3; i++)	
	    ReduceBuff[i] = water_balance[i+1];
	ReduceBuff[3] = nzeros;
	ReduceBuff[4] = num_drops;
	ReduceBuff[5] = liquid_mass;
	ReduceBuff[6] = liquid_mass0;
	pp_gsync();
	pp_global_sum(ReduceBuff,7);
	for (i = 0; i < 3; i++)
	    water_balance[i+1] = ReduceBuff[i];
	nzeros = (int)ReduceBuff[3];
	num_drops = (int)ReduceBuff[4];
	liquid_mass = ReduceBuff[5];
	liquid_mass0 = ReduceBuff[6];
#endif 

	if (!eqn_params->no_droplets && nzeros == 0)
	{
	     printf("No droplets included\n");
	     front->time_limit_reached = YES;
	}
	if(liquid_mass != liquid_mass0)
	{
	    alpha = (log(double(nzeros)/double(num_drops)))
		    /(log(liquid_mass/liquid_mass0));
	}
	else
	{
	    alpha = 0.0;
	}
	if (pp_mynode() == 0)
	{
	    for (i = 0 ; i < 4; i++)
	        fprintf(file,"%20.19f  ",water_balance[i]);
	
	    fprintf(file,"%20.19f\n",alpha);	
	    fclose(file);
	}
}

void VCARTESIAN::setDomain()
{
	static boolean first = YES;
	INTERFACE *grid_intfc;
	Table *T;
	int i;

	if (debugging("trace")) printf("Passed FT_MakeGridIntfc()\n");

	grid_intfc = front->grid_intfc;
	top_grid = &topological_grid(grid_intfc);
	lbuf = front->rect_grid->lbuf;
	ubuf = front->rect_grid->ubuf;
	top_gmax = top_grid->gmax;
	top_L = top_grid->L;
	top_U = top_grid->U;
	top_h = top_grid->h;
	dim = grid_intfc->dim;
	T = table_of_interface(grid_intfc);
	top_comp = T->components;
	eqn_params = (PARAMS*)front->extra2;
	hmin = top_h[0];
	for (i = 1; i < dim; ++i)
	    if (hmin > top_h[i]) hmin = top_h[i];

	if (field == NULL)
	    FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(PHASE_FIELD));

	switch (dim)
	{
	case 1:
            if (first == YES)
            {
		comp_size = top_gmax[0]+1;
                FT_VectorMemoryAlloc((POINTER*)&array,comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&source,comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->vapor,comp_size,
			FLOAT);
		FT_VectorMemoryAlloc((POINTER*)&field->supersat,comp_size,
			FLOAT);
		FT_VectorMemoryAlloc((POINTER*)&field->mrad,comp_size,
			FLOAT);
		FT_VectorMemoryAlloc((POINTER*)&field->temperature,comp_size,
                        FLOAT);
                first = NO;
            }	
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    eqn_params->field = field;
	    break;
	case 2:
	    if (first == YES)
	    {
		comp_size = (top_gmax[0]+1)*(top_gmax[1]+1);
	    	FT_VectorMemoryAlloc((POINTER*)&array,comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&source,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->vapor,
			comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->cloud,
			comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->supersat,comp_size,
                        FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->mrad,comp_size,
                        FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->drops,comp_size,
                        FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->adv,comp_size,
                        FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->adv_old,comp_size,
                        FLOAT);
		FT_VectorMemoryAlloc((POINTER*)&field->sample,comp_size,
			FLOAT);
                FT_MatrixMemoryAlloc((POINTER*)&field->ext_accel,2,comp_size,
                                        FLOAT);
		FT_VectorMemoryAlloc((POINTER*)&field->temperature,
                        comp_size,FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    eqn_params->field = field;
	    break;
	case 3:
	    if (first == YES)
	    {
		comp_size = (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1);
	    	FT_VectorMemoryAlloc((POINTER*)&array,comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&source,comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->vapor,
			comp_size,FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->cloud,
			comp_size,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->supersat,comp_size,
                        FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->mrad,comp_size,
                        FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->drops,comp_size,
                        FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->adv,comp_size,
                        FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->adv_old,comp_size,
                        FLOAT);
		FT_VectorMemoryAlloc((POINTER*)&field->sample,comp_size,
			FLOAT);
                FT_MatrixMemoryAlloc((POINTER*)&field->ext_accel,3,comp_size,
                                        FLOAT);
		FT_VectorMemoryAlloc((POINTER*)&field->temperature,
                        comp_size,FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
	    eqn_params->field = field;
	    break;
	}
	/*save global grid information*/
	eqn_params->global_grid = &(front->pp_grid->Global_grid);
}	/* end setDomain */

void VCARTESIAN::initMovieVariables()
{
	PARAMS *params = (PARAMS*)front->extra2;
	if (dim == 2)
	{	
 	    FT_AddHdfMovieVariable(front,NO,YES,SOLID_COMP,
                                "vapor",0,field->vapor,
				getStateVapor,0,0);
	    FT_AddHdfMovieVariable(front,NO,YES,SOLID_COMP,
                                "supersat",0,field->supersat,
				getStateSuper,0,0);
	    FT_AddHdfMovieVariable(front,NO,YES,SOLID_COMP,
                                "cloud",0,field->cloud,
				getStateSuper,0,0);
 	    FT_AddHdfMovieVariable(front,NO,YES,SOLID_COMP,
                                "temperature",0,field->temperature,
				getStateVapor,0,0);
            if (params->movie_option->plot_particles == YES)
                FT_AddVtkScalarMovieVariable(front,"Cloud",field->cloud);
            if (params->movie_option->plot_vapor == YES)
            	FT_AddVtkScalarMovieVariable(front,"Vapor",field->vapor);

	}
	else
	{
	    /* Added for vtk movie of scalar field */
            if (params->movie_option->plot_vapor == YES)
            	FT_AddVtkScalarMovieVariable(front,"Vapor",field->vapor);
            if (params->movie_option->plot_particles == YES)
                FT_AddVtkScalarMovieVariable(front,"Cloud",field->cloud);
	    if (params->movie_option->plot_temperature)
                FT_AddVtkScalarMovieVariable(front,"Temperature",field->temperature);
	}
        if (debugging("trace"))
            printf("Leaving initMovieVariables()\n");
}	/* end initMovieVariables */

static int find_state_at_crossing(
	Front *front,
	int *icoords,
        GRID_DIRECTION dir,
        int comp,
        POINTER *state,
        HYPER_SURF **hs,
        double *crx_coords)
{
	boolean status;
	INTERFACE *grid_intfc = front->grid_intfc;

	status = FT_StateStructAtGridCrossing(front,grid_intfc,icoords,dir,
				comp,state,hs,crx_coords);
        if (status == NO) return NO_PDE_BOUNDARY;

        if (wave_type(*hs) == FIRST_PHYSICS_WAVE_TYPE) 
	    return NO_PDE_BOUNDARY;
	else if (wave_type(*hs) == DIRICHLET_BOUNDARY)
	{
            return DIRICHLET_PDE_BOUNDARY;
	}
	else if (wave_type(*hs) == GROWING_BODY_BOUNDARY)
	{
            return DIRICHLET_PDE_BOUNDARY;
	}
	else if (wave_type(*hs) == NEUMANN_BOUNDARY)
            return NEUMANN_PDE_BOUNDARY;
}       /* find_state_at_crossing */

void VCARTESIAN::recordLagrangSupersat(const char *out_name)
{
	INTERFACE *intfc = front->interf;
	PARAMS* eqn_params = (PARAMS*)front->extra2;
	char fname[100];
	FILE *file;
	double *PDF;
	double *center, s;
	double min_supersat,max_supersat;
	double bin_size, bin_mid;
	int i,bin_num,index,ic[MAXD];
	static double *supersat_array;
	static int max_array_size = 0;
	
 	if (debugging("trace"))	
	    printf("Entering record Lagrang supersaturation\n");
	if (pp_mynode() == 0)
	{
	    sprintf(fname,"%s/PDF-lagsupersat",out_name);
            if (!create_directory(fname,NO))
            {
                printf("Cannot create directory %s\n",fname);
                clean_up(ERROR);
            }
            sprintf(fname,"%s/supersat-%4.2f",fname,front->time);
            file = fopen(fname,"w");
	}
	if (eqn_params->num_drops > max_array_size)
	{
	    max_array_size = eqn_params->num_drops;
	    free_these(1,supersat_array);
	    FT_VectorMemoryAlloc((POINTER*)&supersat_array,max_array_size,FLOAT);
	}
	for (i = 0; i < eqn_params->num_drops; i++)
	{
	    center = eqn_params->particle_array[i].center;
	    rect_in_which(center,ic,top_grid);
	    index = d_index(ic,top_gmax,dim);
	    s = field->supersat[index];
	    FT_IntrpStateVarAtCoords(front,LIQUID_COMP,center,
                                field->supersat,getStateSuper,&s,&s);
	    supersat_array[i] = s;
	}
	bin_num = 200;
	PDF = ComputePDF(supersat_array,eqn_params->num_drops,bin_size,bin_num,min_supersat,max_supersat);
	if (pp_mynode() == 0)
	{
	    for (i = 0; i < bin_num; i ++)
	    {
	        bin_mid = min_supersat+(0.5+i)*bin_size;
	        fprintf(file,"%15.14f  %15.14f\n",bin_mid,PDF[i]);
	    }
	    fclose(file);
	}
	if (debugging("trace"))
	    printf("Leaving record Lagrang supersaturation\n");
}

void VCARTESIAN::recordRadius(char *out_name)
{
	INTERFACE *intfc = front->interf;
	PARAMS* eqn_params = (PARAMS*)front->extra2;
	char fname[100];
	FILE *file;
	double *PDF;
	double min_radius,max_radius;
	double bin_size, bin_mid;
	int i,bin_num;
	static double *radius_array;
	static int max_array_size = 0;
	
 	if (debugging("trace"))	
	    printf("Entering record radius\n");
	if (pp_mynode() == 0)
	{
	    sprintf(fname,"%s/PDF-radius",out_name);
            if (!create_directory(fname,NO))
            {
                printf("Cannot create directory %s\n",fname);
                clean_up(ERROR);
            }
            sprintf(fname,"%s/radius-%4.2f",fname,front->time);
            file = fopen(fname,"w");
	}
	if (eqn_params->num_drops > max_array_size)
	{
	    max_array_size = eqn_params->num_drops;
	    free_these(1,radius_array);
	    FT_VectorMemoryAlloc((POINTER*)&radius_array,max_array_size,FLOAT);
	}
	for (i = 0; i < eqn_params->num_drops; i++)
	    radius_array[i] = eqn_params->particle_array[i].radius;
	bin_num = 200;
	boolean ignore_zero = YES;
	PDF = ComputePDF(radius_array,eqn_params->num_drops,bin_size,bin_num,min_radius,max_radius,ignore_zero);
	printf("max_radius = %e, min_radius = %e, %d drops contained, %d bins used\n",
		max_radius, min_radius, eqn_params->num_drops,bin_num);
	if (pp_mynode() == 0)
	{
	    for (i = 0; i < bin_num; i ++)
	    {
	        bin_mid = min_radius+(0.5+i)*bin_size;
	        fprintf(file,"%15.14f  %15.14f\n",bin_mid,PDF[i]);
	    }
	    fclose(file);
	}
	if (debugging("trace"))
	    printf("Leaving record radius\n");
}

void VCARTESIAN::initPresetParticles()
{
        int i,j,k,l,id,ic[MAXD];
        int index,count,G_count,G_size = 1;
        double ratio_in_subdomain;
        GAUSS_PARAMS gauss_params;
        UNIFORM_PARAMS uniform_params;
        double r_bar, sigma, x;
        unsigned short int xsubi[3];
        double cL[MAXD]; /*bound of a grid cell*/
        char msg[200];
        char *inname = InName(front);
        FILE *infile = fopen(inname,"r");
        double *supersat = field->supersat;
        PARAMS *eqn_params  = (PARAMS*)front->extra2;
        PARTICLE* particle_array = eqn_params->particle_array;
        double *local_L = front->pp_grid->Zoom_grid.L;
        double *local_U = front->pp_grid->Zoom_grid.U;
	G_size = 1;
	for (i = 0; i < dim; i++)
	    G_size *= front->pp_grid->Global_grid.gmax[i];
        if (particle_array != NULL)
            FT_FreeThese(1,particle_array);

        CursorAfterString(infile,"Enter number of water drops:");
        fscanf(infile,"%d",&eqn_params->num_drops);
        (void) printf("%d\n",eqn_params->num_drops);

        sprintf(msg,"Enter mean radius of water drop:");
        CursorAfterString(infile,msg);
        fscanf(infile,"%lf",&r_bar);
        (void) printf("%f\n",r_bar);
        sprintf(msg,"Enter standard deviation of radius:");
        CursorAfterString(infile,msg);
        fscanf(infile,"%lf",&sigma);
        (void) printf("%f\n",sigma);
        fclose(infile);
        xsubi[0] = 10;
        xsubi[1] = 100;
        xsubi[2] = 1000;

        gauss_params.mu = r_bar;
        gauss_params.sigma = sigma;
        uniform_params.a = 0.0;
        uniform_params.b = 1.0;

        for (i = 0; i < comp_size; i++)
                field->drops[i] = 0;
        count = 0;
        switch (dim)
        {
            case 2:
                for (j = jmin; j <= jmax; j++)
                for (i = imin; i <= imax; i++)
                {
                    index = d_index2d(i,j,top_gmax);
                    if(supersat[index] >= 0)
                        count ++;
                }
                break;
            case 3:
                for (k = kmin; k <= kmax; k++)
                for (j = jmin; j <= jmax; j++)
                for (i = imin; i <= imax; i++)
                {
                    index = d_index3d(i,j,k,top_gmax);
                    if (supersat[index] >= 0)
                        count ++;
                }
                break;
            default:
                printf("Unknow dim = %d\n",dim);
                clean_up(ERROR);
        }

        G_count = count;
#if defined(__MPI__)
        pp_gsync();
        pp_global_isum(&G_count,1);
#endif
        if (G_count != 0)
        {
            ratio_in_subdomain = double(count)/double(G_count);
            eqn_params->num_drops = int(ratio_in_subdomain*eqn_params->num_drops);
            printf("%d droplets in subdomain\n",eqn_params->num_drops);
            printf("%d cells(s>0) in global domain\n",G_count);
	    printf("Fraction of cloudy air = %f\n",(double)G_count/(double)G_size);
        }
        else
        {
            printf("No droplets in the domain!\n");
            clean_up(ERROR);
        }
        FT_VectorMemoryAlloc((POINTER*)&particle_array,
                              eqn_params->num_drops,sizeof(PARTICLE));
        count = 0;
        while(count < eqn_params->num_drops)
        {
            for (i = 0; i < dim; ++i)
            {
                x = dist_uniform((POINTER)&uniform_params,xsubi);
                particle_array[count].center[i]
                        = local_L[i]+x*(local_U[i]-local_L[i]);
                particle_array[count].vel[i] = 0.0;
                ic[i] = floor((particle_array[count].center[i]
                        - top_L[i] + 0.5*top_h[i])/top_h[i]);
            }
            index = d_index(ic,top_gmax,dim);
            if (supersat[index] < 0)
                continue;
            particle_array[count].radius
                        = gauss_center_limit((POINTER)&gauss_params,xsubi);
            particle_array[count].R0
                        = particle_array[count].radius;
            particle_array[count].flag = YES;
            particle_array[count].rho = eqn_params->rho_l;
            count ++;
        }
        count = 1;
        for (l = 0; l < dim; l++) count *= eqn_params->Nclip[l];
        if (count == 1)
            setParticleGlobalIndex(particle_array,eqn_params->num_drops);
        else
            setParticleGroupIndex(particle_array,eqn_params->num_drops,dim,
                                  eqn_params->Nclip,
                                  front->pp_grid->Global_grid.L,
                                  front->pp_grid->Global_grid.U);
        eqn_params->particle_array = particle_array;
	setSamplePoints(field->sample,supersat,comp_size);
        computeSource();
}

static void setSamplePoints(double* sample, double *refs, int size)
{
	for (int i = 0; i < size; i++)
	    sample[i] = (refs[i] >= 0) ? 1 : 0;
}

/*complex operations multiplication*/
static void CXC(fftw_complex a, fftw_complex b, fftw_complex ans)
{
	ans[0] = a[0]*b[0] - a[1]*b[1];
	ans[1] = a[1]*b[0] + a[0]*b[1];
}

/*complex operations division*/
static void CDC(fftw_complex a, fftw_complex b, fftw_complex ans)
{
	double deno = 0.0;
	deno = b[0]*b[0] + b[1]*b[1];
	if (deno == 0.0)
	{
	    ans[0] = ans[1] = 0.0;
	}
	else
	{
	    ans[0] = (a[0]*b[0] + a[1]*b[1])/deno;
	    ans[1] = (a[1]*b[0] - a[0]*b[1])/deno;
	}
}

static void gatherParallelData(
	    Front* front,
	    double* local_var,
	    fftw_complex *global_var)
{

	int i,j,k,l,dim,index0,lmin[MAXD],lmax[MAXD],local_size;
	int i1,j1,k1,index1,N[MAXD],id;
	INTERFACE* grid_intfc = front->grid_intfc;
        RECT_GRID *top_grid = &topological_grid(grid_intfc);
	PP_GRID* pp_grid = front->pp_grid;
	RECT_GRID *global_grid = &(pp_grid->Global_grid);
        int *top_gmax = top_grid->gmax;
	int *global_gmax = global_grid->gmax;
        int *lbuf = front->rect_grid->lbuf;
        int *ubuf = front->rect_grid->ubuf;
	int G_icoords[MAXD], icoords[MAXD], pp_icoords[MAXD];
	int fft_tag = 442;
	static double* local_buff;

	dim = top_grid->dim;
	local_size = 1;
	for (i = 0; i < dim; i++)
	{
	    lmin[i] = (lbuf[i] == 0) ? 1 : lbuf[i];
	    lmax[i] = (ubuf[i] == 0) ? top_gmax[i] - 1 : top_gmax[i] - ubuf[i];
	    N[i] = lmax[i] - lmin[i] + 1;
	    local_size *= (top_gmax[i] + 1); 
	}
	/*make a bufff for communications*/
	if(local_buff == NULL)
	    local_buff = new double[local_size];
	for (i = 0; i < local_size; i++)
	    local_buff[i] = local_var[i];

	pp_gsync();
	switch(dim)
	{
	    case 2:
		if (pp_mynode() != 0)	
		{
		    pp_send(fft_tag,(POINTER)local_buff,local_size*sizeof(double),0);
		}
		else
		{
		    for (id = 0; id < pp_numnodes(); id++)
		    {
			if (id != 0)
			    pp_recv(fft_tag,id,(POINTER)local_buff,local_size*sizeof(double));
			find_Cartesian_coordinates(id,pp_grid,pp_icoords);
			/*send local var to global var for first processor*/
                        for (i = lmin[0]; i <= lmax[0]; i++)
                        for (j = lmin[1]; j <= lmax[1]; j++)
                        {
                            icoords[0] = i; icoords[1] = j;
                            for(l = 0; l < dim; l++)
                                G_icoords[l] = icoords[l] - lmin[l] + pp_icoords[l]*N[l]; 
                            index0 = d_index(icoords,top_gmax,dim);
                            index1 = G_icoords[1] + (global_gmax[1]) * G_icoords[0];
                            global_var[index1][0] = local_buff[index0];
                        }
		    }
		}
		break;
	    case 3:
		if (pp_mynode() != 0)	
		{
		    pp_send(fft_tag,(POINTER)local_buff,local_size*sizeof(double),0);
		}
		else
		{
		    for (id = 0; id < pp_numnodes(); id++)
		    {
			if (id != 0)
		 	    pp_recv(fft_tag,id,(POINTER)local_buff,local_size*sizeof(double));
			/*send local var to global var*/
			find_Cartesian_coordinates(id,pp_grid,pp_icoords);
                        for (i = lmin[0]; i <= lmax[0]; i++)
                        for (j = lmin[1]; j <= lmax[1]; j++)
		        for (k = lmin[2]; k <= lmax[2]; k++)
                        {
			    icoords[0] = i; icoords[1] = j; icoords[2] = k;
                            for(l = 0; l < dim; l++)
                                G_icoords[l] = icoords[l] - lmin[l] + pp_icoords[l]*N[l]; 
                            index0 = d_index(icoords,top_gmax,dim);
                            index1 = G_icoords[2] + (global_gmax[2])
				   *(G_icoords[1] + (global_gmax[1])*G_icoords[0]);
                            global_var[index1][0] = local_buff[index0];
                        }
		    }
		}
		break;
	    default:
		printf("dim can only be 2,3, unknown dim %d\n",dim);
		clean_up(ERROR);
	}	
}

static void scatterParallelData(
	    Front* front,
	    fftw_complex *global_var, 
	    double* local_var)
{
	int i,j,k,l,dim,index0,lmin[MAXD],lmax[MAXD],local_size;
	int i1,j1,k1,index1,N[MAXD],id;
	INTERFACE* grid_intfc = front->grid_intfc;
        RECT_GRID *top_grid = &topological_grid(grid_intfc);
	PP_GRID* pp_grid = front->pp_grid;
	RECT_GRID *global_grid = &(pp_grid->Global_grid);
        int *top_gmax = top_grid->gmax;
	int *global_gmax = global_grid->gmax;
        int *lbuf = front->rect_grid->lbuf;
        int *ubuf = front->rect_grid->ubuf;
	int G_icoords[MAXD], icoords[MAXD], pp_icoords[MAXD];
	int fft_tag = 442;

	dim = top_grid->dim;
	local_size = 1;
	for (i = 0; i < dim; i++)
	{
	    lmin[i] = (lbuf[i] == 0) ? 1 : lbuf[i];
	    lmax[i] = (ubuf[i] == 0) ? top_gmax[i] - 1 : top_gmax[i] - ubuf[i];
	    N[i] = lmax[i] - lmin[i] + 1;
	    local_size *= (top_gmax[i] + 1); 
	}
	pp_gsync();
	switch(dim)
	{
	    case 2:
		if (pp_mynode() == 0)
		{
		    for (id = pp_numnodes()-1; id >= 0; id--)
		    {
			find_Cartesian_coordinates(id,pp_grid,pp_icoords);
			for (i = lmin[0]; i <= lmax[0]; i++)
                    	for (j = lmin[1]; j <= lmax[1]; j++)
			{
			    index0 = d_index2d(i,j,top_gmax);
			    i1 = i - lmin[0] + pp_icoords[0]*N[0];
			    j1 = j - lmin[1] + pp_icoords[1]*N[1];
			    index1 =  j1 + (global_gmax[1])*i1;
			    local_var[index0] = global_var[index1][0];
			}
			if (id != 0)
			     MPI_Send((POINTER)local_var,local_size,MPI_DOUBLE,id,fft_tag,MPI_COMM_WORLD);
			    /*pp_send(fft_tag,(POINTER)local_var,local_size*sizeof(double),id);*/	
		    }
		}
		else
		{
		    pp_recv(fft_tag,0,(POINTER)local_var,local_size*sizeof(double));
		}
		break;
	    case 3:
		if (pp_mynode() == 0)
		{
		    for (id = pp_numnodes()-1; id >= 0; id--)
		    {
			find_Cartesian_coordinates(id,pp_grid,pp_icoords);
			for (i = lmin[0]; i <= lmax[0]; i++)
                    	for (j = lmin[1]; j <= lmax[1]; j++)
                    	for (k = lmin[2]; k <= lmax[2]; k++)
			{
			    index0 = d_index3d(i,j,k,top_gmax);
			    i1 = i - lmin[0] + pp_icoords[0]*N[0];
			    j1 = j - lmin[1] + pp_icoords[1]*N[1];
			    k1 = k - lmin[2] + pp_icoords[2]*N[2];
                     	    index1 = k1 + (global_gmax[2])
				   *(j1 + (global_gmax[1]) * i1);
			    local_var[index0] = global_var[index1][0];
			}
			if (id != 0)
			     MPI_Send((POINTER)local_var,local_size,MPI_DOUBLE,id,fft_tag,MPI_COMM_WORLD);
			    //pp_send(fft_tag,(POINTER)local_var,local_size*sizeof(double),id);	
		    }
		}
		else
		{
		    pp_recv(fft_tag,0,(POINTER)local_var,local_size*sizeof(double));
		}
		break;
	    default:
		printf("dim can only be 2,3, unknown dim %d\n",dim);
		clean_up(ERROR);
	}
}

static double computeInputEnergy
		(double** ext_accel,
		 double** vel, 
		 int dim, 
		 int *imin,
		 int *imax,
		 int *top_gmax)
{
	int i, j, k, l, index, size;
	double eps_in = 0.0;
	double max_a = 0;
	size = 1;
	for (i = 0; i < dim; i++)
	    size *= (imax[i]-imin[i]+1);
	switch(dim)
	{
	    case 2:
		for (i = imin[0]; i <= imax[0]; i++)
		for (j = imin[1]; j <= imax[1]; j++)
		for (l = 0; l < dim; l++)
		{
		    index = d_index2d(i,j,top_gmax);
	    	    eps_in += ext_accel[l][index]*vel[l][index];
		    max_a = FT_Max(max_a,fabs(ext_accel[l][index]));
		}
		break;
	    case 3:
		for (i = imin[0]; i <= imax[0]; i++)
		for (j = imin[1]; j <= imax[1]; j++)
		for (k = imin[2]; k <= imax[2]; k++)
		for (l = 0; l < dim; l++)
		{
		    index = d_index3d(i,j,k,top_gmax);
	    	    eps_in += ext_accel[l][index]*vel[l][index];
		    max_a = FT_Max(max_a,fabs(ext_accel[l][index]));
		}
		break;
	    Default:
		printf("Warning: in computeInputEnergy(), unknown dim = %d\n",dim);
	}	
	printf("max_accel = %f\n",max_a);
	pp_gsync();
	pp_global_isum(&size,1);
	pp_global_sum(&eps_in,1);
	return eps_in/size;
}

static double computeInputEnergy(
		fftw_complex* U, 
		fftw_complex* f, 
		int size)
{
	int i, j;
	double eps_in;
	fftw_complex ans;
	eps_in = 0.0;
	for (i = 0; i < size; i++)
	{
		CXC(U[i],f[i],ans);
		eps_in += ans[0]; 
	}
	return eps_in;
}
/*compute isotropic volume force maintaining turbulence*/
static double computeUrms(double **vel,int dim,int *lmin,int *lmax,int *top_gmax)
{
	int i, j, k, l, index;
	double urms = 0.0;
	int size = 1;
	for (i = 0; i < dim; i++)
	    size *= (lmax[i] - lmin[i] + 1);
	switch(dim)
	{
	    case 2:
		for (i = lmin[0]; i <= lmax[0]; i++)
		for (j = lmin[1]; j <= lmax[1]; j++)
		for (l = 0; l < dim; l++)
		{
		    index = d_index2d(i,j,top_gmax);
		    urms += sqr(vel[l][index]);
		}
		break;
	    case 3:
		for (i = lmin[0]; i <= lmax[0]; i++)
		for (j = lmin[1]; j <= lmax[1]; j++)
		for (k = lmin[2]; k <= lmax[2]; k++)
		for (l = 0; l < dim; l++)
		{
		    index = d_index3d(i,j,k,top_gmax);
		    urms += sqr(vel[l][index]);
		}
		break;
	    default:
		printf("In computeUrms(), Unknown dim = %d\n",dim);
		clean_up(ERROR);
	}
	pp_gsync();
	pp_global_sum(&urms,1);
	pp_global_isum(&size,1);
	urms = (urms/dim)/size;
	return sqrt(urms);
}

void VCARTESIAN::computeVolumeForce()
{
	return computeVolumeForceFourier();
}

void VCARTESIAN::computeVolumeForceLinear()
{
	int i,j,k;
	double urms,A;
	double **ext_accel = field->ext_accel;
	double **vel = field->vel;
	int lmin[MAXD], lmax[MAXD];
	static boolean first = YES;
	static double eps;
	double eps_in;

	if (eqn_params->if_volume_force == NO)
            return;
	printf("computeVolumeForceLinear\n");
	if (first)
	{
	    first = NO;
	    /*compute mean kinetic energy dissipation*/
	    eps = 0.2*computeDspRateLinear(); 	
            eqn_params->disp_rate = eps;
	}
	eqn_params->disp_rate = computeDspRateLinear();
        printf("eps_in = %e, eps_out = %e\n",eps,eqn_params->disp_rate);
	lmin[0] = imin; lmin[1] = jmin;
	lmax[0] = imax; lmax[1] = jmax;
	if (dim == 3)
	{lmin[2] = kmin; lmax[2] = kmax;}
	urms = computeUrms(vel,dim,lmin,lmax,top_gmax);
	printf("urms = %f\n",urms);
	A = eps/(dim*urms*urms);
	printf("A = %f\n",A);
	
	for (i = 0; i < comp_size; i++)
	for (j = 0; j < dim; j++)
	    ext_accel[j][i] = A*vel[j][i];
	for (j = 0; j < dim; j++)
	    FT_ParallelExchGridArrayBuffer(ext_accel[j],front,NULL);
	eps_in = computeInputEnergy(ext_accel,vel,dim,lmin,lmax,top_gmax);
        printf("Phy: eps_in = %e\n",eps_in);
}

void VCARTESIAN::computeVolumeForceFourier()
{
	int i, j, k, l, index0;
	int i1, j1, k1, index1,icrds[MAXD];
	static int Nr, N[MAXD], N0;
	double **ext_accel = field->ext_accel;
	double **vel = field->vel;
	static fftw_complex *U, *V, *W, *fx, *fy, *fz;
	static boolean first = YES;
	static double eps; /*mean kinetic energy dissipation*/
	fftw_complex deno,ans;
	int count = 0;

	FILE* outfile;
	char filename[100];

	if (debugging("volume_force"))
	{
	    sprintf(filename,"%s/vol_forc-%d",
		    OutName(front),front->step);
	    if (pp_mynode() == 0)
	        outfile = fopen(filename,"w");
	}
	if (dim == 2)
	    N0 = 2; /*wave number shell |N|^2 = N0*/
	else
	    N0 = 6; /* 1^2 + 1^2 + 2^2 = 6*/
	
	if (eqn_params->if_volume_force == NO)
	    return;
	printf("computeVolumeForceFourier\n");
	start_clock("volume_force");


	/*for parallelization, global rectangle grid*/
	PP_GRID *pp_grid = front->pp_grid;
	RECT_GRID *global_grid = &(pp_grid->Global_grid);
	
	if (first == YES)
	{
	    first = NO;
	    eps = computeDspRate(); 	    
	    eqn_params->disp_rate = eps;
	    /*compute mean kinetic energy dissipation*/
	    printf("eps_in = %e\n",eps);
	    if (pp_mynode() == 0)
	    {
	        /*allocating memory*/
	        for ( i = 0; i < dim; i++)
	            N[i] = global_grid->gmax[i];
	        /*N should be the global gmax*/
	        Nr = 1;
	        for (i = 0; i < dim; i++)
		    Nr *= N[i]; 
	        U = new fftw_complex[Nr];
	        V = new fftw_complex[Nr];
	        fx = new fftw_complex[Nr];
	        fy = new fftw_complex[Nr];
	        if (dim == 3)
	        {
	            W = new fftw_complex[Nr];
		    fz = new fftw_complex[Nr];
	        }
	    }
	}
	eqn_params->disp_rate = computeDspRate();
	printf("FFT: esp_in = %e, eps_out = %e\n",eps,eqn_params->disp_rate);
	/*collect data in master processor*/
	gatherParallelData(front,vel[0],U);	
	gatherParallelData(front,vel[1],V);	
	if (dim == 3)
	    gatherParallelData(front,vel[2],W);	

	if (pp_mynode() == 0)
	{
	switch (dim)
	{
	    case 2:
		for (i = 0; i < Nr; i++)
		{
		     U[i][1] = V[i][1] = 0.0;
		     fx[i][0] = fx[i][1] = 0.0;
		     fy[i][0] = fy[i][1] = 0.0;
		}
		/*Transform to Fourier space*/
		fftnd(U,dim,N,1);
		fftnd(V,dim,N,1);
		count = 0;
		/*construct forcing term in Fourier space*/
		for (i1 = -N0; i1 <= N0; i1++)
		for (j1 = 0; j1 <= N0; j1++)
		{
		    if ((i1*i1 + j1*j1) != N0)
			continue;
		    count ++;
		    icrds[0] = i1; icrds[1] = j1;
                    for (l = 0; l < dim; l++)
                    {
                        if (icrds[l] < 0)
                            icrds[l] = N[l]+icrds[l];
                    }
                    index1 = icrds[1]+(N[1]/2+1)*icrds[0];
		    deno[0] = sqr(U[index1][0])+sqr(U[index1][1])
		    	    + sqr(V[index1][0])+sqr(V[index1][1]);
		    for (l = 0; l <= 1; l++)
		    {
		       fx[index1][l] = eps * U[index1][l]/deno[0]; 
		       fy[index1][l] = eps * V[index1][l]/deno[0];
		    }
		    
   	            printf("U[%d %d] = [%e %e], fx = [%e %e]\n",
				i1,j1,U[index1][0],U[index1][1],
				fx[index1][0],fx[index1][1]);
		    printf("V[%d %d] = [%e %e], fy = [%e %e]\n",
				i1,j1,V[index1][0],V[index1][1],
				fy[index1][0],fy[index1][1]);
		}
	 	for (i = 0; i < N[0]; i++)
                for (j = 0; j < N[1]/2+1; j++)
		{
		   index1 = j + (N[1]/2+1)*i;
		   for (l = 0; l <= 1; l++)
		   {
		   fx[index1][l] /= (2*count);
		   fy[index1][l] /= (2*count);	
		   }
		}
		/*Transform back to physics space*/
		fftnd(fx,dim,N,-1);
        	fftnd(fy,dim,N,-1);
		if (debugging("volume_force"))
		{
		    printf("N = [%d %d]\n",N[0],N[1]);
		    for (i = 0; i < N[0]; i++)
                    for (j = 0; j < N[1]; j++)
		    {
		        index0 = j + N[1]*i;
		        fprintf(outfile,"%e %e\n",
			fx[index0][0],fy[index0][0]);
		    }
		    fclose(outfile);
		}
		break;
	  case 3:
                /*maping vel field to row major format*/
                for (index1 = 0; index1 < Nr; index1++)
                {
                     U[index1][1] = V[index1][1] = W[index1][1]=0.0;
                     fx[index1][0] = fx[index1][1] = 0.0;
                     fy[index1][0] = fy[index1][1] = 0.0;
                     fz[index1][0] = fz[index1][1] = 0.0;
                }
                /*Transform to Fourier space*/
                fftnd(U,dim,N,1);
                fftnd(V,dim,N,1);
                fftnd(W,dim,N,1);
                /*construct forcing term in Fourier space*/
	 	count = 0;
		for (i1 = -N0; i1 <= N0; i1++)
		for (j1 = -N0; j1 <= N0; j1++)
		for (k1 =   0; k1 <= N0; k1++)
                {
		    if ((i1*i1 + j1*j1 + k1*k1) != N0)
			continue;
		    count ++;
		    icrds[0] = i1; icrds[1] = j1; icrds[2] = k1;
		    for (l = 0; l < dim; l++)
                    {
                        if (icrds[l] < 0)
                            icrds[l] = N[l]+icrds[l];
                    }
	   	    index1 = icrds[2]+(N[2]/2+1)
			    *(icrds[1]+N[1]*icrds[0]);
		    deno[0] = sqr(U[index1][0])+sqr(U[index1][1])
                            + sqr(V[index1][0])+sqr(V[index1][1])
			    + sqr(W[index1][0])+sqr(W[index1][1]);
                    for (l = 0; l <= 1; l++)
                    {
                        fx[index1][l] = eps * U[index1][l]/deno[0];
                        fy[index1][l] = eps * V[index1][l]/deno[0];
                        fz[index1][l] = eps * W[index1][l]/deno[0];
                    }
		   
		    if (debugging("volume_force"))
		    {
   	                printf("U[%d %d %d] = [%e %e], fx = [%e %e]\n",
				i1,j1,k1,U[index1][0],U[index1][1],
				fx[index1][0],fx[index1][1]);
		        printf("V[%d %d %d] = [%e %e], fy = [%e %e]\n",
				i1,j1,k1,V[index1][0],V[index1][1],
				fy[index1][0],fy[index1][1]);
		        printf("W[%d %d %d] = [%e %e], fz = [%e %e]\n",
				i1,j1,k1,W[index1][0],W[index1][1],
				fz[index1][0],fz[index1][1]);
		    }
                }
	        for (i = 0; i < N[0];   i++)
               	for (j = 0; j < N[1];   j++)
               	for (k = 0; k < N[2]/2+1; k++)
                {
                    index1 = k + (N[2]/2+1)*(j + N[1]*i);
                    for (l = 0; l <= 1; l++)
                    {
                        fx[index1][l] /= (2*count);
                        fy[index1][l] /= (2*count);
                   	fz[index1][l] /= (2*count);
                    }
                }

		/*Transform force back to physical space*/
		fftnd(fx,dim,N,-1);
        	fftnd(fy,dim,N,-1);
        	fftnd(fz,dim,N,-1);
		if (debugging("volume_force"))
		{
		    for (k = 0; k < N[2]; k++)
                    for (j = 0; j < N[1]; j++)
		    for (i = 0; i < N[0]; i++)
		    {
		        index0 = k + N[2]*(j + N[1]*i);
		        fprintf(outfile,"%e %e %e\n",
			fx[index0][0],fy[index0][0],fz[index0][0]);
		    }
		    fclose(outfile);
		}

                break;	
	    Default:
		printf("Unknown dim = %d\n",dim);
		clean_up(ERROR);
	}
	}
	scatterParallelData(front,fx,ext_accel[0]);
	scatterParallelData(front,fy,ext_accel[1]);
	if (dim == 3)
	    scatterParallelData(front,fz,ext_accel[2]);
	for (i = 0; i < dim; i++)
	    FT_ParallelExchGridArrayBuffer(ext_accel[i],front,NULL);
	
	/*for test*/
	int lmin[MAXD], lmax[MAXD];
	double eps_in;
	lmin[0] = imin; lmin[1] = jmin; lmin[2] = kmin;
	lmax[0] = imax; lmax[1] = jmax; lmax[2] = kmax;
	eps_in = computeInputEnergy(ext_accel,vel,dim,lmin,lmax,top_gmax);
	printf("Phy: eps_in = %e\n",eps_in);	
	stop_clock("volume_force");
	return;
}

static void computeFluctuation(Front* front, double **ext_accel, int size, int dim)
{
	int i, j;
	FILE* outfile;
	char filename[256]; 
	double mean_buoyancy[MAXD] = {0.0, 0.0, 0.0};
	for (i = 0; i < size; i++)
	for (j = 0; j < dim; j++)
	{
	    mean_buoyancy[j] += ext_accel[j][i];
	}
	
	/*find mean buoyancy*/
	pp_gsync();
	pp_global_sum(mean_buoyancy,3);
	for (j = 0; j < dim; j++)
	    mean_buoyancy[j] /= (pp_numnodes()*size);
	if (pp_mynode() == 0)
        {
            sprintf(filename,"%s/buoyancy",OutName(front));
            outfile = fopen(filename,"a");
            fprintf(outfile,"%f  %15.14f\n",
		    front->time,mean_buoyancy[dim-1]);
            fclose(outfile);
        }
	/*remove mean buoyancy to obtain neutral buoyancy*/
	for (i = 0; i < size; i++)
        {
            ext_accel[0][i] = 0.0;
            ext_accel[dim-2][i] = 0.0;
            ext_accel[dim-1][i] -= mean_buoyancy[dim-1];
        }
}

void VCARTESIAN::computeSource()
{
	computeTemperatureSource();
	computeVaporSource();
}

void VCARTESIAN::computeVaporSource()
{
	/*compute condensation rate*/
	/*this is a source term for vapor mixing ratio equation*/
        int i, j, index, size, num_drops;
        int ic[MAXD];
        PARAMS* eqn_params = (PARAMS*)front->extra2;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        PARTICLE* particle_array = eqn_params->particle_array;
        num_drops = eqn_params->num_drops;
	double *supersat = eqn_params->field->supersat;
	double *T = eqn_params->field->temperature;
	double *qv = eqn_params->field->vapor;
	double *qc = eqn_params->field->cloud;
	double  T0 = eqn_params->T0;
	double  q0 = eqn_params->qv0;
        static boolean first = YES;
        double *coords;
	/*for computing the coefficient*/
	double rho_0 = iFparams->rho2;
	double a3 = 1.0;
	double coeff;

	static double maxsource = -HUGE, minsource = HUGE;

	for(i = 0; i < dim; i++)
	    a3 *= top_h[i];

        for (i = 0; i < comp_size; i++)
        {
	    source[i] = 0.0;
	    field->drops[i] = 0;
	    field->cloud[i] = 0.0;
	    field->mrad[i] = 0.0;
	}
	/*caculate num_drops in each cell: drops[index]*/
	/*compute source term for vapor equation: source[index]*/
	/*compute cloud water mixing ratio: qc[index]*/
        for (i = 0; i < num_drops; i++)
        {
            coords = particle_array[i].center;
            //rect_in_which(coords,ic,top_grid);
            for (j = 0; j < dim; j++)
                ic[j] = floor((coords[j] - top_L[j] + 0.5*top_h[j])/top_h[j]);
            index = d_index(ic,top_gmax,dim);

	    coeff = 4.0*PI*particle_array[i].rho * eqn_params->K
					       / (rho_0 * a3);
	    if (eqn_params->if_condensation)
	    	source[index] += -1000.0 * coeff * supersat[index]
				     	* particle_array[i].radius;
	    if (particle_array[i].radius > 0)
	    {
	        field->drops[index] += 1;
		field->mrad[index] += particle_array[i].radius;
		qc[index] += (4.0/3.0)*PI
				    *   pow(particle_array[i].radius,3)
				    *   particle_array[i].rho
				    /   (a3 * rho_0);
	    }
        }
	FT_ParallelExchGridArrayBuffer(qc,front,NULL);
	FT_ParallelExchGridArrayBuffer(field->drops,front,NULL);
	FT_ParallelExchGridArrayBuffer(field->mrad,front,NULL);
	/*compute mean radius in a cell*/
	for (index = 0; index < comp_size; index++)
	{
		if (field->drops[index] != 0)
		    field->mrad[index] /= field->drops[index];
	}
	/*compute source for Navier Stokes equation:ext_accel[dim][index]*/
	if (eqn_params->if_volume_force)
	    computeVolumeForce();
	else
	{
	    for (index = 0; index < comp_size; index++)
	    for (j = 0; j < dim; j++)
	    {
	        field->ext_accel[j][index] = -iFparams->gravity[j]
	         *((T[index]-T0)/T0 + 0.608 * 0.001 * (qv[index] - q0) - qc[index]);
   	    }
	    for (j = 0; j < dim; j++)
	    FT_ParallelExchGridArrayBuffer(field->ext_accel[j],front,NULL);
	    /*remove mean value to obtain neutral buoyancy*/
	    computeFluctuation(front,field->ext_accel,comp_size,dim);
	}
}

void VCARTESIAN::computeTemperatureSource()
{
        /*compute condensation rate*/
        /*this is a source term for vapor mixing ratio equation*/
        int i, j, index, size, num_drops;
        int ic[MAXD];
        PARAMS* eqn_params = (PARAMS*)front->extra2;
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        PARTICLE* particle_array = eqn_params->particle_array;
        num_drops = eqn_params->num_drops;
        double *supersat = eqn_params->field->supersat;
        static boolean first = YES;
        double *coords;
        /*for computing the coefficient*/
        double rho_0 = iFparams->rho2;
        double a3 = 1.0;
        double coeff;
        double L = eqn_params->Lh;
        double cp = eqn_params->Cp;

        for(i = 0; i < dim; i++)
            a3 *= top_h[i];

        for (i = 0; i < comp_size; i++)
            source[i] = 0.0;

        for (i = 0; i < num_drops; i++)
        {
            coords = particle_array[i].center;
            //rect_in_which(coords,ic,top_grid);
            for (j = 0; j < dim; j++)
                ic[j] = floor((coords[j] - top_L[j] + 0.5*top_h[j])/top_h[j]);
            index = d_index(ic,top_gmax,dim);
            coeff = 4.0*PI*particle_array[i].rho * eqn_params->K
                                               / (rho_0 * a3);
            if (eqn_params->if_condensation)
                source[index] += L/cp * coeff * supersat[index]
                                      * particle_array[i].radius;
        }
        FT_ParallelExchGridArrayBuffer(source,front,NULL);
}

void VCARTESIAN::output()
{
	char* out_name = OutName(front);
	recordRadius(out_name);
        recordPDF(out_name,"all");
        recordCondensationRate(out_name);
        recordWaterBalance();
        recordMixingLine();
        if (eqn_params->Nclip[0] != 1 ||
            eqn_params->Nclip[1] != 1 ||
            eqn_params->Nclip[2] != 1)
            recordGroupParticles();
        if(eqn_params->is_bigdata)
            recordParticles();
}
