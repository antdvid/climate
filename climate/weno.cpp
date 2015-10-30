/* WENO_5 flux for solving scalar equations*/
/* f(T) = u*Tx + v*Ty */
#include "climate.h"

static double linear_flux(
	double wave_speed,
	double u)
{
        return u;
}

static double linear_dF(
	double wave_speed,
	double u)
{
	return 1.0;
}

WENO_FLUX::WENO_FLUX(Front &front):front(&front)
{
}

void WENO_FLUX::Weno5_Get_Flux(
	double *u, /*scalar field*/
	double *v, /*wave speed*/ 
	double *flux, 
	double lambda,
	int n)
{
	int nrad = 3; /* 5th-order WENO */
	int i,j,k;
	int extend_size = n + 2*nrad;
	const double eps = 1.e-8;
	const int p = 2;
	static double *f; /* f(u(x)) */
	static double *fpL, *fpR;
	static int max_n = 0;
	double max_df, norm;
	double aa;

	/* coefficients for 2rd-order ENO interpolation stencil */
	double a[3][3] = {{1.0/3.0, -7.0/6.0, 11.0/6.0}, 
			  {-1.0/6.0, 5.0/6.0, 1.0/3.0}, 
			  {1.0/3.0, 5.0/6.0, -1.0/6.0}};

	/* Optimal weights C_k */
	double c[3] = {0.1,0.6,0.3};

	double is[3]; /* a smoothness measurement of the flux function */
	double alpha[3];
	double omega[3]; /* weights for stencils */
	double q[3]; /* ENO approximation of flux */
	double sum;

	if (max_n < n)
        {
            max_n = n;
            if (f != NULL)
            {
                FT_FreeThese(3,f,fpL,fpR);
            }
            FT_VectorMemoryAlloc((POINTER*)&f,extend_size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&fpL,n+1,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&fpR,n+1,sizeof(double));
        }

	/* Find the maximum of fabs(df(u))*/
	max_df = 0;
	for (i = 0; i < extend_size; ++i) 
	{
	    norm = dF(v[i],u[i]) > 0 ? dF(v[i],u[i]) : -dF(v[i],u[i]);
	    max_df = max_df > norm ? max_df : norm;
	}

	/*compute fL = fpL + fmL, where fmL = 0*/
	/* f[i] = 0.5 * (f_{i - nrad} + max_df * u[i])
	 * ensure that df[i] > 0*/
	for (i = 0; i < extend_size; ++i)
	    f[i] = 0.5 * (flux_func(v[i],u[i]) + max_df * u[i]);

	for (j = 0; j < n + 1; ++j)
	/* To approximate flux at x_{j+1/2}, we use stencil
	 S_k = {x_{j+k-2}, x_{j+k-1}, x_{j+k}} = {x[j+k-2+nrad],
	 x[j+k-1+nrad], x[j+k+nrad]} */
	{
		/* compute the weights omega[] */
	    is[0] = 13.0/12.0*(f[j] - 2.0*f[j+1] + f[j+2])*(f[j] - 2.0
		    *f[j+1] + f[j+2]) + 0.25*(f[j] - 4.0*f[j+1] + 3.0
		    *f[j+2])*(f[j] - 4.0*f[j+1] + 3.0 * f[j+2]);
	    is[1] = 13.0/12.0*(f[j+1] - 2.0*f[j+2] + f[j+3])
		    *(f[j+1] - 2.0*f[j+2] + f[j+3]) + 0.25*(f[j+1]
		    -f[j+3])*(f[j+1] - f[j+3]);
	    is[2] = 13.0/12.0*(f[j+2] - 2.0*f[j+3] + f[j+4])
		    *(f[j+2] - 2.0*f[j+3] + f[j+4]) + 0.25*(3.0*f[j+2] 
		    - 4.0*f[j+3] + f[j+4])*(3.0*f[j+2] - 4.0*f[j+3] + f[j+4]);

	    for (i = 0; i < nrad; ++i)
		alpha[i] = c[i]/pow(eps + is[i],p);

	    sum = alpha[0] + alpha[1] + alpha[2];
	    for (i = 0; i < nrad; ++i)
		omega[i] = alpha[i]/sum;

		/* compute ENO approximation of the flux */
	    for (i = 0; i < nrad; ++i) 
	    {
		q[i] = 0.0;
		for (k = 0; k < nrad; ++k)
		    q[i] += a[i][k] * f[j + i + k];
	    }

	    /* compute the linear combination of the r candidate
	     stencils to get a higher order
	     approximation of the flux */
	    fpL[j] = 0.0;
	    for (i = 0; i < nrad; ++i)
		fpL[j] += omega[i]*q[i];
	}

	/* compute fR = fpR + fmR, where fmR = 0*/
	/* f[i] = 0.5 * (f_{i - nrad} - max_df * u[i])
	 * ensure that df[i] < 0*/
	for (i = 0; i < extend_size; ++i)
	    f[i] = 0.5*(flux_func(v[i],u[i]) + max_df*u[i]);
	for (j = 0; j < n + 1; ++j)
	/* To approximate flux at x_{j+1/2}, we use stencil S_k =
	 {x_{j+1-k+2}, x_{j+1-k+1}, x_{j+1-k}} = {x[j+1-k+2+nrad],
	 x[j+1-k+1+nrad], x[j+1-k+nrad]} */

	{
	    /* compute the weights omega[] */
	    is[0] = 13.0/12.0*(f[j+5] - 2.0*f[j+4] + f[j+3])
		    *(f[j+5] - 2.0*f[j+4] + f[j+3]) + 0.25*(f[j+5]
		    - 4.0*f[j+4] + 3.0*f[j+3])*(f[j+5] - 4.0*f[j+4]
		    + 3.0*f[j+3]);
	    is[1] = 13.0/12.0*(f[j+4] - 2.0*f[j+3] + f[j+2])
		    *(f[j+4] - 2.0*f[j+3] + f[j+2]) + 0.25*(f[j+4]
		    - f[j+2])*(f[j+4] - f[j+2]);
	    is[2] = 13.0/12.0*(f[j+3] - 2.0*f[j+2] + f[j+1])
		    *(f[j+3] - 2.0*f[j+2] + f[j+1]) + 0.25*(3.0*f[j+3] 
		    - 4.0*f[j+2] + f[j+1])*(3.0*f[j+3] - 4.0*f[j+2] + f[j+1]);

	    for (i = 0; i < nrad; ++i)
		alpha[i] = c[i]/pow(eps + is[i], p);

	    sum = alpha[0] + alpha[1] + alpha[2];
	    for (i = 0; i < nrad; ++i)
		omega[i] = alpha[i]/sum;

		/* compute ENO approximation of the flux */
	    for (i = 0; i < nrad; ++i) 
	    {
		q[i] = 0.0;
		for (k = 0; k < nrad; ++k)
		    q[i] += a[i][k]*f[j+5-i-k];
	    }

		/* compute the linear combination of the r candidate stencils
	     to get a higher order approximation of the flux */
	    fpR[j] = 0.0;
	    for (i = 0; i < nrad; ++i)
		fpR[j] += omega[i]*q[i];
	}
	
	/*upwinding strategy*/
	for (j = 0; j < n; ++j)
	{
	    double fL, fR;
	    aa = 0.5*(v[j+nrad] + v[j-1+nrad]);
	    if (aa >= 0)
		fL = fpL[j];
	    else
		fL = fpR[j];

	    aa = 0.5*(v[j+nrad] + v[j+1+nrad]);
	    if (aa >= 0)
		fR = fpL[j+1];
	    else
		fR = fpR[j+1];

	    flux[j+nrad] = aa*lambda*(fR - fL);
	}
}

void WENO_FLUX::setSolverDomain()
{
	static boolean first = YES;
	RECT_GRID *rgr = &topological_grid(front->grid_intfc);
        struct Table *T = table_of_interface(front->grid_intfc);
        int *lbuf = front->rect_grid->lbuf;
        int *ubuf = front->rect_grid->ubuf;
	int i;
	
	nrad = 3;
	dF = linear_dF;
	flux_func = linear_flux;

	dim = Dimension(front->interf);
        top_comp = T->components;
        top_gmax = rgr->gmax;
	top_h = rgr->h;
	top_L = rgr->L;
	if (first)
	{
	    first = NO;
	    size = 1;
	    for (i = 0; i < dim; ++i)
	    	size *= (top_gmax[i] + 1);
	    switch(dim)
	    {
	    case 1:
		imin = (lbuf[0] == 0) ? 1 : lbuf[0];
		imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
		break;
	    case 2:
		imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    	jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
		imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    	jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
		break;
	    case 3:
		imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    	jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
		kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
		imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    	jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
		kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
		break;
	    }
	}

}

void WENO_FLUX::resetAdvection()
{
	int i;
	for (i = 0; i < size; i++)
	    adv[i] = 0.0;
}

void WENO_FLUX::computeAdvection()
{
	int i,l;
	setSolverDomain();
	resetAdvection();
	for (l = 0; l < dim; l++)
	   addFluxInDirection(l,U,T,adv);
	FT_ParallelExchGridArrayBuffer(adv,front,NULL);
}

void WENO_FLUX::printField(double* var,const char *var_name)
{
	/*test advection weno*/
        FILE* outfile;
        char filename[200];
        int ic,i,j;
        sprintf(filename,"%s/%s-%d-nd%d",
        OutName(front),var_name,front->step,pp_mynode());
        outfile = fopen(filename,"w");
        for (j = jmin; j <= jmax; ++j)
        for (i = imin; i <= imax; ++i)
        {
                ic = d_index2d(i,j,top_gmax);
                fprintf(outfile,"%15.14f\n",var[ic]);
        }
        fclose(outfile);
}

void WENO_FLUX::addFluxInDirection(
	int dir,
	double** U, /*wave speed*/
	double*  T, /*temperature*/
	double* flux) /*flux field*/
{
    switch(dim)
    {
	case 1:
	    return addFluxInDirection1d(dir,U,T,flux);
	case 2:
	    return addFluxInDirection2d(dir,U,T,flux);
	case 3:
	    return addFluxInDirection3d(dir,U,T,flux);
	default:
	    printf("Unknown dim %d\n",dim);
	    clean_up(ERROR);
    }
}


void WENO_FLUX::addFluxInDirection1d(
	int dir,
	double **U,
	double *T,
	double *adv)
{
	printf("1d was not implemented!\n");
	return;
}

void WENO_FLUX::addFluxInDirection2d(
	int dir,
	double **U,
	double *T,
	double *adv)
{
	static double *vU;
	static double *vT;
	static double *vflux;
	static boolean first = YES;
	double lambda;
	int icoords[MAXD], n, i, j ,index, array_size;
	
	lambda = -1.0/top_h[dir];

	if (first)
	{
    	    first = NO;
	    array_size = 1;
	    for (i = 0; i < dim; i++)
	    if (array_size < top_gmax[i]+2*nrad+1)
		array_size = top_gmax[i]+2*nrad+1;
	    FT_VectorMemoryAlloc((POINTER*)&vU,array_size,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&vT,array_size,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&vflux,array_size,sizeof(double));
	}
	switch (dir)
	{
	case 0:
	    for (j = jmin; j <= jmax; j++)
	    {
		n = 0;
	        for (i = imin; i <= imax; i++)
	        {
		    index = d_index2d(i,j,top_gmax);
		    vU[n+nrad] = U[dir][index];
		    vT[n+nrad] = T[index];
	            n++;
		}
		icoords[1] = j;
		icoords[0] = imin;
	        appendGhostBuffer(T,vT,n,icoords,0,0);

		icoords[0] = imax;
	        appendGhostBuffer(T,vT,n,icoords,0,1);
	  	Weno5_Get_Flux(vT,vU,vflux,lambda,n);	

		n = 0;
	        for (i = imin; i <= imax; i++)
		{
		    index = d_index2d(i,j,top_gmax);
		    adv[index] += vflux[n+nrad];
		    n++;
		}
	    }
	    break;
	case 1:
	    for (i = imin; i <= imax; i++)
	    {
		n = 0;
	        for (j = jmin; j <= jmax; j++)
	        {
		    index = d_index2d(i,j,top_gmax);
		    vU[n+nrad] = U[dir][index];
		    vT[n+nrad] = T[index];
	            n++;
		}
		icoords[0] = i;
		icoords[1] = jmin;
	        appendGhostBuffer(T,vT,n,icoords,1,0);

		icoords[1] = jmax;
	        appendGhostBuffer(T,vT,n,icoords,1,1);
	  	Weno5_Get_Flux(vT,vU,vflux,lambda,n);	
		n = 0;
	        for (j = jmin; j <= jmax; j++)
		{
		    index = d_index2d(i,j,top_gmax);
		    adv[index] += vflux[n+nrad];
		    n++;
		}
	    }
	    break;
	}
}

void WENO_FLUX::addFluxInDirection3d(
	int dir,
	double **U,
	double *T,
	double *adv)
{
	static double *vU;
	static double *vT;
	static double *vflux;
	static boolean first = YES;
	int icoords[MAXD], n, i, j, k, index, array_size;
	double lambda;

	lambda = -1.0/top_h[dir];
	if (first)
	{
    	    first = NO;
	    array_size = 1;
	    for (i = 0; i < dim; i++)
	    if (array_size < top_gmax[i]+2*nrad)
		array_size = top_gmax[i]+2*nrad+1;
	    FT_VectorMemoryAlloc((POINTER*)&vU,array_size,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&vT,array_size,sizeof(double));
	    FT_VectorMemoryAlloc((POINTER*)&vflux,array_size,sizeof(double));
	}
	n = 0;
	switch (dir)
	{
	case 0:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    {
		n = 0;
	        for (i = imin; i <= imax; i++)
	        {
		    index = d_index3d(i,j,k,top_gmax);
		    vU[n+nrad] = U[dir][index];
		    vT[n+nrad] = T[index];
	            n++;
		}
		icoords[1] = j;
		icoords[2] = k;
		icoords[0] = imin;
	        appendGhostBuffer(T,vT,n,icoords,0,0);

		icoords[0] = imax;
	        appendGhostBuffer(T,vT,n,icoords,0,1);
	  	Weno5_Get_Flux(vT,vU,vflux,lambda,n);	
		n = 0;
	        for (i = imin; i <= imax; i++)
		{
		    index = d_index3d(i,j,k,top_gmax);
		    adv[index] += vflux[n+nrad];
		    n++;
		}
	    }
	    break;
	case 1:
	    for (k = kmin; k <= kmax; k++)
	    for (i = imin; i <= imax; i++)
	    {
		n = 0;
	        for (j = jmin; j <= jmax; j++)
	        {
		    index = d_index3d(i,j,k,top_gmax);
		    vU[n+nrad] = U[dir][index];
		    vT[n+nrad] = T[index];
	            n++;
		}
		icoords[0] = i;
		icoords[2] = k;
		icoords[1] = jmin;
	        appendGhostBuffer(T,vT,n,icoords,1,0);

		icoords[1] = jmax;
	        appendGhostBuffer(T,vT,n,icoords,1,1);
	  	Weno5_Get_Flux(vT,vU,vflux,lambda,n);	
		n = 0;
	        for (j = jmin; j <= jmax; j++)
		{
		    index = d_index3d(i,j,k,top_gmax);
		    adv[index] += vflux[n+nrad];
		    n++;
		}
	    }
	    break;
	case 2:
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		n = 0;
	        for (k = kmin; k <= kmax; k++)
	        {
		    index = d_index3d(i,j,k,top_gmax);
		    vU[n+nrad] = U[dir][index];
		    vT[n+nrad] = T[index];
	            n++;
		}
		icoords[0] = i;
		icoords[1] = j;
		icoords[2] = kmin;
	        appendGhostBuffer(T,vT,n,icoords,1,0);

		icoords[2] = kmax;
	        appendGhostBuffer(T,vT,n,icoords,1,1);
	  	Weno5_Get_Flux(vT,vU,vflux,lambda,n);	
		n = 0;
	        for (k = kmin; k <= kmax; k++)
		{
		    index = d_index3d(i,j,k,top_gmax);
		    adv[index] += vflux[n+nrad];
		    n++;
		}
	    }
	    break;
	}
}
void WENO_FLUX::appendGhostBuffer(
	double *T,
	double *vT,
	int n,
	int *icoords,
	int idir,
	int nb)
{
	int i, index;
	HYPER_SURF *hs;
	double crx_coords[MAXD];
	int ic[MAXD];
	POINTER state;
	int is_crxing;
	COMPONENT comp;
	GRID_DIRECTION ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION rdir[3] = {EAST,NORTH,UPPER};

	index = d_index(icoords,top_gmax,dim);
	comp = top_comp[index];


	for (i = 0; i < dim; i++)
		ic[i] = icoords[i];
	switch (nb)
	{
	case 0:
	    for (i = 1; i <= nrad; i++)
	    {
		ic[idir] = icoords[idir] - i + 1;
		is_crxing = (*findStateAtCrossing)(front,ic,ldir[idir],comp,
                                        &state,&hs,crx_coords);
		if (is_crxing == NO_PDE_BOUNDARY)
		{
		    ic[idir]--;
		    index = d_index(ic,top_gmax,dim);
		    vT[nrad-i] = T[index];
		}
		else
		{
		    printf("BOUNDARY %d is not implemented\n",wave_type(hs));
		    clean_up(ERROR);
		}
	    }
	    break;
	case 1:
	    for (i = 1; i <= nrad; ++i)
	    {
		ic[idir] = icoords[idir] + i -1;
		is_crxing = (*findStateAtCrossing)(front,ic,rdir[idir],comp,
                                        &state,&hs,crx_coords);
                if (is_crxing == NO_PDE_BOUNDARY)
                {
                    ic[idir]++;
                    index = d_index(ic,top_gmax,dim);
                    vT[n+nrad+i-1] = T[index];
                }
                else
                {
		    printf("BOUNDARY %d is not implemented\n",wave_type(hs));
		    clean_up(ERROR);
                }
	    }
	    break;
	}
}



