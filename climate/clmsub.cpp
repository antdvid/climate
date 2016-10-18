/************************************************************************************
FronTier is a set of libraries that implements differnt types of Front Traking algorithms.
Front Tracking is a numerical method for the solution of partial differential equations 
whose solutions have discontinuities.  


Copyright (C) 1999 by The University at Stony Brook. 
 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

******************************************************************************/


#include <iFluid.h>
#include "climate.h"
#include <time.h>
        /*  Function Declarations */
static void zero_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static void rand_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static void Taylor_Green_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static void ABC_state(COMPONENT,double*,IF_FIELD*,int,int,IF_PARAMS*);
static void Fourier_state(COMPONENT,int*,double*,double*,double*,
			RECT_GRID*,IF_PARAMS*);
static void Rogallo_state(COMPONENT,int*,double*,double*,double*,
			RECT_GRID*,IF_PARAMS*);
static double intrp_between(double,double,double,double,double);
static double (*getStateVel[MAXD])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};

extern void melt_flowThroughBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
        Tan_stencil **tsten;
        Nor_stencil *nsten;
        FLOW_THROUGH_PARAMS *ft_params = (FLOW_THROUGH_PARAMS*)params;
        POINT *oldp = ft_params->oldp;
        COMPONENT comp = ft_params->comp;
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        PARAMS *eqn_params = (PARAMS*)front->extra2;
        IF_FIELD *iF_field = iFparams->field;
        PHASE_FIELD *pH_field = eqn_params->field;
        double dir[MAXD];
        double u[3];            /* velocity in the sweeping direction */
        double v[3][MAXD];      /* velocity in the orthogonal direction */
        double vort[3];         /* vorticity stencil */
        double pres[3];         /* pressure stencil */
        double temp[3];         /* temperature stencil */
        double f_u;             /* u flux in the sweeping direction */
        double f_v[MAXD];       /* v flux in the orthogonal direction */
        double f_vort;          /* vort flux */
        double f_pres;          /* pressure flux */
        double f_temp;          /* temperature flux */
        double dn,dt = front->dt;
        STATE *newst = (STATE*)state;
        STATE  **sts;
        int i,j,dim = front->rect_grid->dim;
        int nrad = 2;

        if (debugging("flow_through"))
            printf("Entering melt_flowThroughBoundaryState()\n");
	
	tsten = FrontGetTanStencils(front,oldp,nrad);
        for (i = 0; i < dim; ++i)
            dir[i] = tsten[0]->dir[i];
        dn = FT_GridSizeInDir(dir,front);

        if (comp == negative_component(hs))
            sts = (STATE**)tsten[0]->leftst;
        else
            sts = (STATE**)tsten[0]->rightst;

        if (debugging("flow_through"))
        {
            printf("Ambient component: %d\n",comp);
            printf("hs = %p  oldp->hs = %p\n",(POINTER)hs,(POINTER)oldp->hs);
            printf("Time step = %f  Tangential grid size = %f\n",dt,dn);
            printf("Tangential direction: ");
            for (j = 0; j < dim; ++j)
                printf("%f ",tsten[0]->dir[j]);
            printf("\n");
            printf("Tan_stencil at point p(%f %f)\n",Coords(oldp)[0],
                                Coords(oldp)[1]);
            printf("Left points:\n");
            for (i = 0; i < nrad; ++i)
            {
                for (j = 0; j < dim; ++j)
                    printf("%f ",Coords(tsten[0]->p[-i])[j]);
                printf("\n");
            }
            printf("Right points:\n");
            for (i = 0; i < nrad; ++i)
            {
                for (j = 0; j < dim; ++j)
                    printf("%f ",Coords(tsten[0]->p[i])[j]);
                printf("\n");
            }
        }

	for (j = 0; j < 3; ++j)
            u[j] = 0.0;
        for (j = 0; j < 3; ++j)
        {
            vort[j] = sts[j-1]->vort;
            pres[j] = sts[j-1]->pres;
            temp[j] = sts[j-1]->temperature;
            for (i = 0; i < dim; ++i)
            {
                u[j] += sts[j-1]->vel[i]*dir[i];
                v[j][i] = sts[j-1]->vel[i]*(1.0 - dir[i]);
            }
        }

        f_u = burger_flux(u[0],u[1],u[2]);
        for (i = 0; i < dim; ++i)
            f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
        f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
        f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
        f_temp = linear_flux(u[1],temp[0],temp[1],temp[2]);

        for (i = 0; i < dim; ++i)
            newst->vel[i] = sts[0]->vel[i] - dt/dn*(
                f_u*dir[i] + f_v[i]) ;
        newst->vort = sts[0]->vort - dt/dn*f_vort;
        newst->pres = sts[0]->pres - dt/dn*f_pres;
        newst->temperature = sts[0]->temperature - dt/dn*f_temp;

        nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);
        for (i = 0; i < dim; ++i)
            dir[i] = nsten->nor[i];
        dn = FT_GridSizeInDir(dir,front);

        if (debugging("flow_through"))
        {
            printf("Time step = %f  Normal grid size = %f\n",dt,dn);
            printf("Normal direction: ");
            for (j = 0; j < dim; ++j)
                printf("%f ",nsten->nor[j]);
            printf("\n");
            printf("Nor_stencil at point p(%f %f)\n",Coords(oldp)[0],
                                Coords(oldp)[1]);
            printf("Nor_stencil:\n");
            for (i = 0; i < nrad; ++i)
            {
                for (j = 0; j < dim; ++j)
                    printf("%f ",nsten->pts[i][j]);
                printf("\n");
            }
        }
	
	for (j = 0; j < 3; ++j)
            u[j] = 0.0;
        for (j = 0; j < 2; ++j)
        {
            vort[j] = sts[0]->vort;
            pres[j] = sts[0]->pres;
            temp[j] = sts[0]->temperature;
            for (i = 0; i < dim; ++i)
            {
                u[j] += sts[0]->vel[i]*dir[i];
                v[j][i] = sts[0]->vel[i]*(1.0 - dir[i]);
            }
        }
        for (i = 0; i < dim; ++i)
        {
            double vtmp;
            FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
                        iF_field->vel[i],getStateVel[i],&vtmp,&sts[0]->vel[i]);
            u[2] += vtmp*dir[i];
            v[2][i] = vtmp*(1.0 - dir[i]);
        }
        if (dim == 2)
        {
            FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
                        iF_field->vort,getStateVort,&vort[2],&sts[1]->vort);
        }
        FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],iF_field->pres,
                            getStatePres,&pres[2],&sts[1]->pres);
        FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],pH_field->temperature,
                            getStateTemperature,&temp[2],&sts[1]->temperature);

        f_u = burger_flux(u[0],u[1],u[2]);
        for (i = 0; i < dim; ++i)
            f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
        f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
        f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
        f_temp = linear_flux(u[1],temp[0],temp[1],temp[2]);

        for (i = 0; i < dim; ++i)
            newst->vel[i] += - dt/dn*(f_u*dir[i] + f_v[i]) ;
        newst->vort += - dt/dn*f_vort;
        newst->pres += - dt/dn*f_pres;
        newst->temperature += - dt/dn*f_temp;
	
	if (debugging("flow_through"))
        {
            printf("flow through boundary state:\n");
            print_general_vector("Velocity: ",newst->vel,dim,"\n");
            printf("Pressure: %f\n",newst->pres);
            printf("Vorticity: %f\n",newst->vort);
            printf("Temperature: %f\n",newst->temperature);
        }
}       /*end melt_flowThroughBoundaryState */

static void vFourier_state(
            COMPONENT* c,
            double* coords,
            PHASE_FIELD* field,
            int index,
            int dim,
            PARAMS* eqn_params)
{
	double qs, qe, qmax;
	qs = eqn_params->qs; qe = eqn_params->qe;
        qmax = qs * 1.02;
	if (field->vel[0][index] > 0)
	    field->vapor[index] = qmax;
	else
	    field->vapor[index] = qe;
}

static void lr_state(
	    COMPONENT* c,
	    double* coords,
	    PHASE_FIELD* field,
	    int index,
	    int dim,
	    PARAMS* eqn_params)
{
	double qs, qe, qmax, q0;
	double *L, *U, width, frac, x0;
	double T0,A;
	RECT_GRID* global_grid = eqn_params->global_grid;
	L = global_grid->L; U = global_grid->U;
	frac = eqn_params->frac;
	width = U[0] - L[0];
	x0 = L[0] + 0.5*width;

	qe = eqn_params->qe;
	qs = eqn_params->qs;
	qmax = qs * 1.02;
	if (eqn_params->frac > 0 && eqn_params->frac < 1)
	{
	    if (coords[0] > x0 - width*frac*0.5 &&
	        coords[0] < x0 + width*frac*0.5)
	        field->vapor[index] = qmax;
	    else
	        field->vapor[index] = qe;
	}
	else
	{
	    A = pow(0.512,8.0)*1.0e6;
	    field->vapor[index] = (qmax - qe)
			      * exp(-A*pow(coords[0]/width-0.5,8.0))+qe;
	}
}

static void tb_state(
	    COMPONENT* c,
	    double* coords,
	    PHASE_FIELD* field,
	    int index,
	    int dim,
	    PARAMS* eqn_params)
{
	double qs, qe, qmax;
	double x0, A;
	double* L, *U, width, frac;

	RECT_GRID* global_grid = eqn_params->global_grid;
	L = global_grid->L; U = global_grid->U;
	width = U[dim-1] - L[dim-1];
	x0 = L[dim-1] + 0.5*width;
	frac = eqn_params->frac;

	qs = eqn_params->qs; 
	qe = eqn_params->qe;
	qmax = qs * 1.02;

	if (eqn_params->frac > 0 && eqn_params->frac < 1)
	{
	    if (coords[dim-1] > x0 - width*frac*0.5 &&
		coords[dim-1] < x0 + width*frac*0.5)
	        field->vapor[index] = qmax;
	    else
	        field->vapor[index] = qe;
	}
	else
	{
	    A = pow(0.512,8.0)*1.0e6;
	    field->vapor[index] = (qmax - qe)
			      * exp(-A*pow(coords[dim-1] - x0,6))+qe;
 	}	
}

static void sine_state(
	    COMPONENT* c,
	    double* coords,
	    PHASE_FIELD* field,
	    int index,
	    int dim,
	    PARAMS* eqn_params)
{
	field->vapor[index] = sin(PI*(coords[0]+coords[1]));
}
	
void init_vapor_state_func(
	Front *front,
        VCARTESIAN *vcartesian)
{
        PARAMS *eqn_params = (PARAMS*)front->extra2;
        switch(eqn_params->init_vapor_state)
        {
	    case TB_STATE:
		vcartesian->getInitialVaporState = tb_state;
	        break;
	    case LR_STATE:
		vcartesian->getInitialVaporState = lr_state;
	        break;
	    case FOURIER_STATE:
		vcartesian->getInitialVaporState = vFourier_state;
		break;
	    case RAND_STATE:
	    case CONST_STATE:
		vcartesian->getInitialVaporState = NULL;
                break;
	    case SINE_STATE:
		vcartesian->getInitialVaporState = sine_state;
		break;
            default:
		printf("State %d is not implemented\n",
			eqn_params->init_vapor_state);
		clean_up(ERROR);
		break;
        }
}	/* end init_vapor_state_func */
   
static void nb_state(
	    COMPONENT* c,
	    double* coords,
	    PHASE_FIELD* field,
	    int index,
	    int dim,
	    PARAMS* eqn_params)
{
	double q0, T0;
	
	T0 = eqn_params->T0; q0 = eqn_params->qv0;

	/*lowest temperature is given, transform T0 to mean value*/
	T0 /= 1 - 0.608*0.001*(1.02*eqn_params->qs - q0);
	field->temperature[index] = 
		T0 - 0.608*T0*0.001*(field->vapor[index] - q0);
}/*neutral buoyancy state*/

void init_temp_state_func(
	Front *front,
        VCARTESIAN *vcartesian)
{
        PARAMS *eqn_params = (PARAMS*)front->extra2;
        switch(eqn_params->init_temp_state)
        {
	    case PRESET_STATE:
		vcartesian->getInitialTempState = nb_state;
	        break;
	    case RAND_STATE:
	    case CONST_STATE:
		vcartesian->getInitialTempState = NULL;
		break;
            default:
		printf("State %d is not implemented\n",
			eqn_params->init_temp_state);
		clean_up(ERROR);
        }
}	/* end init_fluid_state_func */

void init_fluid_state_func(
	Front *front,
        Incompress_Solver_Smooth_Basis *l_cartesian)
{
        PARAMS *eqn_params = (PARAMS*)front->extra2;
        switch(eqn_params->init_state)
        {
            case RAND_STATE:
                l_cartesian->getInitialState = rand_state;
                break;
            case TAYLOR_STATE:
                l_cartesian->getInitialState = Taylor_Green_state;
                break;
	    case ZERO_STATE:
		l_cartesian->getInitialState = zero_state;
	        break;
	    case ABC_STATE:
		l_cartesian->getInitialState = ABC_state;
		break;
	    case FOURIER_STATE:
		l_cartesian->setInitialVelocity = Rogallo_state;
		break;
            default:
                l_cartesian->getInitialState = zero_state;
        }
}	/* end init_fluid_state_func */

static void Rogallo_state(
	COMPONENT comp,
	int *N,
	double *vel_x,
	double *vel_y,
	double *vel_z,
	RECT_GRID *gr,
	IF_PARAMS *iFparams)
{
	if (debugging("trace"))
	    printf("Entering Fourier_state()\n");
	int gmax[MAXD],icrds[MAXD]; 
	int i, j, k, l, ll, ld, Nr, index, dim;
	double wn, w0, L[MAXD], w[MAXD];
	double phi, theta, E;
	double u0 = iFparams->Urms;
	if (u0 == 0.0)
	{
	    printf("initial Urms is not given!\n");
	    printf("Enter initial velocity(m/s):\n");
	    clean_up(ERROR);
	}
	fftw_complex *U, *V, *W, *Div, alpha, beta;
	dim = gr->dim;
	Nr = 1; 
	for (i = 0; i < dim; i++)
	{
	    L[i] = gr->GU[i]-gr->GL[i];
	    gmax[i] = N[i];
	    Nr *= N[i];
	}
	w0 = 2.449489743;
	gmax[dim-1] = N[dim-1]/2 + 1; 
	U = new fftw_complex[Nr];
	V = new fftw_complex[Nr];
	if (dim == 3)
	    W = new fftw_complex[Nr];
	Div = new fftw_complex[Nr];
	

	switch (dim)
	{
	    case 2:
		for (i = 0; i < Nr; i++)
		{
	   	    U[i][0] = U[i][1] = 0.0;
	   	    V[i][0] = V[i][1] = 0.0;
		}
		for (i = 0; i < gmax[0]; i++)
		for (j = 0; j < gmax[1]; j++)
		{
		    if (i == 0 && j ==0)
		        continue;
		    index = j + gmax[1] * i; 
		    /*row major format*/
		    icrds[0] = i; icrds[1] = j;
		    for (l = 0; l < dim; l++)
		    {
		        if (icrds[l] > N[l]/2)
		            icrds[l] = (icrds[l]-N[l]);
		    }
		    wn = 0.0;
		    for (l = 0; l < dim; l++)
		    {
			w[l] = (icrds[l]);
			wn += w[l]*w[l];
		    }
		    wn = sqrt(wn);
		    
		    theta  = (double)rand() / (RAND_MAX + 1.0) *(2*M_PI);
		    E = 16.0/sqrt(0.5*M_PI)*u0*u0*pow(wn,4)/pow(w0,5)
			* exp(-2.0*wn*wn/(w0*w0)); 
		    alpha[0] = sqrt(E/(4.0*M_PI*wn*wn))*cos(theta);
		    alpha[1] = sqrt(E/(4.0*M_PI*wn*wn))*sin(theta);
		    U[index][0] = w[1]/wn*alpha[0];
		    U[index][1] = w[1]/wn*alpha[1];
		    V[index][0] = -w[0]/wn*alpha[0];
                    V[index][1] = -w[0]/wn*alpha[1];
		    Div[index][0] = U[index][0]*w[0] + V[index][0]*w[1];
		    Div[index][1] = U[index][0];
		}
		/*iFFT*/
		fftnd(U,dim,N,-1);
		fftnd(V,dim,N,-1);
	   	/*get solution from iFFT*/
		for (i = 0; i < Nr; i++)
		{
	    	    vel_x[i] = U[i][0];
	    	    vel_y[i] = V[i][0];
		}
		delete[] U;
		delete[] V;
	 	break;
	    case 3:
		for (i = 0; i < Nr; i++)
		{
	   	    U[i][0] = U[i][1] = 0.0;
	   	    V[i][0] = V[i][1] = 0.0;
	   	    W[i][0] = W[i][1] = 0.0;
		}
                for (k = 0; k < gmax[2]; k++)
                for (j = 0; j < gmax[1]; j++)
		for (i = 0; i < gmax[0]; i++)
                {
		    if (i == 0 && j == 0 && k == 0)
			continue;
                    index = k + gmax[2]*(j + gmax[1]*i); /*row major format*/ 
		    icrds[0] = i; 
		    icrds[1] = j;
		    icrds[2] = k;
		    /*map indices to wavenumber*/
		    for (l = 0; l < dim; l++)
		    {
		        if (icrds[l] > N[l]/2)
		            icrds[l] = (icrds[l]-N[l]);
		    }
		    wn = 0.0;
		    for (l = 0; l < dim; l++)
		    {
			w[l] = (icrds[l]);
			wn += w[l]*w[l];
		    }
		    wn = sqrt(wn);
		    
		    phi  = (double)rand() / (RAND_MAX + 1.0) * (2*M_PI);
		    E = 16.0/sqrt(0.5*M_PI)*u0*u0*pow(wn,4)/pow(w0,5)
			* exp(-2.0*wn*wn/(w0*w0)); 
		    theta  = (double)rand() / (RAND_MAX + 1.0) *(2*M_PI);
		    alpha[0] = sqrt(E/(4.0*M_PI*wn*wn))*cos(theta)*cos(phi);
		    alpha[1] = sqrt(E/(4.0*M_PI*wn*wn))*sin(theta)*cos(phi);
		    /*use different theta for alpha and beta*/
		    theta  = (double)rand() / (RAND_MAX + 1.0) *(2*M_PI);
		    beta[0] = sqrt(E/(4.0*M_PI*wn*wn))*cos(theta)*sin(phi);
		    beta[1] = sqrt(E/(4.0*M_PI*wn*wn))*sin(theta)*sin(phi);
		    if (w[0] == 0 && w[1] == 0)
			continue;
		    U[index][0] = (alpha[0]*wn*w[1]+beta[0]*w[0]*w[2])
					/ (wn*sqrt(w[0]*w[0]+w[1]*w[1]));
		    U[index][1] = (alpha[1]*wn*w[1]+beta[1]*w[0]*w[2])
                                        / (wn*sqrt(w[0]*w[0]+w[1]*w[1]));
	 	    V[index][0] = (beta[0]*w[1]*w[2]-alpha[0]*wn*w[0])
                                        / (wn*sqrt(w[0]*w[0]+w[1]*w[1]));
                    V[index][1] = (beta[1]*w[1]*w[2]-alpha[0]*wn*w[0])
                                        / (wn*sqrt(w[0]*w[0]+w[1]*w[1]));
		    W[index][0] = -(beta[0]*sqrt(w[0]*w[0]+w[1]*w[1]))
					/ wn;
		    W[index][1] = -(beta[1]*sqrt(w[0]*w[0]+w[1]*w[1]))
                                        / wn;
		    Div[index][0] = U[index][0]*w[0] 
				  + V[index][0]*w[1] + W[index][0]*w[2];
		    Div[index][1] = U[index][0];
                }
		/*iFFT*/
		fftnd(U,dim,N,-1);
		fftnd(V,dim,N,-1);
	    	fftnd(W,dim,N,-1);
		/*get solution from iFFT*/
		for (i = 0; i < Nr; i++)
		{
	    	    vel_x[i] = U[i][0];
	    	    vel_y[i] = V[i][0];
	      	    vel_z[i] = W[i][0];
		}
		delete[] U;
		delete[] V;
	    	delete[] W;
                break;
	    default:
		printf("Dim can only be 2 and 3, unknown dim = %d!\n",dim);
		clean_up(ERROR);
	}

	if (debugging("init_vel"))
	{
	    FILE* file = fopen("init_vel","w");
	    for (i = 0; i < Nr; i++)
	    {
		if (dim == 3)
	    	  fprintf(file,"%9.8f %9.8f %9.8f\n",vel_x[i],vel_y[i],vel_z[i]);
		else if (dim == 2)
	    	  fprintf(file,"%9.8f %9.8f\n",vel_x[i],vel_y[i]);
	    }
	    fclose(file);
	}
	if (debugging("init_vel"))
	{
	    FILE* file = fopen("FFT_div","w");
	    for (i = 0; i < Nr; i++)
	    	  fprintf(file,"%9.8f %9.8f\n",Div[i][0],Div[i][1]);
	    fclose(file);
	}

	if (debugging("trace"))
	    printf("Leaving Fourier_state()\n");
}
/*end Rogallo_state*/

/*Fourier_state*/
/*only used in setParallelVelocity, no for-loops needed*/
/*using FFTW library, real to complex and complex to real functions*/
static void Fourier_state(
	COMPONENT comp,
	int *N,
	double *vel_x,
	double *vel_y,
	double *vel_z,
	RECT_GRID *gr,
	IF_PARAMS *iFparams)
{
	if (debugging("trace"))
	    printf("Entering Fourier_state()\n");
	int gmax[MAXD], i, j, k, Nr, dim, index;
	fftw_complex *U;
	double wn, phi, U_max = 0.0, L[MAXD];

	dim = gr->dim;
	Nr = 1;
	for (i = 0; i < dim; i++)
	{
	    L[i] = gr->GU[i]-gr->GL[i];
	    gmax[i] = N[i];
	    Nr *= N[i];
	}
	printf("grid length = [%f %f]\n",L[0],L[1]);
	gmax[dim-1] = N[dim-1]/2 + 1; 
	U = new fftw_complex[Nr];
	
	switch (dim)
	{
	    case 2:
		for (i = 0; i < gmax[0]; i++)
		for (j = 0; j < gmax[1]; j++)
		{
		    index = j + gmax[1] * i; /*row major format*/
		    if ((i*i + j*j) > 2)
		    {
			U[index][0] = U[index][1] = 0.0;
		 	continue;
		    }
		    wn = (2*M_PI/L[0])*sqrt(i*i+j*j);
		    phi  = (double)rand() / (RAND_MAX + 1.0);
		    U[index][0] = wn*wn*exp(-wn*wn/pow(2*M_PI*4.7568/L[0],2))
                                 * cos(2*M_PI*phi);
                    U[index][1] = wn*wn*exp(-wn*wn/pow(2*M_PI*4.7568/L[0],2))
                                 * sin(2*M_PI*phi);
		}
	 	break;
	    case 3:
                for (k = 0; k < gmax[2]; k++)
                for (j = 0; j < gmax[1]; j++)
		for (i = 0; i < gmax[0]; i++)
                {
                    index = k + gmax[2]*(j + gmax[1]*i); /*row major format*/ 
                    if ((i*i + j*j + k*k) > 6)
                    {
                        U[index][0] = U[index][1] = 0.0;
                        continue;
                    }
                    wn = (2*M_PI/L[0])*sqrt(i*i+j*j+k*k);
                    phi  = (double)rand() / (RAND_MAX + 1.0);
                    U[index][0] = wn*wn*exp(-wn*wn/pow(2*M_PI*4.7568/L[0],2))
                                 * cos(2*M_PI*phi);
                    U[index][1] = wn*wn*exp(-wn*wn/pow(2*M_PI*4.7568/L[0],2))
                                 * sin(2*M_PI*phi);
                }
                break;
	    default:
		printf("Dim can only be 2 and 3, unknown dim = %d!\n",dim);
		clean_up(ERROR);
	}
	fftnd(U,dim,N,-1);
	for (i = 0; i < Nr; i++)
	   U_max = FT_Max(U_max,fabs(U[i][0]));

	for (i = 0; i < Nr; i++)
	{
	    vel_y[i] = vel_x[i] = U[i][0]/U_max * iFparams->Urms;
	    if (dim == 3)
	        vel_z[i] = U[i][0]/U_max * iFparams->Urms;
	}
	delete[] U;
	if (debugging("trace"))
	    printf("Leaving Fourier_state()\n");
}
/*end Fourier_state*/

static void ABC_state(
        COMPONENT comp,
        double *coords,
        IF_FIELD *field,
        int index,
        int dim,
        IF_PARAMS *iFparams)
{
	int i;
	double tcoords[MAXD], a[MAXD] = {2*PI, 2*PI, 2*PI}, A[MAXD] = {1., -1., 0.};
	double **vel = field->vel;
	switch (dim)
	{
	    case 2:
		for (i = 0; i < dim; i++)
	    	    tcoords[i] = coords[i];
	        vel[0][index] = cos(tcoords[1]);
	        vel[1][index] = sin(tcoords[0]);
		break;
	    case 3:
		for (i = 0; i < dim; i++)
	    	    tcoords[i] = coords[i];
		vel[0][index] = cos(tcoords[1])+sin(tcoords[2]);
		vel[1][index] = sin(tcoords[0])+cos(tcoords[2]);
		vel[2][index] = cos(tcoords[0])+sin(tcoords[1]);
		break;
	    default:
		printf("Unknown dim = %d\n",dim);
		clean_up(ERROR);
	}
}       /* end ABC_state */

static void Taylor_Green_state(
        COMPONENT comp,
        double *coords,
        IF_FIELD *field,
        int index,
        int dim,
        IF_PARAMS *iFparams)
{
	int i;
	double tcoords[MAXD], a[MAXD] = {2*PI, 2*PI, 2*PI}, A[MAXD] = {1., -1., 0.};
	double **vel = field->vel;
	switch (dim)
	{
	    case 2:
		for (i = 0; i < dim; i++)
	    	    tcoords[i] = coords[i];
	        vel[0][index] = A[0]*sin(tcoords[0]) * cos(tcoords[1]);
	        vel[1][index] = A[1]*cos(tcoords[0]) * sin(tcoords[1]);
		break;
	    case 3:
		for (i = 0; i < dim; i++)
	    	    tcoords[i] = a[i] * coords[i];
		vel[0][index] = A[0]*cos(tcoords[0])*sin(tcoords[1])*sin(tcoords[2]);
		vel[1][index] = A[1]*sin(tcoords[0])*cos(tcoords[1])*sin(tcoords[2]);
		vel[2][index] = A[2]*sin(tcoords[0])*sin(tcoords[1])*cos(tcoords[2]);
		break;
	    default:
		printf("Unknown dim = %d\n",dim);
		clean_up(ERROR);
	}
}       /* end Taylor_Green_state */

static void rand_state(
        COMPONENT comp,
        double *coords,
        IF_FIELD *field,
        int index,
        int dim,
        IF_PARAMS *iFparams)
{
	short unsigned int seed[3] = {time(NULL)-index,
				      time(NULL),
				      time(NULL)+index};
	//short unsigned int seed[3] = {index,index-10,index+10};
	GAUSS_PARAMS gauss_params;
	double r_bar = 0;
	double sigma = 0.5;
        int i;
	gauss_params.mu = r_bar;
	gauss_params.sigma = sigma;
        for (i = 0; i < dim; ++i)
            field->vel[i][index] = gauss_center_limit((POINTER)&gauss_params,seed);
}       /* end rand_state */

static void zero_state(
        COMPONENT comp,
        double *coords,
        IF_FIELD *field,
        int index,
        int dim,
        IF_PARAMS *iFparams)
{
        int i;
        for (i = 0; i < dim; ++i)
            field->vel[i][index] = 0.0;
        field->pres[index] = iFparams->ref_pres;
}       /* end zero_state */

static double intrp_between(
        double x1,
        double x2,
        double x,
        double y1,
        double y2)
{
        double y;
        if (x1 == x2) return y1;
        y = y1 + (y2 - y1)/(x2 - x1)*(x - x1);
        return y;
}

LOCAL   void Print_particle_array(PARTICLE* particle_array,int num)
{
	int i;
	for (i = 0; i < num; i++)
	    printf("drop[%d]: center=[%f %f %f] flag = %d\n",
		i,particle_array[i].center[0],particle_array[i].center[1],
		particle_array[i].center[2],particle_array[i].flag);
}

LOCAL   PARTICLE* cut_buf_particle(
	PARTICLE* particle_array,
	Front *front,
	int &buf_size,
	int dir,
	int side)
{
	int 	  i,j;
        PARAMS*   eqn_params = (PARAMS*)front->extra2;
        int       num_drops = eqn_params->num_drops;
	double    *ZL,*ZU,*GL,*GU, T;
	static    PARTICLE* buf_particle;
	static    int max_size = 0;
	RECT_GRID *G_gr = &(front->pp_grid->Global_grid);
	RECT_GRID *Z_gr = &(front->pp_grid->Zoom_grid);
	GL = G_gr->L; GU = G_gr->U;
	ZL = Z_gr->L; ZU = Z_gr->U;
	/*count number of particles in buffer*/
	buf_size = 0;
	for (i = 0; i < eqn_params->num_drops; i++)
	{
                if (side == 0 && particle_array[i].center[dir] < ZL[dir])
                     buf_size++;     
                if (side == 1 && particle_array[i].center[dir] > ZU[dir])
                     buf_size++;
	}
	/*allocate memory for buf_particle*/
	if (buf_size > max_size)
	{
		max_size = buf_size;
		free_these(1,buf_particle);
		FT_VectorMemoryAlloc((POINTER*)&buf_particle,max_size,sizeof(PARTICLE));
	}

	j = 0;
	for (i = 0; i < eqn_params->num_drops; i++)
	{
		if ((side == 0 && particle_array[i].center[dir] < ZL[dir]) ||
		    (side == 1 && particle_array[i].center[dir] > ZU[dir]))
		{
		    /*modify coords if position out of global domain*/
		    if (side == 0 && particle_array[i].center[dir] < GL[dir])
		        particle_array[i].center[dir] 
			= GU[dir]+fmod(particle_array[i].center[dir]-GL[dir],
			  GU[dir]-GL[dir]);
		    if (side == 1 && particle_array[i].center[dir] > GU[dir])
		        particle_array[i].center[dir] 
			= GL[dir]+fmod(particle_array[i].center[dir]-GU[dir],
			  GU[dir]-GL[dir]);

		    buf_particle[j] = particle_array[i];
		    particle_array[i].flag = NO;
		    j++;
	      }
	}
	/*IMPORTANT: please donot modify number of drops in domain*/
	/*only set the flag to NO*/
	return buf_particle;
}

LOCAL  void send_particle(PARTICLE* buf_particle,int buf_size,int dst_id)
{
	pp_send(chunk_id(dst_id),(POINTER)&buf_size,sizeof(int),dst_id);
	if (debugging("particles"))
	{
	    printf("In send_particle(),buf_size = %d\n",buf_size);
	    printf("From pp_id[%d] to pp_id[%d]\n",pp_mynode(),dst_id);
	}
	pp_send(array_id(dst_id),(POINTER)buf_particle,
		sizeof(PARTICLE)*buf_size,dst_id);
	return;
}

LOCAL  PARTICLE *receive_particle(int src_id,int &buf_size)
{
	int myid = pp_mynode(), i;
	static int max_buf_size = 0;
	static PARTICLE *adj_particle;

	pp_recv(chunk_id(myid),src_id,(POINTER)&buf_size,sizeof(int));
	if (debugging("particles"))
	{
	    printf("In receive_particle(),buf_size = %d\n",buf_size);
	    printf("From pp_id[%d] to pp_id[%d]\n",src_id,myid);
	}
	if (buf_size > max_buf_size)
	{
	    max_buf_size = buf_size;
	    free(adj_particle);
            FT_VectorMemoryAlloc((POINTER*)&adj_particle,buf_size,sizeof(PARTICLE));	
	}
	pp_recv(array_id(myid),src_id,(POINTER)adj_particle,sizeof(PARTICLE)*buf_size);

	return adj_particle;
}

LOCAL   void merge_particle(PARTICLE** particle_array,
			   PARTICLE* adj_particle,
			   int adj_size,
			   PARAMS* params)
{
	static int max_size = 0;
	int i,j,k,flag,new_size,count_flag;
	int num_drops = params->num_drops;
	static PARTICLE *temp_particle_array;
	new_size = num_drops + adj_size;
	flag = 0;
	if (new_size > max_size)
	{
	    max_size = new_size;
	    free_these(1,temp_particle_array);
	    FT_VectorMemoryAlloc((POINTER*)&temp_particle_array,max_size,sizeof(PARTICLE));
	}
	i = 0;
	count_flag = 0;
	for (k = 0; k < num_drops; k++)	
	{
	    if ((*particle_array)[k].flag == NO)
	    {
		count_flag++;
		continue;
	    }
	    temp_particle_array[i] = (*particle_array)[k];
	    i++;
	}
	new_size -= count_flag;
	for (j = 0; j < adj_size; ++j)
	    temp_particle_array[i+j] = adj_particle[j];

	if (params->num_drops < new_size)
	{
	    free(*particle_array);
            FT_VectorMemoryAlloc((POINTER*)&(*particle_array),new_size,sizeof(PARTICLE));
	}
	for (j = 0; j < new_size; ++j)
	    (*particle_array)[j] = temp_particle_array[j];
	params->num_drops = new_size;
	if (debugging("particles"))
	    printf("After merging %d drops cut, %d drops merged, %d drops contained\n",
		   count_flag,adj_size,new_size);
	return;
}

LOCAL  void ParallelExchParticle(PARTICLE** particle_array,Front *front)
{
	INTERFACE	*intfc = front->interf;
	PP_GRID		*pp_grid = front->pp_grid;
	PARTICLE	*buf_particle, *adj_particle;
	PARAMS		*eqn_params = (PARAMS*)front->extra2;
	int 		*G;
	int 		num_drops = eqn_params->num_drops;
	int 		buf_size, adj_size;
	int 		myid, dst_id;
	int 		me[MAXD], him[MAXD];
	int 		dim = intfc->dim;
	int 		i,j,k;

	myid = pp_mynode();
	G = pp_grid->gmax;
	find_Cartesian_coordinates(myid,pp_grid,me);
	
	
	if (debugging("trace"))
	    printf("Entering ParallelExchParticle()\n");
	for (i = 0; i < dim; i++)
	{
	    for (j = 0; j < 2; j++)
	    {
	    if (debugging("particles"))
	        printf("Exchange particles in direction %d side %d\n",i,j);
		pp_gsync();
		for (k = 0; k < dim; k++)
		    him[k] = me[k];
		
		if (rect_boundary_type(intfc,i,j) == SUBDOMAIN_BOUNDARY)
		{
		    him[i] = me[i] + 2*j -1;
		    him[i] = (him[i]+G[i])%G[i];
		}
		/*Send particles to adjacent domain if necessary*/
		if (rect_boundary_type(intfc,i,j) == SUBDOMAIN_BOUNDARY)
		{
		    dst_id = domain_id(him,G,dim);
		    buf_particle = cut_buf_particle(*particle_array,front,buf_size,i,j); 
		    /*send buffer to corresponding neighbour*/
		    if (me[i] == him[i])
		    {
			adj_particle = buf_particle;
			adj_size = buf_size;
		    }
		    else
		    {
			send_particle(buf_particle,buf_size,dst_id);
		    }
		}
		/*Receive adjacent buffer region if necessary*/
		if (rect_boundary_type(intfc,i,(j+1)%2) == SUBDOMAIN_BOUNDARY)
		{
		    him[i] = me[i] - 2*j + 1;
		    him[i] = (him[i]+G[i])%G[i];
                    dst_id = domain_id(him,G,dim);
                    if (me[i] != him[i])
			adj_particle = receive_particle(dst_id,adj_size);
		}
	    	merge_particle(particle_array,adj_particle,adj_size,eqn_params);	
	    }
	}
	if (debugging("trace"))
            printf("Leaving ParallelExchParticle\n"); 
	return;	
}

extern void ParticlePropagate(Front *fr)
{
	if (debugging("trace"))
	    printf("Entering ParticlePropage()\n");
	start_clock("ParticlePropagate");
        RECT_GRID *gr = FT_GridIntfcTopGrid(fr);
        RECT_GRID *rect_grid = &(fr->pp_grid->Global_grid);
        IF_PARAMS *iFparams = (IF_PARAMS*)fr->extra1;
	PARAMS *eqn_params = (PARAMS*)fr->extra2;
	PARTICLE* particle_array = eqn_params->particle_array;
	double **vel = iFparams->field->vel;
	double *supersat = eqn_params->field->supersat;
	double *gravity = iFparams->gravity;
        int *gmax = FT_GridIntfcTopGmax(fr);
	int i, j, index, dim = gr->dim;
	double T;
	int ic[MAXD];
	double u[MAXD];
	double *center;
	double s; /*restore local supersaturation*/
	double *cvel; /*center velocity for droplets*/
	double a;  /*acceleration*/
	double dt = fr->dt;

        /*computing finite respone time*/
        double rho_0    = iFparams->rho2;/*fluid density*/
        double mu       = iFparams->mu2;/*viscosity*/
	double R, rho, tau_p, delta_R;
	double R_max = 0;
	double R_min = HUGE;
	double w = 2*PI/5.0;
	
	for (i = 0; i < eqn_params->num_drops; i++)
	{
            /*computing finite respone time*/
            R        = particle_array[i].radius;/*droplet radius*/
            rho      = particle_array[i].rho;/*water droplet density*/
            tau_p    = 2 * rho*R*R/(9*rho_0*mu);/*response time*/

	    if (R == 0)
	    {
		R_min = 0;
	        continue;
	    }
	    /*find index at coords*/
	    center = particle_array[i].center;
	    rect_in_which(center,ic,gr);
	    index = d_index(ic,gmax,dim);
	    cvel = particle_array[i].vel;
	    /*compute radius for particle[i]*/
	    s = supersat[index];

	    for (j = 0; j < dim; j++)
	     FT_IntrpStateVarAtCoords(fr,LIQUID_COMP,center,
				vel[j],getStateVel[j],&u[j],&vel[j][index]);
	    FT_IntrpStateVarAtCoords(fr,LIQUID_COMP,center,
				supersat,getStateSuper,&s,&s);

	    if (eqn_params->if_condensation == YES)
	        delta_R = R*R+2*eqn_params->K*s*dt;
	    else
	        delta_R = R*R;

	    if(delta_R < 0)
		R = 0;
	    else
	        R = sqrt(delta_R);

	    particle_array[i].radius = R;
	    /*save max and min radius*/
	    if(R > R_max)
		R_max = R;
	    if(R < R_min)
		R_min = R;
	    /*compute velocity for particle[i] with implicit method*/
	    for(j = 0; j < dim; ++j)
            {
		//update velocity and position of particles
                //according to the ode system
                //dx/dt = v
                //dv/dt = (u-v)/tau + g
                //assume tau and u are constant within one step
                if (eqn_params->if_sedimentation)
                {
                        center[j] += tau_p*(1-exp(-dt/tau_p))*cvel[j]
                                        +(dt-tau_p+tau_p*exp(-dt/tau_p))
                                        *(u[j]+gravity[j]*tau_p);
                        cvel[j] = exp(-dt/tau_p)*cvel[j]
                                + (1-exp(-dt/tau_p))*(u[j]+gravity[j]*tau_p);
                }
                else
                {
                        center[j] += tau_p*(1-exp(-dt/tau_p))*cvel[j]
                                        +(dt-tau_p+tau_p*exp(-dt/tau_p))
                                        *(u[j]);
                        cvel[j] = exp(-dt/tau_p)*cvel[j]
                                + (1-exp(-dt/tau_p))*(u[j]);
                }

		/*compute velocity
		if (eqn_params->if_sedimentation == YES) 
		    cvel[j] += (vel[j][index]/tau_p + gravity[j])*dt;
		else
		    cvel[j] += vel[j][index]/tau_p*dt;

		cvel[j] /= (1+dt/tau_p); 

	  	compute center of particle[i]
		center[j] += cvel[j]*dt;*/

		if(pp_numnodes() > 1)
		    continue;
		/*handle periodic drops for one processor*/
		T = rect_grid->U[j]-rect_grid->L[j];	
		if (center[j] > rect_grid->U[j])
		    center[j] = rect_grid->L[j]+fmod(center[j],T);
		if (center[j] < rect_grid->L[j])
		    center[j] = rect_grid->U[j]+fmod(center[j],T);

		if(isnan(center[j]))
		{
		    printf("center[%d]=nan, T = %f, domain=[%f,%f]\n",
				j,T,rect_grid->L[j],rect_grid->U[j]);
		    clean_up(ERROR);
		}
	    }

	    if (debugging("single_particle"))
	    {
	        printf("\nDrop[%d]:\n",i);
	        printf("Supersat: %f\n",s);
	        printf("Condensation rate: %20.19f\n",eqn_params->K);
	        printf("delta_R = %f\n",delta_R);
	        printf("dt = %f\n",dt);
	        printf("Radius:%15.14f\n",R);
	        printf("center:[%f,%f,%f]\n", center[0],center[1],center[2]);
	        printf("Flow_vel[%f,%f,%f]\nc_vel[%f,%f,%f]\n",
		    vel[0][index],vel[1][index],vel[2][index],
		    cvel[0],cvel[1],cvel[2]);
	        printf("Response time = %f\n",tau_p);
	    }
	}
	if(pp_numnodes() > 1)
	{
	    start_clock("ParallelExchangeParticle");
	    ParallelExchParticle(&particle_array,fr);
	    stop_clock("ParallelExchangeParticle");
	    eqn_params->particle_array = particle_array;
	}
	stop_clock("ParticlePropagate");
	printf("%d droplets in subdomain, max radius = %e, min radius = %e\n",
	eqn_params->num_drops,R_max,R_min);
}

extern void setParticleGroupIndex(
	PARTICLE* particle_array, 
	int size,
	int dim,
	int* Nclip,
	double* L,
	double* U)
{
	int i, j, group_id, group_crds[MAXD], group_gmax[MAXD];
	double *p_cen;
	double  block_size[MAXD];
	for (i = 0; i < dim; i++)
	{
	    block_size[i] = (U[i]-L[i])/Nclip[i];
	    group_gmax[i] = Nclip[i] - 1;
	}

	for (i = 0; i < size; i++)
	{
	    p_cen = particle_array[i].center;
	    for (j = 0; j < dim; j++)
		group_crds[j] = floor(p_cen[j]/block_size[j]);
	    group_id = d_index(group_crds,group_gmax,dim);
	    particle_array[i].Gindex = group_id;
	}
	if (debugging("Gindex"))
	{
	    printf("In setParticleGlobalIndex()\n");
	    printf("%d number of particles in the subdomain %d\n",
		    size, pp_mynode());
	}
}

extern void setParticleGlobalIndex(PARTICLE* particle_array, int size)
{
	int i,Gstart = 0;
	int num_nodes = pp_numnodes();
	int myid = pp_mynode();
	static int *n_drops;
	static boolean first = YES;
	if (debugging("trace"))
	    printf("Entering setParticleGlobalIndex()\n");
	if (first)
	{
	    first = NO;
            FT_VectorMemoryAlloc((POINTER*)&n_drops,num_nodes,sizeof(int));
	}
	for (i = 0; i < num_nodes; i++) n_drops[i] = 0;
	n_drops[myid] = size;
	pp_gsync();
	pp_global_imax(n_drops,num_nodes);
	
	Gstart = 0;
	for (i = 0; i < myid; i++)
	    Gstart += n_drops[i];
	for (i = 0; i < size; i++)
	    particle_array[i].Gindex = (Gstart + i);
	if (debugging("Gindex"))
	{
	    printf("In setParticleGlobalIndex()\n");
	    printf("%d number of particles in the subdomain %d\n",
		    size, myid);
	    printf("Particle index starts and ends= [%d, %d]\n",
		    particle_array[0].Gindex,particle_array[size-1].Gindex);
	    printf("Number of particles in each subdomain:\n");
	    for (i = 0; i < num_nodes;i++)
		printf("%d ",n_drops[i]);
	    printf("\n%d particles before this subdomain\n",Gstart);
	}
}

extern void printDropletsStates(Front* front, char* outname)
{
        char filename[100];
        FILE *outfile;
	PARAMS* eqn_params = (PARAMS*)front->extra2;
	int dim = Dimension(front->interf);
        PARTICLE* particle_array = eqn_params->particle_array;
	int i,j,num_drops = eqn_params->num_drops;
            
        sprintf(filename,"%s/state.ts%s",outname,
                        right_flush(front->step,7));
        sprintf(filename,"%s-drops",filename);
        outfile = fopen(filename,"w");
	for (i = 0; i < num_drops; i++)
	{
	    fprintf(outfile,"%24.18g\n",particle_array[i].radius);
	    for (j = 0; j < dim; j++)
		fprintf(outfile,"%24.18g ",particle_array[i].center[j]);
	    fprintf(outfile,"\n");
	    for (j = 0; j < dim; j++)
		fprintf(outfile,"%24.18g ",particle_array[i].vel[j]);
	    fprintf(outfile,"\n");
	}
	fclose(outfile);
}

/****************PLOT FUNCTIONS****************************************/
extern void gv_plot_scatter(Front* front)
{
	PARAMS* eqn_params = (PARAMS*)front->extra2;
	PARTICLE* particle_array = eqn_params->particle_array;
	int dim = front->rect_grid->dim;
	int i, j, n = eqn_params->num_drops;
	double *coords;
	char *outname = OutName(front);
	char fname[256];
	FILE* file;
	
	sprintf(fname,"%s/particle.off",outname);	
	file = fopen(fname,"w");
	fprintf(file,"appearance {linewidth 3}\n");
	fprintf(file,"VECT\n");
	fprintf(file,"%d %d %d\n\n",n,n,n);
	/*num of polyline, num of point, num of colors*/

        for (i = 0; i < n; i++)
        {
            fprintf(file,"1 ");/*num of points for each polyline*/
        }
	fprintf(file,"\n\n");

        for (i = 0; i < n; i++)
        {
            fprintf(file,"1 ");/*num of colors for each polyline*/
        }
        fprintf(file,"\n\n");
	
	for (i = 0; i < n; i++)
	{
	    coords = particle_array[i].center;
	    if (dim == 2)
	    {
	            fprintf(file,"%f %f 0\n",coords[0],coords[1]);/*coord for each points*/
	    }
	    else if (dim == 3)
	    {
		fprintf(file,"%f %f %f\n",coords[0],coords[1],coords[2]);
	    }
	}
	fprintf(file,"\n");
	for (i = 0; i < n; i++)
	    fprintf(file,"0 0 1 1\n");/*color for each point*/
	fclose(file);
}

extern void vtk_plot_scatter(Front* front)
{
        PARAMS* eqn_params = (PARAMS*)front->extra2;
        PARTICLE* particle_array = eqn_params->particle_array;
        int dim = front->rect_grid->dim;
        int i, j, count, n;
        double *coords;
        char fname[256],dname[256];
        FILE* file;

	count =0;
	for (i = 0; i < eqn_params->num_drops; i++)
	{
	    if (particle_array[i].radius != 0)
		count ++;
	}
	n = count; /*number of droplets with positive radius*/

        sprintf(dname,"%s/vtk",OutName(front));
        if (pp_mynode() == 0)
        {
            if (!create_directory(dname,YES))
            {
                screen("Cannot create directory %s\n",dirname);
                clean_up(ERROR);
            }
        }
        pp_gsync();
        sprintf(dname,"%s/vtk.ts%s",dname,right_flush(front->step,7));
        if (pp_numnodes() > 1)
            sprintf(dname,"%s-nd%s",dname,right_flush(pp_mynode(),4));
        if (!create_directory(dname,YES))
        {
            screen("Cannot create directory %s\n",dname);
            clean_up(ERROR);
        }

	sprintf(fname,"%s/particle.vtk",dname);
        file = fopen(fname,"w");
        fprintf(file,"# vtk DataFile Version 3.0\n");
        fprintf(file,"%s\n","particles");
        fprintf(file,"ASCII\n");

	fprintf(file,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(file,"POINTS %d FLOAT\n",n);

        for (i = 0; i < eqn_params->num_drops; i++)
        {
	    if (particle_array[i].radius == 0)
		continue;
            coords = particle_array[i].center;
            if (dim == 2)
                fprintf(file,"%f %f 0\n",coords[0],coords[1]);
            else if (dim == 3)
                fprintf(file,"%f %f %f\n",coords[0],coords[1],coords[2]);
        }
	fprintf(file,"POINT_DATA %d\n",n);
	fprintf(file,"SCALARS radius FLOAT\n");
	fprintf(file,"LOOKUP_TABLE default\n");
	for (i = 0; i < eqn_params->num_drops; i++)
        {
            if (particle_array[i].radius == 0)
                continue;
            fprintf(file,"%20.14f\n",particle_array[i].radius);
        }
	fclose(file);
}

extern void vtk_plot_sample_traj(Front* front)
{
        PARAMS* eqn_params = (PARAMS*)front->extra2;
	RECT_GRID *rect_grid = front->rect_grid;
        int dim = rect_grid->dim;
        int i, j;
        double *coords,max_dist;
        char *outname = OutName(front);
        char fname[256];
        FILE* file;
	boolean ignore_this;
	/*static array for preserving trajectory*/
	static double traj[MAXD][MAX_STEP];
	static int step = 0;
	int max_step = step;

	if (pp_numnodes() > 1)
	    return;
	for(i = 0; i < dim; i++)
	{
	    traj[i][step] = eqn_params->particle_array[0].center[i];
	}
	if(step == 0)
	{  
	    step++;
	    return;
	}
        sprintf(fname,"%s/vtk/vtk.ts%s/",outname,
				right_flush(front->step,7));
	if (!create_directory(fname,NO))
        {
            printf("Cannot create directory %s\n",fname);
            clean_up(ERROR);
        }
	sprintf(fname,"%s/trajectory.vtk",fname);
        file = fopen(fname,"w");
        fprintf(file,"# vtk DataFile Version 3.0\n");
        fprintf(file,"%s\n","sample_particles");
        fprintf(file,"ASCII\n");

	fprintf(file,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(file,"POINTS %d FLOAT\n",step+1);

	for (i = 0; i < step+1; i++)
	{
            if (dim == 2)
                fprintf(file,"%f %f 0\n",traj[0][i],traj[1][i]);
            else if (dim == 3)
                fprintf(file,"%f %f %f\n",traj[0][i],traj[1][i],traj[2][i]);
	}
        for (i = 0; i < step; i++)
        {
            ignore_this = NO;
            for (j = 0; j < dim; j++)
            {
                max_dist =  sqr(traj[j][i]-traj[j][i+1]);
                max_dist -= 0.5*sqr(rect_grid->U[j]-rect_grid->L[j]);
                if (max_dist > 0)
                    ignore_this = YES;
            }
            if (ignore_this == YES )
                max_step--;
        }
	fprintf(file,"CELLS %d %d\n",max_step,3*max_step);
	for (i = 0; i < step; i++)
	{
	    ignore_this = NO;
	    for (j = 0; j < dim; j++)
	    {
		max_dist =  sqr(traj[j][i]-traj[j][i+1]);
		max_dist -= 0.5*sqr(rect_grid->U[j]-rect_grid->L[j]);
		if (max_dist > 0)
		    ignore_this = YES;	
	    }
	    if (ignore_this == YES )
		continue;
	    fprintf(file,"2 %d %d\n",i,i+1);
	}
	fprintf(file,"CELL_TYPES %d\n",max_step);
	for (i = 0; i < max_step; i++)
        {
            fprintf(file,"3\n");
        }
	step ++;
	fclose(file);
}
/*********************End Plot Function***********************************/

/********************Statistics Functions*******************************/
extern void Deviation(double* array,Front* front,double &Mean,double &Var)
{
	int i,j,k,size=0,index;
	int imin,imax,jmin,jmax,kmin,kmax;
	INTERFACE* grid_intfc = front->grid_intfc;
        RECT_GRID* top_grid = &topological_grid(grid_intfc);
	int dim = top_grid->dim;
	int* lbuf = front->rect_grid->lbuf;
        int* ubuf = front->rect_grid->ubuf;
	int* top_gmax = top_grid->gmax;
	
	/*compute sum using Kahan summation alogrithm*/
	double c = 0.0, t, y;

	Mean = 0.0;
	Var = 0.0;
	switch(dim)
	{
		case 2:
		    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
            	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
            	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
            	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
		    size = (imax-imin+1)*(jmax-jmin+1); 
		    pp_gsync();
		    pp_global_isum(&size,1);
		    for (j = jmin; j <= jmax; ++j)
		    for (i = imin; i <= imax; ++i)
		    {
			index = d_index2d(i,j,top_gmax);
			/*Kahan alogrithm*/
			y = array[index]-c;
			t = Mean+y;
			c = (t-Mean)-y;
			Mean = t;
		    }
        	    pp_gsync();
        	    pp_global_sum(&Mean,1);
		    Mean /= size;
		    c = 0.0;
		    for (j = jmin; j <= jmax; ++j)
                    for (i = imin; i <= imax; ++i)
                    {
                        index = d_index2d(i,j,top_gmax);
			/*Kahan alogrithm*/
			y = sqr(array[index]-Mean)-c;
                        t = Var + y;
			c = (t-Var)-y;
			Var = t; 
                    } 
 		    pp_gsync();
		    pp_global_sum(&Var,1);
		    Var /= size; 
		    break;
		case 3:
		    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	            jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
            	    kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
            	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
            	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
            	    kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
		    size = (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1);
		    pp_global_isum(&size,1);
		    for (k = kmin; k <= kmax; ++k)
		    for (j = jmin; j <= jmax; ++j)
                    for (i = imin; i <= imax; ++i)
		    {
	    		index = d_index3d(i,j,k,top_gmax);
			/*Kahan alogrithm*/
			y = array[index]-c;
			t = Mean+y;
			c = (t-Mean)-y;
			Mean = t;
		    }
		    pp_gsync();
                    pp_global_sum(&Mean,1);
                    Mean /= size;
		    c = 0.0;
		    for (k = kmin; k <= kmax; ++k)
		    for (j = jmin; j <= jmax; ++j)
                    for (i = imin; i <= imax; ++i)
		    {
			index = d_index3d(i,j,k,top_gmax);
			/*Kahan alogrithm*/
			y = sqr(array[index]-Mean)-c;
                        t = Var + y;
			c = (t-Var)-y;
			Var = t; 
		    }
		    pp_gsync();
                    pp_global_sum(&Var,1);
                    Var /= size;
		    break;
		default:
		    printf("Unknown dim = %d\n",dim);
	}
}

extern void Deviation(double* array,int size,double &Mean,double &Var)
{
        int i,nzeros;
        nzeros = 0;
        Mean = 0;
        for (i = 0; i < size; i++)
        {
            if(array[i] != 0)
	    {
                nzeros++;
                Mean += array[i];
	    }
        }
#if defined __MPI__
	pp_gsync();
        pp_global_isum(&nzeros,1);
        pp_global_sum(&Mean,1);
#endif
        Mean /= nzeros;

        Var = 0;
        for (i = 0; i < size; i++)
        {
            if(array[i] != 0)
            {
                Var += sqr(array[i]-Mean);
            }
        }
#if defined __MPI__
	pp_gsync();
        pp_global_sum(&Var,1);
#endif
        Var /= nzeros;
        return;
}

extern double* ComputePDF(
        double *array,
        int size,
        double &bin_size,
        int num_bins,
        double &var_min,
        double &var_max)
{
	return ComputePDF(array,size,bin_size,num_bins,var_min,var_max,NO);
}

extern double* ComputePDF(
	double *array, 
	int size, 
	double &bin_size,
	int num_bins,
	double &var_min,
	double &var_max,
	boolean ignore_zero)
{
	int i, j, total_num;
	int myid = pp_mynode();
	static double *PDF = NULL;
	static int max_bin_num = 0;

	var_min =  HUGE;
	var_max = -HUGE;
	
	if (debugging("trace"))
	    printf("Entering computePDF\n");
	for (i = 0; i < size; i++)
	{
	    if (array[i] < var_min)	    
		var_min = array[i];
	    if (array[i] > var_max)
		var_max = array[i];
	}

#if defined __MPI__
	pp_gsync();
	pp_global_max(&var_max,1);	
	pp_global_min(&var_min,1);	

#endif
	bin_size = (var_max - var_min)/(num_bins-1);
	//num_bins = ceil((var_max - var_min)/bin_size);
	if (num_bins > max_bin_num)
	{
	    max_bin_num = num_bins;
	    if (PDF != NULL)
	       free_these(1,PDF);
	    FT_VectorMemoryAlloc((POINTER*)&PDF,num_bins,FLOAT);
	}
	else if(num_bins == 0)
	{
	    if (PDF == NULL)
	        FT_VectorMemoryAlloc((POINTER*)&PDF,1,FLOAT);
	    PDF[0] = 1.0;
	    num_bins = 1;
	    return PDF;
	}
	for (j = 0; j < num_bins; j++)
		PDF[j] = 0.0;

	for (i = 0; i < size; i++)
	{
	    if (array[i] == 0 && ignore_zero == YES)
		continue;
	    for (j = 0; j < num_bins; j++)
	    {
	        if (array[i] >= var_min+(j-0.5)*bin_size &&
		    array[i] < var_min+(j+0.5)*bin_size)
	        {
		    PDF[j] += 1.0;
		    break;
	        }
	    }
	}
#if defined __MPI__	
	pp_gsync();
	pp_global_sum(PDF,num_bins);
 #endif
	total_num = 0;
	/*normalize PDF*/
	for (j = 0; j < num_bins; j++)
	    total_num += PDF[j];

	for (j = 0; j < num_bins; j++)
	    PDF[j] = double(PDF[j])/(bin_size*double(total_num));
	if (debugging("trace"))
	    printf("Leaving computePDF\n");
	return PDF;
}
