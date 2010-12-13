#include "bem_mpi.h"
//#define _debug
//#define _debug2

#define nodes(i,j)    nodes[(i) + nnod*(j)]
#define elements(i,j) elements[(i) + nelem*(j)]
#define Br(i,j) Br[(i) + tasksize*(j)]
#define Bi(i,j) Bi[(i) + tasksize*(j)]

#define nod_list(i,j) nod_list[(i)+ num_nodes*(j)]
#define elem_list(i,j) elem_list[(i)+num_elem*(j)]

void main_build_matrix(
    double *nodes, double *elements,
    ulong nnod, ulong nelem,
    double omegar, double omegai, double D, ulong mynum_procs,
    double*& Ar, double*& Ai, double*& Br, double*& Bi,
	mpi::environment& env, mpi::communicator& world,
	std::vector<ulong>& tasksizes) {
	
	
	// local variables
	int ngp = 4;
	ulong  i, j, el;

	// Element area, normal and Guassian points for integration.
	double *delG = new double[nelem];
	double *normal_x = new double[nelem];
	double *normal_y = new double[nelem];
	double *normal_z = new double[nelem];
	double *eta1 = new double[ngp];
	double *eta2 = new double[ngp];
	double *w = new double[ngp];
	for(i=0; i<ngp; ++i)
		*(eta1+i) = *(eta2+i) = *(w+i) = 0.0;

#ifdef _debug
	std::cout << "get_gausspts_3d\n";
#endif
	/* getting gauss points for linear triangular element */
	get_gausspts_3d(eta1, eta2, w, ngp);



#ifdef _debug
	std::cout << "get_area_normals\n";
#endif

#if defined(_OPENMP) & defined(_debug)
	double starttime = omp_get_wtime();
#endif	
	get_area_normals(elements, nodes, nelem, nnod, delG, 
					 normal_x, normal_y, normal_z);
#if defined(_OPENMP) & defined(_debug)	
	std::cout << " process " << world.rank() << ": Area calculation took " <<
				 omp_get_wtime()-starttime) " seconds.\n";
#endif

#ifdef _debug
	std::cout << "mxGetPr(A)\n";
#endif
#ifdef _OPENMP
	if (mynum_procs < omp_get_num_procs())
		omp_set_num_threads(mynum_procs);
	else
		omp_set_num_threads(omp_get_num_procs());
#ifdef _debug
	std::cout << " CPUs available: " << omp_get_num_procs() << " , " <<
				 mynum_procs << std::endl;
#endif
#endif

	ulong tasksize = tasksizes[world.rank()];
		
	double xi, yi, zi;
	double xel[3], yel[3], zel[3];
	ulong jel[3];
	double* int_val_real;
	double* int_val_img;
	int_val_real = new double[3];
	int_val_img  = new double[3];
	
	for (ulong i=world.rank(); i<nnod; i+=world.size()) {
		ulong node_id = i+1;
		ulong nidx_local = i / world.size();
		#ifdef _debug2
		std::cout << i << std::endl;
		#endif
		for(el = 0; el < nelem; ++el) {
			jel[0] = (ulong) elements(el,0);	
			jel[1] = (ulong) elements(el,1);
			jel[2] = (ulong) elements(el,2);
			xel[0] = nodes(jel[0]-1,0);
			xel[1] = nodes(jel[1]-1,0);
			xel[2] = nodes(jel[2]-1,0);
			yel[0] = nodes(jel[0]-1,1);
			yel[1] = nodes(jel[1]-1,1);
			yel[2] = nodes(jel[2]-1,1);
			zel[0] = nodes(jel[0]-1,2);
			zel[1] = nodes(jel[1]-1,2);
			zel[2] = nodes(jel[2]-1,2);
			#ifdef _debug2
			std::cout << "el = " << el << std::endl;
			#endif
			#ifdef _debug2
			std::cout << jel[0] << ' ' << jel[1] << ' ' << jel[2] << std::endl;
			#endif
			if(node_id == jel[0]) {
				for(j=0; j<3; *(int_val_real+j)=0.0, *(int_val_img+j)=0.0, ++j);
				singular_integrand_polar(int_val_real, int_val_img, omegar,
				    					omegai, D, delG[el],eta1, eta2, w, ngp);
				Br(nidx_local, jel[0]-1) += int_val_real[0];
				Bi(nidx_local, jel[0]-1) += int_val_img[0];
				Br(nidx_local, jel[1]-1) += int_val_real[1];
				Bi(nidx_local, jel[1]-1) += int_val_img[1];
				Br(nidx_local, jel[2]-1) += int_val_real[2];
				Bi(nidx_local, jel[2]-1) += int_val_img[2];
			}	
		    else if(node_id == jel[1]) {
				for(j=0; j<3; *(int_val_real+j)=0.0, *(int_val_img+j)=0.0, ++j);
				singular_integrand_polar(int_val_real, int_val_img, omegar,
										omegai, D, delG[el],eta1, eta2, w, ngp);
				Br(nidx_local, jel[0]-1) += int_val_real[2];
				Bi(nidx_local, jel[0]-1) += int_val_img[2];
				Br(nidx_local, jel[1]-1) += int_val_real[0];
				Bi(nidx_local, jel[1]-1) += int_val_img[0];
				Br(nidx_local, jel[2]-1) += int_val_real[1];
				Bi(nidx_local, jel[2]-1) += int_val_img[1];
			}	
			else if(node_id == jel[2]) {
				for(j=0; j<3; *(int_val_real+j)=0.0, *(int_val_img+j)=0.0, ++j);
				singular_integrand_polar(int_val_real, int_val_img, omegar,
										omegai, D, delG[el],eta1, eta2, w, ngp);
				Br(nidx_local, jel[0]-1) += int_val_real[1];
				Bi(nidx_local, jel[0]-1)  += int_val_img[1];
				Br(nidx_local, jel[1]-1) += int_val_real[2];
				Bi(nidx_local, jel[1]-1)  += int_val_img[2];
				Br(nidx_local, jel[2]-1) += int_val_real[0];
				Bi(nidx_local, jel[2]-1)  += int_val_img[0];
			}
			else {
				xi = nodes(i,0);
				yi = nodes(i,1);
				zi = nodes(i,2);
				compute_integral_gausspts(Ar, Ai, Br, Bi, ngp, eta1, eta2, w,
										jel, xel, yel, zel,
										tasksize, nidx_local,
										xi, yi, zi,	delG[el],
										normal_x[el], normal_y[el], normal_z[el]
										, omegar, omegai, D);
			}
		}
	}
	delete [] int_val_real;
	delete [] int_val_img;
	delete [] delG;
	delete [] normal_x;
	delete [] normal_y;
	delete [] normal_z;
	delete [] eta1;
	delete [] eta2;
	delete [] w;
	std::cout << "  *Process " << world.rank() << 
		" finished its sub-matrix construction: " << tasksize << ' ' << std::endl;
}

static void compute_integral_gausspts(
	double *AAr, double *AAi, double *BBr, double *BBi,
	int ngpts, 	double *eta1, double *eta2, double *w,
	ulong *jelem, double *xelem, double *yelem, double *zelem,
	ulong tasksize, ulong nidx_local,
	double xx, double yy, double zz, double area,
	double norm_x, 	double norm_y, double norm_z,
	double omegar, double omegai, double diff_coeff) {	

	double phi[4], xs, ys, zs, ri, drdn;
	ulong t_result;
	double z_real, z_img, gi_real, gi_img;
	double dgdr_real, dgdr_img, dgdn_real, dgdn_img;
	double param_real, param_img, *aa_real, *aa_img, *bb_real, *bb_img;
	double temp_result, xss, yss, zss, four_PI_diff;

	param_real = omegar;
	param_img = omegai;
	aa_real = AAr + ((nidx_local)-tasksize);
	aa_img  = AAi + ((nidx_local)-tasksize);
	bb_real = BBr + ((nidx_local)-tasksize);
	bb_img  = BBi + ((nidx_local)-tasksize);
	four_PI_diff = 4*PI*diff_coeff;

	for(int k = ngpts -1; k>=0; --k) {
		phi[0] = eta1[k];
		phi[1] = eta2[k];
		phi[2] = 1 - eta1[k] - eta2[k];
		
		xs = xelem[0]*phi[0] + xelem[1]*phi[1] + xelem[2]*phi[2];
		ys = yelem[0]*phi[0] + yelem[1]*phi[1] + yelem[2]*phi[2];
		zs = zelem[0]*phi[0] + zelem[1]*phi[1] + zelem[2]*phi[2];

		xss = xs-xx;
		yss = ys-yy;
		zss = zs-zz;

		ri = sqrt(xss*xss + yss*yss + zss*zss);
		drdn = (norm_x*xss + norm_y*yss + norm_z*zss)/(area*ri);

		/* z = param*ri */
		z_real = (param_real)*ri;
		z_img = (param_img)*ri;
		/* gi = exp(-z)/(4*PI*ri*diff_coeff) */
		cis((-1.0*(z_real)), (-1.0*(z_img)), gi_real, gi_img);
		(gi_real) /= four_PI_diff*ri;
		(gi_img) /= four_PI_diff*ri;
		/* dgdr = gi*((-1.0/ri) - param) */
		mult_complex(gi_real, gi_img, ((-1.0/ri)-(param_real)), 
				(-1.0*(param_img)), dgdr_real, dgdr_img);
		/* dgdn = dgdr*drdn */
		dgdn_real = (dgdr_real)*drdn;
		dgdn_img = (dgdr_img)*drdn;

		/*make use of dgdn and gi so we only have to multiply by phi[j] in
		 the inner loop...effectively hoisting the multiplications into this loop */
		temp_result = area*w[k]*diff_coeff;
		dgdn_real *= temp_result;
		dgdn_img *= temp_result;
		gi_real *= area*w[k];
		gi_img *= area*w[k];

		for(int j=2; j>=0; --j) {
			t_result = jelem[j]*tasksize;
			
			/* AA[nidx_local][jelem[j]-1] += phi[j]*dgdn*area*w[k]*diff_coeff */
			(*(aa_real + t_result)) +=  phi[j]*dgdn_real;
			(*(aa_img + t_result)) += phi[j]*dgdn_img;

			/* BB[nidx_locale][jelem[j]-1] += phi[j]*gi*area*w[k] */
			(*(bb_real + t_result)) += phi[j]*(gi_real);
			(*(bb_img + t_result)) += phi[j]*(gi_img);
			
		}
		 
	}
	
}

void singular_integrand_polar(double* int_val_real, double* int_val_img,
 				double omegar, double omegai, double diff_coeff,
				double area, double *gp1, double *gp2, double *wgt, int ngpts) {
   
    double basis[3], tau1, tau2;
    int k, j;
    double real_part, img_part;

    double sqrttwo = sqrt(2.0);

    for(k=0; k<ngpts; ++k) {
        tau1 = sqrttwo*gp1[k]*cos(PI*gp2[k]/4.);
        tau2 = sqrttwo*gp1[k]*sin(PI*gp2[k]/4.);
        basis[0] = tau1;
        basis[1] = tau2;
        basis[2] = 1 - tau1 - tau2;

	 	/* do int_value[j] = int_value[j] + 
				exp(-param*gp1[k]*sqrt(2.0))*basis[j]*area*wgt[k]; */
		/* in two states, do the above statements. First do the exponentiation */
		real_part = omegar;
		img_part  = omegai;
		real_part = (real_part)*-1*gp1[k]*sqrttwo;
		img_part  = (img_part )*-1*gp1[k]*sqrttwo;

			/*exp(-param*pg1[k]*sqrt(2.0))*/
		cis(real_part, img_part, real_part, img_part);	
		 
        for(j=0; j<3; ++j) {
			/* then do multiplication and addition */
			(*(int_val_real+j)) += (real_part)*basis[j]*area*wgt[k];
			(*(int_val_img+j))  += ( img_part)*basis[j]*area*wgt[k];
        }
    }
		
    for(j=0; j<3; ++j) {
		/*do int_value[j] = int_value[j]*sqrt(2.0)/(16*diff_coeff) */
		(*(int_val_real+j)) *= sqrttwo/(16.*diff_coeff);
		(*(int_val_img+j))  *= sqrttwo/(16.*diff_coeff);
    }
	
}

void get_gausspts_3d(double *eta1, double *eta2, double *w, int ngp) {
	
	/* getting gass points for 3D triangular element */
	if(ngp == 3) {
		eta1[0] = 1.0/2.0;
		eta1[1] = 0.0;
		eta1[2] = 1.0/2.0;
		eta2[0] = 1.0/2.0;
		eta2[1] = 1.0/2.0;
		eta2[2] = 0.0;
		
		w[0]    = 1.0/3.0;
		w[1]    = 1.0/3.0;
		w[2]    = 1.0/3.0;

	}
	else 
		if(ngp == 4) {
			eta1[0] = 1.0/3.0;
			eta1[1] = 3.0/5.0;
			eta1[2] = 1.0/5.0;
			eta1[3] = 1.0/5.0;
			eta2[0] = 1.0/3.0;
			eta2[1] = 1.0/5.0;
			eta2[2] = 3.0/5.0;
			eta2[3] = 1.0/5.0;
		
			w[0]    = -27.0/48.0;
			w[1]    =  25.0/48.0;
			w[2]    =  25.0/48.0;
			w[3]    =  25.0/48.0;
		}
}

void get_area_normals(double *elem_list, double *nod_list, ulong num_elem,
			ulong num_nodes, double *area_elem_list, double *normal_x_list,
			double *normal_y_list, double *normal_z_list) {

	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	double dx_deta1, dx_deta2, dy_deta1, dy_deta2, dz_deta1, dz_deta2;
	ulong jelem[3], ee;

	for(ee = 0; ee < num_elem; ++ee) {
		jelem[0] = (ulong) elem_list(ee,0);
		jelem[1] = (ulong) elem_list(ee,1);
		jelem[2] = (ulong) elem_list(ee,2);

		x1 = nod_list(jelem[0]-1,0);
		x2 = nod_list(jelem[1]-1,0);
		x3 = nod_list(jelem[2]-1,0);

		y1 = nod_list(jelem[0]-1,1);
		y2 = nod_list(jelem[1]-1,1);
		y3 = nod_list(jelem[2]-1,1);

		z1 = nod_list(jelem[0]-1,2);
		z2 = nod_list(jelem[1]-1,2);
		z3 = nod_list(jelem[2]-1,2);

		dx_deta1 = x1 - x3;
		dx_deta2 = x2 - x3;
		dy_deta1 = y1 - y3;
		dy_deta2 = y2 - y3;
		dz_deta1 = z1 - z3;
		dz_deta2 = z2 - z3; 

		normal_x_list[ee] = dy_deta1*dz_deta2 - dy_deta2*dz_deta1;
		normal_y_list[ee] = dz_deta1*dx_deta2 - dx_deta1*dz_deta2;
		normal_z_list[ee] = dx_deta1*dy_deta2 - dy_deta1*dx_deta2; 
		area_elem_list[ee] = sqrt(pow(normal_x_list[ee],2) +
								  pow(normal_y_list[ee],2) +
								  pow(normal_z_list[ee],2));
	}
}

// Reads the input mesh from 'matfn' .mat Matlab file
// 'matfn' should have following mxArrays names:
// 'nodes', 'elements', 'omega', 'D', 'num_procs'
void get_mesh_from_matlab(const char *matfn, mxArray *&mxNodes,
	mxArray *&mxElements,
	mxArray*& mxOmega, mxArray*& mxD, mxArray*& mxNumProcs) {

	MATFile *pmat;
	
	pmat = matOpen(matfn, "r");
	if (pmat == NULL) {
		std::cerr << " Error opening file " << matfn << '\n';
	    exit(1);
	}
	// read 'nodes'
	mxNodes = matGetVariable(pmat,"nodes");
	if (mxNodes == NULL) {
		std::cerr << " Error reading 'nodes' matrix from .mat file\n";
		exit(1);
	}
	// read 'elements'
	mxElements = matGetVariable(pmat,"elements");
	if (mxElements == NULL) {
		std::cerr << " Error reading 'elements' matrix from .mat file\n";
		exit(1);
	}
	// read 'omega'
	mxOmega = matGetVariable(pmat,"omega");
	if (mxOmega == NULL) {
		std::cerr << " Error reading 'omega' matrix from .mat file\n";
		exit(1);
	}
	// read 'D'
	mxD = matGetVariable(pmat,"D");
	if (mxD == NULL) {
		std::cerr << " Error reading 'D' matrix from .mat file\n";
		exit(1);
	}
	// read 'mynum_procs'
	mxNumProcs = matGetVariable(pmat,"num_procs");
	if (mxNumProcs == NULL) {
		std::cerr << " Error reading 'num_procs' matrix from .mat file\n";
		exit(1);
	}
	if (matClose(pmat) != 0) {
		std::cout << " Could not close the .mat file: " << __FILE__ << ' '
			<< __LINE__ << std::endl;
		exit(1);
	}
}

void write_AB_to_matlab(const char* varname,
		std::vector< std::vector<double> >& All_AB,
		std::vector< ulong >& tasksizes,
		ulong& nnod) {
	
	MATFile *pmat;
	char matfilename[256];
	sprintf(matfilename,"%s.mat",varname);
	ulong worldsize = tasksizes.size();
	
	pmat = matOpen(matfilename, "w");
	if (pmat == NULL) {
		std::cerr << " Error opening file " << matfilename << '\n';
	    exit(1);
	}

	mxArray *mxAB = mxCreateDoubleMatrix(nnod, nnod, mxREAL);
	double *matpointer = mxGetPr(mxAB);
	for (ulong i=0; i<nnod; ++i) {
		ulong procn = i % worldsize;
		ulong rownum = i / worldsize;
		ulong tasksize_local = tasksizes[procn];
		for (ulong j=0; j<nnod; ++j) {
			matpointer[i+nnod*j] = All_AB[procn][rownum+tasksize_local*j];
		}
	}
	int st = matPutVariable(pmat,varname,mxAB);
	if (st!=0) {
		std::cerr << " Could not right variable to .mat file: "
			<< varname << '\n';
	}
	
	if (matClose(pmat) != 0) {
		std::cout << " Could not close the .mat file: " << __FILE__ << ' '
			<< __LINE__ << std::endl;
		exit(1);
	}
	mxDestroyArray(mxAB);
}
int main(int argc, char *argv[]) {
	
	// Boost.MPI
	mpi::environment env(argc,argv);
	mpi::communicator world;
	
	double *nodes, *elements;
	ulong nnod, nelem;
	double omegar, omegai;
	double D;
	ulong mynum_procs;
	double *Arp, *Aip, *Brp, *Bip;
	
	// Read the input mesh from a .mat file
	mxArray *mxnodes, *mxelements, *mxomega, *mxD, *mxmynum_procs;
	get_mesh_from_matlab("bem_mpi_input.mat", mxnodes, mxelements,
		mxomega, mxD, mxmynum_procs);
	nodes = mxGetPr(mxnodes);
	nnod = mxGetM(mxnodes);
	elements = mxGetPr(mxelements);
	nelem = mxGetM(mxelements);
	omegar = *mxGetPr(mxomega);
	if (mxIsComplex(mxomega))
		omegai = *(mxGetPi(mxomega));
	else
		omegai = 0;
	D = *mxGetPr(mxD);
	
	mynum_procs = (ulong) *mxGetPr(mxmynum_procs);
	
	std::vector<ulong> tasksizes(world.size(), nnod / world.size());
	ulong R = nnod % world.size();
	for (ulong i=0; i<tasksizes.size(); ++i) {
		if (i < R) {
			++tasksizes[i];
		}
	}
	ulong tasksize_local = tasksizes[world.rank()];
	
	std::vector<double> VecAr(tasksize_local*nnod);
	std::vector<double> VecAi(tasksize_local*nnod);
	std::vector<double> VecBr(tasksize_local*nnod);
	std::vector<double> VecBi(tasksize_local*nnod);
	Arp = &VecAr[0];
	Aip = &VecAi[0];
	Brp	= &VecBr[0];
	Bip = &VecBi[0];	
	// Build components of A and B matrices on different nodes of a cluster
	main_build_matrix(nodes, elements, nnod, nelem, omegar, omegai,
		D, mynum_procs, Arp, Aip, Brp, Bip, env, world, tasksizes);

	// Gather pieces of A and B from different nodes and assemble a full matrix
	if (world.rank()==0) {
		std::vector< std::vector<double> > All_Ar(world.size());
		std::vector< std::vector<double> > All_Ai(world.size());
		std::vector< std::vector<double> > All_Br(world.size());
		std::vector< std::vector<double> > All_Bi(world.size());
		
		std::cout << " Gathering matrices from different workers...";
		gather(world, VecAr, All_Ar, 0);
		gather(world, VecAi, All_Ai, 0);
		gather(world, VecBr, All_Br, 0);
		gather(world, VecBi, All_Bi, 0);
		std::cout << " done." << std::endl;
		std::cout << " Writing Matlab files...";
		write_AB_to_matlab("Ar",All_Ar,tasksizes,nnod);
		write_AB_to_matlab("Ai",All_Ai,tasksizes,nnod);
		write_AB_to_matlab("Br",All_Br,tasksizes,nnod);
		write_AB_to_matlab("Bi",All_Bi,tasksizes,nnod);
		std::cout << " done." << std::endl;
		// Create a file to tell Matlab it can read the .mat files
		std::ofstream ofs("bem_mpi.lock");
		ofs.close();
	}
	else {
		gather(world, VecAr, 0);
		gather(world, VecAi, 0);
		gather(world, VecBr, 0);
		gather(world, VecBi, 0);
	}
	
	mxDestroyArray(mxnodes);
	mxDestroyArray(mxelements);
	mxDestroyArray(mxomega);
	mxDestroyArray(mxD);
	mxDestroyArray(mxmynum_procs);
	return 0;
}






