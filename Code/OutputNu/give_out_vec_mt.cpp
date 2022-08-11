/**                                                                                                                                
 * GENERAL REMARKS                                                                                                                 
 *                                                                                                                                 
 *  This code is freely available under the following conditions:                                                                  
 *                                                                                                                                 
 *  1) The code is to be used only for non-commercial purposes.                                                                    
 *  2) No changes and modifications to the code without prior permission of the developer.                                         
 *  3) No forwarding the code to a third party without prior permission of the developer.                                          
 *                                                                                                                                 
 *  			MTCalc_with_DFP_COCR                                                                                       
 *  Building a spline on the ground surface to output the value of the required function at an arbitrary point (receiver position) 
 *                                                                                                                                 
 *  Written by Ph.D. Petr A. Domnikov                                                                                              
 *  Novosibirsk State Technical University,                                                                                        
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                                              
 *  p_domnikov@mail.ru                                                                                                             
 *  Version 1.2 March 10, 2021                                                                                                     
*/                                                                                                                                 
#include "stdafx.h"
#include "OutputArbitrary.h"
#include "give_out_vec_mt.h"
#include "GeoPrepDocSettings.h"
extern ofstream logfile;
extern GeoPrepDocSettings GPDocSettings;
//-----------------------------------------------------------------------------
Give_out_vec_mt::Give_out_vec_mt(Vec_Prep_Data *d, T_Mapping_Vec *tmap, double *v3dat)
{
	this->d = d;
	this->tmap = tmap;
	this->v3dat = v3dat;

	this->H = new std::complex<double>[d->n_pointresB][3];
	this->E = new std::complex<double>[d->n_pointresE][3];
	this->impedance = new double[d->n_pointresB];
	this->rho = new double[d->n_pointresB];

	this->H_1d = new std::complex<double>[d->n_pointres][2];
	this->E_1d = new std::complex<double>[d->n_pointres][2];

	if (d->tasktype==0 && !d->fdirect)
		Compute_1d_field();
}
//-----------------------------------------------------------------------------
Give_out_vec_mt::~Give_out_vec_mt()
{
	delete [] H;
	delete [] E;
	delete [] impedance;
	delete [] rho;

	delete [] H_1d;
	delete [] E_1d;

}
//------------------------------------------------------------------------------
//-- Calculates a one-dimensional field (E_1d, H_1d) at a point near the ground (z=0) --
//------------------------------------------------------------------------------
void Give_out_vec_mt::Compute_1d_field()
{
	double uxs,uys,uzs,uxc,uyc,uzc;
	double alpha, w, mu, z;
	int i, j;

	w = d->nu*2.0*PI;
	alpha = d->alfa;
	mu = MU_0;

	ifstream infb,infe;
	infb.open("b2d");
	if(!infb)
	{
		logfile<<"Error in open file "<<"b2d"<<endl;
		exit(1);
	}
	infe.open("e2d");
	if(!infe)
	{
		logfile<<"Error in open file "<<"e2d"<<endl;
		exit(1);
	}
	
	for (j=0; j<d->n_pointres; j++)
	{
		infe>>uxs>>uys>>uzs>>uxc>>uyc>>uzc;
		E_1d[j][0] = std::complex<double>(uxs,uxc);
		E_1d[j][1] = std::complex<double>(uys,uyc);

		infb>>uxs>>uys>>uzs>>uxc>>uyc>>uzc;
		H_1d[j][0] = std::complex<double>(uxs/mu,uxc/mu);
		H_1d[j][1] = std::complex<double>(uys/mu,uyc/mu);
	}

	infb.close();
	infb.clear();
	infe.close();
	infe.clear();
}
//-----------------------------------------------------------------------------
//---Issue E, H, impedance, rok, on a mesh of hexahedrons. 
//-----------------------------------------------------------------------------
void Give_out_vec_mt::Give_out_on_hex(bool for_harm_loop)
{
	int i, j;
	
	double w; // cyclic frequency
	const double mu = MU_0;
	w = d->nu*2.0*PI;

	if (GPDocSettings.vfem_direct_output)
	{
		OutputVect3d output(0,0,0,d->kuzlov, d->kpar, tmap->n, d->n_pointresB, d->pointresB, d->xyz, d->nver, tmap->ed);
		double *vsin=NULL;
		double *vcos=NULL;
		double *resultSin=NULL;
		double *resultCos=NULL;

		if((vsin = new double[tmap->n])==NULL) Memory_allocation_error("vsin", "");
		if((vcos = new double[tmap->n])==NULL) Memory_allocation_error("vcos", "");

		for(i=0; i<tmap->n; i++)
		{
			if(i<tmap->n_c)
			{
				vsin[i] = v3dat[i*2];
				vcos[i] = v3dat[i*2+1];
			}
			else
			{
				vsin[i] = vcos[i] = 0;
				for (j=tmap->ig_t[i]; j<=tmap->ig_t[i+1]-1; j++)
				{
					vsin[i] += v3dat[tmap->jg_t[j]*2]*tmap->gg_t[j];
					vcos[i] += v3dat[tmap->jg_t[j]*2+1]*tmap->gg_t[j];
				}
			}
		}

		if((resultSin = new double[d->n_pointresE])==NULL) Memory_allocation_error("resultSin", "");
		if((resultCos = new double[d->n_pointresE])==NULL) Memory_allocation_error("resultCos", "");

		//Ex
		output.OutputFieldAtReceivers(vsin, resultSin, 0);
		output.OutputFieldAtReceivers(vcos, resultCos, 0);
		for(j=0; j<d->n_pointresE; j++)
		{
			this->E[j][0] = std::complex<double>(w*resultCos[j], -w*resultSin[j]);
			if (!for_harm_loop && d->tasktype) this->E[j][0] += E_1d[j][0];
		}

		//Ey
		output.OutputFieldAtReceivers(vsin, resultSin, 1);
		output.OutputFieldAtReceivers(vcos, resultCos, 1);
		for(j=0; j<d->n_pointresE; j++)
		{
			this->E[j][1] = std::complex<double>(w*resultCos[j], -w*resultSin[j]);
			if (!for_harm_loop && d->tasktype) this->E[j][1] += E_1d[j][1];
		}

		//Ez
		output.OutputFieldAtReceivers(vsin, resultSin, 2);
		output.OutputFieldAtReceivers(vcos, resultCos, 2);
		for(j=0; j<d->n_pointresE; j++)
		{
			this->E[j][2] = std::complex<double>(w*resultCos[j], -w*resultSin[j]);
		}

		delete resultSin; resultSin = NULL;
		delete resultCos; resultCos = NULL;

		if((resultSin = new double[d->n_pointresB])==NULL) Memory_allocation_error("resultSin", "");
		if((resultCos = new double[d->n_pointresB])==NULL) Memory_allocation_error("resultCos", "");

		//Bz
		output.OutputFieldAtReceivers(vsin, resultSin, 3);
		output.OutputFieldAtReceivers(vcos, resultCos, 3);
		for(j=0; j<d->n_pointresB; j++)
		{
			this->H[j][2] = std::complex<double>(resultSin[j]/mu, resultCos[j]/mu);
		}

		//Bx
		output.OutputFieldAtReceivers(vsin, resultSin, 4);
		output.OutputFieldAtReceivers(vcos, resultCos, 4);
		for(j=0; j<d->n_pointresB; j++)
		{
			this->H[j][0] = std::complex<double>(resultSin[j]/mu, resultCos[j]/mu);
		}

		//By
		output.OutputFieldAtReceivers(vsin, resultSin, 5);
		output.OutputFieldAtReceivers(vcos, resultCos, 5);
		for(j=0; j<d->n_pointresB; j++)
		{
			this->H[j][1] = std::complex<double>(resultSin[j]/mu, resultCos[j]/mu);
		}

		delete resultSin; resultSin = NULL;
		delete resultCos; resultCos = NULL;

		delete vsin;
		delete vcos;
	}
	else
	{
		if(d->n_pointresB)
		{
			resultantB->ValueType=vtRotzASin;
			resultantB->Output(0);
			resultantB->ValueType=vtRotzACos;
			resultantB->Output(0);
			resultantB->ValueType=vtRotxASin;
			resultantB->Output(0);
			resultantB->ValueType=vtRotxACos;
			resultantB->Output(0);
			resultantB->ValueType=vtRotyASin;
			resultantB->Output(0);
			resultantB->ValueType=vtRotyACos;
			resultantB->Output(0);
		}
		if (d->n_pointresE)
		{
			resultantA->ValueType=vtAzSin;
			resultantA->Output(0);
			resultantA->ValueType=vtAzCos;
			resultantA->Output(0);
			resultantA->ValueType=vtAxSin;
			resultantA->Output(0);
			resultantA->ValueType=vtAxCos;
			resultantA->Output(0);
			resultantA->ValueType=vtAySin;
			resultantA->Output(0);
			resultantA->ValueType=vtAyCos;
			resultantA->Output(0);
		}
	}

	if(!d->tasktype && for_harm_loop)
	{
		for(j=0; j<d->n_pointresB; j++)
		{

			std::complex<double> Z;
			double e_re_abs, e_im_abs, h_re_abs, h_im_abs;
			std::complex<double> e_abs, h_abs;

			if(!d->fdirect)
			{
				this->H[j][0] += H_1d[j][0];
				this->H[j][1] += H_1d[j][1];
				this->E[j][0] += E_1d[j][0];
				this->E[j][1] += E_1d[j][1];
			}

			e_re_abs = sqrt(pow(real(E[j][0]),2.0) + pow(real(E[j][1]),2.0) + pow(real(E[j][2]),2.0));
			h_re_abs = sqrt(pow(real(H[j][0]),2.0) + pow(real(H[j][1]),2.0) + pow(real(H[j][2]),2.0));
			e_im_abs = sqrt(pow(imag(E[j][0]),2.0) + pow(imag(E[j][1]),2.0) + pow(imag(E[j][2]),2.0));
			h_im_abs = sqrt(pow(imag(H[j][0]),2.0) + pow(imag(H[j][1]),2.0) + pow(imag(H[j][2]),2.0));

			e_abs = std::complex<double>(e_re_abs, e_im_abs);
			h_abs = std::complex<double>(h_re_abs, h_im_abs);

			Z = e_abs/h_abs;

			rho[j] = pow(abs(e_abs)/abs(h_abs), 2.0)/(w*mu);
			impedance[j] = abs(Z);
		}
		Write_result_to_files();
	}
}
//-----------------------------------------------------------------------------
//-------------   writing results to files   --------------------------------
//-----------------------------------------------------------------------------
void Give_out_vec_mt::Write_result_to_files()
{
	long j;
	FILE *ex_c, *ex_s, *ey_c, *ey_s, *ez_c, *ez_s; 
	FILE *hx_c, *hx_s, *hy_c, *hy_s, *hz_c, *hz_s; 
	FILE *rok;
	FILE *f_impedance;


	// Normal field
	ofstream out;

	out.open("normal_field");

	if(!d->fdirect)
	{
		out << "Ex_s:\t" << real(E_1d[0][0]) << '\n';
		out << "Ex_c:\t" << imag(E_1d[0][0]) << '\n';
		out << "Ey_s:\t" << real(E_1d[0][1]) << '\n';
		out << "Ey_c:\t" << imag(E_1d[0][1]) << '\n';

		out << "Hx_s:\t" << real(H_1d[0][0]) << '\n';
		out << "Hx_c:\t" << imag(H_1d[0][0]) << '\n';
		out << "Hy_s:\t" << real(H_1d[0][1]) << '\n';
		out << "Hy_c:\t" << imag(H_1d[0][1]) << '\n';
	}
	else
	{
		out << "Ex_s:\t" << 0.0 << '\n';
		out << "Ex_c:\t" << 0.0 << '\n';
		out << "Ey_s:\t" << 0.0 << '\n';
		out << "Ey_c:\t" << 0.0 << '\n';

		out << "Hx_s:\t" << 0.0 << '\n';
		out << "Hx_c:\t" << 0.0 << '\n';
		out << "Hy_s:\t" << 0.0 << '\n';
		out << "Hy_c:\t" << 0.0 << '\n';
	}

	out.close();
	out.clear();


	///----------

	ex_c = fopen("ex_c", "w"); ex_s = fopen("ex_s", "w");
	ey_c = fopen("ey_c", "w"); ey_s = fopen("ey_s", "w");
	ez_c = fopen("ez_c", "w"); ez_s = fopen("ez_s", "w");
	hx_c = fopen("hx_c", "w"); hx_s = fopen("hx_s", "w");
	hy_c = fopen("hy_c", "w"); hy_s = fopen("hy_s", "w");
	hz_c = fopen("hz_c", "w"); hz_s = fopen("hz_s", "w");
	rok = fopen("rok", "w");
	f_impedance = fopen("impedance", "w");

	fprintf(ex_c, "\n\n");
	fprintf(ex_s, "\n\n");
	fprintf(ey_c, "\n\n");
	fprintf(ey_s, "\n\n");
	fprintf(ez_c, "\n\n");
	fprintf(ez_s, "\n\n");
	fprintf(hx_c, "\n\n");
	fprintf(hx_s, "\n\n");
	fprintf(hy_c, "\n\n");
	fprintf(hy_s, "\n\n");	
	fprintf(hz_c, "\n\n");
	fprintf(hz_s, "\n\n");
	fprintf(rok, "\n\n");
	fprintf(f_impedance, "\n\n");

	for(j=0; j<d->n_pointresB; j++)
	{
		double x;
		x = d->pointresB[j][0];
		fprintf(ex_c, "%g %g\n", x, imag(E[j][0]));
		fprintf(ex_s, "%g %g\n", x, real(E[j][0]));
		fprintf(ey_c, "%g %g\n", x, imag(E[j][1]));
		fprintf(ey_s, "%g %g\n", x, real(E[j][1]));
		fprintf(ez_c, "%g %g\n", x, imag(E[j][2]));
		fprintf(ez_s, "%g %g\n", x, real(E[j][2]));
		fprintf(hx_c, "%g %g\n", x, imag(H[j][0]));
		fprintf(hx_s, "%g %g\n", x, real(H[j][0]));
		fprintf(hy_c, "%g %g\n", x, imag(H[j][1]));
		fprintf(hy_s, "%g %g\n", x, real(H[j][1]));	
		fprintf(hz_c, "%g %g\n", x, imag(H[j][2]));
		fprintf(hz_s, "%g %g\n", x, real(H[j][2]));
		fprintf(rok, "%g %g\n", x, rho[j]);
		fprintf(f_impedance, "%g %g\n", x, impedance[j]);
	}

	fclose(ex_c); fclose(ex_s); fclose(ey_c); fclose(ey_s); fclose(ez_c); fclose(ez_s);
	fclose(hx_c); fclose(hx_s); fclose(hy_c); fclose(hy_s); fclose(hz_c); fclose(hz_s);
	fclose(rok);
	fclose(f_impedance);
}
//-----------------------------------------------------------------------------
void Give_out_vec_mt::Write_B_to_files_for_harm_loop(int StartType,int fdirect)
{
	int i;
	double sin_x, sin_y, sin_z, cos_x, cos_y, cos_z; 
	std::complex <double> (*B2d)[3];
	ifstream inf;
	ofstream outf;

	if(StartType!=2)
	{
		if(!fdirect)
		{
			if((B2d = new std::complex<double>[d->n_pointresB][3])==0)
				Memory_allocation_error("B", "Write_result_to_files_for_harm_loop");

			inf.open("b2d");
			for (i=0; i<d->n_pointresB; i++)
			{
				inf>>sin_x>>sin_y>>sin_z>>cos_x>>cos_y>>cos_z;
				B2d[i][0]=std::complex<double>(sin_x, cos_x); 
				B2d[i][1]=std::complex<double>(sin_y, cos_y); 
				B2d[i][2]=std::complex<double>(sin_z, cos_z); 
			}
			inf.close();
			inf.clear();

			outf.open("b3d");
			outf<<setiosflags(ios_base::scientific)<<setprecision(14);
			for (i=0; i<d->n_pointresB; i++)
				outf<<real(B2d[i][0])+MU_0*real(H[i][0])<<' '
				<<real(B2d[i][1])+MU_0*real(H[i][1])<<' '
				<<real(B2d[i][2])+MU_0*real(H[i][2])<<' '
				<<imag(B2d[i][0])+MU_0*imag(H[i][0])<<' '
				<<imag(B2d[i][1])+MU_0*imag(H[i][1])<<' '
				<<imag(B2d[i][2])+MU_0*imag(H[i][2])<<endl; 
			outf.close();

			if(B2d) {delete [] B2d; B2d=NULL;}
		}
		else
		{
			outf.open("b3d");
			outf<<setiosflags(ios_base::scientific)<<setprecision(14);
			for (i=0; i<d->n_pointresB; i++)
				outf<<MU_0*real(H[i][0])<<' '
				<<MU_0*real(H[i][1])<<' '
				<<MU_0*real(H[i][2])<<' '
				<<MU_0*imag(H[i][0])<<' '
				<<MU_0*imag(H[i][1])<<' '
				<<MU_0*imag(H[i][2])<<endl; 
			outf.close();
		}
	}
	else
	{
		outf.open("b3d_anom");
		outf<<setiosflags(ios_base::scientific)<<setprecision(14);
		for (i=0; i<d->n_pointresB; i++)
			outf<<MU_0*real(H[i][0])<<' '
			<<MU_0*real(H[i][1])<<' '
			<<MU_0*real(H[i][2])<<' '
			<<MU_0*imag(H[i][0])<<' '
			<<MU_0*imag(H[i][1])<<' '
			<<MU_0*imag(H[i][2])<<endl; 
		outf.close();
	}
}
//------------------------------------------------------------------------
void Give_out_vec_mt::Write_E_to_files_for_harm_loop(int StartType,int fdirect)
{
	int i;
	double sin_x, sin_y, sin_z, cos_x, cos_y, cos_z; 
	std::complex <double> (*E2d)[3];
	ifstream inf;
	ofstream outf;

	if((E2d = new std::complex<double>[d->n_pointresE][3])==0)
		Memory_allocation_error("E", "Write_result_to_files_for_harm_loop");

	if(StartType!=2)
	{
		if(!fdirect)
		{
			inf.open("e2d");
			for (i=0; i<d->n_pointresE; i++)
			{
				inf>>sin_x>>sin_y>>sin_z>>cos_x>>cos_y>>cos_z;
				E2d[i][0]=std::complex<double>(sin_x, cos_x); 
				E2d[i][1]=std::complex<double>(sin_y, cos_y); 
				E2d[i][2]=std::complex<double>(sin_z, cos_z); 
			}
			inf.close();
			inf.clear();

			outf.open("e3d");
			outf<<setiosflags(ios_base::scientific)<<setprecision(14);
			for (i=0; i<d->n_pointresE; i++)
				outf<<real(E2d[i][0])+real(E[i][0])<<' '
				<<real(E2d[i][1])+real(E[i][1])<<' '
				<<real(E2d[i][2])+real(E[i][2])<<' '
				<<imag(E2d[i][0])+imag(E[i][0])<<' '
				<<imag(E2d[i][1])+imag(E[i][1])<<' '
				<<imag(E2d[i][2])+imag(E[i][2])<<endl; 
			outf.close();

			if(E2d) {delete [] E2d; E2d=NULL;}
		}
		else
		{
			outf.open("e3d");
			outf<<setiosflags(ios_base::scientific)<<setprecision(14);
			for (i=0; i<d->n_pointresE; i++)
				outf<<real(E2d[i][0])+real(E[i][0])<<' '
				<<real(E[i][1])<<' '
				<<real(E[i][2])<<' '
				<<imag(E[i][0])<<' '
				<<imag(E[i][1])<<' '
				<<imag(E[i][2])<<endl; 
			outf.close();
		}
	}
	else
	{
		outf.open("e3d_anom");
		outf<<setiosflags(ios_base::scientific)<<setprecision(14);
		for (i=0; i<d->n_pointresE; i++)
			outf<<real(E[i][0])<<' '
			<<real(E[i][1])<<' '
			<<real(E[i][2])<<' '
			<<imag(E[i][0])<<' '
			<<imag(E[i][1])<<' '
			<<imag(E[i][2])<<endl; 
		outf.close();
	}
}
//------------------------------------------------------------------------
void Give_out_vec_mt::Write_result_to_files_for_harm_loop()
{
	int i;
	double sin_x, sin_y, sin_z, cos_x, cos_y, cos_z; 
	std::complex <double> (*B2d)[3];
	std::complex <double> (*E2d)[3];
	ifstream inf;
	ofstream outf;

	if((E2d = new std::complex<double>[d->n_pointresE][3])==0)
		Memory_allocation_error("E", "Write_result_to_files_for_harm_loop");

	inf.open("e2d");
	for (i=0; i<d->n_pointresE; i++)
	{
		inf>>sin_x>>sin_y>>sin_z>>cos_x>>cos_y>>cos_z;
		E2d[i][0]=std::complex<double>(sin_x, cos_x); 
		E2d[i][1]=std::complex<double>(sin_y, cos_y); 
		E2d[i][2]=std::complex<double>(sin_z, cos_z); 
	}
	inf.close();
	inf.clear();

	if((B2d = new std::complex<double>[d->n_pointresB][3])==0)
		Memory_allocation_error("B", "Write_result_to_files_for_harm_loop");
		
	inf.open("b2d");
	for (i=0; i<d->n_pointresB; i++)
	{
		inf>>sin_x>>sin_y>>sin_z>>cos_x>>cos_y>>cos_z;
		B2d[i][0]=std::complex<double>(sin_x, cos_x); 
		B2d[i][1]=std::complex<double>(sin_y, cos_y); 
		B2d[i][2]=std::complex<double>(sin_z, cos_z); 
	}
	inf.close();
	inf.clear();	

	outf.open("e3d");
	for (i=0; i<d->n_pointresE; i++)
		outf<<real(E2d[i][0])+real(E[i][0])<<' '
			<<real(E2d[i][1])+real(E[i][1])<<' '
			<<real(E2d[i][2])+real(E[i][2])<<' '
			<<imag(E2d[i][0])+imag(E[i][0])<<' '
			<<imag(E2d[i][1])+imag(E[i][1])<<' '
			<<imag(E2d[i][2])+imag(E[i][2])<<endl; 
	outf.close();
	
	outf.open("b3d");
	for (i=0; i<d->n_pointresB; i++)
		outf<<real(B2d[i][0])+MU_0*real(H[i][0])<<' '
			<<real(B2d[i][1])+MU_0*real(H[i][1])<<' '
			<<real(B2d[i][2])+MU_0*real(H[i][2])<<' '
			<<imag(B2d[i][0])+MU_0*imag(H[i][0])<<' '
			<<imag(B2d[i][1])+MU_0*imag(H[i][1])<<' '
			<<imag(B2d[i][2])+MU_0*imag(H[i][2])<<endl; 
	outf.close();


	if(B2d) {delete [] B2d; B2d=NULL;}
	if(E2d) {delete [] E2d; E2d=NULL;}
}
//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
int Give_out_vec_mt::Read_1d_field(char *fname_in, double &Ex_s, double &Ex_c, double &Ey_s, double &Ey_c, 
				  double &Hx_s, double &Hx_c, double &Hy_s, double &Hy_c)
{
	ifstream f_in;
	char buffer[10];

	f_in.open(fname_in);
	f_in >> buffer >> Ex_s;
	f_in >> buffer >> Ex_c;
	f_in >> buffer >> Ey_s;
	f_in >> buffer >> Ey_c;
	f_in >> buffer >> Hx_s;
	f_in >> buffer >> Hx_c;
	f_in >> buffer >> Hy_s;
	f_in >> buffer >> Hy_c;
	f_in.close();

	return 0;
}
//--------------------------------------------------------------------------------------
//differentiation of the solution in time (according to the 3-layer scheme) at the point t from [t_j2, t_j]
//---------------------------------------------------------------------------------------
double Give_out_vec_mt::dA_dt(double t, double u_j, double u_j1, double u_j2,
			 double dt, double dt0, double dt1, double t_j, double t_j1, double t_j2)
{
	double du_dt;

	if(t < t_j2) t = t_j2;
	if(t > t_j)  t = t_j;

	du_dt = u_j2*(2.0*t - t_j  - t_j1)/(dt1*dt) - 
		u_j1*(2.0*t - t_j  - t_j2)/(dt1*dt0) + 
		u_j*(2.0*t - t_j1 - t_j2)/(dt*dt0);

	return du_dt;
}
//------------------------------------------------------------------------
int Give_out_vec_mt::GetNumberOfNodes()
{
	return d->kuzlov;
}
int Give_out_vec_mt::GetNumberOfElements()
{
	return d->kpar;
}
int Give_out_vec_mt::GetElementNodesNumber()
{
	return 8;
}
const pv::Point3D Give_out_vec_mt::GetNode(const int& i_node)
{
	pv::Point3D Point;
	Point.x()=d->xyz[i_node][0];
	Point.y()=d->xyz[i_node][1];
	Point.z()=d->xyz[i_node][2];
	return Point;
}
const pv::Point3D Give_out_vec_mt::GetNodeTrue(const int& i_node)
{
	pv::Point3D Point;
	Point.x()=d->xyzt[i_node][0];
	Point.y()=d->xyzt[i_node][1];
	Point.z()=d->xyzt[i_node][2];
	return Point;
}
int Give_out_vec_mt::GetNodeNumberOnElement(const int& i_element, const int& i_node)
{
	return d->nver[i_element][i_node];
}
int Give_out_vec_mt::GetElementMaterial(const int& i_element)
{
	return d->nvkat[i_element];
}
int Give_out_vec_mt::GetTypeOfElement(const int& i_element)
{
	return d->nver[i_element][13];
}
double Give_out_vec_mt::GetValueInElementCenter(const int& i_element, const Res3DValueType& r_type)
{
	int i;
	double f[12];
	double x[8],y[8],z[8];
	double in[3],out[3],FieldOnElem;
	int isReal,isNoRot;

	FieldOnElem=0;

	isReal=r_type==vtAxSin ||r_type==vtAySin || r_type==vtAzSin || 
		r_type==vtRotxASin ||r_type==vtRotyASin || r_type==vtRotzASin;

	isNoRot=r_type==vtAxSin || r_type==vtAySin || r_type==vtAzSin || 
		r_type==vtAxCos ||r_type==vtAyCos || r_type==vtAzCos;
	
	for (i=0; i<12; i++)f[i] = v3dat[2*(tmap->ed[i_element][i])+!isReal];


	for (i=0; i<8; i++){
		x[i] = d->xyz[d->nver[i_element][i]][0];
		y[i] = d->xyz[d->nver[i_element][i]][1];
		z[i] = d->xyz[d->nver[i_element][i]][2];
	}
	
	T_Brick L(x, y, z, d->nver[i_element][13]);

	in[0]=in[1]=in[2]=0.0;

	if(isNoRot){
		L.Calc_value_inside_hex(f,in,out);

		if(r_type==vtAxSin || r_type==vtAxCos)FieldOnElem=out[0];
		if(r_type==vtAySin || r_type==vtAyCos)FieldOnElem=out[1];
		if(r_type==vtAzSin || r_type==vtAzCos)FieldOnElem=out[2];
	}
	else{
		L.Calc_rotor_inside_hex(f,in,out);

		if(r_type==vtRotxASin || r_type==vtRotxACos)FieldOnElem=out[0];
		if(r_type==vtRotyASin || r_type==vtRotyACos)FieldOnElem=out[1];
		if(r_type==vtRotzASin || r_type==vtRotzACos)FieldOnElem=out[2];
	}

	return FieldOnElem;
}
int Give_out_vec_mt::GetNumberOfResPoints(const Res3DValueType& r_type)
{
	return (r_type==vtWithDiscontinuity)? d->n_pointresE : d->n_pointresB;
}
pv::Point3D Give_out_vec_mt::GetResPoint(const Res3DValueType& r_type, const int& i_point)
{
	pv::Point3D Point;
	Point.x()=(r_type==vtWithDiscontinuity)? d->pointresE[i_point][0] : d->pointresB[i_point][0];
	Point.y()=(r_type==vtWithDiscontinuity)? d->pointresE[i_point][1] : d->pointresB[i_point][1];
	Point.z()=(r_type==vtWithDiscontinuity)? d->pointresE[i_point][2] : d->pointresB[i_point][2];
	return Point;
}

int * Give_out_vec_mt::GetPointerToRegular()
{
	return &(d->regular[0]);
}

int Give_out_vec_mt::GetXSize()
{
	return d->N_X;
}

int Give_out_vec_mt::GetYSize()
{
	return d->N_Y;
}

int Give_out_vec_mt::GetZSize()
{
	return d->N_Z;
}

double * Give_out_vec_mt::GetPointerToX()
{
	return &(d->Xcrd[0]);
}

double * Give_out_vec_mt::GetPointerToY()
{
	return &(d->Ycrd[0]);
}

double * Give_out_vec_mt::GetPointerToZ()
{
	return &(d->Zcrd[0]);
}

void Give_out_vec_mt::SaveResult(const Res3DValueType& r_type, const double& r_value, const int& j, const int& i_time)
{
	double for_res_mu=MU_0;
	double for_res_w=d->nu*2.0*PI;


	if(r_type==vtRotxASin){
		std::complex <double> ComplexValue(r_value/for_res_mu,imag(H[j][0]));
		H[j][0]=ComplexValue;
	}
	else if(r_type==vtRotxACos){
		std::complex <double> ComplexValue(real(H[j][0]),r_value/for_res_mu);
		H[j][0]=ComplexValue;
	}
	else if(r_type==vtRotyASin){
		std::complex <double> ComplexValue(r_value/for_res_mu,imag(H[j][1]));
		H[j][1]=ComplexValue;
	}
	else if(r_type==vtRotyACos){
		std::complex <double> ComplexValue(real(H[j][1]),r_value/for_res_mu);
		H[j][1]=ComplexValue;
	}
	else if(r_type==vtRotzASin){
		std::complex <double> ComplexValue(r_value/for_res_mu,imag(H[j][2]));
		H[j][2]=ComplexValue;
	}
	else if(r_type==vtRotzACos){
		std::complex <double> ComplexValue(real(H[j][2]),r_value/for_res_mu);
		H[j][2]=ComplexValue;
	}
	else if(r_type==vtAxSin){
		std::complex <double> ComplexValue(real(E[j][0]),-r_value*for_res_w);
		E[j][0]=ComplexValue;
	}
	else if(r_type==vtAxCos){
		std::complex <double> ComplexValue(r_value*for_res_w,imag(E[j][0]));
		E[j][0]=ComplexValue;
	}
	else if(r_type==vtAySin){
		std::complex <double> ComplexValue(real(E[j][1]),-r_value*for_res_w);
		E[j][1]=ComplexValue;
	}
	else if(r_type==vtAyCos){
		std::complex <double> ComplexValue(r_value*for_res_w,imag(E[j][1]));
		E[j][1]=ComplexValue;
	}
	else if(r_type==vtAzSin){
		std::complex <double> ComplexValue(real(E[j][2]),-r_value*for_res_w);
		E[j][2]=ComplexValue;
	}
	else if(r_type==vtAzCos){
		std::complex <double> ComplexValue(r_value*for_res_w,imag(E[j][2]));
		E[j][2]=ComplexValue;
	}
	else
	{ 
		logfile<<"Unknown Res3DValueType in Give_out_vec_mt::SaveResult"<<endl;
	}
}
