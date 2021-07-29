#include "stdafx.h"
#include "open_close.h"
#include "utils.h"
#include "linalg.h"

char str[1024];

const double PI=3.1415926535897932;
const double MU=4*PI*1e-7;

const int size_f=sizeof(float);

char PathToInstallDir[1024];

int CreateProcessForEXE(char *cmdline, char *workdir)
{
	int retp;
	STARTUPINFOA si;
	PROCESS_INFORMATION pi;
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&pi, sizeof(pi));
	// Start the child process. 
	if (!(retp=CreateProcessA(NULL,  // No module name (use command line). 
		(LPSTR)(const char*)cmdline,// Command line. 
		NULL,				// Process handle not inheritable. 
		NULL,				// Thread handle not inheritable. 
		FALSE,				// Set handle inheritance to FALSE. 
		0/*CREATE_NO_WINDOW*/,	// No creation flags. 
		NULL,				// Use parent's environment block. 
		workdir,			// Use parent's starting directory. 
		&si,				// Pointer to STARTUPINFO structure.
		&pi)))				// Pointer to PROCESS_INFORMATION structure.
	{
		sprintf(charbuffer,"Can't create process for %s, error code %d",cmdline,GetLastError());
		cout<<charbuffer<<endl;
		logfile<<charbuffer<<endl;
		exit(retp);
	}
	// Wait until child process exits.
	WaitForSingleObject(pi.hProcess, INFINITE);
	// Get exit code.
	GetExitCodeProcess(pi.hProcess, (LPDWORD)&retp);
	// Close process and thread handles. 
	CloseHandle(pi.hProcess);
	CloseHandle(pi.hThread);
	return retp;
}

void run(char *cmd)
{
	int retp;
	retp=CreateProcessForEXE(cmd,NULL);
	if(retp)
	{
		cout<<"Error: "<<cmd<<" returned "<<retp<<endl;
		logfile<<"Error: "<<cmd<<" returned "<<retp<<endl;
		exit(retp);
	}
}

int SolveGauss(complex<double> *c_A,complex<double> *c_b,int slae_size,complex<double> *(&L),complex<double> *(&f),complex<double> *(&aa))
{
	const int b = 32;
	int blocksize = b;
	int n=slae_size;

	if(!L){L = new complex<double>[(b*b - b)/2];}
	if(!f){f = new complex<double>[n*blocksize];}
	if(!aa){aa = new complex<double>[n*blocksize];}

	getrf(c_A, slae_size, b, blocksize, L, f, aa);
	SolveL1(c_A, c_b, slae_size);
	return SolveU(c_A, c_b, slae_size);
}

const int n=6;
const int nsq=n*n;

void GetMemmory(vector<vector<double>> &Fields,int nrec,int nfrq)
{
	double zero=0.0;
	int ifrq,irec;
	Fields.resize(nrec);
	for(irec=0;irec<nrec;irec++)
	{
		vector<double> &fi=Fields[irec];
		fi.resize(nfrq);
		for(ifrq=0;ifrq<nfrq;ifrq++){fi[ifrq]=zero;}
	}
}


void GetMemmory(vector<vector<complex<double>>> &Fields,int nrec,int nfrq)
{
	complex<double> zero(0.0,0.0);
	int ifrq,irec;
	Fields.resize(nrec);
	for(irec=0;irec<nrec;irec++)
	{
		vector<complex<double>> &fi=Fields[irec];
		fi.resize(nfrq);
		for(ifrq=0;ifrq<nfrq;ifrq++){fi[ifrq]=zero;}
	}
}

void ReadFrqFile(vector<vector<complex<double>>> &Fields,int nrec,int ifrq,char *fnameR,char *fnameI)
{
	ifstream infr,infi;
	int irec;
	double u,v,tmpd;
	char fnr[256],fni[256];
	sprintf(fnr,"%s%d",fnameR,ifrq);
	sprintf(fni,"%s%d",fnameI,ifrq);
	StopIfErrorReturn(OpenInputFile(infr,fnr),"OpenInputFile");
	StopIfErrorReturn(OpenInputFile(infi,fni),"OpenInputFile");
	for(irec=0;irec<nrec;irec++)
	{
		infr>>tmpd>>u;
		infi>>tmpd>>v;
		Fields[irec][ifrq]=complex<double>(u,v);
	}
	StopIfErrorReturn(CloseInputFile(infr,fnr),"CloseInputFile");
	StopIfErrorReturn(CloseInputFile(infi,fni),"CloseInputFile");
}

void ReadAllFile(vector<vector<complex<double>>> &Fields,int nrec,int nfrq,char *fnameR,char *fnameI)
{
	ifstream infr,infi;
	int ifrq,irec;
	double u,v,tmpd;
	StopIfErrorReturn(OpenInputFile(infr,fnameR),"OpenInputFile");
	StopIfErrorReturn(OpenInputFile(infi,fnameI),"OpenInputFile");
	Fields.resize(nrec);
	for(irec=0;irec<nrec;irec++)
	{
		vector<complex<double>> &fi=Fields[irec];
		fi.resize(nfrq);
		infr.ignore(1000,'\n');
		infr.ignore(1000,'\n');
		infr.ignore(1000,'\n');
		infr.ignore(1000,'\n');
		infr.ignore(1000,'\n');
		infr.ignore(1000,'\n');
		infi.ignore(1000,'\n');
		infi.ignore(1000,'\n');
		infi.ignore(1000,'\n');
		infi.ignore(1000,'\n');
		infi.ignore(1000,'\n');
		infi.ignore(1000,'\n');
		for(ifrq=0;ifrq<nfrq;ifrq++)
		{
			infr>>tmpd>>tmpd>>tmpd>>u;
			infi>>tmpd>>tmpd>>tmpd>>v;
			fi[ifrq]=complex<double>(u,v);
		}
		if(nfrq)
		{
			infr.ignore(1000,'\n');
			infi.ignore(1000,'\n');
		}
	}
	StopIfErrorReturn(CloseInputFile(infr,fnameR),"CloseInputFile");
	StopIfErrorReturn(CloseInputFile(infi,fnameI),"CloseInputFile");
}

void ReadAllFile(vector<vector<complex<double>>> &Fx,vector<vector<complex<double>>> &Fy,vector<vector<complex<double>>> &Fz,
				 int nrec,int nfrq,char *fname,char c,int ia,bool fdmu)
{
	ifstream inf;
	int ifrq,irec;
	double u[3],v[3];
	char buf[256];
	
	Fx.resize(nrec);
	Fy.resize(nrec);
	Fz.resize(nrec);
	for(irec=0;irec<nrec;irec++)
	{
		Fx[irec].resize(nfrq);
		Fy[irec].resize(nfrq);
		Fz[irec].resize(nfrq);
	}
	for(ifrq=0;ifrq<nfrq;ifrq++)
	{
		if(!ia)
			sprintf(buf,"%s%c_s0f%d",fname,c,ifrq);
		else
			sprintf(buf,"%s%c_s0f%d_%d",fname,c,ifrq,ia);

		inf.open(buf);
		if(!inf)
		{
			logfile<<"Error in open file "<<buf<<endl;
			exit(1);
		}

		for(irec=0;irec<nrec;irec++)
		{
			inf>>u[0]>>u[1]>>u[2]>>v[0]>>v[1]>>v[2];

			if(fdmu)
			{
				u[0]/=MU;
				u[1]/=MU;
				u[2]/=MU;

				v[0]/=MU;
				v[1]/=MU;
				v[2]/=MU;
			}
			
			Fx[irec][ifrq]=complex<double>(u[0],v[0]);
			Fy[irec][ifrq]=complex<double>(u[1],v[1]);
			Fz[irec][ifrq]=complex<double>(u[2],v[2]);
		}

		inf.close();
		inf.clear();
	}
}

void ReadAllFile(vector<vector<complex<double>>> &Fx,vector<vector<complex<double>>> &Fy,
				 int nrec,int nfrq,char *fname,char c,int ia,bool fdmu)
{
	ifstream inf;
	int ifrq,irec;
	double u[3],v[3];
	char buf[256];
	
	Fx.resize(nrec);
	Fy.resize(nrec);
	for(irec=0;irec<nrec;irec++)
	{
		Fx[irec].resize(nfrq);
		Fy[irec].resize(nfrq);
	}
	for(ifrq=0;ifrq<nfrq;ifrq++)
	{
		if(!ia)
			sprintf(buf,"%s%c_s0f%d",fname,c,ifrq);
		else
			sprintf(buf,"%s%c_s0f%d_%d",fname,c,ifrq,ia);

		inf.open(buf);
		if(!inf)
		{
			logfile<<"Error in open file "<<buf<<endl;
			exit(1);
		}

		for(irec=0;irec<nrec;irec++)
		{
			inf>>u[0]>>u[1]>>u[2]>>v[0]>>v[1]>>v[2];

			if(fdmu)
			{
				u[0]/=MU;
				u[1]/=MU;
				u[2]/=MU;

				v[0]/=MU;
				v[1]/=MU;
				v[2]/=MU;
			}
			
			Fx[irec][ifrq]=complex<double>(u[0],v[0]);
			Fy[irec][ifrq]=complex<double>(u[1],v[1]);
		}

		inf.close();
		inf.clear();
	}
}

void WriteAllFileReal(vector<vector<complex<double>>> &Field,int nrec,int nfrq,vector<double> &sfreq,
				  vector<double> &Xr,vector<double> &Yr,vector<double> &Zr,const char *fname,const char *Leg)
{
	const int bufSize = 10;
	char buffer[bufSize];
	int i,f;
	ofstream outr;
	outr.open(fname);
	outr<<setprecision(7)<<setiosflags(ios_base::scientific);
	for(i=0;i<nrec;i++)
	{
		double r0 = 0.0;
		outr<<(int)Xr[i]<<" "<<(int)Yr[i]<<'\n';
		outr<<Xr[i]<<" "<<Yr[i]<<" "<<Zr[i]<<'\n';
		outr<<"1"<<'\n';
		outr<<"1 1 1 1"<<'\n';
		outr<<"1"<<'\n';
		strcpy_s(buffer, bufSize,Leg);
		outr<<"    v(Hz)    "<<buffer<<'\n';
		for(f=0;f<nfrq;f++)
		{
			outr<<setw(15)<<1e-3/sfreq[f]<<" "
				<<setw(15)<<Field[i][f].real()<<" "
				<<setw(15)<<r0<<" "
				<<setw(15)<<Field[i][f].real()<<" "<<'\n';
		}
	}
	outr.close();
	outr.clear();
}

void WriteAllFileImag(vector<vector<complex<double>>> &Field,int nrec,int nfrq,vector<double> &sfreq,
				  vector<double> &Xr,vector<double> &Yr,vector<double> &Zr,const char *fname,const char *Leg)
{
	const int bufSize = 10;
	char buffer[bufSize];
	int i,f;
	ofstream outr;
	outr.open(fname);
	for(i=0;i<nrec;i++)
	{
		double r0 = 0.0;
		outr<<(int)Xr[i]<<" "<<(int)Yr[i]<<'\n';
		outr<<Xr[i]<<" "<<Yr[i]<<" "<<Zr[i]<<'\n';
		outr<<"1"<<'\n';
		outr<<"1 1 1 1"<<'\n';
		outr<<"1"<<'\n';
		strcpy_s(buffer, bufSize,Leg);
		outr<<"    v(Hz)    "<<buffer<<'\n';
		for(f=0;f<nfrq;f++)
		{
			outr<<setw(15)<<1e-3/sfreq[f]<<" "
				<<setw(15)<<Field[i][f].imag()<<" "
				<<setw(15)<<r0<<" "
				<<setw(15)<<Field[i][f].imag()<<" "<<'\n';
		}
	}
	outr.close();
	outr.clear();
}

void MultMatrixes(complex<double> *A,complex<double> *B,complex<double> *R,int n,int m)
{
	int i,j,k,in,im,kn;
	for(i=0;i<n;i++)
	{
		in=i*n;
		im=i*m;
		for(j=0;j<n;j++)
		{
			R[in+j]=0.0;
			for(k=0;k<m;k++)
			{
				kn=k*n;
				R[in+j]+=A[im+k]*B[kn+j];
			}
		}
	}
}

void MultMatrixVector(complex<double> *A,complex<double> *b,complex<double> *r,int n)
{
	int i,j,in;
	for(i=0;i<n;i++)
	{
		in=i*n;
		r[i]=0.0;
		for(j=0;j<n;j++)
		{
			r[i]+=A[in+j]*b[j];
		}
	}
}

void ReadZTensor(int nrec,int nfreq,vector<vector<double>> &Zxx,vector<vector<double>> &Zxy,
				 vector<vector<double>> &Zyx,vector<vector<double>> &Zyy,
				 char *Fxx,char *Fxy,char *Fyx,char *Fyy)
{
	int i,j;
	float Uxx,Uxy,Uyx,Uyy;
	ifstream infxx,infxy,infyx,infyy;

	StopIfErrorReturn(OpenInputFile(infxx,Fxx,ios::binary),"OpenInputFile");
	StopIfErrorReturn(OpenInputFile(infxy,Fxy,ios::binary),"OpenInputFile");
	StopIfErrorReturn(OpenInputFile(infyx,Fyx,ios::binary),"OpenInputFile");
	StopIfErrorReturn(OpenInputFile(infyy,Fyy,ios::binary),"OpenInputFile");

	for(i=0;i<nrec;i++)
	{
		for(j=0;j<nfreq;j++)
		{
			infxx.read((char *)&Uxx,size_f);
			infxy.read((char *)&Uxy,size_f);
			infyx.read((char *)&Uyx,size_f);
			infyy.read((char *)&Uyy,size_f);
			Zxx[i][j]=Uxx;
			Zxy[i][j]=Uxy;
			Zyx[i][j]=Uyx;
			Zyy[i][j]=Uyy;
		}
	}

	StopIfErrorReturn(CloseInputFile(infxx,Fxx),"CloseInputFile");
	StopIfErrorReturn(CloseInputFile(infxy,Fxy),"CloseInputFile");
	StopIfErrorReturn(CloseInputFile(infyx,Fyx),"CloseInputFile");
	StopIfErrorReturn(CloseInputFile(infyy,Fyy),"CloseInputFile");
}

void ReadZTensor(int nrec,int nfreq,vector<vector<double>> &Zxx,vector<vector<double>> &Zyy,char *Fxx,char *Fyy)
{
	int i,j;
	float Uxx,Uyy;
	ifstream infxx,infyy;

	StopIfErrorReturn(OpenInputFile(infxx,Fxx,ios::binary),"OpenInputFile");
	StopIfErrorReturn(OpenInputFile(infyy,Fyy,ios::binary),"OpenInputFile");

	for(i=0;i<nrec;i++)
	{
		for(j=0;j<nfreq;j++)
		{
			infxx.read((char *)&Uxx,size_f);
			infyy.read((char *)&Uyy,size_f);
			Zxx[i][j]=Uxx;
			Zyy[i][j]=Uyy;
		}
	}

	StopIfErrorReturn(CloseInputFile(infxx,Fxx),"CloseInputFile");
	StopIfErrorReturn(CloseInputFile(infyy,Fyy),"CloseInputFile");
}

void ReadZTensor(int nrec,int nfreq,vector<vector<double>> &Zxx,char *Fxx)
{
	int i,j;
	float Uxx;
	ifstream infxx;

	StopIfErrorReturn(OpenInputFile(infxx,Fxx,ios::binary),"OpenInputFile");

	for(i=0;i<nrec;i++)
	{
		for(j=0;j<nfreq;j++)
		{
			infxx.read((char *)&Uxx,size_f);
			Zxx[i][j]=Uxx;
		}
	}

	StopIfErrorReturn(CloseInputFile(infxx,Fxx),"CloseInputFile");
}

void WriteZTensor(int nrec,vector<double> &Xr,vector<double> &Yr,vector<double> &Zr,int nfreq,vector<double> &sfreq,vector<vector<double>> &Zxx,
				  vector<vector<double>> &Zxy,vector<vector<double>> &Zyx,vector<vector<double>> &Zyy,char *Fname,char *Leg)
{
	int i,j;
	ofstream ofp;
	for(i=0;i<nrec;i++)
	{
		sprintf(str,"%s.%d.txt",Fname,i+1);
		StopIfErrorReturn(OpenOutputFile(ofp,str),"OpenOutputFile");
		ofp<<scientific<<setprecision(7);
		ofp<<0.0<<'\t'<<0.0<<'\t'<<0.0<<'\n';
		ofp<<Xr[i]<<'\t'<<Yr[i]<<'\t'<<Zr[i]<<'\n';
		ofp<<'\n';
		ofp<<Leg<<'\n';
		for(j=0;j<nfreq;j++)
		{
			ofp<<sfreq[j]<<'\t'<<Zxx[i][j]<<'\t'<<Zxy[i][j]<<'\t'<<Zyy[i][j]<<'\t'<<Zyx[i][j]<<'\n';
		}
		StopIfErrorReturn(CloseOutputFile(ofp,str),"CloseOutputFile");
	}
}

void WriteZTensor1(int nrec,vector<double> &Xr,vector<double> &Yr,vector<double> &Zr,int nfreq,vector<double> &sfreq,vector<vector<double>> &Zxx,
				  vector<vector<double>> &Zxy,vector<vector<double>> &Zyx,vector<vector<double>> &Zyy,char *Fname,char *Leg)
{
	int i,j;
	ofstream ofp;
	for(i=0;i<nrec;i++)
	{
		sprintf(str,"%s.%d.txt",Fname,i+1);
		StopIfErrorReturn(OpenOutputFile(ofp,str),"OpenOutputFile");
		ofp<<scientific<<setprecision(7);
		ofp<<0.0<<'\t'<<0.0<<'\t'<<0.0<<'\n';
		ofp<<Xr[i]<<'\t'<<Yr[i]<<'\t'<<Zr[i]<<'\n';
		ofp<<'\n';
		ofp<<Leg<<'\n';
		for(j=0;j<nfreq;j++)
		{
			ofp<<sfreq[j]<<'\t'<<Zxx[i][j]<<'\t'<<Zxy[i][j]-360.0<<'\t'<<Zyy[i][j]<<'\t'<<Zyx[i][j]<<'\n';
		}
		StopIfErrorReturn(CloseOutputFile(ofp,str),"CloseOutputFile");
	}
}

void WriteZTensor(int nrec,vector<double> &Xr,vector<double> &Yr,vector<double> &Zr,int nfreq,vector<double> &sfreq,
				  vector<vector<double>> &Zxx,vector<vector<double>> &Zyy,char *Fname,char *Leg)
{
	int i,j;
	ofstream ofp;
	for(i=0;i<nrec;i++)
	{
		sprintf(str,"%s.%d.txt",Fname,i+1);
		StopIfErrorReturn(OpenOutputFile(ofp,str),"OpenOutputFile");
		ofp<<scientific<<setprecision(7);
		ofp<<0.0<<'\t'<<0.0<<'\t'<<0.0<<'\n';
		ofp<<Xr[i]<<'\t'<<Yr[i]<<'\t'<<Zr[i]<<'\n';
		ofp<<'\n';
		ofp<<Leg<<'\n';
		for(j=0;j<nfreq;j++)
		{
			ofp<<sfreq[j]<<'\t'<<Zxx[i][j]<<'\t'<<Zyy[i][j]<<'\n';
		}
		StopIfErrorReturn(CloseOutputFile(ofp,str),"CloseOutputFile");
	}
}

void WriteZTensor(int nrec,vector<double> &Xr,vector<double> &Yr,vector<double> &Zr,int nfreq,vector<double> &sfreq,
				  vector<vector<double>> &Zxx,char *Fname,char *Leg)
{
	int i,j;
	ofstream ofp;
	for(i=0;i<nrec;i++)
	{
		sprintf(str,"%s.%d.txt",Fname,i+1);
		StopIfErrorReturn(OpenOutputFile(ofp,str),"OpenOutputFile");
		ofp<<scientific<<setprecision(7);
		ofp<<0.0<<'\t'<<0.0<<'\t'<<0.0<<'\n';
		ofp<<Xr[i]<<'\t'<<Yr[i]<<'\t'<<Zr[i]<<'\n';
		ofp<<'\n';
		ofp<<Leg<<'\n';
		for(j=0;j<nfreq;j++)
		{
			ofp<<sfreq[j]<<'\t'<<Zxx[i][j]<<'\n';
		}
		StopIfErrorReturn(CloseOutputFile(ofp,str),"CloseOutputFile");
	}
}

void WriteToEdsAll(ofstream &ofp,ofstream &ofm,int nrec,int nfrq,vector<double> &freq,
				   vector<double> &Xr,vector<double> &Yr,vector<double> &Zr,
				   vector<vector<double>> &Zxx,vector<vector<double>> &Zyy,
				   vector<vector<double>> &Zxy,vector<vector<double>> &Zyx,int type)
{
	int i,irec,ifrq;
	vector<vector<double>> *pZt;

	for(i=0;i<4;i++)
	{
		pZt=(i==3)? &Zyy : (i==2)? &Zyx : (i==1)? &Zxy : &Zxx;

		for(irec=0;irec<nrec;irec++)
		{
			ofm<<type<<'\n';

			ofp<<(int)Xr[irec]<<' '<<(int)Yr[irec]<<'\n';
			ofp<<Xr[irec]<<' '<<Yr[irec]<<' '<<Zr[irec]<<'\n';
			ofp<<1<<'\n';
			ofp<<1<<' '<<1<<' '<<1<<' '<<1<<'\n';
			ofp<<1<<'\n';
			ofp<<'\t'<<"nu"<<'\t'<<"Un"<<'\t'<<"Ua"<<'\t'<<"Us"<<'\n';

			for(ifrq=0;ifrq<nfrq;ifrq++)
			{
				ofp<<freq[nfrq-ifrq-1]<<' '<<0.0<<' '<<(*pZt)[irec][nfrq-ifrq-1]<<' '<<(*pZt)[irec][nfrq-ifrq-1]<<'\n';
			}
		}
	}
}

int main(int argc,char **argv)
{
	int i,j,k;
	ifstream inf;
	ofstream ofp,ofm;
	int ifrq,nfrq,irec,nrec,nrecE,jrec;
	vector<double> freq,sfreq;
	vector<vector<complex<double>>> Hxx,Hyx,Hzx,Exx,Eyx;	// Ix
	vector<vector<complex<double>>> Hxy,Hyy,Hzy,Exy,Eyy;	// Iy
	vector<double> Xr,Yr,Zr;
	int ia;

	logfile.open("LogAssembleMTZ");

	sprintf(PathToInstallDir,".");
	inf.open("PathToInstallDir");
	if(inf)
	{
		inf>>PathToInstallDir;
		inf.close();
	}
	inf.clear();

	_chdir("..\\Results");

	ia=0;
	nfrq=0;
	
	StopIfErrorReturn(OpenInputFile(inf,"frec"),"OpenInputFile");
	inf>>nfrq;
	freq.resize(nfrq);
	sfreq.resize(nfrq);
	for(ifrq=0;ifrq<nfrq;ifrq++)
	{
		inf>>freq[ifrq];	// частоты будут храниться в обратном порядке, а поля в прямом !!!
	}
	for(ifrq=0;ifrq<nfrq;ifrq++)
	{
		sfreq[ifrq]=sqrt(freq[nfrq-ifrq-1]);
	}
	StopIfErrorReturn(CloseInputFile(inf,"frec"),"CloseInputFile");

	StopIfErrorReturn(OpenInputFile(inf,"pointres"),"OpenInputFile");
	inf>>nrec;
	Xr.resize(nrec);
	Yr.resize(nrec);
	Zr.resize(nrec);
	for(irec=0;irec<nrec;irec++){inf>>Xr[irec]>>Yr[irec]>>Zr[irec];}
	StopIfErrorReturn(CloseInputFile(inf,"pointres"),"CloseInputFile");

	nrecE=nrec;

	GetMemmory(Hxx,nrec,nfrq);
	GetMemmory(Hyx,nrec,nfrq);
	GetMemmory(Hzx,nrec,nfrq);
	GetMemmory(Exx,nrecE,nfrq);
	GetMemmory(Eyx,nrecE,nfrq);

	GetMemmory(Hxy,nrec,nfrq);
	GetMemmory(Hyy,nrec,nfrq);
	GetMemmory(Hzy,nrec,nfrq);
	GetMemmory(Exy,nrecE,nfrq);
	GetMemmory(Eyy,nrecE,nfrq);

	for(ifrq=0;ifrq<nfrq;ifrq++)
	{
		ReadFrqFile(Hxx,nrec,ifrq,"hx_sx_s0f","hx_cx_s0f");
		ReadFrqFile(Hyx,nrec,ifrq,"hy_sx_s0f","hy_cx_s0f");
		ReadFrqFile(Hzx,nrec,ifrq,"hz_sx_s0f","hz_cx_s0f");
		ReadFrqFile(Exx,nrecE,ifrq,"ex_sx_s0f","ex_cx_s0f");
		ReadFrqFile(Eyx,nrecE,ifrq,"ey_sx_s0f","ey_cx_s0f");

		ReadFrqFile(Hxy,nrec,ifrq,"hx_sy_s0f","hx_cy_s0f");
		ReadFrqFile(Hyy,nrec,ifrq,"hy_sy_s0f","hy_cy_s0f");
		ReadFrqFile(Hzy,nrec,ifrq,"hz_sy_s0f","hz_cy_s0f");
		ReadFrqFile(Exy,nrecE,ifrq,"ex_sy_s0f","ex_cy_s0f");
		ReadFrqFile(Eyy,nrecE,ifrq,"ey_sy_s0f","ey_cy_s0f");
	}

	for(i=0;i<nrec;i++)
	{
		double ca,sa,xr,yr,xi,yi,vals,valc;
		double rmx,rmy,rmz;
		double rmex,rmey,rmez;
		double rax,ray,raz;
		double raex,raey,raez;

		rmx=rmy=rmz=1.0;
		rax=0.0;
		ray=90.0*PI/180.0;
		raz=0.0;
		rmex=rmey=rmez=1.0;
		raex=0.0;
		raey=90.0*PI/180.0;
		raez=0.0;

		for(j=0;j<nfrq;j++)
		{
			xr=real(Hxx[i][j]);
			yr=real(Hyx[i][j]);
			xi=imag(Hxx[i][j]);
			yi=imag(Hyx[i][j]);

			ca=cos(rax);
			sa=sin(rax);
			vals=ca*xr*rmx+sa*yr*rmy;
			valc=ca*xi*rmx+sa*yi*rmy;
			Hxx[i][j]=complex<double>(vals,valc);

			ca=cos(ray);
			sa=sin(ray);
			vals=ca*xr*rmx+sa*yr*rmy;
			valc=ca*xi*rmx+sa*yi*rmy;
			Hyx[i][j]=complex<double>(vals,valc);

			vals=real(Hzx[i][j]);
			valc=imag(Hzx[i][j]);
			vals*=rmz;
			valc*=rmz;
			Hzx[i][j]=complex<double>(vals,valc);

			xr=real(Hxy[i][j]);
			yr=real(Hyy[i][j]);
			xi=imag(Hxy[i][j]);
			yi=imag(Hyy[i][j]);

			ca=cos(rax);
			sa=sin(rax);
			vals=ca*xr*rmx+sa*yr*rmy;
			valc=ca*xi*rmx+sa*yi*rmy;
			Hxy[i][j]=complex<double>(vals,valc);

			ca=cos(ray);
			sa=sin(ray);
			vals=ca*xr*rmx+sa*yr*rmy;
			valc=ca*xi*rmx+sa*yi*rmy;
			Hyy[i][j]=complex<double>(vals,valc);

			vals=real(Hzy[i][j]);
			valc=imag(Hzy[i][j]);
			vals*=rmz;
			valc*=rmz;
			Hzy[i][j]=complex<double>(vals,valc);
		}
	}
	


	WriteAllFileReal(Hxx,nrec,nfrq,sfreq,Xr,Yr,Zr,"hx_sx_all_s0","Hx(A/m)");
	WriteAllFileReal(Hyx,nrec,nfrq,sfreq,Xr,Yr,Zr,"hy_sx_all_s0","Hy(A/m)");
	WriteAllFileReal(Hzx,nrec,nfrq,sfreq,Xr,Yr,Zr,"hz_sx_all_s0","Hz(A/m)");
	WriteAllFileReal(Hxy,nrec,nfrq,sfreq,Xr,Yr,Zr,"hx_sy_all_s0","Hx(A/m)");
	WriteAllFileReal(Hyy,nrec,nfrq,sfreq,Xr,Yr,Zr,"hy_sy_all_s0","Hy(A/m)");
	WriteAllFileReal(Hzy,nrec,nfrq,sfreq,Xr,Yr,Zr,"hz_sy_all_s0","Hz(A/m)");

	WriteAllFileImag(Hxx,nrec,nfrq,sfreq,Xr,Yr,Zr,"hx_cx_all_s0","Hx(A/m)");
	WriteAllFileImag(Hyx,nrec,nfrq,sfreq,Xr,Yr,Zr,"hy_cx_all_s0","Hy(A/m)");
	WriteAllFileImag(Hzx,nrec,nfrq,sfreq,Xr,Yr,Zr,"hz_cx_all_s0","Hz(A/m)");
	WriteAllFileImag(Hxy,nrec,nfrq,sfreq,Xr,Yr,Zr,"hx_cy_all_s0","Hx(A/m)");
	WriteAllFileImag(Hyy,nrec,nfrq,sfreq,Xr,Yr,Zr,"hy_cy_all_s0","Hy(A/m)");
	WriteAllFileImag(Hzy,nrec,nfrq,sfreq,Xr,Yr,Zr,"hz_cy_all_s0","Hz(A/m)");

	WriteAllFileReal(Exx,nrec,nfrq,sfreq,Xr,Yr,Zr,"ex_sx_all_s0","Hx(A/m)");
	WriteAllFileReal(Eyx,nrec,nfrq,sfreq,Xr,Yr,Zr,"ey_sx_all_s0","Hy(A/m)");
	WriteAllFileReal(Exy,nrec,nfrq,sfreq,Xr,Yr,Zr,"ex_sy_all_s0","Hx(A/m)");
	WriteAllFileReal(Eyy,nrec,nfrq,sfreq,Xr,Yr,Zr,"ey_sy_all_s0","Hy(A/m)");

	WriteAllFileImag(Exx,nrec,nfrq,sfreq,Xr,Yr,Zr,"ex_cx_all_s0","Hx(A/m)");
	WriteAllFileImag(Eyx,nrec,nfrq,sfreq,Xr,Yr,Zr,"ey_cx_all_s0","Hy(A/m)");
	WriteAllFileImag(Exy,nrec,nfrq,sfreq,Xr,Yr,Zr,"ex_cy_all_s0","Hx(A/m)");
	WriteAllFileImag(Eyy,nrec,nfrq,sfreq,Xr,Yr,Zr,"ey_cy_all_s0","Hy(A/m)");

	

	sprintf(str,"%s\\calcimpedance.exe",PathToInstallDir);
	if(IsFileExist(str))
	{
		run(str);
	}
	else
	{
		run("C:\\GI\\calcimpedance.exe");
	}
	
	vector<vector<double>> Zxx,Zxy,Zyx,Zyy;
	vector<vector<double>> Fxx,Fxy,Fyx,Fyy;
	vector<vector<double>> Rxx,Rxy,Ryx,Ryy;
	vector<vector<double>> Fixx,Fixy,Fiyx,Fiyy;
	vector<vector<double>> *pZt;

	GetMemmory(Zxx,nrec,nfrq);
	GetMemmory(Zxy,nrec,nfrq);
	GetMemmory(Zyx,nrec,nfrq);
	GetMemmory(Zyy,nrec,nfrq);

	GetMemmory(Fxx,nrec,nfrq);
	GetMemmory(Fxy,nrec,nfrq);
	GetMemmory(Fyx,nrec,nfrq);
	GetMemmory(Fyy,nrec,nfrq);

	GetMemmory(Rxx,nrec,nfrq);
	GetMemmory(Rxy,nrec,nfrq);
	GetMemmory(Ryx,nrec,nfrq);
	GetMemmory(Ryy,nrec,nfrq);

	GetMemmory(Fixx,nrec,nfrq);
	GetMemmory(Fixy,nrec,nfrq);
	GetMemmory(Fiyx,nrec,nfrq);
	GetMemmory(Fiyy,nrec,nfrq);



	ReadZTensor(nrec,nfrq,Zxx,Zxy,Zyx,Zyy,"Zxx","Zxy","Zyx","Zyy");
	ReadZTensor(nrec,nfrq,Fxx,Fxy,Fyx,Fyy,"Fxx","Fxy","Fyx","Fyy");
	ReadZTensor(nrec,nfrq,Rxx,Rxy,Ryx,Ryy,"Rxx","Rxy","Ryx","Ryy");
	ReadZTensor(nrec,nfrq,Fixx,Fixy,Fiyx,Fiyy,"Fixx","Fixy","Fiyx","Fiyy");

	for(k=0;k<16;k++)
	{
		if((k>=4 && k<=7) || (k>=12 && k<=15))
		{
			switch(k)
			{
				case 4 : {pZt=&Fxx;break;}
				case 5 : {pZt=&Fxy;break;}
				case 6 : {pZt=&Fyx;break;}
				case 7 : {pZt=&Fyy;break;}
				case 12 : {pZt=&Fixx;break;}
				case 13 : {pZt=&Fixy;break;}
				case 14 : {pZt=&Fiyx;break;}
				case 15 : {pZt=&Fiyy;break;}
			}
			
			for(irec=0;irec<nrec;irec++)
			{
				for(ifrq=0;ifrq<nfrq;ifrq++)
				{
					if((*pZt)[irec][ifrq]>180.0)
					{
						(*pZt)[irec][ifrq]-=180.0;
					}
					if((*pZt)[irec][ifrq]<-180.0)
					{
						(*pZt)[irec][ifrq]+=180.0;
					}
				}
			}
		}
	}

	ofm.open("ft0_0");

	sprintf(str,"edsall0_0_imp");
	ofp.open(str);

	WriteToEdsAll(ofp,ofm,nrec,nfrq,freq,Xr,Yr,Zr,Zxx,Zyy,Zxy,Zyx,1);
	WriteToEdsAll(ofp,ofm,nrec,nfrq,freq,Xr,Yr,Zr,Fxx,Fyy,Fxy,Fyx,2);
	WriteToEdsAll(ofp,ofm,nrec,nfrq,freq,Xr,Yr,Zr,Rxx,Ryy,Rxy,Ryx,3);
	WriteToEdsAll(ofp,ofm,nrec,nfrq,freq,Xr,Yr,Zr,Fixx,Fiyy,Fixy,Fiyx,4);

	ofp.close();
	ofp.clear();

	ofm.close();
	ofm.clear();

	logfile.close();
	logfile.clear();

	return 0;
}
