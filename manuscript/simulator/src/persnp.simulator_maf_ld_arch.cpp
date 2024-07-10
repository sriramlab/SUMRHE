#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector> 
//#include <random>

#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/QR>
#include "time.h"

#include "genotype.h"
#include "mailman.h"
#include "arguments.h"
//#include "helper.h"
#include "storage.h"


#include "/usr/include/boost/numeric/ublas/matrix.hpp"
#include "/usr/include/boost/numeric/ublas/io.hpp"
#include "/usr/include/boost/random.hpp"
#include "/usr/include/boost/numeric/ublas/lu.hpp"
#include "/usr/include/boost/numeric/ublas/vector.hpp"
#include "/usr/include/boost/numeric/ublas/vector_proxy.hpp"
#include "/usr/include/boost/numeric/ublas/triangular.hpp"
#include <boost/lexical_cast.hpp>

#if SSE_SUPPORT==1
	#define fastmultiply fastmultiply_sse
	#define fastmultiply_pre fastmultiply_pre_sse
#else
	#define fastmultiply fastmultiply_normal
	#define fastmultiply_pre fastmultiply_pre_normal
#endif

#ifdef _WIN32
#define rdtsc __rdtsc
#else

unsigned long long rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;
}
#endif

using namespace Eigen;
using namespace std;

// Storing in RowMajor Form
typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdr;
//typedef Matrix<int, Dynamic, Dynamic, RowMajor> MatrixXdrInt;

class data {

 public:
     MatrixXdr gen;
     int index;
     MatrixXdr maf;
};


//Intermediate Variables
int blocksize;
int hsegsize;
double *partialsums;
double *sum_op;		
double *yint_e;
double *yint_m;
double **y_e;
double **y_m;


struct timespec t0;

clock_t total_begin = clock();
MatrixXdr pheno;
MatrixXdr mask;
MatrixXdr covariate;  
MatrixXdr Q;
MatrixXdr v1; //W^ty
MatrixXdr v2;            //QW^ty
MatrixXdr v3;    //WQW^ty
MatrixXdr new_pheno;



genotype g;
genotype g1;
genotype g2;
MatrixXdr geno_matrix; //(p,n)
genotype* Geno;
int MAX_ITER;
int k,p,n;
int k_orig;

MatrixXdr c; //(p,k)
MatrixXdr x; //(k,n)
MatrixXdr v; //(p,k)
MatrixXdr means; //(p,1)
MatrixXdr stds; //(p,1)
MatrixXdr sum2;
MatrixXdr sum;  

Eigen::RowVectorXd means_na;
Eigen::RowVectorXd stds_na;
////////
//related to phenotype	
double y_sum; 
double y_mean;

options command_line_opts;

bool debug = false;
bool check_accuracy = false;
bool var_normalize=false;
int accelerated_em=0;
double convergence_limit;
bool memory_efficient = false;
bool missing=false;
bool fast_mode = true;
bool text_version = false;
bool use_cov=false; 


/////
bool mean_impute=true;
std::vector<int> miss_index;
/////
bool dominance=false;
bool only_dominance=false;
int nongen_Nbin;
int total_Nbin;

///
//// jackknife index wich are computed based on annotation file
MatrixXdr dic_index;
MatrixXdr jack_bin_size;
vector<int> len;
vector<int> Annot;
int Njack;
int Nbin;
int Nz=10;
///////

//define random vector z's
MatrixXdr  all_zb;
MatrixXdr  all_Uzb;
MatrixXdr res;
MatrixXdr XXz;
MatrixXdr Xy;
MatrixXdr yXXy;

///// for simulation
MatrixXdr maf_ld;
MatrixXdr A_effect;
MatrixXdr A_effect_multi;
vector<string> fid;
MatrixXdr simul_pheno;
MatrixXdr simul_par;
int num_simul;
/////vars for dominace

MatrixXdr global_MAF;

///
//Matrix<int, Dynamic, Dynamic, RowMajor> gen;
MatrixXdr gen;
bool read_header;
//read variables
unsigned char mask2;
int wordsize;
unsigned int unitsperword;
int unitsize;
int nrow, ncol;
unsigned char *gtype;
int Nsnp;
int Nindv;
bool **bin_annot;
int step_size;
int step_size_rem;
std::vector<std::vector<bool> > annot_bool;
std::vector<std::vector<int> > jack_bin;
vector <data> allgen;
vector <genotype> allgen_mail;
int global_snp_index;
bool use_mailman=true;

///variable for non geno vc

vector <data> nongen_all;
MatrixXdr nongen_vc;

///var for analytic se

MatrixXdr wt;
MatrixXdr vt;

///reading single col annot
vector <int> SNP_annot;
bool use_1col_annot=false;


///Variables for reg out cov on both side of LM
bool both_side_cov=false;
MatrixXdr UXXz;
MatrixXdr XXUz;
MatrixXdr Xz;
MatrixXdr trVK;

std::istream& newline(std::istream& in)
{
    if ((in >> std::ws).peek() != std::char_traits<char>::to_int_type('\n')) {
        in.setstate(std::ios_base::failbit);
    }
    return in.ignore();
}


int read_cov(bool std,int Nind, std::string filename, std::string covname){
	ifstream ifs(filename.c_str(), ios::in); 
	std::string line; 
	std::istringstream in; 
	int covIndex = 0; 
	std::getline(ifs,line); 
	in.str(line); 
	string b;
	vector<vector<int> > missing; 
	int covNum=0;  
	while(in>>b)
	{
		if(b!="FID" && b !="IID"){
		missing.push_back(vector<int>()); //push an empty row  
		if(b==covname && covname!="")
			covIndex=covNum; 
		covNum++; 
		}
	}
	vector<double> cov_sum(covNum, 0); 
	if(covname=="")
	{
		covariate.resize(Nind, covNum); 
		cout<< "Read in "<<covNum << " Covariates.. "<<endl;
	}
	else 
	{
		covariate.resize(Nind, 1); 
		cout<< "Read in covariate "<<covname<<endl;  
	}

	
	int j=0; 
	while(std::getline(ifs, line)){
		in.clear(); 
		in.str(line);
		string temp;
		in>>temp; in>>temp; //FID IID 
		for(int k=0; k<covNum; k++){
			
			in>>temp;
			if(temp=="NA")
			{
				missing[k].push_back(j);
				continue; 
			} 
			double cur = atof(temp.c_str()); 
			if(cur==-9)
			{
				missing[k].push_back(j); 
				continue; 
			}
			if(covname=="")
			{
				cov_sum[k]= cov_sum[k]+ cur; 
				covariate(j,k) = cur; 
			}
			else
				if(k==covIndex)
				{
					covariate(j, 0) = cur;
					cov_sum[k] = cov_sum[k]+cur; 
				}
		}
		//if(j<10) 
		//	cout<<covariate.block(j,0,1, covNum)<<endl; 
		j++;
	}
	//compute cov mean and impute 
	for (int a=0; a<covNum ; a++)
	{
		int missing_num = missing[a].size(); 
		cov_sum[a] = cov_sum[a] / (Nind - missing_num);

		for(int b=0; b<missing_num; b++)
		{
                        int index = missing[a][b];
                        if(covname=="")
                                covariate(index, a) = cov_sum[a];
                        else if (a==covIndex)
                                covariate(index, 0) = cov_sum[a];
                } 
	}
	if(std)
	{
		MatrixXdr cov_std;
		cov_std.resize(1,covNum);  
		MatrixXdr sum = covariate.colwise().sum();
		MatrixXdr sum2 = (covariate.cwiseProduct(covariate)).colwise().sum();
		MatrixXdr temp;
//		temp.resize(Nind, 1); 
//		for(int i=0; i<Nind; i++)
//			temp(i,0)=1;  
		for(int b=0; b<covNum; b++)
		{
			cov_std(0,b) = sum2(0,b) + Nind*cov_sum[b]*cov_sum[b]- 2*cov_sum[b]*sum(0,b);
			cov_std(0,b) =sqrt((Nind- 1)/cov_std(0,b)) ;
			double scalar=cov_std(0,b); 
			for(int j=0; j<Nind; j++)
			{
				covariate(j,b) = covariate(j,b)-cov_sum[b];  
				covariate(j,b) =covariate(j,b)*scalar;
			} 
			//covariate.col(b) = covariate.col(b) -temp*cov_sum[b];
			
		}
	}	
	return covNum; 
}







void multiply_y_pre_fast(MatrixXdr &op, int Ncol_op ,MatrixXdr &res,bool subtract_means){

        for(int k_iter=0;k_iter<Ncol_op;k_iter++){
                sum_op[k_iter]=op.col(k_iter).sum();
        }

                        //cout << "Nops = " << Ncol_op << "\t" <<g.Nsegments_hori << endl;
        #if DEBUG==1
                if(debug){
                        print_time ();
                        cout <<"Starting mailman on premultiply"<<endl;
                        cout << "Nops = " << Ncol_op << "\t" <<g.Nsegments_hori << endl;
                        cout << "Segment size = " << g.segment_size_hori << endl;
                        cout << "Matrix size = " <<g.segment_size_hori<<"\t" <<g.Nindv << endl;
                        cout << "op = " <<  op.rows () << "\t" << op.cols () << endl;
                }
        #endif


        //TODO: Memory Effecient SSE FastMultipy

        for(int seg_iter=0;seg_iter<g.Nsegments_hori-1;seg_iter++){
                mailman::fastmultiply(g.segment_size_hori,g.Nindv,Ncol_op,g.p[seg_iter],op,yint_m,partialsums,y_m);
                int p_base = seg_iter*g.segment_size_hori;
                for(int p_iter=p_base; (p_iter<p_base+g.segment_size_hori) && (p_iter<g.Nsnp) ; p_iter++ ){
                        for(int k_iter=0;k_iter<Ncol_op;k_iter++)
                                res(p_iter,k_iter) = y_m[p_iter-p_base][k_iter];
                }
        }

        int last_seg_size = (g.Nsnp%g.segment_size_hori !=0 ) ? g.Nsnp%g.segment_size_hori : g.segment_size_hori;
        mailman::fastmultiply(last_seg_size,g.Nindv,Ncol_op,g.p[g.Nsegments_hori-1],op,yint_m,partialsums,y_m);
        int p_base = (g.Nsegments_hori-1)*g.segment_size_hori;
        for(int p_iter=p_base; (p_iter<p_base+g.segment_size_hori) && (p_iter<g.Nsnp) ; p_iter++){
                for(int k_iter=0;k_iter<Ncol_op;k_iter++)
                        res(p_iter,k_iter) = y_m[p_iter-p_base][k_iter];
        }

        #if DEBUG==1
                if(debug){
                        print_time ();
                        cout <<"Ending mailman on premultiply"<<endl;
                }
        #endif


        if(!subtract_means)
                return;

        for(int p_iter=0;p_iter<p;p_iter++){
                for(int k_iter=0;k_iter<Ncol_op;k_iter++){
                        res(p_iter,k_iter) = res(p_iter,k_iter) - (g.get_col_mean(p_iter)*sum_op[k_iter]);
                        if(var_normalize)
                                res(p_iter,k_iter) = res(p_iter,k_iter)/(g.get_col_std(p_iter));
                }
        }

}



void multiply_y_post_fast(MatrixXdr &op_orig, int Nrows_op, MatrixXdr &res,bool subtract_means){

        MatrixXdr op;
        op = op_orig.transpose();

        if(var_normalize && subtract_means){
                for(int p_iter=0;p_iter<p;p_iter++){
                        for(int k_iter=0;k_iter<Nrows_op;k_iter++)
                                op(p_iter,k_iter) = op(p_iter,k_iter) / (g.get_col_std(p_iter));
                }
        }

        #if DEBUG==1
                if(debug){
                        print_time ();
                        cout <<"Starting mailman on postmultiply"<<endl;
                }
        #endif

        int Ncol_op = Nrows_op;

        //cout << "ncol_op = " << Ncol_op << endl;

        int seg_iter;
        for(seg_iter=0;seg_iter<g.Nsegments_hori-1;seg_iter++){
mailman::fastmultiply_pre(g.segment_size_hori,g.Nindv,Ncol_op, seg_iter * g.segment_size_hori, g.p[seg_iter],op,yint_e,partialsums,y_e);
        }
        int last_seg_size = (g.Nsnp%g.segment_size_hori !=0 ) ? g.Nsnp%g.segment_size_hori : g.segment_size_hori;
        mailman::fastmultiply_pre(last_seg_size,g.Nindv,Ncol_op, seg_iter * g.segment_size_hori, g.p[seg_iter],op,yint_e,partialsums,y_e);

        for(int n_iter=0; n_iter<n; n_iter++)  {
                for(int k_iter=0;k_iter<Ncol_op;k_iter++) {
                        res(k_iter,n_iter) = y_e[n_iter][k_iter];
                        y_e[n_iter][k_iter] = 0;
                }
        }

        #if DEBUG==1
                if(debug){
                        print_time ();
                        cout <<"Ending mailman on postmultiply"<<endl;
                }
        #endif


        if(!subtract_means)
                return;

        double *sums_elements = new double[Ncol_op];
        memset (sums_elements, 0, Nrows_op * sizeof(int));

        for(int k_iter=0;k_iter<Ncol_op;k_iter++){
                double sum_to_calc=0.0;
                for(int p_iter=0;p_iter<p;p_iter++)
                        sum_to_calc += g.get_col_mean(p_iter)*op(p_iter,k_iter);
                sums_elements[k_iter] = sum_to_calc;
        }
        for(int k_iter=0;k_iter<Ncol_op;k_iter++){
                for(int n_iter=0;n_iter<n;n_iter++)
                        res(k_iter,n_iter) = res(k_iter,n_iter) - sums_elements[k_iter];
        }


}


void initial_var()
{
    /*if(key==1)
        g=g1;
    if(key==2)
        g=g2;*/
   // g=Geno[key];
        

	p = g.Nsnp;
        n = g.Nindv;


        c.resize(p,k);
        x.resize(k,n);
        v.resize(p,k);
        //means.resize(p,1);
        //stds.resize(p,1);
        sum2.resize(p,1);
        sum.resize(p,1);


        if(!fast_mode && !memory_efficient){
                geno_matrix.resize(p,n);
                g.generate_eigen_geno(geno_matrix,var_normalize);
        }

        //TODO: Initialization of c with gaussian distribution
        c = MatrixXdr::Random(p,k);


        // Initial intermediate data structures
        blocksize = k;
         hsegsize = g.segment_size_hori;        // = log_3(n)
        int hsize = pow(3,hsegsize);
        int vsegsize = g.segment_size_ver;              // = log_3(p)
        int vsize = pow(3,vsegsize);

        partialsums = new double [blocksize];
        sum_op = new double[blocksize];
        yint_e = new double [hsize*blocksize];
        yint_m = new double [hsize*blocksize];
        memset (yint_m, 0, hsize*blocksize * sizeof(double));
        memset (yint_e, 0, hsize*blocksize * sizeof(double));

        y_e  = new double*[g.Nindv];
        for (int i = 0 ; i < g.Nindv ; i++) {
                y_e[i] = new double[blocksize];
                memset (y_e[i], 0, blocksize * sizeof(double));
        }

        y_m = new double*[hsegsize];
        for (int i = 0 ; i < hsegsize ; i++)
                y_m[i] = new double[blocksize];
      /*  for(int i=0;i<p;i++){
                means(i,0) = g.get_col_mean(i);
                stds(i,0) =1/g.get_col_std(i);
                //sum2(i,0) =g.get_col_sum2(i); 
                sum(i,0)= g.get_col_sum(i);
        }

*/


}
 



/*void read_cov(int Nind, std::string filename, std::string covname){
	ifstream ifs(filename.c_str(), ios::in); 
	std::string line; 
	std::istringstream in; 
	int covIndex = 0; 
	std::getline(ifs,line); 
	in.str(line); 
	string b;
	vector<vector<int> > missing; 
	int covNum=0;  
	while(in>>b)
	{
		missing.push_back(vector<int>()); //push an empty row  
		if(b==covname && covname!="")
			covIndex=covNum; 
		covNum++; 
	}
	vector<double> cov_sum(covNum, 0); 
	if(covname=="")
	{
		covariate.resize(Nind, covNum); 
		cout<< "Read in "<<covNum << " Covariates.. "<<endl;
	}
	else 
	{
		covariate.resize(Nind, 1); 
		cout<< "Read in covariate "<<covname<<endl;  
	}

	
	int j=0; 
	while(std::getline(ifs, line)){
		in.clear(); 
		in.str(line);
		string temp; 
		for(int k=0; k<covNum; k++){
			in>>temp;
			if(temp=="NA")
			{
				missing[k].push_back(j);
				continue;  
			} 
			int cur = atof(temp.c_str()); 
			if(cur==-9)
			{
				missing[k].push_back(j); 
				continue; 
			}
			if(covname=="")
			{
				cov_sum[k]= cov_sum[k]+ cur; 
				covariate(j,k) = cur; 
			}
			else
				if(k==covIndex)
				{
					covariate(j, 0) = cur;
					cov_sum[k] = cov_sum[k]+cur; 
				}
		} 
		j++;
	}
	//compute cov mean and impute 
	for (int a=0; a<covNum ; a++)
	{
		int missing_num = missing[a].size(); 
		cov_sum[a] = cov_sum[a] / (covNum - missing_num);

		for(int b=0; b<missing_num; b++)
		{
                        int index = missing[a][b];
                        if(covname=="")
                                covariate(index, a) = cov_sum[a];
                        else if (a==covIndex)
                                covariate(index, 0) = cov_sum[a];
                } 
	}
}*/
void read_pheno2(int Nind, std::string filename){
//	pheno.resize(Nind,1); 
	ifstream ifs(filename.c_str(), ios::in); 
	
	std::string line;
	std::istringstream in;  
	int phenocount=0; 
//read header
	std::getline(ifs,line); 
	in.str(line); 
	string b; 
	while(in>>b)
	{
		if(b!="FID" && b !="IID")
			phenocount++; 
	}
	pheno.resize(Nind, phenocount);
	mask.resize(Nind, phenocount);
	int i=0;  
	while(std::getline(ifs, line)){
		in.clear(); 
		in.str(line); 
		string temp;
		//fid,iid
		//todo: fid iid mapping; 
		//todo: handle missing phenotype
		in>>temp; in>>temp; 
		for(int j=0; j<phenocount;j++) {
			in>>temp;
			double cur = atof(temp.c_str());
			if(temp=="NA" || cur==-9){
			pheno(i,j)=0;
			mask(i,j)=0;
			}
			else{
			pheno(i,j)=atof(temp.c_str());
			mask(i,j)=1;

			}

    
		}
		i++;
	}
	//cout<<pheno; 
}

double compute_yXXy(int num_snp){


        MatrixXdr res(num_snp, 1);

	
	if(use_mailman==true)
		multiply_y_pre_fast(pheno,1,res,false);
	else
		 res=gen*pheno;
	

	res = res.cwiseProduct(stds);
        MatrixXdr resid(num_snp, 1);
        resid = means.cwiseProduct(stds);
        resid = resid *y_sum;
        MatrixXdr Xy(num_snp,1);
        Xy = res-resid;
    
        double yXXy = (Xy.array()* Xy.array()).sum();
        

	//cout<<yXXy<<endl;	

        return yXXy;

}

double compute_yVXXVy(int num_snp){
        MatrixXdr new_pheno_sum = new_pheno.colwise().sum();
        MatrixXdr res(num_snp, 1);
        
	



        if(use_mailman==true)
                multiply_y_pre_fast(new_pheno,1,res,false);
        else
                 res=gen*new_pheno;



        res = res.cwiseProduct(stds);
        MatrixXdr resid(num_snp, 1);
        resid = means.cwiseProduct(stds);
        resid = resid *new_pheno_sum;
        MatrixXdr Xy(num_snp,1);
        Xy = res-resid;
        double ytVXXVy = (Xy.array()* Xy.array()).sum();
        return ytVXXVy;

}





	
MatrixXdr  compute_XXz (int num_snp){
	//mask
	/*for (int i=0;i<Nz;i++)
	   for(int j=0;j<Nindv;j++)
		 all_zb(j,i)=all_zb(j,i)*mask(j,0);
*/
         res.resize(num_snp, Nz);

        
	if(use_mailman==true)
		multiply_y_pre_fast(all_zb,Nz,res, false);
	else	
        	res=gen*all_zb;
   

        MatrixXdr zb_sum = all_zb.colwise().sum();
        

	for(int j=0; j<num_snp; j++)
            for(int k=0; k<Nz;k++)
                res(j,k) = res(j,k)*stds(j,0);
            
        MatrixXdr resid(num_snp, Nz);
        MatrixXdr inter = means.cwiseProduct(stds);
        resid = inter * zb_sum;
        MatrixXdr inter_zb = res - resid;
       

	for(int k=0; k<Nz; k++)
            for(int j=0; j<num_snp;j++)
                inter_zb(j,k) =inter_zb(j,k) *stds(j,0);
       MatrixXdr new_zb = inter_zb.transpose();
       MatrixXdr new_res(Nz, Nindv);
       
	
	 if(use_mailman==true)
	    multiply_y_post_fast(new_zb, Nz, new_res, false);
	 else
	    new_res=new_zb*gen;	

       MatrixXdr new_resid(Nz, num_snp);
       MatrixXdr zb_scale_sum = new_zb * means;
       new_resid = zb_scale_sum * MatrixXdr::Constant(1,Nindv, 1);


                      /// new zb 
       MatrixXdr temp=new_res - new_resid;

	for (int i=0;i<Nz;i++)
           for(int j=0;j<Nindv;j++)
                 temp(i,j)=temp(i,j)*mask(j,0);


	return temp.transpose();
       

}






MatrixXdr  compute_XXUz (int num_snp){
        //mask
/*        for (int i=0;i<Nz;i++)
           for(int j=0;j<Nindv;j++)
                 all_Uzb(j,i)=all_Uzb(j,i)*mask(j,0);
*/
         res.resize(num_snp, Nz);


        if(use_mailman==true)
                multiply_y_pre_fast(all_Uzb,Nz,res, false);
        else
                res=gen*all_Uzb;
  

        MatrixXdr zb_sum = all_Uzb.colwise().sum();


        for(int j=0; j<num_snp; j++)
            for(int k=0; k<Nz;k++)
                res(j,k) = res(j,k)*stds(j,0);

        MatrixXdr resid(num_snp, Nz);
        MatrixXdr inter = means.cwiseProduct(stds);
        resid = inter * zb_sum;
        MatrixXdr inter_zb = res - resid;


        for(int k=0; k<Nz; k++)
            for(int j=0; j<num_snp;j++)
                inter_zb(j,k) =inter_zb(j,k) *stds(j,0);
       MatrixXdr new_zb = inter_zb.transpose();
       MatrixXdr new_res(Nz, Nindv);


         if(use_mailman==true)
            multiply_y_post_fast(new_zb, Nz, new_res, false);
         else
            new_res=new_zb*gen;

       MatrixXdr new_resid(Nz, num_snp);
       MatrixXdr zb_scale_sum = new_zb * means;
       new_resid = zb_scale_sum * MatrixXdr::Constant(1,Nindv, 1);


                      /// new zb 
       MatrixXdr temp=new_res - new_resid;

        for (int i=0;i<Nz;i++)
           for(int j=0;j<Nindv;j++)
                 temp(i,j)=temp(i,j)*mask(j,0);


        return temp.transpose();


}



MatrixXdr  compute_Xz (int num_snp){

         MatrixXdr new_zb= MatrixXdr::Random(Nz,num_snp);
         new_zb = new_zb * sqrt(3);

	 MatrixXdr new_res(Nz, Nindv);         

        if(use_mailman==true)
                multiply_y_post_fast(new_zb,Nz,new_res, false);
        else
                new_res=new_zb*gen;


        MatrixXdr new_resid(Nz, num_snp);
       MatrixXdr zb_scale_sum = new_zb * means;
       new_resid = zb_scale_sum * MatrixXdr::Constant(1,Nindv, 1);


                      /// new zb 
       MatrixXdr temp=new_res - new_resid;

        for (int i=0;i<Nz;i++)
           for(int j=0;j<Nindv;j++)
                 temp(i,j)=temp(i,j)*mask(j,0);




	return temp.transpose();
}








void read_annot (string filename){
         
//	int step_size=Nsnp/Njack;
  //      int step_size_rem=Nsnp%Njack;
        vector<bool> snp_annot;
	//jack_bin.resize(Njack, vector<int>(Nbin,0));	
	
//	cout<<step_size<<endl;

	ifstream inp(filename.c_str());
        if (!inp.is_open()){
                cerr << "Error reading file "<< filename <<endl;
                exit(1);
        }
        string line;
        int j = 0 ;
        int linenum = 0 ;
        int num_parti;
        stringstream check1(line);
        string intermediate;
        vector <string> tokens;
        while(std::getline (inp, line)){
                //linenum ++;
                char c = line[0];
                if (c=='#')
                        continue;
                istringstream ss (line);
                if (line.empty())
                        continue;
                j++;
                //cout<<line<<endl;

                stringstream check1(line);
                string intermediate;
                vector <string> tokens;
                // Tokenizing w.r.t. space ' ' 
                while(getline(check1, intermediate, ' '))
                 {
                      tokens.push_back(intermediate);
                 }
                 if(linenum==0){
                 num_parti=tokens.size();
		 Nbin=num_parti;
                  snp_annot.resize(Nbin,0);	 
          	  jack_bin.resize(Njack, vector<int>(Nbin,0));
	         len.resize(num_parti,0);
                }
                int index_annot=0;
                for(int i = 0; i < tokens.size(); i++){
		      snp_annot[i]=0;
		      if (tokens[i]=="1"){
                            len[i]++;
			    snp_annot[i]=1;
 		      }
                }
		annot_bool.push_back(snp_annot);
        linenum++;
       }

	  if(Nsnp!=linenum){
          cout<<"Number of the rows in bim file and annotation file does not match"<<endl;
        }

	Nsnp=linenum;
	//cout<<"Total number of SNPs : "<<Nsnp<<endl;
	int selected_snps=0;
	for (int i=0;i<num_parti;i++){
                cout<<len[i]<<" SNPs in "<<i<<"-th bin"<<endl;
		selected_snps+=len[i];
        }
	
	cout<<" Number of selected SNPs w.r.t  annot file : " <<selected_snps<<endl;


	 step_size=Nsnp/Njack;
         step_size_rem=Nsnp%Njack;
	cout<<"Number of SNPs per block : "<<step_size<<endl;
      //  cout<<"stepsize : "<<step_size_rem<<endl;
	jack_bin.resize(Njack, vector<int>(Nbin,0));
	int temp;
	for (int i=0;i<Nsnp;i++)
	   for(int j=0;j<Nbin;j++)
		 if (annot_bool[i][j]==1){
			temp=i/step_size;
			if (temp==Njack)
				temp--;
			//cout<<i<<"xxx"<<j<<"xxx"<<temp<<endl;
			jack_bin[temp][j]++;	
		 }
/*	
cout<<"jackbin"<<endl;
	for (int i=0;i<Njack;i++){
	   for(int j=0;j<Nbin;j++)
                cout<<jack_bin[i][j]<<" ";
	  cout<<endl;
        }*/
/*
for (int i=0;i<linenum;i++){
  for(int j=0;j<Nbin;j++)
	cout<<annot_bool[i][j]<<" ";
  cout<<endl;
}
*/

}


//vector <int> SNP_annot;

void read_annot_1col (string filename){
	
	
        ifstream ifs(filename.c_str(), ios::in);

        std::string line;
        std::istringstream in;
        
	
	len.resize(Nbin,0);
	 step_size=Nsnp/Njack;
         step_size_rem=Nsnp%Njack;
        cout<<"Number of SNPs per block : "<<step_size<<endl;
	jack_bin.resize(Njack, vector<int>(Nbin,0));
	int i=0;
        while(std::getline(ifs, line)){
                in.clear();
                in.str(line);
                string temp;
                
		in>>temp;        
	        int  cur = atoi(temp.c_str());
		SNP_annot.push_back(cur);
		len[cur-1]++;
	
		int jack_val=i/step_size;
		 if (jack_val==Njack)
                    jack_val--;
		jack_bin[jack_val][SNP_annot[i]-1]++;

               
		 i++;
        }

	if(Nsnp!=i){
	  cout<<"Number of the rows in bim file and annotation file does not match"<<endl;
	}


        cout<<"Total number of SNPs : "<<Nsnp<<endl;
        for (int i=0;i<Nbin;i++){
                cout<<len[i]<<" SNPs in "<<i<<"-th bin"<<endl;
        }

        /*for (int i=0;i<Njack;i++){
           for(int j=0;j<Nbin;j++)
                cout<<jack_bin[i][j]<<" ";
          cout<<endl;
        }*/


}

void read_bim (string filename){
        ifstream inp(filename.c_str());
        if (!inp.is_open()){
                cerr << "Error reading file "<< filename <<endl;
                exit(1);
        }
        string line;
        int j = 0 ;
        int linenum = 0 ;
        while(std::getline (inp, line)){
                linenum ++;
                char c = line[0];
                if (c=='#')
                        continue;
                istringstream ss (line);
                if (line.empty())
                        continue;
                j++;
        }
        Nsnp = j;
        inp.close();
	cout<<"#SNP in bim file "<<Nsnp<<endl;
}





void count_pheno(std::string filename){
        ifstream ifs(filename.c_str(), ios::in);

        std::string line;
        int i=0;
        while(std::getline(ifs, line)){
                i++;
        }
        Nindv=i-1;
}




int  count_fam(std::string filename){
        ifstream ifs(filename.c_str(), ios::in);

        std::string line;
        int i=0;
        while(std::getline(ifs, line)){
                i++;
        }
        return i;
}

void  create_simul_pheno(std::string filename){
	
	//vector<string> fid;
	 ifstream ifs(filename.c_str(), ios::in);

        std::string line;
        std::istringstream in;
        string temp;


        double val;
        int i=0;
        while(std::getline(ifs, line)){
                  in.clear();
                  in.str(line);

                  in >> temp;
                  fid.push_back(temp.c_str());
			//val = atof(temp.c_str());
		//	cout<<val<<endl;
        //          cout<<fid[i]<<endl;
		  i++;
        }
	Nindv=i;
        
}






//// functions related to reading without mailman

template<typename T>
static std::istream & binary_read(std::istream& stream, T& value){
        return stream.read(reinterpret_cast<char*>(&value), sizeof(T));
}

void set_metadata() {
        wordsize = sizeof(char) * 8;
        unitsize = 2;
        unitsperword = wordsize/unitsize;
        mask2 = 0;
        for (int i = 0 ; i < unitsize; i++)
                mask2 = mask2 |(0x1<<i);
    nrow = Nsnp;
    ncol = ceil(1.0*Nindv/unitsperword);
}


int simulate2_geno_from_random(float p_j){
        float rval = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float dist_pj[3] = { (1-p_j)*(1-p_j), 2*p_j*(1-p_j), p_j*p_j };
        if(rval < dist_pj[0] )
                return 0;
        else if( rval >= dist_pj[0] && rval < (dist_pj[0]+dist_pj[1]))
                return 1;
        else
                return 2;
}









float get_observed_pj(const unsigned char* line){
        int y[4];
        int observed_sum=0;
        int observed_ct=0;

        for (int k = 0 ;k < ncol ; k++) {
                unsigned char c = line [k];
                y[0] = (c)&mask2;
                y[1] = (c>>2)&mask2;
                y[2] = (c>>4)&mask2;
                y[3] = (c>>6)&mask2;
                int j0 = k * unitsperword;
                int lmax = 4;
                if (k == ncol - 1)  {
                        lmax = Nindv%4;
                        lmax = (lmax==0)?4:lmax;
                }
                for ( int l = 0 ; l < lmax; l++){
                        int j = j0 + l ;
                        // Extract  PLINK coded genotype and convert into 0/1/2
                        // // PLINK coding: 
                        // // 00->0
                        // // 01->missing
                        // // 10->1
                        // // 11->2
                        int val = y[l];
                        val-- ;
                        if(val != 0){
                                val =  (val < 0 ) ? 0 :val ;
                                observed_sum += val;
                                observed_ct ++;
                        }
                }
        }
        return observed_sum*0.5/observed_ct;

}



void read_bed (std::istream& ifs,bool allow_missing,int num_snp)  {
         //ifstream ifs (filename.c_str(), ios::in|ios::binary);
        char magic[3];
        set_metadata ();

    gtype =  new unsigned char[ncol];

     if(read_header)  
      binary_read(ifs,magic);

        int sum=0;

        // Note that the coding of 0 and 2 can get flipped relative to plink because plink uses allele frequency (minor)
        // allele to code a SNP as 0 or 1.
        // This flipping does not matter for results.
        int y[4];
        
for(int i=0;i<num_snp;i++){
		global_snp_index++;
                ifs.read (reinterpret_cast<char*>(gtype), ncol*sizeof(unsigned char));
                float p_j = get_observed_pj(gtype);
        for (int k = 0 ;k < ncol ; k++) {
                unsigned char c = gtype [k];
                        // Extract PLINK genotypes
                y[0] = (c)&mask2;
                y[1] = (c>>2)&mask2;
                y[2] = (c>>4)&mask2;
                y[3] = (c>>6)&mask2;
                        int j0 = k * unitsperword;
                        // Handle number of individuals not being a multiple of 4
                        int lmax = 4;
                        if (k == ncol - 1)  {
                                lmax = Nindv%4;
                                lmax = (lmax==0)?4:lmax;
                        }
                        for ( int l = 0 ; l < lmax; l++){
                                int j = j0 + l ;
                                // Extract  PLINK coded genotype and convert into 0/1/2
                                // PLINK coding: 
                                // 00->0
                                // 01->missing
                                // 10->1
                                // 11->2
                                int val = y[l];
                                if(val==1 && !allow_missing){
                                        val = simulate2_geno_from_random(p_j);
                                        val++;
                                        val = (val==1) ? 0 : val;
  //                                 val=0;
				 }
                                val-- ;
                                val =  (val < 0 ) ? 0 :val ;
				sum += val;
			   
			    for(int bin_index=0;bin_index<Nbin;bin_index++){
				if(annot_bool[global_snp_index][bin_index]==1){
                        	    
				      int snp_index;
				     //int snp_index=allgen[bin_index].index;
				     //allgen[bin_index].gen(snp_index,j)=val;
						
				     if(use_mailman==true){
				 	 snp_index=allgen_mail[bin_index].index;
					 int horiz_seg_no = snp_index/allgen_mail[bin_index].segment_size_hori; 
					 allgen_mail[bin_index].p[horiz_seg_no][j] = 3 *allgen_mail[bin_index].p[horiz_seg_no][j]  + val;	 
				     // computing sum for every snp to compute mean
				         allgen_mail[bin_index].columnsum[snp_index]+=val;

				      }
				     else{
					 snp_index=allgen[bin_index].index;
                                         allgen[bin_index].gen(snp_index,j)=val;
				     }
				
				}
 			     
			    }

                    }
        }

    for(int bin_index=0;bin_index<Nbin;bin_index++)
       if(annot_bool[global_snp_index][bin_index]==1){
       		//cout<<"global"<<global_snp_index<<endl;
		        //cout<<allgen[bin_index].index<<endl; 
		//       allgen[bin_index].index++;
	
	      if(use_mailman==true)
		  allgen_mail[bin_index].index++;
	       else
	           allgen[bin_index].index++;		
	}

    

   }        
	
	sum = 0 ;
        delete[] gtype;
}




void read_bed2 (std::istream& ifs,bool allow_missing,int num_snp)  {
         //ifstream ifs (filename.c_str(), ios::in|ios::binary);
        char magic[3];
        set_metadata ();

    gtype =  new unsigned char[ncol];

     if(read_header)
      binary_read(ifs,magic);
 
        int sum=0;
	double sum_dom=0;
        // Note that the coding of 0 and 2 can get flipped relative to plink because plink uses allele frequency (minor)
        // allele to code a SNP as 0 or 1.
        // This flipping does not matter for results.
        int y[4];

int bin_pointer;

vector<int> pointer_bins;

for(int i=0;i<num_snp;i++){
                global_snp_index++;
                ifs.read (reinterpret_cast<char*>(gtype), ncol*sizeof(unsigned char));
                float p_j = get_observed_pj(gtype);
	 

	   pointer_bins.clear();      
       	  for(int bin_index=0;bin_index<Nbin;bin_index++)
       		if(annot_bool[global_snp_index][bin_index]==1)
			  pointer_bins.push_back(bin_index);
			//bin_pointer=bin_index;

	  for (int k = 0 ;k < ncol ; k++) {
                unsigned char c = gtype [k];
                        // Extract PLINK genotypes
                y[0] = (c)&mask2;
                y[1] = (c>>2)&mask2;
                y[2] = (c>>4)&mask2;
                y[3] = (c>>6)&mask2;
                        int j0 = k * unitsperword;
                        // Handle number of individuals not being a multiple of 4
                        int lmax = 4;
                        if (k == ncol - 1)  {
                                lmax = Nindv%4;
                                lmax = (lmax==0)?4:lmax;
                        }
                        for ( int l = 0 ; l < lmax; l++){
                                int j = j0 + l ;
                                // Extract  PLINK coded genotype and convert into 0/1/2
                                // PLINK coding: 
                                // 00->0
                                // 01->missing
                                // 10->1
                                // 11->2
                                double val = y[l];
				double val_dom;
                                if(val==1 && !allow_missing){
                                        
					if(use_mailman==false & mean_impute==true){
						miss_index.push_back(j);
						val=0;
					}
					else{
					val = simulate2_geno_from_random(p_j);
                                        val++;
                                        val = (val==1) ? 0 : val;
                                        //val=0;
                                        }

                                 }
                                val-- ;
                                val =  (val < 0 ) ? 0 :val ;
                                
				//	/cout<<val<<endl;
				if(only_dominance==true){
					if(val==1)
					   val=0;
				}
				
			         sum += val;


				if(dominance==true){
				    if (val==1)
					val_dom=2*global_MAF(global_snp_index,0);
				    else if (val==2)
					val_dom=(4*global_MAF(global_snp_index,0))-2;
				    else if (val==0)
					  val_dom=val;

				  sum_dom+=val_dom;
 				}
						
                            for(int bin_index=0;bin_index<pointer_bins.size();bin_index++){
					
				       bin_pointer=pointer_bins[bin_index];
                                      int snp_index;

                                     if(use_mailman==true){
                                         snp_index=allgen_mail[bin_pointer].index;
                                         int horiz_seg_no = snp_index/allgen_mail[bin_pointer].segment_size_hori;
                                         allgen_mail[bin_pointer].p[horiz_seg_no][j] = 3 *allgen_mail[bin_pointer].p[horiz_seg_no][j]  + val;
                                     // computing sum for every snp to compute mean
                                         allgen_mail[bin_pointer].columnsum[snp_index]+=val;

                                      }
                                     else{
                                         snp_index=allgen[bin_pointer].index;
					 if(dominance==false)
                                         	allgen[bin_pointer].gen(snp_index,j)=val;
					 else{
						allgen[bin_pointer].gen(snp_index,j)=val_dom;
						allgen[bin_pointer].maf(snp_index,0)+=(double)val/Nindv/2;
					
					}	

                                     }


                            }

                    }
        }

    
      if(use_mailman==false & mean_impute==true){
                 for(int bin_index=0;bin_index<pointer_bins.size();bin_index++){
			int bin_pointer=pointer_bins[bin_index]; 
			int snp_index=allgen[bin_pointer].index;

			for(int indv_index=0;indv_index<miss_index.size();indv_index++){
				int temp=miss_index[indv_index];
				if(dominance==false) 
					allgen[bin_pointer].gen(snp_index,temp)=(double)sum/Nindv;
        			else
					allgen[bin_pointer].gen(snp_index,temp)=(double)sum_dom/Nindv;

			}
			
	         }
      }

	std::vector<int>().swap(miss_index);


     

     for(int bin_index=0;bin_index<pointer_bins.size();bin_index++){
     		bin_pointer=pointer_bins[bin_index];
	         if(use_mailman==true)
                  allgen_mail[bin_pointer].index++;
               else
                   allgen[bin_pointer].index++;
     }	

   if(dominance==false)
  	global_MAF(global_snp_index,0)=(double)sum/Nindv/2;

   




  sum=0;
  sum_dom=0;
}

      //  sum = 0 ;
        delete[] gtype;
}



void read_bed_1colannot (std::istream& ifs,bool allow_missing,int num_snp)  {
         //ifstream ifs (filename.c_str(), ios::in|ios::binary);
        char magic[3];
        set_metadata ();

    gtype =  new unsigned char[ncol];

     if(read_header)
      binary_read(ifs,magic);

        int sum=0;

        // Note that the coding of 0 and 2 can get flipped relative to plink because plink uses allele frequency (minor)
        // allele to code a SNP as 0 or 1.
        // This flipping does not matter for results.
        int y[4];

int bin_pointer;

//vector<int> pointer_bins;

for(int i=0;i<num_snp;i++){
                global_snp_index++;
                ifs.read (reinterpret_cast<char*>(gtype), ncol*sizeof(unsigned char));
                float p_j = get_observed_pj(gtype);

          for (int k = 0 ;k < ncol ; k++) {
                unsigned char c = gtype [k];
                        // Extract PLINK genotypes
                y[0] = (c)&mask2;
                y[1] = (c>>2)&mask2;
                y[2] = (c>>4)&mask2;
                y[3] = (c>>6)&mask2;
                        int j0 = k * unitsperword;
                        // Handle number of individuals not being a multiple of 4
                        int lmax = 4;
                        if (k == ncol - 1)  {
                                lmax = Nindv%4;
                                lmax = (lmax==0)?4:lmax;
                        }
                        for ( int l = 0 ; l < lmax; l++){
                                int j = j0 + l ;
                                // Extract  PLINK coded genotype and convert into 0/1/2
                                // PLINK coding: 
                                // 00->0
                                // 01->missing
                                // 10->1
                                // 11->2
                                int val = y[l];
                                if(val==1 && !allow_missing){
                                        val = simulate2_geno_from_random(p_j);
                                        val++;
                                        val = (val==1) ? 0 : val;
                                   //val=0;
                                 }
                                val-- ;
                                val =  (val < 0 ) ? 0 :val ;
                                sum += val;
//cout<<"sss"<<endl;
//cout<<global_snp_index<<endl;

                                       bin_pointer=SNP_annot[global_snp_index]-1;
                           
//	cout<<bin_pointer<<endl;
			           int snp_index;
                                     if(use_mailman==true){
                                         snp_index=allgen_mail[bin_pointer].index;
                                         int horiz_seg_no = snp_index/allgen_mail[bin_pointer].segment_size_hori;
                                         allgen_mail[bin_pointer].p[horiz_seg_no][j] = 3 *allgen_mail[bin_pointer].p[horiz_seg_no][j]  + val;
                                     // computing sum for every snp to compute mean
                                         allgen_mail[bin_pointer].columnsum[snp_index]+=val;

                                      }
                                     else{
                                         snp_index=allgen[bin_pointer].index;
                                         allgen[bin_pointer].gen(snp_index,j)=val;
                                     }


                    }
        }

        	bin_pointer=SNP_annot[global_snp_index]-1;
	         if(use_mailman==true)
                  allgen_mail[bin_pointer].index++;
               else
                   allgen[bin_pointer].index++;
}

  sum = 0 ;
  delete[] gtype;



}

MatrixXdr jack_se(MatrixXdr jack){

int nrows=jack.rows();
int ncols=jack.cols();
MatrixXdr sum_row=jack.rowwise().mean();
MatrixXdr SEjack;
SEjack=MatrixXdr::Zero(nrows,1);
double temp_val=0;
for (int i=0;i<nrows;i++){
    for (int j=0;j<ncols;j++){
        temp_val=jack(i,j)-sum_row(i);
        temp_val= temp_val* temp_val;
        SEjack(i,0)+=temp_val;
    }
    SEjack(i,0)=SEjack(i,0)*(Njack-1)/Njack;
    SEjack(i,0)=sqrt(SEjack(i,0));
}


return SEjack;
}


//MatrixXdr maf_ld;
//MatrixXdr A_effect;
void read_mafld(std::string filename){
	
	maf_ld.resize(Nsnp,3);
	ifstream ifs(filename.c_str(), ios::in);

        std::string line;
        std::istringstream in;
	string temp;

	std::getline(ifs,line);
	in.str(line);
	cout<<"reading MAF and LD scores"<<endl;
	while(in>>temp){
		cout<<temp<<" ";
	}
	 
	cout<<endl;
	double val;
	int i=0;
	while(std::getline(ifs, line)){
		  in.clear();
                  in.str(line);
		  
		  in >> temp;  
		  val = atof(temp.c_str());
		  maf_ld(i,0)=val*(1-val);  ///MAF(1-MAF)
		  maf_ld(i,2)=val;		  
  
		  in >> temp; 
		  val = atof(temp.c_str());
                  maf_ld(i,1)=val;

		  i++;
	}
	
}

void read_simul_par(std::string filename){

//double p_casual,double ld_ex,double maf_ex,double min_maf,double max_maf,double h2,int num
        
	int num_par=7;
	int num_vc=2;
	simul_par.resize(num_par,num_vc);
        ifstream ifs(filename.c_str(), ios::in);

        std::string line;
        std::istringstream in;
        string temp;

        std::getline(ifs,line);
        in.str(line);
	cout<<"reading simulation parameters"<<endl;
        while(in>>temp){
                cout<<temp<<" ";
        }
	cout<<endl;
        double val;
        int i=0;
        while(std::getline(ifs, line)){
                  in.clear();
                  in.str(line);

		  for (int j=0;j<num_par;j++){
                  in >> temp;
                  val = atof(temp.c_str());
                  simul_par(j,i)=val;
		  cout<<val<<" ";
		  }
		  cout<<endl;
		  i++;
        }
	cout<<endl;

}



void make_maf_ld_effect(double p_casual,double ld_ex,double maf_ex,double min_maf,double max_maf,double h2,int num){

    std::ofstream outfile;
	if(dominance==true){
	outfile.open("causal_snp_dom.txt", std::ios_base::out);
	}
	else{
        outfile.open(command_line_opts.OUTPUT_FILE_PATH+".causal_snps_annot.txt", std::ios_base::out);
	}
	//MatrixXdr A_effect;
	A_effect.resize(Nsnp,1);
        A_effect_multi.resize(Nsnp,num);

    // method 1: toss a coin (varying # of causals)
	boost::mt19937 seedr;
    //seedr.seed(std::time(0));
    //seed with rdtsc:
    seedr.seed(rdtsc());
    boost::bernoulli_distribution<> b_dist(p_casual);
    boost::variate_generator<boost::mt19937&, boost::bernoulli_distribution<>  > coin(seedr, b_dist);


	MatrixXdr cas;
	cas.resize(Nsnp,1);

    bool fix_ncausals = true;
    if (!fix_ncausals){
	    for (int i=0;i<Nsnp;i++){
	        cas(i,0)=coin();
	        if(maf_ld(i,2)<min_maf || maf_ld(i,2)>max_maf)
		        cas(i,0)=0;
	    }
    }
    else{
        // method 2: random selection (fixed # of causals)
        boost::uniform_int<> idxs(0, Nsnp-1);
        boost::variate_generator<boost::mt19937&, boost::uniform_int<> > get_snp_idx(seedr, idxs);
        int ncausals = p_casual*Nsnp;
        for (int i=0; i<Nsnp; i++)
            cas(i,0)=0;
        while (ncausals != 0){
            int tmp = get_snp_idx();
            if ((cas(tmp, 0) == 0) && (maf_ld(tmp,2) >= min_maf) && (maf_ld(tmp,2) <= max_maf) ){
                cas(tmp,0) = 1;
                ncausals--;
            }
        }
    }

	cout<<cas.sum()<<endl;


	///write causal snps in files 
	outfile<<cas<<endl;
	outfile.close(); 
	///
	
	maf_ld.col(0)=maf_ld.col(0).array().pow(maf_ex);
	maf_ld.col(1)=maf_ld.col(1).array().pow(ld_ex);

	/*cout<<"f(1-f)"<<maf_ld(0,0)<<endl;
	cout<<"ld"<<maf_ld(0,1)<<endl;
	cout<<"c"<<cas(0,0)<<endl;
	*/

	A_effect=maf_ld.col(0).array()*maf_ld.col(1).array()*cas.array();
	
	//cout<<"c_ld_f"<<effect(0,0)<<endl;

	//double cons_fact=(double)h2/A_effect.sum();
    
    // Fix per-SNP heritability instead of total h2
	double cons_fact=(double)h2/Nsnp;
    
    cout << "h2: " << h2 << endl
        << "cons_fact: " << cons_fact << endl
        << "A_effect.sum(): " << A_effect.sum() << endl;

	//cout<<"const"<<cons_fact<<endl;
	
	A_effect=A_effect*cons_fact;



	//cout<<"normalized eff"<<effect(0,0)<<endl;
	cout<<"sum  "<<A_effect.sum()<<endl;

	 for (int i=0;i<Nsnp;i++){
		boost::normal_distribution<> dist(0,sqrt(A_effect(i,0)));
        	boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > effect_size(seedr, dist);
	//	cout<<i<<" "<<sqrt(A_effect(i,0))<<endl;
		for (int j=0;j<num;j++){
		A_effect_multi(i,j)=effect_size();

		}
	}

	
				

       

}

void make_dom_effect(double p_casual,double h2d,int num){

	boost::mt19937 seedr;
        seedr.seed(std::time(0));
        boost::bernoulli_distribution<> b_dist(p_casual);
        boost::variate_generator<boost::mt19937&, boost::bernoulli_distribution<>  > coin(seedr, b_dist);





}

void  read_nongen(std::string filename,int num_indv){
        ifstream ifs(filename.c_str(), ios::in);
        std::string line;
        std::istringstream in;
        
        int j=0;
        while(std::getline(ifs, line)){
                in.clear();
                in.str(line);
                string temp;

		for (int i=0; i<num_indv;i++){
                in>>temp;
                        /*if(temp=="NA")
                        {
                                missing[k].push_back(j);
                                continue;
                        }*/
		double cur = atof(temp.c_str());
		//nongen_vc(j,i)=cur;
		gen(j,i)=cur;
               }
        j++;
        }

}










MatrixXdr  compute_XXy (int num_snp,MatrixXdr pheno){
                
	res.resize(num_snp, 1);
        if(use_mailman==true)
                multiply_y_pre_fast(pheno,1,res, false);
        else
                res=gen*pheno;


        MatrixXdr zb_sum = pheno.colwise().sum();


        for(int j=0; j<num_snp; j++)
            for(int k=0; k<1;k++)
                res(j,k) = res(j,k)*stds(j,0);

        MatrixXdr resid(num_snp, 1);
        MatrixXdr inter = means.cwiseProduct(stds);
        resid = inter * zb_sum;
        MatrixXdr inter_zb = res - resid;


        for(int k=0; k<1; k++)
            for(int j=0; j<num_snp;j++)
                inter_zb(j,k) =inter_zb(j,k) *stds(j,0);
       MatrixXdr new_zb = inter_zb.transpose();
       MatrixXdr new_res(1, Nindv);


         if(use_mailman==true)
            multiply_y_post_fast(new_zb, 1, new_res, false);
         else
            new_res=new_zb*gen;

       MatrixXdr new_resid(1, num_snp);
       MatrixXdr zb_scale_sum = new_zb * means;
       new_resid = zb_scale_sum * MatrixXdr::Constant(1,Nindv, 1);


                      /// new zb 
       MatrixXdr temp=new_res - new_resid;

        for (int i=0;i<1;i++)
           for(int j=0;j<Nindv;j++)
                 temp(i,j)=temp(i,j)*mask(j,0);


        return temp.transpose();


}

MatrixXdr compute_Xtw(int num_snp,MatrixXdr w){

        MatrixXdr res(num_snp, 1);


        if(use_mailman==true)
                multiply_y_pre_fast(w,1,res,false);
        else
                 res=gen*w;

        double w_sum=w.sum();        

        res = res.cwiseProduct(stds);
        MatrixXdr resid(num_snp, 1);
        resid = means.cwiseProduct(stds);
        resid = resid *w_sum;
        MatrixXdr Xy(num_snp,1);
        Xy = res-resid;

        return Xy;

}






MatrixXdr  compute_Xv (int num_snp,MatrixXdr new_zb){

	
	

	Nz=num_simul;
	for (int j=0;j<Nz;j++)
        for (int i=0;i<num_snp;i++)
                new_zb(j,i)= new_zb(j,i)*stds(i,0);


         MatrixXdr new_res(Nz, Nindv);


        if(use_mailman==true)
                multiply_y_post_fast(new_zb,Nz,new_res, false);
        else
                new_res=new_zb*gen;


        MatrixXdr new_resid(Nz, num_snp);
       MatrixXdr zb_scale_sum = new_zb * means;
       new_resid = zb_scale_sum * MatrixXdr::Constant(1,Nindv, 1);


                      /// new zb 
       MatrixXdr temp=new_res - new_resid;




        return temp.transpose();
}










void process(string name, int round){


ifstream ifs (name.c_str(), ios::in|ios::binary);
read_header=true;
global_snp_index=-1;

//////////////////////////


int point_index=0;
for (int jack_index=0;jack_index<Njack;jack_index++){

	//cout<<"jack index"<<jack_index<<endl;
        int read_Nsnp=(jack_index<(Njack-1)) ? (step_size) : (step_size+step_size_rem);

       if(use_mailman==true){
        for (int i=0;i<Nbin;i++){
        allgen_mail[i].segment_size_hori = floor(log(Nindv)/log(3)) - 2 ;
        allgen_mail[i].Nsegments_hori = ceil(jack_bin[jack_index][i]*1.0/(allgen_mail[i].segment_size_hori*1.0));
        allgen_mail[i].p.resize(allgen_mail[i].Nsegments_hori,std::vector<int>(Nindv));
        allgen_mail[i].not_O_i.resize(jack_bin[jack_index][i]);
        allgen_mail[i].not_O_j.resize(Nindv);
        allgen_mail[i].index=0;
        allgen_mail[i].Nsnp=jack_bin[jack_index][i];
        allgen_mail[i].Nindv=Nindv;

         allgen_mail[i].columnsum.resize(jack_bin[jack_index][i],1);
          for (int index_temp=0;index_temp<jack_bin[jack_index][i];index_temp++)
                    allgen_mail[i].columnsum[index_temp]=0;

         }
       }
       else{
           for (int bin_index=0;bin_index<Nbin;bin_index++){
                allgen[bin_index].gen.resize(jack_bin[jack_index][bin_index],Nindv);
                allgen[bin_index].index=0;
		//allgen[bin_index].maf.resize(jack_bin[jack_index][bin_index],1);
		allgen[bin_index].maf=MatrixXdr::Zero(jack_bin[jack_index][bin_index],1);
           }
       }


        if(use_1col_annot==true)
                read_bed_1colannot(ifs,missing,read_Nsnp);
        else
                read_bed2(ifs,missing,read_Nsnp);
       read_header=false;
     for (int bin_index=0;bin_index<Nbin;bin_index++){
          int num_snp;
          if (use_mailman==true)
                num_snp=allgen_mail[bin_index].index;
          else
                num_snp=allgen[bin_index].index;


          if(num_snp!=0){
          stds.resize(num_snp,1);
          means.resize(num_snp,1);
	  Eigen::VectorXd mean;
          if(use_mailman==true){
                for (int i=0;i<num_snp;i++)
                   means(i,0)=(double)allgen_mail[bin_index].columnsum[i]/Nindv;
         }
          else{
		 means=allgen[bin_index].gen.rowwise().mean();
                mean=allgen[bin_index].gen.rowwise().mean();
	  }
	 
         if(mean_impute==true & use_mailman==false){		
		 
		 stds=1/((allgen[bin_index].gen.colwise()-mean).array().square().rowwise().sum()/(allgen[bin_index].gen.cols()-1)).sqrt();
         }
	 else{

		if(dominance==false){
	 		 for (int i=0;i<num_snp;i++)
               			stds(i,0)=1/sqrt((means(i,0)*(1-(0.5*means(i,0)))));
		}  
		else{
			for (int i=0;i<num_snp;i++){
				double mf=allgen[bin_index].maf(i,0);
				stds(i,0)=1/((2*mf)*(1-mf));
			//cout<<stds(i,0)<<endl;
			}
		}
        }
 
	if (use_mailman==true){
                g=allgen_mail[bin_index];
                 g.segment_size_hori = floor(log(Nindv)/log(3)) - 2 ;
                 g.Nsegments_hori = ceil(jack_bin[jack_index][bin_index]*1.0/(g.segment_size_hori*1.0));
                 g.p.resize(g.Nsegments_hori,std::vector<int>(Nindv));
                 g.not_O_i.resize(jack_bin[jack_index][bin_index]);
                 g.not_O_j.resize(Nindv);
                initial_var();
           }
          else{
                 gen=allgen[bin_index].gen;

          }



////////////////////

/// save dominance grm
/*      if(dominance==true){
	MatrixXdr grm= gen*gen.transpose();
	 std::ofstream outfile;
        string add_output;
        add_output="dominance_smallscale.grm.gz";
        outfile.open(add_output.c_str(), std::ios_base::out);
	      for (int i=0;i<Nindv;i++)
                for (int j=0;j<=i;j++){
                        double grmval=(double)grm(i,j)/Nsnp;
                        outfile<<i+1<<" "<<j+1<<" "<<Nsnp<<" "<<grmval<<endl;
                }
     }*/	 
///////




	MatrixXdr vec_eff=A_effect_multi.block(point_index,0,num_snp,num_simul).transpose();
	point_index=point_index+num_snp;
	simul_pheno+=compute_Xv(num_snp,vec_eff);



/////////////////



        //cout<<num_snp<< "SNPs in bin "<<bin_index<<"of jack "<<jack_index<<endl;
        //cout<<" Reading and computing bin "<<bin_index <<"  of "<< jack_index<<"-th is finished"<<endl;

            if(use_mailman==true){
                delete[] sum_op;
                delete[] partialsums;
                 delete[] yint_e;
                delete[] yint_m;
                for (int i  = 0 ; i < hsegsize; i++)
                        delete[] y_m [i];
                delete[] y_m;

                for (int i  = 0 ; i < g.Nindv; i++)
                        delete[] y_e[i];
                delete[] y_e;

                std::vector< std::vector<int> >().swap(g.p);
                std::vector< std::vector<int> >().swap(g.not_O_j);
                std::vector< std::vector<int> >().swap(g.not_O_i);
                std::vector< std::vector<int> >().swap(allgen_mail[bin_index].p);
                std::vector< std::vector<int> >().swap(allgen_mail[bin_index].not_O_j);
                std::vector< std::vector<int> >().swap(allgen_mail[bin_index].not_O_i);
        //g.p.clear();
        //g.not_O_j.clear();
        //g.not_O_i.clear();
                g.columnsum.clear();
                g.columnsum2.clear();
                g.columnmeans.clear();
                 g.columnmeans2.clear();
                 allgen_mail[bin_index].columnsum.clear();
                allgen_mail[bin_index].columnsum2.clear();
                allgen_mail[bin_index].columnmeans.clear();
                 allgen_mail[bin_index].columnmeans2.clear();
            }
        }

     }
//cout<<" Reading and computing  of "<< jack_index<<"-th is finished"<<endl;
}


}







int main(int argc, char const *argv[]){
 

	
	parse_args(argc,argv);
////////////////////////////////////////////
///////////////////////////////////////////
    
        int B = command_line_opts.batchNum;
        k_orig = command_line_opts.num_of_evec ;
        debug = command_line_opts.debugmode ;
        check_accuracy = command_line_opts.getaccuracy;
        var_normalize = false;
        accelerated_em = command_line_opts.accelerated_em;
        k = k_orig + command_line_opts.l;
        k = (int)ceil(k/10.0)*10;
        command_line_opts.l = k - k_orig;
        srand((unsigned int) time(0));
	Nz=command_line_opts.num_of_evec;
        k=Nz;

	Njack=command_line_opts.jack_number;

////
string filename;
//////////////////////////// Read multi genotypes
string line;
int cov_num;
int num_files=0;
string geno_name=command_line_opts.GENOTYPE_FILE_PATH;

/////////Read bim file to count # SNPs
std::stringstream f1;
f1 << geno_name << ".bim";
read_bim (f1.str());

//////Read annotation files
filename=command_line_opts.Annot_PATH;

if(use_1col_annot==true){
read_annot_1col(filename);
}
else{
read_annot(filename);
}
//filename=command_line_opts.Annot_PATH;

///reading phnotype and save the number of indvs
filename=command_line_opts.PHENOTYPE_FILE_PATH;
count_pheno(filename);

std::stringstream f0;
f0 << geno_name << ".fam";
string name_fam=f0.str();
int fam_lines=count_fam(name_fam);

create_simul_pheno(name_fam);


cout<<"Number of Indvs :"<<Nindv<<endl;


///// code for handeling overlapping annotations
std::stringstream f3;
f3 << geno_name << ".bed";
string name=f3.str();
cout<<name<<endl;
ifstream ifs (name.c_str(), ios::in|ios::binary);
read_header=true;
global_snp_index=-1;

////read maf ld file

filename=command_line_opts.maf_ld_PATH;
read_mafld(filename);


global_MAF.resize(Nsnp,1);

allgen_mail.resize(Nbin);  
//allgen.resize(Nbin);

double h2A=0;
double h2D=0;
filename=command_line_opts.simul_par_PATH;
read_simul_par(filename);
make_maf_ld_effect(simul_par(0,0),simul_par(1,0),simul_par(2,0),simul_par(3,0),simul_par(4,0),simul_par(5,0),simul_par(6,0));

num_simul=simul_par(6,0);

Nz=num_simul;
k=Nz;
cout<<"number of simulations: "<<num_simul<<endl;

simul_pheno=MatrixXdr::Zero(Nindv,num_simul);

use_mailman=true;
dominance=false;
process(name,1);
h2A=simul_par(5,0);
cout << "SIMUL PAR, H2A: " << h2A << endl;


bool simul_dom=false;

if(simul_dom==true){
dominance=true;

make_maf_ld_effect(simul_par(0,1),simul_par(1,1),simul_par(2,1),simul_par(3,1),simul_par(4,1),simul_par(5,1),simul_par(6,1));

//dominance=true;
process(name,1);
 h2D=simul_par(5,1);
}





boost::mt19937 seedr;
boost::normal_distribution<> dist(0,sqrt(1-h2A-h2D));
boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > effect_env(seedr, dist);



std::ofstream outfile;                                                                  
string add_output=command_line_opts.OUTPUT_FILE_PATH;                                   

cout<<"writing phenos in "<<add_output<<endl;

for ( int i=0;i<num_simul;i++){
	cout<<"writing "<<i<<"-th simulation"<<endl;

	string index=boost::lexical_cast<string>(i);
	//string name=add_output+index+".phen";
    // for varying causal SNPs:
	string name=add_output+".phen";
	std::ofstream outfile;
	outfile.open(name, std::ios_base::out);
	
	cout<<name<<endl;
	outfile<<"FID IID pheno"<<endl;
	for (int j=0;j<Nindv;j++){
		outfile<<fid[j]<<" "<<fid[j]<<" "<<simul_pheno(j,i)+effect_env()<<endl;		
	}
	outfile.close();
}


/*
cout<<"xxxxxxxxx"<<endl;
use_mailman=false;
dominance=true;
process(name,1);
*/


return 0;

}
