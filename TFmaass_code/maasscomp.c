#include <acb_mat.h>
#include <acb_hypgeom.h>
#include <inttypes.h>
#include <sys/file.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdio.h>
#include <dirent.h>

#define MAXLEVEL 105
#define MAXN 4000000
#define NPROCS 64
#define MATDIM 2000 //Dimension of the matricies used.

static slong prec = 300;
#include "testfunc.c"

int initfactor64(const char *);
int factor64(uint64_t *,int *,uint64_t);
int isprime64(uint64_t);

typedef struct {
	long N; //This denotes the level
	arb_mat_t posterms; //Stores the trace formula values of positive Hecke operators and level N
	arb_mat_t negterms; //Stores the trace formula values of negative Hecke operators and level N
} TFdata;

int FindIndex(const long *a, int size, long value){
    int i = 0;
    while ( i < size && a[i] != value ) i++;
    return ( i == size ? (-1) : (i) );
}

// Read trace formula data values from tracecomp.c
void readtfdata(TFdata tracevals){
    char filename[30],line[3600],gstr0[200],gstr2[200],gstr4[200];;
    sprintf(filename,"./tfdata/%ld.txt",tracevals.N);
    FILE *fp = fopen(filename,"r");
    arb_t gval0,gval2,gval4;
    arb_init(gval0); arb_init(gval2); arb_init(gval4);
    long N,n;
    while(fgets(line,sizeof(line),fp) != NULL){
        sscanf(line,"%li,%li,%100[^,],%100[^,],%100[^\n]",&N,&n,gstr0,gstr2,gstr4);
        arb_load_str(gval0,gstr0);
        arb_load_str(gval2,gstr2);
        arb_load_str(gval4,gstr4);
        if(labs(n) > MAXN){
            continue;
        }
        if(n > 0){
            arb_set(arb_mat_entry(tracevals.posterms,n-1,0),gval0);
            arb_set(arb_mat_entry(tracevals.posterms,n-1,1),gval2);
            arb_set(arb_mat_entry(tracevals.posterms,n-1,2),gval4);
        } else if(n < 0){
            arb_set(arb_mat_entry(tracevals.negterms,labs(n)-1,0),gval0);
            arb_set(arb_mat_entry(tracevals.negterms,labs(n)-1,1),gval2);
            arb_set(arb_mat_entry(tracevals.negterms,labs(n)-1,2),gval4);
        }
        
    }
    fclose(fp);
}

// Convert acb_mat to arb_mat.
void acb_mat_to_arb_mat(arb_mat_t res, acb_mat_t src){
    for(int i = 0; i < acb_mat_nrows(res);i++){
        for(int j = 0; j < acb_mat_ncols(res);j++){
            arb_set(arb_mat_entry(res,i,j),acb_realref(acb_mat_entry(src,i,j)));
        }
    }
}

// Reduce the size of a matrix src and store it in res.
void arb_mat_set_window(arb_mat_t res, arb_mat_t src, slong r1, slong c1, slong r2, slong c2){
    slong i,j;
    if(arb_mat_nrows(res) == r2-r1 && arb_mat_ncols(res) == c2-c1){
        for(i=0; i<r2-r1;i++){
            for(j=0; j<c2-c1;j++){
                arb_set(arb_mat_entry(res,i,j),arb_mat_entry(src,r1+i,c1+j));
            }
        }
    } else {
        printf("Wrong dimensions for res matrix in arb_mat_set_window\n");
        return;
    }
}

// Order the eigenvalues and eigenvectors in place in ascending order.
void ordereigenvalvecasc(acb_ptr E, acb_mat_t P){
    acb_t e;
    int m,n = acb_mat_ncols(P);
    acb_init(e);
    acb_mat_transpose(P,P);

    for(int k = 0; k < n-1; k++){
        m = k;
        for(int l = k+1; l < n; l++){
            if(arb_lt(acb_realref(E+l),acb_realref(E+m))){
                m = l;
            }
        }
        if(k != m){
            acb_swap(E+k,E+m);
            acb_mat_swap_rows(P,NULL,k,m);
        }
    }
    acb_mat_transpose(P,P);
}

// Computes the Kbessel function with rigorous error bound for terms near zero.
void acb_my_Kbessel(acb_t res, const acb_t nu, const acb_t z, slong prec){
    static int init;
	static arb_t gammathird,R,b1,b2,b3,pi2;
	static acb_t numid;
	int k,n;

	if (!init) {
        arb_init(gammathird); arb_init(R); arb_init(pi2);
        arb_init(b1); arb_init(b2); arb_init(b3);
        acb_init(numid);
        arb_one(gammathird);
        arb_div_ui(gammathird, gammathird, 3, prec);
        arb_gamma(gammathird, gammathird,prec);
        arb_const_pi(pi2,prec);
        arb_mul_2exp_si(pi2,pi2,-1);
		init = 1;
	}

    arb_set(R,acb_imagref(nu));
    if(arb_ge(acb_realref(z),R)){

        acb_get_mid(numid,nu);
        acb_hypgeom_bessel_k(res,numid,z,prec);

        // Now compute the error bound, even though we are using acb, the input will always be real.
        arb_mul(b1,R,R,prec);
        arb_mul(b2,acb_realref(z),acb_realref(z),prec);
        arb_sub(b2,b2,b1,prec);
        arb_sqrt(b2,b2,prec);

        arb_mul_ui(b3,pi2,3,prec);
        arb_div(b3,b3,b2,prec);
        arb_sqrt(b3,b3,prec);

        arb_set_ui(b1,3);
        arb_mul_2exp_si(b1,b1,-2);
        arb_div(b1,b1,R,prec);
        arb_root_ui(b1,b1,3,prec);
        arb_mul(b1,b1,gammathird,prec);

        if(arb_lt(b3,b1)){
            arb_set(b1,b3);
        }

        // Still have to compute the exp terms in error term.
        arb_div(b3,R,acb_realref(z),prec);
        arb_acos(b3,b3,prec);
        arb_sub(b3,b3,pi2,prec);
        arb_mul(b3,b3,R,prec);
        arb_sub(b3,b3,b2,prec);
        arb_exp(b3,b3,prec);
        arb_mul(b1,b1,b3,prec);

        //Now unique to this code, we know the worst error is going to occur from the error bound of R. So using
        //the above computed bound and the intermediate value theorem we can get the error radius.
        arb_get_rad_arb(b2,R);
        arb_mul(b1,b1,b2,prec);
        acb_add_error_arf(res,arb_midref(b1));
    } else {
        acb_hypgeom_bessel_k(res,nu,z,prec);
    }
}

//  Performs all the steps to compute and verify the Laplace and Hecke eigenvalues for a given level and parity.
// The parity input is 0 for even forms and 1 for odd forms.
void verifarb(arb_mat_t gvals,int N,int level,int parity){
    arb_t b1,b2,b3,b4,b5;
    acb_t z1,z2,z3,z4;
    arb_mat_t Mg,Mg2,Mg4,D,P,R,Rinv,RinvT,tempvec,tempvecT,mat11,tempmat,Drtinv,ci,R_vals,coeffs;
    int d,i,j,g,m,n,sign,fouriercheck;
    uint64_t pN[15],p[15];
    int eN[15],e[15],kN,k;
    int x,y,omega;
    arb_init(b1); arb_init(b2); arb_init(b3); arb_init(b4); arb_init(b5);
    acb_init(z1); acb_init(z2); acb_init(z3); acb_init(z4);
    arb_mat_init(Mg,N,N);
    arb_mat_init(Mg2,N,N);
    arb_mat_init(Mg4,N,N);

    arb_mat_zero(Mg);
    arb_mat_zero(Mg2);
    arb_mat_zero(Mg4);

    char direc[25],filename[100];
    sprintf(direc,"./maassdata/");
    FILE *fp;

    kN = factor64(pN,eN,level);
    int bindigits[kN+1];

    // Generate the matrices of trace formula values Q and \tilde{Q}.
    i=1;
    j=1;
    for(m=1;m<=N;m++){
        for(n=1;n<=N;n++){
            if(m!=n){
                if(n_gcd(m,level)==1 && n_gcd(n,level)==1){
                    arb_zero(b1);
                    arb_zero(b2);
                    arb_zero(b3);
                    g = n_gcd(m,n);
                    for(d=1;d<=g;d++){
                        if(g % d == 0){
                            arb_add(b1,b1,arb_mat_entry(gvals,m*n/(d*d)-1,0),prec);
                            arb_add(b2,b2,arb_mat_entry(gvals,m*n/(d*d)-1,1),prec);
                            arb_add(b3,b3,arb_mat_entry(gvals,m*n/(d*d)-1,2),prec);
                        }
                    }
                    arb_set(arb_mat_entry(Mg,i-1,j-1),b1);
                    arb_set(arb_mat_entry(Mg2,i-1,j-1),b2);
                    arb_set(arb_mat_entry(Mg4,i-1,j-1),b3);
                    j++;
                }
            }else{
                if(n_gcd(m,level)==1){
                    arb_zero(b1);
                    arb_zero(b2);
                    arb_zero(b3);
                    for(d=1;d<=m;d++){
                        if(m % d == 0){
                            arb_add(b1,b1,arb_mat_entry(gvals,m*m/(d*d)-1,0),prec);
                            arb_add(b2,b2,arb_mat_entry(gvals,m*m/(d*d)-1,1),prec);
                            arb_add(b3,b3,arb_mat_entry(gvals,m*m/(d*d)-1,2),prec);
                        }
                    }
                    arb_set(arb_mat_entry(Mg,i-1,i-1),b1);
                    arb_set(arb_mat_entry(Mg2,i-1,i-1),b2);
                    arb_set(arb_mat_entry(Mg4,i-1,i-1),b3);
                    j++;
                }
            }
        }
        j=1;
        if(n_gcd(m,level)==1){
            i++;
        }
    }
    i=i-1;
    arb_mat_t M1,M2,M4;
    arb_mat_init(M1,i,i);
    arb_mat_init(M2,i,i);
    arb_mat_init(M4,i,i);

    arb_mat_set_window(M1,Mg,0,0,i,i);
    arb_mat_set_window(M2,Mg2,0,0,i,i);
    arb_mat_set_window(M4,Mg4,0,0,i,i);

    arb_mat_clear(Mg);
    arb_mat_clear(Mg2);
    arb_mat_clear(Mg4);

    arb_mat_init(ci,i,1);
    arb_mat_init(D,i,1);
    arb_mat_init(Drtinv,i,i);
    arb_mat_init(P,i,i);
    arb_mat_init(R,i,i);
    arb_mat_init(Rinv,i,i);
    arb_mat_init(RinvT,i,i);
    arb_mat_init(tempvec,i,1);
    arb_mat_init(tempvecT,1,i);
    arb_mat_init(mat11,1,1);
    arb_mat_init(tempmat,i,i);

    // Compute the \tilde{Q} matrix.
    arb_one(b1);
    arb_mul_2exp_si(b1,b1,-2);
    arb_mat_neg(tempmat,M2);
    arb_mat_scalar_addmul_arb(tempmat,M1,b1,prec);
    arb_mat_get_mid(tempmat,tempmat);

    /* Start of the linear algebra to compute the forms */
    acb_mat_t CM1,CP;
    acb_mat_init(CM1,i,i); acb_mat_init(CP,i,i);
    acb_mat_set_arb_mat(CM1,M1);
    acb_mat_get_mid(CM1,CM1);
    acb_ptr E;
    E = _acb_vec_init(i);

    // First diagonalisation of the matrix Q.
    acb_mat_approx_eig_qr(E,NULL,CP,CM1,NULL,0,prec);
    ordereigenvalvecasc(E,CP);
    acb_mat_to_arb_mat(P,CP);

    arb_t tol;
    arb_init(tol);
    arb_zero(tol);

    for(m=0;m<i;m++){
        if(arb_lt(acb_realref(E+m),tol)){
            arb_zero(arb_mat_entry(Drtinv,m,m));
        } else {
            arb_rsqrt(arb_mat_entry(Drtinv,m,m),acb_realref(E+m),prec);
        }
    }
    arb_mat_mul(RinvT,P,Drtinv,prec);
    arb_mat_transpose(Rinv,RinvT);
    arb_mat_get_mid(RinvT,RinvT);

    arb_mat_mul(tempmat,tempmat,RinvT,prec);
    arb_mat_mul(tempmat,Rinv,tempmat,prec);

    arb_mat_get_mid(tempmat,tempmat);
    acb_mat_set_arb_mat(CM1,tempmat);
    acb_mat_zero(CP);
    _acb_vec_zero(E,i);
    // Second diagonalisation of the matrix R^{-1} \tidle{Q} R.
    // This will give us the list of Laplace eigenvalues in E. 
    acb_mat_approx_eig_qr(E,NULL,CP,CM1,NULL,0,prec);

    // b2 = 1/4.
    arb_one(b2);
    arb_mul_2exp_si(b2,b2,-2);

    // Removing eigenvalues that are NaN.
    m=0;
    for(n=0;n<i;n++){
        if(arb_lt(acb_realref(E+n),b2)){
            arb_pos_inf(acb_realref(E+n));
            arb_zero(acb_imagref(E+n));
            continue;
        }
        m++;
    }
    ordereigenvalvecasc(E,CP);
    acb_mat_to_arb_mat(P,CP);
    arb_mat_get_mid(P,P);
    arb_mat_init(R_vals,m,2);
    arb_mat_init(coeffs,N,2);

    // Convert the Laplace eigenvalues lambda = 1/4 + R^2 into the values R.
    for(n=0;n<m;n++){
        arb_set(b1,acb_realref(E+n));  
        arb_sub(b1,b1,b2,prec);
        arb_sqrt(b1,b1,prec);
        mag_zero(arb_radref(b1));
        arb_set(arb_mat_entry(R_vals,n,0),b1);
    }
    /* Start of the linear algebra to compute the forms */

    // This is to ignore forms that we dont really verify.
    arb_set_d(tol,pow(10,-2));

    // Loop over each of the Laplace eigenvalues that we shall verify.
    for(n=0;n<m;n++){
        /* Start of Laplace eigenvalue verification */
        arb_mat_zero(coeffs);
        arb_set(b1,arb_mat_entry(R_vals,n,0));

        arb_mat_window_init(ci,P,0,n,i,n+1);
        arb_mat_mul(ci,RinvT,ci,prec);
        arb_mat_transpose(tempvecT,ci);

        arb_mat_set(tempmat,M4);
        arb_mul(b1,b1,b1,prec);//R^2
        arb_mul_2exp_si(b2,b1,1);
        arb_mat_scalar_addmul_arb(tempmat,M2,b2,prec);
        arb_mul(b1,b1,b1,prec);//R^4
        arb_mat_scalar_addmul_arb(tempmat,M1,b1,prec);

        arb_mat_mul(tempvec,M1,ci,prec);
        arb_mat_mul(mat11,tempvecT,tempvec,prec);

        arb_set(b2,arb_mat_entry(mat11,0,0));

        arb_mat_mul(tempvec,tempmat,ci,prec);
        arb_mat_mul(mat11,tempvecT,tempvec,prec);

        arb_div(b2,arb_mat_entry(mat11,0,0),b2,prec);
        arb_sqrtpos(b2,b2,prec);
        // This checks if error is > 10^(-2), i.e if its too big.
        if(arf_get_d(arb_midref(b2),ARF_RND_NEAR) > arf_get_d(arb_midref(tol),ARF_RND_NEAR) || arf_is_nan(arb_midref(b2))){
            continue;
        }
        arb_add_error(arb_mat_entry(R_vals,n,0),b2);
        arb_set(arb_mat_entry(R_vals,n,1),b2);
        /* End of Laplace eigenvalue verification */

        /* Start of Hecke eigenvalue verification */

        //Finds list of Hecke eigenvalues and stores them in tempvec.
        //Note that this will only include the eigenvalues with gcd(n,N)=1.
        for(j=0;j<i;j++){
            arb_set(arb_mat_entry(tempvecT,0,j),arb_mat_entry(M1,0,j));
        }
        arb_mat_mul(mat11,tempvecT,ci,prec);
        arb_mat_mul(tempvec,M1,ci,prec);
        arb_mat_scalar_div_arb(tempvec,tempvec,arb_mat_entry(mat11,0,0),prec);

        // Find the |lambda_j - \tilde{lambda_i}| >= delta_i = b1.
        if(n==0){
            arb_sub(b1,arb_mat_entry(R_vals,n,0),arb_mat_entry(R_vals,n,1),prec);
            arb_add(b2,arb_mat_entry(R_vals,n+1,0),arb_mat_entry(R_vals,n+1,1),prec);
            arb_sub(b1,b1,b2,prec);
            arb_abs(b1,b1);
        } else if(n==m-1){
            arb_sub(b1,arb_mat_entry(R_vals,n-1,0),arb_mat_entry(R_vals,n-1,1),prec);
            arb_add(b2,arb_mat_entry(R_vals,n,0),arb_mat_entry(R_vals,n,1),prec);
            arb_sub(b1,b1,b2,prec);
            arb_abs(b1,b1);
        } else {
            arb_sub(b1,arb_mat_entry(R_vals,n,0),arb_mat_entry(R_vals,n,1),prec);
            arb_add(b2,arb_mat_entry(R_vals,n+1,0),arb_mat_entry(R_vals,n+1,1),prec);
            arb_sub(b1,b1,b2,prec);
            arb_abs(b1,b1);

            arb_sub(b2,arb_mat_entry(R_vals,n-1,0),arb_mat_entry(R_vals,n-1,1),prec);
            arb_add(b4,arb_mat_entry(R_vals,n,0),arb_mat_entry(R_vals,n,1),prec);
            arb_sub(b2,b2,b4,prec);
            arb_abs(b2,b2);

            if(arb_lt(b2,b1)){
                arb_set(b1,b2);
            }
        }
        j = 0;
        for(g = 0;g<N;g++){
            if(gcd(g+1,level)==1){
                arb_set(arb_mat_entry(coeffs,g,0),arb_mat_entry(tempvec,j,0));
                arb_mul(b2,b3,arb_mat_entry(M1,j,j),prec);
                arb_sqrt(b2,b2,prec);
                arb_mul(b2,b2,arb_mat_entry(R_vals,n,1),prec);
                arb_div(arb_mat_entry(coeffs,g,1),b2,b1,prec);
                arb_add_error(arb_mat_entry(coeffs,g,0),arb_mat_entry(coeffs,g,1));
                j++;
            }
            
        }
        /* End of Hecke eigenvalue verification */

        /* Start of Hecke eigenvalue computation for (n,level) > 1 */
        g = (1 << kN); //Set g to 2^kN+1

        // Here we are trying to calcuate a_p for p|level. We know this is a_p=+-/sqrt(p). To begin we test a_p = +1/sqrt(p) with both w = -1 then w= +1. Then we
        // test a_p = -1/sqrt(p) with both w = -1 then w= +1.
        // If form is even.
	    fouriercheck = 0;
        if(parity==0){
            arb_one(b5); // b5 will be checking which coefficients to choose.
            // Setting up nu = iR
            arb_set(acb_imagref(z1),arb_mat_entry(R_vals,n,0)); arb_zero(acb_realref(z1));
            // Setup b1 = 1/sqrt(level)
            arb_rsqrt_ui(b4,level,prec);
            for(x=0;x<g;x++){
                omega = 1;
                for(y=0;y<kN;y++){
                    arb_rsqrt_ui(arb_mat_entry(coeffs,pN[y]-1,0),pN[y],prec);
                    if(((x >> y) & 1) == 0){
                        arb_neg(arb_mat_entry(coeffs,pN[y]-1,0),arb_mat_entry(coeffs,pN[y]-1,0));
                    } else {
                        omega *= -1;
                    }
                    for(j=2;j<=N/pN[y];j++){
                        arb_mul(arb_mat_entry(coeffs,pN[y]*j-1,0),arb_mat_entry(coeffs,pN[y]-1,0),arb_mat_entry(coeffs,j-1,0),prec);
                    }
                }
                if(omega == -1){
                    // First we try w = -1.
                    arb_zero(b1);
                    for(j=0;j<N;j++){
                        arb_const_pi(b2,prec);
                        arb_mul_ui(b2,b2,2*(j+1),prec);
                        arb_mul(b2,b2,b4,prec);
                        arb_sqrt(b3,b2,prec);
                        acb_set_arb(z2,b2);

                        acb_my_Kbessel(z3,z1,z2,prec);
                        arb_rsqrt_ui(b2,j+1,prec);
                        arb_mul(b2,b2,arb_mat_entry(coeffs,j,0),prec);
                        arb_mul(b2,b2,b3,prec);
                        arb_mul(b2,b2,acb_realref(z3),prec);
                        arb_add(b1,b1,b2,prec);
                    }
                    // Computing error of truncated sum
                    arb_const_pi(b2,prec);
                    arb_mul_2exp_si(b2,b2,1);
                    arb_mul(b2,b2,b4,prec);
                    arb_exp(b3,b2,prec);
                    arb_sub_ui(b3,b3,1,prec);

                    arb_mul_si(b2,b2,-N,prec);
                    arb_exp(b2,b2,prec);
                    arb_div(b2,b2,b3,prec);
                    arb_set_d(b3,1.758); // This is bound of Fourier coefficients
                    arb_mul(b2,b2,b3,prec);
                    arb_const_pi(b3,prec);
                    arb_mul_2exp_si(b3,b3,-1);
                    arb_sqrt(b3,b3,prec);
                    arb_mul(b2,b2,b3,prec);

                    //arb_add_error(b1,b2);
		            arb_abs(b1,b1);
                    if(arb_gt(b1,b2)){
                        continue;
                    } else {
                        if(fouriercheck == 0){
                            sign = x;
                            fouriercheck = 1;
                        }
                        else{
                            fouriercheck = 2;
                            break;
                        }
                    }
                } else if(omega == 1){
                    // Now try w=1.
                    arb_zero(b1);
                    for(j=0;j<N;j++){
                        arb_const_pi(b2,prec);
                        arb_sqrt_ui(b3,2,prec);
                        arb_mul(b2,b2,b3,prec);
                        arb_mul_ui(b2,b2,2*(j+1),prec);
                        arb_mul(b2,b2,b4,prec);
                        arb_sqrt(b3,b2,prec);
                        acb_set_arb(z2,b2);
                        
                        acb_my_Kbessel(z3,z1,z2,prec);
                        acb_mul_arb(z3,z3,b3,prec);

                        arb_mul_2exp_si(b2,b2,-1);
                        arb_sqrt(b3,b2,prec);
                        acb_set_arb(z2,b2);
                        
                        acb_my_Kbessel(z4,z1,z2,prec);
                        acb_mul_arb(z4,z4,b3,prec);
                        acb_sub(z3,z3,z4,prec);
                        arb_rsqrt_ui(b2,j+1,prec);
                        arb_mul(b2,b2,arb_mat_entry(coeffs,j,0),prec);
                        arb_mul(b2,b2,acb_realref(z3),prec);
                        arb_add(b1,b1,b2,prec);
                    }

                    // Computing error of truncated sum
                    arb_const_pi(b2,prec);
                    arb_sqrt_ui(b3,2,prec);
                    arb_mul(b2,b2,b3,prec);
                    arb_mul(b2,b2,b4,prec);
                    arb_exp(b3,b2,prec);
                    arb_sub_ui(b3,b3,1,prec);

                    arb_mul_si(b2,b2,-N,prec);
                    arb_exp(b2,b2,prec);
                    arb_div(b2,b2,b3,prec);
                    arb_set_d(b3,1.758); // This is bound of Fourier coefficients
                    arb_mul(b2,b2,b3,prec);
                    arb_mul_2exp_si(b2,b2,1);
                    arb_const_pi(b3,prec);
                    arb_mul_2exp_si(b3,b3,-1);
                    arb_sqrt(b3,b3,prec);
                    arb_mul(b2,b2,b3,prec);

                    //arb_add_error(b1,b2);
		            arb_abs(b1,b1);
                    if(arb_gt(b1,b2)){
                        continue;
                    } else {
                        if(fouriercheck == 0){
                            sign = x;
                            fouriercheck = 1;
                        }
                        else{
                            fouriercheck = 2;
                            break;
                        }
                    }
                }
            }
        } else if(parity==1){
            arb_one(b5); // b5 will be checking which coefficients to choose.
            // Setting up nu = iR
            arb_set(acb_imagref(z1),arb_mat_entry(R_vals,n,0)); arb_zero(acb_realref(z1));
            // Setup b1 = 1/sqrt(level)
            arb_rsqrt_ui(b4,level,prec);
            for(x=0;x<g;x++){
                omega = 1;
                for(y=0;y<kN;y++){
                    arb_rsqrt_ui(arb_mat_entry(coeffs,pN[y]-1,0),pN[y],prec);
                    if(((x >> y) & 1) == 0){
                        arb_neg(arb_mat_entry(coeffs,pN[y]-1,0),arb_mat_entry(coeffs,pN[y]-1,0));
                    } else {
                        omega *= -1;
                    }
                    for(j=2;j<=N/pN[y];j++){
                        arb_mul(arb_mat_entry(coeffs,pN[y]*j-1,0),arb_mat_entry(coeffs,pN[y]-1,0),arb_mat_entry(coeffs,j-1,0),prec);
                    }
                }
                if(omega == -1){
                    // Now try w=1.
                    arb_zero(b1);
                    for(j=0;j<N;j++){
                        arb_const_pi(b2,prec);
                        arb_sqrt_ui(b3,2,prec);
                        arb_mul(b2,b2,b3,prec);
                        arb_mul_ui(b2,b2,2*(j+1),prec);
                        arb_mul(b2,b2,b4,prec);
                        arb_sqrt(b3,b2,prec);
                        acb_set_arb(z2,b2);
                        
                        acb_my_Kbessel(z3,z1,z2,prec);
                        acb_mul_arb(z3,z3,b3,prec);

                        arb_mul_2exp_si(b2,b2,-1);
                        arb_sqrt(b3,b2,prec);
                        acb_set_arb(z2,b2);
                        
                        acb_my_Kbessel(z4,z1,z2,prec);
                        acb_mul_arb(z4,z4,b3,prec);
                        acb_mul_2exp_si(z4,z4,-1);
                        acb_sub(z3,z3,z4,prec);
                        arb_sqrt_ui(b2,j+1,prec);
                        arb_mul(b2,b2,arb_mat_entry(coeffs,j,0),prec);
                        arb_mul(b2,b2,acb_realref(z3),prec);
                        arb_add(b1,b1,b2,prec);
                    }

                    // Computing error of truncated sum
                    arb_const_pi(b2,prec);
                    arb_sqrt_ui(b3,2,prec);
                    arb_mul(b2,b2,b3,prec);
                    arb_mul(b2,b2,b4,prec);
                    arb_exp(b3,b2,prec);
                    arb_sub_ui(b3,b3,1,prec);
                    arb_mul(b3,b3,b3,prec);

                    arb_mul_si(b2,b2,-N,prec);
                    arb_exp(b2,b2,prec);
                    arb_div(b2,b2,b3,prec);
                    arb_set_d(b3,1.758); // This is bound of Fourier coefficients
                    arb_mul_ui(b3,b3,3,prec);
                    arb_mul_2exp_si(b3,b3,-1);
                    arb_mul(b2,b2,b3,prec);
                    arb_const_pi(b3,prec);
                    arb_mul_2exp_si(b3,b3,-1);
                    arb_sqrt(b3,b3,prec);
                    arb_mul(b2,b2,b3,prec);

                    // sqrt(2*pi^2) = pi * sqrt(2)
                    arb_const_pi(b3,prec);
                    arb_mul(b3,b3,b3,prec);
                    arb_mul_2exp_si(b3,b3,1);
                    arb_sqrt(b3,b3,prec);
                    arb_mul(b3,b3,b4,prec);
                    arb_exp(b3,b3,prec);
                    arb_mul_ui(b3,b3,N+1,prec);
                    arb_sub_ui(b3,b3,N,prec);
                    arb_mul(b2,b2,b3,prec);

                    //arb_add_error(b1,b2);
	  	            arb_abs(b1,b1);
                    if(arb_gt(b1,b2)){
                        continue;
                    } else {
                        if(fouriercheck == 0){
                            sign = x;
                            fouriercheck = 1;
                        }
                        else{
                            fouriercheck = 2;
                            break;
                        }
                    }
                } else if(omega == 1){
                    // First we try w = -1.
                    arb_zero(b1);
                    for(j=0;j<N;j++){
                        arb_const_pi(b2,prec);
                        arb_mul_ui(b2,b2,2*(j+1),prec);
                        arb_mul(b2,b2,b4,prec);
                        arb_sqrt(b3,b2,prec);
                        acb_set_arb(z2,b2);

                        acb_my_Kbessel(z3,z1,z2,prec);
                        arb_sqrt_ui(b2,j+1,prec);
                        arb_mul(b2,b2,arb_mat_entry(coeffs,j,0),prec);
                        arb_mul(b2,b2,b3,prec);
                        arb_mul(b2,b2,acb_realref(z3),prec);
                        arb_add(b1,b1,b2,prec);
                    }

                    // Computing error of truncated sum
                    arb_const_pi(b2,prec);
                    arb_mul_2exp_si(b2,b2,1);
                    arb_mul(b2,b2,b4,prec);
                    arb_exp(b3,b2,prec);
                    arb_sub_ui(b3,b3,1,prec);
                    arb_mul(b3,b3,b3,prec);

                    arb_mul_si(b2,b2,-N,prec);
                    arb_exp(b2,b2,prec);
                    arb_div(b2,b2,b3,prec);
                    arb_set_d(b3,1.758); // This is bound of Fourier coefficients
                    arb_mul(b2,b2,b3,prec);
                    arb_const_pi(b3,prec);
                    arb_mul_2exp_si(b3,b3,-1);
                    arb_sqrt(b3,b3,prec);
                    arb_mul(b2,b2,b3,prec);

                    arb_const_pi(b3,prec);
                    arb_mul_2exp_si(b3,b3,1);
                    arb_mul(b3,b3,b4,prec);
                    arb_exp(b3,b3,prec);
                    arb_mul_ui(b3,b3,N+1,prec);
                    arb_sub_ui(b3,b3,N,prec);
                    arb_mul(b2,b2,b3,prec);

                    //arb_add_error(b1,b2);
		            arb_abs(b1,b1);
                    if(arb_gt(b1,b2)){
                        continue;
                    } else {
                        if(fouriercheck == 0){                    
                            sign = x;
                            fouriercheck = 1;
                        }
                        else{
                            fouriercheck = 2;
                            break;
                        }
                    }
                }
            }
        }
    // Checking to see if our data was computed to enough precision so that we can actually apply this part of the verification. If not we just skip it.
	if(fouriercheck == 2 || fouriercheck == 0){
            omega = 0;
            for(y=0;y<kN;y++){
                arb_zero(arb_mat_entry(coeffs,pN[y]-1,0));
                for(j=2;j<=N/pN[y];j++){
                    arb_mul(arb_mat_entry(coeffs,pN[y]*j-1,0),arb_mat_entry(coeffs,pN[y]-1,0),arb_mat_entry(coeffs,j-1,0),prec);
                }
            }
        } else {
            omega = 1;
            for(y=0;y<kN;y++){
                arb_rsqrt_ui(arb_mat_entry(coeffs,pN[y]-1,0),pN[y],prec);
                if(((sign >> y) & 1) == 0){
                    arb_neg(arb_mat_entry(coeffs,pN[y]-1,0),arb_mat_entry(coeffs,pN[y]-1,0));
                } else {
                    omega *= -1;
                }
                for(j=2;j<=N/pN[y];j++){
                    arb_mul(arb_mat_entry(coeffs,pN[y]*j-1,0),arb_mat_entry(coeffs,pN[y]-1,0),arb_mat_entry(coeffs,j-1,0),prec);
                }
            }
        }

        /* End of Hecke eigenvalue computation for (n,level) > 1 */

        // Prints the data for each Maass form into its own file.
        sprintf(filename,"%s%.10lf.txt",direc,arf_get_d(arb_midref(arb_mat_entry(R_vals,n,0)),ARF_RND_NEAR));
        fp = fopen(filename,"w");
        arf_fprintd(fp,arb_midref(arb_mat_entry(R_vals,n,0)),30);
        fprintf(fp,",");
        arf_fprintd(fp,arb_midref(arb_mat_entry(R_vals,n,1)),30);
        fprintf(fp,"\n");
        fprintf(fp,"%d\n",level);
        fprintf(fp,"%d\n",parity);
        fprintf(fp,"%d\n",omega);
        for(j=0;j<N;j++){
            arf_fprintd(fp,arb_midref(arb_mat_entry(coeffs,j,0)),30);
            fprintf(fp,",");
            arf_fprintd(fp,arb_midref(arb_mat_entry(coeffs,j,1)),30);
            fprintf(fp,"\n");
        }
        fclose(fp);     
    }
    arb_mat_clear(M1);
    arb_mat_clear(M2);
    arb_mat_clear(M4);
    arb_mat_clear(ci);
    arb_mat_clear(D);
    arb_mat_clear(Drtinv);
    arb_mat_clear(P);
    arb_mat_clear(R);
    arb_mat_clear(Rinv);
    arb_mat_clear(RinvT);
    arb_mat_clear(tempvec);
    arb_mat_clear(tempvecT);
    arb_mat_clear(mat11);
    arb_mat_clear(tempmat);

    acb_mat_clear(CM1);
    acb_mat_clear(CP);
}

int main(int argc,char *argv[]) {
    long TFpow,level;
    if (argc < 2 || sscanf(argv[1],"%ld",&TFpow) != 1){
		return 0;
    }
    initfactor64("./factor64-master/factor.bin");
    int kprocs=0,k;
    TFdata tracevals;
    arb_mat_t oddtraces,eventraces;

    // Loop over all levels that we want to perform this procedure for.
    for(level=2;level<=MAXLEVEL;level++){
        if(n_is_squarefree(level)){
            if (!fork()) {
		        printf("Process started for level = %ld\n",level);
                arb_mat_init(oddtraces,MAXN,3);
                arb_mat_init(eventraces,MAXN,3);
                arb_mat_init(tracevals.posterms,MAXN,3);
                arb_mat_init(tracevals.negterms,MAXN,3);
                tracevals.N = level;
                readtfdata(tracevals);
                arb_mat_add(eventraces,tracevals.posterms,tracevals.negterms,prec);
                arb_mat_sub(oddtraces,tracevals.posterms,tracevals.negterms,prec);
                arb_mat_clear(tracevals.posterms);
                arb_mat_clear(tracevals.negterms);               
                arb_mat_scalar_mul_2exp_si(eventraces,eventraces,-1);
                arb_mat_scalar_mul_2exp_si(oddtraces,oddtraces,-1);
                printf("level = %ld, pid = %d\n",level,getpid());
                printf("Starting verification for level = %ld\n",level);
                verifarb(eventraces,MATDIM,level,0);
                arb_mat_clear(eventraces);
                printf("Even terms done for level = %ld\n",level);
                verifarb(oddtraces,MATDIM,level,1);
                arb_mat_clear(oddtraces);
                printf("Odd terms done terms done for level = %ld\n",level);
                printf("Computing forms for level = %ld done.\n",level);\
                _exit(0);
            } else {
                if (kprocs == NPROCS - 1){
                    wait(&k);
                    //CHECK STATUS HERE, CRASH WILL BE NONZERO
                } else{
                    kprocs++;
                }
            }
        }    
    }
    while(wait(&k) > 0);
    return 0;
}
