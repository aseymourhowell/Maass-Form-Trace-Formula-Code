
#include <acb_mat.h>
#include <stdio.h>
#include <sys/file.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <math.h>

#define NPROCS 2
#define NUMSAMPLEPOINTS 20
#define NUMTAYLORTERMS 200

static slong prec = 1000;
#include "testfunc.c"

void printarf(const arf_t x){
    static int init;
    static fmpz_t m,e;

    if(!init){
        fmpz_init(m); fmpz_init(e);
        init = 1;
    }

    arf_get_fmpz_2exp(m,e,x);
    fmpz_print(m); printf(" "); fmpz_print(e);
}

void printarb(const arb_t x){
    static int init;
    static arf_t a;

    if (!init) {
        arf_init(a);
        init = 1;
    }

    printarf(arb_midref(x));
    printf(" ");
    arf_set_mag(a,arb_radref(x));
    printarf(a);
}

//precompute weights
void calcweights(arb_mat_t w,int n){
    int k;
    arb_t x,a,c,s,h;
    arb_init(x); arb_init(a);
    arb_init(c); arb_init(s);
    arb_init(h);

 	arb_log_ui(h,5*n,prec);
	arb_div_ui(h,h,n,prec);

    for (k=-n;k<=n;k++) {
        // a = h*cosh(k*h)/cosh(sinh(k*h))^2, x = tanh(sinh(k*h))
        arb_mul_si(x,h,k,prec);
		arb_sinh_cosh(s,c,x,prec);
		arb_cosh(a,s,prec);
		arb_mul(a,a,a,prec);
		arb_div(a,c,a,prec);
		arb_mul(a,a,h,prec);
		arb_tanh(x,s,prec);

        arb_set(arb_mat_entry(w,0,k+n),x);
        arb_set(arb_mat_entry(w,1,k+n),a);
    }
    arb_clear(x); arb_clear(a);
    arb_clear(c); arb_clear(s);
    arb_clear(h);
}

//Precompute the g values.
void precompg(acb_mat_t vals, tfd g, const arb_mat_t w){
    int k,n,i;
    arb_t a,b,x,c,ab,zeroone;
    acb_t z,y,sh,ch;
    n = acb_mat_ncols(w);

    arb_init(a); arb_init(b);
    arb_init(x); arb_init(c);
    arb_init(ab); arb_init(zeroone);
    acb_init(y); acb_init(z);
    acb_init(sh); acb_init(ch);

    arb_one(zeroone);
    arb_union(zeroone,zeroone,x,53);

    for(i=0;i<g.d;i++){
        arb_set_ui(a,i);
        arb_mul(a,a,g.a,prec);
        arb_add(b,a,g.a,prec);
        arb_sub(c,b,a,prec);
        arb_add(ab,a,b,prec);
        for(k=0;k<n;k++){
		    arb_mul(x,arb_mat_entry(w,0,k),c,prec);
		    arb_add(x,x,ab,prec);
            arb_mul_2exp_si(acb_realref(z),x,-1);
		    arb_zero(acb_imagref(z));
            evaltfd(acb_mat_entry(vals,i,k),&g,i,z);
        }

        //STORE THESE VALUES IN A MATRIX.
        for(k=-n/2;k<n/2;k++){
            arb_add_si(x,zeroone,k,prec);
            arb_div_ui(acb_realref(z),x,n,prec);
            arb_zero(acb_imagref(z));
            acb_exp_pi_i(z,z,prec);
            acb_mul_arb(z,z,c,prec);
            arb_add(acb_realref(z),acb_realref(z),ab,prec);
            evaltfd(y,&g,i,z);
        }
    }
    arb_clear(a); arb_clear(b);
    arb_clear(x); arb_clear(c);
    arb_clear(ab); arb_clear(zeroone);
    acb_clear(y); acb_clear(z);
    acb_clear(sh); acb_clear(ch);
}

//Computes the integrand for the elliptic terms. Here D4n will be |D/4n|.
//This will be used for the split intervals.
void elliptic_integrand(acb_t res, void *param, int j, acb_srcptr x,arb_t D4n){
    static acb_t c,s;
    static int ell_init;
    if(!ell_init){
        acb_init(s);
        acb_init(c);
        ell_init = 1;
    }

    acb_mul_2exp_si(s,x,-1);
    acb_sinh_cosh(s,c,s,prec);
    acb_mul(s,s,s,prec);
    acb_add_arb(s,s,D4n,prec);
    acb_pow_ui(s,s,10,prec);

    evaltfd(res, param, j, x);
    acb_mul(res, res, c, prec);
    acb_div(res, res, s, prec);
}

//compute integral of d^t/dx^t f(x) between a and b. Here D4n = |D/4n| and t = Taylorpow.
void ellintegral(arb_t res,arb_mat_t w,acb_mat_t vals,void *p,int m,arb_t D4n,arb_t a,arb_t b,int Taylorpow,int split){
    int k,n;
    static int init = 0;
    static arb_t x,c,zeroone,h,r,ab;
    static acb_t z,sh,ch,y;

    if(!init){
        arb_init(x); arb_init(c);
        arb_init(zeroone); arb_init(h);
        arb_init(r); arb_init(ab);
        acb_init(z); acb_init(y);
        acb_init(ch); acb_init(sh);
        arb_one(zeroone);
		arb_union(zeroone,zeroone,x,53);
        init = 1;
    }

    n = arb_mat_ncols(w);
    arb_zero(res);
    arb_sub(c,b,a,prec);
    arb_add(ab,a,b,prec);
    for(k=0;k<n;k++){
        arb_mul(x,arb_mat_entry(w,0,k),c,prec);
        arb_add(x,x,ab,prec);
        arb_mul_2exp_si(acb_realref(z),x,-1);
        arb_zero(acb_imagref(z));

        if(split){
            evaltfd(y,p,m,z);
            acb_mul_2exp_si(z,z,-1);
            acb_sinh_cosh(sh,ch,z,prec);
            acb_mul(sh,sh,sh,prec);
            acb_add_arb(sh,sh,D4n,prec);
            acb_pow_ui(sh,sh,Taylorpow,prec);

            acb_mul(y,y,ch,prec);
            acb_div(y,y,sh,prec);
            arb_addmul(res,arb_mat_entry(w,1,k),acb_realref(y),prec);
        }else{
            //Actual elliptic integrand of g(u)cosh(u/2)/(sinh^2(u/2) + |D/4n|)
            acb_mul_2exp_si(z,z,-1);
            acb_sinh_cosh(sh,ch,z,prec);
            acb_mul(sh,sh,sh,prec);
            acb_add_arb(sh,sh,D4n,prec);
            acb_pow_ui(sh,sh,Taylorpow,prec);

            acb_mul(ch,acb_mat_entry(vals,m,k),ch,prec);
            acb_div(ch,ch,sh,prec);

            arb_addmul(res,arb_mat_entry(w,1,k),acb_realref(ch),prec);
        }
    }

    //Compute error term by sampling integrand on the circle
    //of radius b-a around (a+b)/2 using 2n angle intervals.
    arb_zero(r);
    arb_mul_2exp_si(ab,ab,-1); //ab=(a+b)/2
    n = n/2;
    for(k=-n;k<n;k++){
        arb_add_si(x,zeroone,k,prec);
        arb_div_ui(acb_realref(z),x,n,prec);
        arb_zero(acb_imagref(z));
        acb_exp_pi_i(z,z,prec);
        acb_mul_arb(z,z,c,prec);
        arb_add(acb_realref(z),acb_realref(z),ab,prec);

        evaltfd(y,p,m,z);
        acb_mul_2exp_si(z,z,-1);
        acb_sinh_cosh(sh,ch,z,prec);
        acb_mul(sh,sh,sh,prec);
        acb_add_arb(sh,sh,D4n,prec);
        acb_log(sh,sh,prec);
        acb_mul_si(sh,sh,-Taylorpow,prec);
        acb_exp(sh,sh,prec);

        acb_mul(y,acb_mat_entry(vals,m,abs(k)),ch,prec);
        acb_mul(y,y,sh,prec);


        acb_abs(h,y,prec);
        arb_union(r,r,h,prec);
    }
    arb_neg(h,r);
    arb_union(r,r,h,prec);

    // h = exp(4-5*n/log(5*n))
    arb_log_ui(h,5*n,prec);
    arb_ui_div(h,5*n,h,prec);
    arb_sub_ui(h,h,4,prec);
    arb_neg(h,h);
    arb_exp(h,h,prec);

    arb_mul(r,r,h,prec);
    arb_add(res,res,r,prec);
    arb_mul(res, res, c, prec);
    arb_mul_2exp_si(res,res,-1);
}

//Checks to see if singularity wil occur by seeing if the denominator of the integrand has a zero in the disk D(0,2).
//If it does, it will split the interval in 2 and recurisevly call itself until all singularities are missed or
//100 splits have been made.
void splitinterval(arb_mat_t intervals, arb_t D4n, arb_t a, arb_t b, int *i){
    arb_t temp;
    arb_init(temp);
    double d1,d2,u;

    if(*i>100){
        printf("Error: Too many integral splits.\n");
        exit(0);
    }
    u = arf_get_d(arb_midref(D4n),ARF_RND_NEAR);
    if(u < 1){
        u = (4*asin(sqrt(u)));
        u = u*u;
        d1 = arf_get_d(arb_midref(a),ARF_RND_NEAR);
        d2 = arf_get_d(arb_midref(b),ARF_RND_NEAR);

        u=(u+pow(d1+d2,2))/(pow(d2-d1,2));

        if(u<4.5){
            arb_add(temp,a,b,prec);
            arb_mul_2exp_si(temp,temp,-1);
            splitinterval(intervals,D4n,a,temp,i);
            *i=*i+1;
            splitinterval(intervals,D4n,temp,b,i);
            arb_clear(temp);
            return;
        }
    }
    arb_set(arb_mat_entry(intervals,0,*i),b);
    arb_clear(temp);
    return;
}

//Computes the elliptic integral for all 3 test functions in the interval [a,b], then stores the respectives results in res.
//The code will subdivide the interval so that the integrand is holomorphic on the disk D(0,2) and we can apply the rigorous error bound.
void compute_ell_integral(arb_mat_t res, acb_mat_t gvals0, acb_mat_t gvals2, acb_mat_t gvals4, arb_mat_t w, tfd g, tfd gprime2, tfd gprime4, int j, arb_t D4n, arb_t a, arb_t b,int Taylorpow){
    static arb_t inte0, inte2, inte4, temp, b1, b2;
    static int init;
    static arb_mat_t inte100;
    int i;
    int count=1;
    arb_mat_t intervals;
    if(!init){
        arb_init(inte0); arb_init(inte2); arb_init(inte4);
        arb_init(temp); arb_init(b1); arb_init(b2);
        arb_mat_init(inte100,1,100);
        init=1;
    }
    arb_mat_zero(inte100);
    arb_set(arb_mat_entry(inte100,0,0),a);
    splitinterval(inte100,D4n,a,b,&count);
    arb_mat_window_init(intervals,inte100,0,0,1,count+1);

    if(count==1){
        ellintegral(temp, w, gvals0, &g, j, D4n,a,b,Taylorpow,0);
        arb_add(arb_mat_entry(res,0,0), arb_mat_entry(res,0,0), temp, prec);
        ellintegral(temp, w, gvals2, &gprime2, j, D4n,a,b,Taylorpow,0);
        arb_add(arb_mat_entry(res,0,1), arb_mat_entry(res,0,1), temp, prec);
        ellintegral(temp, w, gvals4, &gprime4, j, D4n,a,b,Taylorpow,0);
        arb_add(arb_mat_entry(res,0,2), arb_mat_entry(res,0,2), temp, prec);
    }else if(count>1){
        for(i=0;i<count;i++){
            ellintegral(temp, w, gvals0, &g, j, D4n,arb_mat_entry(intervals,0,i), arb_mat_entry(intervals,0,i+1),Taylorpow,1);
            arb_add(arb_mat_entry(res,0,0), arb_mat_entry(res,0,0), temp, prec);
            ellintegral(temp, w, gvals2, &gprime2, j, D4n,arb_mat_entry(intervals,0,i), arb_mat_entry(intervals,0,i+1),Taylorpow,1);
            arb_add(arb_mat_entry(res,0,1), arb_mat_entry(res,0,1), temp, prec);
            ellintegral(temp, w, gvals4, &gprime4, j, D4n,arb_mat_entry(intervals,0,i), arb_mat_entry(intervals,0,i+1),Taylorpow,1);
            arb_add(arb_mat_entry(res,0,2), arb_mat_entry(res,0,2), temp, prec);
        }
    }
    arb_mat_clear(intervals);
}

void computesamplepoints(arb_mat_t samplepoints, int n){
    arb_t c;
    int j;
    arb_init(c);

    //c = 1+sqrt(2)
    arb_set_ui(c,2);
    arb_sqrt(c,c,prec);
    arb_add_ui(c,c,1,prec);

    for(j=0;j<n;j++){
        arb_pow_ui(arb_mat_entry(samplepoints,0,j),c,j,prec);
        arb_inv(arb_mat_entry(samplepoints,0,j),arb_mat_entry(samplepoints,0,j),prec);
    }
    arb_clear(c);
}

int main(int argc, char const *argv[])
{
    int n,k,kprocs,i,j,sampleindex,Taylorpow;
    kprocs = 0;
    k = 0;

    if (argc < 2 || sscanf(argv[1],"%d",&i) != 1){
		return 0;
    }

    tf h;
    tfd g,gprime,gprime2,gprime3,gprime4;

    hd(&h,i);
    shifttf(&h);

    arb_t a,D4n,sum0,sum2,sum4,b1,b2,temp;
    arb_mat_t res,sample_points,w;
    arb_init(a);
    arb_init(D4n);
    arb_init(sum0); arb_init(temp);
    arb_init(sum2); arb_init(sum4);
    arb_init(b1); arb_init(b2);
    arb_mat_init(res,1,3);

    // Generate test function.
    arb_set_d(a,5.51341248666333);
    arb_div_ui(a,a,i,prec);
    dilateg(&g,&h,a);
    derivgd(&gprime,&g);
    derivgd(&gprime2,&gprime);
    derivgd(&gprime3,&gprime2);
    derivgd(&gprime4,&gprime3);

    //Precompute weights
    for (n=1;5*n/log(5*n)-4<prec*M_LN2;n++);
    //Init matrix for weights. First row is x_k values and second row is a_k values.
    arb_mat_init(w,2,2*n+1);
    calcweights(w,n);

    //Precompute g values
    acb_mat_t gvals0,gvals2,gvals4;
    acb_mat_init(gvals0,g.d,2*n+1);
    acb_mat_init(gvals2,gprime2.d,2*n+1);
    acb_mat_init(gvals4,gprime4.d,2*n+1);
    precompg(gvals0,g,w);
    precompg(gvals2,gprime2,w);
    precompg(gvals4,gprime4,w);

    // Compute the sample points for the Taylor series.
    arb_mat_init(sample_points,1,NUMSAMPLEPOINTS);
    computesamplepoints(sample_points,NUMSAMPLEPOINTS);

    //Computes the elliptic integrals.
    for(sampleindex=0;sampleindex<NUMSAMPLEPOINTS;sampleindex++){
        for(Taylorpow=0;Taylorpow<NUMTAYLORTERMS;Taylorpow++){
            if (!fork()) {
                arb_mat_zero(res);
                for (j = 0; j < g.d; j++){
                    arb_set_ui(b1, j);
                    arb_mul(b1, b1, g.a, prec);
                    arb_add(b2, b1, g.a, prec);
                    compute_ell_integral(res, gvals0, gvals2, gvals4, w, g, gprime2, gprime4, j, arb_mat_entry(sample_points,0,sampleindex), b1, b2,Taylorpow+1);
                }
                //arb_pow_ui(b1,arb_mat_entry(sample_points,0,sampleindex),Taylorpow,prec);
                //arb_mat_scalar_mul_arb(res,res,b1,prec);
                //This is multiplying the integral by (-1)^k.
                if(Taylorpow%2 != 0){
                    arb_mat_neg(res,res);
                }

                flock(1,LOCK_EX);
                //print lines of m,k,elliptic term.
                printf("%d,%d,",sampleindex,Taylorpow);
                arb_dump_file(stdout,arb_mat_entry(res,0,0));
                printf(",");
                arb_dump_file(stdout,arb_mat_entry(res,0,1));
                printf(",");
                arb_dump_file(stdout,arb_mat_entry(res,0,2));
                printf("\n");
                fflush(stdout);
                flock(1,LOCK_UN);
                _exit(0);
            } else {
                if (kprocs == NPROCS-1){
                    wait(&k);
                } else{
                    kprocs++;
                }
            }
        }
    }

    while (wait(&k) > 0);
    return 0;
}
