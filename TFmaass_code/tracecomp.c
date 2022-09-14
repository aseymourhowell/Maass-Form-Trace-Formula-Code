#include <arb_mat.h>
#include <inttypes.h>
#include <sys/file.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdio.h>
#include <dirent.h>

#define NPROCS 2
#define MAXLEVEL 105
#define MAXN 4000000
#define NUMSAMPLEPOINTS 20
#define NUMTAYLORTERMS 200

static slong prec = 300;
#include "testfunc.c"

int initfactor64(const char *);
int factor64(uint64_t *,int *,uint64_t);
int isprime64(uint64_t);

// Code to compute integral with rigorous error bound from Molin-2016.
typedef void (*function)(acb_t,void *,int,acb_srcptr);
void arb_integral(arb_t res,function f,void *p,int m,arb_t a,arb_t b) {
	static int init;
	static arb_t h,x,t,s,c,zeroone;
	static acb_t z,y;
	int k,n;

	if (!init) {
		arb_init(h); arb_init(x);
		arb_init(t); arb_init(zeroone);
		arb_init(s); arb_init(c);
		acb_init(z); acb_init(y);
		arb_one(zeroone);
		arb_union(zeroone,zeroone,x,53);
		init = 1;
	}

	// choose the number of steps so that the error is roughly 2^-prec
	for (n=1;5*n/log(5*n)-4<prec*M_LN2;n++);
	arb_log_ui(h,5*n,prec);
	arb_div_ui(h,h,n,prec);

	arb_zero(res);
	for (k=-n;k<=n;k++) {
		// t = h*cosh(k*h)/cosh(sinh(k*h))^2, x = tanh(sinh(k*h))
		arb_mul_si(x,h,k,prec);
		arb_sinh_cosh(s,c,x,prec);
		arb_cosh(t,s,prec);
		arb_mul(t,t,t,prec);
		arb_div(t,c,t,prec);
		arb_mul(t,t,h,prec);
		arb_tanh(x,s,prec);

		arb_sub(c,b,a,prec);
		arb_mul(x,x,c,prec);
		arb_add(x,x,a,prec);
		arb_add(x,x,b,prec);
		arb_mul_2exp_si(acb_realref(z),x,-1);
		arb_zero(acb_imagref(z));

		f(y,p,m,z);
		arb_addmul(res,t,acb_realref(y),prec);
	}

    //Compute error term by sampling integrand on the circle
    //of radius b-a around (a+b)/2 using 2n angle intervals.
    arb_zero(t);
    arb_add(s,a,b,prec); arb_mul_2exp_si(s,s,-1); // s = (a+b)/2
    for(k=-n;k<n;k++){
        arb_add_si(x,zeroone,k,prec);
		arb_div_ui(acb_realref(z),x,n,prec);
		arb_zero(acb_imagref(z));
		acb_exp_pi_i(z,z,prec);
        acb_mul_arb(z,z,c,prec);
		arb_add(acb_realref(z),acb_realref(z),s,prec);

		f(y,p,m,z);
		acb_abs(x,y,prec);
		arb_union(t,t,x,prec);
    }
    arb_neg(x,t);
	arb_union(t,t,x,prec);

    // h = exp(4-5*n/log(5*n))
	arb_ui_div(h,5,h,prec);
	arb_sub_ui(h,h,4,prec);
	arb_neg(h,h);
	arb_exp(h,h,prec);

    arb_mul(t,t,h,prec);
    arb_add(res,res,t,prec);
	arb_mul(res, res, c, prec);
	arb_mul_2exp_si(res,res,-1);
}

// Code to compute the integrand of g'(x)/sinh(x/2). We treat the points near x=0 seperately
// to ensure the error bounds don't blow up.
void identityintegral(acb_t res, void *param, int j, acb_srcptr x){
    static acb_t s;
    static int init;
    if(!init){
        acb_init(s);
        init = 1;
    }

    if(j==0){
        acb_mul_2exp_si(s,x,-1);
        acb_mul_onei(s,s);
        acb_sinc(s,s,prec);
        acb_mul_2exp_si(s,s,-1);

        evaltfdx(res,param,j,x);
        acb_div(res, res, s, prec);
    }else{
        acb_mul_2exp_si(s,x,-1);
        acb_sinh(s,s,prec);
        evaltfdx(res,param,j,x);
        acb_div(res, res, s, prec);
    }
}

int main(int argc,char *argv[]) {
    double y, scale;
    initfactor64("./factor64-master/factor.bin");
    slong D,t,hecken,d,a,l,nmin,nmax,tstart,sigma1;
    ulong count = 0;
    long i,j,Nprod = 1;
    int kprocs=0,k,N,levelcount=0,proc,*muN;
    char line[1000],Lval[105],inteval0[330],inteval2[330],inteval4[330];
    FILE *fp, *tmp;

    // Variables for parabolic terms
    uint64_t pn[15];
    int en[15],kn;
    long *divisors = (long *)malloc(100 * sizeof(long));

    if (argc < 2 || sscanf(argv[1],"%ld",&i) != 1){
    	return 0;
    }

    arb_t g_0,g2_0,g4_0,inte0, inte2, inte4,b1,b2,b3,D4n,lvalarb,temp,idenint0,idenint2,idenint4;
    acb_t ctemp,gval,gval0,gval2,gval4,hi2;
    arb_init(g_0); arb_init(g2_0); arb_init(g4_0);
    arb_init(inte0); arb_init(inte2); arb_init(inte4);
    arb_init(b1); arb_init(b2); arb_init(b3);
    arb_init(D4n); arb_init(lvalarb); arb_init(temp);
    arb_init(idenint0); arb_init(idenint2); arb_init(idenint4);
    acb_init(ctemp); acb_init(gval); acb_init(hi2);
    acb_init(gval0); acb_init(gval2); acb_init(gval4);

    tf h;
    tfd g,gprime,gprime2,gprime3,gprime4;
    hd(&h,i);
    shifttf(&h);

    // Generate test function.
    arb_t dilfac;
    arb_init(dilfac);
    arb_set_d(dilfac,5.51341248666333);
    arb_div_ui(dilfac,dilfac,i,prec);

    dilateg(&g,&h,dilfac);
    derivgd(&gprime,&g);
    derivgd(&gprime2,&gprime);
    derivgd(&gprime3,&gprime2);
    derivgd(&gprime4,&gprime3);
    scale = arf_get_d(arb_midref(g.a),ARF_RND_NEAR);

    long *levels = (long *)malloc(MAXLEVEL * sizeof(long)); //This is too many terms
    for(i = 2; i <= MAXLEVEL; i++){
        if(n_is_squarefree(i) == 1){
            levels[levelcount] = i;
            levelcount++;
        }
    }
    levels = (long *)realloc(levels, levelcount * sizeof(long)); // Realloc to correct number of terms
    muN = (int *)malloc(levelcount * sizeof(int));

    uint64_t p[levelcount][15];
    int e[levelcount][15],Npcount[levelcount];

    for(i=0;i<levelcount;i++){
        Npcount[i] = factor64(p[i],e[i],levels[i]);
        muN[i] = n_moebius_mu(levels[i]);
    }

    printf("%d,%ld\n",levelcount,levels[levelcount-1]);


    /* Set up for the elliptic terms */
    //Compute g(0), g''(0) and g''''(0) which will be used when computing the error bounds.
    //These will also be used for the parabolic terms.
    acb_zero(ctemp);
    evaltfd(gval,&g,0,ctemp);
    acb_get_real(g_0,gval);
    evaltfd(gval,&gprime2,0,ctemp);
    acb_get_real(g2_0,gval);
    evaltfd(gval,&gprime4,0,ctemp);
    acb_get_real(g4_0,gval);

    //b1 = 1+sqrt(2)
    arb_set_ui(b1,2);
    arb_sqrt(b1,b1,prec);
    arb_add_ui(b1,b1,1,prec);

    // Parameters that can change
    int samplepointindex;
    arb_mat_t samplepoints;
    arb_mat_init(samplepoints,1,NUMSAMPLEPOINTS);

    for(j=0;j<NUMSAMPLEPOINTS;j++){
        arb_pow_ui(arb_mat_entry(samplepoints,0,j),b1,j,prec);
        arb_inv(arb_mat_entry(samplepoints,0,j),arb_mat_entry(samplepoints,0,j),prec);
    }
    
    //Reading and organising the integrals precomputed for the
    //Taylor series for elliptic terms.
    fp = fopen("1000k200j20d13.txt", "r");
    arb_mat_t intevals0, intevals2, intevals4;
    arb_mat_init(intevals0,NUMSAMPLEPOINTS,NUMTAYLORTERMS);
    arb_mat_init(intevals2,NUMSAMPLEPOINTS,NUMTAYLORTERMS);
    arb_mat_init(intevals4,NUMSAMPLEPOINTS,NUMTAYLORTERMS);
    while(fgets(line,sizeof(line),fp) != NULL){
        sscanf(line,"%ld,%ld,%330[^,],%330[^,],%330[^\n]",&i,&j,inteval0,inteval2,inteval4);
        arb_load_str(b1,inteval0);
        arb_load_str(b2,inteval2);
        arb_load_str(b3,inteval4);

        arb_set(arb_mat_entry(intevals0,i,j),b1);
        arb_set(arb_mat_entry(intevals2,i,j),b2);
        arb_set(arb_mat_entry(intevals4,i,j),b3);
    }
    fclose(fp);

    arb_poly_t *T0 = (arb_poly_t *)malloc(NUMSAMPLEPOINTS * sizeof(arb_poly_t));
    arb_poly_t *T2 = (arb_poly_t *)malloc(NUMSAMPLEPOINTS * sizeof(arb_poly_t));
    arb_poly_t *T4 = (arb_poly_t *)malloc(NUMSAMPLEPOINTS * sizeof(arb_poly_t));

    for(i = 0; i < NUMSAMPLEPOINTS; i++){
        arb_poly_init(T0[i]);
        arb_poly_init(T2[i]);
        arb_poly_init(T4[i]);
        arb_poly_fit_length(T0[i],NUMTAYLORTERMS);
        _arb_poly_set_length(T0[i],NUMTAYLORTERMS);
        arb_poly_fit_length(T2[i],NUMTAYLORTERMS);
        _arb_poly_set_length(T2[i],NUMTAYLORTERMS);
        arb_poly_fit_length(T4[i],NUMTAYLORTERMS);
        _arb_poly_set_length(T4[i],NUMTAYLORTERMS);

        for(j=0; j<NUMTAYLORTERMS; j++){
            arb_poly_set_coeff_arb(T0[i],j,arb_mat_entry(intevals0,i,j));
            arb_poly_set_coeff_arb(T2[i],j,arb_mat_entry(intevals2,i,j));
            arb_poly_set_coeff_arb(T4[i],j,arb_mat_entry(intevals4,i,j));
        }
    }
    /* End of setup for elliptic terms*/
    
    /* Set up for the identity terms */
    for(i=0;i<g.d;i++){
        arb_set_ui(b1,i);
        arb_mul(b1,b1,g.a,prec);
        arb_add(b2,b1,g.a,prec);
        arb_integral(temp,identityintegral,&g,i,b1,b2);
        arb_add(idenint0, idenint0, temp, prec);
        arb_integral(temp,identityintegral,&gprime2,i,b1,b2);
        arb_add(idenint2, idenint2, temp, prec);
        arb_integral(temp,identityintegral,&gprime4,i,b1,b2);
        arb_add(idenint4, idenint4, temp, prec);
    }

    arb_div_si(idenint0,idenint0,-6,prec);
    arb_div_si(idenint2,idenint2,-6,prec);
    arb_div_si(idenint4,idenint4,-6,prec);

    // Cpmpute h(i/2)
    acb_onei(ctemp);
    acb_mul_2exp_si(ctemp, ctemp, -1);
    acb_mul_arb(ctemp, ctemp, g.a,prec);
    h_k(hi2, ctemp, g.d);

    /* End of setup for identity terms */

    for(proc = 0; proc < NPROCS; proc++){
        if(!fork()){
            // Assign each process a list of Hecke operators which they will compute the trace formula for.
            nmin = (proc)*(2*MAXN)/NPROCS - MAXN;
            nmax = (proc+1)*(2*MAXN)/NPROCS - MAXN;
            if(proc != NPROCS-1){
                nmax -= 1;
            }

            arb_mat_t *hypervalues = (arb_mat_t *)malloc(levelcount * sizeof(arb_mat_t));
            for(j = 0; j < levelcount; j++){
                arb_mat_init(hypervalues[j],labs(nmax-nmin)+1,3);
            }
            printf("proc = %d, nmax = %ld, nmin = %ld, diff = %ld\n",proc,nmax, nmin, labs(nmax-nmin)+1);
            
            /* Start of hyperbolic terms */
            count = 0;
            // Opens the file containing the L(1,\psi_D) data for 0 < D < 1000000000.
            fp = fopen("0_to_1000000000.out", "r");

            while(fgets(line,sizeof(line),fp) != NULL){
                sscanf(line,"%li %li %li %10000[^\n]",&D,&d,&l,Lval);
                arb_set_str(lvalarb,Lval,300);

                if(D + 4*nmin < 0){
                    tstart = 0;
                } else {
                    tstart = floor(sqrt(D + 4*nmin));
                }
                
                for(t = tstart; t*t <= D + 4*nmax; t++){
                    if((t*t - D) % 4 == 0){
                        
                        hecken = (t*t - D)/4;
                        y = 2*log((t + sqrt(D))/(2*sqrt(labs(hecken))))/scale;
                        j = floor(y);
                        if(j < g.d && (hecken <= nmax && hecken >= nmin)){
                            //b1=log(((|t|+sqrt(D))^2)/4|n|)
                            arb_set_ui(b1, D);
                            arb_sqrt(b1, b1, prec);
                            arb_add_ui(b1, b1, t, prec);
                            arb_mul(b1, b1, b1, prec);
                            arb_div_si(b1, b1, 4 * labs(hecken), prec);
                            arb_log(b1, b1, prec);
                            acb_set_arb(ctemp, b1);

                            // Computing the test functions g(b1),g''(b1),g''''(b1)
                            evaltfd(gval0, &g, j, ctemp);
                            evaltfd(gval2, &gprime2, j, ctemp);
                            evaltfd(gval4, &gprime4, j, ctemp);

                            //Loop over the level
                            for(N=0;N<levelcount;N++){
                                if(gcd(levels[N],labs(hecken)) == 1){
                                    arb_set(temp,lvalarb);
                                    Nprod = 1;
                                    // Calcuate the prod over the primes that divie N
                                    for(a=0;a<Npcount[N];a++){
                                        if(l % p[N][a] == 0){
                                            Nprod *= (kronecker(d,p[N][a])-1);
                                            arb_div_ui(temp,temp,1 + (p[N][a]-kronecker(d,p[N][a]))*((pow(p[N][a],ordp(l,p[N][a])) - 1)/(p[N][a]-1)),prec);
                                        } else {
                                            Nprod *= (kronecker(d,p[N][a])-1);
                                        }
                                    }

                                    if(Nprod == 0){
                                        continue;
                                    }
                                    arb_mul_si(temp,temp,Nprod,prec);
                                    count++;

                                    // Add terms to the trace formula values.
                                    if(hecken > 0){
                                        arb_mul(b1, temp, acb_realref(gval0), prec);
                                        arb_mul_2exp_si(b1, b1, 1);
                                        arb_add(arb_mat_entry(hypervalues[N], hecken - labs(nmin), 0), arb_mat_entry(hypervalues[N], hecken - labs(nmin), 0), b1, prec);

                                        arb_mul(b1, temp, acb_realref(gval2), prec);
                                        arb_mul_2exp_si(b1, b1, 1);
                                        arb_add(arb_mat_entry(hypervalues[N], hecken - labs(nmin), 1), arb_mat_entry(hypervalues[N], hecken - labs(nmin), 1), b1, prec);

                                        arb_mul(b1, temp, acb_realref(gval4), prec);
                                        arb_mul_2exp_si(b1, b1, 1);
                                        arb_add(arb_mat_entry(hypervalues[N], hecken - labs(nmin), 2), arb_mat_entry(hypervalues[N], hecken - labs(nmin), 2), b1, prec);
                                    } else {
                                        arb_mul(b1, temp, acb_realref(gval0), prec);
                                        if (t != 0){
                                            arb_mul_2exp_si(b1, b1, 1);
                                        }
                                        arb_add(arb_mat_entry(hypervalues[N], labs(hecken) - labs(nmax), 0), arb_mat_entry(hypervalues[N], labs(hecken) - labs(nmax), 0), b1, prec);

                                        arb_mul(b1, temp, acb_realref(gval2), prec);
                                        if (t != 0){
                                            arb_mul_2exp_si(b1, b1, 1);
                                        }
                                        arb_add(arb_mat_entry(hypervalues[N], labs(hecken) - labs(nmax), 1), arb_mat_entry(hypervalues[N], labs(hecken) - labs(nmax), 1), b1, prec);

                                        arb_mul(b1, temp, acb_realref(gval4), prec);
                                        if (t != 0){
                                            arb_mul_2exp_si(b1, b1, 1);
                                        }
                                        arb_add(arb_mat_entry(hypervalues[N], labs(hecken) - labs(nmax), 2), arb_mat_entry(hypervalues[N], labs(hecken) - labs(nmax), 2), b1, prec);
                                    }

                                } else {
                                    continue;
                                }   
                            }         
                        } else if(j >= g.d && labs(hecken) > MAXN){
                            goto end_loop;
                        } else {
                            continue;
                        }
                    }
                }
            }
            end_loop:
            printf("count = %ld, proc = %d\n",count,proc);
            fclose(fp);
            /* Hyperbolic terms done */

            /* Start of elliptic terms */
            // Elliptic terms only appear for positive hecken.
            if(nmin >= 0){
                // Opens the file containing the L(1,\psi_D) data for -16000000 < D < 0.
                fp = fopen("-16000000_to_0.out","r");
                while(fgets(line,sizeof(line),fp) != NULL){
                    sscanf(line,"%li %li %li %10000[^\n]",&D,&d,&l,Lval);
                    arb_set_str(lvalarb,Lval,300);

                    if(D + 4*nmin < 0){
                        tstart = 0;
                    } else {
                        tstart = floor(sqrt(D + 4*nmin));
                    }

                    for(t = tstart; t*t <= D + 4*nmax; t++){
                        if((t*t - D) % 4 == 0){
                            hecken = (t*t - D)/4;
                            if((hecken <= nmax && hecken >= nmin)){
                                arb_set_si(D4n,D);
                                arb_div_ui(D4n, D4n, 4*hecken, prec);
                                arb_abs(D4n, D4n);

                                // Calculate which sample point to use for Taylor approximation.
                                samplepointindex = ceil(-log(arf_get_d(arb_midref(D4n),ARF_RND_NEAR)*(2+sqrt(2))/2)/log(1+sqrt(2)));

                                // Compute |D4n - x_i|
                                arb_sub(b1,D4n,arb_mat_entry(samplepoints,0,samplepointindex),prec);
                                arb_poly_evaluate(inte0,T0[samplepointindex],b1,prec);
                                arb_poly_evaluate(inte2,T2[samplepointindex],b1,prec);
                                arb_poly_evaluate(inte4,T4[samplepointindex],b1,prec);

                                // Computing the Taylor error.
                                if(arb_lt(arb_mat_entry(samplepoints,0,samplepointindex),D4n)){
                                    //Computes |1 - x/x_0|^(K+1)
                                    arb_div(b1,D4n,arb_mat_entry(samplepoints,0,samplepointindex),prec);
                                    arb_sub_ui(b1,b1,1,prec);
                                    arb_abs(b1,b1);
                                    arb_pow_ui(b1,b1,NUMTAYLORTERMS,prec);

                                    arb_sqrt(b2,arb_mat_entry(samplepoints,0,samplepointindex),prec);
                                    arb_div(b1,b1,b2,prec);

                                    arb_mul(b2,g_0,b1,prec);
                                    arb_add_error(inte0,b2);
                                    arb_mul(b2,g2_0,b1,prec);
                                    arb_add_error(inte2,b2);
                                    arb_mul(b2,g4_0,b1,prec);
                                    arb_add_error(inte4,b2);
                                }else{
                                    //Computes |1 - x_0/x|^(K+1)
                                    arb_div(b1,arb_mat_entry(samplepoints,0,samplepointindex),D4n,prec);
                                    arb_sub_ui(b1,b1,1,prec);
                                    arb_abs(b1,b1);
                                    arb_pow_ui(b1,b1,NUMTAYLORTERMS,prec);

                                    arb_sqrt(b2,D4n,prec);
                                    arb_div(b1,b1,b2,prec);

                                    arb_mul(b2,g_0,b1,prec);
                                    arb_add_error(inte0,b2);
                                    arb_mul(b2,g2_0,b1,prec);
                                    arb_add_error(inte2,b2);
                                    arb_mul(b2,g4_0,b1,prec);
                                    arb_add_error(inte4,b2);
                                }

                                arb_const_pi(b1, prec);
                                arb_sqrt(b2, D4n, prec);
                                arb_div(b2, b2, b1, prec);
                                arb_mul(inte0, inte0, b2, prec);
                                arb_mul(inte2, inte2, b2, prec);
                                arb_mul(inte4, inte4, b2, prec);

                                for(i=0;i<levelcount;i++){
                                    if(gcd(levels[i],labs(hecken)) == 1){
                                        arb_set(b2,lvalarb);
                                        Nprod = 1;
                                        // Calcuate the prod over the primes that divide N
                                        for(j=0;j<Npcount[i];j++){
                                            if(l % p[i][j] == 0){
                                                Nprod *= (kronecker(d,p[i][j])-1);
                                                arb_div_ui(b2,b2,1 + (p[i][j]-kronecker(d,p[i][j]))*((pow(p[i][j],ordp(l,p[i][j])) - 1)/(p[i][j]-1)),prec);
                                            } else {
                                                Nprod *= (kronecker(d,p[i][j])-1);
                                            }
                                        }

                                        if(Nprod == 0){
                                            continue;
                                        }
                                        arb_mul_si(b2,b2,Nprod,prec);
                                        // Add elliptic terms to the trace formula values.
                                        if(hecken > 0){
                                            arb_mul(b1,b2,inte0,prec);
                                            if (t != 0){
                                                arb_mul_2exp_si(b1, b1, 1);
                                            }
                                            arb_add(arb_mat_entry(hypervalues[i],hecken-nmin,0),arb_mat_entry(hypervalues[i],hecken-nmin,0),b1,prec);

                                            arb_mul(b1,b2,inte2,prec);
                                            if (t != 0){
                                                arb_mul_2exp_si(b1, b1, 1);
                                            }
                                            arb_add(arb_mat_entry(hypervalues[i],hecken-nmin,1),arb_mat_entry(hypervalues[i],hecken-nmin,1),b1,prec);

                                            arb_mul(b1,b2,inte4,prec);
                                            if (t != 0){
                                                arb_mul_2exp_si(b1, b1, 1);
                                            }
                                            arb_add(arb_mat_entry(hypervalues[i],hecken-nmin,2),arb_mat_entry(hypervalues[i],hecken-nmin,2),b1,prec);
                                        }
                                    }
                                }
                            } else {
                                continue;
                            }
                        }
                    }
                }
                fclose(fp);
            }
            /* End of elliptic terms */

            /* Start of parabolic and identity terms */
            for(hecken=nmin;hecken<=nmax;hecken++){
                if(hecken == 0){
                    continue;
                }
                kn = factor64(pn,en,labs(hecken));
                // Count number of divisors.
                N = 1;
                for(i=0;i<kn;i++){
                    N *= (en[i] + 1);
                }

                divisors = (long *)realloc(divisors,N * sizeof(long));
                N = compute_divisors(hecken,pn,en,kn,divisors);
                sigma1 = 0;

                for(i=0;i<N;i++){
                    a = divisors[i];
                    d = labs(hecken)/a;

                    // Compute sigma1(|n|) by adding all the divisors
                    sigma1 += a;

                    if(a != d){
                        //setting up log(a/d), b1 = log(a/d)
                        arb_set_ui(b1,a);
                        arb_div_si(b1,b1,d,prec);
                        arb_log(b1,b1,prec);

                        //computing g(log|a/d|), g''(log|a/d|) and g''''(log|a/d|)
                        j = floor((log(a) - log(d))/scale);
                        if(j >= g.d || j < -g.d){
                            acb_zero(gval0);
                            acb_zero(gval2);
                            acb_zero(gval4);
                        } else if(j >= 0){
                            acb_set_arb(ctemp,b1);
                            evaltfd(gval0, &g, j, ctemp);
                            evaltfd(gval2, &gprime2, j, ctemp);
                            evaltfd(gval4, &gprime4, j, ctemp);
                        } else {
                            acb_set_arb(ctemp,b1);
                            acb_neg(ctemp, ctemp);
                            evaltfd(gval0, &g, abs(j)-1, ctemp);
                            evaltfd(gval2, &gprime2, abs(j)-1, ctemp);
                            evaltfd(gval4, &gprime4, abs(j)-1, ctemp);
                        }

                        for(l=0;l<levelcount;l++){
                            // Parabolic terms only appear if level is prime
                            if(n_is_prime(levels[l])){
                                if(gcd(levels[l],labs(hecken)) == 1){
                                    arb_zero(inte0);
                                    arb_zero(inte2);
                                    arb_zero(inte4);
                                    // Setting temp to Lambda(N) = log(N)
                                    arb_log_ui(temp,levels[l],prec);
                                    D = floor((g.d*scale + log(a) - log(d))/(2*log(levels[l]))) + 1;
                                    for(t = 0; t < D; t++){
                                        j = floor((2*t*log(levels[l]) - log(a) + log(d))/scale);
                                        if(j >= g.d || j < -g.d){
                                            continue;
                                        }else if(j>=0){
                                            //Computing 2*t*log(N) - log(a/d)
                                            arb_mul_ui(b2,temp,2*t,prec);
                                            arb_sub(b2,b2,b1,prec);
                                            acb_set_arb(ctemp,b2);

                                            //Computing N^{t}
                                            arb_ui_pow_ui(b2,levels[l],t,prec);

                                            //Computing N^{-t} * g(2*t*log(N) - log(a/d)) and adding it to g_const sum.
                                            evaltfd(gval,&g,j,ctemp);
                                            arb_div(b3,acb_realref(gval),b2,prec);
                                            arb_add(inte0,inte0,b3,prec);

                                            evaltfd(gval,&gprime2,j,ctemp);
                                            arb_div(b3,acb_realref(gval),b2,prec);
                                            arb_add(inte2,inte2,b3,prec);

                                            evaltfd(gval,&gprime4,j,ctemp);
                                            arb_div(b3,acb_realref(gval),b2,prec);
                                            arb_add(inte4,inte4,b3,prec);
                                        } else {
                                            j = floor((-2*t*log(levels[l]) + log(a) - log(d))/scale);
                                            //Computing log(a/d) - 2*t*log(N)
                                            arb_mul_ui(b2,temp,2*t,prec);
                                            arb_sub(b2,b1,b2,prec);
                                            acb_set_arb(ctemp, b2);

                                            //Computing N^{t}
                                            arb_ui_pow_ui(b2,levels[l],t,prec);

                                            //Computing N^{-t} * g(log(a/d) - 2*t*log(N)) and adding it to g_const sum.
                                            evaltfd(gval,&g,j,ctemp);
                                            arb_div(b3,acb_realref(gval),b2,prec);
                                            arb_add(inte0,inte0,b3,prec);

                                            evaltfd(gval,&gprime2,j,ctemp);
                                            arb_div(b3,acb_realref(gval),b2,prec);
                                            arb_add(inte2,inte2,b3,prec);

                                            evaltfd(gval,&gprime4,j,ctemp);
                                            arb_div(b3,acb_realref(gval),b2,prec);
                                            arb_add(inte4,inte4,b3,prec);
                                        }
                                    }
                                    //Scaling the sum by 2*log(N)
                                    arb_mul(inte0,inte0,temp,prec);
                                    arb_mul_2exp_si(inte0,inte0,1);
                                    arb_mul(inte2,inte2,temp,prec);
                                    arb_mul_2exp_si(inte2,inte2,1);
                                    arb_mul(inte4,inte4,temp,prec);
                                    arb_mul_2exp_si(inte4,inte4,1);

                                    if(hecken > 0){
                                        //Computing gcd(N^{infinity},|a-d|)
                                        t = abs(a-d);
                                        j = n_remove(&t,levels[l]);
                                        t = n_pow(levels[l],j);

                                        arb_mul(b2,acb_realref(gval0),temp,prec);
                                        arb_div_ui(b2,b2,t,prec);
                                        arb_sub(b2,b2,inte0,prec);
                                        arb_add(arb_mat_entry(hypervalues[l],hecken-nmin,0),arb_mat_entry(hypervalues[l],hecken-nmin,0),b2,prec);

                                        arb_mul(b2,acb_realref(gval2),temp,prec);
                                        arb_div_ui(b2,b2,t,prec);
                                        arb_sub(b2,b2,inte2,prec);
                                        arb_add(arb_mat_entry(hypervalues[l],hecken-nmin,1),arb_mat_entry(hypervalues[l],hecken-nmin,1),b2,prec);

                                        arb_mul(b2,acb_realref(gval4),temp,prec);
                                        arb_div_ui(b2,b2,t,prec);
                                        arb_sub(b2,b2,inte4,prec);
                                        arb_add(arb_mat_entry(hypervalues[l],hecken-nmin,2),arb_mat_entry(hypervalues[l],hecken-nmin,2),b2,prec);
                                    } else {
                                        //Computing gcd(N^{infinity},|a-d|)
                                        t = abs(a+d);
                                        j = n_remove(&t,levels[l]);
                                        t = n_pow(levels[l],j);

                                        arb_mul(b2,acb_realref(gval0),temp,prec);
                                        arb_div_ui(b2,b2,t,prec);
                                        arb_sub(b2,b2,inte0,prec);
                                        arb_add(arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 0), arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 0), b2, prec);

                                        arb_mul(b2,acb_realref(gval2),temp,prec);
                                        arb_div_ui(b2,b2,t,prec);
                                        arb_sub(b2,b2,inte2,prec);
                                        arb_add(arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 1), arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 1), b2, prec);

                                        arb_mul(b2,acb_realref(gval4),temp,prec);
                                        arb_div_ui(b2,b2,t,prec);
                                        arb_sub(b2,b2,inte4,prec);
                                        arb_add(arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 2), arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 2), b2, prec);
                                    }
                                }
                            }
                        }
                    } else { // When a==d that also means that hecken is a square.
                        for(l=0;l<levelcount;l++){
                            // Parabolic terms only appear if level is prime
                            if(gcd(levels[l],labs(hecken)) == 1){
                                if(n_is_prime(levels[l])){
                                    arb_zero(inte0);
                                    arb_zero(inte2);
                                    arb_zero(inte4);
                                    // Setting temp to Lambda(N) = log(N)
                                    arb_log_ui(temp,levels[l],prec);
                                    D = floor(g.d*scale/(2*log(levels[l]))) + 1;
                                    for(t = 0; t < D; t++){
                                        j = floor(2*t*log(levels[l])/scale);
                                        //Compute 2*t*log(N)
                                        arb_mul_ui(b2,temp,2*t,prec);
                                        acb_set_arb(ctemp,b2);

                                        //Computing N^{t}
                                        arb_ui_pow_ui(b2,levels[l],t,prec);

                                        //Computing N^{-t} * g(2*t*log(N)) and adding it to g_const sum.
                                        evaltfd(gval,&g,j,ctemp);
                                        arb_div(b3,acb_realref(gval),b2,prec);
                                        arb_add(inte0,inte0,b3,prec);

                                        evaltfd(gval,&gprime2,j,ctemp);
                                        arb_div(b3,acb_realref(gval),b2,prec);
                                        arb_add(inte2,inte2,b3,prec);

                                        evaltfd(gval,&gprime4,j,ctemp);
                                        arb_div(b3,acb_realref(gval),b2,prec);
                                        arb_add(inte4,inte4,b3,prec);
                                    }
                                    //Scaling the sum by 2*log(N)
                                    arb_mul(inte0,inte0,temp,prec);
                                    arb_mul_2exp_si(inte0,inte0,1);
                                    arb_mul(inte2,inte2,temp,prec);
                                    arb_mul_2exp_si(inte2,inte2,1);
                                    arb_mul(inte4,inte4,temp,prec);
                                    arb_mul_2exp_si(inte4,inte4,1);

                                    // Adding parabolic terms to trace formula values.
                                    if(hecken > 0){
                                        arb_sub(arb_mat_entry(hypervalues[l],hecken-nmin,0),arb_mat_entry(hypervalues[l],hecken-nmin,0),inte0,prec);
                                        arb_sub(arb_mat_entry(hypervalues[l],hecken-nmin,1),arb_mat_entry(hypervalues[l],hecken-nmin,1),inte2,prec);
                                        arb_sub(arb_mat_entry(hypervalues[l],hecken-nmin,2),arb_mat_entry(hypervalues[l],hecken-nmin,2),inte4,prec);
                                    } else {
                                        //Computing gcd(N^{infinity},|a-d|)
                                        t = abs(a+d);
                                        j = n_remove(&t,levels[l]);
                                        t = n_pow(levels[l],j);

                                        arb_mul(b2,g_0,temp,prec);
                                        arb_div_ui(b2,b2,t,prec);
                                        arb_sub(b2,b2,inte0,prec);
                                        arb_add(arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 0), arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 0), b2, prec);

                                        arb_mul(b2,g2_0,temp,prec);
                                        arb_div_ui(b2,b2,t,prec);
                                        arb_sub(b2,b2,inte2,prec);
                                        arb_add(arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 1), arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 1), b2, prec);

                                        arb_mul(b2,g4_0,temp,prec);
                                        arb_div_ui(b2,b2,t,prec);
                                        arb_sub(b2,b2,inte4,prec);
                                        arb_add(arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 2), arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 2), b2, prec);
                                    }
                                }
                                //Computing the identity term and adding it to the trace formula values.
                                if(hecken > 0){
                                    Nprod = 1;
                                    for(t=0;t<Npcount[l];t++){
                                        Nprod *= (p[l][t] - 1);
                                    }

                                    arb_mul_ui(b1,idenint0,Nprod,prec);
                                    arb_div_si(b1,b1,a,prec);
                                    arb_add(arb_mat_entry(hypervalues[l],hecken-nmin,0),arb_mat_entry(hypervalues[l],hecken-nmin,0),b1,prec);

                                    arb_mul_ui(b1,idenint2,Nprod,prec);
                                    arb_div_si(b1,b1,a,prec);
                                    arb_add(arb_mat_entry(hypervalues[l],hecken-nmin,1),arb_mat_entry(hypervalues[l],hecken-nmin,1),b1,prec);

                                    arb_mul_ui(b1,idenint4,Nprod,prec);
                                    arb_div_si(b1,b1,a,prec);
                                    arb_add(arb_mat_entry(hypervalues[l],hecken-nmin,2),arb_mat_entry(hypervalues[l],hecken-nmin,2),b1,prec);
                                }
                            }
                        }
                    }
                }
                for(l=0;l<levelcount;l++){
                    if(gcd(levels[l],labs(hecken)) == 1){
                        //Compute mu(N)*sigma_1(|n|)/sqrt(|n|) * h(i/2)
                        D = sigma1 * muN[l];
                        arb_mul_si(b1,acb_realref(hi2),D,prec);
                        arb_rsqrt_ui(b2,labs(hecken),prec);
                        arb_mul(b1,b1,b2,prec);

                        // Subtracting h(i/2) eigenvalue terms from trace formula values.
                        if(hecken > 0){
                            arb_sub(arb_mat_entry(hypervalues[l],hecken-nmin,0),arb_mat_entry(hypervalues[l],hecken-nmin,0),b1,prec);

                            arb_mul_2exp_si(b1,b1,-2);
                            arb_sub(arb_mat_entry(hypervalues[l],hecken-nmin,1),arb_mat_entry(hypervalues[l],hecken-nmin,1),b1,prec);

                            arb_mul_2exp_si(b1,b1,-2);
                            arb_sub(arb_mat_entry(hypervalues[l],hecken-nmin,2),arb_mat_entry(hypervalues[l],hecken-nmin,2),b1,prec);
                        } else {
                            arb_sub(arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 0), arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 0), b1, prec);

                            arb_mul_2exp_si(b1,b1,-2);
                            arb_sub(arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 1), arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 1), b1, prec);

                            arb_mul_2exp_si(b1,b1,-2);
                            arb_sub(arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 2), arb_mat_entry(hypervalues[l], labs(hecken) - labs(nmax), 2), b1, prec);
                        }
                    }
                }
            }
            /* End of parabolic and identity terms */

            // Dumping code into temporary files. We will then use tracecomp.sh to reorganise this data so that each file is for a different level.
            flock(1,LOCK_EX);
            printf("Proc = %d done \n",proc);
            sprintf(line,"./tfdata/%ld.tmp",proc);
            tmp = fopen(line,"w");
            for(N=0;N<levelcount;N++){
                for(j=0;j<labs(nmax-nmin)+1;j++){
                    if(proc < NPROCS/2){
                        fprintf(tmp,"%ld,%ld,",levels[N],-(j + labs(nmax)));
                    } else {
                        fprintf(tmp,"%ld,%ld,",levels[N],j + labs(nmin));
                    }
                    arb_dump_file(tmp,arb_mat_entry(hypervalues[N],j,0));
                    fprintf(tmp,",");
                    arb_dump_file(tmp,arb_mat_entry(hypervalues[N],j,1));
                    fprintf(tmp,",");
                    arb_dump_file(tmp,arb_mat_entry(hypervalues[N],j,2));
                    fprintf(tmp,"\n");
                }
            }
            fclose(tmp);
            //fflush(stdout);
            flock(1,LOCK_UN);
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
    while(wait(&k) > 0);
    
    free(levels);
    return 0;
}
