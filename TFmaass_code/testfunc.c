#include <acb_poly.h>
#include <string.h>
#include <inttypes.h>

// Since this is called by other files, we do not need to define a precision

typedef struct {
	slong d;
	arb_poly_t *A0;
	acb_poly_t *A1;
} tf;
typedef struct {
	slong d;
    arb_t a;
	arb_poly_t *A0;
	acb_poly_t *A1;
} tfd;

/* Start of code to generate the test function via convolutions.*/
/* Code provided by Andrew Booker.*/
static void acb_poly_conj(acb_poly_t res,const acb_poly_t f) {
	slong i;
	acb_ptr p;
	acb_poly_set(res,f);
	for (i=acb_poly_degree(res);i>=0;i--) {
		p = acb_poly_get_coeff_ptr(res,i);
		acb_conj(p,p);
	}
}

// return f(-x)
static void acb_poly_twist(acb_poly_t res,const acb_poly_t f) {
	slong i;
	acb_ptr p;
	acb_poly_set(res,f);
	i=acb_poly_degree(res);
	if (!(i&1)) i--;
	for (;i>=0;i-=2) {
		p = acb_poly_get_coeff_ptr(res,i);
		acb_neg(p,p);
	}
}

static void getamj(acb_poly_t amj,const tf *a,int m,int j) {
	if (j >= 0) {
		if (!m) acb_poly_set_arb_poly(amj,a->A0[j]);
		if (m < 0) acb_poly_conj(amj,a->A1[j]);
		if (m > 0) acb_poly_set(amj,a->A1[j]);
	} else {
		if (!m) acb_poly_set_arb_poly(amj,a->A0[~j]);
		if (m < 0) acb_poly_set(amj,a->A1[~j]);
		if (m > 0) acb_poly_conj(amj,a->A1[~j]);
		acb_poly_twist(amj,amj);
	}
}

// f^{(n)}/n!
static void normalized_deriv(acb_poly_t res,const acb_poly_t f,int n) {
	fmpz_t z;
	int i,l;

	l = acb_poly_length(f)-n;
	if (l <= 0) {
		acb_poly_zero(res);
		return;
	}
	if (res != f) {
		acb_poly_fit_length(res,l);
		_acb_poly_set_length(res,l);
	}

	fmpz_init(z);
	for (i=0;i<l;i++) {
		fmpz_bin_uiui(z,i+n,n);
		acb_mul_fmpz(acb_poly_get_coeff_ptr(res,i),
			acb_poly_get_coeff_ptr(f,i+n),z,prec);
	}
	fmpz_clear(z);
	_acb_poly_set_length(res,l);
}

static void acb_poly_get_arb_poly(arb_poly_t res,const acb_poly_t f) {
	slong i,l;
	l = acb_poly_length(f);
	arb_poly_fit_length(res,l);
	_arb_poly_set_length(res,l);
	for (i=0;i<l;i++)
		arb_poly_set_coeff_arb(res,i,acb_realref(acb_poly_get_coeff_ptr(f,i)));
}

void tfinit(tf *a) {
	a->d = 0;
	a->A0 = (arb_poly_t *)0;
	a->A1 = (acb_poly_t *)0;
}

void tfclear(tf *a) {
	int j;
	if (a->d) {
		for (j=0;j<a->d;j++) {
			arb_poly_clear(a->A0[j]);
			acb_poly_clear(a->A1[j]);
		}
		free(a->A0);
	}
	a->d = 0;
}

void tfeval(arb_t,tf *,arb_srcptr);
void convolution(tf *res,const tf *a,const tf *b) {
	static int init;
	static acb_poly_t amj,bnk,amjr,bnks,ctemp;
	static arb_poly_t rtemp;
	static acb_t t1,t2;
	static fmpz_t z;
	tf c;
	int j,k,m,n,r,s,delta;

	if (!init) {
		acb_poly_init(amj); acb_poly_init(bnk);
		acb_poly_init(amjr); acb_poly_init(bnks);
		acb_poly_init(ctemp); arb_poly_init(rtemp);
		acb_init(t1); acb_init(t2);
		fmpz_init(z);
		init = 1;
	}

	c.d = a->d+b->d;
	c.A0 = (arb_poly_t *)malloc(c.d*(sizeof(c.A0[0])+sizeof(c.A1[0])));
	c.A1 = (acb_poly_t *)(c.A0+c.d);
	for (j=0;j<c.d;j++) {
		arb_poly_init(c.A0[j]);
		acb_poly_init(c.A1[j]);
	}

	for (j=-a->d;j<a->d;j++)
	for (m=-1;m<=1;m++) {
		getamj(amj,a,m,j);
		for (k=-b->d;k<b->d;k++)
		if (j+k >= -1)
		for (n=-(m>=0);n<=1;n++) {
			getamj(bnk,b,n,k);
			if (m == n)
			for (delta=(j+k<0);delta<=1;delta++) {
				acb_poly_set(amjr,amj);
				acb_poly_set(bnks,bnk);
				while (!acb_poly_is_zero(bnks)) {
					acb_set_d(t1,delta-0.5);
					acb_poly_integral(amjr,amjr,prec);
					acb_poly_evaluate(t2,amjr,t1,prec);
					acb_neg(t2,t2);
					acb_poly_set_coeff_acb(amjr,0,t2);
					acb_poly_evaluate(t2,bnks,t1,prec);
					if (delta) acb_neg(t2,t2);
					acb_poly_scalar_mul(ctemp,amjr,t2,prec);
					if (m)
						acb_poly_add(c.A1[j+k+delta],c.A1[j+k+delta],ctemp,prec);
					else {
						acb_poly_get_arb_poly(rtemp,ctemp);
						arb_poly_add(c.A0[j+k+delta],c.A0[j+k+delta],rtemp,prec);
					}
					acb_poly_derivative(bnks,bnks,prec);
				}
			} else {
				for (r=0;r<=acb_poly_degree(amj);r++) {
					normalized_deriv(amjr,amj,r);
					for (s=0;s<=acb_poly_degree(bnk);s++) {
						normalized_deriv(bnks,bnk,s);

						// (-1)^s*(r+s)!/(Pi*I*(n-m))^(r+s+1)
						arb_const_pi(acb_realref(t1),prec);
						arb_mul_si(acb_realref(t1),acb_realref(t1),n-m,prec);
						arb_zero(acb_imagref(t1));
						acb_mul_onei(t1,t1);
						acb_pow_si(t1,t1,-1-r-s,prec);
						fmpz_fac_ui(z,r+s);
						if (s & 1) fmpz_neg(z,z);
						acb_mul_fmpz(t1,t1,z,prec);

						for (delta=(j+k<0);delta<=1;delta++) {
							if (m >= 0) {
								acb_set_d(t2,delta-0.5);
								acb_poly_evaluate(t2,bnks,t2,prec);
								if (((m-n)*(k+delta)+delta+1) & 1)
									acb_neg(t2,t2);
								acb_mul(t2,t2,t1,prec);
								acb_poly_scalar_mul(ctemp,amjr,t2,prec);
								if (m)
									acb_poly_add(c.A1[j+k+delta],c.A1[j+k+delta],ctemp,prec);
								else {
									acb_poly_get_arb_poly(rtemp,ctemp);
									arb_poly_add(c.A0[j+k+delta],c.A0[j+k+delta],rtemp,prec);
								}
							}
							if (n >= 0) {
								acb_set_d(t2,delta-0.5);
								acb_poly_evaluate(t2,amjr,t2,prec);
								if (((m-n)*(j+delta)+delta) & 1)
									acb_neg(t2,t2);
								acb_mul(t2,t2,t1,prec);
								acb_poly_scalar_mul(ctemp,bnks,t2,prec);
								if (n)
									acb_poly_add(c.A1[j+k+delta],c.A1[j+k+delta],ctemp,prec);
								else {
									acb_poly_get_arb_poly(rtemp,ctemp);
									arb_poly_add(c.A0[j+k+delta],c.A0[j+k+delta],rtemp,prec);
								}
							}
						}
					}
				}
			}
		}
	}

	tfclear(res);
	memcpy(res,&c,sizeof(c));
}

static void h1(tf *a) {
	arb_ptr p;

	a->d = 1;
	a->A0 = (arb_poly_t *)malloc(sizeof(a->A0[0])+sizeof(a->A1[0]));
	a->A1 = (acb_poly_t *)(a->A0+1);

	arb_poly_init(a->A0[0]);
	arb_poly_fit_length(a->A0[0],2);
	_arb_poly_set_length(a->A0[0],2);
	p = arb_poly_get_coeff_ptr(a->A0[0],0);

	// Pi^2/(Pi^2+4)
	arb_const_pi(p,prec);
	arb_inv(p,p,prec);
	arb_mul_2exp_si(p,p,1);
	arb_sqr(p,p,prec);
	arb_add_ui(p,p,1,prec);
	arb_inv(p,p,prec);

	arb_neg(arb_poly_get_coeff_ptr(a->A0[0],1),p);
	arb_mul_2exp_si(p,p,-1);

	acb_poly_init(a->A1[0]);
	acb_poly_set_arb_poly(a->A1[0],a->A0[0]);
	acb_poly_scalar_mul_2exp_si(a->A1[0],a->A1[0],-1);
}

void hd(tf *res,int d) {
	tf a;
	h1(res);
	if (d > 1) {
		h1(&a);
		while (--d)
			convolution(res,res,&a);
		tfclear(&a);
	}
}

/* End of code to generate the test function via convolutions. */

// Computes the test function h_1(t)^k
void h_k(acb_t res, acb_t z, int k){
    static arb_t piconst, pi2;
    static acb_t t;
    static int init;
    if(!init){
        arb_init(piconst);
        arb_init(pi2);
        acb_init(t);
    }

    arb_const_pi(pi2, prec);
    arb_mul(pi2, pi2, pi2, prec);
    arb_set(piconst, pi2);
    arb_add_ui(piconst, piconst, 4, prec);
    arb_div(piconst, pi2, piconst, prec);

    acb_set(t, z);
    acb_mul_2exp_si(t,t,-1);
    acb_sinc(t, t, prec);
    acb_mul(res,t,t,prec);

    acb_const_pi(t,prec);
    acb_sub(t, z, t, prec);
    acb_mul_2exp_si(t,t,-1);
    acb_sinc(t, t, prec);
    acb_mul(t,t,t,prec);
    acb_mul_2exp_si(t,t,-1);
    acb_add(res, res, t, prec);

    acb_const_pi(t,prec);
    acb_add(t, z, t, prec);
    acb_mul_2exp_si(t,t,-1);
    acb_sinc(t, t, prec);
    acb_mul(t,t,t,prec);
    acb_mul_2exp_si(t,t,-1);
    acb_add(res, res, t, prec);

    acb_mul_arb(res, res, piconst,prec);
    acb_pow_ui(res, res, k, prec);
}

//Shift test function so it agress with g_d(x)
void shifttf(void *param){
	int j;
    acb_t temp;
    acb_init(temp);
	tf *a = (tf *)param;
    for(j = 0; j< a->d; j++){
        acb_set_si(temp,-1-2*j);
        acb_mul_2exp_si(temp, temp, -1); //maybe need some more prec
        arb_poly_taylor_shift(a->A0[j], a->A0[j],acb_realref(temp), prec);
		acb_poly_taylor_shift(a->A1[j], a->A1[j],temp, prec);
    }
    acb_clear(temp);
}

//Finds the derivative of h and stores it as a tf
void hderiv(tf *hprime, tf *h){
	acb_t z1;
	arb_poly_t a0;
	acb_poly_t a1,a2;
	acb_init(z1);
	int j;

	arb_poly_init(a0);
	acb_poly_init(a1);
	acb_poly_init(a2);
    hprime->d = h->d;
	hprime->A0 = (arb_poly_t *)malloc((h->d)*sizeof(h->A0[0]));
	hprime->A1 = (acb_poly_t *)malloc((h->d)*sizeof(h->A1[0]));


	for(j=0;j<h->d;j++){
		arb_poly_set(a0, h->A0[j]);
		acb_poly_set(a1, h->A1[j]);

		arb_poly_init(hprime->A0[j]);
		acb_poly_init(hprime->A1[j]);

		arb_poly_derivative(hprime->A0[j], a0 ,prec);
		acb_poly_derivative(a2, a1, prec);

		acb_const_pi(z1, prec);
		acb_mul_onei(z1,z1);

		acb_poly_scalar_mul(a1, a1, z1, prec);
		acb_poly_add(hprime->A1[j], a2, a1, prec);
	}
}

//Takes test function g(x) and return test function g(x/a)/a.
void dilateg(tfd *g,tf *h,arb_t a){
    arb_poly_t a0;
	acb_poly_t a1;
	acb_t z1;
	arb_t b1;
    int j;

	arb_poly_init(a0);
	acb_poly_init(a1);
	acb_init(z1);
	arb_init(b1);

    g->A0 = (arb_poly_t *)malloc((h->d)*sizeof(h->A0[0]));
	g->A1 = (acb_poly_t *)malloc((h->d)*sizeof(h->A1[0]));

	arb_inv(b1,a,prec);
    acb_set_arb(z1,b1);
	arb_poly_zero(a0);
	arb_poly_set_coeff_arb(a0,1,b1);
    acb_poly_set_arb_poly(a1,a0);

    for(j=0;j<h->d;j++){
        arb_poly_init(g->A0[j]);
		acb_poly_init(g->A1[j]);

        arb_poly_compose(g->A0[j],h->A0[j],a0,prec);
        acb_poly_compose(g->A1[j],h->A1[j],a1,prec);

        arb_poly_scalar_mul(g->A0[j],g->A0[j],b1,prec);
        acb_poly_scalar_mul(g->A1[j],g->A1[j],z1,prec);
    }
    arb_init(g->a);
    arb_set(g->a,a);
    g->d = h->d;
}

// Computes the test function g_d(x) with x in [j,j+1).
void evaltfd(acb_t res, void *param, int j, acb_srcptr x){
	static int init;
	static acb_t z1,z2;
	tfd *h = (tfd *)param;
	if(!init){
		acb_init(z1); acb_init(z2);
		init = 1;
	}
	arb_poly_evaluate_acb(res, h->A0[j], x, prec);
	acb_poly_evaluate(z1, h->A1[j], x, prec);
	acb_div_arb(z2,x,h->a,prec);
	acb_exp_pi_i(z2,z2,prec);
	acb_mul(z1, z1, z2, prec);
	acb_add(res, res, z1, prec);

	acb_conj(z1, z1);
	acb_add(res, res, z1, prec);
}

// Prints parts of the test funcion g_d.
void printingdg(void *param){
	tfd *a = (tfd *)param;
    flint_printf("d = %d\n", a->d);
    int j;
    for(j = 0; j< a->d; j++){
        flint_printf("Polynomial of A_0 for j = %d is ",j);
        arb_poly_printd(a->A0[j],10);
        flint_printf("\n");
		flint_printf("Polynomial of A_1 for j = %d is ",j);
        acb_poly_printd(a->A1[j],5);
        flint_printf("\n");
    }
}

//Calculates the derivative of g_d and stores it in hprime.
void derivgd(tfd *hprime, tfd *h){
	acb_t z1;
	arb_poly_t a0;
	acb_poly_t a1,a2;
	acb_init(z1);
	int j;

	arb_poly_init(a0);
	acb_poly_init(a1);
	acb_poly_init(a2);
    hprime->d = h->d;
	arb_init(hprime->a);
	arb_set(hprime->a,h->a);
	hprime->A0 = (arb_poly_t *)malloc((h->d)*sizeof(h->A0[0]));
	hprime->A1 = (acb_poly_t *)malloc((h->d)*sizeof(h->A1[0]));


	for(j=0;j<h->d;j++){
		arb_poly_set(a0, h->A0[j]);
		acb_poly_set(a1, h->A1[j]);

		arb_poly_init(hprime->A0[j]);
		acb_poly_init(hprime->A1[j]);

		arb_poly_derivative(hprime->A0[j], a0 ,prec);
		acb_poly_derivative(a2, a1, prec);

		acb_const_pi(z1, prec);
		acb_mul_onei(z1,z1);
		acb_div_arb(z1,z1,h->a,prec);

		acb_poly_scalar_mul(a1, a1, z1, prec);
		acb_poly_add(hprime->A1[j], a2, a1, prec);
	}
}

//Function that agrees with g'(x) on [aj,a(j+1)) , j>0 and g'(x)/x on [0,a) for dilated g
//where a is the dilation factor.
void evaltfdx(acb_t res, void *param, int j, acb_srcptr x){
	static int init;
	static acb_t z1,z2;
	static acb_poly_t a0,a1;
	tfd *h = (tfd *)param;
	if(!init){
		acb_init(z1); acb_init(z2);
		acb_poly_init(a0); acb_poly_init(a1);
		init = 1;
	}

	acb_poly_derivative(a1,h->A1[j],prec);
	acb_const_pi(z1, prec);
	acb_mul_onei(z1,z1);
	acb_div_arb(z1,z1,h->a,prec);
	acb_poly_scalar_mul(a0,h->A1[j],z1,prec);
	acb_poly_add(a1,a1,a0,prec);

	acb_poly_set_arb_poly(a0,h->A0[j]);
	acb_poly_derivative(a0,a0,prec);

	if (j) {
		acb_poly_evaluate(res,a0,x,prec);

		acb_div_arb(z2,x,h->a,prec);
		acb_exp_pi_i(z2,z2,prec);
		acb_poly_evaluate(z1,a1,x,prec);
		acb_mul(z1,z1,z2,prec);
		acb_add(res,res,z1,prec);

		acb_conj(z1,z1);
		acb_add(res,res,z1,prec);
		return;
	}
	acb_poly_shift_right(a0,a0,1); //Dividing by x
	acb_poly_evaluate(res,a0,x,prec);

	acb_poly_conj(a0,a1);
	acb_poly_add(a0,a0,a1,prec);
	acb_poly_get_coeff_acb(z2,a0,0);

	acb_div_arb(z1,x,h->a,prec);
	acb_mul_2exp_si(z1,z1,-1);
	acb_sinc_pi(z1,z1,prec);
	acb_mul(z2,z2,z1,prec);
	acb_div_arb(z1,x,h->a,prec);
	acb_mul_2exp_si(z1,z1,-1);
	acb_sin_pi(z1,z1,prec);
	acb_mul(z2,z2,z1,prec);
	acb_const_pi(z1,prec);
	acb_mul(z2,z2,z1,prec);
	acb_div_arb(z2,z2,h->a,prec);
	acb_sub(res,res,z2,prec);

	acb_poly_shift_right(a0,a0,1);
	acb_poly_evaluate(z1,a0,x,prec);
	acb_div_arb(z2,x,h->a,prec);
	acb_cos_pi(z2,z2,prec);
	acb_mul(z2,z2,z1,prec);
	acb_add(res,res,z2,prec);

	acb_poly_conj(a0,a1);
	acb_poly_sub(a0,a1,a0,prec);
	acb_poly_evaluate(z1,a0,x,prec);
	acb_const_pi(z2,prec);
	acb_div_arb(z2,z2,h->a,prec);
	acb_mul(z1,z1,z2,prec);
	acb_mul(z2,z2,x,prec);
	acb_sinc(z2,z2,prec);
	acb_mul(z1,z1,z2,prec);
	acb_mul_onei(z1,z1);
	acb_add(res,res,z1,prec);
}

//Function that agrees with g(x+y)-g(y) on [aj-y,a(j+1)-y),j>0 and (g(x+y)-g(y))/x on [0,a) for dilated g
//where a is the dilation factor and y is a real number.
void evaltfdxy(acb_t res, void *param, int j, acb_srcptr x, acb_t y){
	static int init,n,k,m;
	static slong jfloor;
	static acb_t z1,z2,s,c,s2,c2,snc,ukprod;
	static arb_t ak,bk,ck;
	static arb_poly_t a0;
	static acb_poly_t a1,a2;
	tfd *h = (tfd *)param;
	if(!init){
		acb_init(z1); acb_init(z2);
		acb_init(s); acb_init(c);
		acb_init(s2); acb_init(c2);
		acb_init(snc); acb_init(ukprod);
		arb_init(ak); arb_init(bk);
		arb_init(ck);
		arb_poly_init(a0); acb_poly_init(a1); acb_poly_init(a2);
		init = 1;
	}

	if (j) {
		acb_add(z1,x,y,prec);
		return;
	}

	arb_floor(ak,acb_realref(y),prec);
	jfloor = arf_get_si(arb_midref(ak),ARF_RND_FLOOR);
	printf("jfloor = %ld\n",jfloor);
	if(jfloor>h->d){
		return;
	}

	arb_poly_set(a0,h->A0[jfloor]);
	acb_poly_set(a1,h->A1[jfloor]);

	n = acb_poly_degree(a1);

	//create 2*cos((pi*(u+y))/a) and 2*sin((pi*(u+y))/a)
	acb_add(z1,x,y,prec);
	acb_div_arb(z1,z1,h->a,prec);
	acb_sin_cos_pi(s,c,z1,prec);
	acb_mul_2exp_si(s,s,1);
	acb_mul_2exp_si(c,c,1);

	//create cos((pi*(u+2*y))/2a) and sin((pi*(u+2*y))/2a)
	acb_mul_2exp_si(z1,x,-1);
	acb_add(z1,z1,y,prec);
	acb_div_arb(z1,z1,h->a,prec);
	acb_sin_cos_pi(s2,c2,z1,prec);

	//create (2*pi/a)*sinc(u*pi/2a).
	acb_div_arb(z1,x,h->a,prec);
	acb_mul_2exp_si(z1,z1,-1);
	acb_sinc_pi(snc,z1,prec);
	acb_const_pi(z1, prec);
	acb_mul_2exp_si(z1,z1,1);
	acb_div_arb(z1,z1,h->a,prec);
	acb_mul(snc,snc,z1,prec);

	//a2 poly is y^k poly.
	acb_poly_zero(a2);
	acb_zero(res);
	for(k=0;k<=n;k++){
		acb_poly_get_coeff_acb(z1, a1, k);
		acb_get_real(ak, z1);
		acb_get_imag(bk, z1);
		acb_mul_arb(z1,s2,ak,prec);
		acb_mul_arb(z2,c2,bk,prec);
		acb_add(z1,z1,z2,prec);
		acb_poly_set_coeff_acb(a2, k, z1);

		if(k==1){
			arb_poly_get_coeff_arb(ck,a0,k);
			acb_mul_arb(z1,c,ak,prec);
			acb_mul_arb(z2,s,bk,prec);
			acb_sub(z1,z1,z2,prec);
			acb_add_arb(z1,z1,ck,prec);
			acb_set(res,z1);
		}else if (k>1){
			acb_zero(ukprod);
			for(m=0;m<=k-1;m++){
				acb_pow_ui(z1,y,m,prec);
				acb_pow_ui(z2,x,k-m-1,prec);
				acb_mul(z1,z1,z2,prec);
				acb_zero(z2);
				arb_bin_uiui(acb_realref(z2),k,m,prec);
				acb_mul(z1,z1,z2,prec);
				acb_add(ukprod,ukprod,z1,prec);
			}

			arb_poly_get_coeff_arb(ck,a0,k);
			acb_mul_arb(z1,c,ak,prec);
			acb_mul_arb(z2,s,bk,prec);
			acb_sub(z1,z1,z2,prec);
			acb_add_arb(z1,z1,ck,prec);

			acb_mul(z1,z1,ukprod,prec);
			acb_add(res,res,z1,prec);
		}
	}

	acb_poly_evaluate(z1,a2,y,prec);
	acb_mul(z1,z1,snc,prec);
	acb_sub(res,res,z1,prec);
}

// These integer functions are taken from Andy's opt.c
static long gcd(long x,long y) {
	long t;
	while (x)
		t = y % x, y = x, x = t;
	return y;
}

static inline int ordp(unsigned long x,unsigned long p) {
	int k;
	for (k=0;!(x%p);k++)
		x /= p;
	return k;
}

static inline int modmul(int x,int y,int p) {
	return (int)((long)x*y % p);
}

static int modpow(int x,int n,int p) {
	int a = 1;
	while (n > 0) {
		if (n & 1) a = modmul(a,x,p);
		x = modmul(x,x,p);
		n >>= 1;
	}
	return a;
}

static int chi8[] = {0,1,0,-1,0,-1,0,1};
static int kronecker(long int d,int p) {
	int dmodp;
	if (p == 2) return chi8[d & 7];
	if ( !(dmodp = d % p) ) return 0;
	if (dmodp < 0) dmodp += p;
	if (modpow(dmodp,(p-1)>>1,p) == 1)
		return 1;
	return -1;
}

// Recusively generate divisors from prime factorisation.
void generate_divisors(uint64_t *p, int *e, int k, int curindex, uint64_t curdiv, uint64_t *d, int *divindex){
    if(curindex == k){
        d[*divindex] = curdiv;
        *divindex += 1;
        return;
    }

    for(int i=0; i <= e[curindex]; i++){
        generate_divisors(p,e,k,curindex+1,curdiv,d,divindex);
        curdiv *= p[curindex]; 
    }
}

int compare (const void * a, const void * b)
{
  if ( *(uint64_t*)a <  *(uint64_t*)b ) return -1;
  if ( *(uint64_t*)a == *(uint64_t*)b ) return 0;
  if ( *(uint64_t*)a >  *(uint64_t*)b ) return 1;
  return 0;
}

// Compute divisors of number n from known prime factorisation.
// Returns the number of divisors
int compute_divisors(uint64_t n,uint64_t *p,int *e, int k, uint64_t *d){
    int i,divisor_count=1;
    
    int *divindex = (int *)malloc(sizeof(int));  
    *divindex=0;
    
    // Compute the number of divisors to malloc the d.
    for(i=0;i<k;i++){
        divisor_count *= (e[i] + 1);
    }
    //d = (uint64_t *)malloc(divisor_count * sizeof(uint64_t));
    generate_divisors(p,e,k,0,1,d,divindex);
    qsort(d,divisor_count,sizeof(uint64_t),compare);
    
    return divisor_count;
}
