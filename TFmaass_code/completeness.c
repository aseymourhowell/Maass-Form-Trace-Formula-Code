#include <acb_poly.h>
#include <arb_mat.h>
#include <time.h>
#include <dirent.h>
#include <string.h>

#define MAXLEVEL 105

static slong prec = 300;
#include "testfunc.c"

int descalphasort(const struct dirent **a, const struct dirent **b)
{
    return alphasort(b, a);
}

//Uses the bisect method to find x when h(x) = y and y known, in the interval [a,b]. f(x)=h(x)-y.
void hbisect(arb_t x,arb_t a,arb_t b,arb_t y,int TFpow,arb_t gscale) {
    static arb_t a1,b1,c1,f_a,f_b,f_c;
    static acb_t z1,z2;
    static int init;
    int n,nmax;
    if(!init){
        arb_init(a1); arb_init(b1); arb_init(c1);
        arb_init(f_a); arb_init(f_b); arb_init(f_c);
        acb_init(z1); acb_init(z2);
        init=1;
    }

    if(arb_lt(b,a)){
        printf("Error: a < b for bisection\n");
        return;
    }

    //Find f(a)
    acb_set_arb(z1,a);
    acb_mul_arb(z1,z1,gscale,prec);
    h_k(z2,z1,TFpow);
    arb_sub(f_a,acb_realref(z2),y,prec);
    //Find f(b)
    acb_set_arb(z1,b);
    acb_mul_arb(z1,z1,gscale,prec);
    h_k(z2,z1,TFpow);
    arb_sub(f_b,acb_realref(z2),y,prec);

    arb_set(a1,a);
    arb_set(b1,b);

    nmax=1000;
    n=0;
    while(n<nmax){
        //c=(a+b)/2
        arb_add(c1,a1,b1,prec);
        arb_mul_2exp_si(c1,c1,-1);

        //Find f(c)
        acb_zero(z2);
        acb_set_arb(z1,c1);
        acb_mul_arb(z1,z1,gscale,prec);
        h_k(z2,z1,TFpow);
        arb_sub(f_c,acb_realref(z2),y,prec);

        if(arb_is_zero(f_c)){
            arb_set(x,c1);
            return;
        }else if(arb_sgn_nonzero(f_a)==arb_sgn_nonzero(f_c)){
            arb_set(a1,c1);
            arb_set(f_a,f_c);
        }else{
            arb_set(b1,c1);
            arb_set(f_b,f_c);
        }
        n++;
    }
    arb_set(x,c1);
}

// Read trace formula values for T_{\pm 1} operators.
void readT1value(arb_t res,long level){
    char filename[30],line[3600],gstr0[200],gstr2[200],gstr4[200];;
    sprintf(filename,"./tfdata/t1_13/0.tmp",level); // REMEMBER TO CHANGE THIS
    FILE *fp = fopen(filename,"r");
    long N,n;
    while(fgets(line,sizeof(line),fp) != NULL){
        sscanf(line,"%li,%li,%100[^,],%100[^,],%100[^\n]",&N,&n,gstr0,gstr2,gstr4);
        if(N==level){
            if(n==-1){
            arb_load_str(res,gstr0);
            break;
            }
        } 
    }
    fclose(fp);
}

// Read trace formula values for T_{\pm 1} operators for a given parity. 0 for even forms, 1 for odd forms.
void readT1value_parity(arb_t res,long level,int parity){
    arb_t tval;
    arb_init(tval);
    char filename[30],line[3600],gstr0[200],gstr2[200],gstr4[200];;
    sprintf(filename,"./tfdata/t1_13/0.tmp",level); // REMEMBER TO CHANGE THIS
    FILE *fp = fopen(filename,"r");
    long N,n;
    while(fgets(line,sizeof(line),fp) != NULL){
        sscanf(line,"%li,%li,%100[^,],%100[^,],%100[^\n]",&N,&n,gstr0,gstr2,gstr4);
        if(N==level){
            if(n==-1){
            arb_load_str(tval,gstr0);
            break;
            }
        } 
    }
    fclose(fp);

    sprintf(filename,"./tfdata/t1_13/1.tmp",level);
    fp = fopen(filename,"r");
    while(fgets(line,sizeof(line),fp) != NULL){
        sscanf(line,"%li,%li,%100[^,],%100[^,],%100[^\n]",&N,&n,gstr0,gstr2,gstr4);
        if(N==level){
            if(n==1){
            arb_load_str(res,gstr0);
            if(parity == 0){
                arb_add(res,res,tval,prec);
                arb_mul_2exp_si(res,res,-1);
            } else if(parity == 1){
                arb_sub(res,tval,res,prec);
                arb_mul_2exp_si(res,res,-1);
            }
            break;
            }
        } 
    }
    fclose(fp);
    arb_clear(tval);
}

// Computes the range of R that we computed completeness for.
void completeness_fullspec(){
    arb_t Rval,eps;
    arb_init(Rval); arb_init(eps);
    char full[200],line[200], direc[30] = "./maassdata/";
    FILE *fp;
    struct dirent **namelist;
    int i,n,j;
    char *token;
    size_t last_idx;

    long level,levelcount=0,*levels = (long *)malloc(MAXLEVEL * sizeof(long)); //This is too many terms
    for(i = 2; i <= MAXLEVEL; i++){
        if(n_is_squarefree(i) == 1){
            levels[levelcount] = i;
            levelcount++;
        }
    }
    levels = (long *)realloc(levels, levelcount * sizeof(long)); // Realloc to correct number of terms
    long levelindexcounters[levelcount];
    n = scandir(direc, &namelist, 0,descalphasort);
    arb_mat_t Rvalsover[levelcount],Rvals[levelcount];
    for(j = 0; j < levelcount; j++){
        arb_mat_init(Rvalsover[j],n,3);
        levelindexcounters[j]=0;
    }

    // Loop through Maass form files.
    if (n < 0)
        perror("scandir");
    else{
        //n = n-2;
        while (n--) {
            strcpy(full,direc);
            strcat(full,namelist[n]->d_name);
            fp = fopen(full,"r");
            // Get R value
            if(fgets(line,sizeof(line),fp) != NULL){
                last_idx = strlen(line) - 1;

                if( line[last_idx] == '\n' ) {
                    line[last_idx] = '\0';
                }
                token = strtok(line, ",");
                arb_set_str(Rval,token,prec);
                while(token != NULL){
                    arb_set_str(eps,token,prec);
                    token = strtok(NULL, ",");
                }
                if(arf_is_nan(arb_midref(eps))){
                    fclose(fp);
                    free(namelist[n]);
                    continue;
                }
                arb_add_error(Rval,eps);
            }
            // Get level
            if(fgets(line,sizeof(line),fp) != NULL){
                sscanf(line,"%ld\n",&level);
                for(j=0;j<levelcount;j++){
                    if(levels[j]==level){
                        arb_set(arb_mat_entry(Rvalsover[j],levelindexcounters[j],0),Rval);
                        arb_set(arb_mat_entry(Rvalsover[j],levelindexcounters[j],1),eps);
                        levelindexcounters[j]++;
                        break;
                    }
                }
            }
            fclose(fp);
            free(namelist[n]);
        }
    free(namelist);
    }
    acb_t z1,z2;
    arb_t b1,S,gscale,t1val;
    acb_init(z1); acb_init(z2);
    arb_init(b1); arb_init(S); arb_init(gscale); arb_init(t1val);
    int TFpow = 13;
    arb_set_d(gscale,5.51341248666333);
    arb_div_ui(gscale,gscale,TFpow,prec);
    arb_t tol;
    arb_init(tol);
    arb_set_d(tol,pow(10,-2));
    for(j=0;j<levelcount;j++){
        arb_mat_window_init(Rvals[j],Rvalsover[j],0,0,levelindexcounters[j],2);
        arb_zero(b1);
        arb_zero(S);
        for(i=0;i<arb_mat_nrows(Rvals[j]);i++){
            if(arf_is_nan(arb_midref(arb_mat_entry(Rvals[j],i,1)))){
                continue;
            } else if(arb_lt(tol,arb_mat_entry(Rvals[j],i,1)) || arb_is_zero(arb_mat_entry(Rvals[j],i,1))){
                continue;
            }
            //arb_printd(arb_mat_entry(Rvals[j],i,0),16);
            //printf("\n");
            arb_add(b1, arb_mat_entry(Rvals[j],i,0),arb_mat_entry(Rvals[j],i,1),prec);
            arb_mul(b1,b1,gscale,prec);
            acb_set_arb(z1,b1);
            h_k(z2,z1,TFpow);
            arb_add(S,S,acb_realref(z2),prec);
        }
        readT1value(t1val,levels[j]);

        printf("TF1 = ");
        arb_printd(t1val,16);
        printf("\n");

        printf("level = %ld, S = ",levels[j]);
        arb_printd(S,16);
        printf("\n");

        arb_sub(S,t1val,S,prec);
        printf("diff = ");
        arb_printd(S,16);
        printf("\n");
        acb_set_d_d(z1,20,0);

        acb_printd(z1,16);
        printf("\n");

        hbisect(b1,acb_imagref(z1),acb_realref(z1),S,TFpow,gscale);

        printf("List is complete for all lambda < ");
        arb_printd(b1,16);
        printf("\n");
    }
}

// Computes the range of R that we computed completeness for, with a given parity. 0 for even forms, 1 for odd forms.
void completeness_parity(int sign_wanted){
    arb_t Rval,eps;
    arb_init(Rval); arb_init(eps);
    char full[200],line[200], direc[30] = "./maassdata/";
    FILE *fp;
    struct dirent **namelist;
    int i,n,j,sign;
    char *token;
    size_t last_idx;

    double test;

    long level,levelcount=0,*levels = (long *)malloc(MAXLEVEL * sizeof(long)); //This is too many terms
    for(i = 2; i <= MAXLEVEL; i++){
        if(n_is_squarefree(i) == 1){
            levels[levelcount] = i;
            levelcount++;
        }
    }
    levels = (long *)realloc(levels, levelcount * sizeof(long)); // Realloc to correct number of terms
    long levelindexcounters[levelcount];
    n = scandir(direc, &namelist, 0,descalphasort);
    arb_mat_t Rvalsover[levelcount],Rvals[levelcount];
    for(j = 0; j < levelcount; j++){
        arb_mat_init(Rvalsover[j],n,3);
        levelindexcounters[j]=0;
    }
    if (n < 0)
        perror("scandir");
    else{
        //n = n-2;
        while (n--) {
            strcpy(full,direc);
            strcat(full,namelist[n]->d_name);
            fp = fopen(full,"r");
            // Get R value
            if(fgets(line,sizeof(line),fp) != NULL){
                last_idx = strlen(line) - 1;

                if( line[last_idx] == '\n' ) {
                    line[last_idx] = '\0';
                }
                token = strtok(line, ",");
                arb_set_str(Rval,token,prec);
                while(token != NULL){
                    arb_set_str(eps,token,prec);
                    token = strtok(NULL, ",");
                }
                if(arf_is_nan(arb_midref(eps))){
                    fclose(fp);
                    free(namelist[n]);
                    continue;
                }
                arb_add_error(Rval,eps);
            }
            // Get level
            if(fgets(line,sizeof(line),fp) != NULL){
                sscanf(line,"%ld\n",&level);
            }

            if(fgets(line,sizeof(line),fp) != NULL){
                sscanf(line,"%d\n",&sign);
                if(sign != sign_wanted){
                    fclose(fp);
                    free(namelist[n]);
                    continue;
                }
                for(j=0;j<levelcount;j++){
                    if(levels[j]==level){
                        arb_set(arb_mat_entry(Rvalsover[j],levelindexcounters[j],0),Rval);
                        arb_set(arb_mat_entry(Rvalsover[j],levelindexcounters[j],1),eps);
                        levelindexcounters[j]++;
                        break;
                    }
                }
            }
            fclose(fp);
            free(namelist[n]);
        }
    free(namelist);
    }
    acb_t z1,z2;
    arb_t b1,S,gscale,t1val;
    acb_init(z1); acb_init(z2);
    arb_init(b1); arb_init(S); arb_init(gscale); arb_init(t1val);
    int TFpow = 13;
    arb_set_d(gscale,5.51341248666333);
    arb_div_ui(gscale,gscale,TFpow,prec);
    arb_t tol;
    arb_init(tol);
    arb_set_d(tol,pow(10,-2));

    if(sign_wanted == 0){
        sprintf(full,"./tfdata/completedata_even.txt");
        fp = fopen(full,"w");
    } else if(sign_wanted == 1){
        sprintf(full,"./tfdata/completedata_odd.txt");
        fp = fopen(full,"w");
    }
    

    for(j=0;j<levelcount;j++){
        arb_mat_window_init(Rvals[j],Rvalsover[j],0,0,levelindexcounters[j],2);
        arb_zero(b1);
        arb_zero(S);
        for(i=0;i<arb_mat_nrows(Rvals[j]);i++){
            if(arf_is_nan(arb_midref(arb_mat_entry(Rvals[j],i,1)))){
                continue;
            } else if(arb_lt(tol,arb_mat_entry(Rvals[j],i,1)) || arb_is_zero(arb_mat_entry(Rvals[j],i,1))){
                continue;
            }
            //arb_printd(arb_mat_entry(Rvals[j],i,0),16);
            //printf("\n");
            arb_add(b1, arb_mat_entry(Rvals[j],i,0),arb_mat_entry(Rvals[j],i,1),prec);
            arb_mul(b1,b1,gscale,prec);
            acb_set_arb(z1,b1);
            h_k(z2,z1,TFpow);
            arb_add(S,S,acb_realref(z2),prec);
        }
        readT1value_parity(t1val,levels[j],sign_wanted);

        printf("TF1 = ");
        arb_printd(t1val,16);
        printf("\n");

        printf("level = %ld, S = ",levels[j]);
        arb_printd(S,16);
        printf("\n");

        arb_sub(S,t1val,S,prec);
        printf("diff = ");
        arb_printd(S,16);
        printf("\n");

        acb_set_d_d(z1,20,0);

        acb_printd(z1,16);
        printf("\n");

        hbisect(b1,acb_imagref(z1),acb_realref(z1),S,TFpow,gscale);

        printf("List is complete for all lambda < ");
        arb_printd(b1,16);
        printf("\n");

        test = arf_get_d(arb_midref(b1),ARF_RND_NEAR);
        if(test == 20){
            test = 0;
        }

        fprintf(fp,"%ld,%.30f\n",levels[j],test);
    }
    //arb_mat_printd(Rvals[0],10);
    fclose(fp);
}

int main(int argc, char const *argv[])
{
    completeness_parity(1);
    return 0;
}