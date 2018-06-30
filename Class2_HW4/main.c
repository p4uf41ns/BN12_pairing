#include "bn12.h"

/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
int main(void){
    BN12_init();
    BN12_print_parameters();
    
    //test_Field();
    //test_Frobenius_map();
    //test_skew_frobenius_map();
    //test_rational_point();
    //test_twist();
    //test_plain_ate_pairing();
    //test_opt_ate_pairing();
    //test_x_ate_pairing();
    //test_G1_SCM();
    //test_G2_SCM();
    //test_G3_EXP();
    compare_pairings();
    
    BN12_clear();
    return 0;
}

/*============================================================================*/
/* Field                                                                      */
/*============================================================================*/
/*----------------------------------------------------------------------------*/
//Fp
void Fp_init(Fp *A){
    mpz_init(A->x0);
}

void Fp_clear(Fp *A){
    mpz_clear(A->x0);
}

void Fp_printf(Fp *A,char *str){
    gmp_printf("%s%Zd",str,A->x0);
}

void Fp_set(Fp *ANS,Fp *A){
    mpz_set(ANS->x0,A->x0);
}

void Fp_set_ui(Fp *ANS,unsigned long int UI){
    mpz_set_ui(ANS->x0,UI);
}

void Fp_set_mpz(Fp *ANS,mpz_t A){
    mpz_set(ANS->x0,A);
}

void Fp_set_neg(Fp *ANS,Fp *A){
    //mpz_cost.mpz_add++;
    //Fp_cost.Fp_neg++;
    mpz_sub(ANS->x0,prime,A->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_set_random(Fp *ANS,gmp_randstate_t state){
    mpz_urandomm(ANS->x0,state,prime);
}

void Fp_mul(Fp *ANS,Fp *A,Fp *B){
    /*if(mpz_cmp(A->x0,B->x0)==0){
        mpz_cost.mpz_sqr++;
        Fp_cost.Fp_sqr++;
    }else{
        mpz_cost.mpz_mul++;
        Fp_cost.Fp_mul++;
    }*/
    mpz_mul(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_mul_ui(Fp *ANS,Fp *A,unsigned long int UI){
    //Fp_cost.Fp_mul_ui++;
    //mpz_cost.mpz_mul_ui++;
    mpz_mul_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_mul_mpz(Fp *ANS,Fp *A,mpz_t B){
    /*if(mpz_cmp(A->x0,B)==0){
        mpz_cost.mpz_sqr++;
        Fp_cost.Fp_sqr++;
    }else{
        mpz_cost.mpz_mul++;
        Fp_cost.Fp_mul_mpz++;
    }*/
    mpz_mul(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_add(Fp *ANS,Fp *A,Fp *B){
    //Fp_cost.Fp_add++;
    //mpz_cost.mpz_add++;
    mpz_add(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_add_ui(Fp *ANS,Fp *A,unsigned long int UI){
    //Fp_cost.Fp_add_ui++;
    //mpz_cost.mpz_add_ui++;
    mpz_add_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_add_mpz(Fp *ANS,Fp *A,mpz_t B){
    //Fp_cost.Fp_add_mpz++;
    //mpz_cost.mpz_add++;
    mpz_add(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_sub(Fp *ANS,Fp *A,Fp *B){
    //Fp_cost.Fp_add++;
    //mpz_cost.mpz_add++;
    mpz_sub(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_sub_ui(Fp *ANS,Fp *A,unsigned long int UI){
    //Fp_cost.Fp_add_ui++;
    //mpz_cost.mpz_add_ui++;
    mpz_sub_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_sub_mpz(Fp *ANS,Fp *A,mpz_t B){
    //Fp_cost.Fp_add_mpz++;
    //mpz_cost.mpz_add++;
    mpz_sub(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_inv(Fp *ANS,Fp *A){
    //Fp_cost.Fp_inv++;
    //mpz_cost.mpz_invert++;
    mpz_invert(ANS->x0,A->x0,prime);
}

int  Fp_legendre(Fp *A){
    return mpz_legendre(A->x0,prime);
}

int  Fp_isCNR(Fp *A){
    Fp tmp;
    Fp_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_sub_ui(exp,prime,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp_pow(&tmp,A,exp);
    
    if(Fp_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp_clear(&tmp);
        return -1;
    }
}

void Fp_sqrt(Fp *ANS,Fp *A){
    Fp x,y,t,k,n,tmp;
    Fp_init(&x);
    Fp_init(&y);
    Fp_init(&t);
    Fp_init(&k);
    Fp_init(&n);
    Fp_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    Fp_set_random(&n,state);
    
    while(Fp_legendre(&n)!=-1){
        Fp_set_random(&n,state);
    }
    mpz_sub_ui(q,prime,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp_pow(&x,A,exp);
    Fp_mul(&tmp,&x,&x);
    Fp_mul(&k,&tmp,A);
    Fp_mul(&x,&x,A);
    while(mpz_cmp_ui(k.x0,1)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp_pow(&tmp,&k,exp);
        while(mpz_cmp_ui(tmp.x0,1)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp_pow(&t,&y,result);
        Fp_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp_mul(&x,&x,&t);
        Fp_mul(&k,&k,&y);
    }
    Fp_set(ANS,&x);
    
    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
    Fp_clear(&x);
    Fp_clear(&y);
    Fp_clear(&t);
    Fp_clear(&k);
    Fp_clear(&n);
    Fp_clear(&tmp);
}

void Fp_pow(Fp *ANS,Fp *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp tmp;
    Fp_init(&tmp);
    
    Fp_set(&tmp,A);
    
    for(i=1; binary[i]!='\0'; i++){
        Fp_mul(&tmp,&tmp,&tmp);
        if(binary[i]=='1'){
            Fp_mul(&tmp,A,&tmp);
        }
    }
    Fp_set(ANS,&tmp);
    
    Fp_clear(&tmp);
}


int  Fp_cmp(Fp *A,Fp *B){
    if(mpz_cmp(A->x0,B->x0)==0){
        return 0;   
    }
    return 1;
}

int  Fp_cmp_ui(Fp *A,unsigned long int UI){
    if(mpz_cmp_ui(A->x0,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp_cmp_mpz(Fp *A,mpz_t B){
    if(mpz_cmp(A->x0,B)==0){
        return 0;
    }
    return 1;
}

int  Fp_cmp_zero(Fp *A){
    if(mpz_cmp_ui(A->x0,0)==0){
        return 0;
    }
    return 1;
}

int  Fp_cmp_one(Fp *A){
    if(mpz_cmp_ui(A->x0,1)==0){
        return 0;
    }
    return 1;
}
/*----------------------------------------------------------------------------*/
//Fp2
void Fp2_init(Fp2 *A){
    Fp_init(&A->x0);
    Fp_init(&A->x1);
}

void Fp2_clear(Fp2 *A){
    Fp_clear(&A->x0);
    Fp_clear(&A->x1);
}

void Fp2_printf(Fp2 *A,char *str){
    gmp_printf("%s(",str);
    Fp_printf(&A->x0,"");
    gmp_printf(",");
    Fp_printf(&A->x1,"");
    gmp_printf(")");
}

void Fp2_set(Fp2 *ANS,Fp2 *A){
    Fp_set(&ANS->x0,&A->x0);
    Fp_set(&ANS->x1,&A->x1);
}

void Fp2_set_ui(Fp2 *ANS,unsigned long int UI){
    Fp_set_ui(&ANS->x0,UI);
    Fp_set_ui(&ANS->x1,UI);
}    

void Fp2_set_mpz(Fp2 *ANS,mpz_t A){
    Fp_set_mpz(&ANS->x0,A);
    Fp_set_mpz(&ANS->x1,A);
}

void Fp2_set_neg(Fp2 *ANS,Fp2 *A){
    Fp_set_neg(&ANS->x0,&A->x0);
    Fp_set_neg(&ANS->x1,&A->x1);
}

void Fp2_set_random(Fp2 *ANS,gmp_randstate_t state){
    Fp_set_random(&ANS->x0,state);
    Fp_set_random(&ANS->x1,state);
}

void Fp2_mul(Fp2 *ANS,Fp2 *A,Fp2 *B){
    Fp tmp1,tmp2,tmp3,tmp4;
	Fp_init(&tmp1);
	Fp_init(&tmp2);
	Fp_init(&tmp3);
	Fp_init(&tmp4);
	
	//set
	Fp_mul(&tmp1,&A->x0,&B->x0);//a*c	
	Fp_mul(&tmp2,&A->x1,&B->x1);//b*d
	Fp_add(&tmp3,&A->x0,&A->x1);//a+b
	Fp_add(&tmp4,&B->x0,&B->x1);//c+d
	//x0
	Fp_sub(&ANS->x0,&tmp1,&tmp2);//a*c+b*d*v
	//x1
	Fp_mul(&ANS->x1,&tmp3,&tmp4);//(a+b)(c+d)
	Fp_sub(&ANS->x1,&ANS->x1,&tmp1);
	Fp_sub(&ANS->x1,&ANS->x1,&tmp2);
	
	//clear
	Fp_clear(&tmp1);
	Fp_clear(&tmp2);
	Fp_clear(&tmp3);
	Fp_clear(&tmp4);
}

void Fp2_mul_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI){
    Fp_mul_ui(&ANS->x0,&A->x0,UI);
    Fp_mul_ui(&ANS->x1,&A->x1,UI);
}

void Fp2_mul_mpz(Fp2 *ANS,Fp2 *A,mpz_t B){
    Fp_mul_mpz(&ANS->x0,&A->x0,B);
    Fp_mul_mpz(&ANS->x1,&A->x1,B);
}

void Fp2_mul_basis(Fp2 *ANS,Fp2 *A){
    Fp tmp;
    Fp_init(&tmp);
    Fp_set(&tmp,&A->x0);
    
    Fp_sub(&ANS->x0,&tmp,&A->x1);
    Fp_add(&ANS->x1,&tmp,&A->x1);
    
    Fp_clear(&tmp);
}

void Fp2_inv_basis(Fp2 *ANS,Fp2 *A){
    Fp2 tmp;
    Fp2_init(&tmp);
    Fp2_set(&tmp,A);
    
    Fp_add(&ANS->x0,&tmp.x0,&tmp.x1);
    Fp_mul_mpz(&ANS->x0,&ANS->x0,Alpha_1_inv.x0.x0);
    Fp_sub(&ANS->x1,&tmp.x1,&tmp.x0);
    Fp_mul_mpz(&ANS->x1,&ANS->x1,Alpha_1_inv.x0.x0);
    
    Fp2_clear(&tmp);
}

void Fp2_sqr(Fp2 *ANS,Fp2 *A){
    Fp tmp1,tmp2;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    
    Fp_add(&tmp1,&A->x0,&A->x1);
	Fp_sub(&tmp2,&A->x0,&A->x1);
	//x1
	Fp_mul(&ANS->x1,&A->x0,&A->x1);
	Fp_add(&ANS->x1,&ANS->x1,&ANS->x1);
	//x0
	Fp_mul(&ANS->x0,&tmp1,&tmp2);
    
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
}

void Fp2_add(Fp2 *ANS,Fp2 *A,Fp2 *B){
    Fp_add(&ANS->x0,&A->x0,&B->x0);
    Fp_add(&ANS->x1,&A->x1,&B->x1);
}

void Fp2_add_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI){
    Fp_add_ui(&ANS->x0,&A->x0,UI);
    Fp_add_ui(&ANS->x1,&A->x1,UI);
}

void Fp2_add_mpz(Fp2 *ANS,Fp2 *A,mpz_t B){
    Fp_add_mpz(&ANS->x0,&A->x0,B);
    Fp_add_mpz(&ANS->x1,&A->x1,B);
}

void Fp2_sub(Fp2 *ANS,Fp2 *A,Fp2 *B){
    Fp_sub(&ANS->x0,&A->x0,&B->x0);
    Fp_sub(&ANS->x1,&A->x1,&B->x1);
}

void Fp2_sub_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI){
    Fp_sub_ui(&ANS->x0,&A->x0,UI);
    Fp_sub_ui(&ANS->x1,&A->x1,UI);
}

void Fp2_sub_mpz(Fp2 *ANS,Fp2 *A,mpz_t B){
    Fp_sub_mpz(&ANS->x0,&A->x0,B);
    Fp_sub_mpz(&ANS->x1,&A->x1,B);
}

void Fp2_inv(Fp2 *ANS,Fp2 *A){
    Fp c_x0,c_x1,t0,t1;
    Fp_init(&c_x0);
    Fp_init(&c_x1);
    Fp_init(&t0);
    Fp_init(&t1);
    
    Fp_set(&c_x0,&A->x0);
    Fp_set_neg(&c_x1,&A->x1);
    
    Fp_mul(&t0,&c_x0,&A->x0);
    Fp_mul(&t1,&c_x1,&A->x1);
    Fp_sub(&t0,&t0,&t1);
    Fp_inv(&t0,&t0);
    Fp_mul(&ANS->x0,&c_x0,&t0);
    Fp_mul(&ANS->x1,&c_x1,&t0);
    
    Fp_clear(&c_x0);
    Fp_clear(&c_x1);
    Fp_clear(&t0);
    Fp_clear(&t1);
}

int  Fp2_legendre(Fp2 *A){
    Fp2 tmp;
    Fp2_init(&tmp);
    
    mpz_t exp;
    mpz_init(exp);
    mpz_pow_ui(exp,prime,2);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp,A,exp);
    
    if(Fp2_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp2_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp2_clear(&tmp);
        return -1;
    }
}

int  Fp2_isCNR(Fp2 *A){
    Fp2 tmp;
    Fp2_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_pow_ui(exp,prime,2);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp2_pow(&tmp,A,exp);
    
    if(Fp2_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp2_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp2_clear(&tmp);
        return -1;
    }

}

void Fp2_sqrt(Fp2 *ANS,Fp2 *A){
    Fp2 x,y,t,k,n,tmp;
    Fp2_init(&x);
    Fp2_init(&y);
    Fp2_init(&t);
    Fp2_init(&k);
    Fp2_init(&n);
    Fp2_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    Fp2_set_random(&n,state);
    while(Fp2_legendre(&n)!=-1){
        Fp2_set_random(&n,state);
    }
    mpz_pow_ui(q,prime,2);
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp2_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&x,A,exp);
    Fp2_mul(&tmp,&x,&x);
    Fp2_mul(&k,&tmp,A);
    Fp2_mul(&x,&x,A);
    while(Fp2_cmp_one(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp2_pow(&tmp,&k,exp);
        while(Fp2_cmp_one(&tmp)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp2_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp2_pow(&t,&y,result);
        Fp2_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp2_mul(&x,&x,&t);
        Fp2_mul(&k,&k,&y);
    }
    Fp2_set(ANS,&x);
    
    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
    Fp2_clear(&x);
    Fp2_clear(&y);
    Fp2_clear(&t);
    Fp2_clear(&k);
    Fp2_clear(&n);
    Fp2_clear(&tmp);
}

void Fp2_pow(Fp2 *ANS,Fp2 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp2 tmp;
    Fp2_init(&tmp);
    Fp2_set(&tmp,A);
    
    for(i=1; binary[i]!='\0'; i++){
        Fp2_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp2_mul(&tmp,A,&tmp);
        }
    }
    
    Fp2_set(ANS,&tmp);
    Fp2_clear(&tmp);
}


int  Fp2_cmp(Fp2 *A,Fp2 *B){
    if(Fp_cmp(&A->x0,&B->x0)==0 && Fp_cmp(&A->x1,&B->x1)==0){
        return 0;   
    }
    return 1;
}

int  Fp2_cmp_ui(Fp2 *A,unsigned long int UI){
    if(Fp_cmp_ui(&A->x0,UI)==0 && Fp_cmp_ui(&A->x1,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp2_cmp_mpz(Fp2 *A,mpz_t B){
    if(Fp_cmp_mpz(&A->x0,B)==0 && Fp_cmp_mpz(&A->x1,B)==0){
        return 0;
    }
    return 1;
}

int  Fp2_cmp_zero(Fp2 *A){
    if(Fp_cmp_zero(&A->x0)==0 && Fp_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

int  Fp2_cmp_one(Fp2 *A){
    if(Fp_cmp_one(&A->x0)==0 && Fp_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}
/*----------------------------------------------------------------------------*/
//Fp6
void Fp6_init(Fp6 *A){
    Fp2_init(&A->x0);
    Fp2_init(&A->x1);
    Fp2_init(&A->x2);
}

void Fp6_clear(Fp6 *A){
    Fp2_clear(&A->x0);
    Fp2_clear(&A->x1);
    Fp2_clear(&A->x2);
}

void Fp6_printf(Fp6 *A,char *str){
    gmp_printf("%s(",str);
    Fp2_printf(&A->x0,"");
    gmp_printf(",");
    Fp2_printf(&A->x1,"");
    gmp_printf(",");
    Fp2_printf(&A->x2,"");
    gmp_printf(")");
}

void Fp6_set(Fp6 *ANS,Fp6 *A){
    Fp2_set(&ANS->x0,&A->x0);
    Fp2_set(&ANS->x1,&A->x1);
    Fp2_set(&ANS->x2,&A->x2);
}

void Fp6_set_ui(Fp6 *ANS,unsigned long int UI){
    Fp2_set_ui(&ANS->x0,UI);
    Fp2_set_ui(&ANS->x1,UI);
    Fp2_set_ui(&ANS->x2,UI);
}

void Fp6_set_mpz(Fp6 *ANS,mpz_t A){
    Fp2_set_mpz(&ANS->x0,A);
    Fp2_set_mpz(&ANS->x1,A);
    Fp2_set_mpz(&ANS->x2,A);
}

void Fp6_set_neg(Fp6 *ANS,Fp6 *A){
    Fp2_set_neg(&ANS->x0,&A->x0);
    Fp2_set_neg(&ANS->x1,&A->x1);
    Fp2_set_neg(&ANS->x2,&A->x2);
}

void Fp6_set_random(Fp6 *ANS,gmp_randstate_t state){
    Fp2_set_random(&ANS->x0,state);
    Fp2_set_random(&ANS->x1,state);
    Fp2_set_random(&ANS->x2,state);
}

void Fp6_mul(Fp6 *ANS,Fp6 *A,Fp6 *B){
    Fp2 tmp00,tmp11,tmp22,tmp,t0,t1,t2;
    Fp2_init(&tmp00);
    Fp2_init(&tmp11);
    Fp2_init(&tmp22);
    Fp2_init(&tmp);
    Fp2_init(&t0);
    Fp2_init(&t1);
    Fp2_init(&t2);
    
    //set
    Fp2_mul(&tmp00,&A->x0,&B->x0);//x0*y0
    Fp2_mul(&tmp11,&A->x1,&B->x1);//x1*y1
    Fp2_mul(&tmp22,&A->x2,&B->x2);//x2*y2
    
    Fp2_add(&t0,&A->x0,&A->x1);//x0+x1
    Fp2_add(&tmp,&B->x0,&B->x1);//y0+y1
    Fp2_mul(&t0,&t0,&tmp);//(x0+x1)(y0+y1)
    
    Fp2_add(&t1,&A->x1,&A->x2);//x1+x2
    Fp2_add(&tmp,&B->x1,&B->x2);//y1+y2
    Fp2_mul(&t1,&t1,&tmp);//(x1+x2)(y1+y2)
    
    Fp2_add(&t2,&B->x0,&B->x2);//y2+y0
    Fp2_add(&tmp,&A->x0,&A->x2);//x2+x0
    Fp2_mul(&t2,&t2,&tmp);//(x2+x0)(y2+y0)
    //x0
    Fp2_sub(&t1,&t1,&tmp11);
    Fp2_sub(&t1,&t1,&tmp22);//(x1+x2)(y1+y2)-x1y1-x2y2
    Fp2_mul_basis(&tmp,&t1);
    Fp2_add(&ANS->x0,&tmp00,&tmp);
    //x1
    Fp2_sub(&t0,&t0,&tmp00);
    Fp2_sub(&t0,&t0,&tmp11);
    Fp2_mul_basis(&tmp,&tmp22);
    Fp2_add(&ANS->x1,&tmp,&t0);
    //x2
    Fp2_sub(&t2,&t2,&tmp00);
    Fp2_sub(&t2,&t2,&tmp22);
    Fp2_add(&ANS->x2,&tmp11,&t2);
    
    //clear
    Fp2_clear(&tmp00);
    Fp2_clear(&tmp11);
    Fp2_clear(&tmp22);
    Fp2_clear(&tmp);
    Fp2_clear(&t0);
    Fp2_clear(&t1);
    Fp2_clear(&t2);
}

void Fp6_mul_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI){
    Fp2_mul_ui(&ANS->x0,&A->x0,UI);
    Fp2_mul_ui(&ANS->x1,&A->x1,UI);
    Fp2_mul_ui(&ANS->x2,&A->x2,UI);
}

void Fp6_mul_mpz(Fp6 *ANS,Fp6 *A,mpz_t B){
    Fp2_mul_mpz(&ANS->x0,&A->x0,B);
    Fp2_mul_mpz(&ANS->x1,&A->x1,B);
    Fp2_mul_mpz(&ANS->x2,&A->x2,B);
}

void Fp6_mul_basis(Fp6 *ANS,Fp6 *A){
    Fp6 tmp;
	Fp6_init(&tmp);
	Fp6_set(&tmp,A);
	
	Fp_sub(&ANS->x0.x0,&tmp.x2.x0,&tmp.x2.x1);
	Fp_add(&ANS->x0.x1,&tmp.x2.x0,&tmp.x2.x1);
	Fp_set(&ANS->x1.x0,&tmp.x0.x0);
	Fp_set(&ANS->x1.x1,&tmp.x0.x1);
	Fp_set(&ANS->x2.x0,&tmp.x1.x0);
	Fp_set(&ANS->x2.x1,&tmp.x1.x1);
	
	Fp6_clear(&tmp);
}

void Fp6_sqr(Fp6 *ANS,Fp6 *A){
    Fp2 tmp00,tmp12_2,tmp01_2,tmp22,tmp;
    Fp2_init(&tmp00);
    Fp2_init(&tmp22);
    Fp2_init(&tmp12_2);
    Fp2_init(&tmp01_2);
    Fp2_init(&tmp);
    
    Fp2_sqr(&tmp00,&A->x0);        //x0^2
    Fp2_sqr(&tmp22,&A->x2);        //x2^2
    Fp2_add(&tmp,&A->x1,&A->x1);        //2x1
    Fp2_mul(&tmp12_2,&tmp,&A->x2);  //2x1x2
    Fp2_mul(&tmp01_2,&A->x0,&tmp);  //2x0x1
    Fp2_add(&tmp,&A->x0,&A->x1);        //x0+x1+x2
    Fp2_add(&tmp,&tmp,&A->x2);
    
    //x0
    Fp2_mul_basis(&ANS->x0,&tmp12_2);
    Fp2_add(&ANS->x0,&ANS->x0,&tmp00);
    //x1
    Fp2_mul_basis(&ANS->x1,&tmp22);
    Fp2_add(&ANS->x1,&ANS->x1,&tmp01_2);
    //x2
    Fp2_sqr(&ANS->x2,&tmp);
    Fp2_add(&tmp,&tmp00,&tmp22);
    Fp2_add(&tmp,&tmp,&tmp12_2);
    Fp2_add(&tmp,&tmp,&tmp01_2);
    Fp2_sub(&ANS->x2,&ANS->x2,&tmp);
    
    Fp2_clear(&tmp00);
    Fp2_clear(&tmp22);
    Fp2_clear(&tmp12_2);
    Fp2_clear(&tmp01_2);
    Fp2_clear(&tmp);
}

void Fp6_add(Fp6 *ANS,Fp6 *A,Fp6 *B){
    Fp2_add(&ANS->x0,&A->x0,&B->x0);
    Fp2_add(&ANS->x1,&A->x1,&B->x1);
    Fp2_add(&ANS->x2,&A->x2,&B->x2);
}

void Fp6_add_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI){
    Fp2_add_ui(&ANS->x0,&A->x0,UI);
    Fp2_add_ui(&ANS->x1,&A->x1,UI);
    Fp2_add_ui(&ANS->x2,&A->x2,UI);
}

void Fp6_add_mpz(Fp6 *ANS,Fp6 *A,mpz_t B){
    Fp2_add_mpz(&ANS->x0,&A->x0,B);
    Fp2_add_mpz(&ANS->x1,&A->x1,B);
    Fp2_add_mpz(&ANS->x2,&A->x2,B);
}

void Fp6_sub(Fp6 *ANS,Fp6 *A,Fp6 *B){
    Fp2_sub(&ANS->x0,&A->x0,&B->x0);
    Fp2_sub(&ANS->x1,&A->x1,&B->x1);
    Fp2_sub(&ANS->x2,&A->x2,&B->x2);
}

void Fp6_sub_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI){
    Fp2_sub_ui(&ANS->x0,&A->x0,UI);
    Fp2_sub_ui(&ANS->x1,&A->x1,UI);
    Fp2_sub_ui(&ANS->x2,&A->x2,UI);
}

void Fp6_sub_mpz(Fp6 *ANS,Fp6 *A,mpz_t B){
    Fp2_sub_mpz(&ANS->x0,&A->x0,B);
    Fp2_sub_mpz(&ANS->x1,&A->x1,B);
    Fp2_sub_mpz(&ANS->x2,&A->x2,B);
}

void Fp6_inv(Fp6 *ANS,Fp6 *A){
    Fp6 barA;
    Fp6_init(&barA);
    Fp2 s,t0,t1,t2,t3;
    Fp2_init(&s);
    Fp2_init(&t0);
    Fp2_init(&t1);
    Fp2_init(&t2);
    Fp2_init(&t3);
    
    Fp2_sqr(&t0,&A->x0); //t0=a0^2
    Fp2_sqr(&t1,&A->x1); //t1=a1^2
    Fp2_sqr(&t2,&A->x2); //t2=a2^2
    
    Fp2_mul(&t3,&A->x1,&A->x2);
    Fp2_mul_basis(&t3,&t3);
    Fp2_sub(&barA.x0,&t0,&t3);
    
    Fp2_mul(&t3,&A->x0,&A->x1);
    Fp2_mul_basis(&barA.x1,&t2);
    Fp2_sub(&barA.x1,&barA.x1,&t3);
    
    Fp2_mul(&t3,&A->x0,&A->x2);
    Fp2_sub(&barA.x2,&t1,&t3);
    
    Fp2_mul(&t0,&t0,&A->x0);
    Fp2_mul(&t2,&t2,&A->x2);
    Fp2_mul_basis(&t2,&t2);
    
    Fp2_add(&s,&t3,&t3);
    Fp2_add(&s,&s,&t3);
    Fp2_sub(&s,&t1,&s);
    Fp2_mul(&s,&s,&A->x1);
    Fp2_add(&s,&s,&t2);
    Fp2_mul_basis(&s,&s);
    Fp2_add(&s,&s,&t0);
    
    Fp2_inv(&s,&s);
    
    Fp2_mul(&ANS->x0,&barA.x0,&s);
    Fp2_mul(&ANS->x1,&barA.x1,&s);
    Fp2_mul(&ANS->x2,&barA.x2,&s);
    
    Fp6_clear(&barA);
    Fp2_clear(&s);
    Fp2_clear(&t0);
    Fp2_clear(&t1);
    Fp2_clear(&t2);
    Fp2_clear(&t3);
}

int  Fp6_legendre(Fp6 *A){
    mpz_t exp;
    mpz_init(exp);
    Fp6 tmp;
    Fp6_init(&tmp);
    
    mpz_pow_ui(exp,prime,6);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp6_pow(&tmp,A,exp);
    
    if(Fp6_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp6_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp6_clear(&tmp);
        return -1;
    }
}

int  Fp6_isCNR(Fp6 *A){
    Fp6 tmp;
    Fp6_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_pow_ui(exp,prime,6);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp6_pow(&tmp,A,exp);
    
    if(Fp6_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp6_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp6_clear(&tmp);
        return -1;
    }
}

void Fp6_sqrt(Fp6 *ANS,Fp6 *A){
    Fp6 tmp1,tmp2;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    mpz_t exp,buf;
    mpz_init(exp);
    mpz_init(buf);
    
    Fp2_set(&tmp1.x0,&A->x0);
    Fp2_mul_mpz(&tmp1.x1,&A->x1,frobenius_constant[f_p4][1].x0.x0);
    Fp2_mul_mpz(&tmp1.x2,&A->x2,frobenius_constant[f_p4][2].x0.x0);
    
    Fp2_set(&tmp2.x0,&A->x0);
    Fp2_mul_mpz(&tmp2.x1,&A->x1,frobenius_constant[f_p2][1].x0.x0);
    Fp2_mul_mpz(&tmp2.x2,&A->x2,frobenius_constant[f_p2][2].x0.x0);
    
    Fp6_mul(&tmp1,&tmp1,&tmp2);
    Fp6_mul(&tmp1,&tmp1,A);
    Fp6_set_ui(&tmp2,0);
    Fp2_sqrt(&tmp2.x0,&tmp1.x0);
    Fp2_inv(&tmp2.x0,&tmp2.x0);
    Fp2_set(&tmp2.x0,&tmp2.x0);
    mpz_pow_ui(exp,prime,8);
    mpz_pow_ui(buf,prime,4);
    mpz_add(exp,exp,buf);
    mpz_add_ui(exp,exp,2);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp6_pow(&tmp1,A,exp);
    Fp6_mul(&tmp1,&tmp1,&tmp2);
    Fp6_set(ANS,&tmp1);
    
    mpz_clear(exp);
    mpz_clear(buf);
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
}

void Fp6_pow(Fp6 *ANS,Fp6 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp6 tmp;
    Fp6_init(&tmp);
    Fp6_set(&tmp,A);
    
    for(i=1; binary[i]!='\0'; i++){
        Fp6_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp6_mul(&tmp,A,&tmp);
        }
    }
    
    Fp6_set(ANS,&tmp);
    Fp6_clear(&tmp);
}

int  Fp6_cmp(Fp6 *A,Fp6 *B){
    if(Fp2_cmp(&A->x0,&B->x0)==0 && Fp2_cmp(&A->x1,&B->x1)==0 && Fp2_cmp(&A->x2,&B->x2)==0){
        return 0;   
    }
    return 1;
}

int  Fp6_cmp_ui(Fp6 *A,unsigned long int UI){
    if(Fp2_cmp_ui(&A->x0,UI)==0 && Fp2_cmp_ui(&A->x1,UI)==0 && Fp2_cmp_ui(&A->x2,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp6_cmp_mpz(Fp6 *A,mpz_t B){
    if(Fp2_cmp_mpz(&A->x0,B)==0 && Fp2_cmp_mpz(&A->x1,B)==0 && Fp2_cmp_mpz(&A->x2,B)==0){
        return 0;
    }
    return 1;
}

int  Fp6_cmp_zero(Fp6 *A){
    if(Fp2_cmp_zero(&A->x0)==0 && Fp2_cmp_zero(&A->x1)==0 && Fp2_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}

int  Fp6_cmp_one(Fp6 *A){
    if(Fp2_cmp_one(&A->x0)==0 && Fp2_cmp_zero(&A->x1)==0 && Fp2_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}
/*----------------------------------------------------------------------------*/
//Fp12
void Fp12_init(Fp12 *A){
    Fp6_init(&A->x0);
    Fp6_init(&A->x1);
}

void Fp12_clear(Fp12 *A){
    Fp6_clear(&A->x0);
    Fp6_clear(&A->x1);
}

void Fp12_printf(Fp12 *A,char *str){
    gmp_printf("%s(",str);
    Fp6_printf(&A->x0,"");
    gmp_printf(",");
    Fp6_printf(&A->x1,"");
    gmp_printf(")");
}

void Fp12_set(Fp12 *ANS,Fp12 *A){
    Fp6_set(&ANS->x0,&A->x0);
    Fp6_set(&ANS->x1,&A->x1);
}

void Fp12_set_ui(Fp12 *ANS,unsigned long int UI){
    Fp6_set_ui(&ANS->x0,UI);
    Fp6_set_ui(&ANS->x1,UI);
}

void Fp12_set_mpz(Fp12 *ANS,mpz_t A){
    Fp6_set_mpz(&ANS->x0,A);
    Fp6_set_mpz(&ANS->x1,A);    
}

void Fp12_set_neg(Fp12 *ANS,Fp12 *A){
    Fp6_set_neg(&ANS->x0,&A->x0);
    Fp6_set_neg(&ANS->x1,&A->x1);
}

void Fp12_set_random(Fp12 *ANS,gmp_randstate_t state){
    Fp6_set_random(&ANS->x0,state);
    Fp6_set_random(&ANS->x1,state);
}

void Fp12_mul(Fp12 *ANS,Fp12 *A,Fp12 *B){
    Fp6 tmp1,tmp2;
	Fp6_init(&tmp1);
	Fp6_init(&tmp2);
	
	//set
	Fp6_mul(&tmp2,&A->x1,&B->x1);//b*d
	Fp6_add(&tmp1,&A->x0,&A->x1);//a+b
	Fp6_add(&ANS->x1,&B->x0,&B->x1);//c+d
	Fp6_mul(&ANS->x1,&tmp1,&ANS->x1);//(a+b)(c+d)
	Fp6_mul(&tmp1,&A->x0,&B->x0);//a*c
	//x0
	Fp6_mul_basis(&ANS->x0,&tmp2);//b*d*v
	Fp6_add(&ANS->x0,&ANS->x0,&tmp1);//a*c+b*d*v
	//x1
	Fp6_sub(&ANS->x1,&ANS->x1,&tmp1);
	Fp6_sub(&ANS->x1,&ANS->x1,&tmp2);
	
	//clear
	Fp6_clear(&tmp1);
	Fp6_clear(&tmp2);
}

void Fp12_mul_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_mul_ui(&ANS->x0,&A->x0,UI);
    Fp6_mul_ui(&ANS->x1,&A->x1,UI);
}

void Fp12_mul_mpz(Fp12 *ANS,Fp12 *A,mpz_t B){
    Fp6_mul_mpz(&ANS->x0,&A->x0,B);
    Fp6_mul_mpz(&ANS->x1,&A->x1,B);
}

void Fp12_sqr(Fp12 *ANS,Fp12 *A){
    Fp6 tmp1,tmp2,tmp3;
	Fp6_init(&tmp1);
	Fp6_init(&tmp2);
	Fp6_init(&tmp3);
	
	Fp6_add(&tmp1,&A->x0,&A->x1);
	Fp6_mul_basis(&tmp2,&A->x1);
	Fp6_add(&tmp2,&tmp2,&A->x0);
	Fp6_mul(&tmp3,&A->x0,&A->x1);
	
	//x0
	Fp6_mul(&ANS->x0,&tmp1,&tmp2);
	Fp6_sub(&ANS->x0,&ANS->x0,&tmp3);
	Fp6_mul_basis(&tmp1,&tmp3);
	Fp6_sub(&ANS->x0,&ANS->x0,&tmp1);
	//x1
	Fp6_add(&ANS->x1,&tmp3,&tmp3);
	
	Fp6_clear(&tmp1);
	Fp6_clear(&tmp2);
	Fp6_clear(&tmp3);
}

void Fp12_sqr_cyclotomic(Fp12 *ANS,Fp12 *A){
    //A=a+b*gamma in G3
    //A^2=(1+2b^2*beta)+((a+b)^2-1-b^2*beta-b^2)
    
    Fp6 tmp1,tmp2;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    
    Fp6_add(&tmp1,&A->x0,&A->x1);
    Fp6_sqr(&tmp1,&tmp1);           //(a+b)^2
    Fp6_sqr(&tmp2,&A->x1);          //b^2
    Fp6_mul_basis(&ANS->x1,&tmp2);  //b^2*beta
    Fp6_add(&ANS->x0,&ANS->x1,&ANS->x1);
    Fp_add_ui(&ANS->x0.x0.x0,&ANS->x0.x0.x0,1);   //1+2b^2*beta
    
    Fp6_sub(&ANS->x1,&tmp1,&ANS->x1);
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp2);
    Fp_sub_ui(&ANS->x1.x0.x0,&ANS->x1.x0.x0,1);   //(a+b)^2-1-b^2*beta-b^2
    
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
}

void Fp12_add(Fp12 *ANS,Fp12 *A,Fp12 *B){
    Fp6_add(&ANS->x0,&A->x0,&B->x0);
    Fp6_add(&ANS->x1,&A->x1,&B->x1);
}

void Fp12_add_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_add_ui(&ANS->x0,&A->x0,UI);
    Fp6_add_ui(&ANS->x1,&A->x1,UI);
}

void Fp12_add_mpz(Fp12 *ANS,Fp12 *A,mpz_t B){
    Fp6_add_mpz(&ANS->x0,&ANS->x0,B);
    Fp6_add_mpz(&ANS->x1,&ANS->x1,B);
}

void Fp12_sub(Fp12 *ANS,Fp12 *A,Fp12 *B){
    Fp6_sub(&ANS->x0,&A->x0,&B->x0);
    Fp6_sub(&ANS->x1,&A->x1,&B->x1);
}

void Fp12_sub_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_sub_ui(&ANS->x0,&ANS->x0,UI);
    Fp6_sub_ui(&ANS->x1,&ANS->x1,UI);
}

void Fp12_sub_mpz(Fp12 *ANS,Fp12 *A,mpz_t B){
    Fp6_sub_mpz(&ANS->x0,&ANS->x0,B);
    Fp6_sub_mpz(&ANS->x1,&ANS->x1,B);
}

void Fp12_inv(Fp12 *ANS,Fp12 *A){
    Fp6 c_x0,c_x1,t0,t1;
    Fp6_init(&c_x0);
    Fp6_init(&c_x1);
    Fp6_init(&t0);
    Fp6_init(&t1);
    
    Fp6_set(&c_x0,&A->x0);
    Fp6_set_neg(&c_x1,&A->x1);
    
    Fp6_mul(&t0,&c_x0,&A->x0);
    Fp6_mul(&t1,&c_x1,&A->x1);
    Fp6_mul_basis(&t1,&t1);
    Fp6_add(&t0,&t0,&t1);
    Fp6_inv(&t0,&t0);
    Fp6_mul(&ANS->x0,&c_x0,&t0);
    Fp6_mul(&ANS->x1,&c_x1,&t0);
    
    Fp6_clear(&c_x0);
    Fp6_clear(&c_x1);
    Fp6_clear(&t0);
    Fp6_clear(&t1);
}

int  Fp12_legendre(Fp12 *A){
    mpz_t exp;
    mpz_init(exp);
    Fp12 tmp;
    Fp12_init(&tmp);
    
    mpz_pow_ui(exp,prime,12);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp12_pow(&tmp,A,exp);
    
    if(Fp12_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp12_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp12_clear(&tmp);
        return -1;
    }
}

void Fp12_sqrt(Fp12 *ANS,Fp12 *A){
    Fp12 x,y,t,k,n,tmp;
    Fp12_init(&x);
    Fp12_init(&y);
    Fp12_init(&t);
    Fp12_init(&k);
    Fp12_init(&n);
    Fp12_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
	gmp_randinit_default (state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    Fp12_set_random(&n,state);
    while(Fp12_legendre(&n)!=-1){
        Fp12_set_random(&n,state);
    }
    mpz_pow_ui(q,prime,12);
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp12_pow(&y,&n,q);
    mpz_set_ui(z,e);    
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp12_pow(&x,A,exp);
    Fp12_mul(&tmp,&x,&x);
    Fp12_mul(&k,&tmp,A);
    Fp12_mul(&x,&x,A);
    while(Fp12_cmp_one(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp12_pow(&tmp,&k,exp);
        while(Fp12_cmp_one(&tmp)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp12_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp12_pow(&t,&y,result);
        Fp12_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp12_mul(&x,&x,&t); 
        Fp12_mul(&k,&k,&y);
    }
    Fp12_set(ANS,&x);
    
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(exp);
    mpz_clear(result);
    Fp12_clear(&x);
    Fp12_clear(&y);
    Fp12_clear(&t);
    Fp12_clear(&k);
    Fp12_clear(&n);
    Fp12_clear(&tmp);
}

void Fp12_pow(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp12 tmp;
    Fp12_init(&tmp);
    Fp12_set(&tmp,A);
    
    for(i=1; binary[i]!='\0'; i++){
        Fp12_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp12_mul(&tmp,A,&tmp);
        }
    }
    
    Fp12_set(ANS,&tmp);
    Fp12_clear(&tmp);
}

int  Fp12_cmp(Fp12 *A,Fp12 *B){
    if(Fp6_cmp(&A->x0,&B->x0)==0 && Fp6_cmp(&A->x1,&B->x1)==0){
        return 0;   
    }
    return 1;
}

int  Fp12_cmp_ui(Fp12 *A,unsigned long int UI){
    if(Fp6_cmp_ui(&A->x0,UI)==0 && Fp6_cmp_ui(&A->x1,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp12_cmp_mpz(Fp12 *A,mpz_t B){
    if(Fp6_cmp_mpz(&A->x0,B)==0 && Fp6_cmp_mpz(&A->x1,B)==0){
        return 0;
    }
    return 1;
}

int  Fp12_cmp_zero(Fp12 *A){
    if(Fp6_cmp_zero(&A->x0)==0 && Fp6_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

int  Fp12_cmp_one(Fp12 *A){
    if(Fp6_cmp_one(&A->x0)==0 && Fp6_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

void Fp12_frobenius_map_p1(Fp12 *ANS,Fp12 *A){
	Fp tmp;
    Fp_init(&tmp);
    
    //x0
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    Fp_set(&tmp,&A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1,&tmp);
    Fp2_mul_mpz(&ANS->x0.x1,&ANS->x0.x1,frobenius_constant[f_p1][1].x1.x0);
    Fp_set(&ANS->x0.x2.x0,&A->x0.x2.x0);
    Fp_set_neg(&ANS->x0.x2.x1,&A->x0.x2.x1);
    Fp2_mul_mpz(&ANS->x0.x2,&ANS->x0.x2,frobenius_constant[f_p1][2].x0.x0);
    //x1
    Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    Fp2_mul(&ANS->x1.x0,&ANS->x1.x0,&frobenius_constant[f_p1][3]);
    Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    Fp2_mul_mpz(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant[f_p1][4].x0.x0);
    Fp_add(&tmp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_set(&ANS->x1.x1.x1,&tmp);
    
    Fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    Fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    Fp2_mul(&ANS->x1.x2,&ANS->x1.x2,&frobenius_constant[f_p1][5]);
    
    Fp_clear(&tmp);
}

void Fp12_frobenius_map_p2(Fp12 *ANS,Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpz(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p2][1].x0.x0);
    Fp2_mul_mpz(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p2][2].x0.x0);
    //x1
    Fp2_mul_mpz(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p2][3].x0.x0);
    Fp2_set_neg(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpz(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p2][5].x0.x0);
}

void Fp12_frobenius_map_p3(Fp12 *ANS,Fp12 *A){
    Fp tmp;
    Fp_init(&tmp);
    
    //x0
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    Fp_set(&tmp,&A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1,&tmp);
    Fp_set_neg(&ANS->x0.x2.x0,&A->x0.x2.x0);
    Fp_set(&ANS->x0.x2.x1,&A->x0.x2.x1);
    //x1
    Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    Fp2_mul(&ANS->x1.x0,&ANS->x1.x0,&frobenius_constant[f_p3][3]);
    Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    Fp2_mul_mpz(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant[f_p3][4].x0.x0);
    Fp_add(&tmp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_set(&ANS->x1.x1.x1,&tmp);
    Fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    Fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    Fp2_mul(&ANS->x1.x2,&ANS->x1.x2,&frobenius_constant[f_p3][5]);
    
    Fp_clear(&tmp); 
}

void Fp12_frobenius_map_p4(Fp12 *ANS,Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpz(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p4][1].x0.x0);
    Fp2_mul_mpz(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p4][2].x0.x0);
    //x1
    Fp2_mul_mpz(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p4][3].x0.x0);
    Fp2_set(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpz(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p4][5].x0.x0);
}

void Fp12_frobenius_map_p6(Fp12 *ANS,Fp12 *A){
    //x0
	Fp6_set(&ANS->x0,&A->x0);
	//x1
	Fp6_set_neg(&ANS->x1,&A->x1);
}

void Fp12_frobenius_map_p8(Fp12 *ANS,Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpz(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p8][1].x0.x0);
    Fp2_mul_mpz(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p8][2].x0.x0);
    //x1
    Fp2_mul_mpz(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p8][3].x0.x0);
    Fp2_set(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpz(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p8][5].x0.x0);
}

void Fp12_frobenius_map_p10(Fp12 *ANS,Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpz(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p10][1].x0.x0);
    Fp2_mul_mpz(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p10][2].x0.x0);
    //x1
    Fp2_mul_mpz(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p10][3].x0.x0);
    Fp2_set_neg(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpz(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p10][5].x0.x0);
}

/*============================================================================*/
/* Elliptic curve                                                             */
/*============================================================================*/
/*----------------------------------------------------------------------------*/
//EFp
void EFp_init(EFp *P){
    Fp_init(&P->x);
    Fp_init(&P->y);
    P->infinity=0;
}

void EFp_set(EFp *ANS,EFp *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp_set_ui(EFp *ANS,unsigned long int UI){
    Fp_set_ui(&ANS->x,UI);
    Fp_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp_set_mpz(EFp *ANS,mpz_t A){
    Fp_set_mpz(&ANS->x,A);
    Fp_set_mpz(&ANS->y,A);
    ANS->infinity=0;
}

void EFp_set_neg(EFp *ANS,EFp *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp_clear(EFp *P){
    Fp_clear(&P->x);
    Fp_clear(&P->y);
}

void EFp_printf(EFp *P,char *str){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp_printf(&P->x,"");
        printf(",");
        Fp_printf(&P->y,"");
        printf(")");
    }else{
        printf("0");
    }
}

void EFp_rational_point(EFp *P){
    Fp tmp1,tmp2,tmp_x;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&tmp_x);
	gmp_randinit_default (state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
	
    while(1){
        Fp_set_random(&P->x,state);
        Fp_mul(&tmp1,&P->x,&P->x);
        Fp_mul(&tmp2,&tmp1,&P->x);
        Fp_add_mpz(&tmp_x,&tmp2,curve_b);
        if(Fp_legendre(&tmp_x)==1){
            Fp_sqrt(&P->y,&tmp_x);
            break;
        }
    }
    
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    Fp_clear(&tmp_x);
}

void EFp_ECD(EFp *ANS,EFp *P){
    if(Fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp Tmp_P;
    EFp_init(&Tmp_P);
    EFp_set(&Tmp_P,P);
    Fp tmp1,tmp2,lambda;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&lambda);
    
    Fp_add(&tmp1,&Tmp_P.y,&Tmp_P.y);
    Fp_inv(&tmp1,&tmp1);
    Fp_mul(&tmp2,&Tmp_P.x,&Tmp_P.x);
    Fp_add(&lambda,&tmp2,&tmp2);
    Fp_add(&tmp2,&tmp2,&lambda);
    Fp_mul(&lambda,&tmp1,&tmp2);
    Fp_mul(&tmp1,&lambda,&lambda);
    Fp_add(&tmp2,&Tmp_P.x,&Tmp_P.x);
    Fp_sub(&ANS->x,&tmp1,&tmp2);
    Fp_sub(&tmp1,&Tmp_P.x,&ANS->x);
    Fp_mul(&tmp2,&lambda,&tmp1);
    Fp_sub(&ANS->y,&tmp2,&Tmp_P.y);
    
    //clear
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    Fp_clear(&lambda);
    EFp_clear(&Tmp_P);
}

void EFp_ECA(EFp *ANS,EFp *P1,EFp *P2){
    if(P1->infinity==1){
        EFp_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp_set(ANS,P1);
        return;
    }else if(Fp_cmp(&P1->x,&P2->x)==0){
        if(Fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp_ECD(ANS,P1);
            return;
        }
    }
    
    EFp Tmp_P1,Tmp_P2;
    EFp_init(&Tmp_P1);
    EFp_set(&Tmp_P1,P1);
    EFp_init(&Tmp_P2);
    EFp_set(&Tmp_P2,P2);
    Fp tmp1,tmp2,lambda;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&lambda);
    
    Fp_sub(&tmp1,&Tmp_P2.x,&Tmp_P1.x);
    Fp_inv(&tmp1,&tmp1);
    Fp_sub(&tmp2,&Tmp_P2.y,&Tmp_P1.y);
    Fp_mul(&lambda,&tmp1,&tmp2);
    Fp_mul(&tmp1,&lambda,&lambda);
    Fp_sub(&tmp2,&tmp1,&Tmp_P1.x);
    Fp_sub(&ANS->x,&tmp2,&Tmp_P2.x);
    Fp_sub(&tmp1,&Tmp_P1.x,&ANS->x);
    Fp_mul(&tmp2,&lambda,&tmp1);
    Fp_sub(&ANS->y,&tmp2,&Tmp_P1.y);
        
    //clear 
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    Fp_clear(&lambda);
    EFp_clear(&Tmp_P1);
    EFp_clear(&Tmp_P2);
}

void EFp_SCM(EFp *ANS,EFp *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp_set(ANS,P);
        return;
    }
    
    EFp Tmp_P,Next_P;
    EFp_init(&Tmp_P);
    EFp_set(&Tmp_P,P);
    EFp_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp_set(&Next_P,&Tmp_P);
    for(i=1; binary[i]!='\0'; i++){
        EFp_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    
    EFp_set(ANS,&Next_P);
    
    EFp_clear(&Next_P);
    EFp_clear(&Tmp_P);
}

//skew frobenius map
void EFp_skew_frobenius_map_p2(EFp *ANS,EFp *A){
    Fp_mul_mpz(&ANS->x,&A->x,epsilon1);
    Fp_set_neg(&ANS->y,&A->y);
}

/*----------------------------------------------------------------------------*/
//EFp2
void EFp2_init(EFp2 *P){
    Fp2_init(&P->x);
    Fp2_init(&P->y);
    P->infinity=0;
}

void EFp2_set(EFp2 *ANS,EFp2 *A){
    Fp2_set(&ANS->x,&A->x);
    Fp2_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp2_set_ui(EFp2 *ANS,unsigned long int UI){
    Fp2_set_ui(&ANS->x,UI);
    Fp2_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp2_set_mpz(EFp2 *ANS,mpz_t A){
    Fp2_set_mpz(&ANS->x,A);
    Fp2_set_mpz(&ANS->y,A);
    ANS->infinity=0;
}

void EFp2_set_neg(EFp2 *ANS,EFp2 *A){
    Fp2_set(&ANS->x,&A->x);
    Fp2_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp2_clear(EFp2 *P){
    Fp2_clear(&P->x);
    Fp2_clear(&P->y);
}

void EFp2_printf(EFp2 *P,char *str){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp2_printf(&P->x,"");
        printf(",");
        Fp2_printf(&P->y,"");
        printf(")");
    }else{
        printf("0");
    }
}

void EFp2_rational_point(EFp2 *P){
    Fp2 tmp1,tmp2;
    Fp2_init(&tmp1);
    Fp2_init(&tmp2);
	gmp_randinit_default (state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp2_set_random(&P->x,state);
        Fp2_sqr(&tmp1,&P->x);
        Fp2_mul(&tmp2,&tmp1,&P->x);
        Fp_add_mpz(&tmp2.x0,&tmp2.x0,curve_b);
        if(Fp2_legendre(&tmp2)==1){
            Fp2_sqrt(&P->y,&tmp2);
            break;
        }
    }
    
    Fp2_clear(&tmp1);
    Fp2_clear(&tmp2);
}

void EFp2_ECD(EFp2 *ANS,EFp2 *P){
    if(Fp2_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp2 Tmp_P;
    EFp2_init(&Tmp_P);
    EFp2_set(&Tmp_P,P);
    Fp2 tmp1,tmp2,lambda;
    Fp2_init(&tmp1);
    Fp2_init(&tmp2);
    Fp2_init(&lambda);
    
    Fp2_add(&tmp1,&Tmp_P.y,&Tmp_P.y);
    
    Fp2_inv(&tmp1,&tmp1);
    Fp2_sqr(&tmp2,&Tmp_P.x);
    Fp2_add(&lambda,&tmp2,&tmp2);
    Fp2_add(&tmp2,&tmp2,&lambda);
    Fp2_mul(&lambda,&tmp1,&tmp2);
    
    Fp2_sqr(&tmp1,&lambda);
    Fp2_add(&tmp2,&Tmp_P.x,&Tmp_P.x);
    Fp2_sub(&ANS->x,&tmp1,&tmp2);
    
    Fp2_sub(&tmp1,&Tmp_P.x,&ANS->x);
    Fp2_mul(&tmp2,&lambda,&tmp1);
    Fp2_sub(&ANS->y,&tmp2,&Tmp_P.y);
    
    Fp2_clear(&tmp1);
    Fp2_clear(&tmp2);
    Fp2_clear(&lambda);
    EFp2_clear(&Tmp_P);
}

void EFp2_ECA(EFp2 *ANS,EFp2 *P1,EFp2 *P2){
    if(P1->infinity==1){
        EFp2_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp2_set(ANS,P1);
        return;
    }else if(Fp2_cmp(&P1->x,&P2->x)==0){
        if(Fp2_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp2_ECD(ANS,P1);
            return;
        }
    }
    
    EFp2 Tmp_P1,Tmp_P2;
    EFp2_init(&Tmp_P1);
    EFp2_set(&Tmp_P1,P1);
    EFp2_init(&Tmp_P2);
    EFp2_set(&Tmp_P2,P2);
    Fp2 tmp1,tmp2,lambda;
    Fp2_init(&tmp1);
    Fp2_init(&tmp2);
    Fp2_init(&lambda);
    
    Fp2_sub(&tmp1,&Tmp_P2.x,&Tmp_P1.x);
    Fp2_inv(&tmp1,&tmp1);
    Fp2_sub(&tmp2,&Tmp_P2.y,&Tmp_P1.y);
    Fp2_mul(&lambda,&tmp1,&tmp2);
    Fp2_sqr(&tmp1,&lambda);
    Fp2_sub(&tmp2,&tmp1,&Tmp_P1.x);
    Fp2_sub(&ANS->x,&tmp2,&Tmp_P2.x);
    Fp2_sub(&tmp1,&Tmp_P1.x,&ANS->x);
    Fp2_mul(&tmp2,&lambda,&tmp1);
    Fp2_sub(&ANS->y,&tmp2,&Tmp_P1.y);
        
    //clear 
    Fp2_clear(&tmp1);
    Fp2_clear(&tmp2);
    Fp2_clear(&lambda);
    EFp2_clear(&Tmp_P1);
    EFp2_clear(&Tmp_P2);
}

void EFp2_SCM(EFp2 *ANS,EFp2 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp2_set(ANS,P);
        return;
    }
    
    EFp2 Tmp_P,Next_P;
    EFp2_init(&Tmp_P);
    EFp2_set(&Tmp_P,P);
    EFp2_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp2_set(&Next_P,&Tmp_P);
    for(i=1; binary[i]!='\0'; i++){
        EFp2_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp2_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp2_set(ANS,&Next_P);
    
    EFp2_clear(&Next_P);
    EFp2_clear(&Tmp_P);
}

//skew_frobenius_map
void EFp2_skew_frobenius_map_p1(EFp2 *ANS,EFp2 *A){
    //x
    Fp_set(&ANS->x.x0,&A->x.x0);
    Fp_set_neg(&ANS->x.x1,&A->x.x1);
    Fp2_mul(&ANS->x,&ANS->x,&frobenius_constant[f_p1][1]);
    //y
    Fp_set(&ANS->y.x0,&A->y.x0);
    Fp_set_neg(&ANS->y.x1,&A->y.x1);
    Fp2_mul(&ANS->y,&ANS->y,&frobenius_constant[f_p1][4]);
}

void EFp2_skew_frobenius_map_p2(EFp2 *ANS,EFp2 *A){
    //x
	Fp2_mul(&ANS->x,&A->x,&frobenius_constant[f_p2][1]);
	//y
	Fp2_mul(&ANS->y,&A->y,&frobenius_constant[f_p2][4]);
}

void EFp2_skew_frobenius_map_p3(EFp2 *ANS,EFp2 *A){
    //x
    Fp_set(&ANS->x.x0,&A->x.x0);
    Fp_set_neg(&ANS->x.x1,&A->x.x1);
	Fp2_mul(&ANS->x,&ANS->x,&frobenius_constant[f_p3][1]);
	//y
	Fp_set(&ANS->y.x0,&A->y.x0);
    Fp_set_neg(&ANS->y.x1,&A->y.x1);
	Fp2_mul(&ANS->y,&ANS->y,&frobenius_constant[f_p3][4]);
}

void EFp2_skew_frobenius_map_p10(EFp2 *ANS,EFp2 *A){
    //x
	Fp2_mul(&ANS->x,&A->x,&frobenius_constant[f_p10][1]);
	//y
	Fp2_mul(&ANS->y,&A->y,&frobenius_constant[f_p10][4]);
}

/*----------------------------------------------------------------------------*/
//EFp6
void EFp6_init(EFp6 *P){
    Fp6_init(&P->x);
    Fp6_init(&P->y);
    P->infinity=0;
}

void EFp6_set(EFp6 *ANS,EFp6 *A){
    Fp6_set(&ANS->x,&A->x);
    Fp6_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp6_set_ui(EFp6 *ANS,unsigned long int UI){
    Fp6_set_ui(&ANS->x,UI);
    Fp6_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp6_set_mpz(EFp6 *ANS,mpz_t A){
    Fp6_set_mpz(&ANS->x,A);
    Fp6_set_mpz(&ANS->y,A);
    ANS->infinity=0;
}

void EFp6_set_neg(EFp6 *ANS,EFp6 *A){
    Fp6_set(&ANS->x,&A->x);
    Fp6_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp6_clear(EFp6 *P){
    Fp6_clear(&P->x);
    Fp6_clear(&P->y);
}

void EFp6_printf(EFp6 *P,char *str){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp6_printf(&P->x,"");
        printf(",");
        Fp6_printf(&P->y,"");
        printf(")");
    }else{
        printf("0");
    }
}

void EFp6_rational_point(EFp6 *P){
    Fp6 tmp1,tmp2;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
	gmp_randinit_default (state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp6_set_random(&P->x,state);
        Fp6_sqr(&tmp1,&P->x);
        Fp6_mul(&tmp2,&tmp1,&P->x);
        Fp_add_mpz(&tmp2.x0.x0,&tmp2.x0.x0,curve_b);
        if(Fp6_legendre(&tmp2)==1){
            Fp6_sqrt(&P->y,&tmp2);
            break;
        }
    }
    
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
}

void EFp6_ECD(EFp6 *ANS,EFp6 *P){
    if(Fp6_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp6 Tmp_P;
    EFp6_init(&Tmp_P);
    EFp6_set(&Tmp_P,P);
    Fp6 tmp1,tmp2,lambda;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    Fp6_init(&lambda);
    
    Fp6_add(&tmp1,&Tmp_P.y,&Tmp_P.y);
    
    Fp6_inv(&tmp1,&tmp1);
    Fp6_sqr(&tmp2,&Tmp_P.x);
    Fp6_add(&lambda,&tmp2,&tmp2);
    Fp6_add(&tmp2,&tmp2,&lambda);
    Fp6_mul(&lambda,&tmp1,&tmp2);
    
    Fp6_sqr(&tmp1,&lambda);
    Fp6_add(&tmp2,&Tmp_P.x,&Tmp_P.x);
    Fp6_sub(&ANS->x,&tmp1,&tmp2);
    
    Fp6_sub(&tmp1,&Tmp_P.x,&ANS->x);
    Fp6_mul(&tmp2,&lambda,&tmp1);
    Fp6_sub(&ANS->y,&tmp2,&Tmp_P.y);
    
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
    Fp6_clear(&lambda);
    EFp6_clear(&Tmp_P);
}

void EFp6_ECA(EFp6 *ANS,EFp6 *P1,EFp6 *P2){
    if(P1->infinity==1){
        EFp6_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp6_set(ANS,P1);
        return;
    }else if(Fp6_cmp(&P1->x,&P2->x)==0){
        if(Fp6_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp6_ECD(ANS,P1);
            return;
        }
    }
    
    EFp6 Tmp_P1,Tmp_P2;
    EFp6_init(&Tmp_P1);
    EFp6_set(&Tmp_P1,P1);
    EFp6_init(&Tmp_P2);
    EFp6_set(&Tmp_P2,P2);
    Fp6 tmp1,tmp2,lambda;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    Fp6_init(&lambda);
    
    Fp6_sub(&tmp1,&Tmp_P2.x,&Tmp_P1.x);
    Fp6_inv(&tmp1,&tmp1);
    Fp6_sub(&tmp2,&Tmp_P2.y,&Tmp_P1.y);
    Fp6_mul(&lambda,&tmp1,&tmp2);
    Fp6_sqr(&tmp1,&lambda);
    Fp6_sub(&tmp2,&tmp1,&Tmp_P1.x);
    Fp6_sub(&ANS->x,&tmp2,&Tmp_P2.x);
    Fp6_sub(&tmp1,&Tmp_P1.x,&ANS->x);
    Fp6_mul(&tmp2,&lambda,&tmp1);
    Fp6_sub(&ANS->y,&tmp2,&Tmp_P1.y);
        
    //clear 
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
    Fp6_clear(&lambda);
    EFp6_clear(&Tmp_P1);
    EFp6_clear(&Tmp_P2);
}

void EFp6_SCM(EFp6 *ANS,EFp6 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp6_set(ANS,P);
        return;
    }
    
    EFp6 Tmp_P,Next_P;
    EFp6_init(&Tmp_P);
    EFp6_set(&Tmp_P,P);
    EFp6_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp6_set(&Next_P,&Tmp_P);
    for(i=1; binary[i]!='\0'; i++){
        EFp6_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp6_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp6_set(ANS,&Next_P);
    
    EFp6_clear(&Next_P);
    EFp6_clear(&Tmp_P);
}
/*----------------------------------------------------------------------------*/
//EFp12
void EFp12_init(EFp12 *P){
    Fp12_init(&P->x);
    Fp12_init(&P->y);
    P->infinity=0;
}

void EFp12_set(EFp12 *ANS,EFp12 *A){
    Fp12_set(&ANS->x,&A->x);
    Fp12_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp12_set_ui(EFp12 *ANS,unsigned long int UI){
    Fp12_set_ui(&ANS->x,UI);
    Fp12_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp12_set_mpz(EFp12 *ANS,mpz_t A){
    Fp12_set_mpz(&ANS->x,A);
    Fp12_set_mpz(&ANS->y,A);
    ANS->infinity=0;
}

void EFp12_set_neg(EFp12 *ANS,EFp12 *A){
    Fp12_set(&ANS->x,&A->x);
    Fp12_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp12_clear(EFp12 *P){
    Fp12_clear(&P->x);
    Fp12_clear(&P->y);
}

void EFp12_printf(EFp12 *P,char *str){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp12_printf(&P->x,"");
        printf(",");
        Fp12_printf(&P->y,"");
        printf(")");
    }else{
        printf("0");
    }
}

void EFp12_rational_point(EFp12 *P){
    Fp12 tmp1,tmp2;
    Fp12_init(&tmp1);
    Fp12_init(&tmp2);
	gmp_randinit_default (state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp12_set_random(&P->x,state);
        Fp12_sqr(&tmp1,&P->x);
        Fp12_mul(&tmp2,&tmp1,&P->x);
        Fp_add_mpz(&tmp2.x0.x0.x0,&tmp2.x0.x0.x0,curve_b);
        if(Fp12_legendre(&tmp2)==1){
            Fp12_sqrt(&P->y,&tmp2);
            break;
        }
    }
    
    Fp12_clear(&tmp1);
    Fp12_clear(&tmp2);
}

void EFp12_generate_G1(EFp12 *P){
    EFp tmp_P;
    EFp_init(&tmp_P);
    
    EFp_rational_point(&tmp_P);
    EFp12_set_ui(P,0);
    Fp_set(&P->x.x0.x0.x0,&tmp_P.x);
    Fp_set(&P->y.x0.x0.x0,&tmp_P.y);
    P->infinity=tmp_P.infinity;
    
    EFp_clear(&tmp_P);
}

void EFp12_generate_G2(EFp12 *Q){
    EFp12 random_P,P,frobenius_P;
    EFp12_init(&random_P);
    EFp12_init(&P);
    EFp12_init(&frobenius_P);
    mpz_t exp;
    mpz_init(exp);
    
    EFp12_rational_point(&random_P);
    mpz_pow_ui(exp,order,2);
    mpz_tdiv_q(exp,EFp12_total,exp);
    EFp12_SCM(&P,&random_P,exp);
    Fp12_frobenius_map_p1(&frobenius_P.x,&P.x);
    Fp12_frobenius_map_p1(&frobenius_P.y,&P.y);
    EFp12_set_neg(&P,&P);
    EFp12_ECA(Q,&P,&frobenius_P);
    
    mpz_clear(exp);
    EFp12_clear(&random_P);
    EFp12_clear(&P);
    EFp12_clear(&frobenius_P);
}

void EFp12_ECD(EFp12 *ANS,EFp12 *P){
    if(Fp12_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp12 Tmp_P;
    EFp12_init(&Tmp_P);
    EFp12_set(&Tmp_P,P);
    Fp12 tmp1,tmp2,lambda;
    Fp12_init(&tmp1);
    Fp12_init(&tmp2);
    Fp12_init(&lambda);
    
    Fp12_add(&tmp1,&Tmp_P.y,&Tmp_P.y);
    
    Fp12_inv(&tmp1,&tmp1);
    Fp12_sqr(&tmp2,&Tmp_P.x);
    Fp12_add(&lambda,&tmp2,&tmp2);
    Fp12_add(&tmp2,&tmp2,&lambda);
    Fp12_mul(&lambda,&tmp1,&tmp2);
    
    Fp12_sqr(&tmp1,&lambda);
    Fp12_add(&tmp2,&Tmp_P.x,&Tmp_P.x);
    Fp12_sub(&ANS->x,&tmp1,&tmp2);
    
    Fp12_sub(&tmp1,&Tmp_P.x,&ANS->x);
    Fp12_mul(&tmp2,&lambda,&tmp1);
    Fp12_sub(&ANS->y,&tmp2,&Tmp_P.y);
    
    Fp12_clear(&tmp1);
    Fp12_clear(&tmp2);
    Fp12_clear(&lambda);
    EFp12_clear(&Tmp_P);
}

void EFp12_ECA(EFp12 *ANS,EFp12 *P1,EFp12 *P2){
    if(P1->infinity==1){
        EFp12_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp12_set(ANS,P1);
        return;
    }else if(Fp12_cmp(&P1->x,&P2->x)==0){
        if(Fp12_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp12_ECD(ANS,P1);
            return;
        }
    }
    
    EFp12 Tmp_P1,Tmp_P2;
    EFp12_init(&Tmp_P1);
    EFp12_set(&Tmp_P1,P1);
    EFp12_init(&Tmp_P2);
    EFp12_set(&Tmp_P2,P2);
    Fp12 tmp1,tmp2,lambda;
    Fp12_init(&tmp1);
    Fp12_init(&tmp2);
    Fp12_init(&lambda);
    
    Fp12_sub(&tmp1,&Tmp_P2.x,&Tmp_P1.x);
    Fp12_inv(&tmp1,&tmp1);
    Fp12_sub(&tmp2,&Tmp_P2.y,&Tmp_P1.y);
    Fp12_mul(&lambda,&tmp1,&tmp2);
    Fp12_sqr(&tmp1,&lambda);
    Fp12_sub(&tmp2,&tmp1,&Tmp_P1.x);
    Fp12_sub(&ANS->x,&tmp2,&Tmp_P2.x);
    Fp12_sub(&tmp1,&Tmp_P1.x,&ANS->x);
    Fp12_mul(&tmp2,&lambda,&tmp1);
    Fp12_sub(&ANS->y,&tmp2,&Tmp_P1.y);
        
    //clear 
    Fp12_clear(&tmp1);
    Fp12_clear(&tmp2);
    Fp12_clear(&lambda);
    EFp12_clear(&Tmp_P1);
    EFp12_clear(&Tmp_P2);
}

void EFp12_SCM(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp12_set(ANS,P);
        return;
    }
    
    EFp12 Tmp_P,Next_P;
    EFp12_init(&Tmp_P);
    EFp12_set(&Tmp_P,P);
    EFp12_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp12_set(&Next_P,&Tmp_P);
    for(i=1; binary[i]!='\0'; i++){
        EFp12_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp12_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp12_set(ANS,&Next_P);
    
    EFp12_clear(&Next_P);
    EFp12_clear(&Tmp_P);
}

/*============================================================================*/
/* Pairing functions                                                          */
/*============================================================================*/
/*----------------------------------------------------------------------------*/
//protptype
void EFp12_lTP(Fp12 *ANS,EFp12 *T,EFp12 *P,EFp12 *Q,Fp *L){
	EFp12 Tmp1,Tmp2,Tmp3;
	EFp12_init(&Tmp1);
	EFp12_set(&Tmp1,T);
	EFp12_init(&Tmp2);
	EFp12_set(&Tmp2,P);
	EFp12_init(&Tmp3);
	EFp12_set(&Tmp3,Q);
	Fp12 Buf1,Buf2,Buf3,ZERO;
	Fp12_init(&Buf1);
	Fp12_init(&Buf2);
	Fp12_init(&Buf3);
	Fp12_init(&ZERO);
	Fp12_set_ui(&ZERO,0);
	
	Fp12_sub(&Buf1,&Tmp2.y,&Tmp1.y);
	Fp12_sub(&Buf3,&Tmp2.x,&Tmp1.x);
	if(Fp12_cmp(&Buf3,&ZERO)==0){
		Fp12_sub(ANS,&Tmp3.x,&Tmp1.x);
	}else{
		Fp12_inv(&Buf2,&Buf3);
		Fp12_mul(&Buf3,&Buf1,&Buf2);
		Fp12_sub(&Buf1,&Tmp3.x,&Tmp2.x);
		Fp12_mul(&Buf2,&Buf1,&Buf3);
		Fp12_sub(&Buf1,&Tmp3.y,&Tmp2.y);
		Fp12_sub(ANS,&Buf1,&Buf2);
		
		Fp12_mul_mpz(ANS,ANS,L->x0);
		//Fp12_set_neg(ANS,ANS);
	}
	
	Fp12_clear(&Buf1);
	Fp12_clear(&Buf2);
	Fp12_clear(&Buf3);
	Fp12_clear(&ZERO);
	EFp12_clear(&Tmp1);
	EFp12_clear(&Tmp2);
	EFp12_clear(&Tmp3);
}

void EFp12_lTT(Fp12 *ANS,EFp12 *T,EFp12 *Q,Fp *L){
	EFp12 Tmp1,Tmp2;
	EFp12_init(&Tmp1);
	EFp12_set(&Tmp1,T);
	EFp12_init(&Tmp2);
	EFp12_set(&Tmp2,Q);
	Fp12 Buf1,Buf2,Buf3,ZERO;
	Fp12_init(&Buf1);
	Fp12_init(&Buf2);
	Fp12_init(&Buf3);
	Fp12_init(&ZERO);
	Fp12_set_ui(&ZERO,0);
	mpz_t k;
	mpz_init(k);
	
	if(Fp12_cmp(&Tmp1.y,&ZERO)==0){
		Fp12_sub(ANS,&Tmp2.x,&Tmp1.x);
	}else{
		Fp12_sqr(&Buf1,&Tmp1.x);
		Fp12_mul_ui(&Buf1,&Buf1,3);
		Fp12_mul_ui(&Buf2,&Tmp1.y,2);
		Fp12_inv(&Buf3,&Buf2);
		Fp12_mul(&Buf2,&Buf1,&Buf3);
		
		Fp12_sub(&Buf1,&Tmp2.x,&Tmp1.x);
		Fp12_mul(&Buf3,&Buf1,&Buf2);
		Fp12_sub(&Buf1,&Tmp2.y,&Tmp1.y);
		Fp12_sub(ANS,&Buf1,&Buf3);
		
		Fp12_mul_mpz(ANS,ANS,L->x0);
		//Fp12_set_neg(ANS,ANS);
	}
	
	mpz_clear(k);
	Fp12_clear(&Buf1);
	Fp12_clear(&Buf2);
	Fp12_clear(&Buf3);
	Fp12_clear(&ZERO);
	EFp12_clear(&Tmp1);
	EFp12_clear(&Tmp2);
}

void EFp12_vTP(Fp12 *ANS,EFp12 *T,EFp12 *P,EFp12 *Q){
	EFp12 Tmp1,Tmp2,Tmp3;
	EFp12_init(&Tmp1);
	EFp12_set(&Tmp1,T);
	EFp12_init(&Tmp2);
	EFp12_set(&Tmp2,P);
	EFp12_init(&Tmp3);
	EFp12_set(&Tmp3,Q);
	
	EFp12_ECA(&Tmp1,&Tmp1,&Tmp2);
	Fp12_sub(ANS,&Tmp3.x,&Tmp1.x);
	
	EFp12_clear(&Tmp1);
	EFp12_clear(&Tmp2);
	EFp12_clear(&Tmp3);
}

void EFp12_vTT(Fp12 *ANS,EFp12 *T,EFp12 *Q){
	EFp12 Tmp1,Tmp2;
	EFp12_init(&Tmp1);
	EFp12_set(&Tmp1,T);
	EFp12_init(&Tmp2);
	EFp12_set(&Tmp2,Q);
	
	EFp12_ECD(&Tmp1,&Tmp1);
	Fp12_sub(ANS,&Tmp2.x,&Tmp1.x);
	
	EFp12_clear(&Tmp1);
	EFp12_clear(&Tmp2);
}

void Miller_prototype(Fp12 *ANS,EFp12 *P,EFp12 *Q){
	EFp12 Tmp1,Tmp2,T,N;
	EFp12_init(&Tmp1);
	EFp12_init(&Tmp2);
	EFp12_init(&T);
	EFp12_init(&N);
	EFp twist_P;
	EFp_init(&twist_P);
	EFp2 twist_Q;
	EFp2_init(&twist_Q);
	Fp12 f,Buf1,Buf2,Buf3,Buf4;
	Fp12_init(&f);
	Fp12_init(&Buf1);
	Fp12_init(&Buf2);
	Fp12_init(&Buf3);
	Fp12_init(&Buf4);
	Fp L;
	Fp_init(&L);
	mpz_t buf;
	mpz_init(buf);
	mpz_sub_ui(buf,trace,1);
	int i,length;
	length=(int)mpz_sizeinbase(buf,2);
	char binary[length];
	mpz_get_str(binary,2,buf);
	
	EFp12_set(&Tmp1,P);
	EFp12_set(&Tmp2,Q);
	EFp12_to_EFp(&twist_P,&Tmp2);
	EFp12_to_EFp2(&twist_Q,&Tmp1);
	Pseudo_8_sparse_mapping(&twist_P,&twist_Q,&L);
	EFp_to_EFp12(&Tmp2,&twist_P);
	EFp2_to_EFp12(&Tmp1,&twist_Q);
	
	EFp12_set(&T,&Tmp1);
	Fp_set_ui(&f.x0.x0.x0,1);
	for(i=1; binary[i]!='\0'; i++){
		EFp12_ECD(&N,&T);
		if(N.infinity==1){
			Fp12_sqr(&Buf1,&f);
			EFp12_lTT(&Buf2,&T,&Tmp2,&L);
			Fp12_mul(&f,&Buf1,&Buf2);
		}else{
			Fp12_sqr(&Buf1,&f);
			EFp12_lTT(&Buf2,&T,&Tmp2,&L);
			//EFp12_vTT(&Buf3,&T,&Tmp2);
			//Fp12_inv(&Buf4,&Buf3);
			//Fp12_mul(&Buf3,&Buf1,&Buf2);
			//Fp12_mul(&f,&Buf3,&Buf4);
			Fp12_mul(&f,&Buf1,&Buf2);
		}
		EFp12_set(&T,&N);
		if(binary[i]=='1'){
			EFp12_ECA(&N,&T,&Tmp1);
			if(N.infinity==1){
				EFp12_lTP(&Buf1,&T,&Tmp1,&Tmp2,&L);
				Fp12_mul(&f,&f,&Buf1);
			}else{
				EFp12_lTP(&Buf1,&T,&Tmp1,&Tmp2,&L);
				//EFp12_vTP(&Buf2,&T,&Tmp1,&Tmp2);
				//Fp12_inv(&Buf3,&Buf2);
				//Fp12_mul(&Buf2,&f,&Buf1);
				//Fp12_mul(&f,&Buf2,&Buf3);
				Fp12_mul(&f,&f,&Buf1);
			}
			EFp12_set(&T,&N);	
		}
	}
	
	Fp12_set(ANS,&f);
	
	EFp_clear(&twist_P);
	EFp2_clear(&twist_Q);
	mpz_clear(buf);
	Fp12_clear(&f);
	Fp12_clear(&Buf1);
	Fp12_clear(&Buf2);
	Fp12_clear(&Buf3);
	Fp12_clear(&Buf4);
	EFp12_clear(&T);
	EFp12_clear(&N);
	EFp12_clear(&Tmp1);
	EFp12_clear(&Tmp2);
	Fp_clear(&L);
}

void Prototype_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    Miller_prototype(ANS,P,Q);
    Final_exp_optimal(ANS,ANS);
}

/*----------------------------------------------------------------------------*/
//twist
void EFp12_to_EFp2(EFp2 *ANS,EFp12 *A){
    Fp2_set(&ANS->x,&A->x.x0.x1);
    Fp2_set(&ANS->y,&A->y.x1.x1);
    ANS->infinity=A->infinity;
}

void EFp2_to_EFp12(EFp12 *ANS,EFp2 *A){
    Fp12_set_ui(&ANS->x,0);
    Fp12_set_ui(&ANS->y,0);
    Fp2_set(&ANS->x.x0.x1,&A->x);
    Fp2_set(&ANS->y.x1.x1,&A->y);
    ANS->infinity=A->infinity;
}

void EFp12_to_EFp(EFp *ANS,EFp12 *A){
    Fp_set(&ANS->x,&A->x.x0.x0.x0);
    Fp_set(&ANS->y,&A->y.x0.x0.x0);
    ANS->infinity=A->infinity;
}

void EFp_to_EFp12(EFp12 *ANS,EFp *A){
    Fp12_set_ui(&ANS->x,0);
    Fp12_set_ui(&ANS->y,0);
    Fp_set(&ANS->x.x0.x0.x0,&A->x);
    Fp_set(&ANS->y.x0.x0.x0,&A->y);
    ANS->infinity=A->infinity;
}

/*----------------------------------------------------------------------------*/
//Pseudo 8-sparse
void Pseudo_8_sparse_mapping(EFp *P,EFp2 *Q,Fp *L){
    EFp2 Tmp_Q;
	EFp2_init(&Tmp_Q);
	EFp Tmp_P;
	EFp_init(&Tmp_P);
	Fp A,B,C,D,c;
	Fp_init(&A);
	Fp_init(&B);
	Fp_init(&C);
	Fp_init(&D);
	Fp_init(&c);
	
	EFp_set(&Tmp_P,P);
	EFp2_set(&Tmp_Q,Q);
	
	Fp_mul(&A,&Tmp_P.x,&Tmp_P.y);
	Fp_inv(&A,&A);
	Fp_mul(&B,&Tmp_P.x,&Tmp_P.x);
	Fp_mul(&B,&B,&A);
	Fp_mul(&C,&Tmp_P.y,&A);
	Fp_mul(&D,&B,&B);
	
	Fp2_mul_mpz(&Q->x,&Tmp_Q.x,D.x0);
	Fp_mul(&c,&B,&D);
	Fp2_mul_mpz(&Q->y,&Tmp_Q.y,c.x0);
	
	Fp_mul(&P->x,&D,&Tmp_P.x);
	Fp_set(&P->y,&P->x);
	
	Fp_mul(L,&C,&Tmp_P.y);
	Fp_mul(L,L,L);
	Fp_mul(L,L,&C);
	
	
	EFp2_clear(&Tmp_Q);
	EFp_clear(&Tmp_P);
	Fp_clear(&A);
	Fp_clear(&B);
	Fp_clear(&C);
	Fp_clear(&D);
	Fp_clear(&c);
}

void Pseudo_8_sparse_mul(Fp12 *ANS,Fp12 *A,Fp12 *B){
    //A= f0 + f1γ^2 + f2γ^4 + f3γ　+ f4γ^3 + f5γ^5
	//B= 1                  +  aγ +  bγ^3
	// x0.x0  x0.x1  x0.x2  x1.x0   x1.x1   x1.x2
	Fp12 ans;
	Fp12_init(&ans);
	Fp2 tmp0,tmp1,tmp2,tmp3;
	Fp2_init(&tmp0);
	Fp2_init(&tmp1);
	Fp2_init(&tmp2);
	Fp2_init(&tmp3);
	
	Fp2_mul(&tmp0,&A->x0.x0,&B->x1.x0);		//tmp0=b3*f0
	Fp2_mul(&tmp1,&A->x0.x1,&B->x1.x1);		//tmp1=b4*f1
	Fp2_add(&tmp2,&A->x0.x0,&A->x0.x1);		//tmp2=f0+f1
	Fp2_add(&tmp3,&B->x1.x0,&B->x1.x1);		//tmp3=b3+b4
	Fp2_mul(&tmp2,&tmp2,&tmp3);			//tmp2=tmp2*tmp3
	Fp2_sub(&tmp2,&tmp2,&tmp0);			//tmp2=tmp2-tmp0
	Fp2_sub(&tmp2,&tmp2,&tmp1);			//tmp2=tmp2-tmp1
	Fp2_add(&ans.x1.x1,&tmp2,&A->x1.x1);	//ans[γ^3]=tmp2+f4
	Fp2_add(&tmp0,&tmp0,&A->x1.x0);		//tmp0=tmp0+f3
	Fp2_mul(&tmp2,&A->x0.x2,&B->x1.x1);		//tmp2=b4*f2
	Fp2_mul_basis(&tmp2,&tmp2);			//tmp2=tmp2*(α+1)
	Fp2_add(&ans.x1.x0,&tmp0,&tmp2);		//ans[γ]=tmp0+tmp2
	Fp2_add(&tmp0,&tmp1,&A->x1.x2);		//tmp0=tmp1+f5
	Fp2_mul(&tmp1,&A->x0.x2,&B->x1.x0);		//tmp1=b3*f2
	Fp2_add(&ans.x1.x2,&tmp0,&tmp1);		//ans[γ^5]=tmp0+tmp1
	Fp2_mul(&tmp0,&A->x1.x0,&B->x1.x0);		//tmp0=b3*f3
	Fp2_mul(&tmp1,&A->x1.x1,&B->x1.x1);		//tmp1=b4*f4
	Fp2_add(&tmp2,&A->x1.x0,&A->x1.x1);		//tmp2=f3+f4
	Fp2_mul(&tmp2,&tmp2,&tmp3);			//tmp2=tmp2*tmp3
	Fp2_sub(&tmp2,&tmp2,&tmp0);			//tmp2=tmp2-tmp0
	Fp2_sub(&tmp2,&tmp2,&tmp1);			//tmp2=tmp2-tmp1
	Fp2_add(&ans.x0.x2,&tmp2,&A->x0.x2);	//ans[γ^4]=tmp2+f4
	Fp2_add(&tmp0,&tmp0,&A->x0.x1);		//tmp0=tmp0+f1
	Fp2_mul(&tmp2,&A->x1.x2,&B->x1.x1);		//tmp2=b4*f5
	Fp2_mul_basis(&tmp2,&tmp2);			//tmp2=tmp2*(α+1)
	Fp2_add(&ans.x0.x1,&tmp0,&tmp2);		//ans[γ^2]=tmp0+tmp2
	Fp2_mul(&tmp0,&A->x1.x2,&B->x1.x0);		//tmp0=b3*f5
	Fp2_add(&tmp0,&tmp0,&tmp1);			//tmp0=tmp0+tmp1
	Fp2_mul_basis(&tmp0,&tmp0);			//tmp0=tmp0*(α+1)
	Fp2_add(&ans.x0.x0,&tmp0,&A->x0.x0);	//ans[1]=tmp0+f0
	Fp12_set(ANS,&ans);
	
	Fp12_clear(&ans);
	Fp2_clear(&tmp0);
	Fp2_clear(&tmp1);
	Fp2_clear(&tmp2);
	Fp2_clear(&tmp3);
}

void ff_ltt(Fp12 *f,EFp2 *T,EFp *P,Fp *L){
    EFp2 Tmp_T;
	EFp2_init(&Tmp_T);
	Fp12 ff,ltt;
	Fp12_init(&ff);
	Fp12_init(&ltt);
	Fp2 A,B,C,D,E;
	Fp2_init(&A);
	Fp2_init(&B);
	Fp2_init(&C);
	Fp2_init(&D);
	Fp2_init(&E);
	EFp2_set(&Tmp_T,T);
	
	Fp12_sqr(&ff,f);
	
	//ltt
	Fp2_add(&A,&Tmp_T.y,&Tmp_T.y);		//A=1/(2*T.y)
	Fp2_inv(&A,&A);
	Fp2_sqr(&B,&Tmp_T.x);			//B=3(T.x)^2
	Fp2_add(&C,&B,&B);
	Fp2_add(&B,&C,&B);
	Fp2_mul(&C,&A,&B);				//C=A*B
	
	Fp2_add(&D,&Tmp_T.x,&Tmp_T.x);		//D=2T.x
	Fp2_sqr(&T->x,&C);				//next_T.x=C^2-D
	Fp2_sub(&T->x,&T->x,&D);
	Fp2_mul(&E,&C,&Tmp_T.x);			//E=C*T.x-T.y
	Fp2_sub(&E,&E,&Tmp_T.y);
	Fp2_mul(&T->y,&C,&T->x);			//next_T.y=E-C*next_T.x
	Fp2_sub(&T->y,&E,&T->y);
	
	//set ltt
	Fp_set_ui(&ltt.x0.x0.x0,1);
	Fp2_set_neg(&ltt.x1.x0,&C);
	Fp2_mul_mpz(&ltt.x1.x1,&E,L->x0);
	
	Pseudo_8_sparse_mul(f,&ff,&ltt);
	
	EFp2_clear(&Tmp_T);
	Fp2_clear(&A);
	Fp2_clear(&B);
	Fp2_clear(&C);
	Fp2_clear(&D);
	Fp2_clear(&E);
	Fp12_clear(&ff);
	Fp12_clear(&ltt);
}

void f_ltq(Fp12 *f,EFp2 *T,EFp2 *Q,EFp *P,Fp *L){
    EFp2 Tmp_T;
	EFp2_init(&Tmp_T);
	Fp12 ltq;
	Fp12_init(&ltq);
	Fp2 A,B,C,D,E;
	Fp2_init(&A);
	Fp2_init(&B);
	Fp2_init(&C);
	Fp2_init(&D);
	Fp2_init(&E);
	EFp2_set(&Tmp_T,T);
		
	//ltq
	Fp2_sub(&A,&Q->x,&Tmp_T.x);		//A=(Q->x-T.x)^-1
	Fp2_inv(&A,&A);
	Fp2_sub(&B,&Q->y,&Tmp_T.y);		//B=(Q->y-T.y)
	Fp2_mul(&C,&A,&B);			//C=A*B
	Fp2_add(&D,&Tmp_T.x,&Q->x);		//D=Q->x+T.x
	Fp2_sqr(&T->x,&C);			//next_T.x=C^2-D
	Fp2_sub(&T->x,&T->x,&D);
	Fp2_mul(&E,&C,&Tmp_T.x);		//E=C*T.x-T.y
	Fp2_sub(&E,&E,&Tmp_T.y);
	Fp2_mul(&T->y,&C,&T->x);		//next_T.y=E-C*next_T.x
	Fp2_sub(&T->y,&E,&T->y);
	
	//set ltq
	Fp_set_ui(&ltq.x0.x0.x0,1);
	Fp2_set_neg(&ltq.x1.x0,&C);
	Fp2_mul_mpz(&ltq.x1.x1,&E,L->x0);
	
	Pseudo_8_sparse_mul(f,f,&ltq);
	
	EFp2_clear(&Tmp_T);
	Fp12_clear(&ltq);
	Fp2_clear(&A);
	Fp2_clear(&B);
	Fp2_clear(&C);
	Fp2_clear(&D);
}

/*----------------------------------------------------------------------------*/
//miller
void Miller_algo_for_plain_ate(Fp12 *ANS,EFp12 *Q,EFp12 *P){
    EFp12 Test;
	EFp12_init(&Test);
	EFp2 T;
	EFp2_init(&T);
	EFp2 mapped_Q;
	EFp2_init(&mapped_Q);
	EFp mapped_P;
	EFp_init(&mapped_P);
	Fp12 f;
	Fp12_init(&f);
	Fp L;
	Fp_init(&L);
	mpz_t loop;
	mpz_init(loop);
	mpz_sub_ui(loop,trace,1);
	int i,length;
	length=(int)mpz_sizeinbase(loop,2);
	char binary[length];
	mpz_get_str(binary,2,loop);
	
	EFp12_to_EFp(&mapped_P,P);
	EFp12_to_EFp2(&mapped_Q,Q);
	Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
	EFp2_set(&T,&mapped_Q);     //set T
	Fp_set_ui(&f.x0.x0.x0,1);   //set f
	
	//miller
    for(i=1; binary[i]!='\0'; i++){
		ff_ltt(&f,&T,&mapped_P,&L);
		if(binary[i]=='1'){
			f_ltq(&f,&T,&mapped_Q,&mapped_P,&L);
		}
	}
	
	Fp12_set(ANS,&f);
	
	Fp12_clear(&f);
	EFp2_clear(&T);
	EFp2_clear(&mapped_Q);
	EFp_clear(&mapped_P);
	Fp_clear(&L);
	EFp12_clear(&Test);
	mpz_clear(loop);
}

void Miller_algo_for_opt_ate(Fp12 *ANS,EFp12 *Q,EFp12 *P){
    EFp12 Buf;
    EFp12_init(&Buf);
    EFp2 T;
    EFp2_init(&T);
    EFp2 mapped_Q,mapped_Q_neg,mapped_Q1,mapped_Q2_neg;
    EFp2_init(&mapped_Q);
    EFp2_init(&mapped_Q_neg);
    EFp2_init(&mapped_Q1);
    EFp2_init(&mapped_Q2_neg);
    EFp mapped_P;
    EFp_init(&mapped_P);
    Fp12 f;
    Fp12_init(&f);
    Fp L;
    Fp_init(&L);
    int i;
    
    //set
    EFp12_to_EFp(&mapped_P,P);//set mapped_P
    EFp12_to_EFp2(&mapped_Q,Q);//set mapped_Q
    Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
    EFp2_set(&mapped_Q_neg,&mapped_Q);//set mapped_Q_neg
    Fp2_set_neg(&mapped_Q_neg.y,&mapped_Q_neg.y);
    EFp2_set(&T,&mapped_Q);     //set T
    Fp_set_ui(&f.x0.x0.x0,1);
    //miller
    for(i=X6_2_length-1; i>=0; i--){
        switch(X6_2_binary[i]){
            case 0:
                ff_ltt(&f,&T,&mapped_P,&L);
                break;
            case 1:
                ff_ltt(&f,&T,&mapped_P,&L);
                f_ltq(&f,&T,&mapped_Q,&mapped_P,&L);
                break;
            case -1:
                ff_ltt(&f,&T,&mapped_P,&L);
                f_ltq(&f,&T,&mapped_Q_neg,&mapped_P,&L);
                break;
            default:
                break;
        }
        
    }
    
    EFp2_skew_frobenius_map_p1(&mapped_Q1,&mapped_Q);//Q^p
    EFp2_skew_frobenius_map_p2(&mapped_Q2_neg,&mapped_Q);//Q^(p^2)
    EFp2_set_neg(&mapped_Q2_neg,&mapped_Q2_neg);
    f_ltq(&f,&T,&mapped_Q1,&mapped_P,&L);
    f_ltq(&f,&T,&mapped_Q2_neg,&mapped_P,&L);
    
    Fp12_set(ANS,&f);
    
    EFp12_clear(&Buf);
    Fp12_clear(&f);
    EFp2_clear(&T);
    EFp2_clear(&mapped_Q);
    EFp2_clear(&mapped_Q_neg);
    EFp2_clear(&mapped_Q1);
    EFp2_clear(&mapped_Q2_neg);
    EFp_clear(&mapped_P);
    Fp_clear(&L);
}

void Miller_algo_for_x_ate(Fp12 *ANS,EFp12 *Q,EFp12 *P){
    EFp12 Buf;
    EFp12_init(&Buf);
    EFp2 T;
    EFp2_init(&T);
    EFp2 mapped_Q,mapped_Q_neg,mapped_Q1,mapped_Q2;
    EFp2_init(&mapped_Q);
    EFp2_init(&mapped_Q_neg);
    EFp2_init(&mapped_Q1);
    EFp2_init(&mapped_Q2);
    EFp mapped_P;
    EFp_init(&mapped_P);
    Fp12 f;
    Fp12_init(&f);
    Fp L;
    Fp_init(&L);
    int i;
    
    //set
    EFp12_to_EFp(&mapped_P,P);//set mapped_P
    EFp12_to_EFp2(&mapped_Q,Q);//set mapped_Q
    Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
    EFp2_set(&mapped_Q_neg,&mapped_Q);//set mapped_Q_neg
    Fp2_set_neg(&mapped_Q_neg.y,&mapped_Q_neg.y);
    EFp2_set(&T,&mapped_Q);     //set T
    Fp_set_ui(&f.x0.x0.x0,1);
    //miller
    for(i=X_length-1; i>=0; i--){
        switch(X_binary[i]){
            case 0:
                ff_ltt(&f,&T,&mapped_P,&L);
                break;
            case 1:
                ff_ltt(&f,&T,&mapped_P,&L);
                f_ltq(&f,&T,&mapped_Q,&mapped_P,&L);
                break;
            case -1:
                ff_ltt(&f,&T,&mapped_P,&L);
                f_ltq(&f,&T,&mapped_Q_neg,&mapped_P,&L);
                break;
            default:
                break;
        }
    }
    
    Fp12_frobenius_map_p3(&Buf.x,&f);       //f←f*f^(p^3)
    Fp12_mul(&f,&Buf.x,&f);
    EFp2_skew_frobenius_map_p3(&mapped_Q1,&T);//mapped_Q1←T^(p^3)
    EFp2_set(&mapped_Q2,&T);
    f_ltq(&f,&mapped_Q2,&mapped_Q1,&mapped_P,&L);
    Fp12_frobenius_map_p10(&Buf.x,&f);      //f←f*f^(p^10)
    Fp12_mul(&f,&Buf.x,&f);
    EFp2_skew_frobenius_map_p10(&T,&mapped_Q2);//T←Q2^(p^10)
    f_ltq(&f,&T,&mapped_Q2,&mapped_P,&L);
    
    Fp12_set(ANS,&f);
    
    EFp12_clear(&Buf);
    Fp12_clear(&f);
    EFp2_clear(&T);
    EFp2_clear(&mapped_Q);
    EFp2_clear(&mapped_Q_neg);
    EFp2_clear(&mapped_Q1);
    EFp2_clear(&mapped_Q2);
    EFp_clear(&mapped_P);
    Fp_clear(&L);
}

/*----------------------------------------------------------------------------*/
//final exp
void Final_exp_plain(Fp12 *ANS,Fp12 *A){
    Fp12 Tmp,Buf1,Buf2;
	Fp12_init(&Tmp);
	Fp12_set(&Tmp,A);
	Fp12_init(&Buf1);
	Fp12_init(&Buf2);
	mpz_t exp,buf;
	mpz_init(exp);
	mpz_init(buf);
	
	Fp12_frobenius_map_p6(&Buf1,&Tmp);
	Fp12_inv(&Buf2,&Tmp);
	Fp12_mul(&Tmp,&Buf1,&Buf2);
	
	Fp12_frobenius_map_p2(&Buf1,&Tmp);
	Fp12_mul(&Tmp,&Buf1,&Tmp);
	
	mpz_pow_ui(exp,prime,4);
	mpz_pow_ui(buf,prime,2);
	mpz_sub(exp,exp,buf);
	mpz_add_ui(exp,exp,1);
	mpz_tdiv_q(exp,exp,order);
	Fp12_pow(ANS,&Tmp,exp);
	
	mpz_clear(exp);
	mpz_clear(buf);
	Fp12_clear(&Tmp);
	Fp12_clear(&Buf1);
	Fp12_clear(&Buf2);
}

void Fp12_pow_X(Fp12 *ANS,Fp12 *A){
    int i;
    Fp12 tmp,A_inv;
	Fp12_init(&tmp);
	Fp12_init(&A_inv);
	Fp12_frobenius_map_p6(&A_inv,A);
	
	Fp12_set(&tmp,A);
	for(i=X_length-1; i>=0; i--){
		switch(X_binary[i]){
			case 0:
				Fp12_sqr_cyclotomic(&tmp,&tmp);
				break;
			case 1:
				Fp12_sqr_cyclotomic(&tmp,&tmp);
				Fp12_mul(&tmp,&tmp,A);
				break;
			case -1:
				Fp12_sqr_cyclotomic(&tmp,&tmp);
				Fp12_mul(&tmp,&tmp,&A_inv);
				break;
			default:
				break;
		}
	}
	Fp12_set(ANS,&tmp);
	
	Fp12_clear(&tmp);
	Fp12_clear(&A_inv);
}

void Final_exp_optimal(Fp12 *ANS,Fp12 *A){
    Fp12 t0,t1,t2,t3,t4;
    Fp12_init(&t0);
    Fp12_init(&t1);
    Fp12_init(&t2);
    Fp12_init(&t3);
    Fp12_init(&t4);
    
    gettimeofday(&tv_start,NULL);
    //f←f^(p^6)*f^-1
    Fp12_frobenius_map_p6(&t0,A);//f^(p^6)
    Fp12_inv(&t1,A);//f^-1
    Fp12_mul(A,&t0,&t1);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    Fp12_frobenius_map_p2(&t0,A);//f^(p^2)
    Fp12_mul(A,&t0,A);//f^(p^2)*f
    
    gettimeofday(&tv_end,NULL);
    FINALEXP_OPT_EASY=timedifference_msec(tv_start,tv_end);
    
    gettimeofday(&tv_start,NULL);
    
    Fp12_pow_X(&t0,A);   //t0←f^(-u)
    Fp12_frobenius_map_p6(&t0,&t0);
    Fp12_sqr_cyclotomic(&t0,&t0);              //t0←t0^2
    Fp12_sqr_cyclotomic(&t1,&t0);              //t1←t0^2
    Fp12_mul(&t1,&t0,&t1);              //t1←t0*t1
    Fp12_pow_X(&t2,&t1);        //t2←t1^(-u)
    Fp12_frobenius_map_p6(&t2,&t2);
    Fp12_frobenius_map_p6(&t3,&t1);         //t3←t1^-1
    Fp12_mul(&t1,&t2,&t3);              //t1←t2*t3
    Fp12_sqr_cyclotomic(&t3,&t2);              //t3←t2^2
    Fp12_pow_X(&t4,&t3);        //t4←t3^(-u)
    Fp12_frobenius_map_p6(&t4,&t4);
    Fp12_frobenius_map_p6(&t4,&t4);         //t4←t4^(-1)
    Fp12_mul(&t4,&t4,&t1);              //t4←t4*t1
    Fp12_mul(&t3,&t4,&t0);              //t3←t4*t0
    Fp12_mul(&t0,&t2,&t4);              //t0←t2*t4
    Fp12_mul(&t0,&t0,A);             //t0←t0*f
    Fp12_frobenius_map_p1(&t2,&t3);         //t2←t3^p
    Fp12_mul(&t0,&t2,&t0);              //t0←t2*t0
    Fp12_frobenius_map_p2(&t2,&t4);         //t2←t4^(p^2)
    Fp12_mul(&t0,&t2,&t0);              //t0←t2*t0
    Fp12_frobenius_map_p6(&t2,A);            //t2←f^(-1)
    Fp12_mul(&t2,&t2,&t3);              //t2←t2*t3
    Fp12_frobenius_map_p3(&t2,&t2);         //t2←t2^(p^3)
    Fp12_mul(ANS,&t2,&t0);              //t0←t2*t0
    
    gettimeofday(&tv_end,NULL);
    FINALEXP_OPT_HARD=timedifference_msec(tv_start,tv_end);
    //Print_mpz_Cost(&mpz_cost,"Final Exp. (hard part)\n");
    
    Fp12_clear(&t0);
    Fp12_clear(&t1);
    Fp12_clear(&t2);
    Fp12_clear(&t3);
    Fp12_clear(&t4);
}

/*----------------------------------------------------------------------------*/
//pairing
void Plain_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    //Miller's Algo.
    //Init_mpz_Cost(&mpz_cost);
    gettimeofday(&tv_start,NULL);
    Miller_algo_for_plain_ate(ANS,P,Q);
    gettimeofday(&tv_end,NULL);
    MILLER_PLAINATE=timedifference_msec(tv_start,tv_end);
    //Print_mpz_Cost(&mpz_cost,"Miller's Algo.\n");
    
    //Final Exp.
    Final_exp_optimal(ANS,ANS);
}

void Opt_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    //Miller's Algo.
    //Init_mpz_Cost(&mpz_cost);
    gettimeofday(&tv_start,NULL);
    Miller_algo_for_opt_ate(ANS,P,Q);
    gettimeofday(&tv_end,NULL);
    MILLER_OPTATE=timedifference_msec(tv_start,tv_end);
    //Print_mpz_Cost(&mpz_cost,"Miller's Algo.\n");
    
    //Final Exp.
    Final_exp_optimal(ANS,ANS);
}

void X_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    //Miller's Algo.
    //Init_mpz_Cost(&mpz_cost);
    gettimeofday(&tv_start,NULL);
    Miller_algo_for_x_ate(ANS,P,Q);
    gettimeofday(&tv_end,NULL);
    MILLER_XATE=timedifference_msec(tv_start,tv_end);
    //Print_mpz_Cost(&mpz_cost,"Miller's Algo.\n");
    
    Final_exp_optimal(ANS,ANS);
}

/*----------------------------------------------------------------------------*/
//JSF
void Joint_sparse_form(int **binary,mpz_t scalar[2],int *loop_length){
	int i,j;
	unsigned long int u;
	mpz_t mod_2,mod_4,mod_8;
	mpz_init(mod_2);
	mpz_init(mod_4);
	mpz_init(mod_8);
	
	mpz_t k[2];
	mpz_init(k[0]);
	mpz_init(k[1]);
	//set
	j=0;
	mpz_set(k[0],scalar[0]);
	mpz_set(k[1],scalar[1]);
	
	while(mpz_cmp_ui(k[0],0)>0 || mpz_cmp_ui(k[1],0)>0){
		for(i=0; i<2; i++){
			mpz_mod_ui(mod_2,k[i],2);
			if(mpz_cmp_ui(mod_2,0)==0){
				u=0;
			}else{
				mpz_mod_ui(mod_4,k[i],4);
				u=mpz_get_ui(mod_4);
				if(u==3){
					u=-1;
				}
				mpz_mod_ui(mod_8,k[i],8);
				mpz_mod_ui(mod_4,k[1-i],4);
				if((mpz_cmp_ui(mod_8,3)==0 || mpz_cmp_ui(mod_8,5)==0) && mpz_cmp_ui(mod_4,2)==0){
					u=-u;
				}
			}
			binary[i][j]=u;
		}
		for(i=0; i<2; i++){
			u=binary[i][j];
			switch (u){
				case 1:
					mpz_sub_ui(k[i],k[i],1);
					break;
				case -1:
					mpz_add_ui(k[i],k[i],1);
					break;
				default:
					break;
			}
			mpz_tdiv_q_ui(k[i],k[i],2);
		}
		j=j+1;
	}
	*loop_length=j-1;
	
	mpz_clear(mod_2);
	mpz_clear(mod_4);
	mpz_clear(mod_8);
	mpz_clear(k[0]);
	mpz_clear(k[1]);	
}

/*----------------------------------------------------------------------------*/
//G1 SCM
void EFp12_G1_SCM_plain(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    //Init_mpz_Cost(&mpz_cost);
    gettimeofday(&tv_start,NULL);
    
    EFp tmp_P;
	EFp_init(&tmp_P);
	
	EFp12_to_EFp(&tmp_P,P);
	EFp_SCM(&tmp_P,&tmp_P,scalar);
	EFp_to_EFp12(ANS,&tmp_P);
	
	EFp_clear(&tmp_P);
	
	gettimeofday(&tv_end,NULL);
	G1SCM_PLAIN=timedifference_msec(tv_start,tv_end);
	//Print_mpz_Cost(&mpz_cost,"G1 SCM (plain)\n");
}

void EFp12_G1_SCM_2split(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    //Init_mpz_Cost(&mpz_cost);
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[2],loop_length;
	EFp next_tmp_P,tmp_P,skew_P;
	EFp_init(&next_tmp_P);
	EFp_init(&tmp_P);
	EFp_init(&skew_P);
	mpz_t s[2],buf,V1,V2,s1,s2,s3,s4,s5,CHECK;
	mpz_init(buf);
	mpz_init(V1);
	mpz_init(V2);
	mpz_init(s1);
	mpz_init(s2);
	mpz_init(s3);
	mpz_init(s4);
	mpz_init(s5);
	mpz_init(CHECK);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	EFp table[4];
	for(i=0; i<4; i++){
		EFp_init(&table[i]);
	}
	
	//set V1
	mpz_mul(V1,X,X);
	mpz_mul_ui(V1,V1,6);
	mpz_mul_ui(buf,X,4);
	mpz_add(V1,V1,buf);
	mpz_add_ui(V1,V1,1);
	//set V2
	mpz_add(V2,X,X);
	mpz_add_ui(V2,V2,1);
	mpz_tdiv_qr(s1,s2,scalar,V1);	//s1,s2
	mpz_mul(buf,V2,s1);		//s3,s4
	mpz_tdiv_qr(s3,s4,buf,V1);
	mpz_mul(s5,V2,s3);		//s5
	mpz_add(s[1],s4,s5);			//s[1]
	mpz_mod(s[1],s[1],order);
	mpz_sub(s[0],s2,s5);			//s[0]
	mpz_mod(s[0],s[0],order);
	//set CHECK
	mpz_sub_ui(CHECK,order,1);
	mpz_tdiv_q_ui(CHECK,CHECK,2);
	
	//set
	EFp12_to_EFp(&tmp_P,P);				//tmp_P
	EFp_skew_frobenius_map_p2(&skew_P,&tmp_P);	//skew_P
	if(mpz_cmp(s[0],CHECK)>0){
		mpz_sub(s[0],order,s[0]);
		EFp_set_neg(&tmp_P,&tmp_P);
	}
	//set table
	table[0].infinity=1;					//00
	EFp_set(&table[1],&tmp_P);			//01
	EFp_set(&table[2],&skew_P);			//10
	EFp_ECA(&table[3],&tmp_P,&skew_P);		//11
	
	//binary
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	//set binary
	char binary_s[2][loop_length+1];
	char str[5],*e;
	int binary[loop_length+1];
	for(i=0; i<2; i++){
		if(length_s[i]==loop_length){
			mpz_get_str(binary_s[i],2,s[i]);
		}else{
			char binary_buf[loop_length+1];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c",binary_s[1][i],binary_s[0][i]);
		binary[i]=(int)strtol(str,&e,2);
	}
	EFp_set(&next_tmp_P,&table[binary[0]]);
	
	//SCM
	for(i=1; i<loop_length; i++){
		EFp_ECD(&next_tmp_P,&next_tmp_P);
		EFp_ECA(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
	}
	
	EFp_to_EFp12(ANS,&next_tmp_P);
	
	mpz_clear(buf);
	mpz_clear(V1);
	mpz_clear(V2);
	mpz_clear(s1);
	mpz_clear(s2);
	mpz_clear(s3);
	mpz_clear(s4);
	mpz_clear(s5);
	mpz_clear(CHECK);
	EFp_clear(&next_tmp_P);
	EFp_clear(&tmp_P);
	EFp_clear(&skew_P);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<4; i++){
		EFp_clear(&table[i]);
	}
	
	gettimeofday(&tv_end,NULL);
	G1SCM_2SPLIT=timedifference_msec(tv_start,tv_end);
	//Print_mpz_Cost(&mpz_cost,"G1 SCM (2split)\n");
}

void EFp12_G1_SCM_2split_JSF(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    //Init_mpz_Cost(&mpz_cost);
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[2],loop_length;
	EFp next_tmp_P,tmp_P,tmp_P_neg,skew_P,skew_P_neg;
	EFp_init(&next_tmp_P);
	EFp_init(&tmp_P);
	EFp_init(&tmp_P_neg);
	EFp_init(&skew_P);
	EFp_init(&skew_P_neg);
	mpz_t s[2],buf,V1,V2,s1,s2,s3,s4,s5,CHECK;
	mpz_init(buf);
	mpz_init(V1);
	mpz_init(V2);
	mpz_init(s1);
	mpz_init(s2);
	mpz_init(s3);
	mpz_init(s4);
	mpz_init(s5);
	mpz_init(CHECK);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	EFp table[9];
	for(i=0; i<9; i++){
		EFp_init(&table[i]);
	}
	
	//set V1
	mpz_mul(V1,X,X);
	mpz_mul_ui(V1,V1,6);
	mpz_mul_ui(buf,X,4);
	mpz_add(V1,V1,buf);
	mpz_add_ui(V1,V1,1);
	//set V2
	mpz_add(V2,X,X);
	mpz_add_ui(V2,V2,1);
	mpz_tdiv_qr(s1,s2,scalar,V1);	//s1,s2
	mpz_mul(buf,V2,s1);		//s3,s4
	mpz_tdiv_qr(s3,s4,buf,V1);
	mpz_mul(s5,V2,s3);		//s5
	mpz_add(s[1],s4,s5);			//s[1]
	mpz_mod(s[1],s[1],order);
	mpz_sub(s[0],s2,s5);			//s[0]
	mpz_mod(s[0],s[0],order);
	//set CHECK
	mpz_sub_ui(CHECK,order,1);
	mpz_tdiv_q_ui(CHECK,CHECK,2);
	
	//set
	EFp12_to_EFp(&tmp_P,P);					//tmp_P
	EFp_set_neg(&tmp_P_neg,&tmp_P);			//tmp_P_neg
	EFp_skew_frobenius_map_p2(&skew_P,&tmp_P);		//skew_P
	EFp_set_neg(&skew_P_neg,&skew_P);			//skew_P_neg
	
	if(mpz_cmp(s[0],CHECK)>0){
		mpz_sub(s[0],order,s[0]);
		EFp_set_neg(&tmp_P,&tmp_P);
		EFp_set_neg(&tmp_P_neg,&tmp_P_neg);
	}
	
	//set table
	table[0].infinity=1;						//00
	EFp_set(&table[1],&tmp_P);				//01
	EFp_set(&table[2],&skew_P);				//10
	EFp_ECA(&table[3],&skew_P,&tmp_P);			//11
	EFp_set(&table[4],&tmp_P_neg);			//0-1
	EFp_set(&table[5],&skew_P_neg);			//-10
	EFp_ECA(&table[6],&skew_P_neg,&tmp_P_neg);	//-1-1
	EFp_ECA(&table[7],&skew_P,&tmp_P_neg);		//1-1
	EFp_ECA(&table[8],&skew_P_neg,&tmp_P);		//-11
	
	//get loop_length
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	//JSF
	int JSF_length;
	int JSF_binary[2][loop_length+1];
	for(i=0; i<loop_length; i++){
		JSF_binary[0][i]=0;
		JSF_binary[1][i]=0;
	}
	int *JSF_pointer[2];
	JSF_pointer[0]=JSF_binary[0];
	JSF_pointer[1]=JSF_binary[1];
	Joint_sparse_form(JSF_pointer,s,&JSF_length);
	int binary[JSF_length+1];
	for(i=JSF_length; i>=0; i--){
		if(JSF_binary[1][i]==0 && JSF_binary[0][i]==0) 		binary[i]=0;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==1) 	binary[i]=1;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==0) 	binary[i]=2;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==1)	binary[i]=3;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==-1)	binary[i]=4;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==0)	binary[i]=5;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==-1)	binary[i]=6;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==-1)	binary[i]=7;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==1)	binary[i]=8;
	}
	EFp_set(&next_tmp_P,&table[binary[JSF_length]]);
	//SCM
	for(i=JSF_length-1; i>=0; i--){
		EFp_ECD(&next_tmp_P,&next_tmp_P);
		EFp_ECA(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
	}
	EFp_to_EFp12(ANS,&next_tmp_P);
	
	mpz_clear(buf);
	mpz_clear(V1);
	mpz_clear(V2);
	mpz_clear(s1);
	mpz_clear(s2);
	mpz_clear(s3);
	mpz_clear(s4);
	mpz_clear(s5);
	mpz_clear(CHECK);
	EFp_clear(&next_tmp_P);
	EFp_clear(&tmp_P);
	EFp_clear(&tmp_P_neg);
	EFp_clear(&skew_P);
	EFp_clear(&skew_P_neg);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<9; i++){
		EFp_clear(&table[i]);
	}
	
	gettimeofday(&tv_end,NULL);
	G1SCM_2SPLIT_JSF=timedifference_msec(tv_start,tv_end);
	//Print_mpz_Cost(&mpz_cost,"G1 SCM (2split-JSF)\n");
}

/*----------------------------------------------------------------------------*/
//G2 SCM
void EFp12_G2_SCM_plain(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
    //Init_mpz_Cost(&mpz_cost);
    gettimeofday(&tv_start,NULL);
    
    EFp2 tmp_Q;
	EFp2_init(&tmp_Q);
	
	EFp12_to_EFp2(&tmp_Q,Q);
	EFp2_SCM(&tmp_Q,&tmp_Q,scalar);
	EFp2_to_EFp12(ANS,&tmp_Q);
	
	EFp2_clear(&tmp_Q);
	
	gettimeofday(&tv_end,NULL);
	G2SCM_PLAIN=timedifference_msec(tv_start,tv_end);
	//Print_mpz_Cost(&mpz_cost,"G2 SCM (plain)\n");
}

void EFp12_G2_SCM_2split(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
    //Init_mpz_Cost(&mpz_cost);
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[2],loop_length;
	EFp2 next_twisted_Q,twisted_Q,skew_Q;
	EFp2_init(&next_twisted_Q);
	EFp2_init(&twisted_Q);
	EFp2_init(&skew_Q);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	EFp2 table[4];
	for(i=0; i<4; i++){
		EFp2_init(&table[i]);
	}
	
	//set
	EFp12_to_EFp2(&twisted_Q,Q);				//twisted_Q
	EFp2_skew_frobenius_map_p1(&skew_Q,&twisted_Q);//skew_Q
	
	//set table
	table[0].infinity=1;						//00
	EFp2_set(&table[1],&twisted_Q);			//01
	EFp2_set(&table[2],&skew_Q);				//10
	EFp2_ECA(&table[3],&twisted_Q,&skew_Q);		//11
	
	//s0,s1
	mpz_sub_ui(buf,trace,1);
	mpz_tdiv_qr(s[1],s[0],scalar,buf);
	
	//binary
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	//set binary
	char binary_s[2][loop_length+1];
	char str[5],*e;
	int binary[loop_length+1];
	for(i=0; i<2; i++){
		if(length_s[i]==loop_length){
			mpz_get_str(binary_s[i],2,s[i]);
		}else{
			char binary_buf[loop_length+1];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c",binary_s[1][i],binary_s[0][i]);
		binary[i]=(int)strtol(str,&e,2);
	}
	EFp2_set(&next_twisted_Q,&table[binary[0]]);
	
	//SCM
	for(i=1; i<loop_length; i++){
		EFp2_ECD(&next_twisted_Q,&next_twisted_Q);
		EFp2_ECA(&next_twisted_Q,&next_twisted_Q,&table[binary[i]]);
	}
	
	EFp2_to_EFp12(ANS,&next_twisted_Q);
	ANS->infinity=next_twisted_Q.infinity;
	
	mpz_clear(buf);
	EFp2_clear(&next_twisted_Q);
	EFp2_clear(&twisted_Q);
	EFp2_clear(&skew_Q);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<4; i++){
		EFp2_clear(&table[i]);
	}
	
	gettimeofday(&tv_end,NULL);
	G2SCM_2SPLIT=timedifference_msec(tv_start,tv_end);
	//Print_mpz_Cost(&mpz_cost,"G2 SCM (2split)\n");
}

void EFp12_G2_SCM_2split_JSF(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
    //Init_mpz_Cost(&mpz_cost);
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[2],loop_length;
	EFp2 next_tmp_Q,tmp_Q,tmp_Q_neg,skew_Q,skew_Q_neg;
	EFp2_init(&next_tmp_Q);
	EFp2_init(&tmp_Q);
	EFp2_init(&tmp_Q_neg);
	EFp2_init(&skew_Q);
	EFp2_init(&skew_Q_neg);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	EFp2 table[9];
	for(i=0; i<9; i++){
		EFp2_init(&table[i]);
	}
	
	//set
	EFp12_to_EFp2(&tmp_Q,Q);					//tmp_Q
	EFp2_set_neg(&tmp_Q_neg,&tmp_Q);			//tmp_Q_neg
	EFp2_skew_frobenius_map_p1(&skew_Q,&tmp_Q);		//skew_Q
	EFp2_set_neg(&skew_Q_neg,&skew_Q);			//skew_Q_neg
	
	//set table
	table[0].infinity=1;						//00
	EFp2_set(&table[1],&tmp_Q);				//01
	EFp2_set(&table[2],&skew_Q);				//10
	EFp2_ECA(&table[3],&skew_Q,&tmp_Q);		//11
	EFp2_set(&table[4],&tmp_Q_neg);			//0-1
	EFp2_set(&table[5],&skew_Q_neg);			//-10
	EFp2_ECA(&table[6],&skew_Q_neg,&tmp_Q_neg);	//-1-1
	EFp2_ECA(&table[7],&skew_Q,&tmp_Q_neg);		//1-1
	EFp2_ECA(&table[8],&skew_Q_neg,&tmp_Q);		//-11
	
	//s0,s1
	mpz_sub_ui(buf,trace,1);
	mpz_tdiv_qr(s[1],s[0],scalar,buf);
	
	//get loop_length
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	//JSF
	int JSF_length;
	int JSF_binary[2][loop_length+1];
	for(i=0; i<loop_length; i++){
		JSF_binary[0][i]=0;
		JSF_binary[1][i]=0;
	}
	int *JSF_pointer[2];
	JSF_pointer[0]=JSF_binary[0];
	JSF_pointer[1]=JSF_binary[1];
	Joint_sparse_form(JSF_pointer,s,&JSF_length);
	int binary[JSF_length+1];
	for(i=JSF_length; i>=0; i--){
		if(JSF_binary[1][i]==0 && JSF_binary[0][i]==0) 		binary[i]=0;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==1) 	binary[i]=1;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==0) 	binary[i]=2;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==1)	binary[i]=3;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==-1)	binary[i]=4;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==0)	binary[i]=5;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==-1)	binary[i]=6;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==-1)	binary[i]=7;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==1)	binary[i]=8;
	}
	EFp2_set(&next_tmp_Q,&table[binary[JSF_length]]);
	//SCM
	for(i=JSF_length-1; i>=0; i--){
		EFp2_ECD(&next_tmp_Q,&next_tmp_Q);
		EFp2_ECA(&next_tmp_Q,&next_tmp_Q,&table[binary[i]]);
	}
	EFp2_to_EFp12(ANS,&next_tmp_Q);
	
	mpz_clear(buf);
	EFp2_clear(&next_tmp_Q);
	EFp2_clear(&tmp_Q);
	EFp2_clear(&tmp_Q_neg);
	EFp2_clear(&skew_Q);
	EFp2_clear(&skew_Q_neg);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<9; i++){
		EFp2_clear(&table[i]);
	}
	
	gettimeofday(&tv_end,NULL);
	G2SCM_2SPLIT_JSF=timedifference_msec(tv_start,tv_end);
	//Print_mpz_Cost(&mpz_cost,"G2 SCM (2split-JSF)\n");
}

void EFp12_G2_SCM_4split(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
    //Init_mpz_Cost(&mpz_cost);
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[4],loop_length;
	EFp2 next_twisted_Q,twisted_Q,twisted_Q_6x,twisted_Q_6xx,twisted_Q_36xxx,skew_Q,skew_Q_neg,skew_Q_puls1,minus_skew_Q_puls1;
	EFp2_init(&next_twisted_Q);
	EFp2_init(&twisted_Q);
	EFp2_init(&twisted_Q_6x);
	EFp2_init(&twisted_Q_6xx);
	EFp2_init(&twisted_Q_36xxx);
	EFp2_init(&skew_Q);
	EFp2_init(&skew_Q_neg);
	EFp2_init(&skew_Q_puls1);
	EFp2_init(&minus_skew_Q_puls1);
	
	mpz_t buf,A,B,s[4];
	mpz_init(buf);
	mpz_init(A);
	mpz_init(B);
	for(i=0; i<4; i++){
		mpz_init(s[i]);
	}
	//table
	EFp2 table[16];
	for(i=0; i<16; i++){
		EFp2_init(&table[i]);
	}
	
	//twisted_Q
	EFp12_to_EFp2(&twisted_Q,Q);
	//twisted_Q_6xx
	EFp2_skew_frobenius_map_p1(&skew_Q,&twisted_Q);
	EFp2_set_neg(&skew_Q_neg,&skew_Q);
	EFp2_set(&twisted_Q_6xx,&skew_Q);
	//twisted_Q_6x
	EFp2_ECA(&skew_Q_puls1,&skew_Q,&twisted_Q);
	EFp2_ECA(&minus_skew_Q_puls1,&skew_Q_neg,&twisted_Q);
	EFp2_skew_frobenius_map_p3(&minus_skew_Q_puls1,&minus_skew_Q_puls1);
	EFp2_ECA(&twisted_Q_6x,&skew_Q_puls1,&minus_skew_Q_puls1);
	EFp2_set_neg(&twisted_Q_6x,&twisted_Q_6x);
	//twisted_Q_36xxx
	EFp2_skew_frobenius_map_p1(&twisted_Q_36xxx,&twisted_Q_6x);
	
	//set table
	table[0].infinity=1;								//0000
	EFp2_set(&table[1],&twisted_Q);					//0001
	EFp2_set(&table[2],&twisted_Q_6x);					//0010
	EFp2_ECA(&table[3],&twisted_Q_6x,&twisted_Q);		//0011
	EFp2_set(&table[4],&twisted_Q_6xx);				//0100
	EFp2_ECA(&table[5],&twisted_Q_6xx,&twisted_Q);		//0101
	EFp2_ECA(&table[6],&twisted_Q_6xx,&twisted_Q_6x);		//0110
	EFp2_ECA(&table[7],&table[6],&twisted_Q);			//0111
	EFp2_set(&table[8],&twisted_Q_36xxx);				//1000
	EFp2_ECA(&table[9],&twisted_Q_36xxx,&twisted_Q);		//1001
	EFp2_ECA(&table[10],&twisted_Q_36xxx,&twisted_Q_6x);	//1010
	EFp2_ECA(&table[11],&twisted_Q_36xxx,&table[3]);		//1011
	EFp2_ECA(&table[12],&twisted_Q_36xxx,&twisted_Q_6xx);	//1100
	EFp2_ECA(&table[13],&table[12],&twisted_Q);			//1101
	EFp2_ECA(&table[14],&table[12],&twisted_Q_6x);		//1110
	EFp2_ECA(&table[15],&table[14],&twisted_Q);			//1111
	
	//set
	//s0,s1,s2,s3
	mpz_sub_ui(buf,trace,1);
	mpz_tdiv_qr(B,A,scalar,buf);
	mpz_mul_ui(buf,X,6);
	mpz_tdiv_qr(s[1],s[0],A,buf);
	mpz_tdiv_qr(s[3],s[2],B,buf);
	
	//binary
	loop_length=0;
	for(i=0; i<4; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	//set binary
	char binary_s[4][loop_length+1];
	char str[5],*e;
	int binary[loop_length+1];
	for(i=0; i<4; i++){
		if(length_s[i]==loop_length){
			mpz_get_str(binary_s[i],2,s[i]);
		}else{
			char binary_buf[loop_length+1];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c%c%c",binary_s[3][i],binary_s[2][i],binary_s[1][i],binary_s[0][i]);
		binary[i]=(int)strtol(str,&e,2);
	}
	
	EFp2_set(&next_twisted_Q,&table[binary[0]]);
	
	//SCM
	for(i=1; i<loop_length; i++){
		EFp2_ECD(&next_twisted_Q,&next_twisted_Q);
		EFp2_ECA(&next_twisted_Q,&next_twisted_Q,&table[binary[i]]);
	}
	
	EFp2_to_EFp12(ANS,&next_twisted_Q);
	
	mpz_clear(buf);
	mpz_clear(A);
	mpz_clear(B);
	EFp2_clear(&next_twisted_Q);
	EFp2_clear(&twisted_Q);
	EFp2_clear(&twisted_Q_6x);
	EFp2_clear(&twisted_Q_6xx);
	EFp2_clear(&twisted_Q_36xxx);
	EFp2_clear(&skew_Q);
	EFp2_clear(&skew_Q_neg);
	EFp2_clear(&skew_Q_puls1);
	EFp2_clear(&minus_skew_Q_puls1);
	for(i=0; i<4; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<16; i++){
		EFp2_clear(&table[i]);
	}
	
	gettimeofday(&tv_end,NULL);
	G2SCM_4SPLIT=timedifference_msec(tv_start,tv_end);
	//Print_mpz_Cost(&mpz_cost,"G2 SCM (4split)\n");
}

/*----------------------------------------------------------------------------*/
//G3 EXP
void Fp12_G3_EXP_plain(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    //Init_mpz_Cost(&mpz_cost);
    gettimeofday(&tv_start,NULL);
    
    int i,length;
	length=(int)mpz_sizeinbase(scalar,2);
	char binary[length];
	mpz_get_str(binary,2,scalar);
	Fp12 buf;
	Fp12_init(&buf);
	Fp12_set(&buf,A);
	
	for(i=1; binary[i]!='\0'; i++){
		Fp12_sqr_cyclotomic(&buf,&buf);
		if(binary[i]=='1'){
			Fp12_mul(&buf,A,&buf);
		}
	}
	
	Fp12_set(ANS,&buf);
	Fp12_clear(&buf);
	
	gettimeofday(&tv_end,NULL);
	G3EXP_PLAIN=timedifference_msec(tv_start,tv_end);
	//Print_mpz_Cost(&mpz_cost,"G3 EXP (plain)\n");
}

void Fp12_G3_EXP_2split(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    //Init_mpz_Cost(&mpz_cost);
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[2],loop_length;
	Fp12 Buf,next_f,f,frobenius_f;
	Fp12_init(&Buf);
	Fp12_init(&next_f);
	Fp12_init(&f);
	Fp12_init(&frobenius_f);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	Fp12 table[4];
	for(i=0; i<4; i++){
		Fp12_init(&table[i]);
	}
	
	//set
	Fp12_set(&f,A);						//f
	Fp12_frobenius_map_p1(&frobenius_f,&f);			//frobenius_f
	
	//set table
	Fp_set_ui(&table[0].x0.x0.x0,1);			//00
	Fp12_set(&table[1],&f);					//01
	Fp12_set(&table[2],&frobenius_f);			//10
	Fp12_mul(&table[3],&f,&frobenius_f);		//11
	
	//s0,s1
	mpz_sub_ui(buf,trace,1);
	mpz_tdiv_qr(s[1],s[0],scalar,buf);
	
	//binary
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	//set binary
	char binary_s[2][loop_length+1];
	char str[5],*e;
	int binary[loop_length+1];
	for(i=0; i<2; i++){
		if(length_s[i]==loop_length){
			mpz_get_str(binary_s[i],2,s[i]);
		}else{
			char binary_buf[loop_length+1];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c",binary_s[1][i],binary_s[0][i]);
		binary[i]=(int)strtol(str,&e,2);
	}
	Fp12_set(&next_f,&table[binary[0]]);
	
	//EXP
	for(i=1; i<loop_length; i++){
		Fp12_sqr_cyclotomic(&next_f,&next_f);
		Fp12_mul(&next_f,&next_f,&table[binary[i]]);
	}
	
	Fp12_set(ANS,&next_f);
	
	mpz_clear(buf);
	Fp12_clear(&Buf);
	Fp12_clear(&next_f);
	Fp12_clear(&f);
	Fp12_clear(&frobenius_f);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<4; i++){
		Fp12_clear(&table[i]);
	}
	
	gettimeofday(&tv_end,NULL);
	G3EXP_2SPLIT=timedifference_msec(tv_start,tv_end);
	//Print_mpz_Cost(&mpz_cost,"G3 EXP (2split)\n");
}

void Fp12_G3_EXP_2split_JSF(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    //Init_mpz_Cost(&mpz_cost);
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[2],loop_length;
	Fp12 next_f,f,f_inv,frobenius_f,frobenius_f_inv;
	Fp12_init(&next_f);
	Fp12_init(&f);
	Fp12_init(&f_inv);
	Fp12_init(&frobenius_f);
	Fp12_init(&frobenius_f_inv);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	Fp12 table[9];
	for(i=0; i<9; i++){
		Fp12_init(&table[i]);
	}
	
	//set
	Fp12_set(&f,A);							//f
	Fp12_frobenius_map_p6(&f_inv,&f);						//f_inv
	Fp12_frobenius_map_p1(&frobenius_f,&f);				//frobenius_f
	Fp12_frobenius_map_p6(&frobenius_f_inv,&frobenius_f);		//frobenius_f_inv
	
	//set table
	Fp_set_ui(&table[0].x0.x0.x0,1);			//00
	Fp12_set(&table[1],&f);					//01
	Fp12_set(&table[2],&frobenius_f);			//10
	Fp12_mul(&table[3],&frobenius_f,&f);		//11
	Fp12_set(&table[4],&f_inv);				//0-1
	Fp12_set(&table[5],&frobenius_f_inv);		//-10
	Fp12_mul(&table[6],&frobenius_f_inv,&f_inv);	//-1-1
	Fp12_mul(&table[7],&frobenius_f,&f_inv);	//1-1
	Fp12_mul(&table[8],&frobenius_f_inv,&f);	//-11
	
	//s0,s1
	mpz_sub_ui(buf,trace,1);
	mpz_tdiv_qr(s[1],s[0],scalar,buf);
	
	//get loop_length
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	//JSF
	int JSF_length;
	int JSF_binary[2][loop_length+1];
	for(i=0; i<loop_length; i++){
		JSF_binary[0][i]=0;
		JSF_binary[1][i]=0;
	}
	int *JSF_pointer[2];
	JSF_pointer[0]=JSF_binary[0];
	JSF_pointer[1]=JSF_binary[1];
	Joint_sparse_form(JSF_pointer,s,&JSF_length);
	int binary[JSF_length+1];
	for(i=JSF_length; i>=0; i--){
		if(JSF_binary[1][i]==0 && JSF_binary[0][i]==0) 		binary[i]=0;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==1) 	binary[i]=1;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==0) 	binary[i]=2;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==1)	binary[i]=3;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==-1)	binary[i]=4;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==0)	binary[i]=5;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==-1)	binary[i]=6;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==-1)	binary[i]=7;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==1)	binary[i]=8;
	}
	Fp12_set(&next_f,&table[binary[JSF_length]]);
	//SCM
	for(i=JSF_length-1; i>=0; i--){
		Fp12_sqr_cyclotomic(&next_f,&next_f);
		Fp12_mul(&next_f,&next_f,&table[binary[i]]);
	}
	Fp12_set(ANS,&next_f);
	
	mpz_clear(buf);
	Fp12_clear(&next_f);
	Fp12_clear(&f);
	Fp12_clear(&f_inv);
	Fp12_clear(&frobenius_f);
	Fp12_clear(&frobenius_f_inv);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<9; i++){
		Fp12_clear(&table[i]);
	}
	
	gettimeofday(&tv_end,NULL);
	G3EXP_2SPLIT_JSF=timedifference_msec(tv_start,tv_end);
	//Print_mpz_Cost(&mpz_cost,"G3 EXP (2split-JSF)\n");
}

void Fp12_G3_EXP_4split(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    //Init_mpz_Cost(&mpz_cost);
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[4],loop_length;
	Fp12 Buf;
	Fp12_init(&Buf);
	Fp12 next_f,f,f_6x,f_6xx,f_36xxx,frobenius_f,frobenius_f_inv,frobenius_f_f,frobenius_f_inv_f;
	Fp12_init(&next_f);
	Fp12_init(&f);
	Fp12_init(&f_6x);
	Fp12_init(&f_6xx);
	Fp12_init(&f_36xxx);
	Fp12_init(&frobenius_f);
	Fp12_init(&frobenius_f_inv);
	Fp12_init(&frobenius_f_f);
	Fp12_init(&frobenius_f_inv_f);
	mpz_t buf,a,b,s[4];
	mpz_init(buf);
	mpz_init(a);
	mpz_init(b);
	for(i=0; i<4; i++){
		mpz_init(s[i]);
	}
	//table
	Fp12 table[16];
	for(i=0; i<16; i++){
		Fp12_init(&table[i]);
	}
	
	//f
	Fp12_set(&f,A);
	//f_6xx
	Fp12_frobenius_map_p1(&frobenius_f,&f);
	Fp12_set(&f_6xx,&frobenius_f);
	Fp12_frobenius_map_p6(&frobenius_f_inv,&frobenius_f);
	//f_6x
	Fp12_mul(&frobenius_f_f,&frobenius_f,&f);
	Fp12_mul(&frobenius_f_inv_f,&frobenius_f_inv,&f);
	Fp12_frobenius_map_p3(&frobenius_f_inv_f,&frobenius_f_inv_f);
	Fp12_mul(&f_6x,&frobenius_f_f,&frobenius_f_inv_f);
	Fp12_frobenius_map_p6(&f_6x,&f_6x);
	//f_36xxx
	Fp12_frobenius_map_p1(&f_36xxx,&f_6x);
	
	//set table
	Fp_set_ui(&table[0].x0.x0.x0,1);			//0000
	Fp12_set(&table[1],&f);						//0001
	Fp12_set(&table[2],&f_6x);					//0010
	Fp12_mul(&table[3],&f_6x,&f);					//0011
	Fp12_set(&table[4],&f_6xx);					//0100
	Fp12_mul(&table[5],&f_6xx,&f);				//0101
	Fp12_mul(&table[6],&f_6xx,&f_6x);				//0110
	Fp12_mul(&table[7],&table[6],&f);				//0111
	Fp12_set(&table[8],&f_36xxx);					//1000
	Fp12_mul(&table[9],&f_36xxx,&f);				//1001
	Fp12_mul(&table[10],&f_36xxx,&f_6x);			//1010
	Fp12_mul(&table[11],&f_36xxx,&table[3]);		//1011
	Fp12_mul(&table[12],&f_36xxx,&f_6xx);			//1100
	Fp12_mul(&table[13],&table[12],&f);			//1101
	Fp12_mul(&table[14],&table[12],&f_6x);			//1110
	Fp12_mul(&table[15],&table[14],&f);			//1111
	
	//set
	//s0,s1,s2,s3
	mpz_sub_ui(buf,trace,1);
	mpz_tdiv_qr(b,a,scalar,buf);
	mpz_mul_ui(buf,X,6);
	mpz_tdiv_qr(s[1],s[0],a,buf);
	mpz_tdiv_qr(s[3],s[2],b,buf);
	
	//binary
	loop_length=0;
	for(i=0; i<4; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	//set binary
	char binary_s[4][loop_length+1];
	char str[5],*e;
	int binary[loop_length+1];
	for(i=0; i<4; i++){
		if(length_s[i]==loop_length){
			mpz_get_str(binary_s[i],2,s[i]);
		}else{
			char binary_buf[loop_length+1];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c%c%c",binary_s[3][i],binary_s[2][i],binary_s[1][i],binary_s[0][i]);
		binary[i]=(int)strtol(str,&e,2);
	}
	
	Fp12_set(&next_f,&table[binary[0]]);
	
	//SCM
	for(i=1; i<loop_length; i++){
		Fp12_sqr_cyclotomic(&next_f,&next_f);
		Fp12_mul(&next_f,&next_f,&table[binary[i]]);
	}
	
	Fp12_set(ANS,&next_f);
	
	mpz_clear(buf);
	mpz_clear(a);
	mpz_clear(b);
	Fp12_clear(&Buf);
	Fp12_clear(&next_f);
	Fp12_clear(&f);
	Fp12_clear(&f_6x);
	Fp12_clear(&f_6xx);
	Fp12_clear(&f_36xxx);
	Fp12_clear(&frobenius_f);
	Fp12_clear(&frobenius_f_inv);
	Fp12_clear(&frobenius_f_f);
	Fp12_clear(&frobenius_f_inv_f);
	for(i=0; i<4; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<16; i++){
		Fp12_clear(&table[i]);
	}
	
	gettimeofday(&tv_end,NULL);
	G3EXP_4SPLIT=timedifference_msec(tv_start,tv_end);
	//Print_mpz_Cost(&mpz_cost,"G3 EXP (4split)\n");
}

/*----------------------------------------------------------------------------*/
//init/set/clear
void BN12_init(){
    init_parameters();
    generate_X();
    if(generate_prime()==1 && generate_order()==1){
        generate_trace();
        weil();
        get_epsilon();
        get_Two_inv();
        set_basis();
        set_frobenius_constant();
        set_curve_parameter();
    }else{
        BN12_clear();
        printf("error : prime\nexit\n");
    }
}

void init_parameters(){
    int i,j;
    
    mpz_init(X);
    mpz_init(prime);
    mpz_init(order);
    mpz_init(trace);
    
    mpz_init(EFp_total);
    mpz_init(EFp12_total);
    mpz_init(curve_b);
    
    for(i=0; i<X_length+1; i++){
        X_binary[i]=0;
    }
    for(i=0; i<X6_2_length+1; i++){
        X6_2_binary[i]=0;
    }
    
    mpz_init(epsilon1);
    mpz_init(epsilon2);
    mpz_init(Two_inv);
    Fp2_init(&Alpha_1);
    Fp2_init(&Alpha_1_inv);
    
    for(i=0; i<12; i++){
        for(j=0; j<6; j++){
            Fp2_init(&frobenius_constant[i][j]);
        }
        for(j=0; j<2; j++){
            Fp2_init(&skew_frobenius_constant[i][j]);
        }
    }
}

void generate_X(){
    int i;
    mpz_t buf;
    mpz_init(buf);
    
    //X_binary
    X_binary[114]=1;
    X_binary[84]=1;
    X_binary[53]=-1;
    X_binary[0]=-1;
    
    //X_binary_opt
    X6_2_binary[116]=1;
    X6_2_binary[115]=1;
    X6_2_binary[86]=1;
    X6_2_binary[85]=1;
    X6_2_binary[55]=-1;
    X6_2_binary[54]=-1;
    X6_2_binary[2]=-1;
    
    //BN12.X
    mpz_set_ui(X,0);
    for(i=X_length; i>=0; i--){
        if(X_binary[i]==1){
            mpz_ui_pow_ui(buf,2,i);
            mpz_add(X,X,buf);
        }else if(X_binary[i]==-1){
            mpz_ui_pow_ui(buf,2,i);
            mpz_sub(X,X,buf);
        }
    }
    
    mpz_clear(buf);
}

int generate_prime(){
    mpz_t buf,result;
    mpz_init(buf);
    mpz_init(result);
    
    //prime
    mpz_pow_ui(buf,X,4);
    mpz_mul_ui(buf,buf,36);
    mpz_set(result,buf);
    mpz_pow_ui(buf,X,3);
    mpz_mul_ui(buf,buf,36);
    mpz_add(result,result,buf);
    mpz_pow_ui(buf,X,2);
    mpz_mul_ui(buf,buf,24);
    mpz_add(result,result,buf);
    mpz_mul_ui(buf,X,6);
    mpz_add(result,result,buf);
    mpz_add_ui(result,result,1);
    
    //isprime
    if(mpz_probab_prime_p(result,25)==0){
        mpz_clear(buf);
        mpz_clear(result);
        return 0;
    }else{
        mpz_set(prime,result);
        mpz_clear(buf);
        mpz_clear(result);
        return 1;
    }
}

int generate_order(){
    mpz_t buf,result;
    mpz_init(buf);
    mpz_init(result);
    
    //prime
    mpz_pow_ui(buf,X,4);
    mpz_mul_ui(buf,buf,36);
    mpz_set(result,buf);
    mpz_pow_ui(buf,X,3);
    mpz_mul_ui(buf,buf,36);
    mpz_add(result,result,buf);
    mpz_pow_ui(buf,X,2);
    mpz_mul_ui(buf,buf,18);
    mpz_add(result,result,buf);
    mpz_mul_ui(buf,X,6);
    mpz_add(result,result,buf);
    mpz_add_ui(result,result,1);
    
    //isprime
    if(mpz_probab_prime_p(result,25)==0){
        mpz_clear(buf);
        mpz_clear(result);
        return 0;
    }else{
        mpz_set(order,result);
        mpz_clear(buf);
        mpz_clear(result);
        return 1;
    }
}

void generate_trace(){
    mpz_t buf;
    mpz_init(buf);
    
    mpz_pow_ui(buf,X,2);
    mpz_mul_ui(buf,buf,6);
    mpz_add_ui(trace,buf,1);
    
    mpz_clear(buf);
}

void weil(){
    mpz_t t2,t6,t12,p2,p6,buf;
    mpz_init(t2);
    mpz_init(t6);
    mpz_init(t12);
    mpz_init(p2);
    mpz_init(p6);
    mpz_init(buf);
    
    //EFp_total
    mpz_add_ui(buf,prime,1);
    mpz_sub(EFp_total,buf,trace);
    
    //t2←α^2+β^2
    mpz_pow_ui(t2,trace,2);
    mpz_mul_ui(buf,prime,2);
    mpz_sub(t2,t2,buf);
    mpz_pow_ui(p2,prime,2);
    
    //α^6+β^6
    mpz_pow_ui(t6,t2,3);
    mpz_mul(buf,t2,p2);
    mpz_mul_ui(buf,buf,3);
    mpz_sub(t6,t6,buf);
    mpz_pow_ui(p6,p2,3);
    
    //α^12+β^12
    mpz_pow_ui(t12,t6,2);
    mpz_mul_ui(buf,p6,2);
    mpz_sub(t12,t12,buf);
    
    //EFp12_232_total
    mpz_pow_ui(buf,p6,2);
    mpz_sub(buf,buf,t12);
    mpz_add_ui(EFp12_total,buf,1);
    
    mpz_clear(t2);
    mpz_clear(t6);
    mpz_clear(t12);
    mpz_clear(p2);
    mpz_clear(p6);
    mpz_clear(buf);
}

void get_epsilon(){
    Fp inv,buf,result1,result2;
    Fp_init(&inv);
    Fp_init(&buf);
    Fp_init(&result1);
    Fp_init(&result2);
    
    Fp_set_ui(&buf,2);
    Fp_inv(&inv,&buf);
    mpz_sub_ui(buf.x0,prime,3);
    
    Fp_sqrt(&buf,&buf);
    Fp_sub_ui(&buf,&buf,1);
    Fp_mul(&result1,&buf,&inv);
    Fp_mul(&result2,&result1,&result1);
    
    mpz_set(epsilon1,result1.x0);
    mpz_set(epsilon2,result2.x0);
    
    Fp_clear(&inv);
    Fp_clear(&buf);
    Fp_clear(&result1);
    Fp_clear(&result2);
}

void get_Two_inv(){
    mpz_set_ui(Two_inv,2);
    mpz_invert(Two_inv,Two_inv,prime);
}

void set_basis(){
    Fp2_set_ui(&Alpha_1,1);
    Fp2_inv(&Alpha_1_inv,&Alpha_1);
}

void set_frobenius_constant(){
	Fp2 tmp1,tmp2,tmp3;
	Fp2_init(&tmp1);
	Fp2_init(&tmp2);
	Fp2_init(&tmp3);
	
	mpz_t exp,buf,p2,p3,p4,p6,p8,p10;
	mpz_init(exp);
	mpz_init(buf);
	mpz_init(p2);
	mpz_init(p3);
	mpz_init(p4);
	mpz_init(p6);
	mpz_init(p8);
	mpz_init(p10);
	
	mpz_mul(p2,prime,prime);
	mpz_mul(p3,p2,prime);
	mpz_mul(p4,p3,prime);
	mpz_mul(p6,p4,p2);
	mpz_mul(p8,p6,p2);
	mpz_mul(p10,p8,p2);
	
	//frobenius_1
	mpz_sub_ui(exp,prime,1);
	mpz_tdiv_q_ui(exp,exp,3);
	Fp2_pow(&tmp1,&Alpha_1,exp);
	Fp2_mul(&tmp2,&tmp1,&tmp1);
	mpz_tdiv_q_ui(exp,exp,2);
	Fp2_pow(&tmp3,&Alpha_1,exp);
	//set f_p1
	Fp_set_ui(&frobenius_constant[f_p1][0].x0,1);
	Fp2_set(&frobenius_constant[f_p1][1],&tmp1);
	Fp2_set(&frobenius_constant[f_p1][2],&tmp2);
	Fp2_set(&frobenius_constant[f_p1][3],&tmp3);
	Fp2_mul(&frobenius_constant[f_p1][4],&tmp1,&tmp3);
	Fp2_mul(&frobenius_constant[f_p1][5],&tmp2,&tmp3);
	
	//set skew_f_p1
	Fp2_inv(&tmp1,&tmp1);
	mpz_sub_ui(exp,prime,1);
	mpz_tdiv_q_ui(exp,exp,2);
	Fp2_pow(&tmp2,&Alpha_1,exp);
	Fp2_inv(&tmp2,&tmp2);
	Fp2_set(&skew_frobenius_constant[f_p1][0],&tmp1);
	Fp2_set(&skew_frobenius_constant[f_p1][1],&tmp2);
	
	//frobenius_2
	mpz_sub_ui(exp,p2,1);
	mpz_tdiv_q_ui(exp,exp,3);
	Fp2_pow(&tmp1,&Alpha_1,exp);
	Fp2_mul(&tmp2,&tmp1,&tmp1);
	mpz_tdiv_q_ui(exp,exp,2);
	Fp2_pow(&tmp3,&Alpha_1,exp);
	//set f_p2
	Fp_set_ui(&frobenius_constant[f_p2][0].x0,1);
	Fp2_set(&frobenius_constant[f_p2][1],&tmp1);
	Fp2_set(&frobenius_constant[f_p2][2],&tmp2);
	Fp2_set(&frobenius_constant[f_p2][3],&tmp3);
	Fp2_mul(&frobenius_constant[f_p2][4],&tmp1,&tmp3);
	Fp2_mul(&frobenius_constant[f_p2][5],&tmp2,&tmp3);
	//set skew_f_p2
	Fp2_inv(&tmp1,&tmp1);
	mpz_sub_ui(exp,p2,1);
	mpz_tdiv_q_ui(exp,exp,2);
	Fp2_pow(&tmp2,&Alpha_1,exp);
	Fp2_inv(&tmp2,&tmp2);
	Fp2_set(&skew_frobenius_constant[f_p2][0],&tmp1);
	Fp2_set(&skew_frobenius_constant[f_p2][1],&tmp2);
	
	//frobenius_3
	mpz_sub_ui(exp,p3,1);
	mpz_tdiv_q_ui(exp,exp,3);
	Fp2_pow(&tmp1,&Alpha_1,exp);
	Fp2_mul(&tmp2,&tmp1,&tmp1);
	mpz_tdiv_q_ui(exp,exp,2);
	Fp2_pow(&tmp3,&Alpha_1,exp);
	//set f_p3
	Fp_set_ui(&frobenius_constant[f_p3][0].x0,1);
	Fp2_set(&frobenius_constant[f_p3][1],&tmp1);
	Fp2_set(&frobenius_constant[f_p3][2],&tmp2);
	Fp2_set(&frobenius_constant[f_p3][3],&tmp3);
	Fp2_mul(&frobenius_constant[f_p3][4],&tmp1,&tmp3);
	Fp2_mul(&frobenius_constant[f_p3][5],&tmp2,&tmp3);
	//set skew_f_p3
	Fp2_inv(&tmp1,&tmp1);
	mpz_sub_ui(exp,p3,1);
	mpz_tdiv_q_ui(exp,exp,2);
	Fp2_pow(&tmp2,&Alpha_1,exp);
	Fp2_inv(&tmp2,&tmp2);
	Fp2_set(&skew_frobenius_constant[f_p3][0],&tmp1);
	Fp2_set(&skew_frobenius_constant[f_p3][1],&tmp2);
	
	//frobenius_constant[f_p4]
	mpz_sub_ui(exp,p4,1);
	mpz_tdiv_q_ui(exp,exp,3);
	Fp2_pow(&tmp1,&Alpha_1,exp);
	Fp2_mul(&tmp2,&tmp1,&tmp1);
	mpz_tdiv_q_ui(exp,exp,2);
	Fp2_pow(&tmp3,&Alpha_1,exp);
	//set frobenius_constant[f_p4]
	Fp_set_ui(&frobenius_constant[f_p4][0].x0,1);
	Fp2_set(&frobenius_constant[f_p4][1],&tmp1);
	Fp2_set(&frobenius_constant[f_p4][2],&tmp2);
	Fp2_set(&frobenius_constant[f_p4][3],&tmp3);
	Fp2_mul(&frobenius_constant[f_p4][4],&tmp1,&tmp3);
	Fp2_mul(&frobenius_constant[f_p4][5],&tmp2,&tmp3);
	
	//frobenius_constant[f_p8]
	mpz_sub_ui(exp,p8,1);
	mpz_tdiv_q_ui(exp,exp,3);
	Fp2_pow(&tmp1,&Alpha_1,exp);
	Fp2_mul(&tmp2,&tmp1,&tmp1);
	mpz_tdiv_q_ui(exp,exp,2);
	Fp2_pow(&tmp3,&Alpha_1,exp);
	//set frobenius_constant[f_p8]
	Fp_set_ui(&frobenius_constant[f_p8][0].x0,1);
	Fp2_set(&frobenius_constant[f_p8][1],&tmp1);
	Fp2_set(&frobenius_constant[f_p8][2],&tmp2);
	Fp2_set(&frobenius_constant[f_p8][3],&tmp3);
	Fp2_mul(&frobenius_constant[f_p8][4],&tmp1,&tmp3);
	Fp2_mul(&frobenius_constant[f_p8][5],&tmp2,&tmp3);
	
	//frobenius_10
	mpz_sub_ui(exp,p10,1);
	mpz_tdiv_q_ui(exp,exp,3);
	Fp2_pow(&tmp1,&Alpha_1,exp);
	Fp2_mul(&tmp2,&tmp1,&tmp1);
	mpz_tdiv_q_ui(exp,exp,2);
	Fp2_pow(&tmp3,&Alpha_1,exp);
	//set frobenius_10
	Fp_set_ui(&frobenius_constant[f_p10][0].x0,1);
	Fp2_set(&frobenius_constant[f_p10][1],&tmp1);
	Fp2_set(&frobenius_constant[f_p10][2],&tmp2);
	Fp2_set(&frobenius_constant[f_p10][3],&tmp3);
	Fp2_mul(&frobenius_constant[f_p10][4],&tmp1,&tmp3);
	Fp2_mul(&frobenius_constant[f_p10][5],&tmp2,&tmp3);
	//set skew_f_10
	Fp2_inv(&tmp1,&tmp1);
	mpz_sub_ui(exp,p10,1);
	mpz_tdiv_q_ui(exp,exp,2);
	Fp2_pow(&tmp2,&Alpha_1,exp);
	Fp2_inv(&tmp2,&tmp2);
	Fp2_set(&skew_frobenius_constant[f_p10][0],&tmp1);
	Fp2_set(&skew_frobenius_constant[f_p10][1],&tmp2);
	
	
	//skew_frobenius_1
	mpz_sub_ui(exp,prime,1);
	mpz_mul_ui(exp,exp,2);
	Fp2_pow(&skew_frobenius_constant[f_p1][0],&Alpha_1,exp);
	Fp2_inv(&skew_frobenius_constant[f_p1][0],&skew_frobenius_constant[f_p1][0]);
	Fp2_mul(&skew_frobenius_constant[f_p1][0],&skew_frobenius_constant[f_p1][0],&frobenius_constant[f_p1][1]);
	mpz_sub_ui(exp,prime,1);
	mpz_mul_ui(exp,exp,3);
	Fp2_pow(&skew_frobenius_constant[f_p1][1],&Alpha_1,exp);
	Fp2_inv(&skew_frobenius_constant[f_p1][1],&skew_frobenius_constant[f_p1][1]);
	Fp2_mul(&skew_frobenius_constant[f_p1][1],&skew_frobenius_constant[f_p1][1],&frobenius_constant[f_p1][4]);
	
	//skew_frobenius_2
	mpz_sub_ui(exp,p2,1);
	mpz_mul_ui(exp,exp,2);
	Fp2_pow(&skew_frobenius_constant[f_p2][0],&Alpha_1,exp);
	Fp2_inv(&skew_frobenius_constant[f_p2][0],&skew_frobenius_constant[f_p2][0]);
	Fp2_mul(&skew_frobenius_constant[f_p2][0],&skew_frobenius_constant[f_p2][0],&frobenius_constant[f_p2][1]);
	mpz_sub_ui(exp,p2,1);
	mpz_mul_ui(exp,exp,3);
	Fp2_pow(&skew_frobenius_constant[f_p2][1],&Alpha_1,exp);
	Fp2_inv(&skew_frobenius_constant[f_p2][1],&skew_frobenius_constant[f_p2][1]);
	Fp2_mul(&skew_frobenius_constant[f_p2][1],&skew_frobenius_constant[f_p2][1],&frobenius_constant[f_p2][4]);
	
	//skew_frobenius_3
	mpz_sub_ui(exp,p3,1);
	mpz_mul_ui(exp,exp,2);
	Fp2_pow(&skew_frobenius_constant[f_p3][0],&Alpha_1,exp);
	Fp2_inv(&skew_frobenius_constant[f_p3][0],&skew_frobenius_constant[f_p3][0]);
	Fp2_mul(&skew_frobenius_constant[f_p3][0],&skew_frobenius_constant[f_p3][0],&frobenius_constant[f_p3][1]);
	mpz_sub_ui(exp,p3,1);
	mpz_mul_ui(exp,exp,3);
	Fp2_pow(&skew_frobenius_constant[f_p3][1],&Alpha_1,exp);
	Fp2_inv(&skew_frobenius_constant[f_p3][1],&skew_frobenius_constant[f_p3][1]);
	Fp2_mul(&skew_frobenius_constant[f_p3][1],&skew_frobenius_constant[f_p3][1],&frobenius_constant[f_p3][4]);
	
	//skew_frobenius_10
	mpz_sub_ui(exp,p10,1);
	mpz_mul_ui(exp,exp,2);
	Fp2_pow(&skew_frobenius_constant[f_p10][0],&Alpha_1,exp);
	Fp2_inv(&skew_frobenius_constant[f_p10][0],&skew_frobenius_constant[f_p10][0]);
	Fp2_mul(&skew_frobenius_constant[f_p10][0],&skew_frobenius_constant[f_p10][0],&frobenius_constant[f_p10][1]);
	mpz_sub_ui(exp,p10,1);
	mpz_mul_ui(exp,exp,3);
	Fp2_pow(&skew_frobenius_constant[f_p10][1],&Alpha_1,exp);
	Fp2_inv(&skew_frobenius_constant[f_p10][1],&skew_frobenius_constant[f_p10][1]);
	Fp2_mul(&skew_frobenius_constant[f_p10][1],&skew_frobenius_constant[f_p10][1],&frobenius_constant[f_p10][4]);
	
	Fp2_clear(&tmp1);
	Fp2_clear(&tmp2);
	Fp2_clear(&tmp3);
	mpz_clear(exp);
	mpz_clear(buf);
	mpz_clear(p2);
	mpz_clear(p3);
	mpz_clear(p4);
	mpz_clear(p6);
	mpz_clear(p8);
	mpz_clear(p10);
}

void set_curve_parameter(){
    mpz_set_ui(curve_b,2);
}

void BN12_print_parameters(){
    mpz_t mod_12;
    mpz_init(mod_12);
    printf("====================================================================================\n");
    printf("BN12 Class 2 (HW=4)\n\n");
    
    gmp_printf("Parameters\n");
    
    gmp_printf("X = 2^{114}+2^{84}-2^{53}-1\n");
    
    mpz_mod_ui(mod_12,X,12);
    gmp_printf("X mod 12 = %Zd\n\n",mod_12);
    
    gmp_printf("X     (%dbit length) : %Zd \n",(int)mpz_sizeinbase(X,2),X);
    gmp_printf("prime (%dbit length) : %Zd \n",(int)mpz_sizeinbase(prime,2),prime);
    gmp_printf("order (%dbit length) : %Zd \n",(int)mpz_sizeinbase(order,2),order);
    gmp_printf("trace (%dbit length) : %Zd \n\n",(int)mpz_sizeinbase(trace,2),trace);
    
    gmp_printf("BN curve\n");
    gmp_printf("E:y^2=x^3+2\n\n");
    
    gmp_printf("Twisted curve\n");
    gmp_printf("E':y^2=x^3+2(alpha+1)^{-1}\n\n");
    
    gmp_printf("Extension field\n");
    gmp_printf("Fp2 = Fp[alpha]/(alpha^2+1)\n");
    gmp_printf("Fp6 = Fp2[beta]/(beta^3-(alpha+1))\n");
    gmp_printf("Fp12= Fp6[gamma]/(gamma^2-beta)\n");
    
    mpz_clear(mod_12);
}

void BN12_clear(){
    int i,j;
    
    mpz_clear(X);
    mpz_clear(prime);
    mpz_clear(order);
    mpz_clear(trace);
    
    mpz_clear(EFp_total);
    mpz_clear(EFp12_total);
    mpz_clear(curve_b);
    
    mpz_clear(epsilon1);
    mpz_clear(epsilon2);
    mpz_clear(Two_inv);
    Fp2_clear(&Alpha_1);
    Fp2_clear(&Alpha_1_inv);
    
    for(i=0; i<12; i++){
        for(j=0; j<6; j++){
            Fp2_clear(&frobenius_constant[i][j]);
        }
        for(j=0; j<2; j++){
            Fp2_clear(&skew_frobenius_constant[i][j]);
        }
    }
}

/*============================================================================*/
/* Test functions                                                             */
/*============================================================================*/
/*----------------------------------------------------------------------------*/
//time
float timedifference_msec(struct timeval tv_start, struct timeval tv_end){
    return (tv_end.tv_sec - tv_start.tv_sec) * 1000.0f + (tv_end.tv_usec - tv_start.tv_usec) / 1000.0f;
}

float timedifference_usec(struct timeval tv_start, struct timeval tv_end){
    return (tv_end.tv_sec - tv_start.tv_sec) + (tv_end.tv_usec - tv_start.tv_usec);
}

/*----------------------------------------------------------------------------*/
//cost
void Init_mpz_Cost(struct mpz_Cost *cost){
    cost->mpz_mul=0;
    cost->mpz_mul_ui=0;
    cost->mpz_sqr=0;
    cost->mpz_add=0;
    cost->mpz_add_ui=0;
    cost->mpz_invert=0;
}

void Print_mpz_Cost(struct mpz_Cost *cost,char *str){
    printf("%s",str);
    printf("mpz_mul,mpz_mul_ui,mpz_sqr,mpz_add,mpz_add_ui,mpz_invert\n");
    printf("%ld,",cost->mpz_mul);
    printf("%ld,",cost->mpz_mul_ui);
    printf("%ld,",cost->mpz_sqr);
    printf("%ld,",cost->mpz_add);
    printf("%ld,",cost->mpz_add_ui);
    printf("%ld",cost->mpz_invert);
    printf("\n");
}

void Init_Fp_Cost(struct Fp_Cost *cost){
    cost->Fp_mul=0;
    cost->Fp_mul_mpz=0;
    cost->Fp_mul_ui=0;
    cost->Fp_sqr=0;
    cost->Fp_basis=0;
    cost->Fp_add=0;
    cost->Fp_add_mpz=0;
    cost->Fp_add_ui=0;
    cost->Fp_inv=0;
    cost->Fp_neg=0;
}

void Print_Fp_Cost(struct Fp_Cost *cost,char *str){
    printf("%s",str);
    printf("Fp_mul,Fp_mul_mpz,Fp_mul_ui,Fp_sqr,Fp_basis,Fp_add,Fp_add_mpz,Fp_add_ui,Fp_inv,Fp_neg\n");
    printf("%ld,",cost->Fp_mul);
    printf("%ld,",cost->Fp_mul_mpz);
    printf("%ld,",cost->Fp_mul_ui);
    printf("%ld,",cost->Fp_sqr);
    printf("%ld,",cost->Fp_basis);
    printf("%ld,",cost->Fp_add);
    printf("%ld,",cost->Fp_add_mpz);
    printf("%ld,",cost->Fp_add_ui);
    printf("%ld,",cost->Fp_inv);
    printf("%ld",cost->Fp_neg);
    printf("\n");
}

/*----------------------------------------------------------------------------*/
//test
void test_Field(){
    printf("====================================================================================\n");
    Fp6 tmp_Fp6,test1,test2;
    Fp6_init(&tmp_Fp6);
    Fp6_init(&test1);
    Fp6_init(&test2);
    mpz_t exp;
    mpz_init(exp);
    
    gmp_randinit_default (state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    Fp6_set_random(&tmp_Fp6,state);
    Fp6_printf(&tmp_Fp6,"");
    printf("\n\n");
    
    printf("mul/sqr\n");
    Fp6_mul(&test1,&tmp_Fp6,&tmp_Fp6);
    Fp6_printf(&test1,"");
    printf("\n");
    
    Fp6_sqr(&test2,&tmp_Fp6);
    Fp6_printf(&test2,"");
    printf("\n");
    
    if(Fp6_cmp(&test1,&test2)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    printf("pow/inv\n");
    mpz_pow_ui(exp,prime,12);
    mpz_sub_ui(exp,exp,2);
    Fp6_pow(&test1,&tmp_Fp6,exp);
    Fp6_printf(&test1,"");
    printf("\n");
    
    Fp6_inv(&test2,&tmp_Fp6);
    Fp6_printf(&test2,"");
    printf("\n");
    
    if(Fp6_cmp(&test1,&test2)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(exp);
    Fp6_clear(&tmp_Fp6);
    Fp6_clear(&test1);
    Fp6_clear(&test2);
}

void test_Frobenius_map(){
    printf("====================================================================================\n");
     Fp12 tmp_Fp12,test1,test2;
    Fp12_init(&tmp_Fp12);
    Fp12_init(&test1);
    Fp12_init(&test2);
    mpz_t exp;
    mpz_init(exp);
    
    gmp_randinit_default (state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    Fp12_set_random(&tmp_Fp12,state);
    Fp12_printf(&tmp_Fp12,"");
    printf("\n\n");
    
    printf("frobenius\n");
    Fp12_frobenius_map_p10(&test1,&tmp_Fp12);
    Fp12_printf(&test1,"");
    printf("\n\n");
    
    mpz_pow_ui(exp,prime,10);
    Fp12_pow(&test2,&tmp_Fp12,exp);
    Fp12_printf(&test2,"");
    printf("\n");
    
    if(Fp12_cmp(&test1,&test2)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(exp);
    Fp12_clear(&tmp_Fp12);
    Fp12_clear(&test1);
    Fp12_clear(&test2);
}

void test_skew_frobenius_map(){
    printf("====================================================================================\n");
    EFp12 Q,test1,test2;
    EFp12_init(&Q);
    EFp12_init(&test1);
    EFp12_init(&test2);
    EFp2 twisted_Q;
    EFp2_init(&twisted_Q);
    
    EFp12_generate_G2(&Q);
    EFp12_to_EFp2(&twisted_Q,&Q);
    
    Fp12_frobenius_map_p10(&test1.x,&Q.x);
    Fp12_frobenius_map_p10(&test1.y,&Q.y);
    EFp12_printf(&test1,"");
    printf("\n\n");
    
    EFp2_skew_frobenius_map_p10(&twisted_Q,&twisted_Q);
    EFp2_to_EFp12(&test2,&twisted_Q);
    EFp12_printf(&test2,"");
    printf("\n\n");
    
    if(Fp12_cmp(&test1.x,&test2.x)==0 && Fp12_cmp(&test1.y,&test2.y)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    EFp12_clear(&Q);
    EFp12_clear(&test1);
    EFp12_clear(&test2);
    EFp2_clear(&twisted_Q);
}

void test_rational_point(){
    printf("====================================================================================\n");
    EFp12 test_G1,test_G2;
    EFp12_init(&test_G1);
    EFp12_init(&test_G2);
    
    EFp12_generate_G1(&test_G1);
    EFp12_printf(&test_G1,"G1\n");
    printf("\n");
    EFp12_SCM(&test_G1,&test_G1,order);
    EFp12_printf(&test_G1,"G1 test\n");
    printf("\n");
    
    EFp12_generate_G2(&test_G2);
    EFp12_printf(&test_G2,"G2\n");
    printf("\n");
    EFp12_SCM(&test_G2,&test_G2,order);
    EFp12_printf(&test_G2,"G2 test\n");
    printf("\n");
    
    EFp12_clear(&test_G1);
    EFp12_clear(&test_G2);
}

void test_twist(){
    printf("====================================================================================\n");
    EFp12 Q,test1,test2;
    EFp12_init(&Q);
    EFp12_init(&test1);
    EFp12_init(&test2);
    EFp2 twist_Q;
    EFp2_init(&twist_Q);
    
    
    EFp12_generate_G2(&Q);
    EFp12_printf(&Q,"Q\n");
    printf("\n\n");
    
    EFp12_to_EFp2(&twist_Q,&Q);
    EFp2_ECD(&twist_Q,&twist_Q);
    EFp2_to_EFp12(&test1,&twist_Q);
    EFp12_printf(&test1,"");
    printf("\n\n");
    
    EFp12_ECD(&test2,&Q);
    EFp12_printf(&test2,"");
    printf("\n\n");
    
    if(Fp12_cmp(&test1.x,&test2.x)==0 && Fp12_cmp(&test1.y,&test2.y)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    EFp12_clear(&Q);
    EFp12_clear(&test1);
    EFp12_clear(&test2);
    EFp2_clear(&twist_Q);
}

void test_prototype_pairing(){
    printf("====================================================================================\n");
    EFp12 P,Q,s1P,s2Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&s1P);
    EFp12_init(&s2Q);
    Fp12 Z,test1,test2;
    Fp12_init(&Z);
    Fp12_init(&test1);
    Fp12_init(&test2);
    mpz_t s12,s1,s2;
    mpz_init(s12);
    mpz_init(s1);
    mpz_init(s2);
    
    mpz_set_str(s1,"723875134982754847463759",10);
    mpz_set_str(s2,"127468259486348759437251",10);
    mpz_mul(s12,s1,s2);
    
    EFp12_generate_G1(&P);
    EFp12_generate_G2(&Q);
    EFp12_SCM(&s1P,&P,s1);
    EFp12_SCM(&s2Q,&Q,s2);
    
    Prototype_pairing(&Z,&Q,&P);
    Fp12_pow(&test1,&Z,s12);
    Fp12_printf(&test1,"");
    printf("\n\n");
    
    Prototype_pairing(&test2,&s2Q,&s1P);
    Fp12_printf(&test2,"");
    printf("\n\n");
    
    if(Fp12_cmp(&test1,&test2)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    
    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);
    EFp12_clear(&P);
    EFp12_clear(&Q);
    EFp12_clear(&s1P);
    EFp12_clear(&s2Q);
    Fp12_clear(&Z);
    Fp12_clear(&test1);
    Fp12_clear(&test2);
}

void test_plain_ate_pairing(){
    printf("====================================================================================\n");
    printf("Plain-ate pairing\n\n");
    EFp12 P,Q,s1P,s2P,s1Q,s2Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&s1P);
    EFp12_init(&s2P);
    EFp12_init(&s1Q);
    EFp12_init(&s2Q);
    Fp12 Z,test1,test2,test3;
    Fp12_init(&Z);
    Fp12_init(&test1);
    Fp12_init(&test2);
    Fp12_init(&test3);
    mpz_t s12,s1,s2;
    mpz_init(s12);
    mpz_init(s1);
    mpz_init(s2);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(s1,state,order);
	mpz_urandomm(s2,state,order);
    mpz_mul(s12,s1,s2);
    mpz_mod(s12,s12,order);
    
    EFp12_generate_G1(&P);
    EFp12_printf(&P,"P\n");
    printf("\n\n");
    EFp12_generate_G2(&Q);
    EFp12_printf(&Q,"Q\n");
    printf("\n\n");
    EFp12_SCM(&s1P,&P,s1);
    EFp12_SCM(&s2P,&P,s2);
    EFp12_SCM(&s1Q,&Q,s1);
    EFp12_SCM(&s2Q,&Q,s2);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain_ate(Q,P)^s1*s2\n");
    Plain_ate_pairing(&Z,&Q,&P);
    Fp12_pow(&test1,&Z,s12);
    printf("Miller's Algo. (Plain-ate) : %.2f[ms]\n",MILLER_PLAINATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("plain_ate([s2]Q,[s1]P)\n");
    Plain_ate_pairing(&test2,&s2Q,&s1P);
    printf("Miller's Algo. (Plain-ate) : %.2f[ms]\n",MILLER_PLAINATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("plain_ate([s1]Q,[s2]P)\n");
    Plain_ate_pairing(&test3,&s1Q,&s2P);
    printf("Miller's Algo. (Plain-ate) : %.2f[ms]\n",MILLER_PLAINATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test3,"");
    printf("\n\n");
    
    printf("bilinear test\n");
    if(Fp12_cmp(&test1,&test2)==0 && Fp12_cmp(&test1,&test3)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);
    EFp12_clear(&P);
    EFp12_clear(&Q);
    EFp12_clear(&s1P);
    EFp12_clear(&s2P);
    EFp12_clear(&s1Q);
    EFp12_clear(&s2Q);
    Fp12_clear(&Z);
    Fp12_clear(&test1);
    Fp12_clear(&test2);
    Fp12_clear(&test3);
}

void test_opt_ate_pairing(){
    printf("====================================================================================\n");
    printf("Opt-ate pairing\n\n");
    EFp12 P,Q,s1P,s2P,s1Q,s2Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&s1P);
    EFp12_init(&s2P);
    EFp12_init(&s1Q);
    EFp12_init(&s2Q);
    Fp12 Z,test1,test2,test3;
    Fp12_init(&Z);
    Fp12_init(&test1);
    Fp12_init(&test2);
    Fp12_init(&test3);
    mpz_t s12,s1,s2;
    mpz_init(s12);
    mpz_init(s1);
    mpz_init(s2);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(s1,state,order);
	mpz_urandomm(s2,state,order);
    mpz_mul(s12,s1,s2);
    mpz_mod(s12,s12,order);
    
    EFp12_generate_G1(&P);
    EFp12_printf(&P,"P\n");
    printf("\n\n");
    EFp12_generate_G2(&Q);
    EFp12_printf(&Q,"Q\n");
    printf("\n\n");
    EFp12_SCM(&s1P,&P,s1);
    EFp12_SCM(&s2P,&P,s2);
    EFp12_SCM(&s1Q,&Q,s1);
    EFp12_SCM(&s2Q,&Q,s2);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("opt_ate(Q,P)^s1*s2\n");
    Opt_ate_pairing(&Z,&Q,&P);
    Fp12_pow(&test1,&Z,s12);
    printf("Miller's Algo. (Opt-ate) : %.2f[ms]\n",MILLER_OPTATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("opt_ate([s2]Q,[s1]P)\n");
    Opt_ate_pairing(&test2,&s2Q,&s1P);
    printf("Miller's Algo. (Opt-ate) : %.2f[ms]\n",MILLER_OPTATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("opt_ate([s1]Q,[s2]P)\n");
    Opt_ate_pairing(&test3,&s1Q,&s2P);
    printf("Miller's Algo. (Opt-ate) : %.2f[ms]\n",MILLER_OPTATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test3,"");
    printf("\n\n");
    
    printf("bilinear test\n");
    if(Fp12_cmp(&test1,&test2)==0 && Fp12_cmp(&test1,&test3)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);
    EFp12_clear(&P);
    EFp12_clear(&Q);
    EFp12_clear(&s1P);
    EFp12_clear(&s2P);
    EFp12_clear(&s1Q);
    EFp12_clear(&s2Q);
    Fp12_clear(&Z);
    Fp12_clear(&test1);
    Fp12_clear(&test2);
    Fp12_clear(&test3);
}

void test_x_ate_pairing(){
    printf("====================================================================================\n");
    printf("X-ate pairing\n\n");
    EFp12 P,Q,s1P,s2P,s1Q,s2Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&s1P);
    EFp12_init(&s2P);
    EFp12_init(&s1Q);
    EFp12_init(&s2Q);
    Fp12 Z,test1,test2,test3;
    Fp12_init(&Z);
    Fp12_init(&test1);
    Fp12_init(&test2);
    Fp12_init(&test3);
    mpz_t s12,s1,s2;
    mpz_init(s12);
    mpz_init(s1);
    mpz_init(s2);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(s1,state,order);
	mpz_urandomm(s2,state,order);
    mpz_mul(s12,s1,s2);
    mpz_mod(s12,s12,order);
    
    EFp12_generate_G1(&P);
    EFp12_printf(&P,"P\n");
    printf("\n\n");
    EFp12_generate_G2(&Q);
    EFp12_printf(&Q,"Q\n");
    printf("\n\n");
    EFp12_SCM(&s1P,&P,s1);
    EFp12_SCM(&s2P,&P,s2);
    EFp12_SCM(&s1Q,&Q,s1);
    EFp12_SCM(&s2Q,&Q,s2);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("x_ate(Q,P)^s1*s2\n");
    X_ate_pairing(&Z,&Q,&P);
    Fp12_pow(&test1,&Z,s12);
    printf("Miller's Algo. (X-ate) : %.2f[ms]\n",MILLER_XATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("x_ate([s2]Q,[s1]P)\n");
    X_ate_pairing(&test2,&s2Q,&s1P);
    printf("Miller's Algo. (X-ate) : %.2f[ms]\n",MILLER_XATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("x_ate([s1]Q,[s2]P)\n");
    X_ate_pairing(&test3,&s1Q,&s2P);
    printf("Miller's Algo. (X-ate) : %.2f[ms]\n",MILLER_XATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test3,"");
    printf("\n\n");
    
    printf("bilinear test\n");
    if(Fp12_cmp(&test1,&test2)==0 && Fp12_cmp(&test1,&test3)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);
    EFp12_clear(&P);
    EFp12_clear(&Q);
    EFp12_clear(&s1P);
    EFp12_clear(&s2P);
    EFp12_clear(&s1Q);
    EFp12_clear(&s2Q);
    Fp12_clear(&Z);
    Fp12_clear(&test1);
    Fp12_clear(&test2);
    Fp12_clear(&test3);
}

void test_G1_SCM(){
    printf("====================================================================================\n");
    printf("G1 SCM\n\n");
    EFp12 P,test1,test2,test3;
    EFp12_init(&P);
    EFp12_init(&test1);
    EFp12_init(&test2);
    EFp12_init(&test3);
    mpz_t scalar;
    mpz_init(scalar);
    
    //gmp_randinit_default (state);
    //gmp_randseed_ui(state,(unsigned long)time(NULL));
    //mpz_urandomm(scalar,state,order);
    
    mpz_set_str(scalar,"6433987872172996767370742739789192855206084819751810317146156999801714004530205817441825065255166612237926117694268872766909697908042480267",10);
    printf("scalar : %d-bit\n",(int)mpz_sizeinbase(scalar,2));
    
    printf("generating rational point P in G1\n\n");
    EFp12_generate_G1(&P);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain\n");
    EFp12_G1_SCM_plain(&test1,&P,scalar);
    printf("G1 SCM (plain) : %.2f[ms]\n",G1SCM_PLAIN);
    EFp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("2split\n");
    EFp12_G1_SCM_2split(&test2,&P,scalar);
    printf("G1 SCM (2split) : %.2f[ms]\n",G1SCM_2SPLIT);
    EFp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("2split-JSF\n");
    EFp12_G1_SCM_2split_JSF(&test3,&P,scalar);
    printf("G1 SCM (2split-JSF) : %.2f[ms]\n",G1SCM_2SPLIT_JSF);
    EFp12_printf(&test3,"");
    printf("\n\n");
    
    if(Fp12_cmp(&test1.x,&test2.x)==0 && Fp12_cmp(&test1.y,&test2.y)==0
    && Fp12_cmp(&test1.x,&test3.x)==0 && Fp12_cmp(&test1.y,&test3.y)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(scalar);
    EFp12_clear(&P);
    EFp12_clear(&test1);
    EFp12_clear(&test2);
    EFp12_clear(&test3);
}

void test_G2_SCM(){
    printf("====================================================================================\n");
    printf("G2 SCM\n\n");
    EFp12 Q,test1,test2,test3,test4;
    EFp12_init(&Q);
    EFp12_init(&test1);
    EFp12_init(&test2);
    EFp12_init(&test3);
    EFp12_init(&test4);
    mpz_t scalar;
    mpz_init(scalar);
    
    //gmp_randinit_default (state);
    //gmp_randseed_ui(state,(unsigned long)time(NULL));
    //mpz_urandomm(scalar,state,order);
    
    mpz_set_str(scalar,"6433987872172996767370742739789192855206084819751810317146156999801714004530205817441825065255166612237926117694268872766909697908042480267",10);
    printf("scalar : %d-bit\n",(int)mpz_sizeinbase(scalar,2));
    
    printf("generating rational point Q in G2\n\n");
    EFp12_generate_G2(&Q);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain\n");
    EFp12_G2_SCM_plain(&test1,&Q,scalar);
    printf("G2 SCM (plain) : %.2f[ms]\n",G2SCM_PLAIN);
    EFp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("2split\n");
    EFp12_G2_SCM_2split(&test2,&Q,scalar);
    printf("G2 SCM (2split) : %.2f[ms]\n",G2SCM_2SPLIT);
    EFp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("2split-JSF\n");
    EFp12_G2_SCM_2split_JSF(&test3,&Q,scalar);
    printf("G2 SCM (2split-JSF) : %.2f[ms]\n",G2SCM_2SPLIT_JSF);
    EFp12_printf(&test3,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test4\n\n");
    printf("4split\n");
    EFp12_G2_SCM_4split(&test4,&Q,scalar);
    printf("G2 SCM (4split) : %.2f[ms]\n",G2SCM_4SPLIT);
    EFp12_printf(&test4,"");
    printf("\n\n");
    
    if(Fp12_cmp(&test1.x,&test2.x)==0 && Fp12_cmp(&test1.y,&test2.y)==0
    && Fp12_cmp(&test1.x,&test3.x)==0 && Fp12_cmp(&test1.y,&test3.y)==0
    && Fp12_cmp(&test1.x,&test4.x)==0 && Fp12_cmp(&test1.y,&test4.y)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(scalar);
    EFp12_clear(&Q);
    EFp12_clear(&test1);
    EFp12_clear(&test2);
    EFp12_clear(&test3);
    EFp12_clear(&test4);
}

void test_G3_EXP(){
    printf("====================================================================================\n");
    printf("G3 Exp.\n\n");
    EFp12 P,Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    Fp12 Z,test1,test2,test3,test4;
    Fp12_init(&Z);
    Fp12_init(&test1);
    Fp12_init(&test2);
    Fp12_init(&test3);
    Fp12_init(&test4);
    mpz_t scalar;
    mpz_init(scalar);
    
    //gmp_randinit_default (state);
    //gmp_randseed_ui(state,(unsigned long)time(NULL));
    //mpz_urandomm(scalar,state,order);
    
    mpz_set_str(scalar,"6433987872172996767370742739789192855206084819751810317146156999801714004530205817441825065255166612237926117694268872766909697908042480267",10);
    printf("scalar : %d-bit\n",(int)mpz_sizeinbase(scalar,2));
    
    printf("generating rational point P in G1\n\n");
    EFp12_generate_G1(&P);
    printf("generating rational point Q in G2\n\n");
    EFp12_generate_G2(&Q);
    printf("x-ate(Q,P)\n");
    Opt_ate_pairing(&Z,&Q,&P);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain\n");
    Fp12_G3_EXP_plain(&test1,&Z,scalar);
    printf("G3 SCM (plain) : %.2f[ms]\n",G3EXP_PLAIN);
    Fp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("2split\n");
    Fp12_G3_EXP_2split(&test2,&Z,scalar);
    printf("G3 SCM (2split) : %.2f[ms]\n",G3EXP_2SPLIT);
    Fp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("2split-JSF\n");
    Fp12_G3_EXP_2split_JSF(&test3,&Z,scalar);
    printf("G3 SCM (2split-JSF) : %.2f[ms]\n",G3EXP_2SPLIT_JSF);
    Fp12_printf(&test3,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test4\n\n");
    printf("4split\n");
    Fp12_G3_EXP_4split(&test4,&Z,scalar);
    printf("G3 SCM (4split) : %.2f[ms]\n",G3EXP_4SPLIT);
    Fp12_printf(&test4,"");
    printf("\n\n");
    
    if(Fp12_cmp(&test1,&test2)==0 && Fp12_cmp(&test1,&test3)==0 && Fp12_cmp(&test1,&test4)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(scalar);
    EFp12_clear(&P);
    EFp12_clear(&Q);
    Fp12_clear(&Z);
    Fp12_clear(&test1);
    Fp12_clear(&test2);
    Fp12_clear(&test3);
    Fp12_clear(&test4);
}

void compare_pairings(){
    printf("====================================================================================\n");
    printf("Ate-based pairing\n\n");
    EFp12 P,Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    Fp12 Z;
    Fp12_init(&Z);
    EFp2 twisted_Q;
    EFp2_init(&twisted_Q);
    
    printf("generating rational point P in G1...\n\n");
    EFp12_generate_G1(&P);
    //EFp12_printf(&P,"P\n");
    //printf("\n\n");
    printf("generating rational point Q in G2...\n\n");
    EFp12_generate_G2(&Q);
    //EFp12_printf(&Q,"Q\n");
    //printf("\n\n");
    
    /*printf("------------------------------------------------------------------------------------\n");
    printf("Plain-ate pairing\n\n");
    Plain_ate_pairing(&Z,&Q,&P);
    printf("Miller's Algo. (Plain-ate) : %.2f[ms]\n",MILLER_PLAINATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Pairing (total) : %.2f[ms]\n",MILLER_PLAINATE+FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    printf("\n");*/
    
    printf("------------------------------------------------------------------------------------\n");
    printf("Opt-ate pairing\n\n");
    Opt_ate_pairing(&Z,&Q,&P);
    printf("Miller's Algo. (Opt-ate)   : %.2f[ms]\n",MILLER_OPTATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Pairing (total) : %.2f[ms]\n",MILLER_OPTATE+FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    printf("\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("X-ate pairing\n\n");
    X_ate_pairing(&Z,&Q,&P);
    printf("Miller's Algo. (X-ate)     : %.2f[ms]\n",MILLER_XATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Pairing (total) : %.2f[ms]\n",MILLER_XATE+FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    printf("\n");
    
    /*printf("------------------------------------------------------------------------------------\n");
    printf("Sextic twist\n\n");
    gettimeofday(&tv_start,NULL);
    EFp12_to_EFp2(&twisted_Q,&Q);
    gettimeofday(&tv_end,NULL);
    printf("EFp12 to EFp2     : %.2f[us]\n",timedifference_usec(tv_start,tv_end));
    
    gettimeofday(&tv_start,NULL);
    EFp2_to_EFp12(&Q,&twisted_Q);
    gettimeofday(&tv_end,NULL);
    printf("EFp2 to EFp12     : %.2f[us]\n",timedifference_usec(tv_start,tv_end));
    */
    
    EFp12_clear(&P);
    EFp12_clear(&Q);
    Fp12_clear(&Z);
}

void operation_cost(){
    printf("====================================================================================\n");
    printf("Operation cost\n\n");
    Fp2 tmp1_Fp2,tmp2_Fp2;
    Fp2_init(&tmp1_Fp2);
    Fp2_init(&tmp2_Fp2);
    Fp6 tmp1_Fp6,tmp2_Fp6;
    Fp6_init(&tmp1_Fp6);
    Fp6_init(&tmp2_Fp6);
    Fp12 tmp1_Fp12,tmp2_Fp12;
    Fp12_init(&tmp1_Fp12);
    Fp12_init(&tmp2_Fp12);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    Fp2_set_random(&tmp1_Fp2,state);
    Fp6_set_random(&tmp1_Fp6,state);
    Fp12_set_random(&tmp1_Fp12,state);
    Fp2_set_random(&tmp2_Fp2,state);
    Fp6_set_random(&tmp2_Fp6,state);
    Fp12_set_random(&tmp2_Fp12,state);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("Extension field Fp2\n\n");
    Init_Fp_Cost(&Fp_cost);
    Fp2_mul(&tmp1_Fp2,&tmp1_Fp2,&tmp2_Fp2);
    Print_Fp_Cost(&Fp_cost,"Fp2_mul\n");
    
    Init_Fp_Cost(&Fp_cost);
    Fp2_sqr(&tmp1_Fp2,&tmp2_Fp2);
    Print_Fp_Cost(&Fp_cost,"Fp2_sqr\n");
    
    Init_Fp_Cost(&Fp_cost);
    Fp2_inv(&tmp1_Fp2,&tmp2_Fp2);
    Print_Fp_Cost(&Fp_cost,"Fp2_inv\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("Extension field Fp6\n\n");
    Init_Fp_Cost(&Fp_cost);
    Fp6_mul(&tmp1_Fp6,&tmp1_Fp6,&tmp2_Fp6);
    Print_Fp_Cost(&Fp_cost,"Fp6_mul\n");
    
    Init_Fp_Cost(&Fp_cost);
    Fp6_sqr(&tmp1_Fp6,&tmp2_Fp6);
    Print_Fp_Cost(&Fp_cost,"Fp6_sqr\n");
    
    Init_Fp_Cost(&Fp_cost);
    Fp6_inv(&tmp1_Fp6,&tmp2_Fp6);
    Print_Fp_Cost(&Fp_cost,"Fp6_inv\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("Extension field Fp12\n\n");
    Init_Fp_Cost(&Fp_cost);
    Fp12_mul(&tmp1_Fp12,&tmp1_Fp12,&tmp2_Fp12);
    Print_Fp_Cost(&Fp_cost,"Fp12_mul\n");
    
    Init_Fp_Cost(&Fp_cost);
    Fp12_sqr(&tmp1_Fp12,&tmp2_Fp12);
    Print_Fp_Cost(&Fp_cost,"Fp12_sqr\n");
    
    Init_Fp_Cost(&Fp_cost);
    Fp12_inv(&tmp1_Fp12,&tmp2_Fp12);
    Print_Fp_Cost(&Fp_cost,"Fp12_inv\n");
    
    Init_Fp_Cost(&Fp_cost);
    Fp12_frobenius_map_p1(&tmp1_Fp12,&tmp2_Fp12);
    Print_Fp_Cost(&Fp_cost,"Fp12_frobenius_map_p1\n");
    
    Init_Fp_Cost(&Fp_cost);
    Fp12_frobenius_map_p2(&tmp1_Fp12,&tmp2_Fp12);
    Print_Fp_Cost(&Fp_cost,"Fp12_frobenius_map_p2\n");
    
    Init_Fp_Cost(&Fp_cost);
    Fp12_frobenius_map_p3(&tmp1_Fp12,&tmp2_Fp12);
    Print_Fp_Cost(&Fp_cost,"Fp12_frobenius_map_p3\n");
    
    Init_Fp_Cost(&Fp_cost);
    Fp12_frobenius_map_p4(&tmp1_Fp12,&tmp2_Fp12);
    Print_Fp_Cost(&Fp_cost,"Fp12_frobenius_map_p4\n");
    
    Init_Fp_Cost(&Fp_cost);
    Fp12_frobenius_map_p6(&tmp1_Fp12,&tmp2_Fp12);
    Print_Fp_Cost(&Fp_cost,"Fp12_frobenius_map_p6\n");
    
    Init_Fp_Cost(&Fp_cost);
    Fp12_frobenius_map_p8(&tmp1_Fp12,&tmp2_Fp12);
    Print_Fp_Cost(&Fp_cost,"Fp12_frobenius_map_p8\n");
    
    Init_Fp_Cost(&Fp_cost);
    Fp12_frobenius_map_p10(&tmp1_Fp12,&tmp2_Fp12);
    Print_Fp_Cost(&Fp_cost,"Fp12_frobenius_map_p10\n");
    
    
    Fp2_clear(&tmp1_Fp2);
    Fp6_clear(&tmp1_Fp6);
    Fp12_clear(&tmp1_Fp12);
    Fp2_clear(&tmp2_Fp2);
    Fp6_clear(&tmp2_Fp6);
    Fp12_clear(&tmp2_Fp12);
}

