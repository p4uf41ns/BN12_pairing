#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#define X_length 114
#define X6_2_length 116

typedef struct{
	mpz_t x0;
}Fp;

typedef struct{
	Fp x0,x1;
}Fp2;

typedef struct{
    Fp2 x0,x1,x2;
}Fp6;

typedef struct{
    Fp6 x0,x1;
}Fp12;
/*============================================================================*/
/* Field                                                                      */
/*============================================================================*/
//Fp
void Fp_init(Fp *A);
void Fp_clear(Fp *A);
void Fp_printf(Fp *A,char *str);
void Fp_set(Fp *ANS,Fp *A);
void Fp_set_ui(Fp *ANS,unsigned long int UI);
void Fp_set_mpz(Fp *ANS,mpz_t A);
void Fp_set_neg(Fp *ANS,Fp *A);
void Fp_set_random(Fp *ANS,gmp_randstate_t state);
void Fp_mul(Fp *ANS,Fp *A,Fp *B);
void Fp_mul_ui(Fp *ANS,Fp *A,unsigned long int UI);
void Fp_mul_mpz(Fp *ANS,Fp *A,mpz_t B);
void Fp_add(Fp *ANS,Fp *A,Fp *B);
void Fp_add_ui(Fp *ANS,Fp *A,unsigned long int UI);
void Fp_add_mpz(Fp *ANS,Fp *A,mpz_t B);
void Fp_sub(Fp *ANS,Fp *A,Fp *B);
void Fp_sub_ui(Fp *ANS,Fp *A,unsigned long int UI);
void Fp_sub_mpz(Fp *ANS,Fp *A,mpz_t B);
void Fp_inv(Fp *ANS,Fp *A);
int  Fp_legendre(Fp *A);
int  Fp_isCNR(Fp *A);
void Fp_sqrt(Fp *ANS,Fp *A);
void Fp_pow(Fp *ANS,Fp *A,mpz_t scalar);
int  Fp_cmp(Fp *A,Fp *B);
int  Fp_cmp_ui(Fp *A,unsigned long int UI);
int  Fp_cmp_mpz(Fp *A,mpz_t B);
int  Fp_cmp_zero(Fp *A);
int  Fp_cmp_one(Fp *A);
/*----------------------------------------------------------------------------*/
//Fp2
void Fp2_init(Fp2 *A);
void Fp2_clear(Fp2 *A);
void Fp2_printf(Fp2 *A,char *str);
void Fp2_set(Fp2 *ANS,Fp2 *A);
void Fp2_set_ui(Fp2 *ANS,unsigned long int UI);
void Fp2_set_mpz(Fp2 *ANS,mpz_t A);
void Fp2_set_neg(Fp2 *ANS,Fp2 *A);
void Fp2_set_random(Fp2 *ANS,gmp_randstate_t state);
void Fp2_mul(Fp2 *ANS,Fp2 *A,Fp2 *B);
void Fp2_mul_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI);
void Fp2_mul_mpz(Fp2 *ANS,Fp2 *A,mpz_t B);
void Fp2_mul_basis(Fp2 *ANS,Fp2 *A);
void Fp2_inv_basis(Fp2 *ANS,Fp2 *A);
void Fp2_sqr(Fp2 *ANS,Fp2 *A);
void Fp2_add(Fp2 *ANS,Fp2 *A,Fp2 *B);
void Fp2_add_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI);
void Fp2_add_mpz(Fp2 *ANS,Fp2 *A,mpz_t B);
void Fp2_sub(Fp2 *ANS,Fp2 *A,Fp2 *B);
void Fp2_sub_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI);
void Fp2_sub_mpz(Fp2 *ANS,Fp2 *A,mpz_t B);
void Fp2_inv(Fp2 *ANS,Fp2 *A);
int  Fp2_legendre(Fp2 *A);
int  Fp2_isCNR(Fp2 *A);
void Fp2_sqrt(Fp2 *ANS,Fp2 *A);
void Fp2_pow(Fp2 *ANS,Fp2 *A,mpz_t scalar);
int  Fp2_cmp(Fp2 *A,Fp2 *B);
int  Fp2_cmp_ui(Fp2 *A,unsigned long int UI);
int  Fp2_cmp_mpz(Fp2 *A,mpz_t B);
int  Fp2_cmp_zero(Fp2 *A);
int  Fp2_cmp_one(Fp2 *A);
/*----------------------------------------------------------------------------*/
//Fp6
void Fp6_init(Fp6 *A);
void Fp6_clear(Fp6 *A);
void Fp6_printf(Fp6 *A,char *str);
void Fp6_set(Fp6 *ANS,Fp6 *A);
void Fp6_set_ui(Fp6 *ANS,unsigned long int UI);
void Fp6_set_mpz(Fp6 *ANS,mpz_t A);
void Fp6_set_neg(Fp6 *ANS,Fp6 *A);
void Fp6_set_random(Fp6 *ANS,gmp_randstate_t state);
void Fp6_mul(Fp6 *ANS,Fp6 *A,Fp6 *B);
void Fp6_mul_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI);
void Fp6_mul_mpz(Fp6 *ANS,Fp6 *A,mpz_t B);
void Fp6_mul_basis(Fp6 *ANS,Fp6 *A);
void Fp6_sqr(Fp6 *ANS,Fp6 *A);
void Fp6_add(Fp6 *ANS,Fp6 *A,Fp6 *B);
void Fp6_add_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI);
void Fp6_add_mpz(Fp6 *ANS,Fp6 *A,mpz_t B);
void Fp6_sub(Fp6 *ANS,Fp6 *A,Fp6 *B);
void Fp6_sub_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI);
void Fp6_sub_mpz(Fp6 *ANS,Fp6 *A,mpz_t B);
void Fp6_inv(Fp6 *ANS,Fp6 *A);
int  Fp6_legendre(Fp6 *A);
int  Fp6_isCNR(Fp6 *A);
void Fp6_sqrt(Fp6 *ANS,Fp6 *A);
void Fp6_pow(Fp6 *ANS,Fp6 *A,mpz_t scalar);
int  Fp6_cmp(Fp6 *A,Fp6 *B);
int  Fp6_cmp_ui(Fp6 *A,unsigned long int UI);
int  Fp6_cmp_mpz(Fp6 *A,mpz_t B);
int  Fp6_cmp_zero(Fp6 *A);
int  Fp6_cmp_one(Fp6 *A);
/*----------------------------------------------------------------------------*/
//Fp12
void Fp12_init(Fp12 *A);
void Fp12_clear(Fp12 *A);
void Fp12_printf(Fp12 *A,char *str);
void Fp12_set(Fp12 *ANS,Fp12 *A);
void Fp12_set_ui(Fp12 *ANS,unsigned long int UI);
void Fp12_set_mpz(Fp12 *ANS,mpz_t A);
void Fp12_set_neg(Fp12 *ANS,Fp12 *A);
void Fp12_set_random(Fp12 *ANS,gmp_randstate_t state);
void Fp12_mul(Fp12 *ANS,Fp12 *A,Fp12 *B);
void Fp12_mul_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI);
void Fp12_mul_mpz(Fp12 *ANS,Fp12 *A,mpz_t B);
void Fp12_sqr(Fp12 *ANS,Fp12 *A);
void Fp12_sqr_cyclic(Fp12 *ANS,Fp12 *A);
void Fp12_add(Fp12 *ANS,Fp12 *A,Fp12 *B);
void Fp12_add_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI);
void Fp12_add_mpz(Fp12 *ANS,Fp12 *A,mpz_t B);
void Fp12_sub(Fp12 *ANS,Fp12 *A,Fp12 *B);
void Fp12_sub_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI);
void Fp12_sub_mpz(Fp12 *ANS,Fp12 *A,mpz_t B);
void Fp12_inv(Fp12 *ANS,Fp12 *A);
int  Fp12_legendre(Fp12 *A);
int  Fp12_isCNR(Fp12 *A);
void Fp12_sqrt(Fp12 *ANS,Fp12 *A);
void Fp12_pow(Fp12 *ANS,Fp12 *A,mpz_t scalar);
int  Fp12_cmp(Fp12 *A,Fp12 *B);
int  Fp12_cmp_ui(Fp12 *A,unsigned long int UI);
int  Fp12_cmp_mpz(Fp12 *A,mpz_t B);
int  Fp12_cmp_zero(Fp12 *A);
int  Fp12_cmp_one(Fp12 *A);
void Fp12_frobenius_map_p1(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p2(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p3(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p4(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p6(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p8(Fp12 *ANS,Fp12 *A);
void Fp12_frobenius_map_p10(Fp12 *ANS,Fp12 *A);

/*============================================================================*/
/* Elliptic Curve                                                             */
/*============================================================================*/
typedef struct{
	Fp x,y;
	int infinity;
}EFp;

typedef struct{
    Fp2 x,y;
	int infinity;
}EFp2;

typedef struct{
    Fp6 x,y;
	int infinity;
}EFp6;

typedef struct{
    Fp12 x,y;
	int infinity;
}EFp12;
/*----------------------------------------------------------------------------*/
//EFp
void EFp_init(EFp *P);
void EFp_set(EFp *P,EFp *A);
void EFp_set_ui(EFp *ANS,unsigned long int UI);
void EFp_set_mpz(EFp *ANS,mpz_t A);
void EFp_set_neg(EFp *ANS,EFp *A);
void EFp_clear(EFp *P);
void EFp_printf(EFp *P,char *str);
void EFp_rational_point(EFp *P);
void EFp_ECD(EFp *ANS,EFp *P);
void EFp_ECA(EFp *ANS,EFp *P1,EFp *P2);
void EFp_SCM(EFp *ANS,EFp *P,mpz_t scalar);
//skew_frobenius_map
void EFp_skew_frobenius_map_p2(EFp *ANS,EFp *A);
/*----------------------------------------------------------------------------*/
//EFp2
void EFp2_init(EFp2 *P);
void EFp2_set(EFp2 *ANS,EFp2 *A);
void EFp2_set_ui(EFp2 *ANS,unsigned long int UI);
void EFp2_set_mpz(EFp2 *ANS,mpz_t A);
void EFp2_set_neg(EFp2 *ANS,EFp2 *A);
void EFp2_clear(EFp2 *P);
void EFp2_printf(EFp2 *P,char *str);
void EFp2_rational_point(EFp2 *P);
void EFp2_ECD(EFp2 *ANS,EFp2 *P);
void EFp2_ECA(EFp2 *ANS,EFp2 *P1,EFp2 *P2);
void EFp2_SCM(EFp2 *ANS,EFp2 *P,mpz_t scalar);
//skew_frobenius_map
void EFp2_skew_frobenius_map_p1(EFp2 *ANS,EFp2 *A);
void EFp2_skew_frobenius_map_p2(EFp2 *ANS,EFp2 *A);
void EFp2_skew_frobenius_map_p3(EFp2 *ANS,EFp2 *A);
void EFp2_skew_frobenius_map_p10(EFp2 *ANS,EFp2 *A);
/*----------------------------------------------------------------------------*/
//EFp6
void EFp6_init(EFp6 *P);
void EFp6_set(EFp6 *ANS,EFp6 *A);
void EFp6_set_ui(EFp6 *ANS,unsigned long int UI);
void EFp6_set_mpz(EFp6 *ANS,mpz_t A);
void EFp6_set_neg(EFp6 *ANS,EFp6 *A);
void EFp6_clear(EFp6 *P);
void EFp6_printf(EFp6 *P,char *str);
void EFp6_rational_point(EFp6 *P);
void EFp6_ECD(EFp6 *ANS,EFp6 *P);
void EFp6_ECA(EFp6 *ANS,EFp6 *P1,EFp6 *P2);
void EFp6_SCM(EFp6 *ANS,EFp6 *P,mpz_t scalar);
/*----------------------------------------------------------------------------*/
//EFp12
void EFp12_init(EFp12 *P);
void EFp12_set(EFp12 *ANS,EFp12 *A);
void EFp12_set_ui(EFp12 *ANS,unsigned long int UI);
void EFp12_set_mpz(EFp12 *ANS,mpz_t A);
void EFp12_set_neg(EFp12 *ANS,EFp12 *A);
void EFp12_clear(EFp12 *P);
void EFp12_printf(EFp12 *P,char *str);
void EFp12_rational_point(EFp12 *P);
void EFp12_generate_G1(EFp12 *P);
void EFp12_generate_G2(EFp12 *Q);
void EFp12_ECD(EFp12 *ANS,EFp12 *P);
void EFp12_ECA(EFp12 *ANS,EFp12 *P1,EFp12 *P2);
void EFp12_SCM(EFp12 *ANS,EFp12 *P,mpz_t scalar);

/*============================================================================*/
/* Pairing functions                                                          */
/*============================================================================*/
enum f_state{
	f_p1,f_p2,f_p3,f_p4,f_p5,f_p6,f_p7,f_p8,f_p9,f_p10,f_p11,f_p12
};
gmp_randstate_t state;
int X_binary[X_length+1];
int X6_2_binary[X6_2_length+1];
mpz_t X,prime,order,trace;
mpz_t EFp_total,EFp12_total;
mpz_t curve_b;
Fp2 frobenius_constant[12][6];
Fp2 skew_frobenius_constant[12][2];
Fp2 Alpha_1,Alpha_1_inv;
mpz_t epsilon1,epsilon2;
mpz_t Two_inv;

/*----------------------------------------------------------------------------*/
//prototype
void EFp12_lTP(Fp12 *ANS,EFp12 *T,EFp12 *P,EFp12 *Q,Fp *L);
void EFp12_lTT(Fp12 *ANS,EFp12 *T,EFp12 *Q,Fp *L);
void EFp12_vTP(Fp12 *ANS,EFp12 *T,EFp12 *P,EFp12 *Q);
void EFp12_vTT(Fp12 *ANS,EFp12 *T,EFp12 *Q);
void Miller_prototype(Fp12 *ANS,EFp12 *P,EFp12 *Q);
void Prototype_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q);
/*----------------------------------------------------------------------------*/
//twist
void EFp12_to_EFp2(EFp2 *ANS,EFp12 *A);
void EFp2_to_EFp12(EFp12 *ANS,EFp2 *A);
void EFp12_to_EFp(EFp *ANS,EFp12 *A);
void EFp_to_EFp12(EFp12 *ANS,EFp *A);
/*----------------------------------------------------------------------------*/
//Pseudo 8-sparse
void Pseudo_8_sparse_mapping(EFp *P,EFp2 *Q,Fp *L);
void Pseudo_8_sparse_mul(Fp12 *ANS,Fp12 *A,Fp12 *B);
void ff_ltt(Fp12 *f,EFp2 *T,EFp *P,Fp *L);
void f_ltq(Fp12 *f,EFp2 *T,EFp2 *Q,EFp *P,Fp *L);
/*----------------------------------------------------------------------------*/
//miller
void Miller_algo_for_plain_ate(Fp12 *ANS,EFp12 *Q,EFp12 *P);
void Miller_algo_for_opt_ate(Fp12 *ANS,EFp12 *Q,EFp12 *P);
void Miller_algo_for_x_ate(Fp12 *ANS,EFp12 *Q,EFp12 *P);
/*----------------------------------------------------------------------------*/
//final exp
void Final_exp_plain(Fp12 *ANS,Fp12 *A);
void Fp12_pow_X(Fp12 *ANS,Fp12 *A);
void Final_exp_optimal(Fp12 *ANS,Fp12 *A);
/*----------------------------------------------------------------------------*/
//pairing
void Plain_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q);
void Opt_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q);
void X_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q);
/*----------------------------------------------------------------------------*/
//JSF
void Joint_sparse_form(int **binary,mpz_t S[2],int *loop_length);
//G1 SCM
void EFp12_G1_SCM_plain(EFp12 *ANS,EFp12 *P,mpz_t scalar);
void EFp12_G1_SCM_2split(EFp12 *ANS,EFp12 *P,mpz_t scalar);
void EFp12_G1_SCM_2split_JSF(EFp12 *ANS,EFp12 *P,mpz_t scalar);
/*----------------------------------------------------------------------------*/
//G2 SCM
void EFp12_G2_SCM_plain(EFp12 *ANS,EFp12 *Q,mpz_t scalar);
void EFp12_G2_SCM_2split(EFp12 *ANS,EFp12 *Q,mpz_t scalar);
void EFp12_G2_SCM_2split_JSF(EFp12 *ANS,EFp12 *Q,mpz_t scalar);
void EFp12_G2_SCM_4split(EFp12 *ANS,EFp12 *Q,mpz_t scalar);
/*----------------------------------------------------------------------------*/
//G3 EXP
void Fp12_G3_EXP_plain(Fp12 *ANS,Fp12 *A,mpz_t scalar);
void Fp12_G3_EXP_2split(Fp12 *ANS,Fp12 *A,mpz_t scalar);
void Fp12_G3_EXP_2split_JSF(Fp12 *ANS,Fp12 *A,mpz_t scalar);
void Fp12_G3_EXP_4split(Fp12 *ANS,Fp12 *A,mpz_t scalar);
/*----------------------------------------------------------------------------*/
//init/set/clear
void BN12_init();
void BN12_clear();
void BN12_print_parameters();
void init_parameters();
void generate_X();
int  generate_prime();
int  generate_order();
void generate_trace();
void weil();
void get_epsilon();
void get_Two_inv();
void set_basis();
void set_frobenius_constant();
void set_curve_parameter();

/*============================================================================*/
/* Test functions                                                             */
/*============================================================================*/
struct timeval tv_start,tv_end;
float MILLER_TATE,MILLER_PLAINATE,MILLER_OPTATE,MILLER_XATE;
float FINALEXP_PLAIN,FINALEXP_OPT_EASY,FINALEXP_OPT_HARD;
float G1SCM_PLAIN,G1SCM_2SPLIT,G1SCM_2SPLIT_JSF;
float G2SCM_PLAIN,G2SCM_2SPLIT,G2SCM_2SPLIT_JSF,G2SCM_4SPLIT;
float G3EXP_PLAIN,G3EXP_2SPLIT,G3EXP_2SPLIT_JSF,G3EXP_4SPLIT;
struct mpz_Cost{
    unsigned long int mpz_mul;
    unsigned long int mpz_mul_ui;
    unsigned long int mpz_sqr;
    unsigned long int mpz_add;
    unsigned long int mpz_add_ui;
    unsigned long int mpz_invert;
};
struct Fp_Cost{
    unsigned long int Fp_mul;
    unsigned long int Fp_mul_mpz;
    unsigned long int Fp_mul_ui;
    unsigned long int Fp_sqr;
    unsigned long int Fp_basis;
    unsigned long int Fp_add;
    unsigned long int Fp_add_mpz;
    unsigned long int Fp_add_ui;
    unsigned long int Fp_inv;
    unsigned long int Fp_neg;
};
struct mpz_Cost mpz_cost;
struct Fp_Cost Fp_cost;

/*----------------------------------------------------------------------------*/
//time
float timedifference_msec(struct timeval tv_start, struct timeval tv_end);
float timedifference_usec(struct timeval tv_start, struct timeval tv_end);
/*----------------------------------------------------------------------------*/
//cost
void Init_mpz_Cost(struct mpz_Cost *cost);
void Print_mpz_Cost(struct mpz_Cost *cost,char *str);
void Init_Fp_Cost(struct Fp_Cost *cost);
void Print_Fp_Cost(struct Fp_Cost *cost,char *str);
/*----------------------------------------------------------------------------*/
//test
void test_Field();
void test_Frobenius_map();
void test_skew_frobenius_map();
void test_rational_point();
void test_twist();
void test_prototype_pairing();
void test_plain_ate_pairing();
void test_opt_ate_pairing();
void test_x_ate_pairing();
void test_G1_SCM();
void test_G2_SCM();
void test_G3_EXP();
void compare_pairings();
void operation_cost();
