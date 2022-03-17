#ifndef ISC2019_ARB_VHSS_H
#define ISC2019_ARB_VHSS_H

#endif //ISC2019_ARB_VHSS_H

#include "PKE.h"
#include "GenData.h"

typedef struct {
    fmpz_poly_t alpha;
}VHSS_Para;

typedef struct {
    fmpz_poly_t s1;
    fmpz_poly_t s2;
}t_b;

void HSS_Gen(VHSS_Para* vhss_para, PKE_Para* pke_para, PKE_PK* pke_pk, PKE_SK* pke_sk,
             t_b* ek1_s,  t_b* ek1_as,t_b *ek2_s, t_b *ek2_as,int l_bit,int d)
{
    PKE_Gen(pke_para,pke_pk,pke_sk,l_bit,d);
    //ek1 s
    flint_randinit(pke_para->state);
    fmpz_poly_randtest_not_zero(ek1_s->s1,pke_para->state,pke_para->N,pke_para->q_bit);
    flint_randinit(pke_para->state);
    fmpz_poly_randtest_not_zero(ek1_s->s2,pke_para->state,pke_para->N,pke_para->q_bit);
    //ek2 s
    fmpz_poly_sub(ek2_s->s1,pke_sk->s1,ek1_s->s1);
    fmpz_poly_scalar_smod_fmpz(ek2_s->s1,ek2_s->s1,pke_para->q);
    fmpz_poly_sub(ek2_s->s2,pke_sk->s2,ek1_s->s2);
    fmpz_poly_scalar_smod_fmpz(ek2_s->s2,ek2_s->s2,pke_para->q);

    flint_randinit(pke_para->state);
    fmpz_poly_randtest_not_zero(vhss_para->alpha,pke_para->state,pke_para->N,1);
    fmpz_poly_t as;
    fmpz_poly_init(as);
    fmpz_poly_mul(as,vhss_para->alpha,pke_sk->s2);
    //ek1 alpha*s
    flint_randinit(pke_para->state);
    fmpz_poly_randtest_not_zero(ek1_as->s1,pke_para->state,pke_para->N,pke_para->q_bit);
    flint_randinit(pke_para->state);
    fmpz_poly_randtest_not_zero(ek1_as->s2,pke_para->state,pke_para->N,pke_para->q_bit);
    //ek2 alpha*s
    fmpz_poly_sub(ek2_as->s1,vhss_para->alpha,ek1_as->s1);
    fmpz_poly_scalar_smod_fmpz(ek2_as->s1,ek2_as->s1,pke_para->q);
    fmpz_poly_sub(ek2_as->s2,as,ek1_as->s2);
    fmpz_poly_scalar_smod_fmpz(ek2_as->s2,ek2_as->s2,pke_para->q);

}

void VHSS_Enc(fmpz_poly_t* x_c, fmpz_poly_t* xs_c, PKE_Para* pke_para,PKE_PK* pke_pk, fmpz_t x)
{
    PKE_OKDM(x_c, xs_c, pke_para,pke_pk, x);
}

void VHSS_Mult(t_b* tb, PKE_Para* pke_para, t_b* ek, fmpz_poly_t* x_c, fmpz_poly_t* xs_c)
{
    PKE_DDec(tb->s1, tb->s2, pke_para,ek->s1,ek->s2, x_c, xs_c);
}




void VHSS_Eval(int b, fmpz_poly_t tb_y, PKE_Para* pke_para,
              t_b* ek, fmpz_t coeff, fmpz_poly_t** X_c, fmpz_poly_t** Xs_c, fmpz_poly_t *PRF){

    t_b* tb=malloc(sizeof(*tb));
    fmpz_poly_init(tb->s1);
    fmpz_poly_init(tb->s2);
    t_b* tb2=malloc(sizeof(*tb2));
    fmpz_poly_init(tb2->s1);
    fmpz_poly_init(tb2->s2);


    fmpz_poly_scalar_mul_fmpz(tb2->s1,ek->s1,coeff);
    fmpz_poly_scalar_mul_fmpz(tb2->s2,ek->s2,coeff);
    VHSS_Mult(tb,pke_para,tb2,X_c[0],Xs_c[0]);
    if(b==1)
    {
        fmpz_poly_add(tb->s1,tb->s1,PRF[0]);
    } else{
        fmpz_poly_sub(tb->s1,tb->s1,PRF[0]);
    }
    for(int i=1;i<pke_para->d;i++)
    {
        VHSS_Mult(tb2,pke_para,tb,X_c[i],Xs_c[i]);
        fmpz_poly_set(tb->s1,tb2->s1);
        fmpz_poly_set(tb->s2,tb2->s2);
        if(b==1)
        {
            fmpz_poly_add(tb->s1,tb->s1,PRF[1]);
        } else{
            fmpz_poly_sub(tb->s1,tb->s1,PRF[1]);
        }
    }
    fmpz_poly_set(tb_y,tb->s1);
}

void VHSS_Ver(VHSS_Para* vhss_para, fmpz_poly_t t1_y,fmpz_poly_t t2_y,fmpz_poly_t tau1_y,fmpz_poly_t tau2_y)
{
    fmpz_t y;
    fmpz_init(y);
    fmpz_poly_t t_y,tau_y,alpha_y;
    fmpz_poly_init(t_y);
    fmpz_poly_init(tau_y);
    fmpz_poly_init(alpha_y);
    fmpz_poly_add(t_y,t1_y,t2_y);
    fmpz_poly_get_coeff_fmpz(y,t_y,0);
    fmpz_poly_add(tau_y,tau1_y,tau2_y);
    fmpz_poly_scalar_mul_fmpz(alpha_y,vhss_para->alpha,y);
    if(fmpz_poly_equal(alpha_y,tau_y)==1)
    {
        printf("Verification passed!\n");
    } else{
        printf("Verification failed!\n");
    }
}

