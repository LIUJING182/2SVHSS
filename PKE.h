#ifndef ISC2019_ARB_PKE_H
#define ISC2019_ARB_PKE_H

#endif //ISC2019_ARB_PKE_H

#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

#include "GaussRand.h"
typedef struct {
    int N;
    int msg_bit;
    fmpz_t msg;
    int l_bit;
    fmpz_t r;
    int p_bit;
    fmpz_t p;
    int q_bit;
    fmpz_t q;
    fmpz_poly_t xN;

    flint_rand_t state;
    int d;
    int num_item;
    int num_data;
}PKE_Para;

typedef struct {
    fmpz_poly_t a;
    fmpz_poly_t b;
}PKE_PK;

typedef struct {
    fmpz_poly_t s1;
    fmpz_poly_t s2;
}PKE_SK;

void Set_Para(PKE_Para* pke_para);

void PKE_Gen(PKE_Para* pke_para,PKE_PK* pke_pk,PKE_SK* pke_sk,int l_bit,int d)
{
    //initialize
    fmpz_init(pke_para->p);
    fmpz_init(pke_para->q);
    fmpz_init(pke_para->msg);
    fmpz_init(pke_para->r);
    fmpz_poly_init(pke_para->xN);

    //assignment
    pke_para->l_bit=l_bit;
    pke_para->d=d;
    pke_para->num_data=pke_para->d;
    pke_para->msg_bit=pke_para->l_bit*(pke_para->d+1);
    Set_Para(pke_para);

    fmpz_set_d_2exp(pke_para->msg,1,pke_para->msg_bit);
    fmpz_set_d_2exp(pke_para->r,1,pke_para->l_bit);
    fmpz_set_d_2exp(pke_para->p,1,pke_para->p_bit);
    fmpz_set_d_2exp(pke_para->q,1,pke_para->q_bit);
    fmpz_poly_set_coeff_ui(pke_para->xN,pke_para->N,1);
    fmpz_poly_set_coeff_ui(pke_para->xN,0,1);

    //gen pk,sk
    fmpz_poly_t e,hat_s;
    fmpz_poly_init(pke_pk->a);
    fmpz_poly_init(pke_pk->b);
    fmpz_poly_init(e);
    fmpz_poly_init(hat_s);
    fmpz_poly_init(pke_sk->s1);
    fmpz_poly_init(pke_sk->s2);
    //pk
    flint_randinit(pke_para->state);
    fmpz_poly_randtest_not_zero(pke_pk->a,pke_para->state,pke_para->N,pke_para->q_bit);
    flint_randinit(pke_para->state);
    fmpz_poly_randtest_unsigned(hat_s,pke_para->state,pke_para->N,1);
    gaussrand(e);
    fmpz_poly_mul(pke_pk->b,pke_pk->a,hat_s);
    fmpz_poly_add(pke_pk->b,pke_pk->b,e);
    fmpz_poly_rem(pke_pk->b,pke_pk->b,pke_para->xN);
    fmpz_poly_scalar_smod_fmpz(pke_pk->b,pke_pk->b,pke_para->q);
    //sk
    fmpz_poly_set_coeff_ui(pke_sk->s1,0,1);
    fmpz_poly_set(pke_sk->s2,hat_s);
}

void PKE_Enc(fmpz_poly_t* x_c, PKE_Para* pke_para,PKE_PK* pke_pk, fmpz_t x)
{
    fmpz_poly_t v,e1,e2;
    fmpz_poly_init(v);
    fmpz_poly_init(e1);
    fmpz_poly_init(e2);
    flint_randinit(pke_para->state);
    fmpz_poly_randtest_not_zero(v,pke_para->state,pke_para->N,1);
    gaussrand(e1);
    gaussrand(e2);

    fmpz_t q_p;
    fmpz_init(q_p);
    fmpz_divexact(q_p,pke_para->q,pke_para->p);
    fmpz_mul(q_p,q_p,x);

    fmpz_poly_mul(x_c[0],pke_pk->b,v);
    fmpz_poly_add(x_c[0],x_c[0],e1);
    fmpz_poly_add_fmpz(x_c[0],x_c[0],q_p);

    fmpz_poly_mul(x_c[1],pke_pk->a,v);
    fmpz_poly_sub(x_c[1],e2,x_c[1]);

    fmpz_poly_rem(x_c[0],x_c[0],pke_para->xN);
    fmpz_poly_scalar_smod_fmpz(x_c[0],x_c[0],pke_para->q);
    fmpz_poly_rem(x_c[1],x_c[1],pke_para->xN);
    fmpz_poly_scalar_smod_fmpz(x_c[1],x_c[1],pke_para->q);
}

void PKE_OKDM(fmpz_poly_t* x_c, fmpz_poly_t* xs_c, PKE_Para* pke_para,PKE_PK* pke_pk, fmpz_t x)
{
    fmpz_t zero;
    fmpz_init_set_ui(zero,0);
    PKE_Enc(x_c,pke_para,pke_pk,x);
    PKE_Enc(xs_c,pke_para,pke_pk,zero);

    fmpz_t q_p;
    fmpz_init(q_p);
    fmpz_divexact(q_p,pke_para->q,pke_para->p);
    fmpz_mul(q_p,q_p,x);
    fmpz_poly_add_fmpz(xs_c[1],xs_c[1],q_p);
    fmpz_poly_scalar_smod_fmpz(xs_c[1],xs_c[1],pke_para->q);
}

void Round(PKE_Para* pke_para,fmpz_poly_t db,fmpq_t p_q)
{
    fmpz_t coff_fmpz;
    fmpq_t coff_fmpq;
    mpz_t coff_mpz;
    mpq_t coff_mpq;
    mpq_t one_div_two;

    fmpz_init(coff_fmpz);
    fmpq_init(coff_fmpq);
    mpz_init(coff_mpz);
    mpq_init(coff_mpq);
    mpq_init(one_div_two);

    mpq_set_d(one_div_two,0.5);

    for(int i=0;i<pke_para->N;i++)
    {
        fmpz_poly_get_coeff_fmpz(coff_fmpz,db,i);
        fmpq_mul_fmpz(coff_fmpq,p_q,coff_fmpz);
        fmpq_get_mpq(coff_mpq,coff_fmpq);
        if(mpq_sgn(coff_mpq)>0)
        {
            mpq_add(coff_mpq,coff_mpq,one_div_two);
        }
        else{
            mpq_sub(coff_mpq,coff_mpq,one_div_two);
        }
        mpz_set_q(coff_mpz,coff_mpq);
        fmpz_set_mpz(coff_fmpz,coff_mpz);
        fmpz_poly_set_coeff_fmpz(db,i,coff_fmpz);
    }

}

void PKE_DDec(fmpz_poly_t db_1, fmpz_poly_t db_2, PKE_Para* pke_para,fmpz_poly_t tb_1,fmpz_poly_t tb_2,
              fmpz_poly_t* x_c, fmpz_poly_t* xs_c)
{
    fmpz_poly_t temp1,temp2;
    fmpz_poly_init(temp1);fmpz_poly_init(temp2);
    fmpq_t p_q;
    fmpq_init(p_q);
    fmpq_set_fmpz_frac(p_q,pke_para->p,pke_para->q);

    fmpz_poly_mul(temp1,x_c[0],tb_1);
    fmpz_poly_mul(temp2,x_c[1],tb_2);
    fmpz_poly_add(db_1,temp1,temp2);
    fmpz_poly_rem(db_1,db_1,pke_para->xN);
    fmpz_poly_scalar_smod_fmpz(db_1,db_1,pke_para->q);
    Round(pke_para,db_1,p_q);
    fmpz_poly_scalar_smod_fmpz(db_1,db_1,pke_para->p);

    fmpz_poly_mul(temp1,xs_c[0],tb_1);
    fmpz_poly_mul(temp2,xs_c[1],tb_2);
    fmpz_poly_add(db_2,temp1,temp2);
    fmpz_poly_rem(db_2,db_2,pke_para->xN);
    fmpz_poly_scalar_smod_fmpz(db_2,db_2,pke_para->q);
    Round(pke_para,db_2,p_q);
    fmpz_poly_scalar_smod_fmpz(db_2,db_2,pke_para->p);

    fmpz_poly_clear(temp1);
    fmpz_poly_clear(temp2);
    fmpq_clear(p_q);
}


void Set_Para(PKE_Para* pke_para)
{
    switch (pke_para->l_bit) {
        case 16:
            pke_para->N = 32768;
            pke_para->p_bit = 319;
            pke_para->q_bit = 662;
            break;
        case 32:
            pke_para->N = 65536;
            pke_para->p_bit = 576;
            pke_para->q_bit = 1177;
            break;
        case 64:
            pke_para->N = 131072;
            pke_para->p_bit = 1088;
            pke_para->q_bit = 2204;
            break;
        default:
            printf("error\n");
            break;
    }
}
