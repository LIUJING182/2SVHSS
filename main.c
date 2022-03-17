#include <stdio.h>
#include <flint/fmpz_mod_poly.h>
#include "VHSS.h"
void Eval_Poly(int l_bit,int d);
int main() {
    int l_bit;
    int d;
    printf("Enter the bits (16, 32, 64) of the input:\n");
    scanf(" %d", &l_bit);
    printf("Enter the degree (2,3,...,10) of term: \n");
    scanf("%d", &d);
    printf("*****************running*****************\n");
    Eval_Poly(l_bit,d);
    return 0;
}

void Eval_Poly(int l_bit,int d)
{
    VHSS_Para* vhss_para = malloc(sizeof(*vhss_para));
    PKE_Para* pke_para = malloc(sizeof(*pke_para));
    PKE_PK* pke_pk=malloc(sizeof(*pke_pk));
    PKE_SK* pke_sk=malloc(sizeof(*pke_sk));
    t_b* ek1_s=malloc(sizeof(*ek1_s));
    t_b* ek1_as=malloc(sizeof(*ek1_as));
    t_b* ek2_s=malloc(sizeof(*ek2_s));
    t_b* ek2_as=malloc(sizeof(*ek2_as));
    fmpz_poly_init(ek1_s->s1);fmpz_poly_init(ek1_s->s2);
    fmpz_poly_init(ek1_as->s1);fmpz_poly_init(ek1_as->s2);
    fmpz_poly_init(ek2_s->s1);fmpz_poly_init(ek2_s->s2);
    fmpz_poly_init(ek2_as->s1);fmpz_poly_init(ek2_as->s2);
    fmpz_poly_t t1_y,t2_y,tau1_y,tau2_y;
    fmpz_poly_init(t1_y);
    fmpz_poly_init(t2_y);
    fmpz_poly_init(tau1_y);
    fmpz_poly_init(tau2_y);
    Data* data=malloc(sizeof(*data));
    HSS_Gen(vhss_para,pke_para,pke_pk,pke_sk,ek1_s,ek1_as,ek2_s,ek2_as,l_bit,d);
    Gen_Data(data,pke_para,pke_pk);
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_randbits(coeff,pke_para->state,pke_para->l_bit);
    clock_t S1_start, S1_finish, S2_start, S2_finish, C_start, C_finish;
    double S1_time, S2_time, C_time;
    S1_start = clock();
    VHSS_Eval(1,t1_y,pke_para,ek1_s,coeff,data->X_c,data->Xs_c,data->PRF);
    VHSS_Eval(1,tau1_y,pke_para,ek1_as,coeff,data->X_c,data->Xs_c,data->PRF);
    S1_finish = clock();
    S1_time = (double)(S1_finish - S1_start) / CLOCKS_PER_SEC;
    S2_start = clock();
    VHSS_Eval(2,t2_y,pke_para,ek2_s,coeff,data->X_c,data->Xs_c,data->PRF);
    VHSS_Eval(2,tau2_y,pke_para,ek2_as,coeff,data->X_c,data->Xs_c,data->PRF);
    S2_finish = clock();
    S2_time = (double)(S2_finish - S2_start) / CLOCKS_PER_SEC;
    printf("all server's running time: %f ms\n", (S1_time+S2_time)*1000);
    C_start = clock();
    VHSS_Ver(vhss_para,t1_y,t2_y,tau1_y,tau2_y);
    C_finish = clock();
    C_time = (double)(C_finish - C_start) / CLOCKS_PER_SEC*1000;
    printf("the client's running time: %f ms\n", C_time);
}
