#ifndef ISC2019_ARB_GENDATA_H
#define ISC2019_ARB_GENDATA_H

#endif //ISC2019_ARB_GENDATA_H

typedef struct {
    fmpz_t* X;
    fmpz_t* coeff;
    fmpz_poly_t* X_bin;
    fmpz_poly_t* PRF;
    fmpz_poly_t** X_c;
    fmpz_poly_t** Xs_c;
}Data;

void Gen_Data(Data* data, PKE_Para* pke_para, PKE_PK* pke_pk)
{
    int x;
    data->X=malloc(sizeof(fmpz_t)*pke_para->num_data);
    data->X_bin=malloc(sizeof(fmpz_poly_t)*pke_para->num_data);
    data->PRF=malloc(sizeof(fmpz_poly_t)*10);
    data->X_c=malloc(sizeof(fmpz_poly_t)*pke_para->num_data);
    data->Xs_c=malloc(sizeof(fmpz_poly_t)*pke_para->num_data);
    for(int i=0;i<pke_para->num_data;i++)
    {
        data->X_c[i]=malloc(sizeof(fmpz_poly_t)*2);
        data->Xs_c[i]=malloc(sizeof(fmpz_poly_t)*2);

        fmpz_init(data->X[i]);
        fmpz_poly_init(data->X_bin[i]);
        fmpz_poly_init(data->X_c[i][0]);
        fmpz_poly_init(data->X_c[i][1]);
        fmpz_poly_init(data->Xs_c[i][0]);
        fmpz_poly_init(data->Xs_c[i][1]);
        x=rand()% (int)pow(2,pke_para->l_bit);
        fmpz_set_ui(data->X[i],x);
        PKE_OKDM(data->X_c[i],data->Xs_c[i],pke_para,pke_pk,data->X[i]);
    }
    for(int i=0;i<10;i++)
    {
        fmpz_poly_init(data->PRF[i]);
        fmpz_poly_randtest_not_zero(data->PRF[i],pke_para->state,pke_para->N,pke_para->q_bit);
        fmpz_poly_scalar_smod_fmpz(data->PRF[i],data->PRF[i],pke_para->q);
    }
}