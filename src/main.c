#include "Matrix.c"
int main(){ 
    int en;
    printf("Enter the number of equations : ");
    scanf("%d",&en);
    Matrix A;
    Matrix F;
    init_Matrix(en , en , &A);
    init_Matrix(en , 1 , &F);
    for(int i = 1 ; i<= en ; i++){
        printf("### For Equation number %d \n",i);
        double item = 0;
        for(int j = 1 ; j<= en ; j++){
            printf("Enter the term of X%d : ",j);
            scanf("%lf",&item);
            setItem_Matrix(i , j , item , &A);
        }
        printf("For equation number %d enter the free term : ",i);
        scanf("%lf",&item);
        setItem_Matrix(i , 1 , item , &F);
    }
    if(determinant_Matrix(&A)==0){
        printf("There equations don't have solution");
    }
    else{
        Matrix X = solveEquations_Matrix(&A , &F);
        for(int i = 1 ; i<= en ; i++){
            printf("X%d = %0.2lf \n",i, getItem_Matrix(i , 1 , &X));
        }
    }
    system("pause");
    return 0;
}