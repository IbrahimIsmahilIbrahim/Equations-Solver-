#ifndef _Matrix_H_
#define _Matrix_H_

#include <stdio.h>
#include <stdlib.h>
typedef struct Matrix{
    double** data;
    int numOfRows;
    int numOfColumn;
}Matrix;

void init_Matrix(int nr , int nc , Matrix* m){
    m->data = (double**) malloc(nr*sizeof(double));
    m->numOfRows = nr;
    m->numOfColumn = nc;
    for(int i = 0 ; i<nr ; i++){
        m->data[i] = (double*) malloc(nc*sizeof(double));
    } 
    //assign values
    for(int i = 0 ; i<nr ; i++){
        for(int j = 0  ; j<nc ; j++){
            m->data[i][j] = 0.0;
        }
    }
}

void delete_Matrix(Matrix* m){
    free(m->data);
}

void setItem_Matrix(int nr , int nc , double newValue , Matrix* m){
    if(m->numOfRows >= nr && m->numOfColumn >= nc){
        m->data[nr-1][nc-1] = newValue;
    }
}

double getItem_Matrix(int nr , int nc , Matrix* m){
    return (m->numOfRows >= nr && m->numOfColumn >= nc) ? m->data[nr-1][nc-1] : -1;
}

void copy_Matrix(Matrix* m1 , Matrix* m2){
    init_Matrix(m2->numOfRows , m2->numOfColumn , m1);
    for(int i = 1 ; i<=m2->numOfRows ; i++){
        for(int j = 1 ; j<=m2->numOfColumn ; j++){
            setItem_Matrix(i , j , getItem_Matrix(i , j , m2) , m1);
        }
    }
}



void show_Matrix(Matrix* m){
    for(int i = 0 ; i<m->numOfRows ; i++){
        for(int j = 0 ; j<m->numOfColumn ; j++){
            printf("%0.2lf \t" , m->data[i][j]);
        }
        printf("\n");
    }
}

void deleteRow_Matrix(int nr , Matrix* m){
    if(m->numOfRows >= nr && nr>=1){
        Matrix copy;
        copy_Matrix(&copy , m);
        //rebuild to m
        delete_Matrix(m);
        init_Matrix(copy.numOfRows-1 , copy.numOfColumn , m);
        int id = 1;
        for(int i = 1 ; i<= copy.numOfRows ; i++){
            if(i==nr) continue;
            for(int j = 1 ; j<= copy.numOfColumn ; j++){
                setItem_Matrix(id , j , getItem_Matrix(i , j , &copy) , m );
            }
            id++;
        }
        delete_Matrix(&copy);
    }
}

void deleteColumn_Matrix(int nc , Matrix* m){
    if(m->numOfColumn>=nc && nc>=1){
        Matrix copy;
        copy_Matrix(&copy , m);
        delete_Matrix(m);
        init_Matrix(copy.numOfRows , copy.numOfColumn-1 , m);
        int id = 1;
        for(int i = 1 ; i<=copy.numOfRows ; i++){
            id=1;
            for(int j = 1 ; j<=copy.numOfColumn ; j++){
                if(nc==j) continue;
                setItem_Matrix(i , id , getItem_Matrix(i , j , &copy) , m);
                id++;
            }
        }
        delete_Matrix(&copy);
    }
}

double determinant_Matrix(Matrix* m){
    if(m->numOfRows == 1) return getItem_Matrix(1,1,m);
    else{
        double result = 0.0;
        for(int nc = 1 ; nc<= m->numOfColumn ; nc++){
            Matrix copy;
            copy_Matrix(&copy , m);
            deleteRow_Matrix(1 , &copy);
            deleteColumn_Matrix(nc , &copy);
            int signal = (nc%2==1) ? 1 : -1;
            result += (getItem_Matrix(1 , nc , m) * signal * determinant_Matrix(&copy));
            delete_Matrix(&copy);
        }
        return result;
    }
}

Matrix transpose_Matrix(Matrix* m){
    Matrix result;
    init_Matrix(m->numOfColumn , m->numOfRows , &result);
    for(int nr = 1 ; nr<= m->numOfRows ; nr++){
        for(int nc = 1 ; nc<= m->numOfColumn ; nc++){
            setItem_Matrix(nc , nr , getItem_Matrix(nr , nc , m) , &result);
        }
    }
    return result;
}

Matrix adj_Matrix(Matrix* m){
    Matrix result;
    init_Matrix(m->numOfRows , m->numOfColumn , &result);
    for(int nr = 1 ; nr<=m->numOfRows ; nr++){
        for(int nc=1 ; nc<=m->numOfColumn ; nc++){
            Matrix copy;
            copy_Matrix(&copy , m);
            deleteRow_Matrix(nr , &copy);
            deleteColumn_Matrix(nc , &copy);
            int signal = ((nr%2==1 && nc%2==1)||(nr%2==0 && nc%2==0))?1:-1;
            setItem_Matrix(nr , nc , signal*determinant_Matrix(&copy) , &result);
            delete_Matrix(&copy);
        }
    }
    return transpose_Matrix(&result);
}

Matrix multMatrix2Number(Matrix* m , double num){
    Matrix re;
    init_Matrix(m->numOfRows , m->numOfColumn , &re);
    for(int nr = 1 ; nr<= m->numOfRows ; nr++){
        for(int nc = 1 ; nc<= m->numOfColumn ; nc++){
            setItem_Matrix(nr , nc , getItem_Matrix(nr , nc , m)*num , &re);
        }
    }
    return re;
}
Matrix inverse_Matrix(Matrix* m){
    Matrix c = adj_Matrix(m);
    return multMatrix2Number(&c,(1/determinant_Matrix(m)) );
}
Matrix multMatrix2Matrix(Matrix* m1 , Matrix* m2){
    Matrix result;
    init_Matrix(m1->numOfRows , m2->numOfColumn , &result);
    for(int nr = 1 ; nr<=m1->numOfRows ; nr++){
        for(int nc = 1 ; nc<=m2->numOfColumn ; nc++){
            double item = 0.0;
            for(int i = 1 ; i<= m1->numOfColumn ; i++){
                item += (getItem_Matrix(nr , i , m1)*getItem_Matrix(i , nc , m2));
            }
            setItem_Matrix(nr , nc , item , &result);
        }
    }
    return result;
}
Matrix solveEquations_Matrix(Matrix* A , Matrix* F){
    Matrix c = inverse_Matrix(A);
    return multMatrix2Matrix(&c , F);
}
#endif