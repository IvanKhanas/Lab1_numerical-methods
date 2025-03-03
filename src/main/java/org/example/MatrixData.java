package org.example;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;

public abstract class MatrixData {
    public RealMatrix A = new Array2DRowRealMatrix(new double[][]{{12,7,4,1}, {5,20,8,7}, {1,4,14,8}, {0,1,4,5}}, false);

    public double[] Accurence = {0.01, 0.001, 0.0001};


    public RealMatrix E = new Array2DRowRealMatrix(new double[][]{{1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}}, false);

     public RealMatrix AInverse = new LUDecomposition(A).getSolver().getInverse(); //получаем обратную матрицу для A

     DoubleMatrix Ain = new DoubleMatrix(AInverse.getData());

     public double[] getAEigenvalueModules(){ //вычисление собственных значений матрицы A
         DoubleMatrix Aother = new DoubleMatrix(A.getData());
         ComplexDoubleMatrix eigenvalues = Eigen.eigenvalues(Aother);

         double[] eigenvalueModules = new double[eigenvalues.length];

         for (int i = 0; i < eigenvalues.length; i++){
             eigenvalueModules[i] = eigenvalues.get(i).abs();
         }
         return eigenvalueModules;
     }

    public double[] getEigenvalueModules() {   // Метод получения собственных значений обратной матрицы
        ComplexDoubleMatrix eigenvalues = Eigen.eigenvalues(Ain);
        double[] eigenvalueModules = new double[eigenvalues.length];

        for (int i = 0; i < eigenvalues.length; i++) {
            eigenvalueModules[i] = eigenvalues.get(i).abs(); // получаем модуль комплексного числа
        }

        return eigenvalueModules; // возвращаем массив модулей
    }


}
