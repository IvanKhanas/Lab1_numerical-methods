package org.example;

import org.apache.commons.math3.linear.*;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;

import javax.swing.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public abstract class MatrixData {


/// ВПИШИТЕ СВОИ ЗНАЧЕНИЯ
    public RealMatrix A = new Array2DRowRealMatrix(new double[][]{{12,7,4,1}, {5,20,8,7}, {1,4,14,8}, {0,1,4,5}}, false); //Матрица A


    public RealVector XAccurate  = new ArrayRealVector(new double[]{111.0/155, 119.0/155,173.0/310, 9.0/5});//Значения X точное
    public RealVector b = new ArrayRealVector(new double[] {18.0, 36.0, 26.0, 12.0});//Вектор значений после =

    public RealVector x0_1 = new ArrayRealVector(new double[] {0.0, 10.0, 3.0, 0.0});
    public RealVector x0_2 = new ArrayRealVector(new double[] {0.0, 8.0, 1.0, 4.0});
    public RealVector x0_3 = new ArrayRealVector(new double[] {1.0, 6.0, 7.0, 9.0}); //Выбранные вектора
    public RealVector x0_4 = new ArrayRealVector(new double[] {0.0, 2.0, 5.0, 11.0});
    public RealVector x0_5 = new ArrayRealVector(new double[] {1.0, 2.0, 5.0, 14.0});


    /// ПОСЛЕ ТОГО КАК ВПИСАЛИ СВОИ ЗНАЧЕНИЯ МОЖЕТЕ ЗАПУСКАТЬ




    public List<RealVector> VectorList = new ArrayList<>();


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
    JFrame frame = new JFrame("Graphs");




}
