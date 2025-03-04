package org.example.SLAEsolution;

import org.apache.commons.math3.linear.*;
import org.example.MatrixData;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import org.knowm.xchart.XChartPanel;
import org.knowm.xchart.XYChart;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static java.lang.Math.abs;

public class Jacobi extends MatrixData { //Создание и заполнение матриц L,U и D

    int n = A.getRowDimension();

    RealMatrix L = new Array2DRowRealMatrix(n, n);
    RealMatrix U = new Array2DRowRealMatrix(n, n);
    RealMatrix D = new Array2DRowRealMatrix(n, n);

    public List<Double> NormsAnswers [] = new List[3];
    public List<Double> IterationsAnswers[] =new List[3];




    private static void printMatrix(RealMatrix matrix) { //Вывод матриц

        int rows = matrix.getRowDimension();
        int cols = matrix.getColumnDimension();

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                System.out.printf("%8.2f ", matrix.getEntry(i, j));
            }
            System.out.println();
        }
        System.out.println();
    }

        public void fillMatrix(){

        // Заполнение матриц L и U
        for (int i = 0; i < n; i++) {
            // Заполнение верхней треугольной матрицы U
            for (int j = i; j < n; j++) {
                U.setEntry(i, j, A.getEntry(i, j) );
            }

            // Заполнение нижней треугольной матрицы L
            for (int j = i; j < n; j++) {
                if (i == j) {
                    L.setEntry(i, j, 0); // Главная диагональ равна 1
                } else {

                    L.setEntry(j, i, (A.getEntry(j, i)));
                }
            }

        }

        // Заполнение диагональной матрицы D
        for (int i = 0; i < n; i++) {
            D.setEntry(i, i, A.getEntry(i, i)); // Диагональные элементы из U
            U.setEntry(i, i, 0); // Обнуляем диагональ в U
            L.setEntry(i, i, 0);
        }

            System.out.println("Матрица A:");

        printMatrix(A);

        System.out.println("Матрица L:");

        printMatrix(L);

        System.out.println("Матрица U:");

        printMatrix(U);

        System.out.println("Матрица D:");

        printMatrix(D);


    }



    public void convergence() { //Проверка по достаточному условию сходимости
        fillMatrix();
        System.out.println("\n" + "Проверка по достаточному условию сходимости: ");
        int points = 0;

        for (int i = 0; i < A.getRowDimension(); i++) {
            double elemsum = 0.0;

            for (int j = 0; j < A.getColumnDimension(); j++) {
                if (i != j) {
                    elemsum += abs(A.getEntry(i, j));
                }
            }

            if (abs(A.getEntry(i, i)) > elemsum) {
                points++;
                System.out.println(abs(A.getEntry(i, i)) + " > " + elemsum);

            } else if (abs(A.getEntry(i, i)) < elemsum) {
                System.out.println(abs(A.getEntry(i, i)) + " < " + elemsum);
            } else if (abs(A.getEntry(i, i)) == elemsum) {
                System.out.println(abs(A.getEntry(i, i)) + " = " + elemsum);

            }

        }

        if (points == A.getRowDimension()) {
            System.out.println("True! Матрица обладает диагональным преобладанием.");
        } else {
            System.out.println("False! Матрица не обладает диагональным преобладанием.");
        }
    }

    public void CriterionOfConvergence() { //Проверка по критерию сходимости
        convergence();
        RealMatrix Dinv = new LUDecomposition(D).getSolver().getInverse();
        RealMatrix FinalMatrix = Dinv.multiply(L.add(U));


        DoubleMatrix Aother = new DoubleMatrix(FinalMatrix.getData());
        ComplexDoubleMatrix eigenvalues = Eigen.eigenvalues(Aother);

        double[] eigenvalueModules = new double[eigenvalues.length];

        int points = 0;

        for (int i = 0; i < eigenvalues.length; i++){
            eigenvalueModules[i] = eigenvalues.get(i).abs();
        }
        System.out.println("\n Проверка по критерию сходимости: \n");
        System.out.println("(Все собственные значения D^-1 * (L+U) < 1)");

        for (int i = 0; i < eigenvalueModules.length; i++) {
            if (eigenvalueModules[i]<1.0){System.out.println(eigenvalueModules[i]+" < "+1);points++;}
            else if (eigenvalueModules[i]>1.0) {System.out.println(eigenvalueModules[i]+" > "+1);

            } else if (eigenvalueModules[i]==1.0) {
                System.out.println(eigenvalueModules[i]+" = "+1);
            }
        }
        if (points == eigenvalueModules.length) {System.out.println("True! Метод сходится!\n");}
        else{System.out.println("False! Метод не сходится!\n");}
    }



    public void Iteration(){ //Нахождение решений
        System.out.println("\n\n"+"                    Решение методом Якоби:"+ "\n");
        CriterionOfConvergence();
        RealMatrix Dinv = new LUDecomposition(D).getSolver().getInverse();


        List<RealVector> VectorList = new ArrayList<>();
        Collections.addAll(VectorList, x0_1, x0_2, x0_3,x0_4,x0_5);

        List<RealVector> Answers = new ArrayList<>();

        System.out.println("\nРешения по этому методу:");
        System.out.println("Значение x точное : \t" + XAccurate);

        System.out.println("Выбранные вектора: \t" + VectorList+ "\n\n");

        for (int j=0; j<Accurence.length;j++) {

            List<Double> NormsList = new ArrayList<>();
            List<Double> NumsOfIterations = new ArrayList<>();

            for (int i = 0; i < VectorList.size(); i++) {
                RealVector x1 = Dinv.operate(b.subtract((L.add(U)).operate(VectorList.get(i))));

                RealVector x2 = Dinv.operate(b.subtract((L.add(U)).operate(x1)));

                double NumOfIterations = 0;

                while (true) {
                    RealVector x3 = Dinv.operate(b.subtract((L.add(U)).operate(x2)));
                    if ((x3.subtract(x2)).getNorm() < Accurence[j]) {
                        Answers.add(x3);

                        NumsOfIterations.add(NumOfIterations);

                        NormsList.add((x3.subtract(VectorList.get(i))).getNorm());

                        break;
                    }

                    x2 = x3;
                    NumOfIterations++;


                }
                System.out.println("Вектор номер " +(i+1)+" : "+"\t" +VectorList.get(i)+"\n"  + "Точность :" + Accurence[j]);
                System.out.println("Решение "  + " : " + Answers + "\tКоличество итераций: " + NumOfIterations);
                System.out.println("||x_3-x_0|| : " + NormsList.get(i) + "\n");


                NumOfIterations = 0;
                Answers.clear();


            }
            Collections.sort(NormsList);

            IterationsAnswers[j] = new ArrayList<>(NumsOfIterations);

            NormsAnswers[j] =(new ArrayList<>(NormsList));
            System.out.println("\n");
        }



    }

    public void PrintGRaphs(){ //Рисуем графики

        JFrame frame = new JFrame("Graphs");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setLayout(new GridLayout(1, 3));




        XYChart chart1 = new XYChart(600, 600 );
        chart1.setTitle("Jacobi E=0.01");
        chart1.setXAxisTitle("||x_точн-x_0||");
        chart1.setYAxisTitle("Number of Iterations");

        double YData [] = (IterationsAnswers[0]).stream().mapToDouble(Double::doubleValue).toArray();
        double XData [] = (NormsAnswers[0]).stream().mapToDouble(Double::doubleValue).toArray();

        chart1.addSeries("Зависимость числа итераций от ||x_точн-x_0||", XData, YData);
        frame.add(new XChartPanel(chart1));


        XYChart chart2 = new XYChart(600, 600 );
        chart2.setTitle("Jacobi E=0.001");
        chart2.setXAxisTitle("||x_точн-x_0||");
        chart2.setYAxisTitle("Number of Iterations");

        double YData2 [] = (IterationsAnswers[1]).stream().mapToDouble(Double::doubleValue).toArray();
        double XData2 [] = (NormsAnswers[1]).stream().mapToDouble(Double::doubleValue).toArray();

        chart2.addSeries("Зависимость числа итераций от ||x_точн-x_0||", XData2, YData2);
        frame.add(new XChartPanel<>(chart2));

        XYChart chart3 = new XYChart(600, 600 );
        chart3.setTitle("Jacobi E=0.0001");
        chart3.setXAxisTitle("||x_точн-x_0||");
        chart3.setYAxisTitle("Number of Iterations");

        double YData3 [] = (IterationsAnswers[2]).stream().mapToDouble(Double::doubleValue).toArray();

        double XData3 [] = (NormsAnswers[2]).stream().mapToDouble(Double::doubleValue).toArray();

        chart3.addSeries("Зависимость числа итераций от ||x_точн-x_0||", XData3, YData3);
        frame.add(new XChartPanel<>(chart3));

        frame.pack();
        frame.setVisible(true);


    }

    }






