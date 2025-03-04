package org.example.SLAEsolution;


import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
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

public class GaussSeidel extends Jacobi{



    public List<Double> NormsAnswers [] = new List[3];
    public List<Double> IterationsAnswers[] =new List[3];


    @Override
    public void convergence(){ //Проверка по критерию сходимости
        System.out.println("                    Решение методом Гаусса-Зейделя :\n");
        fillMatrix();

        System.out.println("Проверка по критерию сходимости : ");
        System.out.println("метод сходится тогда и только тогда, когда собственные значения\n" +
                "матрицы − (L+D)^-1 * U   по модулю не превосходят 1.\n");


        RealMatrix AlmostFinalMAtrix = new LUDecomposition(L.add(D)).getSolver().getInverse();
        RealMatrix FinalMatrix = AlmostFinalMAtrix.multiply(U);


        DoubleMatrix Aother = new DoubleMatrix(FinalMatrix.getData());
        ComplexDoubleMatrix eigenvalues = Eigen.eigenvalues(Aother);

        double[] eigenvalueModules = new double[eigenvalues.length];

        int points = 0;

        for (int i = 0; i < eigenvalues.length; i++){
            eigenvalueModules[i] = eigenvalues.get(i).abs();
        }

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

    @Override
    public void Iteration(){ //Нахождение решений
        convergence();

        RealMatrix Inv = new LUDecomposition(L.add(D)).getSolver().getInverse();



        Collections.addAll(VectorList, x0_1, x0_2, x0_3,x0_4,x0_5);

        List<RealVector> Answers = new ArrayList<>();

        System.out.println("\nРешения по этому методу:");
        System.out.println("Значение x точное : \t" + XAccurate);

        System.out.println("Выбранные вектора: \t" + VectorList+ "\n\n");

        for (int j=0; j<Accurence.length;j++) {
            List<Double> NormsList = new ArrayList<>();
            List<Double> NumsOfIterations = new ArrayList<>();

            for (int i = 0; i < VectorList.size(); i++) {
                RealVector x1 = Inv.operate(b).subtract((Inv.multiply(U)).operate(VectorList.get(i)));

                RealVector x2 = Inv.operate(b).subtract((Inv.multiply(U)).operate(x1));

                double NumOfIterations = 0;

                while (true) {
                    RealVector x3 = Inv.operate(b).subtract((Inv.multiply(U)).operate(x2));
                    if ((x3.subtract(x2)).getNorm() < Accurence[j]) {
                        Answers.add(x3);

                        NumsOfIterations.add(NumOfIterations);

                        NormsList.add((x3.subtract(VectorList.get(i))).getNorm());

                        break;
                    }

                    x2 = x3;
                    NumOfIterations++;


                }
                System.out.println("Вектор номер " +(i+1)+" : "+"\t" +VectorList.get(j)+"\n"  + "Точность :" + Accurence[j]);
                System.out.println("Решение "  + " : " + Answers + "\tКоличество итераций: " + NumOfIterations);
                System.out.println("||x_3-x_0|| : " + NormsList.get(i) + "\n");


                NumOfIterations = 0;
                Answers.clear();


            }
            IterationsAnswers[j] = new ArrayList<>(NumsOfIterations);

            Collections.sort(NormsList);

            NormsAnswers[j] =(new ArrayList<>(NormsList));
            System.out.println("\n");

    }



    }
    @Override
    public void PrintGRaphs(){ //Рисуем графики

        JFrame frame = new JFrame("Graphs");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setLayout(new GridLayout(1, 3));




        XYChart chart1 = new XYChart(600, 600 );
        chart1.setTitle("GaussSeidel E=0.01");
        chart1.setXAxisTitle("||x_точн-x_0||");
        chart1.setYAxisTitle("Number of Iterations");

        double YData [] = (IterationsAnswers[0]).stream().mapToDouble(Double::doubleValue).toArray();
        double XData [] = (NormsAnswers[0]).stream().mapToDouble(Double::doubleValue).toArray();

        chart1.addSeries("Зависимость числа итераций от ||x_точн-x_0||", XData, YData);
        frame.add(new XChartPanel(chart1));


        XYChart chart2 = new XYChart(600, 600 );
        chart2.setTitle("GaussSeidel E=0.001");
        chart2.setXAxisTitle("||x_точн-x_0||");
        chart2.setYAxisTitle("Number of Iterations");

        double YData2 [] = (IterationsAnswers[1]).stream().mapToDouble(Double::doubleValue).toArray();
        double XData2 [] = (NormsAnswers[1]).stream().mapToDouble(Double::doubleValue).toArray();

        chart2.addSeries("Зависимость числа итераций от ||x_точн-x_0||", XData2, YData2);
        frame.add(new XChartPanel<>(chart2));

        XYChart chart3 = new XYChart(600, 600 );
        chart3.setTitle("GaussSeidel E=0.0001");
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
