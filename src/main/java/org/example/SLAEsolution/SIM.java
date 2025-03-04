package org.example.SLAEsolution;


import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
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

import static java.lang.Math.*;
import static org.apache.commons.math3.stat.StatUtils.max;

public class SIM extends MatrixData {

    double Tao;



    public List<Double> NormsAnswers [] = new List[3];
    public List<Double> IterationsAnswers[] =new List[3];


    int NumOFx_0;



    public void convergence () { //Проверка критерия сходимости

        System.out.println("\n\n" + "                    Решение методом простой итерации: ");

        System.out.println("Проверка по критерию сходимости :");
        for ( double t=1.00; t>=0.00; t-=0.01){
            RealMatrix result= E.subtract(A.scalarMultiply(t));

            DoubleMatrix res1 = new DoubleMatrix(result.getData());
            ComplexDoubleMatrix lambdas = Eigen.eigenvalues(res1);
            double[] lambdasModules = new double[lambdas.length];

            for (int i = 0; i < lambdas.length; i++) {
                lambdasModules[i] = lambdas.get(i).abs(); // получаем модуль комплексного числа
            }

            double norm = sqrt(max(lambdasModules));

          if (norm<1.00) { Tao = t;
                System.out.println("True! " + "Значение параметра тао : "+ t + "\n");
                System.out.println("Значение выражения ||E-Tao*A|| : "+norm); break;}

        }


    }

    public void Iteration(){ //Нахождение решений
        convergence();

        List<RealVector> VectorList = new ArrayList<>();
        Collections.addAll(VectorList, x0_1, x0_2, x0_3, x0_4, x0_5);

        NumOFx_0 = VectorList.size();

        List<RealVector> Answers = new ArrayList<>();

        System.out.println("\nРешения по этому методу:");
        System.out.println("Значение x точное : \t" + XAccurate);



        System.out.println("Выбранные вектора: \t" + VectorList+ "\n\n");

        for (int j=0; j<Accurence.length; j++){

            List<Double> NormsList = new ArrayList<>();
            List<Double> NumsOfIterations = new ArrayList<>();

        for (int i=0; i<VectorList.size(); i++){
            RealVector x1 = ((E.subtract(A.scalarMultiply(Tao))).operate(VectorList.get(i))).add(b.mapMultiply(Tao));

            RealVector x2 = ((E.subtract(A.scalarMultiply(Tao))).operate(x1).add(b.mapMultiply(Tao)));

            double NumOfIterations = 0.0;

            double norm = 0.0;

            while(true){
                RealVector x3 = ((E.subtract(A.scalarMultiply(Tao))).operate(x2).add(b.mapMultiply(Tao)));
                if((x3.subtract(x2)).getNorm() < Accurence[j]){
                    Answers.add(x3);

                    NumsOfIterations.add(NumOfIterations);

                    NormsList.add((XAccurate.subtract(VectorList.get(i))).getNorm());
                    norm = (XAccurate.subtract(VectorList.get(i))).getNorm();

                    break;
                }

                x2 = x3;
                NumOfIterations++;

            }


            System.out.println("Вектор номер " +(i+1)+" : "+"\t" +VectorList.get(i)+"\n"  + "Точность :" + Accurence[j]);
            System.out.println("Решение "  + " : " + Answers + "\tКоличество итераций: " + NumOfIterations );
            System.out.println("||x_точн-x_0|| : " + (norm) + "\n");
            NumOfIterations = 0;
            norm = 0.0;
            Answers.clear();

        }
            Collections.sort(NormsList);

            IterationsAnswers[j] = new ArrayList<>(NumsOfIterations);

            NormsAnswers[j] =(new ArrayList<>(NormsList));

    }

    }

    public void PrintGRaphs(){ //Рисуем графики

        JFrame frame = new JFrame("Graphs");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setLayout(new GridLayout(1, 3));




        XYChart chart1 = new XYChart(600, 600 );
        chart1.setTitle("SIM E=0.01");
        chart1.setXAxisTitle("Number of Iterations");
        chart1.setYAxisTitle("||x_точн-x_0||");

        double  YData [] = (IterationsAnswers[0]).stream().mapToDouble(Double::doubleValue).toArray();

        double XData [] = (NormsAnswers[0]).stream().mapToDouble(Double::doubleValue).toArray();

        chart1.addSeries("Зависимость ||x_точн-x_0|| от числа итераций", XData, YData);
        frame.add(new XChartPanel(chart1));


        XYChart chart2 = new XYChart(600, 600 );
        chart2.setTitle("SIM E=0.001");
        chart2.setXAxisTitle("Number of Iterations");
        chart2.setYAxisTitle("||x_точн-x_0||");

        double YData2 [] = (IterationsAnswers[1]).stream().mapToDouble(Double::doubleValue).toArray();

        double XData2 [] = (NormsAnswers[1]).stream().mapToDouble(Double::doubleValue).toArray();

        chart2.addSeries("Зависимость ||x_точн-x_0|| от числа итераций", XData2, YData2);
        frame.add(new XChartPanel<>(chart2));

        XYChart chart3 = new XYChart(600, 600 );
        chart3.setTitle("SIM E=0.0001");
        chart3.setXAxisTitle("Number of Iterations");
        chart3.setYAxisTitle("||x_точн-x_0||");

        double YData3 [] = (IterationsAnswers[2]).stream().mapToDouble(Double::doubleValue).toArray();

        double XData3 [] = (NormsAnswers[2]).stream().mapToDouble(Double::doubleValue).toArray();

        chart3.addSeries("Зависимость ||x_точн-x_0|| от числа итераций", XData3, YData3);
        frame.add(new XChartPanel<>(chart3));

        frame.pack();
        frame.setVisible(true);


    }


}
