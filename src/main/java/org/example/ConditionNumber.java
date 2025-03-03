package org.example;


import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static java.lang.Math.sqrt;
import static org.apache.commons.math3.stat.StatUtils.max;
import static org.apache.commons.math3.stat.StatUtils.sum;

public class ConditionNumber extends  MatrixData{

    List<Double> M = new ArrayList<>(); // Числа обусловленности

    double A1, A2, A3;

    double IA1, IA2, IA3;


    public void ANormas(){

        A1 = sqrt(max(getAEigenvalueModules()));

        IA1 = sqrt(max(getEigenvalueModules()));


        //тут будет A1 инвертированной

        List<Double> A2Rows = new ArrayList<>(); //Хранит суммы элементов по строкам для матрицы A

        List<Double> IA2Rows = new ArrayList<>(); //Хранит суммы элементов по строкам для Inverse A

        List<Double> IA3Cols = new ArrayList<>();//Хранит суммы элементов по столбцам для Inverse A

        List<Double> A3Cols = new ArrayList<>(); //Хранит суммы элементов по столбцам для A

        for(int i=0;i<4;i++){
              A2Rows.add(sum(A.getRow(i))); // добавляется макс значение каждой строки для матрицы A

              IA2Rows.add(sum(AInverse.getRow(i))); // добавляется макс значение каждой строки для обратной матрицы A

              IA3Cols.add(sum(A.getColumn(i))); // добавляется макс значение каждого столбца для матрицы A

              A3Cols.add(sum(AInverse.getColumn(i))); // добавляется макс значение каждого столбца для обратной матрицы A

          }

        // присваеваем максимальные значения листов

        A2 = Collections.max(A2Rows);

        IA2 = Collections.max(IA2Rows);

        A3 = Collections.max(A3Cols);

        IA3 = Collections.max(IA3Cols);

    }

    void PrintResults () {
        ANormas();

        List<Double> As = new ArrayList<>();

        List<Double> IAs = new ArrayList<>();

        Collections.addAll(As, A1, A2, A3);
        Collections.addAll(IAs, IA1, IA2, IA3);
        Collections.addAll(M, IA1*A1, A2*IA2,A3*IA3);

        System.out.println("\n" +"Числа обусловленности используя все 3 типа норм: \n" + M + "\n");

        System.out.println("                    Нормы матриц : " + "\n");

        System.out.println(" Нормы матрицы A: \n" + As +"\n"+"\n"+ "Нормы для обратной матрицы A: \n"+ IAs);
    }

    }
