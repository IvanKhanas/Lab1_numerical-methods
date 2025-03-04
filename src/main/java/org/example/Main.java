package org.example;


import org.example.SLAEsolution.GaussSeidel;
import org.example.SLAEsolution.Jacobi;
import org.example.SLAEsolution.SIM;



public class Main {
    public static void main(String[] args) {




        ConditionNumber cond  = new ConditionNumber();
        cond.PrintResults();



        SIM sim = new SIM();
        sim.Iteration();
        sim.PrintGRaphs();


        Jacobi jacobi = new Jacobi();
        jacobi.Iteration();
        jacobi.PrintGRaphs();


        GaussSeidel gaussSeidel = new GaussSeidel();
        gaussSeidel.Iteration();
        gaussSeidel.PrintGRaphs();








    }
}