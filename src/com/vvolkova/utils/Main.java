package com.vvolkova.utils;

import com.vvolkova.impedance.ImpedanceMatrix;
import com.vvolkova.rwg.GeometryAnalyzer;

import java.io.File;

public class Main {

    public static void main(String[] args) throws Exception {
        //PolygoneReader poly = new PolygoneReader(new File("Points.txt"), new File("Triangles.txt"));
        PolygoneReader poly = new PolygoneReader(new File("Points.txt"), new File("Triangles.txt"), new File("Feeding Triangle.txt"),
                new File("H.txt"));
        GeometryAnalyzer geometryAnalyzer = new GeometryAnalyzer(poly.getPolygone());
        ImpedanceMatrix impedanceMatrix = new ImpedanceMatrix(geometryAnalyzer);
        //MomSolver momSolver = new MomSolver(impedanceMatrix);





    }


}
