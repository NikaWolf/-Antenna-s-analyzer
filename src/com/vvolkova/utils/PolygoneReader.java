package com.vvolkova.utils;

import com.vvolkova.geometric.Edge;
import com.vvolkova.geometric.Point;
import com.vvolkova.geometric.Polygone;
import com.vvolkova.geometric.Triangle;

import java.io.File;
import java.util.ArrayList;

public class PolygoneReader {
    private Scanner polygoneScanner;
    private Scanner meshScanner;
    private Scanner feedScanner;
    private Scanner hScanner;
    private Polygone polygone;

    public Polygone getPolygone() {
        return polygone;
    }

    public PolygoneReader(File poly, File ele) {
        polygoneScanner = new Scanner(poly);
        meshScanner = new Scanner(ele);
        polygone = new Polygone();
        readPoints();
        //readBoundaryEdges();
        readMeshTriangle();


    }

    public PolygoneReader(File poly, File ele, File feed, File h) {
        this(poly,ele);
        feedScanner = new Scanner(feed);
        hScanner = new Scanner(h);

        readFeed();
        readH();
        polygone.check();
    }



    private void readPoints() {
        int n = polygoneScanner.nextInt();
        int n2 = polygoneScanner.nextInt();
        //int n3 = polygoneScanner.nextInt();
        //int n4 = polygoneScanner.nextInt();
        ArrayList<Point> points = new ArrayList<Point>();
        //Point p = new Point(100,100,100,100);
        //points.add(p);
        for (int i = 0; i < n; i++) {
            //int num = polygoneScanner.nextInt();
            double x = polygoneScanner.nextDouble();
            double y = polygoneScanner.nextDouble();
            double z = polygoneScanner.nextDouble();
            //if(num == i) {
                Point point = new Point(x, y, z, i);
                points.add(point);
            //}
        }
        polygone.setPoints(points);
    }

    private void readBoundaryEdges() {
        int n = polygoneScanner.nextInt();
        int n2 = polygoneScanner.nextInt();
        ArrayList<Edge> boundaryEdges = new ArrayList<Edge>();
        for (int i = 0; i < n; i++) {
            //int num = polygoneScanner.nextInt();
            int a = polygoneScanner.nextInt();
            int b = polygoneScanner.nextInt();
            int area = polygoneScanner.nextInt();
            //if(num == i) {
                Edge edge = new Edge(a, b, i);
                boundaryEdges.add(edge);
            //}
        }
        polygone.setBoundaryEdges(boundaryEdges);
    }

    private void readMeshTriangle() {
        int n = meshScanner.nextInt();
        //int n2 = meshScanner.nextInt();
        //int n3 = meshScanner.nextInt();
        ArrayList<Triangle> mesh = new ArrayList<Triangle>();
        //Triangle tr = new Triangle(-100,-100,-100,-100);
        //mesh.add(tr);
        for (int i = 0; i < n; i++) {
            //int num = meshScanner.nextInt();
            int a = meshScanner.nextInt();
            int b = meshScanner.nextInt();
            int c = meshScanner.nextInt();
            int num = meshScanner.nextInt();
            if(num < 3) {
                Triangle triangle = new Triangle(a, b, c, num);
                mesh.add(triangle);
            }
        }
        polygone.setMesh(mesh);
    }

    private void readFeed(){
        int n = feedScanner.nextInt();
        ArrayList<Integer> feedTriangles = new ArrayList<Integer>();
        for (int i = 0; i < n; i++) {
            int a = feedScanner.nextInt();
            feedTriangles.add(a);
        }
        polygone.setFeedTriangles(feedTriangles);
    }

    private void readH() {
        double h = hScanner.nextDouble();
        polygone.setH(h);
    }


}
