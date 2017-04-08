package com.vvolkova.geometric;

import java.util.ArrayList;

public class Polygone {
    private ArrayList<Point> points;
    private ArrayList<Edge> boundaryEdges;
    private ArrayList<Triangle> mesh;
    private ArrayList<Integer> feedTriangles;
    private double h;

    public void setPoints (ArrayList<Point> points) {
        if(points != null) {
            this.points = points;
        }
    }

    public void setBoundaryEdges (ArrayList<Edge> boundaryEdges) {
        if(boundaryEdges != null) {
            this.boundaryEdges = boundaryEdges;
        }
    }

    public void setMesh (ArrayList<Triangle> mesh) {
        if (mesh != null ) {
            this.mesh = mesh;
        }
    }

    public void setFeedTriangles(ArrayList<Integer> feedTriangles) {
        this.feedTriangles = feedTriangles;
    }

    public void setH(double h) {
        this.h = h;
    }

    public ArrayList<Triangle> getMesh() {
        return mesh;
    }

    public ArrayList<Point> getPoints() {
        return points;
    }

    public ArrayList<Integer> getFeedTriangles() {
        return feedTriangles;
    }

    public double getH() {
        return h;
    }

    public void check() {
        System.out.println("POINTS:");
        for (int i = 0; i < points.size(); i++) {
            Point p = points.get(i);
            System.out.println(i + ": " + p.getX() + " " + p.getY()+ " " + p.getZ());
        }
        System.out.println("BOUNDARY:");
        /*for (int i = 0; i < boundaryEdges.size(); i++) {
            Edge e = boundaryEdges.get(i);
            System.out.println(i + ": " + e.getA() + " " + e.getB());
        }*/
        System.out.println("MESH TRIANGLES:");
        for (int i = 0; i < mesh.size(); i++) {
            Triangle t = mesh.get(i);
            System.out.println(i + ": " + t.getA() + " " + t.getB() + " " + t.getC() + " " + t.getNum());
        }
        System.out.println("FEED TRIANGLES:");
        for (int i = 0; i < feedTriangles.size(); i++) {
            System.out.println(feedTriangles.get(i));
        }
        System.out.println("H:" + h);

    }
}
