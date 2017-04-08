package com.vvolkova.rwg;

import com.vvolkova.geometric.Edge;
import com.vvolkova.geometric.Triangle;

import java.util.ArrayList;

public class RWGstructure {
    private ArrayList<Edge> innerEdges;
    private ArrayList<Double> edgesLenght;
    private ArrayList<RWGelement> RWGElements;
    private ArrayList<Triangle> trianglesPlus;
    private ArrayList<Triangle> trianglesMinus;


    public void setInnerEdges(ArrayList<Edge> innerEdges) {
        if(innerEdges != null) {
            this.innerEdges = innerEdges;
        }
    }

    public void setRWGElements(ArrayList<RWGelement> RWGElements) {
        if(RWGElements != null) {
            this.RWGElements = RWGElements;
        }
    }

    public void setEdgesLenght(ArrayList<Double> edgesLenght) {
        if(edgesLenght != null) {
            this.edgesLenght = edgesLenght;
        }
    }

    public void setTrianglesPlus(ArrayList<Triangle> trianglesPlus) {
        if(trianglesPlus != null) {
            this.trianglesPlus = trianglesPlus;
        }
    }

    public void setTrianglesMinus(ArrayList<Triangle> trianglesMinus) {
        if(trianglesMinus != null) {
            this.trianglesMinus = trianglesMinus;
        }
    }

    public ArrayList<Edge> getInnerEdges() {
        return innerEdges;
    }

    public ArrayList<RWGelement> getRWGElements() {
        return RWGElements;
    }

    public ArrayList<Double> getEdgesLenght() {
        return edgesLenght;
    }

    public ArrayList<Triangle> getTrianglesPlus() {
        return trianglesPlus;
    }

    public ArrayList<Triangle> getTrianglesMinus() {
        return trianglesMinus;
    }

    public void check() {
        System.out.println("INNER EDGES");
        for (int i = 0; i < innerEdges.size(); i++) {
            Edge edge = innerEdges.get(i);
            System.out.println(i + ": " + edge.getA() + " " + edge.getB());
        }

        System.out.println("RWG ELEMENTS");
        for (int i = 0; i < RWGElements.size(); i++) {
            RWGelement rwgElem = RWGElements.get(i);
            System.out.println(i + ": minus_" + rwgElem.getMinusTriangle() + " plus_" + rwgElem.getPlusTriangle()
                    + " innerEdge_" + rwgElem.getInnerEdge());
        }
    }
}
