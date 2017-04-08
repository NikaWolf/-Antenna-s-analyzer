package com.vvolkova.rwg;

import com.vvolkova.geometric.Edge;
import com.vvolkova.geometric.Point;
import com.vvolkova.geometric.Polygone;
import com.vvolkova.geometric.Triangle;

import java.util.ArrayList;

public class GeometryAnalyzer {
    private Polygone polygone;
    private RWGstructure rwGstructure;

    ArrayList<Point> trianglesCenters;
    ArrayList<ArrayList<Point>> trianglesSubCenters;
    ArrayList<Integer> edgeIndecator;

    public GeometryAnalyzer(Polygone polygone){
        if(polygone != null) {
            this.polygone = polygone;
            rwGstructure = new RWGstructure();
            fillRWGStructure();
            rwGstructure.check();
            calcualteTrianglesSubcentres();
            calculateEdgeIndecators();
        }
    }

    public ArrayList< ArrayList<Point>> getTrianglesSubCenters() {
        return trianglesSubCenters;
    }

    public ArrayList<Point> getTrianglesCenters() {
        return trianglesCenters;
    }

    public RWGstructure getRwGstructure() {
        return rwGstructure;
    }

    public Polygone getPolygone() {
        return polygone;
    }

    private void fillRWGStructure() {
        ArrayList<Triangle> mesh = polygone.getMesh();
        ArrayList<Edge> innerEdges = new ArrayList<Edge>();
        ArrayList<RWGelement> rwGElements = new ArrayList<RWGelement>();
        ArrayList<Triangle> trianglesPlus = new ArrayList<Triangle>();
        ArrayList<Triangle> trianglesMinus = new ArrayList<Triangle>();
        int innerEdgesNum = 0;
        for (int i = 0; i < mesh.size(); i++) {
            Triangle t1 = mesh.get(i);
            for (int j = i+1; j < mesh.size(); j++){
                Triangle t2 = mesh.get(j);
                Edge edge = foundCommonEdge(t1, t2, innerEdgesNum);
                if (edge != null) {
                    innerEdges.add(innerEdgesNum, edge);
                    RWGelement rwgElement = new RWGelement(i,j,innerEdgesNum, innerEdgesNum);
                    innerEdgesNum++;
                    rwGElements.add(rwgElement);
                    trianglesPlus.add(mesh.get(i));
                    trianglesMinus.add(mesh.get(j));
                }
            }
        }
        rwGstructure.setInnerEdges(innerEdges);
        rwGstructure.setRWGElements(rwGElements);
        rwGstructure.setEdgesLenght(calculateEdgesLenght());
        rwGstructure.setTrianglesPlus(trianglesPlus);
        rwGstructure.setTrianglesMinus(trianglesMinus);
    }

    private  ArrayList<Double> calculateEdgesLenght(){
        ArrayList<Point> points = polygone.getPoints();
        ArrayList<Double> edgesLengh = new ArrayList<Double>();
        for (Edge e: rwGstructure.getInnerEdges()) {
            Point a = points.get(e.getA());
            Point b = points.get(e.getB());
            double res = Math.sqrt(((b.getX() - a.getX())*(b.getX() - a.getX())) +
                    ((b.getY() - a.getY())*(b.getY() - a.getY())) +
                    ((b.getZ() - a.getZ())*(b.getZ() - a.getZ())));
            edgesLengh.add(res);
        }
        return edgesLengh;

    }

    private Edge foundCommonEdge(Triangle t1, Triangle t2, int num) {
        Edge edge = null;
        int n = 0;
        int a = 0;
        int b = 0;
        int[] t1t2 = {t1.getA(), t1.getB(), t1.getC(), t2.getA(), t2.getB(), t2.getC()};
        for (int i = 0; i < 3; i++) {
            for (int j = 3; j < t1t2.length; j++) {
                if(t1t2[i] == t1t2[j]){
                    n++;
                    if(n==1)
                        a = t1t2[i];
                    else
                        b = t1t2[i];
                }
            }
        }
        if(n==2)
            edge = new Edge(a, b, num);
        return edge;
    }

    private void calcualteTrianglesSubcentres() {
        ArrayList<Triangle> triangles= polygone.getMesh();
        ArrayList<Point> points = polygone.getPoints();

        trianglesCenters = new ArrayList<Point>();
        trianglesSubCenters = new ArrayList< ArrayList<Point>>();

        for (int i = 0; i < triangles.size(); i++) {
            Triangle t = triangles.get(i);

            Point a = points.get(t.getA());
            Point b = points.get(t.getB());
            Point c = points.get(t.getC());

            double x = (a.getX() + b.getX() + c.getX())/3;
            double y = (a.getY() + b.getY() + c.getY())/3;
            double z = (a.getZ() + b.getZ() + c.getZ())/3;

            Point center = new Point(x, y, z, -1);
            trianglesCenters.add(center);
            trianglesSubCenters.add(calculate_9_SubCentres(a, b,c,center));
        }
    }

    private ArrayList<Point> calculate_9_SubCentres(Point a, Point b, Point c, Point center) {
        ArrayList<Point> subCenters = new  ArrayList<Point>();
        Point ab = b.minus(a);
        Point bc = c.minus(b);
        Point ac = c.minus(a);

        double s = 1.0/3;
        double s2 = 2.0/3;

        Point c1 = a.plus(ab.mult(s ));
        Point c2 = a.plus(ab.mult(s2));
        Point c3 = b.plus(bc.mult(s ));
        Point c4 = b.plus(bc.mult(s2));
        Point c5 = a.plus(ac.mult(s ));
        Point c6 = a.plus(ac.mult(s2));

        subCenters.add(c1.plus(c5).plus(a).mult(s));
        subCenters.add(c1.plus(c2).plus(center).mult(s));
        subCenters.add( c2.plus(c3).plus(b).mult(s));
        subCenters.add(c2.plus(c3).plus(center).mult(s));
        subCenters.add(c3.plus(c4).plus(center).mult(s));
        subCenters.add(c1.plus(c5).plus(center).mult(s));
        subCenters.add(c5.plus(c6).plus(center).mult(s));
        subCenters.add(c4.plus(c6).plus(center).mult(s));
        subCenters.add(c4.plus(c6).plus(c).mult(s));
        return subCenters;
    }


    private void calculateEdgeIndecators() {
        edgeIndecator = new ArrayList<Integer>();
        for (int i = 0; i < rwGstructure.getRWGElements().size(); i++) {
            RWGelement el = rwGstructure.getRWGElements().get(i);
            Triangle tp = polygone.getMesh().get(el.getPlusTriangle());
            Triangle tm = polygone.getMesh().get(el.getMinusTriangle());
            edgeIndecator.add(tp.getNum()+tm.getNum());
        }
        System.out.println();
    }
}
