package com.vvolkova.geometric;

public class Triangle{
    private int a;
    private int b;
    private int c;
    private int num;

    public Triangle(int a, int b, int c, int num) {
        this.num = num;
        this.a = a;
        this.b = b;
        this.c = c;
    }

    public int getC() {
        return c;
    }

    public void setC(int c) {
        this.c = c;
    }

    public int getB() {
        return b;
    }

    public void setB(int b) {
        this.b = b;
    }

    public int getA() {
        return a;
    }

    public void setA(int a) {
        this.a = a;
    }

    public int getNum() {
        return num;
    }

    public void setNum(int num) {
        this.num = num;
    }

    /*public double area() {
        //S = ((x1 - x3)(y2 - y3) - (x2 - x3)(y1 - y3)) / 2
        if(area == -1) {
            double m1 = (p1.getX() - p3.getX())*(p2.getY() - p3.getX());
            double m2 = (p2.getX() - p3.getX())*(p1.getY() - p3.getX());
            area = Math.abs(m1 - m2)/2;
        }
        return area;
    }

    public Point centre() {
        if(centre == null) {
            double x = (p1.getX() + p2.getX() + p3.getX())/3;
            double y = (p1.getY() + p2.getY() + p3.getY())/3;
            centre = new Point(-1, x, y, 0);
        }
        return centre;
    }*/
}