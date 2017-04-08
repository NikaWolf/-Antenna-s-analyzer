package com.vvolkova.geometric;

public class Point {
    private double x,y,z;
    private int num;

    public Point(double x, double y, double z, int num) {
        this.num = num;
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public double getX() {
        return x;
    }

    public void setX(double x) {
        this.x = x;
    }

    public double getY() {
        return y;
    }

    public void setY(double y) {
        this.y = y;
    }

    public double getZ() {
        return z;
    }

    public void setZ(double z) {
        this.z = z;
    }

    public int getNum() {
        return num;
    }

    public void setNum(int num) {
        this.num = num;
    }

    public boolean equals(Point obj) {
        return (this.getX()==obj.getX() &&
                this.getY() == obj.getY() &&
                this.getZ() == obj.getZ());
    }

    public Point minus(Point p){
        double x = this.x - p.getX();
        double y = this.y - p.getY();
        double z = this.z - p.getZ();
        return new Point(x, y, z, -1);
    }

    public Point plus(Point p) {
        double x = this.x + p.getX();
        double y = this.y + p.getY();
        double z = this.z + p.getZ();
        return new Point(x, y, z, -1);
    }

    public Point mult(double d) {
        double x = this.x * d;
        double y = this.y * d;
        double z = this.z * d;
        return new Point(x, y, z, -1);
    }

    public Point mult(Point p) {
        double x = this.x * p.getX();
        double y = this.y * p.getY();
        double z = this.z * p.getZ();
        return new Point(x,y,z,-1);
    }

    public Point divide(double d) {
        double x = this.x / d;
        double y = this.y / d;
        double z = this.z / d;
        return new Point(x, y, z, -1);
    }

    @Override
    public String toString() {
        return "Point{" +
                "x=" + String.format("%.4f", x)+
                ", y=" + String.format("%.4f", y) +
                ", z=" + String.format("%.4f", z) +
                '}';
    }
}
