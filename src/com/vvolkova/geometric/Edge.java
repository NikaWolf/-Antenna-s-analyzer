package com.vvolkova.geometric;

public class Edge {
    private int a;
    private int b;
    private int num = -1;

    public Edge(int a, int b) {
        this.a = a;
        this.b = b;
    }

    public Edge(int a, int b, int num) {
        this.a = a;
        this.b = b;
        this.num = num;
    }

    public int getA() {
        return a;
    }

    public int getB() {
        return b;
    }

    public int getNum() {
        return num;
    }

    @Override
    public String toString() {
        return "Edge{" +
                "a=" + a +
                ", b=" + b +
                ", num=" + num +
                '}' + "\n";
    }
}
