package com.vvolkova.rwg;

public class RWGelement {
    private int plusTriangle;
    private int minusTriangle;
    private int innerEdge;
    private int num;

    public RWGelement(int plus, int minus, int innerEdge, int num) {
        this.plusTriangle = plus;
        this.minusTriangle = minus;
        this.innerEdge = innerEdge;
        this.num = num;
    }

    public int getPlusTriangle() {
        return plusTriangle;
    }

    public int getInnerEdge() {
        return innerEdge;
    }

    public int getMinusTriangle() {
        return minusTriangle;
    }
}
