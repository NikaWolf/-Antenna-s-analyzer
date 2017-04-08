package com.vvolkova.impedance;

import com.vvolkova.geometric.Edge;
import com.vvolkova.geometric.Point;
import com.vvolkova.geometric.Polygone;
import com.vvolkova.geometric.Triangle;
import com.vvolkova.rwg.GeometryAnalyzer;
import com.vvolkova.rwg.RWGelement;
import com.vvolkova.rwg.RWGstructure;

import java.util.ArrayList;

public class ImpedanceMatrix {
    private final double f           = 3e8;
    private final int NumberOfSteps  = 21;
    private final double FreqStart   = 1.0e9;
    private final double FreqStop    = 6.0e9;
    private final double step         = (FreqStop-FreqStart)/(NumberOfSteps-1);
    private final double epsilon_    = 8.854e-012;
    private final double mu_         = 1.257e-006;

    private final double c_          = 1.0 / Math.sqrt(epsilon_ * mu_);


    private final double eta_        = Math.sqrt(mu_ / epsilon_);


    private final double pi          = Math.PI;
    private double omega       = 2.0 * pi * f;
    private double k           = omega / c_;
    private Complex K          = new Complex(0,k);
    private double Constant1   = mu_ / (4 * pi);
    private final Complex oneComp    = new Complex(1, 0);
    private Complex Constant2  = oneComp.divides(new Complex(0,4.0*pi*omega*epsilon_));// = 1.0 / (j * 4.0 * pi * omega * epsilon_);
    private double Factor      = 1.0 / 9.0;

    private ArrayList<Complex> FactorA;
    private ArrayList<Complex> FactorFi;

    private ArrayList<Edge> edges;
    private ArrayList<RWGelement> rwgElements;
    private ArrayList<Triangle> triangles;
    private ArrayList<Triangle> trianglesPlus;
    private ArrayList<Triangle> trianglesMinus;
    private ArrayList<Point> points;
    private ArrayList<Point> trianglesCenters;
    private ArrayList<ArrayList<Point>> trianglesSubCenters;
    private ArrayList<Double> edgesLenght;

    private ArrayList<Point> rhoPlus;
    private ArrayList<Point> rhoMinus;
    private ArrayList<ArrayList<Point>> rhoSubPlus;
    private ArrayList<ArrayList<Point>> rhoSubMinus;

    private ArrayList<Triangle> DP;
    private int M;
    private Point delta;
    private double h;
    private Point middle;
    private ArrayList<Integer> feedingTreangle;
    private ArrayList<Point> point;
    private ArrayList<ArrayList<Point>> IMT;



    public ImpedanceMatrix(Polygone polygone, RWGstructure rwgStruct, ArrayList<Point> trianglesCenters,
                           ArrayList<ArrayList<Point>> trianglesSubCenters) {
        this.edges = rwgStruct.getInnerEdges();
        this.rwgElements = rwgStruct.getRWGElements();
        this.triangles = polygone.getMesh();
        this.points = polygone.getPoints();
        this.trianglesCenters = trianglesCenters;
        this.trianglesSubCenters = trianglesSubCenters;
        this.edgesLenght = rwgStruct.getEdgesLenght();
        this.trianglesPlus = rwgStruct.getTrianglesPlus();
        this.trianglesMinus = rwgStruct.getTrianglesMinus();
        this.h = polygone.getH();
        this.feedingTreangle = polygone.getFeedTriangles();


        calculateRho();
        initConstants();
        //fillImpedanceMatrixWithZero();

        fillDP();
        FrequencySeries();


        //calculateImpedanceMatrix();
        //checkMatrixZ();
    }

    private void fillDP() {
        DP = new ArrayList<Triangle>();
        for (int i = 0; i < triangles.size(); i++) {
            Triangle tr = triangles.get(i);
            if(tr.getNum()==0)
                DP.add(tr);
        }
        M = DP.size();
        delta = new Point(0,0,1e-12,-1);
        middle = new Point(0, 0, h/2, -1);
        point = new ArrayList<Point>();
        IMT = new ArrayList<ArrayList<Point>>();


        for (int i = 0; i < M; i++) {
            point.add(trianglesCenters.get(i).plus(middle));
            ArrayList<Point> pp = new ArrayList<Point>();
            for (int j = 0; j < trianglesSubCenters.get(i).size(); j++) {
                Point p = trianglesSubCenters.get(i).get(j).plus(middle);
                pp.add(p);
            }
            IMT.add(pp);
        }

    }


    private void FrequencySeries() {
        ArrayList<Double> f = new ArrayList<Double>();
        for (int FF = 0; FF < NumberOfSteps; FF++) {
            f.add(FreqStart + step*(FF));//FF-1
            omega = 2*pi*f.get(FF);
            k =omega/c_;
            K = new Complex(0,k);
            Constant1   = mu_/(4.0*pi);
            Constant2   = oneComp.divides(new Complex(0,4.0*pi*omega*epsilon_));
            Factor      = 1.0 / 9.0;
            initConstants();
            //Metal-to-metal impedance matrix (SS)
            ArrayList<ArrayList<Complex>> ZSS = new ArrayList<ArrayList<Complex>>();
            ZSS = fillImpedanceMatrixWithZero(ZSS);
            ZSS = calculateImpedanceMatrix(ZSS);
            //Dielectric-to-dielectric impedance matrix (DD)
            ArrayList<ArrayList<Complex>> ZDD = new ArrayList<ArrayList<Complex>>();
            for (int i = 0; i < M; i++) {
                ArrayList<Complex> zdd = new ArrayList<Complex>();
                for (int j = 0; j < M; j++) {
                    zdd.add(new Complex(0,0));
                }
            }
            ArrayList<Point> E = new ArrayList<Point>();
            for (int i = 0; i < M; i++) {
                Point OP = point.get(i).plus(delta);
                //Point e =

            }



        }
        System.out.println();
    }

    public ImpedanceMatrix(GeometryAnalyzer geomAn) {
        this(geomAn.getPolygone(), geomAn.getRwGstructure(), geomAn.getTrianglesCenters(), geomAn.getTrianglesSubCenters());
    }

    private void initConstants() {
        FactorA = new ArrayList<Complex>(); //= Factor * (j * omega * EdgeLength / 4) * Constant1;
        FactorFi = new ArrayList<Complex>();//= Factor * EdgeLength * Constant2;
        for (int i = 0; i < edges.size(); i++) {
            Edge e = edges.get(i);
            double lenght = edgesLenght.get(i);
            FactorA.add( new Complex(0, (Factor * omega * lenght * Constant1) / 4.0));
            FactorFi.add(Constant2.times(Factor*lenght));
        }
    }


    private void calculateRho() {
        rhoPlus = new ArrayList<Point>();
        rhoMinus = new ArrayList<Point>();
        rhoSubPlus = new ArrayList<ArrayList<Point>>();
        rhoSubMinus = new ArrayList<ArrayList<Point>>();
        calculateRhoPlus();
        calculateRhoMinus();
    }

    private void calculateRhoPlus() {
        for (int i = 0; i < edges.size(); i++) {
            int NoPlus = rwgElements.get(i).getPlusTriangle();
            Triangle tr = triangles.get(NoPlus);
            Edge e = edges.get(i);

            int node;
            if (tr.getA() != e.getA() && tr.getA() != e.getB())
                node = tr.getA();
            else if (tr.getB() != e.getA() && tr.getB() != e.getB())
                node = tr.getB();
            else
                node = tr.getC();

            Point freeVertex = points.get(node);
            Point trCenter = trianglesCenters.get(NoPlus);
            Point pho = trCenter.minus(freeVertex);
            rhoPlus.add(pho);

            ArrayList<Point> phoSub = new ArrayList<Point>();
            ArrayList<Point> subPoints = trianglesSubCenters.get(NoPlus);
            for (int j = 0; j < 9; j++) {
                phoSub.add(subPoints.get(j).minus(freeVertex));
            }
            rhoSubPlus.add(i,phoSub);
        }
    }

    private void calculateRhoMinus() {
        for (int i = 0; i < edges.size(); i++) {
            int trNum = rwgElements.get(i).getMinusTriangle();
            Triangle tr = triangles.get(trNum);
            Edge e = edges.get(i);

            int node;
            if (tr.getA() != e.getA() && tr.getA() != e.getB())
                node = tr.getA();
            else if (tr.getB() != e.getA() && tr.getB() != e.getB())
                node = tr.getB();
            else
                node = tr.getC();

            Point freeVertex = points.get(node);
            Point trCenter = trianglesCenters.get(trNum);
            Point pho = freeVertex.minus(trCenter);
            rhoMinus.add(pho);

            ArrayList<Point> phoSub = new ArrayList<Point>();
            ArrayList<Point> subPoints = trianglesSubCenters.get(trNum);
            for (int j = 0; j < 9; j++) {
                //phoSub.add(subPoints.get(j).minus(freeVertex));
                phoSub.add(freeVertex.minus(subPoints.get(j)));
            }
            rhoSubMinus.add(phoSub);
        }
    }
    
    private ArrayList<ArrayList<Complex>> fillImpedanceMatrixWithZero(ArrayList<ArrayList<Complex>> Z){
        Z = new ArrayList<ArrayList<Complex>>();
        for (int i = 0; i < edges.size(); i++) {
            ArrayList<Complex> z = new ArrayList<Complex>();
            for (int j = 0; j < edges.size(); j++) {
                z.add(new Complex(0,0));
            }
            Z.add(z);
        }
        return Z;
    }

    private ArrayList<ArrayList<Complex>> calculateImpedanceMatrix(ArrayList<ArrayList<Complex>> Z) {
        //Z = new ArrayList<ArrayList<Complex>>();
        for (int i = 0; i < triangles.size(); i++) {
            Point triangleCenter = trianglesCenters.get(i);
            ArrayList<double[]> R = new ArrayList<double[]>();
            ArrayList<Complex[]> G = new ArrayList<Complex[]>();

            for (int j = 0; j < trianglesSubCenters.size(); j++) {
                ArrayList<Point> subCentres = trianglesSubCenters.get(j);
                double[] r = new double[9];
                Complex[] g = new Complex[9];
                for (int l = 0; l < subCentres.size(); l++) {
                    Point d = subCentres.get(l).minus(triangleCenter);
                    d = d.mult(d);
                    r[l] = Math.sqrt(d.getX() + d.getY() + d.getZ());//R=sqrt(sum(D.*D));
                    g[l] = K.times(-1).times(r[l]);
                    g[l] = g[l].exp().divides(new Complex(r[l],0));//g=exp(-K*R)./R;

                }
                R.add(r);
                G.add(g);
            }
            ArrayList<Complex[]> GP = new ArrayList<Complex[]>();
            ArrayList<Complex[]> GM = new ArrayList<Complex[]>();
            ArrayList<Complex> Fi = new ArrayList<Complex>();
            ArrayList<Complex> ZF = new ArrayList<Complex>();

            for (int j = 0; j <edges.size(); j++) {
                RWGelement rwgEl = rwgElements.get(j);
                Complex[] gp = new Complex[9];
                Complex[] gm = new Complex[9];
                Complex fi = new Complex(0, 0);
                for (int f = 0; f < 9; f++) {
                    gp[f] = G.get(rwgEl.getPlusTriangle())[f].clone();
                    gm[f] = G.get(rwgEl.getMinusTriangle())[f].clone();
                    fi = fi.plus(gp[f]).minus(gm[f]);
                }
                GM.add(gm);
                GP.add(gp);
                Fi.add(fi);
                ZF.add(fi.times(FactorFi.get(j)));
            }
            calculateZ_matrix(GP, GM, ZF, i, Z);
        }
        return Z;
    }

    private ArrayList<ArrayList<Complex>> calculateZ_matrix(ArrayList<Complex[]> GP, ArrayList<Complex[]> GM, ArrayList<Complex> ZF, int n, ArrayList<ArrayList<Complex>> Z) {
        ArrayList<Integer> PLUS = new ArrayList<Integer>();
        ArrayList<Integer> MINUS = new ArrayList<Integer>();
        for (int i = 0; i < edges.size(); i++) {
            RWGelement rwgEl = rwgElements.get(i);
            if(n == rwgEl.getPlusTriangle())
                PLUS.add(i);
            if(n == rwgEl.getMinusTriangle())
                MINUS.add(i);
        }

        for (int i = 0; i < PLUS.size(); i++) {
            int m = PLUS.get(i);
            ArrayList<ArrayList<Point>> RP = new ArrayList<ArrayList<Point>>();
            ArrayList<ArrayList<Complex>> RP_2 = new ArrayList<ArrayList<Complex>>();
            ArrayList<ArrayList<Complex>> RP_3 = new ArrayList<ArrayList<Complex>>();
            ArrayList<Complex> RP_4 = new ArrayList<Complex>();

            ArrayList<ArrayList<Point>> RM = new ArrayList<ArrayList<Point>>();
            ArrayList<ArrayList<Complex>> RM_2 = new ArrayList<ArrayList<Complex>>();
            ArrayList<ArrayList<Complex>> RM_3 = new ArrayList<ArrayList<Complex>>();
            ArrayList<Complex> RM_4 = new ArrayList<Complex>();

            ArrayList<Complex> A = new ArrayList<Complex>();
            ArrayList<Complex> Z1 = new ArrayList<Complex>();

            for (int j = 0; j < edges.size(); j++) {
                RP.add(rhoSubPlus.get(m));
                RM.add(rhoSubPlus.get(m));
                ArrayList<Point> RP1 = new ArrayList<Point>();
                ArrayList<Complex> RP2 = new ArrayList<Complex>();
                ArrayList<Complex> RP3 = new ArrayList<Complex>();
                Complex RP4 = new Complex(0,0);

                ArrayList<Point> RM1 = new ArrayList<Point>();
                ArrayList<Complex> RM2 = new ArrayList<Complex>();
                ArrayList<Complex> RM3 = new ArrayList<Complex>();
                Complex RM4 = new Complex(0,0);

                for (int q = 0; q < 9; q++) {
                    RP1.add(RP.get(j).get(q).mult(rhoPlus.get(j)));
                    RM1.add(RP.get(j).get(q).mult(rhoMinus.get(j)));
                }
                RP.set(j, RP1);
                RM.set(j, RP1);

                for (int q = 0; q < 9; q++) {
                    Complex cp = new Complex(0,0);
                    cp = cp.plus(RP.get(j).get(q).getX());
                    cp = cp.plus(RP.get(j).get(q).getY());
                    cp = cp.plus(RP.get(j).get(q).getZ());
                    RP2.add(cp);

                    Complex cm = new Complex(0,0);
                    cm = cm.plus(RM.get(j).get(q).getX());
                    cm = cm.plus(RM.get(j).get(q).getY());
                    cm = cm.plus(RM.get(j).get(q).getZ());
                    RM2.add(cm);
                }
                RP_2.add(j, RP2);
                RM_2.add(j, RP2);

                for (int q = 0; q < 9; q++) {
                    Complex cp = GP.get(j)[q].times(RP_2.get(j).get(q));
                    RP3.add(cp);

                    Complex cm = GM.get(j)[q].times(RM_2.get(j).get(q));
                    RM3.add(cm);
                }
                RP_3.add(RP3);
                RM_3.add(RM3);

                for (int q = 0; q < 9; q++) {
                    RP4 = RP4.plus(RP_3.get(j).get(q));
                    RM4 = RM4.plus(RM_3.get(j).get(q));
                }
                RP_4.add(RP4);
                RM_4.add(RM4);
                A.add(RP4.plus(RM4));
                Z1.add(FactorA.get(j).times(A.get(j)));
                Complex z = Z.get(j).get(m);
                z = z.plus((Z1.get(j).plus(ZF.get(j)).times(edgesLenght.get(m))));

                Z.get(j).set(m,z);
            }

        }

        for (int i = 0; i < MINUS.size(); i++) {
            int m = MINUS.get(i);
            ArrayList<ArrayList<Point>> RP = new ArrayList<ArrayList<Point>>();
            ArrayList<ArrayList<Complex>> RP_2 = new ArrayList<ArrayList<Complex>>();
            ArrayList<ArrayList<Complex>> RP_3 = new ArrayList<ArrayList<Complex>>();
            ArrayList<Complex> RP_4 = new ArrayList<Complex>();

            ArrayList<ArrayList<Point>> RM = new ArrayList<ArrayList<Point>>();
            ArrayList<ArrayList<Complex>> RM_2 = new ArrayList<ArrayList<Complex>>();
            ArrayList<ArrayList<Complex>> RM_3 = new ArrayList<ArrayList<Complex>>();
            ArrayList<Complex> RM_4 = new ArrayList<Complex>();

            ArrayList<Complex> A = new ArrayList<Complex>();
            ArrayList<Complex> Z1 = new ArrayList<Complex>();

            for (int j = 0; j < edges.size(); j++) {
                RP.add(rhoSubMinus.get(m));
                RM.add(rhoSubMinus.get(m));
                ArrayList<Point> RP1 = new ArrayList<Point>();
                ArrayList<Complex> RP2 = new ArrayList<Complex>();
                ArrayList<Complex> RP3 = new ArrayList<Complex>();
                Complex RP4 = new Complex(0,0);

                ArrayList<Point> RM1 = new ArrayList<Point>();
                ArrayList<Complex> RM2 = new ArrayList<Complex>();
                ArrayList<Complex> RM3 = new ArrayList<Complex>();
                Complex RM4 = new Complex(0,0);

                for (int q = 0; q < 9; q++) {
                    RP1.add(RP.get(j).get(q).mult(rhoPlus.get(j)));
                    RM1.add(RP.get(j).get(q).mult(rhoMinus.get(j)));
                }
                RP.set(j, RP1);
                RM.set(j, RP1);

                for (int q = 0; q < 9; q++) {
                    Complex cp = new Complex(0,0);
                    cp = cp.plus(RP.get(j).get(q).getX());
                    cp = cp.plus(RP.get(j).get(q).getY());
                    cp = cp.plus(RP.get(j).get(q).getZ());
                    RP2.add(cp);

                    Complex cm = new Complex(0,0);
                    cm = cm.plus(RM.get(j).get(q).getX());
                    cm = cm.plus(RM.get(j).get(q).getY());
                    cm = cm.plus(RM.get(j).get(q).getZ());
                    RM2.add(cm);
                }
                RP_2.add(j, RP2);
                RM_2.add(j, RP2);

                for (int q = 0; q < 9; q++) {
                    Complex cp = GP.get(j)[q].times(RP_2.get(j).get(q));
                    RP3.add(cp);

                    Complex cm = GM.get(j)[q].times(RM_2.get(j).get(q));
                    RM3.add(cm);
                }
                RP_3.add(RP3);
                RM_3.add(RM3);

                for (int q = 0; q < 9; q++) {
                    RP4 = RP4.plus(RP_3.get(j).get(q));
                    RM4 = RM4.plus(RM_3.get(j).get(q));
                }
                RP_4.add(RP4);
                RM_4.add(RM4);
                A.add(RP4.plus(RM4));
                Z1.add(FactorA.get(j).times(A.get(j)));
                //Complex z = Z.get(j).get(m);
                Complex z = (Z1.get(j).minus(ZF.get(j)).times(edgesLenght.get(m)));
                z = z.plus(Z.get(j).get(m));

                Z.get(j).set(m,z);
            }

        }
        return Z;

    }

    private void checkMatrixZ(ArrayList<ArrayList<Complex>> Z) {
        for (int i = 0; i < Z.size(); i++) {
            ArrayList<Complex> z = Z.get(i);
            for (int j = 0; j < z.size(); j++) {
                Complex z1 = z.get(j);
                System.out.print(z1 + " ");
            }
            System.out.println();
        }
    }

    public double getF() {
        return f;
    }

    public double getOmega() {
        return omega;
    }

    public double getMu_() {
        return mu_;
    }

    public double getEpsilon_() {
        return epsilon_;
    }

    public double getC_() {
        return c_;
    }

    public double getEta_() {
        return eta_;
    }

    public ArrayList<RWGelement> getRwgElements() {
        return rwgElements;
    }

    public ArrayList<Triangle> getTriangles() {
        return triangles;
    }

    public ArrayList<Triangle> getTrianglesPlus() {
        return trianglesPlus;
    }

    public ArrayList<Triangle> getTrianglesMinus() {
        return trianglesMinus;
    }

    public ArrayList<Point> getTrianglesCenters() {
        return trianglesCenters;
    }

    public ArrayList<ArrayList<Point>> getTrianglesSubCenters() {
        return trianglesSubCenters;
    }

    public ArrayList<Double> getEdgesLenght() {
        return edgesLenght;
    }

    public ArrayList<Point> getRhoPlus() {
        return rhoPlus;
    }

    public ArrayList<Point> getRhoMinus() {
        return rhoMinus;
    }
}
