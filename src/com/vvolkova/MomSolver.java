package com.vvolkova;

import com.vvolkova.geometric.Point;
import com.vvolkova.impedance.Complex;
import com.vvolkova.impedance.ImpedanceMatrix;
import com.vvolkova.rwg.RWGelement;

import java.util.ArrayList;

public class MomSolver {
    private Complex f;
    private Complex epsilon_;
    private Complex mu_;
    private Complex c_ ;
    private Complex eta_;
    private Complex pi;
    private Complex omega;
    private ArrayList<RWGelement> rwgElements;
    private ArrayList<Point> trianglesCenters;
    private ArrayList<Point> rhoPlus;
    private ArrayList<Point> rhoMinus;
    private ArrayList<Double> edgesLenght;

    private Complex[] d;
    private Complex[] pol;
    private Complex k;
    private Complex kv[];

    private ArrayList<ArrayList<Complex>> Z;
    private ArrayList<Complex> V;
    private ArrayList<Complex> I;

    public MomSolver(ImpedanceMatrix m) {
        this.f = new Complex( m.getF(),0);
        this.epsilon_ = new Complex(m.getEpsilon_(),0);
        this.mu_ = new Complex(m.getMu_(), 0);
        this.c_ = new Complex(m.getC_(), 0);
        this.eta_ = new Complex(m.getEta_(), 0);
        this.omega = new Complex(m.getOmega(), 0);
        this.rwgElements = m.getRwgElements();
        this.trianglesCenters = m.getTrianglesCenters();
        this.rhoPlus = m.getRhoPlus();
        this.rhoMinus = m.getRhoMinus();
        this.edgesLenght = m.getEdgesLenght();

        //this.Z = m.getZ();!!!!!!!

        initConstants();
        calculateV();
        foundSystemSolution();
        checkI();
    }

    private void initConstants(){
        d = new Complex[] {new Complex(0, 0), new Complex(0, 0), new Complex(-1, 0)};
        pol = new Complex[] {new Complex(1, 0), new Complex(0, 0), new Complex(0, 0)};
        k = omega.divides(c_);
        kv = new Complex[3];
        V = new ArrayList<Complex>();
        I = new ArrayList<Complex>();

        for (int i = 0; i < 3; i++) {
            kv[i] = k.times(d[i]);
        }

        for (int i = 0; i < Z.size(); i++) {
            I.add(new Complex(0,0));
        }
    }

    private void calculateV() {
        for (int m = 0; m < rwgElements.size(); m++) {
            RWGelement rwGelement = rwgElements.get(m);
            Complex scalarProduct = new Complex(0,0);
            scalarProduct = scalarProduct.plus(kv[0].times(trianglesCenters.get(rwGelement.getPlusTriangle()).getX()));
            scalarProduct = scalarProduct.plus(kv[1].times(trianglesCenters.get(rwGelement.getPlusTriangle()).getY()));
            scalarProduct = scalarProduct.plus(kv[2].times(trianglesCenters.get(rwGelement.getPlusTriangle()).getZ()));

            Complex[] emPlus = new Complex[3];
            emPlus[0] = scalarProduct.times(new Complex(0,-1)).exp().times(pol[0]);
            emPlus[1] = scalarProduct.times(new Complex(0,-1)).exp().times(pol[1]);
            emPlus[2] = scalarProduct.times(new Complex(0,-1)).exp().times(pol[2]);

            scalarProduct = new Complex(0,0);
            scalarProduct = scalarProduct.plus(kv[0].times(trianglesCenters.get(rwGelement.getMinusTriangle()).getX()));
            scalarProduct = scalarProduct.plus(kv[1].times(trianglesCenters.get(rwGelement.getMinusTriangle()).getY()));
            scalarProduct = scalarProduct.plus(kv[2].times(trianglesCenters.get(rwGelement.getMinusTriangle()).getZ()));

            Complex[] emMinus = new Complex[3];
            emMinus[0] = scalarProduct.times(new Complex(0,-1)).exp().times(pol[0]);
            emMinus[1] = scalarProduct.times(new Complex(0,-1)).exp().times(pol[1]);
            emMinus[2] = scalarProduct.times(new Complex(0,-1)).exp().times(pol[2]);


            Point rp = rhoPlus.get(m);
            Complex scalarPlus = new Complex(0,0);
            scalarPlus = scalarPlus.plus(emPlus[0].times(rp.getX()));
            scalarPlus = scalarPlus.plus(emPlus[1].times(rp.getY()));
            scalarPlus = scalarPlus.plus(emPlus[2].times(rp.getZ()));
            Point rm = rhoMinus.get(m);
            Complex scalarMinus = new Complex(0,0);
            scalarMinus = scalarMinus.plus(emMinus[0].times(rm.getX()));
            scalarMinus = scalarMinus.plus(emMinus[1].times(rm.getY()));
            scalarMinus = scalarMinus.plus(emMinus[2].times(rm.getZ()));

            Complex v = (scalarPlus.divides(new Complex(2,0)).plus(scalarMinus.divides(new Complex(2,0)))).times(edgesLenght.get(m));
            V.add(v);
        }
    }

    /*private void foundSystemSolution() {
        //СЛАУ
        // Считываем размер вводимой матрицы
        int size = Z.size();
        Complex[][] matrix = new Complex[size][size+1];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size+1; j++) {
                matrix[i][j] = new Complex(0,0);
            }
        }


        // Матрица будет иметь размер (size) x (size + 1),
        // c учетом столбца свободных членов
        for (int i = 0; i < size; i++) {
            ArrayList<Complex> z = Z.get(i);
            for (int j = 0; j < size; j++) {
                matrix[i][j] = z.get(j).clone();
            }
            matrix[i][size] = V.get(i).clone();
        }

        // Считываем необходимую точность решения
        double eps = 0.1;

        // Введем вектор значений неизвестных на предыдущей итерации,
        // размер которого равен числу строк в матрице, т.е. size,
        // причем согласно методу изначально заполняем его нулями
        Complex[] previousVariableValues = new Complex[size];
        for (int i = 0; i < size; i++) {
            previousVariableValues[i] = new Complex(0,0);
        }

        // Будем выполнять итерационный процесс до тех пор,
        // пока не будет достигнута необходимая точность
        while (true) {
            // Введем вектор значений неизвестных на текущем шаге
            Complex[] currentVariableValues = new Complex[size];
            for (int i = 0; i < size; i++) {
                currentVariableValues[i] = new Complex(0,0);
            }

            // Посчитаем значения неизвестных на текущей итерации
            // в соответствии с теоретическими формулами
            for (int i = 0; i < size; i++) {
                // Инициализируем i-ую неизвестную значением
                // свободного члена i-ой строки матрицы
                currentVariableValues[i] = matrix[i][size].clone();

                // Вычитаем сумму по всем отличным от i-ой неизвестным
                for (int j = 0; j < size; j++) {
                    // При j < i можем использовать уже посчитанные
                    // на этой итерации значения неизвестных
                    if (j < i) {
                        //currentVariableValues[i] -= matrix[i][j] * currentVariableValues[j];
                        currentVariableValues[i] = currentVariableValues[i].minus(matrix[i][j].times(currentVariableValues[j]));
                    }

                    // При j > i используем значения с прошлой итерации
                    if (j > i) {
                        //currentVariableValues[i] -= matrix[i][j] * previousVariableValues[j];
                        currentVariableValues[i] = currentVariableValues[i].minus(matrix[i][j].times(previousVariableValues[j]));
                    }
                }

                // Делим на коэффициент при i-ой неизвестной
                //currentVariableValues[i] /= matrix[i][i];
                currentVariableValues[i] = currentVariableValues[i].divides(matrix[i][i]);
            }

            // Посчитаем текущую погрешность относительно предыдущей итерации
            double error = 0.0;

            for (int i = 0; i < size; i++) {
                //error += Math.abs(currentVariableValues[i] - previousVariableValues[i]);
                error += (currentVariableValues[i].minus(previousVariableValues[i])).abs();
            }

            // Если необходимая точность достигнута, то завершаем процесс
            if (error < eps) {
                break;
            }

            // Переходим к следующей итерации, так
            // что текущие значения неизвестных
            // становятся значениями на предыдущей итерации

            for (int i = 0; i < previousVariableValues.length; i++) {
                previousVariableValues[i] = currentVariableValues[i].clone();
            }
        }

        for (int i = 0; i < previousVariableValues.length; i++) {
            I.add(previousVariableValues[i]);
        }
    }*/

    private void foundSystemSolution(){
        //СЛАУ метод Гаусса
        int size = Z.size();
        Complex[][] a = new Complex[size][size+1];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size+1; j++) {
                a[i][j] = new Complex(0,0);
            }
        }
        for (int i = 0; i < size; i++) {
            ArrayList<Complex> z = Z.get(i);
            for (int j = 0; j < size; j++) {
                a[i][j] = z.get(j).clone();
            }
            a[i][size] = V.get(i).clone();
        }
        Complex x[]=new Complex[a.length];
        for (int i = 0; i < x.length; i++) {
            x[i] = a[i][a[i].length - 1];
        }
        Complex m;
        for (int k = 1; k < a.length; k++) {
            for (int j = k; j < a.length; j++) {
                m = a[j][k - 1].divides( a[k - 1][k - 1]);
                for (int i = 0; i < a[j].length; i++) {
                    a[j][i] = a[j][i].minus( m.times( a[k - 1][i]));
                }
                x[j] = x[j].minus( m.times(x[k - 1]));
            }
        }
        for (int i=a.length-1;i>=0;i--) {
            for (int j=i+1;j<a.length;j++) {
                x[i]= x[i].minus(a[i][j].times(x[j]));
            }
            x[i] = x[i].divides(a[i][i]);
        }
        for (int i = 0; i < x.length; i++) {
            I.set(i,x[i]);
        }
    }

    private void checkI() {
        for (int i = 0; i < I.size(); i++) {
            System.out.println((i+1) + ") " + I.get(i));
        }
    }

}
