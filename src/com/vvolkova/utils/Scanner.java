package com.vvolkova.utils;

import java.io.*;
import java.util.StringTokenizer;

public class Scanner {
    BufferedReader br;
    StringTokenizer st;

    public Scanner(File f) {
        try {
            br = new BufferedReader(new FileReader(f));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public Scanner(InputStream f) {
        br = new BufferedReader(new InputStreamReader(f));
    }

    public String next() {
        while (st == null || !st.hasMoreTokens()) {
            try {
                st = new StringTokenizer(br.readLine());
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return st.nextToken();
    }

    public int nextInt() {
        return Integer.parseInt(next());
    }

    public double nextDouble() {
        return Double.parseDouble(next());
    }
}