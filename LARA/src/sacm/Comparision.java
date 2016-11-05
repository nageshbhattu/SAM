/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sacm;

import aspectSegmenter.Reviews2Vectors;
import java.util.Arrays;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

/**
 *
 * @author agnibesh
 */
public class Comparision {
    
   // SACModel sacm = new SACModel();
    PearsonsCorrelation pc = new PearsonsCorrelation();
    int numAspects = Reviews2Vectors.numAspects;
    int numDocs = Reviews2Vectors.dlist.size();
    public double final_mse;
    
    public void mse(double[][] x,double[][] y){
        double[][] predicted = x;
        double[][] groundtruth = y;
        double sum = 0.0;
        double mse_sum = 0.0;
        
        for(int di=0;di<numDocs;di++){
            for(int ai=0;ai<numAspects;ai++){
            double numerator = (predicted[di][ai]-groundtruth[di][ai]*(predicted[di][ai]-groundtruth[di][ai]));
            mse_sum += numerator;
           }   
        }
        final_mse = mse_sum/(numAspects*numDocs);
        System.out.println(final_mse);
    }
    
    public void p_aspect(double[][] x, double[][] y){
        double[][] predicted = x;
        double[][] groundTruth = y;
        double psc = 0.0;
        for(int di=0; di<numDocs; di++){
        double ps = pc.correlation(predicted[di],groundTruth[di]);
        psc += ps;
        }
        double p_aspect_final = psc/numDocs;
        System.out.println(p_aspect_final);
    }
    
    public void p_review(double[][] x, double[][] y){
        double[][] predicted = x;
        double[][] groundTruth = y;
        double p_re = 0.0;
        for(int ai=0; ai<numAspects; ai++){
        double p_r = pc.correlation(predicted[ai],groundTruth[ai]);
        p_re += p_r;
        }
        double p_rev_final = p_re/numAspects;
        System.out.println(p_rev_final);
    }
    
    public void map(double[][] x, double[][] y){
        double[][] predicted = x;
        double[][] groundTruth = y;
        int rel = 0;
        double map = 0.0;
        double AveP;
        double[] p_rev = new double[numAspects];
        int Q = Math.floorDiv(numDocs, 10);
        for(int ai=0; ai<numAspects; ai++){
           //Arrays.sort(predicted[ai]);
           //Arrays.sort(groundTruth[ai]);
           bubbleSort(predicted[ai]);
           bubbleSort(groundTruth[ai]);
           for(int di=0;di<Q;di++){
               for(int aii=0;aii<numAspects;aii++){
               double match = predicted[di][ai];
                    for(int dii=0;dii<Q;dii++){
                    if(groundTruth[dii][aii]==match){
                    rel ++;
                    }
                    else{
                    rel += 0;
                    }
                    }
               }
                double pre = rel/Q;
                AveP = (pre*rel)/Q;
                map += AveP;
           }
           
        
        }
        
       
       
    }
    public void nDCG(double[][] x, double[][] y){
        double[][] predicted = x;
        double[][] groundTruth = y;
        double p_sum = 0.0;
        double g_sum = 0.0;
        double DCG = 0;
        double IDCG = 0;
        
       for(int di=0;di<numDocs;di++) {
        for(int ai=1;ai<predicted.length;ai++){
            p_sum += predicted[di][ai]/(Math.log(ai)/Math.log(2));
            g_sum += groundTruth[di][ai]/(Math.log(ai)/Math.log(2));
        }
         DCG += predicted[di][0]+p_sum;
         IDCG += groundTruth[di][0]+g_sum;
        
       }
       double nDCG = DCG/IDCG;
       System.out.println(nDCG);
    }
      private static void bubbleSort(double[] intArray) {
               int n = intArray.length;
                double temp = 0;
               
                for(int i=0; i < n; i++){
                        for(int j=1; j < (n-i); j++){
                               
                                if(intArray[j-1] < intArray[j]){
                                        //swap the elements!
                                        temp = intArray[j-1];
                                        intArray[j-1] = intArray[j];
                                        intArray[j] = temp;
                                }
                        }
                }
        }
}
