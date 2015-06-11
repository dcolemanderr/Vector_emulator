/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.classifier.train.validation;

import java.util.Comparator;
import java.util.TreeSet;

/**
 *
 * @author wangqion
 */
public class SortedROCSet extends TreeSet< SortedROCSet.PredictionCount>{
    
    public SortedROCSet(){
        super(new ResultComparator());
    }
     
        
    public static class ResultComparator implements Comparator<PredictionCount>{
        public int compare(PredictionCount lhs, PredictionCount rhs){
            if ( Double.isNaN(lhs.fpr) || Double.isNaN(lhs.se) || Double.isNaN(rhs.fpr) || Double.isNaN(rhs.se))
                return 0;
            
            if ( lhs.fpr == rhs.fpr){
                if ( lhs.se == rhs.se){
                    return 0;
                }
                return (lhs.se > rhs.se ? 1: -1);
            }
            if ( lhs.fpr > rhs.fpr){
                return 1;
            }
            return -1;
        }
    }
    
    
    public static class PredictionCount{
        double se;
        double fpr;
        
        public PredictionCount(double f, double s){
            se = s;
            fpr = f;
        }
        
        public double getSe(){
            return se;
        }
        
        public double getFPR(){
            return fpr;
        }
    }
}
