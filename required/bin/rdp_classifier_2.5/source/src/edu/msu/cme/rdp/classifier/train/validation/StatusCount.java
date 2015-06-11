/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.classifier.train.validation;


/**
 *
 * @author wangqion
 */
public class StatusCount {
    private int numTP = 0;
    private int numFP = 0;
    private int numTN = 0;
    private int numFN = 0;

    public void incNumTP(int i){
        numTP += i;
    }

    public void incNumFP(int i){
        numFP += i;
    }

    public void incNumTN(int i){
        numTN += i;
    }

    public void incNumFN(int i){
        numFN += i;
    }

    public int getNumTP(){
        return numTP;
    }

    public int getNumTN(){
        return numTN;
    }

    public int getNumFP(){
        return numFP;
    }

    public int getNumFN(){
        return numFN;
    }

    /**
     * sensitivity = #TP / (#TP + #FN)
     */
    public double calSensitivity(){
        return (double)numTP / ( (double) (numTP +numFN));
    }

    /**
     * specificity = #TN / (#TN + #FP)
     * @return
     */
    public double calSpecificity(){
        return (double)numTN / ( (double) (numTN +numFP));
    }
    
    /**
     * false positive rate = #FP / (#TN + #FP)
     * @return
     */
    public double calFPR(){
        return (double)numFP / ( (double) (numTN +numFP));
    }
    
    /**
     * F1 score = 2#TP / (2#TP + #FP + #FN)
     * @return
     */
    public double calF1score(){
        return 2.0 * (double)numTP / ( (double) (2*numTP + +numFP + numFN));
    }
    
    
}
