/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.classifier.io;

import edu.msu.cme.rdp.classifier.ClassificationResult;
import edu.msu.cme.rdp.classifier.RankAssignment;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

/**
 *
 * @author wangqion
 */
public class ClassificationResultFormatter {
    // list of major rankd

    public static String[] RANKS = { "domain", "phylum", "class", "order", "family", "genus"};

    public enum FORMAT {

        allRank, fixRank, dbformat;
    }

    public static String getOutput(ClassificationResult result, FORMAT format) {
        switch (format) {
            case allRank:
                return getAllRankOutput(result);
            case fixRank:
                return getFixRankOutput(result);
            case dbformat:
                return getDBOutput(result);
            default:
                getAllRankOutput(result);
        }
        return null;
    }

    public static String getAllRankOutput(ClassificationResult result) {
        StringBuilder assignmentStr = new StringBuilder(result.getSequence().getSeqName() + "\t");
        if (result.isReverse()) {
            assignmentStr.append("-");
        }
        for (RankAssignment assignment : (List<RankAssignment>) result.getAssignments()) {
            assignmentStr.append("\t").append(assignment.getName()).append("\t").append(assignment.getRank()).append("\t").append(assignment.getConfidence());
        }
        assignmentStr.append("\n");
        return assignmentStr.toString();
    }

    public static String getAllRankOutput(ClassificationResult result, double conf) {
        StringBuilder assignmentStr = new StringBuilder(result.getSequence().getSeqName() + "\t");
        if (result.isReverse()) {
            assignmentStr.append("-");
        }
        for (RankAssignment assignment : (List<RankAssignment>) result.getAssignments()) {

            if (assignment.getConfidence() >= conf) {
                assignmentStr.append("\t").append(assignment.getName()).append("\t").append(assignment.getRank()).append("\t").append(assignment.getConfidence());
            }

        }
        assignmentStr.append("\n");
        return assignmentStr.toString();
    }

    public static String getFixRankOutput(ClassificationResult result) {
        return getFixRankOutput(RANKS, result);
    }

    public static String getFixRankOutput(String[] ranks, ClassificationResult result) {
        StringBuilder assignmentStr = new StringBuilder();

        HashMap<String, RankAssignment> rankMap = new HashMap<String, RankAssignment>();
        for (RankAssignment assignment : (List<RankAssignment>) result.getAssignments()) {
            rankMap.put(assignment.getRank(), assignment);
        }
        
        // if the score is missing for the rank, report the conf and name from the lower rank
        RankAssignment prevAssign = null;
        for (int i = ranks.length -1; i>=0; i--) {
            RankAssignment assign = rankMap.get(ranks[i]);
            if (assign != null) {
                assignmentStr.insert(0, "\t" + assign.getName() +"\t" + assign.getRank() + "\t" + assign.getConfidence());
                prevAssign = assign;
            } else {
                assignmentStr.insert(0, "\t" + prevAssign.getName() +"\t" + ranks[i] + "\t" + prevAssign.getConfidence());
            }
            
        }
        if (result.isReverse()) {
            assignmentStr.insert(0,"-");
        } else {
            assignmentStr.insert(0, "");
        }
        assignmentStr.insert(0, result.getSequence().getSeqName() + "\t");
        assignmentStr.append("\n");
        
        return assignmentStr.toString();

    }

    public static String getDBOutput(ClassificationResult result) {
        StringBuilder assignmentStr = new StringBuilder();
        Iterator<RankAssignment> it = result.getAssignments().iterator();

        while (it.hasNext()) {
            RankAssignment assign = it.next();
            assignmentStr.append(result.getSequence().getSeqName()).append("\t").append(result.getTrainsetNo()).append("\t").append(assign.getTaxid()).append("\t").append(assign.getConfidence()).append("\n");
        }

        return assignmentStr.toString();
    }
}
