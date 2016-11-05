/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aspectSegmenter;

import cern.colt.Arrays;
import java.util.ArrayList;

/**
 *
 * @author agnibesh
 */
public class Document{
    public int[] indices;
    public int [] counts;
    public double[][] aspectRatings;
    public double[] ratings;
    public int userID;
    public int itemID;
   
    public Document(int[] wordarr,int[] countarr, double[] ratings, int userID, int itemID){
        this.userID = userID;
        this.itemID = itemID;
        this.ratings = ratings;
       
        counts = countarr;
        indices = wordarr; 
    }
    public Document(ArrayList<Integer> countarr, ArrayList<Integer> wordarr)
    {
        counts = countarr.stream().mapToInt(i -> i).toArray();
        indices = wordarr.stream().mapToInt(i -> i).toArray(); 
    }
    public Document(int[] ind, int[] cnts){
        counts = cnts;
        indices = ind;
    }
    
    
    public String toString(){
        String str = "";
        
        str+=Arrays.toString(ratings)+" "+ userID + " " + itemID+ " ";
        for(int wi = 0;wi<indices.length;wi++)
            str += indices[wi]+ ":" + counts[wi] + " ";
        return str + "\n";
    } 
    public void setRating(double[] r){
        ratings = r;
    }
    public double[] getRatings(){
        return ratings;
    }

    public double getOverallRating() {
        return ratings[0];
    }

}
