/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sacm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/**
 *
 * @author agnibesh
 */
public class testing {
    static String tmpTxt;
    static ArrayList<String> ss = new ArrayList<>();
    static HashMap<String,Integer> snh = new HashMap<>();
    static HashMap<String,Integer> sih = new HashMap<>();
    static ArrayList<Integer> ssh = new ArrayList<>();
    static ArrayList<String> tot = new ArrayList<>();
    static int[] s;
    static String[] total;
    static String t;
    static String[] arr;
    static int info=0;
    static int warn=0;
    public static void load() throws UnsupportedEncodingException, FileNotFoundException, IOException{
        Scanner input = new Scanner(System.in);
    //BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(f), "UTF-8"));
    while ((tmpTxt = input.nextLine()) != null) {
      total = tmpTxt.split(" ");
      tot.addAll(Arrays.asList(total));
    if (tmpTxt.startsWith("INFO")) {
        
        
        Pattern p = Pattern.compile("\"([^\"]*)\"");
        Matcher m = p.matcher(tmpTxt);
        while (m.find()) {
        ss.add(m.group(1));
       
        }
                      info++; 
                      arr = new String[ss.size()];
                      for(int i=0;i<ss.size();i++){
                          arr[i]=ss.get(i);
                      }
    }
    else if (tmpTxt.startsWith("WARN")) {
        
        
        Pattern p = Pattern.compile("\"([^\"]*)\"");
        Matcher m = p.matcher(tmpTxt);
        while (m.find()) {
        ss.add(m.group(1));
       
        }
                      warn++; 
                      arr = new String[ss.size()];
                      for(int i=0;i<ss.size();i++){
                          arr[i]=ss.get(i);
                      }
    }
    }
    
    
    
    }
    static void run(){
    
    
    for(int i=0;i<arr.length;i++){
        if("sn".equals(arr[i])){
        t=arr[i+1];
    snh.putIfAbsent(t, i);
    }
        else if("si".equals(arr[i])){
        t=arr[i+1];
    sih.putIfAbsent(t, i);
    }
       /* else if("ss".equals(arr[i])){
        int t=Integer.parseInt(arr[i+1]);
        ssh.add(t);
        
        }
        */
       
    }
    
    
    System.out.println(info);
    System.out.println(warn);
    System.out.println(snh.size());
    System.out.println(sih.size());
    
    }
    
    public static void len(){
         String str="\"ss\":";
         int abc = tot.size();
         String stem="";
    for(int j=0;j<total.length;j++){
       
        if(str.equals(tot.get(j))){
        stem = tot.get(j+1).replace(",","");
        int t=Integer.parseInt(stem);
        ssh.add(t);
    
    }
       }
    for(int i=0;i<ssh.size();i++){
    s = new int[ssh.size()];
    s[i]=ssh.get(i);
    Arrays.sort(s);
    }
    System.out.println(s[s.length-1]);
    }
                      
    public static void main(String[] args) throws FileNotFoundException, IOException{
    load();
    run();
    len();
    }
    
}
