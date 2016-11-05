/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package utilities;

/**
 *
 * @author nageshbhattu
 */
public class test {
    public static void main(String[] args){
        String str = "1 2:10 4:5   ";
        String[] strs = str.split(" ");
        System.out.println(" Length is "+ strs.length);
        for(int i = 0;i<strs.length;i++)
            System.out.println(strs[i]);
    }
}
