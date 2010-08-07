package io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;

import utils.Collects;
import utils.Function;

public class IOTools {
	private IOTools() {
	}
	
	public static BufferedReader fileBufferedReader(File f) 
		throws Exception
	{
		return new BufferedReader(fileReader(f));
	}
	
	public static Reader fileReader(File f)
		throws Exception
	{
		if (f.getName().endsWith(".gz")) {
			return new InputStreamReader(new GZIPInputStream(new FileInputStream(f)));
		} else {
			return new FileReader(f);
		}
	}

	public static Reader nameReader(String arg)
		throws Exception
	{
		if ("-".equals(arg)) {
			return new InputStreamReader(System.in);
		} else {
			return fileReader(new File(arg));
		}
	}
	
	public static BufferedReader nameBufferedReader(String arg) 
		throws Exception
	{
		return new BufferedReader(nameReader(arg));
	}
	
	public static Reader inputReader(String[] args) 
		throws Exception
	{
		if (args.length == 0) {
			return nameReader("-");
		} else if (args.length == 1) {
			return nameReader(args[0]);
		} else {
			return new SequenceReader(
					Collects.map(
							new Function<String,Reader>() {
									public Reader apply(String name) {
										try {
											return nameReader(name);
										} catch (Exception ex) {
											throw new RuntimeException(ex);
										}
									}
								
							},
							Arrays.asList(args)
					)
			);
		}
	}
	
	public static BufferedReader inputBufferedReader(String[] args) 
		throws Exception
	{
		return new BufferedReader(inputReader(args));
	}
}
