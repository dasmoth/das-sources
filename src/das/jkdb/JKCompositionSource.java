/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */

package das.jkdb;

import java.util.*;
import java.io.*;

import javax.servlet.*;

import org.biojava.bio.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.servlets.dazzle.datasource.AbstractDataSource;
import org.biojava.servlets.dazzle.datasource.DataSourceException;
import org.biojava.servlets.dazzle.datasource.DazzleReferenceSource;
import org.biojava.servlets.dazzle.datasource.TilingFeatureSource;
import org.biojava.utils.SmallSet;

/**
 * Simple example datasource backed by an EMBL file.
 *
 * @author Thomas Down
 * @version 1.00
 */

public class JKCompositionSource extends AbstractDataSource implements TilingFeatureSource, DazzleReferenceSource {
    private String fileName;
    private Map<String,Sequence> seqs = new HashMap<String, Sequence>();
    private SequenceDB db;
    
	private int minTile = 10;
	private int minWindow = 500;
	private int defaultMaxBins = 500;
    
    public String getDataSourceType() {
        return "2bit-comp";
    }
    
    public String getDataSourceVersion() {
        return "1.00";
    }

    public void setFileName(String s) {
        fileName = s;
    }

    public String getMapMaster() {
        return null;
    }

    public void init(ServletContext ctx) 
        throws DataSourceException
    {
        super.init(ctx);
        try {
        	db = new JKSequenceDB(new File(fileName), true);  // elide repeats.
        } catch (Exception ex) {
            throw new DataSourceException(ex, "Couldn't load sequence file");
        }
    }


    public String getLandmarkVersion(String ref)
        throws DataSourceException, NoSuchElementException
    {
        return getVersion();
    }

    public Sequence getSequence(String ref)
        throws NoSuchElementException, DataSourceException
    {
        Sequence seq = (Sequence) seqs.get(ref);
        if (seq == null) {
        	if (db.ids().contains(ref)) {
        		try {
        			seq = db.getSequence(ref);
        		} catch (Exception ex) {
        			throw new DataSourceException(ex);
        		}
        	} else if (!ref.startsWith("chr") && db.ids().contains("chr" + ref)) {
        		try {
        			seq = db.getSequence("chr" + ref);
        		} catch (Exception ex) {
        			throw new DataSourceException(ex);
        		}
        	} else {
        		throw new NoSuchElementException("No sequence " + ref);
        	}
        	seqs.put(ref, seq);
        }
        return seq;
    }

    private Set ids;
    
    public Set getEntryPoints() {
    	if (ids == null) {
    		Set _ids = new HashSet();
    		for (Iterator<?> i = db.ids().iterator(); i.hasNext(); ) {
    			String id = i.next().toString();
    			if (id.startsWith("chr")) {
    				id = id.substring(3);
    			}
    			
    			if (id.endsWith("_random")) {
    				continue;
    			}
    			
    			if (id.equals("M")) {
    				continue;
    			}
    			
    			_ids.add(id);
    		}
    		ids = Collections.unmodifiableSet(_ids);
    	}
        return ids;
    }

    public FeatureHolder getFeatures(String ref) throws DataSourceException, NoSuchElementException {
    	return new Seq(getSequence(ref));
    }
    
	public FeatureHolder getFeatures(String ref, int maxbins) throws DataSourceException, NoSuchElementException {
		return new Seq(getSequence(ref), maxbins);
	}
	
	private class Seq extends SimpleSequence {
		private int maxbins = -1;
		
		public Seq(Sequence p) {
			super(
					p, 
					p.getName(), 
					p.getURN(), 
					Annotation.EMPTY_ANNOTATION
			);
		}
		
		public Seq(Sequence p, int maxbins) {
			this(p);
			this.maxbins = maxbins;
		}
		
		public FeatureHolder filter(FeatureFilter ff) {
			try {
				Location loc = extractShadowOverlappingLocation(ff);
				if (loc == null) {
					loc = new RangeLocation(1, length());
				}
				if (maxbins < 0) {
					maxbins = defaultMaxBins;
				}
				
				int tileSize = Math.max(minTile, (loc.getMax() - loc.getMin() + 1) / maxbins);
				int minTile = (int) Math.floor((1.0 * loc.getMin()) / tileSize);
				int maxTile = (int) Math.ceil((1.0 * loc.getMax()) / tileSize);
				
				Feature.Template templ = new Feature.Template();
				templ.source = "comp";
				templ.type = "gc";
				templ.annotation = new SmallAnnotation();
				SimpleFeatureHolder result = new SimpleFeatureHolder();
				for (int t = minTile; t <= maxTile; ++t) {
					int tmin = (t * tileSize) + 1;
					int tmax = (t+ 1) * tileSize;
					
					int cmin = tmin, cmax = tmax;
					if (tileSize < minWindow) {
						cmin = (tmin + tmax - minWindow) / 2;
						cmax = cmin + minWindow - 1;
					}
					if (cmin < 1) {
						cmin = 1;
					}
					if (cmax > length()) {
						cmax = length();
					}
					if ((cmax - cmin + 1) < ((int) (0.4*minWindow))) {
						continue;
					}
					
					Count c = count(this, new RangeLocation(cmin, cmax));
					
					templ.location = new RangeLocation(tmin, tmax);
					
					{
						templ.type = "gc";
						templ.annotation.setProperty("score", new Double((1.0 * c.gc) / (cmax - cmin + 1)));
						Feature f = this.createFeature(templ);
						if (ff.accept(f)) {
							result.addFeature(f);
						}
					}
					{
						templ.type = "cpgoe";
						templ.annotation.setProperty("score", new Double(c.oe));
						Feature f = this.createFeature(templ);
						if (ff.accept(f)) {
							result.addFeature(f);
						}
					}
				}
				
				return result;
			} catch (Exception ex) {
				throw new RuntimeException(ex);
			}
		}
		
		public FeatureHolder filter(FeatureFilter ff, boolean rec) {
			return filter(ff);
		}
		
		public Iterator features() {
			return filter(FeatureFilter.all).features();
		}
		
		public int countFeatures() {
			return filter(FeatureFilter.all).countFeatures();
		}
	}
	
	private Count count(SymbolList seq, Location block) {
		SymbolList sl = seq.subList(block.getMin(), block.getMax());
		int gCount = 0, cCount = 0;
		for (int p = 1; p <= sl.length(); ++p) {
			Symbol s = sl.symbolAt(p);
			if (s == DNATools.c()) {
				++cCount;
			} else if (s == DNATools.g()) {
				++gCount;
			}
		}
		int gcCount = gCount + cCount;
		int cpgCount = 0;
		for (int p = 1; p < sl.length(); ++p) {
			if (sl.symbolAt(p) == DNATools.c() && sl.symbolAt(p + 1) == DNATools.g()) {
				++cpgCount;
			}
		}
		double cpgRat = (1.0 * cpgCount * sl.length()) / (cCount * gCount);
		
		return new Count(gcCount, cpgCount, cpgRat);
		
	}
	
	private static class Count {
		final int gc;
		final int cpg;
		final double oe;
		
		public Count(int gc, int cpg, double oe) {
			this.gc = gc;
			this.cpg = cpg;
			this.oe = oe;
		}
	}
	
	public String getScore(Feature f) {
		return f.getAnnotation().getProperty("score").toString();
	}
	
	public Set getAllTypes() {
		Set<String> s = new SmallSet();
		s.add("cpgoe");
		s.add("gc");
		return s;
	}
	
	   public static Location extractShadowOverlappingLocation(FeatureFilter ff) {
	    	if (ff instanceof FeatureFilter.OverlapsLocation) {
	    	    return ((FeatureFilter.OverlapsLocation) ff).getLocation();
	    	} else if (ff instanceof FeatureFilter.ContainedByLocation) {
	    	    return ((FeatureFilter.ContainedByLocation) ff).getLocation();
	    	} if (ff instanceof FeatureFilter.ShadowOverlapsLocation) {
	    	    return ((FeatureFilter.ShadowOverlapsLocation) ff).getLocation();
	    	} else if (ff instanceof FeatureFilter.And) {
	    	    FeatureFilter.And ffa = (FeatureFilter.And) ff;
	    	    Location l1 = extractShadowOverlappingLocation(ffa.getChild1());
	    	    Location l2 = extractShadowOverlappingLocation(ffa.getChild2());

	    	    if (l1 != null) {
	    		if (l2 != null) {
	    		    return l1.intersection(l2);
	    		} else {
	    		    return l1;
	    		}
	    	    } else {
	    		if (l2 != null) {
	    		    return l2;
	    		} else {
	    		    return null;
	    		}
	    	    }
	    	} else if (ff instanceof FeatureFilter.Or) {
	    	    FeatureFilter.Or ffo = (FeatureFilter.Or) ff;
	    	    Location l1 = extractShadowOverlappingLocation(ffo.getChild1());
	    	    Location l2 = extractShadowOverlappingLocation(ffo.getChild2());

	    	    if (l1 != null && l2 != null) {
	    		return LocationTools.union(l1, l2);
	    	    }
	    	}

	    	return null;
	    }
}
