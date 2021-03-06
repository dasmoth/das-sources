/* -*- mode: java; c-basic-offset: 4; indent-tabs-mode: nil -*- */

package das.bam;

import java.io.File;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.util.*;

import javax.servlet.ServletContext;
import javax.sql.DataSource;

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;

import org.biojava.bio.Annotation;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SimpleFeatureHolder;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.DummySymbolList;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.servlets.dazzle.datasource.*;
import org.biojava.utils.JDBCPooledDataSource;
import org.biojava.utils.SmallSet;

import utils.Collects;

/**
 * DAS source backed by an (indexed) BAM file.
 *
 * @author Thomas Down
 */
public class BAMMappingFeatureSource extends AbstractDataSource implements TilingFeatureSource {
    private String bamPath;
    private String bamIndexPath;
    private int minTile = 100;
    private int defaultMaxBins = 500;
    private int qualityThreshold = -1;
    private boolean groupPairs = false;
    private Set<String> seqNames;

    private SAMFileReader db;
	
    public void setGroupPairs(boolean b) {
	this.groupPairs = b;
    }

    public void setQualityThreshold(int i) {
	this.qualityThreshold = i;
    }
    
    public void setMinTile(int i) {
	this.minTile = i;
    }
	
    public void setBamPath(String s) {
	    this.bamPath = s;
    }
	
    public void setBamIndexPath(String s) {
	this.bamIndexPath = s;
    }
	

    public void init(ServletContext context)
    	throws DataSourceException
    {
	super.init(context);
    	try {
	    if (bamIndexPath == null) {
		bamIndexPath = bamPath + ".bai";
	    }
	    db = new SAMFileReader(new File(bamPath), new File(bamIndexPath));
	    db.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);

            seqNames = new HashSet<String>();
            for (SAMSequenceRecord ssr : db.getFileHeader().getSequenceDictionary().getSequences()) {
                seqNames.add(ssr.getSequenceName());
            }
    	} catch (Exception ex) {
	    throw new DataSourceException(ex);
    	}
    }
	
	@Override
	public Sequence getSequence(String ref) throws DataSourceException, NoSuchElementException {
		return makeSeq(ref);
	}
	
	public FeatureHolder getFeatures(String ref, int maxbins) {
		return makeSeq(ref, maxbins);
	}

	public Set getAllTypes() {
		Set<String> s = new SmallSet();
		s.add("mapping");
		s.add("density");
		return s;
	}

	public String getDataSourceType() {
		return "bam";
	}

	public String getDataSourceVersion() {
		return "0.1.0";
	}

	public String getLandmarkVersion(String ref) throws DataSourceException, NoSuchElementException {
		return "";
	}

	public String getMapMaster() {
		return "";
	}
	
	public String getScore(Feature f) {
		Annotation fa = f.getAnnotation();
		if (fa.containsProperty("score")) {
			return fa.getProperty("score").toString();
		} else {
			return super.getScore(f);
		}
	} 
	
	
	private double mean(Collection<? extends Number> l) {
		double tot = 0;
		for (Number n : l) {
			tot += n.doubleValue();
		}
		return tot / l.size();
	}
	
    private Seq makeSeq(String name, int maxbins) {
        if (!seqNames.contains(name)) {
            if (name.startsWith("chr")) {
                name = name.substring(3);
            } else {
                name = "chr" + name;
            }
        }
        if (!seqNames.contains(name)) {
            throw new RuntimeException("No sequence " + name);
        }
        return new Seq(name, maxbins);
    }
    
    private Seq makeSeq(String name) {
        return new Seq(name, -1);
    }

	private class Seq extends SimpleSequence {
		private int maxbins = -1;
		
            public Seq(String name, int maxbins) {
			super(
					new DummySymbolList(DNATools.getDNA(), Integer.MAX_VALUE), 
					name, 
					name, 
					Annotation.EMPTY_ANNOTATION
			);
                        this.maxbins = maxbins;
            }
		
		public FeatureHolder filter(FeatureFilter ff) {
			try {
				synchronized (db) {
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
					double[] tileCounts = new double[maxTile - minTile + 1];
					
					StrandedFeature.Template templ = new StrandedFeature.Template();
					templ.source = "sam";
					templ.type = "mapping";
					templ.annotation = new SmallAnnotation();
					SimpleFeatureHolder result = new SimpleFeatureHolder();
					
					CloseableIterator<SAMRecord> i = db.query(getName() , loc.getMin(), loc.getMax(), false);
					try {
						for (; i.hasNext(); ) {
							SAMRecord r = i.next();
							
							if (r.getAlignmentStart() <= 0 || r.getAlignmentEnd() <= 0) {
								continue;
							}
							
							if (r.getMappingQuality() < qualityThreshold) {
								continue;
							}
							
							templ.location = new RangeLocation(r.getAlignmentStart(), r.getAlignmentEnd());
							templ.strand = r.getReadNegativeStrandFlag() ? StrandedFeature.NEGATIVE : StrandedFeature.POSITIVE;
							// templ.annotation.setProperty("score", new Double(score));
							if (groupPairs) {
							    templ.annotation.setProperty("pair", r.getReadName());
							}
							Feature f = this.createFeature(templ);
							if (ff.accept(f)) {
								result.addFeature(f);
							}
							
							if (!r.getReadPairedFlag() || (r.getFirstOfPairFlag() && r.getProperPairFlag())) {
							    int minPos = templ.location.getMin();
							    int maxPos = templ.location.getMax();
							    if (r.getReadPairedFlag()) {
								int as = r.getAlignmentStart();
								int ae = r.getAlignmentEnd();
								int ms = r.getMateAlignmentStart();
								if (as < ms) {
								    minPos = as;
								    maxPos = ms + (ae - as);
								} else {
								    minPos = ms;
								    maxPos = ae;
								}
							    }

							    int minReadTile = Math.max(minTile, minPos / tileSize);
							    int maxReadTile = Math.min(maxTile, maxPos / tileSize);
							    
							    for (int t = minReadTile; t <= maxReadTile; ++t) {
								int tileStart = (t) * tileSize + 1;
								int tileEnd = (t + 1) * tileSize;
								int lapStart = Math.max(tileStart, minPos);
								int lapEnd = Math.min(tileEnd, maxPos);
								
								tileCounts[t - minTile] += ((1.0 * (lapEnd - lapStart + 1)) / (maxPos - minPos + 1));
							    }
							}
						}
					} finally {
						i.close();
					}
					
					templ.source = "sam";
					templ.type = "density";
					templ.strand = StrandedFeature.UNKNOWN;
					templ.annotation = new SmallAnnotation();
					for (int t = 0; t < tileCounts.length; ++t) {
						templ.location = new RangeLocation((minTile + t) * tileSize + 1, (minTile + t + 1) * tileSize);
						templ.annotation.setProperty("score", new Double((1000.0 * tileCounts[t]) / tileSize));
						Feature f = this.createFeature(templ);
						if (ff.accept(f)) {
							result.addFeature(f);
						}
					}
					
					return result;
				}
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

    public List getGroups(Feature f) {
	Annotation a = f.getAnnotation();
	if (a.containsProperty("pair")) {
	    return Collections.singletonList(new DASGFFGroup(a.getProperty("pair").toString(), "pair"));
	} else {
	    return Collections.emptyList();
	}
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
