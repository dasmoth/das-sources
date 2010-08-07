package das.jkdb;

import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FeatureRealizer;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SimpleFeatureHolder;
import org.biojava.bio.seq.SimpleFeatureRealizer;
import org.biojava.bio.seq.Feature.Template;
import org.biojava.bio.seq.impl.FeatureImpl;
import org.biojava.bio.symbol.AbstractSymbolList;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.AssertionFailure;
import org.biojava.utils.ChangeVetoException;

class JKSequence extends AbstractSymbolList implements Sequence {
	private final static Symbol[] LUT = new Symbol[] {DNATools.t(), DNATools.c(), DNATools.a(), DNATools.g()};
	
	private final String name;
	private final ByteBuffer buffer;
	private final int length;
	private Location nLocation;
	private Location maskLocation;
	private final int seqStart;
	private int dnaStart;
	private int longestMaskBlock = 0;
	private boolean elideRepeats;
	
	JKSequence(String name, ByteBuffer buffer, int start, boolean elideRepeats) {
		this.name = name;
		this.buffer = buffer;
		this.seqStart = start;
		this.elideRepeats = elideRepeats;
		
		synchronized (buffer) {
			buffer.position(seqStart);
			length = buffer.getInt();
		}
	}
	
	private void init() {
		synchronized (buffer) {
			buffer.position(seqStart + 4);
			{
				int nBlockCnt = buffer.getInt();
				int[] nBlockStarts = new int[nBlockCnt];
				for (int i = 0; i < nBlockCnt; ++i) {
					nBlockStarts[i] = buffer.getInt();
				}
				int[] nBlockSizes = new int[nBlockCnt];
				for (int i = 0; i < nBlockCnt; ++i) {
					nBlockSizes[i] = buffer.getInt();
				}
				List<Location> nBlocks = new ArrayList<Location>();
				for (int i = 0; i < nBlockCnt; ++i) {
					nBlocks.add(new RangeLocation(nBlockStarts[i], nBlockStarts[i] + nBlockSizes[i] - 1));
				}
				nLocation = LocationTools.union(nBlocks);
			}
			int mBlockCnt = buffer.getInt();
			if (!elideRepeats) {
				int[] mBlockStarts = new int[mBlockCnt];
				for (int i = 0; i < mBlockCnt; ++i) {
					mBlockStarts[i] = buffer.getInt();
				}
				int[] mBlockSizes = new int[mBlockCnt];
				for (int i = 0; i < mBlockCnt; ++i) {
					mBlockSizes[i] = buffer.getInt();
					longestMaskBlock = Math.max(longestMaskBlock, mBlockSizes[i]);
				}
				List<Location> mBlocks = new ArrayList<Location>();
				for (int i = 0; i < mBlockCnt; ++i) {
					mBlocks.add(new RangeLocation(mBlockStarts[i], mBlockStarts[i] + mBlockSizes[i] - 1));
				}
				maskLocation = LocationTools.union(mBlocks);
			} else {
				buffer.position(buffer.position() + (mBlockCnt*8));
				maskLocation = Location.empty;
			}
			buffer.getInt();
			dnaStart = buffer.position();
		}
	}
	
	public String getName() {
		return name;
	}

	public String getURN() {
		return getName();
	}

	public Alphabet getAlphabet() {
		return DNATools.getDNA();
	}

	public int length() {
		return length;
	}

	public Symbol symbolAt(int index) throws IndexOutOfBoundsException 
	{
		if (nLocation == null) {
			init();
		}
		
		if (index < 1 || index > length) {
			throw new IndexOutOfBoundsException(String.format("%d is outside 1:%d", index, length));
		}
		
		index -= 1;
		if (nLocation.contains(index)) {
			return DNATools.n();
		}
		
		int major = index >> 2;
		int minor = index & 0x3;
		
		int b = buffer.get(dnaStart + major);
		switch (minor) {
		case 0:
			return LUT[b>>6 & 0x3];
		case 1:
			return LUT[b>>4 & 0x3];
		case 2:
			return LUT[b>>2 & 0x3];
		case 3:
			return LUT[b & 0x3];
		}
		throw new AssertionFailure("Fall-though in 2bit unpacker");
	}

	public boolean containsFeature(Feature f) {
		return f.getSequence() == this;
	}

	public int countFeatures() {
		return filter(FeatureFilter.all).countFeatures(); // EEP this is inefficient!
	}

	public Feature createFeature(Template ft) throws BioException, ChangeVetoException {
		throw new ChangeVetoException();
	}

	public Iterator features() {
		return filter(FeatureFilter.all).features();
	}

	public FeatureHolder filter(FeatureFilter fc, boolean recurse) {
		return filter(fc);
	}

	public FeatureHolder filter(FeatureFilter filter) {
		if (maskLocation == null) {
			init();
		}
		
		Location ol = FilterUtils.extractOverlappingLocation(filter);
		Location l;
		if (ol == null) {
			l = maskLocation;
		} else {
			l = LocationTools.intersection(maskLocation, new RangeLocation(ol.getMin() - longestMaskBlock, ol.getMax() + longestMaskBlock));
		}
		SimpleFeatureHolder fh = new SimpleFeatureHolder();
		Feature.Template temp = new Feature.Template();
		temp.type = "repeat";
		temp.source = "ucsc";
		temp.annotation = new SmallAnnotation();
		for (Iterator<?> i = l.blockIterator(); i.hasNext(); ) {
			Location bloc = (Location) i.next();
			temp.location = bloc;
			try {
				Feature f = FeatureImpl.DEFAULT.realizeFeature(this, this, temp);
				if (filter.accept(f)) {
					fh.addFeature(f);
				}
			} catch (Exception ex) {
				throw new BioError(ex);
			}
		}
		return fh;
	}

	public FeatureFilter getSchema() {
		return new FeatureFilter.ByType("repeat");
	}

	public void removeFeature(Feature f) throws ChangeVetoException, BioException {
		throw new ChangeVetoException();
	}

	public Annotation getAnnotation() {
		return Annotation.EMPTY_ANNOTATION;
	}
}
