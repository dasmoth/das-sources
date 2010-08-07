package das.jkdb;

import java.io.File;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Unchangeable;


/**
 * SequenceDB backed by Jim Kent's 2bit format.
 * 
 * @author thomasdown
 */
public class JKSequenceDB extends Unchangeable implements SequenceDB {
	private static final int TWOBIT_SIGNATURE = 0x1a412743;
	
	private ByteBuffer buffer;
	private Map<String,Integer> seqOffsets = new HashMap<String, Integer>();
	private boolean elideRepeats;
	
	public JKSequenceDB(File f)
		throws Exception
	{
		this(f, false);
	}
	
	public JKSequenceDB(File f, boolean elideRepeats)
		throws Exception
	{
		this.elideRepeats = elideRepeats;
		
		long size = f.length();
		buffer = new FileInputStream(f).getChannel().map(FileChannel.MapMode.READ_ONLY, 0, (int) size);
		int sig = buffer.getInt();
		if (sig != TWOBIT_SIGNATURE) {
			buffer.order(buffer.order() == ByteOrder.LITTLE_ENDIAN ? ByteOrder.BIG_ENDIAN : ByteOrder.LITTLE_ENDIAN);
			buffer.rewind();
			sig = buffer.getInt();
		}
		if (sig != TWOBIT_SIGNATURE) {
			throw new BioException("Bad signature");
		}
		
		int version = buffer.getInt();
		if (version != 0) {
			throw new BioException(String.format("Unsupported 2bit version %d", version));
		}
		
		int seqCnt = buffer.getInt();
		buffer.getInt();  // reserved word.
		
		for (int s = 0; s < seqCnt; ++s) {
			int ns = buffer.get();
			byte[] nameBuffer = new byte[ns];
			buffer.get(nameBuffer);
			int offset = buffer.getInt();
			String name = new String(nameBuffer);
			seqOffsets.put(name, offset);
		}
	}

	public FeatureHolder filter(FeatureFilter filter) {
		return FeatureHolder.EMPTY_FEATURE_HOLDER;
	}

	public Set ids() {
		return seqOffsets.keySet();
	}

	public SequenceIterator sequenceIterator() {
		// TODO Auto-generated method stub
		return null;
	}

	public void addSequence(Sequence seq) throws IllegalIDException,
			BioException, ChangeVetoException 
	{
		throw new ChangeVetoException();
	}

	public String getName() {
		// TODO Auto-generated method stub
		return null;
	}

	public Sequence getSequence(String id) throws IllegalIDException,
			BioException 
	{
		Integer idx = seqOffsets.get(id);
		if (idx == null) {
			throw new IllegalIDException(String.format("Can't find %s", id));
		}

		return new JKSequence(id, buffer, idx.intValue(), elideRepeats);
	}

	public void removeSequence(String id) throws IllegalIDException,
			BioException, ChangeVetoException 
	{
		throw new ChangeVetoException();
	}
}
