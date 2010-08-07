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

import javax.servlet.ServletContext;


import org.biojava.bio.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.*;
import org.biojava.servlets.dazzle.datasource.AbstractDataSource;
import org.biojava.servlets.dazzle.datasource.DataSourceException;
import org.biojava.servlets.dazzle.datasource.DazzleReferenceSource;
import org.biojava.utils.xml.*;

/**
 * Simple example datasource backed by an EMBL file.
 *
 * @author Thomas Down
 * @version 1.00
 */

public class JKSequenceSource extends AbstractDataSource implements DazzleReferenceSource {
    private String fileName;
    private Map<String,Sequence> seqs = new HashMap<String, Sequence>();
    private SequenceDB db;
    
    public String getDataSourceType() {
        return "2bit";
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
        	db = new JKSequenceDB(new File(fileName));
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

    public Set getAllTypes() {
        return Collections.emptySet();
    }
    
    public Set getEntryPoints() {
        return db.ids();
    }
}
