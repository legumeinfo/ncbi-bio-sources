package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2016 NCGR
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.BufferedReader;
import java.io.Reader;

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Import chromosome and supercontig mappings from an NCBI assembly report.
 *
 * #Sequence-Name Sequence-Role      Assigned-Molecule Assigned-Molecule-Location/Type GenBank-Accn    Relationship     RefSeq-Accn     Assembly-Unit    Sequence-Length UCSC-name
 * Ca1            assembled-molecule Ca1               Chromosome                      CM001764.1      =                NC_021160.1     Primary Assembly 48359943        na
 * scaffold1      unplaced-scaffold  na                na                              KB210354.1      =                NW_004515636.1  Primary Assembly 1037    na
 *
 * @author Sam Hokin
 */
public class NCBIAssemblyReportFileConverter extends BioFileConverter {

    private static final Logger LOG = Logger.getLogger(NCBIAssemblyReportFileConverter.class);

    /**
     * Create a new NCBIAssemblyReportFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public NCBIAssemblyReportFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * Read in the assembly report and do the work.
     *
     * seqName   seqRole            molecule moleculeLocation genBankAcc   relationship  refSecAcc      assemblyUnit     length    ucscName
     * 0         1                  2        3                4            5             6              7                8         9
     * Ca1       assembled-molecule Ca1      Chromosome       CM001764.1   =             NC_021160.1    Primary Assembly 48359943  na
     * scaffold1 unplaced-scaffold  na       na               KB210354.1   =             NW_004515636.1 Primary Assembly 1037      na
     */
    @Override
    public void process(Reader reader) throws Exception {

        // only process *_assembly_report.txt
        if (!getCurrentFile().getName().endsWith("assembly_report.txt")) return;

        LOG.info("Processing NCBI assembly report "+getCurrentFile().getName()+"...");

        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) {
                // comment or blank
                continue;
            }

            String[] parts = line.split("\t");
            if (!(parts.length==10)) {
                continue; // shouldn't happen
            }
            
            String seqName = parts[0];
            String seqRole = parts[1];
            String molecule = parts[2];
            String moleculeLocation = parts[3];
            String genBankAcc = parts[4];
            String relationship = parts[5];
            String refSecAcc = parts[6];
            String assemblyUnit = parts[7];
            int length = Integer.parseInt(parts[8]);
            String ucscName = parts[9];

            if (relationship.equals("=")) {
                if (refSecAcc==null || refSecAcc.trim().length()==0) {
                    throw new RuntimeException("Assembly report line is missing refSecAcc:"+line);
                }
                if (moleculeLocation.equals("Chromosome")) {
                    Item chromosome = createItem("Chromosome");
                    chromosome.setAttribute("primaryIdentifier", refSecAcc);
                    chromosome.setAttribute("secondaryIdentifier", seqName);
                    store(chromosome);
                } else if (seqRole.contains("scaffold")) {
                    Item supercontig = createItem("Supercontig");
                    supercontig.setAttribute("primaryIdentifier", refSecAcc);
                    supercontig.setAttribute("secondaryIdentifier", seqName);
                    store(supercontig);
                }
            }
        }

        // wrap up this file
        br.close();
    }
}
