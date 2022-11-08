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
 * Import organism data from an NCBI taxonomy dump.
 * Taxonomy names file (names.dmp):
 *	tax_id					-- the id of node associated with this name
 *	name_txt				-- name itself
 *	unique name				-- the unique variant of this name if name not unique
 *	name class				-- (synonym, common name, ...)
 * Field terminator is "\t|\t"
 * Row terminator is "\t|\n"
 * <pre>
 * 3827	| chickpea	     |   | genbank common name |
 * 3827	| Cicer arietinum L. |   | authority |
 * 3827	| Cicer arietinum    |   | scientific name |
 * 3827	| garbanzo	     |   | common name |
 * </pre>
 *
 * @author Sam Hokin
 */
public class NCBITaxonomyFileConverter extends BioFileConverter {

    private static final Logger LOG = Logger.getLogger(NCBITaxonomyFileConverter.class);

    // could have multiple strains per organism; only store organism once
    Map<String,Item> organismMap = new HashMap<>();

    // the list of organisms to load
    List<String> taxonIds = new ArrayList<>();
    
    /**
     * Create a new NCBITaxonomyFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public NCBITaxonomyFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * Set the comma-separated organism taxon IDs in project.xml, e.g. "3827,3847,3885"
     */
    public void setTaxonIds(String taxonIdString) {
        String[] ids = taxonIdString.split(",");
        for (String id : ids) {
            taxonIds.add(id);
        }
    }

    /**
     * {@inheritDoc}
     * Read in the names.dmp file and store the organism info.
     */
    @Override
    public void process(Reader reader) throws Exception {

        // only process names.dmp
        if (!getCurrentFile().getName().equals("names.dmp")) return;

        LOG.info("Processing NCBI taxonomy dump "+getCurrentFile().getName()+"...");

        // tax_id	-- the id of node associated with this name
        // name_txt	-- name itself
        // unique name	-- the unique variant of this name if name not unique
        // name class	-- (synonym, common name, ...)
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            line = line.replace("\t", "");        // get rid of the tab parts
            String[] parts = line.split("\\|");
            String taxId = parts[0];
            String nameTxt = parts[1];
            String uniqueName = parts[2];
            String nameClass = parts[3];
            if (taxonIds.contains(taxId)) {
                // DEBUG
                System.out.println(taxId+"|"+nameTxt+"|"+uniqueName+"|"+nameClass);
                //
                Item organism = organismMap.get(taxId);
                if (organism==null) {
                    organism = createItem("Organism");
                    organism.setAttribute("taxonId", taxId);
                    organismMap.put(taxId, organism);
                }
                if (nameClass.equals("genbank common name")) {
                    organism.setAttribute("commonName", nameTxt);
                } else if (nameClass.equals("scientific name")) {
                    organism.setAttribute("name", nameTxt);
                    String[] genusSpecies = nameTxt.split(" ");
                    organism.setAttribute("genus", genusSpecies[0]);
                    organism.setAttribute("species", genusSpecies[1]);
                    organism.setAttribute("shortName", genusSpecies[0].charAt(0)+". "+genusSpecies[1]);
                }
            }
        }

        // wrap up this file
        br.close();

        // store the items
        store(organismMap.values());
    }
}
