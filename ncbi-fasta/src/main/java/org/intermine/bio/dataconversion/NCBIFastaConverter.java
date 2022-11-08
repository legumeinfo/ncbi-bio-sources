package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2019 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.Reader;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NoSuchElementException;

import org.intermine.dataconversion.FileConverter;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.metadata.Util;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;
import org.intermine.objectstore.query.PendingClob;

import org.biojava.nbio.core.exceptions.ParserException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.DNASequenceCreator;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.io.PlainFastaHeaderParser;
import org.biojava.nbio.core.sequence.template.Sequence;

import org.apache.commons.lang.StringUtils;
import org.apache.tools.ant.BuildException;
import org.apache.log4j.Logger;

/**
 * Extends FileConverter so we can do more things than the FastaLoaderTask accomplishes when loading NCBI FASTAs.
 *
 * @author Sam Hokin
 */

public class NCBIFastaConverter extends FileConverter {
    static final Logger LOG = Logger.getLogger(NCBIFastaConverter.class);

    String className = "org.intermine.model.bio.Chromosome";
    String classAttribute = "primaryIdentifier";
    
    String taxonId;
    String strainIdentifier;
    String assemblyVersion;   // for proteins, CDS
    String annotationVersion; // for proteins, CDS
    String dataSourceName = null;
    String dataSourceUrl = null;
    String dataSourceDescription = null;
    String dataSetTitle;
    String dataSetUrl;
    String dataSetDescription; 
    String dataSetVersion; // from file name

    ItemWriter writer;
    Model model;

    Item organism;
    Item strain;

    

    /**
     * Constructor
     * @param writer ItemWriter
     * @param model the model to create items in
     * @throws ObjectStoreException if something goes wrong
     */
    public NCBIFastaConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
        this.writer = writer;
        this.model = model;
    }

    /**
     * Set the organism taxon ID.
     */
    public void setTaxonId(String taxonId) {
        this.taxonId = taxonId;
    }

    /**
     * Set the strain identifier.
     */
    public void setStrainIdentifier(String strainIdentifier) {
        this.strainIdentifier = strainIdentifier;
    }

    /**
     * DataSource.name for any bioentities created
     * @param dataSourceName name of datasource for items created
     */
    public void setDataSourceName(String dataSourceName) {
        this.dataSourceName = dataSourceName;
    }

    /**
     * DataSource.url for any bioentities created
     * @param dataSourceUrl url of datasource for items created
     */
    public void setDataSourceUrl(String dataSourceUrl) {
        this.dataSourceUrl = dataSourceUrl;
    }

    /**
     * DataSource.description for any bioentities created
     * @param dataSourceDescription description of datasource for items created
     */
    public void setDataSourceDescription(String dataSourceDescription) {
        this.dataSourceDescription = dataSourceDescription;
    }

    /**
     * If a value is specified this title will used when a DataSet is created.
     * @param dataSetTitle the title of the DataSets of any new features
     */
    public void setDataSetTitle(String dataSetTitle) {
        this.dataSetTitle = dataSetTitle;
    }

    /**
     * If a value is specified this url will used when a DataSet is created.
     * @param dataSetUrl the title of the DataSets of any new features
     */
    public void setDataSetUrl(String dataSetUrl) {
        this.dataSetUrl = dataSetUrl;
    }

    /**
     * If a value is specified this description will used when a DataSet is created.
     * @param dataSetDescription the description of the DataSets of any new features
     */
    public void setDataSetDescription(String dataSetDescription) {
        this.dataSetDescription = dataSetDescription;
    }

    /**
     * Set the assembly version string in project.xml.
     * @param assemblyVersion the version of this assembly
     */
    public void setAssemblyVersion(String assemblyVersion) {
        this.assemblyVersion = assemblyVersion;
    }

    /**
     * Set the annotation version string in project.xml.
     * @param annotationVersion the version of this annotation
     */
    public void setAnnotationVersion(String annotationVersion) {
        this.annotationVersion = annotationVersion;
    }

    /**
     * Process a BufferedReader
     * @param reader the Reader
     * @throws java.io.IOException if an error occurs reading GFF
     * @throws ObjectStoreException if an error occurs storing items
     */
    public void process(Reader reader) throws IOException, ObjectStoreException {
        // only process FASTA
        String sequenceType = null;
        if (getCurrentFile().getName().endsWith("fna")) {
            sequenceType = "dna";
        } else if (getCurrentFile().getName().endsWith("faa")) {
            sequenceType = "protein";
        }

        // bail if we don't know what's in the FASTA
        if (sequenceType==null) return;

        // be sure we've set the required project.xml parameters
        if (taxonId==null) throw new RuntimeException("taxonId is not provided in project.xml.");
        if (strainIdentifier==null) throw new RuntimeException("strainIdentifier is not provided in project.xml.");

        System.out.println("##############################################################################################################################");
        System.out.println("Reading "+sequenceType+" sequences from: "+getCurrentFile().getName());
        System.out.println("##############################################################################################################################");
        LOG.info("NCBIFastaLoaderTask loading file "+getCurrentFile().getName());
        
        // process FASTA file
        if (sequenceType.equalsIgnoreCase("dna")) {
            LinkedHashMap<String,DNASequence> sequenceMap = FastaReaderHelper.readFastaDNASequence(getCurrentFile());
            for (String sequenceId : sequenceMap.keySet()) {
                processDNASequence(sequenceId, sequenceMap.get(sequenceId));
            }
        } else if (sequenceType.equalsIgnoreCase("protein")) {
            LinkedHashMap<String,ProteinSequence> sequenceMap = FastaReaderHelper.readFastaProteinSequence(getCurrentFile());
            for (String sequenceId : sequenceMap.keySet()) {
                processProteinSequence(sequenceId, sequenceMap.get(sequenceId));
            }
        } else {
            throw new RuntimeException("Sequence type set in project.xml is neither dna nor protein.");
        }
    }

    /**
     * Be sure to close the data loader so the last batch gets stored. only needed for tests
     * since the data loading task usually does that for hte live builds.
     * @throws ObjectStoreException if we can't store to db
     */
    public void close() throws ObjectStoreException {
    }

    /**
     * Get (and store) the organism Item to reference when creating new objects.
     * @param taxonId the taxon ID of the desired organism
     * @throws ObjectStoreException if there is a problem
     * @return the Organism Item
     */
    Item getOrganism(String taxonId) throws ObjectStoreException {
        if (organism!=null) {
            return organism;
        } else {
            organism = createItem("Organism");
            organism.setAttribute("taxonId", taxonId);
            store(organism);
            return organism;
        }
    }

    /**
     * Get (and store) the strain Item to reference when creating new objects.
     * @param taxonId the taxon ID of the desired strain's organism
     * @param strainIdentifier the identifier of the desired strain
     * @throws ObjectStoreException if there is a problem
     * @return the Strain Item
     */
    Item getStrain(String taxonId, String strainIdentifier) throws ObjectStoreException {
        if (strain!=null) {
            return strain;
        } else {
            strain = createItem("Strain");
            strain.setReference("organism", getOrganism(taxonId));
            strain.setAttribute("identifier", strainIdentifier);
            store(strain);
            return strain;
        }
    }

    /**
     * Create a Sequence and referencing object for the given DNASequence.
     * @param sequenceId the sequence ID
     * @param sequence the DNA sequence
     * @throws ObjectStoreException if store() fails
     */
    void processDNASequence(String sequenceId, DNASequence dnaSequence) throws ObjectStoreException {
        String sequence = dnaSequence.getSequenceAsString();
        String md5checksum = Util.getMd5checksum(sequence);

        System.out.println(getIdentifier(sequenceId));
        
    //     bioSequence.setResidues(new PendingClob(sequence));
    //     bioSequence.setLength(bioJavaSequence.getLength());
    //     bioSequence.setMd5checksum(md5checksum);

    //     // the identifier and name
    //     String identifier = getIdentifier(bioJavaSequence);
    //     String name = getName(bioJavaSequence); // may be null
        
    //     // HACK: don't allow spaces or tabs in sequence primary identifiers; set symbol=extra part
    //     String symbol = null;
    //     String[] spaceChunks = identifier.split(" ");
    //     if (spaceChunks.length>1) {
    //         identifier = spaceChunks[0];
    //         symbol = spaceChunks[1];
    //     }
    //     String[] tabChunks = identifier.split("\t");
    //     if (tabChunks.length>1) {
    //         identifier = tabChunks[0];
    //         symbol = tabChunks[1];
    //     }

    //     // HACK: toggle the className between "Chromosome" and "Supercontig" based on identifier.
    //     // NC_021160.1 is chromosome
    //     // NW_004522746.1 is scaffold/supercontig
    //     if (className.equals("org.intermine.model.bio.Chromosome") || className.equals("org.intermine.model.bio.Supercontig")) {
    //         if (identifier.startsWith("NC")) {
    //             className = "org.intermine.model.bio.Chromosome";
    //         } else if (identifier.startsWith("NW")) {
    //             className = "org.intermine.model.bio.Supercontig";
    //         } else {
    //             throw new RuntimeException("Cannot determine whether ID="+identifier+" is Chromosome or Supercontig.");
    //         }
    //     }

    //     Class<? extends InterMineObject> imClass;
    //     Class<?> c;
    //     try {
    //         c = Class.forName(className);
    //         if (InterMineObject.class.isAssignableFrom(c)) {
    //             imClass = (Class<? extends InterMineObject>) c;
    //         } else {
    //             throw new RuntimeException("Feature className must be a valid class in the model that inherits from InterMineObject, but was: " + className);
    //         }
    //     } catch (ClassNotFoundException e1) {
    //         throw new RuntimeException("unknown class: " + className + " while creating new Sequence object");
    //     }

    //     // create the object that has the sequence
    //     BioEntity imo = (BioEntity) getDirectDataLoader().createObject(imClass);

    //     try {
    //         imo.setFieldValue(classAttribute, identifier);
    //     } catch (Exception e) {
    //         throw new IllegalArgumentException("Error setting: "+className+"."+classAttribute+" to: "+identifier+". Does the attribute exist?");
    //     }
        
    //     try {
    //         imo.setFieldValue("sequence", bioSequence);
    //     } catch (Exception e) {
    //         throw new IllegalArgumentException("Error setting: "+className+".sequence to: "+identifier+". Does the sequence attribute exist?");
    //     }

    //     imo.setOrganism(organism);
    //     if (strain!=null) imo.setStrain(strain);
    //     if (assemblyVersion!=null) imo.setAssemblyVersion(assemblyVersion);
    //     if (annotationVersion!=null) imo.setAnnotationVersion(annotationVersion);
        

    //     try {
    //         imo.setFieldValue("length", new Integer(bioSequence.getLength()));
    //     } catch (Exception e) {
    //         throw new IllegalArgumentException("Error setting: "+className+".length to: "+bioSequence.getLength()+". Does the attribute exist?");
    //     }

    //     try {
    //         imo.setFieldValue("md5checksum", md5checksum);
    //     } catch (Exception e) {
    //         // Ignore - we don't care if the field doesn't exist.
    //     }

    //     try {
    //         if (symbol!=null) imo.setFieldValue("symbol", symbol);
    //     } catch (Exception e) {
    //         throw new IllegalArgumentException("Error setting: "+className+".symbol to: "+identifier+". Does the symbol attribute exist?");
    //     }

    //     try {
    //         if (name!=null) imo.setFieldValue("name", name);
    //     } catch (Exception e) {
    //         throw new IllegalArgumentException("Error setting: "+className+".name to: "+name+". Does the name attribute exist?");
    //     }

    //     if (StringUtils.isEmpty(dataSetTitle)) {
    //         throw new RuntimeException("DataSet title (ncbi-fasta.dataSetTitle) not set.");
    //     }

    //     extraProcessing(bioJavaSequence, bioSequence, imo, organism, strain, getDataSet());

    //     DataSet dataSet = getDataSet();
    //     imo.addDataSets(dataSet);
    //     try {
    //         getDirectDataLoader().store(bioSequence);
    //         getDirectDataLoader().store(imo);
    //         storeCount += 2;
    //     } catch (ObjectStoreException e) {
    //         throw new BuildException("store failed", e);
    //     }
    }

    /**
     * Create a Sequence and referencing object for the given ProteinSequence.
     * @param sequenceId the sequence ID
     * @param sequence the protein sequence
     * @throws ObjectStoreException if store() fails
     */
    void processProteinSequence(String sequenceId, ProteinSequence proteinSequence) throws ObjectStoreException {
        String sequence = proteinSequence.getSequenceAsString();
        String md5checksum = Util.getMd5checksum(sequence);

        System.out.println(getIdentifier(sequenceId));
    }

    /**
     * Return the DataSet to add to each object.
     * @return the DataSet
     * @throws ObjectStoreException if there is an ObjectStore problem
     */
    // public DataSet getDataSet() throws ObjectStoreException {
    //     if (dataSets.containsKey(dataSetTitle)) {
    //         return dataSets.get(dataSetTitle);
    //     }
    //     DataSet dataSet = getDirectDataLoader().createObject(DataSet.class);
    //     dataSet.setName(dataSetTitle);
    //     if (dataSetUrl!=null) dataSet.setUrl(dataSetUrl);
    //     if (dataSetDescription!=null) dataSet.setDescription(dataSetDescription);
    //     if (dataSetVersion!=null) dataSet.setVersion(dataSetVersion);
    //     if (dataSourceName!=null) dataSet.setDataSource(getDataSource());
    //     getDirectDataLoader().store(dataSet);
    //     dataSets.put(dataSetTitle, dataSet);
    //     return dataSet;
    // }


    /**
     * Do any extra processing needed for this record (extra attributes, objects, references, etc.).
     * This method is called before the new objects are stored.
     * @param bioJavaSequence the BioJava Sequence
     * @param imSequence the IntermMine Sequence
     * @param bioEntity the object that references the sequence
     * @param organism the Organism object for the new InterMineObject
     * @param dataSet the DataSet object
     * @throws ObjectStoreException if a store() fails during processing
     *
     * lcl|NC_021160.1_cds_XP_004485403.1_1 [gene=LOC101488545] [db_xref=GeneID:101488545]
                                            [protein=protein phosphatase 1 regulatory subunit INH3-like] [protein_id=XP_004485403.1]
                                            [location=37410..37766] [gbkey=CDS]
     */
    // protected void extraProcessing(Sequence bioJavaSequence, org.intermine.model.bio.Sequence imSequence,
    //                                BioEntity bioEntity, Organism organism, Strain strain, DataSet dataSet) throws ObjectStoreException {
    //     // NOTHING SO FAR
    //     // Can't associate a bioEntity with a gene or protein since it doesn't have the references
    // }

    /**
     * For the given BioJava Sequence object, return an identifier to be used when creating the corresponding BioEntity.
     * 0=identifier   1=name                             N=ignore
     * NP_001265926.1 rab-type small GTP-binding protein [Cicer arietinum]
     * 0.0 0.1=identifier
     * lcl|NC_021160.1_cds_XP_004485403.1_1 [gene=LOC101488545] [db_xref=GeneID:101488545]
                                            [protein=protein phosphatase 1 regulatory subunit INH3-like] [protein_id=XP_004485403.1]
                                            [location=37410..37766] [gbkey=CDS]
     * @param bioJavaSequence the Sequence
     * @return an identifier
     */
    protected String getIdentifier(String header) {
        String[] bits = header.split(" ");
        if (bits[0].contains("|")) {
            String[] subbits = bits[0].split("\\|");
            return subbits[1];
        } else {
            return bits[0];
        }
    }

    /**
     * For the given BioJava Sequence object, return the name attribute.
     * 0=identifier   1=name                             N=ignore
     * NP_001265926.1 rab-type small GTP-binding protein [Cicer arietinum]
     * @param bioJavaSequence the Sequence
     * @return a name
     */
    protected String getName(String header) {
        String[] bits = header.split(" ");
        if (bits.length==1) {
            return null;
        } else {
            String name = "";
            for (int i=1; i<bits.length; i++) {
                if (bits[i].contains("[")) {
                    break;
                } else {
                    name += bits[i]+" ";
                }
            }
            name = name.trim();
            if (name.length()==0) {
                return null;
            } else {
                return name;
            }
        }
    }

    /**
     * Store and/or return the DataSource set in project.xml.
     */
    // DataSource getDataSource() throws ObjectStoreException {
    //     if (StringUtils.isEmpty(dataSourceName)) {
    //         throw new RuntimeException("dataSourceName not set");
    //     }
    //     if (dataSource==null) {
    //         dataSource = getDirectDataLoader().createObject(DataSource.class);
    //         dataSource.setName(dataSourceName);
    //         if (dataSourceUrl!=null) dataSource.setUrl(dataSourceUrl);
    //         if (dataSourceDescription!=null) dataSource.setDescription(dataSourceDescription);
    //         getDirectDataLoader().store(dataSource);
    //         storeCount += 1;
    //     }
    //     return dataSource;
    // }

}
