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

import java.util.Map;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.NoSuchElementException;

import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;

import org.biojava.nbio.core.exceptions.ParserException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.template.AbstractSequence;

import org.intermine.metadata.Util;
import org.intermine.model.InterMineObject;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.objectstore.query.PendingClob;
import org.intermine.task.FileDirectDataLoaderTask;

import org.intermine.model.bio.BioEntity;
import org.intermine.model.bio.CDS;
import org.intermine.model.bio.Chromosome;
import org.intermine.model.bio.DataSet;
import org.intermine.model.bio.DataSource;
import org.intermine.model.bio.Gene;
import org.intermine.model.bio.Organism;
import org.intermine.model.bio.Protein;
import org.intermine.model.bio.Strain;
import org.intermine.model.bio.Sequence;
import org.intermine.model.bio.Supercontig;

/**
 * A task that can read a set of FASTA files and create the corresponding Sequence objects in an ObjectStore.
 *
 * @author Kim Rutherford
 * @author Peter Mclaren
 * @author Sam Hokin
 */

public class NCBIFastaLoaderTask extends FileDirectDataLoaderTask {
    static final Logger LOG = Logger.getLogger(NCBIFastaLoaderTask.class);

    String sequenceType = "dna";

    String taxonId = null;
    String strainIdentifier;
    String assemblyVersion;   // for proteins, CDS
    String annotationVersion; // for proteins, CDS
    String className;

    String dataSourceName;
    String dataSourceUrl;
    String dataSourceDescription;

    String dataSetName;
    String dataSetUrl;
    String dataSetDescription; 
    String dataSetVersion;

    Organism organism;
    Strain strain;
    DataSource dataSource;
    DataSet dataSet;

    Map<String,Gene> genes = new HashMap<>();
    Map<String,Protein> proteins = new HashMap<>();

    //Set this if we want to do some testing...
    File[] files = null;

    /**
     * Set the sequence type to be passed to the FASTA parser.  The default is "dna".
     * @param sequenceType the sequence type
     */
    public void setSequenceType(String sequenceType) {
        if ("${fasta.sequenceType}".equals(sequenceType)) {
            this.sequenceType = "dna";
        } else {
            this.sequenceType = sequenceType;
        }
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
     * The default class name to use for objects created during load.  Generally this is
     * "org.intermine.model.bio.Chromosome" or "org.intermine.model.bio.Protein"
     * @param className the class name
     */
    public void setClassName(String className) {
        this.className = className;
    }

    /**
     * Return the class name set with setClassName().
     * @return the class name
     */
    public String getClassName() {
        return className;
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
     * @param dataSetName the title of the DataSet of any new features
     */
    public void setDataSetName(String dataSetName) {
        this.dataSetName = dataSetName;
    }

    /**
     * If a value is specified this url will used when a DataSet is created.
     * @param dataSetUrl the title of the DataSet of any new features
     */
    public void setDataSetUrl(String dataSetUrl) {
        this.dataSetUrl = dataSetUrl;
    }

    /**
     * If a value is specified this description will used when a DataSet is created.
     * @param dataSetDescription the description of the DataSet of any new features
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
     * Directly set the array of files to read from.  Use this for testing with junit.
     * @param files the File objects
     */
    protected void setFileArray(File[] files) {
        this.files = files;
    }

    /**
     * Process and load all of the fasta files.
     */
    @Override
    public void process() {
        try {
            super.process();
            getIntegrationWriter().commitTransaction();
            getIntegrationWriter().beginTransaction();
            getDirectDataLoader().close();
        } catch (ObjectStoreException e) {
            throw new BuildException("Failed to store object", e);
        }
    }

    /**
     * Be sure to close the data loader so the last batch gets stored. only needed for tests
     * since the data loading task usually does that for hte live builds.
     * @throws ObjectStoreException if we can't store to db
     */
    public void close() throws ObjectStoreException {
        // store any data left over
        getDirectDataLoader().close();
    }

    /**
     * @throws BuildException if an ObjectStore method fails
     */
    @Override
    public void execute() {
        // don't configure dynamic attributes if this is a unit test!
        if (getProject()!=null) {
            configureDynamicAttributes(this);
        }
        // required project.xml parameters
        if (className==null || className.trim().length()==0) {
            throw new RuntimeException("className must be set in project.xml.");
        }
        if (dataSourceName==null || dataSourceName.trim().length()==0) {
            throw new RuntimeException("dataSourceName must be set in project.xml.");
        }
        if (dataSetName==null || dataSetName.trim().length()==0) {
            throw new RuntimeException("dataSetName must be set in project.xml.");
        }
        // optional but recommended project.xml parameters
        if (strainIdentifier==null || strainIdentifier.trim().length()==0) {
            System.out.println("NOTE: strainIdentifier is missing in project.xml.");
        }
        if (assemblyVersion==null || assemblyVersion.trim().length()==0) {
            System.out.println("NOTE: assemblyVersion is missing in projext.xml.");
        }
        if (annotationVersion==null || annotationVersion.trim().length()==0) {
            System.out.println("NOTE: annotationVersion is missing in project.xml.");
        }
        // this will call processFile() for each file
        super.execute();
    }

    /**
     * Handles each fasta file. Factored out so we can supply files for testing.
     *
     * @param file the File to process.
     * @throws BuildException if the is a problem
     */
    @Override
    public void processFile(File file) {
        System.out.println("##################################################################################################################################################");
        System.out.println("Reading "+sequenceType+" sequences from: "+file);
        System.out.println("##################################################################################################################################################");
        LOG.info("NCBIFastaLoaderTask loading file "+file.getName());
        try {
            // process FASTA file
            if (sequenceType.equalsIgnoreCase("dna")) {
                LinkedHashMap<String,DNASequence> sequenceMap = FastaReaderHelper.readFastaDNASequence(file);
                for (DNASequence sequence : sequenceMap.values()) {
                    processSequence(sequence);
                }
            } else if (sequenceType.equalsIgnoreCase("protein")) {
                LinkedHashMap<String,ProteinSequence> sequenceMap = FastaReaderHelper.readFastaProteinSequence(file);
                for (ProteinSequence sequence : sequenceMap.values()) {
                    processSequence(sequence);
                }
            } else {
                throw new RuntimeException("Sequence type set in project.xml is neither dna nor protein.");
            }
        } catch (ParserException e) {
            throw new BuildException("Sequence not in FASTA format or wrong alphabet for: "+file, e);
        } catch (NoSuchElementException e) {
            throw new BuildException("No FASTA sequences in: "+file, e);
        } catch (FileNotFoundException e) {
            throw new BuildException("Problem reading file - file not found: "+file, e);
        } catch (ObjectStoreException e) {
            throw new BuildException("ObjectStore problem while processing: "+file, e);
        } catch (IOException e) {
            throw new BuildException("Error while closing FileReader for: "+file, e);
        }
    }

    /**
     * Get and store() the Organism object to reference when creating new objects.
     * @param bioJavaSequence the biojava sequence to be parsed
     * @throws ObjectStoreException if there is a problem
     * @return the new Organism
     */
    protected Organism getOrganism() throws ObjectStoreException {
        if (organism==null) {
            organism = getDirectDataLoader().createObject(Organism.class);
            organism.setTaxonId(taxonId);
            getDirectDataLoader().store(organism);
        }
        return organism;
    }

    /**
     * Get and store() the Strain object to reference when creating new objects.
     * @param bioJavaSequence the biojava sequence to be parsed (not used)
     * @throws ObjectStoreException if there is a problem
     * @return the new Strain
     */
    protected Strain getStrain() throws ObjectStoreException {
        if (strain==null && strainIdentifier!=null) {
            strain = getDirectDataLoader().createObject(Strain.class);
            strain.setIdentifier(strainIdentifier);
            strain.setOrganism(getOrganism());
            getDirectDataLoader().store(strain);
        }
        return strain;
    }

    /**
     * Store and/or return the DataSource set in project.xml.
     * @return the DataSet
     * @throws ObjectStoreException if there is an ObjectStore problem
     */
    DataSet getDataSet() throws ObjectStoreException {
        if (dataSet==null) {
            dataSet = getDirectDataLoader().createObject(DataSet.class);
            dataSet.setName(dataSetName);
            dataSet.setDataSource(getDataSource());
            if (dataSetUrl!=null) dataSet.setUrl(dataSetUrl);
            if (dataSetDescription!=null) dataSet.setDescription(dataSetDescription);
            if (dataSetVersion!=null) dataSet.setVersion(dataSetVersion);
            getDirectDataLoader().store(dataSet);
        }
        return dataSet;
    }

    /**
     * Store and/or return the DataSource set in project.xml.
     * @return the DataSource
     * @throws ObjectStoreException if there is a problem storing it.
     */
    DataSource getDataSource() throws ObjectStoreException {
        if (dataSource==null) {
            dataSource = getDirectDataLoader().createObject(DataSource.class);
            dataSource.setName(dataSourceName);
            if (dataSourceUrl!=null) dataSource.setUrl(dataSourceUrl);
            if (dataSourceDescription!=null) dataSource.setDescription(dataSourceDescription);
            getDirectDataLoader().store(dataSource);
        }
        return dataSource;
    }

    /**
     * Store and/or return the Gene given by the identifier.
     * @return the Gene
     * @throws ObjectStoreException
     */
    Gene getGene(String identifier) throws ObjectStoreException {
        if (genes.containsKey(identifier)) {
            return genes.get(identifier);
        } else {
            Gene gene = getDirectDataLoader().createObject(Gene.class);
            gene.setPrimaryIdentifier(identifier);
            getDirectDataLoader().store(gene);
            genes.put(identifier, gene);
            return gene;
        }
    }

    /**
     * Store and/or return the Gene given by the identifier, plus assign the name to the protein name and add the given Protein to the genesproteins collection.
     * @return the Gene
     * @throws ObjectStoreException
     */
    Gene getGene(String identifier, Protein protein, String proteinName) throws ObjectStoreException {
        if (genes.containsKey(identifier)) {
            return genes.get(identifier);
        } else {
            Gene gene = getDirectDataLoader().createObject(Gene.class);
            gene.setPrimaryIdentifier(identifier);
            gene.setName(proteinName);
            gene.addProteins(protein);
            getDirectDataLoader().store(gene);
            genes.put(identifier, gene);
            return gene;
        }
    }

    /**
     * Store and/or return the Protein given by the identifier.
     * @return the Protein
     * @throws ObjectStoreException
     */
    Protein getProtein(String identifier) throws ObjectStoreException {
        if (proteins.containsKey(identifier)) {
            return proteins.get(identifier);
        } else {
            Protein protein = getDirectDataLoader().createObject(Protein.class);
            protein.setPrimaryIdentifier(identifier);
            getDirectDataLoader().store(protein);
            proteins.put(identifier, protein);
            return protein;
        }
    }

    /**
     * Create a Sequence and an object of type className for the given BioJava Sequence.
     * @param bioJavaSequence the AbstractSequence object, either DNASequence or ProteinSequence
     * @throws ObjectStoreException if store() fails
     */
    void processSequence(AbstractSequence bioJavaSequence) throws ObjectStoreException {
        // the sequence string and MD5
        String sequence = bioJavaSequence.getSequenceAsString();
        String md5checksum = Util.getMd5checksum(sequence);
        
        // an InterMine Sequence
        Sequence bioSequence = getDirectDataLoader().createObject(org.intermine.model.bio.Sequence.class);
        bioSequence.setResidues(new PendingClob(sequence));
        bioSequence.setLength(bioJavaSequence.getLength());
        bioSequence.setMd5checksum(md5checksum);

        // the feature identifier and name
        String identifier = getIdentifier(bioJavaSequence);
        String name = getName(bioJavaSequence); // may be null
        
        // HACK: don't allow spaces or tabs in primary identifiers; set symbol=extra part
        String symbol = null;
        String[] spaceChunks = identifier.split(" ");
        if (spaceChunks.length>1) {
            identifier = spaceChunks[0];
            symbol = spaceChunks[1];
        }
        String[] tabChunks = identifier.split("\t");
        if (tabChunks.length>1) {
            identifier = tabChunks[0];
            symbol = tabChunks[1];
        }

        // HACK: set the className to "Chromosome" or "Supercontig" based on identifier.
        // NC_021160.1    is chromosome
        // NW_004522746.1 is scaffold=supercontig
        if (className.equals("org.intermine.model.bio.Chromosome") || className.equals("org.intermine.model.bio.Supercontig")) {
            if (identifier.startsWith("NC")) {
                className = "org.intermine.model.bio.Chromosome";
            } else if (identifier.startsWith("NW")) {
                className = "org.intermine.model.bio.Supercontig";
            } else {
                throw new RuntimeException("Cannot determine whether ID="+identifier+" is Chromosome or Supercontig.");
            }
        }

        // create the feature class
        Class<? extends InterMineObject> imClass;
        Class<?> c;
        try {
            c = Class.forName(className);
            if (InterMineObject.class.isAssignableFrom(c)) {
                imClass = (Class<? extends InterMineObject>) c;
            } else {
                throw new RuntimeException("Feature className must be a valid class in the model that inherits from InterMineObject, but was "+className);
            }
        } catch (ClassNotFoundException e1) {
            throw new RuntimeException("Unknown class: "+className+" while creating new Sequence object");
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // do the work separately for each class since some objects have special attributes to assign from header data //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if (className.equals("org.intermine.model.bio.Chromosome")) {

            // NC_021160.1 Cicer arietinum cultivar CDC Frontier chromosome Ca1, ASM33114v1, whole genome shotgun sequence
            Chromosome feature = (Chromosome) getDirectDataLoader().createObject(imClass);
            feature.setPrimaryIdentifier(identifier);
            setCommonAttributes(feature);
            if (name!=null) feature.setName(name);
            store(feature, bioSequence);

        } else if (className.equals("org.intermine.model.bio.Supercontig")) {

            // NW_004515636.1 Cicer arietinum cultivar CDC Frontier unplaced genomic scaffold, ASM33114v1 scaffold1, whole genome shotgun sequence
            Supercontig feature = (Supercontig) getDirectDataLoader().createObject(imClass);
            feature.setPrimaryIdentifier(identifier);
            setCommonAttributes(feature);
            if (name!=null) feature.setName(name);
            store(feature, bioSequence);

        } else if (className.equals("org.intermine.model.bio.CDS")) {

            // lcl|NW_004522122.1_cds_XP_004517134.1_35487 [gene=LOC101515228] [db_xref=GeneID:101515228]
            //                                             [protein=sucrose transport protein-like] [protein_id=XP_004517134.1]
            //                                             [location=join(17..211,356..419,536..579,665..844)] [gbkey=CDS]
            CDS feature = (CDS) getDirectDataLoader().createObject(imClass);
            feature.setPrimaryIdentifier(identifier);
            setCommonAttributes(feature);
            // gene, protein, gene-protein
            Protein protein = null;
            String proteinId = getBracketIdentifier(bioJavaSequence, "protein_id");
            if (proteinId!=null) {
                protein = getProtein(proteinId);
                feature.setProtein(protein);
            }
            String geneId = getBracketIdentifier(bioJavaSequence, "gene");
            if (geneId!=null) {
                Gene gene = null;
                if (protein==null) {
                    gene = getGene(geneId);
                } else {
                    String proteinName = getBracketIdentifier(bioJavaSequence, "protein");
                    gene = getGene(geneId, protein, proteinName);
                }
                feature.setGene(gene);
            }
            store(feature, bioSequence);

        } else if (className.equals("org.intermine.model.bio.Protein")) {

            // XP_027192614.1 SWI/SNF-related matrix-associated actin-dependent regulator of chromatin subfamily A-like protein 1 isoform X1 [Cicer arietinum]
            Protein feature = (Protein) getDirectDataLoader().createObject(imClass);
            feature.setPrimaryIdentifier(identifier);
            setCommonAttributes(feature);
            if (name!=null) feature.setName(name);
            if (symbol!=null) feature.setSymbol(symbol);
            store(feature, bioSequence);

        } else {
            throw new RuntimeException("Loading of "+className+" from FASTA isn't currently supported.");
        }
    }

    /**
     * Store a feature and its sequence.
     */
    void store(BioEntity feature, Sequence bioSequence) throws ObjectStoreException {
        feature.setFieldValue("sequence", bioSequence);
        feature.setFieldValue("length", new Integer(bioSequence.getLength()));
        getDirectDataLoader().store(bioSequence);
        getDirectDataLoader().store(feature);
    }
    
    /**
     * Set the attributes and references common to all BioEntity objects.
     */
    void setCommonAttributes(BioEntity feature) throws ObjectStoreException {
        feature.addDataSets(getDataSet());
        feature.setOrganism(getOrganism());
        if (strainIdentifier!=null) feature.setStrain(getStrain());
        if (assemblyVersion!=null) feature.setAssemblyVersion(assemblyVersion);
        if (annotationVersion!=null) feature.setAnnotationVersion(annotationVersion);
    }

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
    protected String getIdentifier(AbstractSequence bioJavaSequence) {
        String identifier = null;
        String header = bioJavaSequence.getAccession().getID();
        String[] bits = header.split(" ");
        if (bits[0].contains("|")) {
            String[] subbits = bits[0].split("\\|");
            identifier = subbits[1];
        } else {
            identifier = bits[0];
        }
        return identifier;
    }

    /**
     * Return the identifier stored in brackets with the given field name; else null.
     *
     * lcl|NW_004522122.1_cds_XP_004517134.1_35487 [gene=LOC101515228] [db_xref=GeneID:101515228]
     *                                             [protein=sucrose transport protein-like] [protein_id=XP_004517134.1]
     *                                             [location=join(17..211,356..419,536..579,665..844)] [gbkey=CDS]
     * @param bioJavaSequence the AbstractSequence
     * @param field the name of the desired field, e.g. protein_id
     * @return an identifier
     */
    protected String getBracketIdentifier(AbstractSequence bioJavaSequence, String field) {
        String header = bioJavaSequence.getAccession().getID();
        if (!header.contains(field)) return null;
        String[] split = header.split(field+"=");
        if (split.length==1) return null; // could be doesn't have protein_id= even though has protein_id.
        String secondHalf = split[1];
        String[] parts = secondHalf.split("\\]");
        return parts[0];
    }

    /**
     * For the given BioJava Sequence object, return the name attribute.
     * 0=identifier   1=name                             N=ignore
     * NP_001265926.1 rab-type small GTP-binding protein [Cicer arietinum]
     * @param bioJavaSequence the Sequence
     * @return a name
     */
    protected String getName(AbstractSequence bioJavaSequence) {
        String header = bioJavaSequence.getAccession().getID();
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
}
