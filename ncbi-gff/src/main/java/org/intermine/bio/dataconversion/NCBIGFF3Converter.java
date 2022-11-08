package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.Set;

import org.intermine.bio.io.gff3.GFF3Record;
import org.intermine.dataconversion.FileConverter;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.ClassDescriptor;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

import org.apache.log4j.Logger;

/**
 * Class to load Items from an NCBI GFF3 file.
 *
 * @author Sam Hokin
 */
public class NCBIGFF3Converter extends FileConverter {
    private static final Logger LOG = Logger.getLogger(NCBIGFF3Converter.class);

    // put all items in a single map since they're supposed to have distinct identifiers
    Map<String,Item> featureMap = new HashMap<>();

    // use a counter to suffix IDs
    Map<String,Integer> suffixCounter = new HashMap<>();

    // for ensuring unique primary identifiers
    List<String> genePrimaryIdentifiers = new ArrayList<>();

    ItemWriter writer;
    Model model;

    // set in project.xml
    String taxonId;                          // organism
    String strainIdentifier;                 // strain

    // from GFF header
    String dataSourceName;                                     // dataSource
    String dataSetName, dataSetDescription, dataSetLicence;    // dataSet
    String annotationVersion, assemblyVersion;                 // attributes

    // constant per invocation
    Item organism;
    Item strain;
    Item dataSource;
    Item dataSet;
    
    /**
     * Constructor
     * @param writer ItemWriter
     * @param model the model to create items in
     * @throws ObjectStoreException if something goes wrong
     */
    public NCBIGFF3Converter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
        this.writer = writer;
        this.model = model;
    }

    /**
     * Set the data source name.
     */
    public void setDataSourceName(String dataSourceName) {
        this.dataSourceName = dataSourceName;
    }

    /**
     * Set the data set name.
     */
    public void setDataSetName(String dataSetName) {
        this.dataSetName = dataSetName;
    }

    /**
     * Set the data set description.
     */
    public void setDataSetDescription(String dataSetDescription) {
        this.dataSetDescription = dataSetDescription;
    }

    /**
     * Set the data set licence.
     */
    public void setDataSetLicence(String dataSetLicence) {
        this.dataSetLicence = dataSetLicence;
    }

    /**
     * Set the strain identifier.
     */
    public void setStrainIdentifier(String strainIdentifier) {
        this.strainIdentifier = strainIdentifier;
    }

    /**
     * Set the organism taxon ID.
     */
    public void setTaxonId(String taxonId) {
        this.taxonId = taxonId;
    }

    /**
     * Set the annotation version.
     */
    public void setAnnotationVersion(String annotationVersion) {
        this.annotationVersion = annotationVersion;
    }

    /**
     * Set the assembly version.
     */
    public void setAssemblyVersion(String assemblyVersion) {
        this.assemblyVersion = assemblyVersion;
    }
    
    /**
     * Process a BufferedReader
     * @param reader the Reader
     * @throws java.io.IOException if an error occurs reading GFF
     * @throws ObjectStoreException if an error occurs storing items
     */
    public void process(Reader reader) throws IOException, ObjectStoreException {
        // only process GFF files
        if (!getCurrentFile().getName().endsWith("gff3") && !getCurrentFile().getName().endsWith("gff")) return;

        // be sure we've set the required project.xml parameters
        if (taxonId==null) throw new RuntimeException("taxonId is not provided in project.xml.");
        if (strainIdentifier==null) throw new RuntimeException("strainIdentifier is not provided in project.xml.");

        // process GFF lines
        Map<String,Integer> typeCount = new HashMap<>();
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#!")) {
                // #!gff-spec-version 1.21
                // #!processor NCBI annotwriter
                // #!genome-build ASM33114v1
                // #!genome-build-accession NCBI_Assembly:GCF_000331145.1
                // #!annotation-source NCBI Cicer arietinum Annotation Release 102
                if (line.startsWith("#!genome-build-accession")) {
                    dataSetName = line.replace("#!genome-build-accession ", "");
                    assemblyVersion = dataSetName;
                } else if (line.startsWith("#!annotation-source")) {
                    dataSetDescription = line.replace("#!annotation-source ", "");
                    annotationVersion = dataSetDescription;
                }
            } else if (line.startsWith("#")) {
                // comment
            } else {
                GFF3Record record = new GFF3Record(line);
                // keep count of types
                String recordType = record.getType();
                if (typeCount.containsKey(recordType)) {
                    int count = typeCount.get(recordType) + 1;
                    typeCount.put(recordType, count);
                } else {
                    int count = 1;
                    typeCount.put(recordType, count);
                }
                // process the GFF record
                process(record);
            }
        }
        for (String recordType : typeCount.keySet()) {
            System.out.println(typeCount.get(recordType)+" "+recordType+" records");
        }
    }

    /**
     * Process a GFF3 record and give a xml presentation
     * @param record GFF3Record
     * @throws ObjectStoreException if an error occurs storing items
     * @throws IOException
     * String 	                getSequenceID() Return the sequenceID field of this record.
     * String 	                getType()Return the type field of this record.
     * String 	                getId() Return the first value of the Id field from the attributes of this record.
     * int 	                getStart() Return the start field of this record.
     * int 	                getEnd() Return the end field of this record.
     * List<String>             getNames() Return the list of the Name field from the attributes of this record.
     * String 	                getGap() Return the first value of the Gap field from the attributes of this record.
     * String 	                getNote() Return the first value of the Note field from the attributes of this record.
     * String 	                getOntologyTerm() Return the first value of the OntologyTerm field from the attributes of this record.
     * String 	                getFirstAlias() Return the first value of the Alias field from the attributes of this record.
     * String 	                getPhase() Return the phase field of this record.
     * Double 	                getScore() Return the score field of this record.
     * String 	                getSource() Return the source field of this record.
     * String 	                getStrand() Return the strand field of this record.
     * String 	                getTarget() Return the first value of the Target field from the attributes of this record.
     * Map<String,List<String>> getAttributes() Return the attributes of this record as a Map from attribute key to Lists of attribute values.
     * List<String>             getAliases() Return all values of the Alias field from the attributes of this record.
     * List<String>             getDbxrefs() Return the first value of the Dbxref field from the attributes of this record.
     * List<String>             getParents() Return the list of the Parent field from the attributes of this record.
     */
    public void process(GFF3Record record) throws ObjectStoreException {
        if (record.getType().equals("region")) {
            Item feature = getChromosomeOrSupercontig(record);
        } else if (record.getType().equals("gene")) {
            Item feature = getGene(record);
        } else if (record.getType().equals("mRNA")) {
            Item feature = getMRNA(record);
        } else if (record.getType().equals("exon")) {
            Item feature = getExon(record);
        } else if (record.getType().contains("RNA")) {
            Item feature = getNcRNA(record);
        } else if (record.getType().equals("transcript")) {
            Item feature = getTranscript(record);
        }
    }

    /**
     * Return true if the given id identifies a chromosome.
     */
    public static boolean isChromosome(String sequenceId) {
        return sequenceId.startsWith("NC_");
    }

    /**
     * Return true if the given id identifies a supercontig.
     */
    public static boolean isSupercontig(String sequenceId) {
        return sequenceId.startsWith("NW_");
    }

    /**
     * Create and/or return the given Chromosome or Supercontig Item given only the sequence ID.
     */
    public Item getChromosomeOrSupercontig(String sequenceId) throws ObjectStoreException {
        if (featureMap.containsKey(sequenceId)) {
            return featureMap.get(sequenceId);
        } else {
            Item feature = null;
            if (isChromosome(sequenceId)) {
                feature = createItem("Chromosome");
            } else if (isSupercontig(sequenceId)) {
                feature = createItem("Supercontig");
            } else {
                throw new RuntimeException("Sequence ID is neither a chromosome nor supercontig:"+sequenceId);
            }
            feature.setAttribute("primaryIdentifier", sequenceId);
            feature.setReference("organism", getOrganism());
            feature.setReference("strain", getStrain());
            feature.addToCollection("dataSets", getDataSet());
            feature.setAttribute("assemblyVersion", assemblyVersion);
            feature.setAttribute("annotationVersion", annotationVersion);
            featureMap.put(sequenceId, feature);
            return feature;
        }
    }
    
    /**
     * Create or update and return the given Chromosome or Supercontig Item given a GFF3Record.
     * NC_021160.1 RefSeq region 1 48359943 . + . ID=NC_021160.1:1..48359943;Dbxref=taxon:3827;Name=Ca1;chromosome=Ca1;
     *                                            collection-date=Oct-2011;country=India;cultivar=CDC Frontier;gbkey=Src;genome=chromosome;
     *                                            mol_type=genomic DNA;tissue-type=etiolated seedlings
     *
     * NC_023759.1 RefSeq region 1 52205531 . + . ID=id0;Dbxref=taxon:3885;Name=1;chromosome=1;cultivar=G19833;gbkey=Src;genome=chromosome;mol_type=genomic DNA
     *
     * NOTE: use the sequence ID, not the ID from the attributes, which varies unpredictably!!
     */
    public Item getChromosomeOrSupercontig(GFF3Record record) throws ObjectStoreException {
        Item feature = null;
        String primaryIdentifier = record.getSequenceID();
        if (featureMap.containsKey(primaryIdentifier)) {
            feature = featureMap.get(primaryIdentifier);
        } else {
            if (isChromosome(primaryIdentifier)) {
                feature = createItem("Chromosome");
            } else if (isSupercontig(primaryIdentifier)) {
                feature = createItem("Supercontig");
            } else {
                throw new RuntimeException("GFF record is neither a chromosome nor supercontig:"+record.toString());
            }
        }
        feature.setAttribute("primaryIdentifier", primaryIdentifier);
        feature.setReference("organism", getOrganism());
        feature.setReference("strain", getStrain());
        feature.addToCollection("dataSets", getDataSet());
        feature.setAttribute("assemblyVersion", assemblyVersion);
        feature.setAttribute("annotationVersion", annotationVersion);
        feature.setAttribute("length", String.valueOf(record.getEnd()));
        if (record.getNames().size()>0) feature.setAttribute("secondaryIdentifier", record.getNames().get(0));
        if (record.getScore()!=null) feature.setAttribute("score", record.getScore().toString());
        featureMap.put(primaryIdentifier, feature); // update
        return feature;
    }

    /**
     * Generic routine for getting a feature item from a GFF record. Updates an existing record.
     * NOTE: the featureMap uses the ID (e.g. gene-FOO, gene1) to enforce uniqueness.
     *
     * sequenceID: NC_021167.1 source: Gnomon type: gene start: 16256160 end: 16258858 score: null strand: - phase: null
     * attributes: ID=gene-LOC101506378
     *             Dbxref=GeneID:101506378
     *             Name=LOC101506378
     *             gbkey=Gene
     *             gene=LOC101506378
     *             gene_biotype=protein_coding
     *
     * sequenceID: NC_021160.1 source: Gnomon type: exon start: 834 end: 982 score: null strand: - phase: null 
     * attributes: ID=[exon-XR_001143330.2-1]
     *             Parent=[rna-XR_001143330.2]
     *             Dbxref=[GeneID:101497325, Genbank:XR_001143330.2]
     *             gbkey=[ncRNA]
     *             gene=[LOC101497325]
     *             product=[uncharacterized LOC101497325]
     *             transcript_id=[XR_001143330.2]
     *
     * sequenceID: NC_023759.1 source: RefSeq type: exon start: 4553 end: 4760 score: null strand: - phase: null
     * attributes: ID=[id1]
     *             Parent=[rna0]
     *             Dbxref=[Genbank:XM_007160519.1, GeneID:18642687]
     *             Note=[transcript PHAVU_001G00010012m]
     *             end_range=[4760, .]
     *             gbkey=[mRNA]
     *             partial=[true]
     *             product=[hypothetical protein]
     *             transcript_id=[XM_007160519.1]
     */
    Item getFeatureItem(String itemClass, GFF3Record record) throws ObjectStoreException {
        Item feature = null;
        if (featureMap.containsKey(record.getId())) {
            feature = featureMap.get(record.getId());
        } else {
            feature = createItem(itemClass);
            featureMap.put(record.getId(), feature);
        }
        // use the Name attribute for primaryIdentifier, unless it's already been used for a gene, in which case use ID
        if (record.getNames()!=null && record.getNames().size()>0) {
            String name = record.getNames().get(0);
            if (record.getType().equals("gene")) {
                if (genePrimaryIdentifiers.contains(name)) {
                    feature.setAttribute("primaryIdentifier", record.getId());
                } else {
                    feature.setAttribute("primaryIdentifier", name);
                    genePrimaryIdentifiers.add(name);
                }
            } else {
                // presuming Name uniqueness for non-genes
                feature.setAttribute("primaryIdentifier", name);
            }
        } else {
            // append counter to gene or transcript_id
            String counterKey = null;
            Map<String,List<String>> attributes = record.getAttributes();
            if (attributes.containsKey("gene")) {
                // use gene
                counterKey = attributes.get("gene").get(0);
            } else if (attributes.containsKey("transcript_id")) {
                // use transcript
                counterKey = attributes.get("transcript_id").get(0);
            }
            if (counterKey!=null) {
                String primaryIdentifier = counterKey;
                int count = 1;
                if (suffixCounter.containsKey(counterKey)) {
                    count = suffixCounter.get(counterKey) + 1;
                }
                suffixCounter.put(counterKey, count);
                primaryIdentifier += "."+count;
                feature.setAttribute("primaryIdentifier", primaryIdentifier);
            } else {
                // last resort, use ID
                feature.setAttribute("primaryIdentifier", record.getId());
            }                
        }
        feature.setReference("organism", getOrganism());
        feature.setReference("strain", getStrain());
        feature.addToCollection("dataSets", getDataSet());
        feature.setAttribute("assemblyVersion", assemblyVersion);
        feature.setAttribute("annotationVersion", annotationVersion);
        feature.setAttribute("length", String.valueOf(record.getEnd()-record.getStart()));
        setChromosomeOrSupercontigReference(feature, record.getSequenceID());
        setChromosomeOrSupercontigLocation(feature, record);
        if (record.getScore()!=null) feature.setAttribute("score", String.valueOf(record.getScore()));
        return feature;
    }

    /**
     * Generic routine for getting a feature item that may or may not already exists with the given ID (map key).
     */
    Item getFeatureItem(String itemClass, String id, String sequenceId) throws ObjectStoreException {
        if (featureMap.containsKey(id)) {
            return featureMap.get(id);
        } else {
            // create a stub feature with the given ID as the map key AND primaryIdentifier (which will hopefully be updated)
            Item feature = createItem(itemClass);
            featureMap.put(id, feature);
            feature.setAttribute("primaryIdentifier", id);
            feature.setReference("organism", getOrganism());
            feature.setReference("strain", getStrain());
            feature.addToCollection("dataSets", getDataSet());
            feature.setAttribute("assemblyVersion", assemblyVersion);
            feature.setAttribute("annotationVersion", annotationVersion);
            setChromosomeOrSupercontigReference(feature, sequenceId);
            return feature;
        }
    }

    /**
     * Return or create and return the given Gene Item.
     * NC_021160.1 Gnomon gene 91418 106326 . + . ID=gene-LOC101490835;Dbxref=GeneID:101490835;
     *                                            Name=LOC101490835;gbkey=Gene;gene=LOC101490835;gene_biotype=protein_coding
     */
    public Item getGene(GFF3Record record) throws ObjectStoreException {
        Item feature = getFeatureItem("Gene", record);
        return feature;
    }
    
    /**
     * Return or create and return the given mRNA Item.
     * NC_021160.1 Gnomon mRNA 91418 106326 . + . ID=rna-XM_004485354.3;Parent=gene-LOC101490835;Dbxref=GeneID:101490835,Genbank:XM_004485354.3;
     *                                            Name=XM_004485354.3;gbkey=mRNA;gene=LOC101490835;
     *                                            model_evidence=Supporting evidence includes similarity to: 100%25 coverage of the annotated...;
     *                                            product=squamous cell carcinoma antigen recognized by T-cells 3%2C transcript variant X1;
     *                                            transcript_id=XM_004485354.3
     */
    Item getMRNA(GFF3Record record) throws ObjectStoreException {
        Item feature = getFeatureItem("MRNA", record);
        List<String> parents = record.getAttributes().get("Parent");
        if (parents.size()>0) {
            String parent = parents.get(0);
            feature.setReference("gene", getFeatureItem("Gene", parent, record.getSequenceID()));
        }
        return feature;
    }

    /**
     * Return or create and return the given ncRNA Item.
     */
    Item getNcRNA(GFF3Record record) throws ObjectStoreException {
        Item feature = getFeatureItem("NcRNA", record);
        List<String> parents = record.getAttributes().get("Parent");
        if (parents.size()>0) {
            String parent = parents.get(0);
            feature.setReference("gene", getFeatureItem("Gene", parent, record.getSequenceID()));
        }
        return feature;
    }

    /**
     * Return or create and return the given Transcript Item.
     */
    Item getTranscript(GFF3Record record) throws ObjectStoreException {
        Item feature = getFeatureItem("Transcript", record);
        List<String> parents = record.getAttributes().get("Parent");
        if (parents.size()>0) {
            String parent = parents.get(0);
            feature.setReference("gene", getFeatureItem("Gene", parent, record.getSequenceID()));
        }
        return feature;
    }

    /**
     * Return or create and return the given Exon Item.
     * ID=exon-XM_004485349.3-8;Parent=rna-XM_004485349.3
     */
    Item getExon(GFF3Record record) throws ObjectStoreException {
        Item feature = getFeatureItem("Exon", record);
        String gbKey = record.getAttributes().get("gbkey").get(0);
        String parent = record.getAttributes().get("Parent").get(0);
        if (gbKey.equals("mRNA")) {
            feature.setReference("transcript", getFeatureItem("MRNA", parent, record.getSequenceID()));
        } else if (gbKey.equals("ncRNA")) {
            feature.setReference("transcript", getFeatureItem("NcRNA", parent, record.getSequenceID()));
        } else if (gbKey.equals("tRNA")) {
            feature.setReference("transcript", getFeatureItem("TRNA", parent, record.getSequenceID()));
        } else if (gbKey.equals("rRNA")) {
            feature.setReference("transcript", getFeatureItem("RRNA", parent, record.getSequenceID()));
        } else if (gbKey.equals("misc_RNA")) {
            feature.setReference("transcript", getFeatureItem("Transcript", parent, record.getSequenceID()));
        } else if (gbKey.equals("exon")) {
            // do nothing
        } else {
            System.out.println("Unusual exon gbkey="+gbKey);
        }
        return feature;
    }

    /**
     * Perform any necessary clean-up after post-conversion
     * @throws Exception if an error occurs
     */
    @Override
    public void close() throws Exception {
        store(featureMap.values());
    }

    /**
     * Return the organism Item created for this NCBIGFF3Converter from the organism abbreviation passed
     * to the constructor.
     * @return the organism item
     * @throws ObjectStoreException if the Organism item can't be stored
     */
    public Item getOrganism() throws ObjectStoreException {
        if (organism==null && taxonId!=null) {
            organism = createItem("Organism");
            organism.setAttribute("taxonId", taxonId);
            store(organism);
        }
        return organism;
    }

    /**
     * Return the strain Item created for this NCBIGFF3Converter from the strain passed
     * to the constructor.
     * @param organism the Organism that this Strain references
     * @return the strain item
     * @throws ObjectStoreException if the Strain item can't be stored
     */
    public Item getStrain() throws ObjectStoreException {
        if (strain==null && strainIdentifier!=null && organism!=null) {
            strain = createItem("Strain");
            strain.setAttribute("identifier", strainIdentifier);
	    strain.setReference("organism", organism);
            store(strain);
        }
        return strain;
    }

    /**
     * Create and add a synonym Item from the given information.
     * @param subject the subject of the new Synonym
     * @param value the Synonym value
     * @return the new Synonym Item
     */
    public Item getSynonym(Item subject, String value) {
        Item synonym = createItem("Synonym");
        synonym.setAttribute("value", value);
        synonym.setReference("subject", subject.getIdentifier());
        return synonym;
    }

    /**
     * Return a DataSource item for the given name
     * @return the DataSource Item
     * @throws ObjectStoreException
     */
    public Item getDataSource() throws ObjectStoreException, RuntimeException {
        if (dataSourceName==null || dataSourceName.trim().length()==0) {
            throw new RuntimeException("dataSourceName must be set in project.xml.");
        }
        if (dataSource==null) {
            dataSource = createItem("DataSource");
            dataSource.setAttribute("name", dataSourceName);
            store(dataSource);
        }
        return dataSource;
    }

    /**
     * Return a DataSet item with the given details.
     * @return the DataSet Item
     * @throws ObjectStoreException
     */
    public Item getDataSet() throws ObjectStoreException {
        if (dataSetName==null || dataSetName.trim().length()==0) {
            throw new RuntimeException("dataSetName must be set in project.xml.");
        }
        if (dataSet==null) {
            dataSet = createItem("DataSet");
            dataSet.setAttribute("name", dataSetName);
	    dataSet.setReference("dataSource", getDataSource());
            if (dataSetDescription!=null) dataSet.setAttribute("description", dataSetDescription);
            if (dataSetLicence!=null) dataSet.setAttribute("licence", dataSetLicence);
	    // if (url!=null) dataSet.setAttribute("url", url);
	    // if (version!=null) dataSet.setAttribute("version", version);
            store(dataSet);
        }
        return dataSet;
    }

    /**
     * Set a reference to either a Chromosome or Supercontig identified by sequenceId.
     */
    void setChromosomeOrSupercontigReference(Item feature, String sequenceId) throws ObjectStoreException {
        if (isChromosome(sequenceId)) {
            feature.setReference("chromosome", getChromosomeOrSupercontig(sequenceId));
        } else if (isSupercontig(sequenceId)) {
            feature.setReference("supercontig", getChromosomeOrSupercontig(sequenceId));
        } else {
            throw new RuntimeException("Sequence ID "+sequenceId+" is neither a Chromosome nor Supercontig.");
        }
    }

    /**
     * Create a Location to a Chromosome or Supercontig and set the reference for the given feature.
     */
    void setChromosomeOrSupercontigLocation(Item feature, GFF3Record record) throws ObjectStoreException {
        Item location = createItem("Location");
        location.setReference("locatedOn", getChromosomeOrSupercontig(record.getSequenceID()));
        location.setAttribute("start", String.valueOf(record.getStart()));
        location.setAttribute("end", String.valueOf(record.getEnd()));
        location.setAttribute("strand", record.getStrand());
        location.setReference("feature", feature);
        store(location);
        if (isChromosome(record.getSequenceID())) {
            feature.setReference("chromosomeLocation", location);
        } else if (isSupercontig(record.getSequenceID())) {
            feature.setReference("supercontigLocation", location);
        } else {
            throw new RuntimeException("Sequence ID "+record.getSequenceID()+" neither a Chromosome nor Supercontig.");
        }
    }
}
