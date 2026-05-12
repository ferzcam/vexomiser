/**
 * Semantic similarity baseline (Resnik/Lin + BMA) on the PAVS benchmark.
 *
 * Based on multihop-gda/semantic_similarity.groovy. Key changes for PAVS:
 *   - Loads upheno.owl (RDF_XML, same as original)
 *   - Gene phenotypes: training-case HPO terms (instead of MGI MP terms)
 *   - "Disease" phenotypes: patient HPO terms from test cases directly
 *     (instead of phenotype.hpoa lookup)
 *   - Output TSV matches eval/metrics.py format
 *
 * Usage:
 *   groovy eval/semsim_baseline.groovy \
 *     --data-dir data \
 *     --upheno /path/to/upheno.owl \
 *     --ic resnik --pw resnik --gw bma \
 *     --split test \
 *     --out data/results/semsim_resnik_bma_test.tsv
 */

@Grab(group='com.github.sharispe', module='slib-sml',            version='0.9.1')
@Grab(group='net.sourceforge.owlapi', module='owlapi-api',        version='4.2.5')
@Grab(group='net.sourceforge.owlapi', module='owlapi-apibinding', version='4.2.5')
@Grab(group='net.sourceforge.owlapi', module='owlapi-impl',       version='4.2.5')
@Grab(group='ch.qos.logback',  module='logback-classic',          version='1.2.3')
@Grab(group='org.slf4j',       module='slf4j-api',                version='1.7.30')
@Grab(group='org.codehaus.gpars', module='gpars',                 version='1.1.0')

import org.semanticweb.owlapi.model.*
import org.semanticweb.owlapi.apibinding.OWLManager

import slib.sml.sm.core.engine.SM_Engine
import slib.sml.sm.core.metrics.ic.utils.*
import slib.sml.sm.core.utils.SMConstants
import slib.sml.sm.core.utils.SMconf
import slib.graph.model.impl.graph.memory.GraphMemory
import slib.graph.io.conf.GDataConf
import slib.graph.io.util.GFormat
import slib.graph.io.loader.GraphLoaderGeneric
import slib.graph.model.impl.repo.URIFactoryMemory
import slib.graph.model.impl.graph.elements.Edge
import slib.graph.algo.utils.*

import org.openrdf.model.vocabulary.RDF

import groovyx.gpars.GParsPool

import groovy.cli.commons.CliBuilder
import java.nio.file.Paths

import java.util.logging.*

Logger log = Logger.getLogger("semsim_baseline")
ConsoleHandler ch = new ConsoleHandler()
ch.setLevel(Level.ALL)
ch.setFormatter(new SimpleFormatter())
log.addHandler(ch)
log.setLevel(Level.ALL)
log.setUseParentHandlers(false)

// ── CLI ─────────────────────────────────────────────────────────────────────
def cli = new CliBuilder(usage: 'semsim_baseline.groovy [options]')
cli.d(longOpt: 'data-dir',  args: 1, required: true,  'Path to vexomiser/data/')
cli.u(longOpt: 'upheno',    args: 1, required: true,  'Path to upheno.owl')
cli.ic(longOpt: 'ic',       args: 1, defaultValue: 'resnik', 'IC measure (resnik)')
cli.pw(longOpt: 'pw',       args: 1, defaultValue: 'resnik', 'Pairwise measure (resnik|lin)')
cli.gw(longOpt: 'gw',       args: 1, defaultValue: 'bma',    'Groupwise measure (bma|bmm)')
cli.s(longOpt:  'split',    args: 1, defaultValue: 'test',   'Split (test|val|train)')
cli.o(longOpt:  'out',      args: 1, required: true,  'Output TSV file')

def opts = cli.parse(args)
if (!opts) return

String dataDir   = opts.d
String uphenoOwl = opts.u
String icMeasure = opts.ic
String pwMeasure = opts.pw
String gwMeasure = opts.gw
String split     = opts.s
String outPath   = opts.o

log.info("ic=${icMeasure} pw=${pwMeasure} gw=${gwMeasure} split=${split}")

// ── 1. Load UPheno OWL and extract HP term URIs ─────────────────────────────
log.info("Loading UPheno OWL: ${uphenoOwl}")
def manager  = OWLManager.createOWLOntologyManager()
def ontology = manager.loadOntologyFromOntologyDocument(new File(uphenoOwl))
def classes  = ontology.getClassesInSignature().collect { it.toStringID() }

def existingHpUris = new HashSet<String>()
classes.each { cls ->
    if (cls.contains("HP_")) existingHpUris.add(cls)
}
log.info("HP terms in UPheno: ${existingHpUris.size()}")

def hpToUri = { String hp -> "http://purl.obolibrary.org/obo/" + hp.trim().replace(":", "_") }

// ── 2. Build slib graph from UPheno OWL ─────────────────────────────────────
def factory  = URIFactoryMemory.getSingleton()
def graphUri = factory.getURI("http://purl.obolibrary.org/obo/GDA_")
factory.loadNamespacePrefix("GDA", graphUri.toString())
def graph    = new GraphMemory(graphUri)

def goConf = new GDataConf(GFormat.RDF_XML, Paths.get(uphenoOwl).toString())
GraphLoaderGeneric.populate(goConf, graph)

def virtualRoot = factory.getURI("http://purl.obolibrary.org/obo/GDA_virtual_root")
def rooting = new GAction(GActionType.REROOTING)
rooting.addParameter("root_uri", virtualRoot.stringValue())
GraphActionExecutor.applyAction(factory, rooting, graph)

// ── 3. Eval gene pool (all unique genes in dataset, sorted) ─────────────────
log.info("Loading eval gene pool")
def evalGenes = new File("${dataDir}/track1_cases.tsv").readLines().tail()
    .collect { it.split('\t')[5] }
    .unique()
    .sort()
log.info("Eval genes: ${evalGenes.size()}")

// ── 4. Build gene → HP map from training split ──────────────────────────────
log.info("Building gene phenotype map from training split")
def gene2hpo = new HashMap<String, Set<String>>()
new File("${dataDir}/splits/train.tsv").readLines().tail().each { line ->
    def parts = line.split('\t')
    if (parts.length < 6) return
    def gene   = parts[5]
    def hpoStr = parts[4]
    hpoStr.split(';').each { entry ->
        def uri = hpToUri(entry.split('\\|')[0])
        if (existingHpUris.contains(uri)) {
            gene2hpo.computeIfAbsent(gene, { new HashSet<>() }).add(uri)
        }
    }
}
log.info("Genes with HP annotations: ${gene2hpo.size()} / ${evalGenes.size()}")

// ── 5. Add gene annotations to graph for corpus-based IC ───────────────────
gene2hpo.each { gene, hpUris ->
    def geneUri = factory.getURI("http://vexomiser.org/gene/" + gene.replaceAll('[^A-Za-z0-9]', '_'))
    hpUris.each { hp ->
        graph.addE(new Edge(geneUri, RDF.TYPE, factory.getURI(hp)))
    }
}

// ── 6. Build SM_Engine ──────────────────────────────────────────────────────
log.info("Building SM_Engine")
def engine   = new SM_Engine(graph)
def icConf   = new IC_Conf_Corpus(icMeasureResolver(icMeasure))
def smConfPW = new SMconf(pwMeasureResolver(pwMeasure))
smConfPW.setICconf(icConf)
def smConfGW = new SMconf(gwMeasureResolver(gwMeasure))

// ── 7. Load test cases ──────────────────────────────────────────────────────
log.info("Loading ${split} split")
def testCases = new File("${dataDir}/splits/${split}.tsv").readLines().tail()
log.info("Cases to score: ${testCases.size()}")

// ── 8. Score in parallel ────────────────────────────────────────────────────
log.info("Scoring…")
def allResults = GParsPool.withPool {
    testCases.collectParallel { line ->
        try {
            def parts      = line.split('\t')
            def caseId     = parts[0]
            def causalGene = parts[5]
            def hpoStr     = parts[4]

            def patientHPs = hpoStr.split(';')
                .collect  { hpToUri(it.split('\\|')[0]) }
                .findAll  { existingHpUris.contains(it) }
                .collect  { factory.getURI(it) }
                .toSet()

            def geneIdx = evalGenes.indexOf(causalGene)

            def scores = evalGenes.collect { g ->
                def geneHPs = (gene2hpo.get(g) ?: [] as Set)
                    .collect { factory.getURI(it) }.toSet()
                if (geneHPs.isEmpty() || patientHPs.isEmpty()) return 0.0d
                (double) engine.compare(smConfGW, smConfPW, geneHPs, patientHPs)
            }

            [causalGene, caseId, geneIdx, scores]
        } catch (Exception e) {
            log.warning("Error: ${e.message}")
            null
        }
    }.findAll { it != null }
}

// ── 9. Write output ─────────────────────────────────────────────────────────
log.info("Writing ${allResults.size()} results to ${outPath}")
new File(outPath).parentFile?.mkdirs()
new File(outPath).withWriter { w ->
    allResults.each { r ->
        w.write("${r[0]}\t${r[1]}\t${r[2]}\t${r[3].join('\t')}\n")
    }
}
log.info("Done.")

// ── Measure resolvers ────────────────────────────────────────────────────────
static icMeasureResolver(String m) {
    switch (m.toLowerCase()) {
        case 'resnik':  return SMConstants.FLAG_IC_ANNOT_RESNIK_1995_NORMALIZED
        case 'sanchez': return SMConstants.FLAG_ICI_SANCHEZ_2011
        default: throw new IllegalArgumentException("Unknown IC measure: ${m}")
    }
}

static pwMeasureResolver(String m) {
    switch (m.toLowerCase()) {
        case 'resnik': return SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_RESNIK_1995
        case 'lin':    return SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_LIN_1998
        default: throw new IllegalArgumentException("Unknown pairwise measure: ${m}")
    }
}

static gwMeasureResolver(String m) {
    switch (m.toLowerCase()) {
        case 'bma': return SMConstants.FLAG_SIM_GROUPWISE_BMA
        case 'bmm': return SMConstants.FLAG_SIM_GROUPWISE_BMM
        default: throw new IllegalArgumentException("Unknown groupwise measure: ${m}")
    }
}
