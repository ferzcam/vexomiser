/**
 * Semantic similarity baseline (Resnik/Lin + BMA) on the PAVS benchmark.
 *
 * Adapted from multihop-gda/semantic_similarity.groovy for the vexomiser eval setup:
 *   - Loads HP OBO (not UPheno)
 *   - Gene phenotypes come from training-case HPO terms (not MGI MP)
 *   - "Disease" phenotypes are the patient HPO terms from each test case directly
 *   - Output TSV matches the format expected by eval/metrics.py
 *
 * Usage:
 *   groovy eval/semsim_baseline.groovy \
 *     --data-dir data \
 *     --hp-obo exomiser-data/2406_phenotype/hp.obo \
 *     --ic resnik --pw resnik --gw bma \
 *     --split test \
 *     --out data/results/semsim_resnik_bma_test.tsv
 */

@Grab(group='com.github.sharispe', module='slib-sml',         version='0.9.1')
@Grab(group='net.sourceforge.owlapi', module='owlapi-api',     version='4.2.5')
@Grab(group='net.sourceforge.owlapi', module='owlapi-apibinding', version='4.2.5')
@Grab(group='net.sourceforge.owlapi', module='owlapi-impl',    version='4.2.5')
@Grab(group='ch.qos.logback',  module='logback-classic',       version='1.2.3')
@Grab(group='org.slf4j',       module='slf4j-api',             version='1.7.30')

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

import groovy.cli.commons.CliBuilder

import java.util.concurrent.*
import java.util.logging.*

// ── Logger ─────────────────────────────────────────────────────────────────
Logger log = Logger.getLogger("semsim_baseline")
ConsoleHandler ch = new ConsoleHandler()
ch.setLevel(Level.INFO)
ch.setFormatter(new SimpleFormatter())
log.addHandler(ch)
log.setLevel(Level.INFO)
log.setUseParentHandlers(false)

// ── CLI ─────────────────────────────────────────────────────────────────────
def cli = new CliBuilder(usage: 'semsim_baseline.groovy [options]')
cli.d(longOpt: 'data-dir', args: 1, required: true,  'Path to vexomiser/data/')
cli.hp(longOpt: 'hp-obo',  args: 1, required: true,  'Path to hp.obo')
cli.ic(longOpt: 'ic',      args: 1, defaultValue: 'resnik', 'IC measure (resnik)')
cli.pw(longOpt: 'pw',      args: 1, defaultValue: 'resnik', 'Pairwise measure (resnik|lin)')
cli.gw(longOpt: 'gw',      args: 1, defaultValue: 'bma',    'Groupwise measure (bma|bmm)')
cli.s(longOpt:  'split',   args: 1, defaultValue: 'test',   'Split (test|val|train)')
cli.o(longOpt:  'out',     args: 1, required: true,  'Output TSV file')
cli.t(longOpt:  'threads', args: 1, defaultValue: '0', 'Worker threads (0 = #CPUs)')

def opts = cli.parse(args)
if (!opts) return

String dataDir   = opts.d
String hpObo     = opts.hp
String icMeasure = opts.ic
String pwMeasure = opts.pw
String gwMeasure = opts.gw
String split     = opts.s
String outPath   = opts.o
int    nThreads  = opts.t.toInteger() ?: Runtime.getRuntime().availableProcessors()

log.info("ic=${icMeasure} pw=${pwMeasure} gw=${gwMeasure} split=${split} threads=${nThreads}")

// ── 1. Load HP OBO ──────────────────────────────────────────────────────────
log.info("Loading HP ontology: ${hpObo}")
def factory  = URIFactoryMemory.getSingleton()
def graphUri = factory.getURI("http://purl.obolibrary.org/obo/")
def graph    = new GraphMemory(graphUri)
GraphLoaderGeneric.populate(new GDataConf(GFormat.OBO, hpObo), graph)

// Add a virtual root (required by slib for IC computation)
def virtualRoot = factory.getURI("http://vexomiser.org/semsim_virtual_root")
def rooting = new GAction(GActionType.REROOTING)
rooting.addParameter("root_uri", virtualRoot.stringValue())
GraphActionExecutor.applyAction(factory, rooting, graph)

def knownUris = graph.getV().collect { it.stringValue() }.toSet()
log.info("Ontology nodes after loading: ${knownUris.size()}")

// HP:0001263 → http://purl.obolibrary.org/obo/HP_0001263
def hpToUri = { String hp -> "http://purl.obolibrary.org/obo/" + hp.trim().replace(":", "_") }

// ── 2. Eval gene pool (all unique genes in dataset, sorted) ─────────────────
log.info("Loading eval gene pool from track1_cases.tsv")
def evalGenes = new File("${dataDir}/track1_cases.tsv").readLines().tail()
    .collect { it.split('\t')[5] }
    .unique()
    .sort()
log.info("Eval genes: ${evalGenes.size()}")

// ── 3. Build gene→HP map from training split only ───────────────────────────
log.info("Building gene phenotype map from training split")
def gene2hpo = new HashMap<String, Set<String>>()
new File("${dataDir}/splits/train.tsv").readLines().tail().each { line ->
    def parts = line.split('\t')
    if (parts.length < 6) return
    def gene   = parts[5]
    def hpoStr = parts[4]
    hpoStr.split(';').each { entry ->
        def hp  = entry.split('\\|')[0].trim()
        def uri = hpToUri(hp)
        if (knownUris.contains(uri)) {
            gene2hpo.computeIfAbsent(gene, { new HashSet<>() }).add(uri)
        }
    }
}
log.info("Genes with HP annotations: ${gene2hpo.size()} / ${evalGenes.size()}")

// ── 4. Add gene annotations to graph (for corpus-based IC) ─────────────────
gene2hpo.each { gene, hpUris ->
    def geneUri = factory.getURI("http://vexomiser.org/gene/" + gene.replaceAll('[^A-Za-z0-9]', '_'))
    hpUris.each { hp ->
        graph.addE(new Edge(geneUri, RDF.TYPE, factory.getURI(hp)))
    }
}

// ── 5. Build SM_Engine ──────────────────────────────────────────────────────
log.info("Building SM_Engine")
def engine = new SM_Engine(graph)

def icConf = new IC_Conf_Corpus(icMeasureResolver(icMeasure))
def smConfPW = new SMconf(pwMeasureResolver(pwMeasure))
smConfPW.setICconf(icConf)
def smConfGW = new SMconf(gwMeasureResolver(gwMeasure))

// ── 6. Load test cases ──────────────────────────────────────────────────────
log.info("Loading ${split} split")
def testCases = new File("${dataDir}/splits/${split}.tsv").readLines().tail()
log.info("Cases to score: ${testCases.size()}")

// ── 7. Score in parallel ────────────────────────────────────────────────────
log.info("Scoring with ${nThreads} threads…")
def pool    = Executors.newFixedThreadPool(nThreads)
def futures = testCases.collect { line ->
    pool.submit({
        try {
            def parts      = line.split('\t')
            def caseId     = parts[0]
            def causalGene = parts[5]
            def hpoStr     = parts[4]

            def patientHPs = hpoStr.split(';')
                .collect { hpToUri(it.split('\\|')[0]) }
                .findAll  { knownUris.contains(it) }
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
    } as Callable)
}

def allResults = futures.collect { it.get() }.findAll { it != null }
pool.shutdown()

// ── 8. Write output ─────────────────────────────────────────────────────────
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
