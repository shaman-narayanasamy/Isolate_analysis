## Input and output parameters
#

PWD = os.getcwd()

SAMPLE = os.environ.get("SAMPLE")
INPUTDIR = os.environ.get("INPUTDIR")
OUTPUTDIR = os.environ.get("OUTPUTDIR")
TMPDIR = os.environ.get("TMPDIR", "/tmp")
SRCDIR = os.environ.get("SRCDIR", "%s/src" % PWD)
CONFIG = os.environ.get("CONFIG", "conf/config_normalNode.json")
DBPATH = os.environ.get("DBPATH", "/mnt/nfs/projects/ecosystem_biology/local_tools/IMP/dependencies/prokka/db")

configfile: CONFIG

MEMCORE = os.environ.get("MEMCORE", config['memory_per_core_gb'])
THREADS = os.environ.get("THREADS", config['threads'])
MEMTOTAL = os.environ.get("MEMTOTAL", config['memory_total_gb'])
FILTER = os.environ.get("FILTER", "hg38")

workdir:
    OUTPUTDIR

include:
    "workflows/Preprocessing"

include:
    "workflows/Assembly"

include:
    "workflows/Analysis"

rule ALL:
    input:
        'Preprocessing/preprocessing.done',
        'Assembly/assembly.done',
        'Analysis/analysis.done'
    output:
        touch('new_assembly.done')
