from inferelator import workflow
from inferelator.distributed.inferelator_mp import MPControl
from inferelator.crossvalidation_workflow import CrossValidationManager

INPUT_DIR = "/scratch/cjs9968/project/dan_danr/clusterin_atacseq"
OUTPUT_DIR = "/scratch/cjs9968/project/dan_danr/clusterin_atacseq/try"
TF_NAMES = "FlyBaseTFs.txt"
PRIOR_FILE = "inferelator_prior.tsv_unfiltered_matrix.tsv.gz"
EXPRESSION_FILE = "expression_matrix_flybase.tsv"
GOLD_STANDARD = None

MPControl.set_multiprocess_engine("local")
MPControl.connect()

cv_wrap = CrossValidationManager()
cv_wrap.add_gridsearch_parameter('random_seed', [42])

# Initialize workflow
worker = workflow.inferelator_workflow(regression="bbsr", workflow="tfa")

# Set file paths
worker.set_file_paths(
    input_dir=INPUT_DIR,
    output_dir=OUTPUT_DIR,
    tf_names_file=TF_NAMES,
    priors_file=PRIOR_FILE,
    expression_matrix_file=EXPRESSION_FILE,
    gold_standard_file=GOLD_STANDARD
)

worker.set_file_properties(expression_matrix_columns_are_genes=True)
worker.set_run_parameters(num_bootstraps=5)

# Configure priors and network flags
worker.set_network_data_flags(
    use_no_gold_standard=True,
    use_no_prior=False  # Set to True if you want to ignore priors
)

cv_wrap.workflow = worker
cv_wrap.run()

