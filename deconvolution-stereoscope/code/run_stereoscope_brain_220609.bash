# Run Stereoscope with bulk samples

# do before running bash script instead?
conda activate scanpy_scvi

############################################################
# folders
OUTDIR="../../../data/deconv/brain/results"
INDIR="../../../data/deconv/brain/inputs"
SC_DATA="$INDIR/snRNAseq_subsampled_220609.h5ad"
#SC_LIST="$INDIR/top50_degs_220609.csv"
#SC_LIST="$INDIR/top100_degs_220609.csv"
SC_LIST="$INDIR/top20_degs_220609.csv"
#SC_LIST="$INDIR/hvg2000_220609.csv"

# clustering column
CELLTYPE="wsnn_cc_res.0.2"

# Epochs for training/testing, check the plots training_epochs.png in each folder for saturation and decide.
sc_epochs=500
st_epochs=5000

############################################################
# create directories
mkdir -p $OUTDIR

SC_NAME="$(basename $SC_LIST .csv)"
OUTDIR2="$OUTDIR/$SC_NAME"
mkdir -p $OUTDIR2


############################################################
# create SC ref
outfile_ref="$OUTDIR2/sc_ref/model/attr.pkl"
echo $oufile_ref
if test -f "$outfile_ref"; then
    echo "SC ref already created..."
else
    echo "Creating SC ref to $OUTDIR2/sc_ref"
    python stereoscope_create_sc_ref.py -i $SC_DATA -o $OUTDIR2/sc_ref -a $CELLTYPE -g $SC_LIST -e $sc_epochs
fi

############################################################
# run deconv
FILES="$INDIR/*h5ad"
for f in $FILES
do
    if [ "$f" = "$SC_DATA" ]; then
	continue
    fi
    sample="$(basename $f .h5ad)"
    outdir="$OUTDIR2/$sample"
    mkdir -p $outdir
    outfile="$outdir/proportions.csv"
    if test -f "$outfile"; then
	echo "Already done $outfile..."
    else
	echo "Processing $f to $outdir..."
	python stereoscope_run_deconv.py -i $f -o $outdir -r $OUTDIR2/sc_ref -e $st_epochs
    fi
    
done

