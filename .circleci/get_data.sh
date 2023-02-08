
if [[ "$SHELL" == zsh ]]; then
  setopt SH_WORD_SPLIT
fi

# Edit these for project-wide testing
WGET="wget --retry-connrefused --waitretry=5 --read-timeout=20 --timeout=15 -t 0 -q"
IMAGE="pennlinc/aslprep:unstable"

# Determine if we're in a CI test
if [[ "${CIRCLECI}" = "true" ]]; then
  IN_CI="true"
  NTHREADS=2
  OMP_NTHREADS=2

  if [[ -n "${CIRCLE_CPUS}" ]]; then
    NTHREADS=${CIRCLE_CPUS}
    OMP_NTHREADS=$(expr $NTHREADS - 1)
  fi

else
  IN_CI="false"
  NTHREADS=2
  OMP_NTHREADS=2

  LOCAL_PATCH_FILE="local_xcpd_path.txt"

  # check that the patch file exists
  if [ ! -f $LOCAL_PATCH_FILE ]
  then
    echo "File $LOCAL_PATCH_FILE DNE"
    exit 1
  fi

  LOCAL_PATCH="$( cat ${LOCAL_PATCH_FILE} )"  # Load path from file

  # check that the local xcp_d path exists
  if [ ! -d $LOCAL_PATCH ]
  then
    echo "Path $LOCAL_PATCH DNE"
    exit 1
  fi

fi
export IN_CI NTHREADS OMP_NTHREADS


get_bids_data() {
    WORKDIR=$1
    DS=$2
    echo "working dir: ${WORKDIR}"
    echo "fetching dataset: ${DS}"
    ENTRYDIR=`pwd`
    TEST_DATA_DIR="${WORKDIR}/data"
    mkdir -p $TEST_DATA_DIR
    cd $TEST_DATA_DIR

    # without freesurfer, sub-01
    if [[ ${DS} = downsampled ]]
    then
      dataset_dir="$TEST_DATA_DIR/$DS"
      # Do not re-download if the folder exists
      if [ ! -d $dataset_dir ]
      then
        echo "Downloading ${DS} data to $dataset_dir"

        # Download the raw BIDS dataset
        ${WGET} \
          -O downsampled.tar.xz \
        "https://upenn.box.com/shared/static/og1ixccv5v8eir76emii6rrgnwu4thad.xz"
        tar xvfJ downsampled.tar.xz -C $TEST_DATA_DIR
        mkdir dset
        mv testingbids dset
        rm downsampled.tar.xz

        # Download the pre-generated smriprep derivatives
        ${WGET} \
          -O smriprepx.tar.xz \
        "https://upenn.box.com/shared/static/i64rbrpzfinpej0vct96mjw6ve7yycov.xz"
        tar xvfJ smriprepx.tar.xz -C $TEST_DATA_DIR
        mkdir dset/derivatives
        mv smriprep dset/derivatives/smriprep
        rm smriprepx.tar.xz

      else
        echo "Data directory ($dataset_dir) already exists. If you need to re-download the data, remove the data folder."
      fi
    fi

    cd ${ENTRYDIR}
}
