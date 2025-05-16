#!/bin/sh

#-----------
version=3.0
#-----------

data_dir=xfiles_v${version}
cd /p/lustre2/pruett6/xfiles/xfiles_sim/

echo "======================="
echo "Packaging '${data_dir}'"
echo "======================="

if [ -d xfiles_v${version} ]
then
    rm -rf xfiles_v${version}
fi

if [ -d xfiles_v* ]
then
    # Move last version to devel if it exists
    echo "Moving previous data package to /devel"
    mv xfiles_v*/ devel/
fi

# Make new version directory
mkdir -p ${data_dir}

# Copy the readme and changelog over
cp change_log.md ${data_dir}/

#-----------
# IOD branch
#-----------
echo "*** TRANSFERRING IOD BRANCH ***"
mkdir -p ${data_dir}/iod/iod_datasets
rsync --info=progress2 branches/iod/public/iod_datasets/*.csv ${data_dir}/iod/iod_datasets/
wait
cp branches/iod/public/truth_5.fits ${data_dir}/iod/
cp branches/iod/public/sample_submission*.yaml ${data_dir}/iod/
cp branches/iod/score_submission.py ${data_dir}/iod/
cp branches/iod/README.md ${data_dir}/iod/
echo "*** ZIPPING IOD BRANCH ***"
cd /p/lustre2/pruett6/xfiles/xfiles_sim/${data_dir}/iod/
# zip -qr ../${data_dir}_web/iod.zip .
7z a ../iod.zip .
cd /p/lustre2/pruett6/xfiles/xfiles_sim/

#----------------
# Sidereal Branch
#----------------
echo "*** TRANSFERRING SIDEREAL BRANCH ***"
mkdir -p ${data_dir}/sidereal/images
rsync --info=progress2 branches/tracking/sidereal/public/*.fits ${data_dir}/sidereal/images/
wait
mv ${data_dir}/sidereal/images/sky_flat.fits ${data_dir}/sidereal/
mv ${data_dir}/sidereal/images/truth_10.fits ${data_dir}/sidereal/
cp branches/tracking/sidereal/public/sample_submission_10.yaml ${data_dir}/sidereal/
cp branches/tracking/README.md ${data_dir}/sidereal/
cp branches/tracking/score_submission.py ${data_dir}/sidereal/
echo "*** ZIPPING SIDEREAL BRANCH ***"
cd /p/lustre2/pruett6/xfiles/xfiles_sim/${data_dir}/sidereal/
# zip -qr ../${data_dir}_web/sidereal.zip .
7z a ../sidereal.zip .
cd /p/lustre2/pruett6/xfiles/xfiles_sim/

#--------------
# Target Branch
#--------------
echo "*** TRANSFERRING TARGET BRANCH ***"
mkdir -p ${data_dir}/target/images
rsync --info=progress2 branches/tracking/target/public/*.fits ${data_dir}/target/images/
wait
mv ${data_dir}/target/images/sky_flat.fits ${data_dir}/target/
mv ${data_dir}/target/images/truth_10.fits ${data_dir}/target/
cp branches/tracking/target/public/sample_submission_10.yaml ${data_dir}/target/
cp branches/tracking/README.md ${data_dir}/target/
cp branches/tracking/score_submission.py ${data_dir}/target/
echo "*** ZIPPING TARGET BRANCH ***"
cd /p/lustre2/pruett6/xfiles/xfiles_sim/${data_dir}/target/
# zip -qr ../${data_dir}_web/target.zip .
7z a ../target.zip .
cd /p/lustre2/pruett6/xfiles/xfiles_sim/

#----------------
# Sidereal Track Images Branch
#----------------
echo "*** TRANSFERRING SIDEREAL TRACK IMAGES BRANCH ***"
mkdir -p ${data_dir}/sidereal_track/images
rsync --info=progress2 branches/track_images/sidereal_track/public/*.fits ${data_dir}/sidereal_track/images/
wait
mv ${data_dir}/sidereal_track/images/sky_flat.fits ${data_dir}/sidereal_track/
mv ${data_dir}/sidereal_track/images/truth_5.fits ${data_dir}/sidereal_track/
cp branches/track_images/sidereal_track/public/sample_submission_5.yaml ${data_dir}/sidereal_track/
cp branches/track_images/README.md ${data_dir}/sidereal_track/
cp branches/track_images/score_submission.py ${data_dir}/sidereal_track/
echo "*** ZIPPING SIDEREAL TRACK IMAGES BRANCH ***"
cd /p/lustre2/pruett6/xfiles/xfiles_sim/${data_dir}/sidereal_track/
# zip -qr ../${data_dir}_web/sidereal_track.zip .
7z a ../sidereal_track.zip .
cd /p/lustre2/pruett6/xfiles/xfiles_sim/

#--------------
# Target Track Images Branch
#--------------
echo "*** TRANSFERRING TARGET TRACK IMAGES BRANCH ***"
mkdir -p ${data_dir}/target_track/images
rsync --info=progress2 branches/track_images/target_track/public/*.fits ${data_dir}/target_track/images/
wait
mv ${data_dir}/target_track/images/sky_flat.fits ${data_dir}/target_track/
mv ${data_dir}/target_track/images/truth_5.fits ${data_dir}/target_track/
cp branches/track_images/target_track/public/sample_submission_5.yaml ${data_dir}/target_track/
cp branches/track_images/README.md ${data_dir}/target_track/
cp branches/track_images/score_submission.py ${data_dir}/target_track/
echo "*** ZIPPING TARGET TRACK IMAGES BRANCH ***"
cd /p/lustre2/pruett6/xfiles/xfiles_sim/${data_dir}/target_track/
# zip -qr ../${data_dir}_web/target_track.zip .
7z a ../target_track.zip .
cd /p/lustre2/pruett6/xfiles/xfiles_sim/

echo "*** REMOVING FILES ***"
rm -rf ${data_dir}/target
rm -rf ${data_dir}/sidereal
rm -rf ${data_dir}/iod
rm -rf ${data_dir}/target_track
rm -rf ${data_dir}/sidereal_track

echo "${data_dir} package is complete!" 

