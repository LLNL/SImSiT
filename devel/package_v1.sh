#!/bin/tcsh
#
# package_v1
#
#  Package the X-Files v1.0 files for sharing
#
## Check usage
if ($#argv != 0) then
    echo "Usage: packagev1.sh"
    echo "Package the X-Files v1.0 files for sharing"
    goto done
endif

# rm -rf xfiles_v1.tar.gz

# mkdir -p xfiles_v1
# mkdir -p xfiles_v1/detect/sidereal/detect_sidereal_images
# mkdir -p xfiles_v1/detect/target/detect_target_images
# mkdir -p xfiles_v1/iod/target/iod_catalogs

# cp README.md xfiles_v1
# # Detect / Target
# cp detect/target/check_submission_detect_target.py xfiles_v1/detect/target
# cp detect/target/README_detect_target.md xfiles_v1/detect/target
# cp detect/target/score_submission_detect_target.py xfiles_v1/detect/target
# cp detect/target/public/sample_submission_5_detect_target.yaml xfiles_v1/detect/target
# cp detect/target/public/truth_5_detect_target.fits xfiles_v1/detect/target
# # Detect / Sidereal
# cp detect/sidereal/check_submission_detect_sidereal.py xfiles_v1/detect/sidereal
# cp detect/sidereal/README_detect_sidereal.md xfiles_v1/detect/sidereal
# cp detect/sidereal/score_submission_detect_sidereal.py xfiles_v1/detect/sidereal
# cp detect/sidereal/public/sample_submission_5_detect_sidereal.yaml xfiles_v1/detect/sidereal
# cp detect/sidereal/public/truth_5_detect_sidereal.fits xfiles_v1/detect/sidereal
# # IOD / Target
# cp iod/target/check_submission_iod_target.py xfiles_v1/iod/target
# cp iod/target/README_iod_target.md xfiles_v1/iod/target
# cp iod/target/score_submission_iod_target.py xfiles_v1/iod/target
# cp iod/target/public/sample_submission_5_iod_target.yaml xfiles_v1/iod/target
# cp iod/target/public/truth_5_iod_target.fits xfiles_v1/iod/target

# #
# # Copy the test data files
# #
# # cp -r detect/target/public/*.fits xfiles_v1/detect/target/detect_target_images/
# # cp -r detect/sidereal/public/*.fits xfiles_v1/detect/sidereal/detect_sidereal_images/
# cp -r iod/target/public/iod_datasets/* xfiles_v1/iod/target/iod_catalogs/

# tar -czvf xfiles_v1.tar.gz xfiles_v1

# tar -czvf xfiles_v1_detect_target_0.tar.gz detect/target/public/0*.fits
# tar -czvf xfiles_v1_detect_target_1.tar.gz detect/target/public/1*.fits
# tar -czvf xfiles_v1_detect_target_2.tar.gz detect/target/public/2*.fits
# tar -czvf xfiles_v1_detect_target_3.tar.gz detect/target/public/3*.fits
# tar -czvf xfiles_v1_detect_target_4.tar.gz detect/target/public/4*.fits

tar -czvf xfiles_v1_detect_sidereal_0.tar.gz detect/sidereal/public/0*.fits
tar -czvf xfiles_v1_detect_sidereal_1.tar.gz detect/sidereal/public/1*.fits
tar -czvf xfiles_v1_detect_sidereal_2.tar.gz detect/sidereal/public/2*.fits
tar -czvf xfiles_v1_detect_sidereal_3.tar.gz detect/sidereal/public/3*.fits
tar -czvf xfiles_v1_detect_sidereal_4.tar.gz detect/sidereal/public/4*.fits

###
# cp mbi-theory_arxiv.tar.gz ~/tmp
# cd ~/tmp
# tar -zxvf mbi-theory_arxiv.tar.gz
# pdflatex mbi-theory.tex

## labels for proper exit or error
done:
	echo "Created xfiles_v1.tar.gz"
    exit 0
error:
    exit 1
