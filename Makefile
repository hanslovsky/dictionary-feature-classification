# dist:    recursive

MCC = /misc/local/matlab-2013b/bin/mcc

dist:
	$(MCC)  -o runObjMethod -W main:runObjMethod -T link:exe -d bin/ -w enable:specified_file_mismatch -w enable:repeated_file -w enable:switch_ignored -w enable:missing_lib_sentinel -w enable:demo_license -R -nodesktop -R -singleCompThread -v ~/dev/dawmr/dawmr_lib/janelia_common/utils/dfeval_gbh_worker.m -a ~/dev/dawmr/dawmr_lib/janelia_common/utils -a ./startup.m -a ../hhmi-exp/matlab/imgproc -a ../hhmi-exp/matlab/dictionary -a ../hhmi-exp/matlab/patches -a ../hhmi-exp/matlab/util -a matlab/dictionary -a ~/dev/dawmr/dawmr_lib/janelia_common -a /groups/saalfeld/home/bogovicj/libraries/patchSearchStandalone_20141021.jar

