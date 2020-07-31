#   usage : sh sh_handeye_calib_and_show_poses.sh 0.027 image-d%d.png 8 pose_fPe_ pose_cPo_ eMc
if [ "$#" -eq 7 ]; then
echo "# of arguments is 7"
METER_SQUARE_SIZE=$1
FORMAT_IMG_PATH=$2
ID_START=$3
N_IMG=$4
PREFIX_BASE2GRIPPER=$5
PREFIX_CAM2BOARD=$6
FN_GRIPPER2CAMERA=$7
./tutorial-chessboard-pose --square_size $METER_SQUARE_SIZE --input $FORMAT_IMG_PATH
./tutorial-hand-eye-calibration --ndata $N_IMG --prefix_yaml_base2gripper $PREFIX_BASE2GRIPPER --prefix_yaml_camera2board $PREFIX_CAM2BOARD --fn_gripper2camera $FN_GRIPPER2CAMERA 
#python hand_eye_calibration_show_extrinsics_with_robot_base.py --ndata $N_IMG --eMc_yaml $FN_GRIPPER2CAMERA.yaml
python hand_eye_calibration_show_extrinsics_with_robot_base.py --start_index $ID_START --ndata $N_IMG --fPe_file_pattern $PREFIX_BASE2GRIPPER%d.yaml --cPo_file_pattern $PREFIX_CAM2BOARD%d.yaml --eMc_yaml $FN_GRIPPER2CAMERA.yaml
elif [ "$#" -eq 3 ]; then
#   usage : sh sh_handeye_calib_and_show_poses.sh 0 10 floyd_200730_
echo "# of arguments is 3"
ID_START=$1
N_IMG=$2
PREFIX_BASE2GRIPPER_CAM2BOARD=$3
PREFIX_BASE2GRIPPER=$PREFIX_BASE2GRIPPER_CAM2BOARD"fPe_"
PREFIX_CAM2BOARD=$PREFIX_BASE2GRIPPER_CAM2BOARD"cPo_"
FN_GRIPPER2CAMERA=$PREFIX_BASE2GRIPPER_CAM2BOARD"eMc"
#echo "FN_GRIPPER2CAMERA : "$FN_GRIPPER2CAMERA
#exit
./tutorial-hand-eye-calibration --id_start $ID_START --ndata $N_IMG --floyd_txt $PREFIX_BASE2GRIPPER_CAM2BOARD --prefix_yaml_base2gripper $PREFIX_BASE2GRIPPER --prefix_yaml_camera2board $PREFIX_CAM2BOARD --fn_gripper2camera $FN_GRIPPER2CAMERA 
python hand_eye_calibration_show_extrinsics_with_robot_base.py --start_index $ID_START --ndata $N_IMG --fPe_file_pattern $PREFIX_BASE2GRIPPER%d.yaml --cPo_file_pattern $PREFIX_CAM2BOARD%d.yaml --eMc_yaml $FN_GRIPPER2CAMERA.yaml
elif [ "$#" -eq 5 ]; then
echo "# of arguments is 5"
ID_START=$1
N_IMG=$2
PREFIX_BASE2GRIPPER=$3
PREFIX_CAM2BOARD=$4
FN_GRIPPER2CAMERA=$5
./tutorial-hand-eye-calibration --id_start $ID_START --ndata $N_IMG --prefix_yaml_base2gripper $PREFIX_BASE2GRIPPER --prefix_yaml_camera2board $PREFIX_CAM2BOARD --fn_gripper2camera $FN_GRIPPER2CAMERA 
python hand_eye_calibration_show_extrinsics_with_robot_base.py --start_index $ID_START --ndata $N_IMG --fPe_file_pattern $PREFIX_BASE2GRIPPER%d.yaml --cPo_file_pattern $PREFIX_CAM2BOARD%d.yaml --eMc_yaml $FN_GRIPPER2CAMERA.yaml
else
ID_START=$1
N_IMG=$2
FN_GRIPPER2CAMERA=$3
python hand_eye_calibration_show_extrinsics_with_robot_base.py --start_index $ID_START --ndata $N_IMG --eMc_yaml $FN_GRIPPER2CAMERA.yaml
fi
