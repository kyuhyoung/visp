#   usage : sh sh_handeye_calib_and_show_poses.sh 0.027 image-d%d.png 8 pose_fPe_ pose_cPo_ eMc
if [ "$#" -ge 6 ]; then
METER_SQUARE_SIZE=$1
FORMAT_IMG_PATH=$2
N_IMG=$3
PREFIX_BASE2GRIPPER=$4
PREFIX_CAM2BOARD=$5
FN_GRIPPER2CAMERA=$6
./tutorial-chessboard-pose --square_size $METER_SQUARE_SIZE --input $FORMAT_IMG_PATH
./tutorial-hand-eye-calibration --ndata $N_IMG --prefix_yaml_base2gripper $PREFIX_BASE2GRIPPER --prefix_yaml_camera2board $PREFIX_CAM2BOARD --fn_gripper2camera $FN_GRIPPER2CAMERA 
#python hand_eye_calibration_show_extrinsics_with_robot_base.py --ndata $N_IMG --eMc_yaml $FN_GRIPPER2CAMERA.yaml
python hand_eye_calibration_show_extrinsics_with_robot_base.py --ndata $N_IMG --fPe_file_pattern $PREFIX_BASE2GRIPPER%d.yaml --cPo_file_pattern $PREFIX_CAM2BOARD%d.yaml --eMc_yaml $FN_GRIPPER2CAMERA.yaml
elif [ "$#" -ge 4 ]; then
N_IMG=$1
PREFIX_BASE2GRIPPER=$2
PREFIX_CAM2BOARD=$3
FN_GRIPPER2CAMERA=$4
./tutorial-hand-eye-calibration --ndata $N_IMG --prefix_yaml_base2gripper $PREFIX_BASE2GRIPPER --prefix_yaml_camera2board $PREFIX_CAM2BOARD --fn_gripper2camera $FN_GRIPPER2CAMERA 
python hand_eye_calibration_show_extrinsics_with_robot_base.py --ndata $N_IMG --fPe_file_pattern $PREFIX_BASE2GRIPPER%d.yaml --cPo_file_pattern $PREFIX_CAM2BOARD%d.yaml --eMc_yaml $FN_GRIPPER2CAMERA.yaml
else
N_IMG=$1
FN_GRIPPER2CAMERA=$2
python hand_eye_calibration_show_extrinsics_with_robot_base.py --ndata $N_IMG --eMc_yaml $FN_GRIPPER2CAMERA.yaml
fi
