//! \example tutorial-chessboard-pose.cpp
#include <iostream>

#include <visp3/core/vpConfig.h>

#if VISP_HAVE_OPENCV_VERSION >= 0x020300

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <visp3/gui/vpDisplayX.h>
#include <visp3/gui/vpDisplayGDI.h>
#include <visp3/gui/vpDisplayOpenCV.h>
#include <visp3/gui/vpDisplayD3D.h>
#include <visp3/core/vpIoTools.h>
#include <visp3/core/vpPoint.h>
#include <visp3/core/vpPixelMeterConversion.h>
#include <visp3/core/vpXmlParserCamera.h>
#include <visp3/io/vpVideoReader.h>
#include <visp3/vision/vpPose.h>

using namespace cv;

namespace {
void calcChessboardCorners(int width, int height, double squareSize, std::vector<vpPoint> &corners) {
  corners.resize(0);

  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      vpPoint pt;
      pt.set_oX(j*squareSize);
      pt.set_oY(i*squareSize);
      pt.set_oZ(0.0);
      corners.push_back(pt);
    }
  }
}
} //namespace


// Checks if a matrix is a valid rotation matrix.
bool is_rotation_matrix(const Mat &R)
{
    Mat Rt;
        transpose(R, Rt);
            Mat shouldBeIdentity = Rt * R;
                Mat I = Mat::eye(3,3, shouldBeIdentity.type());
                    
                        return  norm(I, shouldBeIdentity) < 1e-6;
                            
                            }


cv::Mat combine_rotation_translation_into_homogeneous_matrix(const cv::Mat& rot, const cv::Mat& tra)
{
    cv::Mat rot_mat;
    bool is_rot_valid = true;
    int n_elem_rot = rot.rows * rot.cols;
    if(3 == n_elem_rot)
    {   
        cv::Rodrigues(rot, rot_mat);    
    }
    else if(9 == n_elem_rot)
    {
        if(is_rotation_matrix(rot)) rot_mat = rot;
        else is_rot_valid = false;
    }
    else is_rot_valid = false;
    if(!is_rot_valid) { std::cout << "The given rotation matirx is NOT a real rotation matrix." << std::endl; exit(0); }
    cv::Mat mat_homo = cv::Mat::eye(4, 4, CV_32F);
    rot_mat.copyTo(mat_homo(cv::Rect(0, 0, 3, 3)));
    tra.copyTo(mat_homo(cv::Rect(3, 0, 1, 3)));
    return mat_homo;
}

template <typename T> 
Mat init_mat_with_array_of_values(int n_r, int n_c, T li_val[]/*, int n_sp*/)
{
    std::cout << "AAA init" << std::endl;
    Mat_<T> mat(n_r, n_c, li_val);
    std::cout << "BBB init" << std::endl;
    return mat;
#if 0    
    //cout_indented(n_sp, "init_mat START");
        int iR, iE, iC;
            Mat_<T> mat(n_r, n_c);
                for(iR = 0, iE = 0; iR < n_r; iR++)
                    {
                            for(iC = 0; iC < n_c; iC++, iE++)
                                    {
                                                mat(iR, iC) = li_val[iE];
                                                        }      
                                                            }
                                                                //cout_indented(n_sp, "init_mat END");
                                                                    return mat;
#endif //   0                                                                    
                                                                    }

std::pair<cv::Mat, cv::Mat> split_homogeneous_transform_matrix_into_rotation_and_translation(const Mat& mat_homo)
{
    if(!(mat_homo.rows >= 3 && mat_homo.cols >= 3 && mat_homo.cols == mat_homo.rows))
    {
        std::cout << "The given matrix is NOT homogeneous." << std::endl; exit(0);
    }
    int dim_rot = mat_homo.rows - 1;
    cv::Mat mat_rot = mat_homo(cv::Rect(0, 0, dim_rot, dim_rot)).clone();
    cv::Mat vec_tra = mat_homo(cv::Rect(dim_rot, 0, 1, dim_rot)).clone();   
    return std::pair<cv::Mat, cv::Mat>(mat_rot, vec_tra);
}

int main(int argc, const char ** argv) 
{
#if 0
    float data_rad[] = {0.1, 0.2, 0.3}, data_tra[] = {4, 8, 12};
    cv::Mat rot_mat, rot_vec = init_mat_with_array_of_values<float>(3, 1, data_rad), tra_vec = init_mat_with_array_of_values<float>(3, 1, data_tra); 
    std::cout << "rot_vec : " << std::endl << rot_vec << std::endl;   
    std::cout << "tra_vec : " << std::endl << tra_vec << std::endl;   
    cv::Mat mat_homo_1 = combine_rotation_translation_into_homogeneous_matrix(rot_vec, tra_vec);
    std::cout << "mat_homo_1 : " << std::endl << mat_homo_1 << std::endl;    
    cv::Rodrigues(rot_vec, rot_mat);
    std::cout << "rot_mat : " << std::endl << rot_mat << std::endl;   
    cv::Mat mat_homo_2 = combine_rotation_translation_into_homogeneous_matrix(rot_mat, tra_vec);
    std::cout << "mat_homo_2 : " << std::endl << mat_homo_2 << std::endl;
    std::pair<cv::Mat, cv::Mat> mat_rot_and_vec_tra = split_homogeneous_transform_matrix_into_rotation_and_translation(mat_homo_2);
    cv::Mat mat_rot = mat_rot_and_vec_tra.first, vec_tra = mat_rot_and_vec_tra.second;
    std::cout << "mat_rot : " << std::endl << mat_rot << std::endl;
    std::cout << "vec_tra : " << std::endl << vec_tra << std::endl;
    exit(0);
#endif
  int chessboard_width = 9, chessboard_height = 6;
  double chessboard_square_size = 0.03;
  std::string input_filename = "";
  std::string intrinsic_file = "camera.xml";
  std::string camera_name = "Camera";

  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "-w" && i+1 < argc) {
      chessboard_width = atoi(argv[i+1]);
    } else if (std::string(argv[i]) == "-h" && i+1 < argc) {
      chessboard_height = atoi(argv[i+1]);
    } else if (std::string(argv[i]) == "--square_size" && i+1 < argc) {
      chessboard_square_size = atof(argv[i+1]);
    } else if (std::string(argv[i]) == "--input" && i+1 < argc) {
      input_filename = std::string(argv[i+1]);
    } else if (std::string(argv[i]) == "--intrinsic" && i+1 < argc) {
      intrinsic_file = std::string(argv[i+1]);
    } else if (std::string(argv[i]) == "--camera_name" && i+1 < argc) {
      camera_name = std::string(argv[i+1]);
    }
    else if (std::string(argv[i]) == "--help") {
      std::cout << argv[0] << " [-w <chessboard width>] [-w <chessboard height>] [--square_size <square size in meter>] [--input <input images path>] [--intrinsic <Camera intrinsic parameters xml file>] [--camera_name <Camera name in the xml intrinsic file>]" << std::endl;
      return EXIT_SUCCESS;
    }
  }

  std::cout << "Parameters:" << std::endl;
  std::cout << "chessboard_width=" << chessboard_width << std::endl;
  std::cout << "chessboard_height=" << chessboard_height << std::endl;
  std::cout << "chessboard_square_size=" << chessboard_square_size << std::endl;
  std::cout << "input_filename=" << input_filename << std::endl;
  std::cout << "intrinsic_file=" << intrinsic_file << std::endl;
  std::cout << "camera_name=" << camera_name << std::endl;

  vpVideoReader reader;
  if (input_filename.empty()) {
    std::cerr << "input_filename.empty()" << std::endl;
    return EXIT_FAILURE;
  }
  reader.setFileName(input_filename);

  vpImage<vpRGBa> I;
  reader.open(I);
    //std::cout << "AAA chess" << std::endl;
#ifdef VISP_HAVE_X11
    std::cout << "BBB chess" << std::endl;
  vpDisplayX d(I);
#elif defined VISP_HAVE_GDI
    std::cout << "CCC chess" << std::endl;
  vpDisplayGDI d(I);
#elif defined VISP_HAVE_OPENCV
    std::cout << "DDD chess" << std::endl;
  vpDisplayOpenCV d(I);
#endif

  std::vector<vpPoint> corners_pts;
  calcChessboardCorners(chessboard_width, chessboard_height, chessboard_square_size, corners_pts);

  vpCameraParameters cam;
#ifdef VISP_HAVE_PUGIXML
  vpXmlParserCamera parser;
  if (!intrinsic_file.empty() && !camera_name.empty()) {
    parser.parse(cam, intrinsic_file, camera_name, vpCameraParameters::perspectiveProjWithDistortion);
  }
#endif
  std::cout << "cam:\n" << cam << std::endl;

  bool quit = false;
  do {
    reader.acquire(I);

    cv::Mat matImg;
    vpImageConvert::convert(I, matImg);

    vpDisplay::displayText(I, 20, 20, "Right click to quit.", vpColor::red);

    cv::Size chessboardSize(chessboard_width, chessboard_height);
    std::vector<cv::Point2f> corners2D;
    bool found = cv::findChessboardCorners(matImg, chessboardSize, corners2D,
#if (VISP_HAVE_OPENCV_VERSION >= 0x030000)
                                   cv::CALIB_CB_ADAPTIVE_THRESH | cv::CALIB_CB_FAST_CHECK | cv::CALIB_CB_NORMALIZE_IMAGE);
#else
                                   CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_FAST_CHECK | CV_CALIB_CB_NORMALIZE_IMAGE);
#endif

    vpHomogeneousMatrix cMo;
    if (found) {
      cv::Mat matImg_gray;
      cv::cvtColor(matImg, matImg_gray, cv::COLOR_BGR2GRAY);
      cv::cornerSubPix(matImg_gray, corners2D, cv::Size(11,11),
                    cv::Size(-1,-1),
#if (VISP_HAVE_OPENCV_VERSION >= 0x030000)
                    cv::TermCriteria( cv::TermCriteria::EPS+cv::TermCriteria::COUNT, 30, 0.1 ));
#else
                    cv::TermCriteria( CV_TERMCRIT_EPS+CV_TERMCRIT_ITER, 30, 0.1 ));
#endif

      for (size_t i = 0; i < corners_pts.size(); i++) {
        vpImagePoint imPt(corners2D[i].y, corners2D[i].x);
        double x = 0.0, y = 0.0;
        vpPixelMeterConversion::convertPoint(cam, imPt, x, y);
        corners_pts[i].set_x(x);
        corners_pts[i].set_y(y);
      }

      vpPose pose;
      pose.addPoints(corners_pts);
      vpHomogeneousMatrix cMo_dementhon, cMo_lagrange;
      double r_dementhon = std::numeric_limits<double>::max(), r_lagrange = std::numeric_limits<double>::max();
      bool pose_dementhon = pose.computePose(vpPose::DEMENTHON, cMo_dementhon);
      if (pose_dementhon)
        r_dementhon = pose.computeResidual(cMo_dementhon);

      bool pose_lagrange = pose.computePose(vpPose::LAGRANGE, cMo_lagrange);
      if (pose_lagrange)
        r_lagrange = pose.computeResidual(cMo_lagrange);

      cMo = (r_dementhon < r_lagrange) ? cMo_dementhon : cMo_lagrange;
      if (!pose.computePose(vpPose::VIRTUAL_VS, cMo)) {
        std::cerr << "Problem when computing final pose using VVS" << std::endl;
        return EXIT_FAILURE;
      }

      cv::drawChessboardCorners(matImg, chessboardSize, corners2D, found);
      vpImageConvert::convert(matImg, I);
    }

    vpDisplay::display(I);

    vpDisplay::displayText(I, 20, 20, "Left click for the next image, right click to quit.", vpColor::red);
    if (found)
      vpDisplay::displayFrame(I, cMo, cam, 0.05, vpColor::none, 3);

    vpDisplay::flush(I);

    if (found) {
      vpPoseVector pose_vec(cMo);
      std::stringstream ss;
      ss << "pose_cPo_" << reader.getFrameIndex() << ".yaml";
      std::cout << "Save " << ss.str() << std::endl;
      pose_vec.saveYAML(ss.str(), pose_vec);
    }

    vpMouseButton::vpMouseButtonType button;
    if (vpDisplay::getClick(I, button, true)) {
      switch (button) {
        case vpMouseButton::button3:
          quit = true;
          break;

        default:
          break;
      }
    }
  } while (!quit && !reader.end());

  return EXIT_SUCCESS;
}
#else
int main() {
  std::cerr << "Current OpenCV version is " << VISP_HAVE_OPENCV_VERSION << std::endl;
  std::cerr << "OpenCV 2.3.0 or higher is requested to run the calibration." << std::endl;
  return EXIT_SUCCESS;
}
#endif

