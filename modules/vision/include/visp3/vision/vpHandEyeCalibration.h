/****************************************************************************
 *
 * ViSP, open source Visual Servoing Platform software.
 * Copyright (C) 2005 - 2019 by Inria. All rights reserved.
 *
 * This software is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * See the file LICENSE.txt at the root directory of this source
 * distribution for additional information about the GNU GPL.
 *
 * For using ViSP with software that can not be combined with the GNU
 * GPL, please contact Inria about acquiring a ViSP Professional
 * Edition License.
 *
 * See http://visp.inria.fr for more information.
 *
 * This software was developed at:
 * Inria Rennes - Bretagne Atlantique
 * Campus Universitaire de Beaulieu
 * 35042 Rennes Cedex
 * France
 *
 * If you have questions regarding the use of this file, please contact
 * Inria at visp@inria.fr
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Description:
 * Hand-eye calibration.
 *
 * Authors:
 * Francois Chaumette
 * Fabien Spindler
 *
 *****************************************************************************/

/*!
  \file vpHandEyeCalibration.h
  \brief Tools for hand-eye calibration.

  \sa The example in calibrate-hand-eye.cpp
*/
#ifndef _vpHandEyeCalibration_h_
#define _vpHandEyeCalibration_h_

#include <vector>
#include <visp3/core/vpExponentialMap.h>
#include <visp3/core/vpHomogeneousMatrix.h>
#include <visp3/core/vpMath.h>
#include <visp3/core/vpMatrix.h>

#if (VISP_HAVE_OPENCV_VERSION >= 0x020101) // Require opencv >= 2.1.1
#define HE_WO_VISP
#endif //   VISP_HAVE_OPENCV_VERSION

#ifdef HE_WO_VISP
#include <opencv2/core.hpp>
#include <opencv2/calib3d.hpp>
//#include <opencv2/sfm.hpp>
using namespace cv;
#endif  //  HE_WO_VISP

/*!
  \class vpHandEyeCalibration

  \ingroup group_vision_calib

  \brief Tool for hand-eye calibration.

*/
class VISP_EXPORT vpHandEyeCalibration
{
public:

  static int calibrate(const std::vector<vpHomogeneousMatrix> &cMo, const std::vector<vpHomogeneousMatrix> &rMe,
                       vpHomogeneousMatrix &eMc);
#ifdef HE_WO_VISP
    static int calibrate_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, Mat &eMc, int n_sp);
    static double rad2deg(double radian);
    static std::pair<cv::Mat, cv::Mat> split_homogeneous_transform_matrix_into_rotation_and_translation(const Mat& mat_homo, int n_sp);
#endif  //  HE_WO_VISP

private:
  static void calibrationVerifrMo(const std::vector<vpHomogeneousMatrix> &cMo, const std::vector<vpHomogeneousMatrix> &rMe,
                                  const vpHomogeneousMatrix &eMc);
  static int calibrationRotationTsai(const std::vector<vpHomogeneousMatrix> &cMo, const std::vector<vpHomogeneousMatrix> &rMe, vpRotationMatrix &eRc);
  static int calibrationRotationTsaiOld(const std::vector<vpHomogeneousMatrix> &cMo, const std::vector<vpHomogeneousMatrix> &rMe,
                                        vpRotationMatrix &eRc);
  static int calibrationRotationProcrustes(const std::vector<vpHomogeneousMatrix> &cMo, const std::vector<vpHomogeneousMatrix> &rMe, vpRotationMatrix &eRc);
  static int calibrationTranslation(const std::vector<vpHomogeneousMatrix> &cMo, const std::vector<vpHomogeneousMatrix> &rMe, vpRotationMatrix &eRc, vpTranslationVector &eTc);
  static int calibrationTranslationOld(const std::vector<vpHomogeneousMatrix> &cMo, const std::vector<vpHomogeneousMatrix> &rMe,
                                       vpRotationMatrix &eRc, vpTranslationVector &eTc);
  static double calibrationErrVVS(const std::vector<vpHomogeneousMatrix> &cMo, const std::vector<vpHomogeneousMatrix> &rMe,
                                  const vpHomogeneousMatrix &eMc, vpColVector &errVVS);
  static int calibrationVVS(const std::vector<vpHomogeneousMatrix> &cMo, const std::vector<vpHomogeneousMatrix> &rMe, vpHomogeneousMatrix &eMc);

#ifdef HE_WO_VISP
    static void calibrationVerifrMo_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, const Mat &eMc, int n_sp);
    static double calibrationErrVVS_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, const Mat &eMc, Mat& errVVS, int n_sp);
    static int calibrationRotationTsaiOld_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, Mat& eRc, int n_sp);
    static int calibrationTranslationOld_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, Mat &eRc, Mat &eTc, int n_sp);
    static int calibrationTranslation_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, Mat &eRc, Mat &eTc, int n_sp);
    static int calibrationRotationTsai_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, Mat &eRc, int n_sp);
    static int calibrationRotationProcrustes_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, Mat &eRc, int n_sp);
    static int calibrationVVS_wo_visp(const std::vector<Mat> &cMo, const std::vector<Mat> &rMe, Mat &eMc, int n_sp);
#endif  //  HE_WO_wo_visp_VISP  

};

#endif
