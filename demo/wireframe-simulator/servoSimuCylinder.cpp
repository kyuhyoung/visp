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
 * Demonstration of the wireframe simulator with a simple visual servoing
 *
 * Authors:
 * Nicolas Melchior
 *
 *****************************************************************************/

/*!
  \example servoSimuCylinder.cpp

  Demonstration of the wireframe simulator with a simple visual servoing.
*/

#include <stdlib.h>

#include <visp3/core/vpCameraParameters.h>
#include <visp3/core/vpCylinder.h>
#include <visp3/core/vpHomogeneousMatrix.h>
#include <visp3/core/vpImage.h>
#include <visp3/core/vpIoTools.h>
#include <visp3/core/vpMath.h>
#include <visp3/core/vpTime.h>
#include <visp3/core/vpVelocityTwistMatrix.h>
#include <visp3/gui/vpDisplayD3D.h>
#include <visp3/gui/vpDisplayGDI.h>
#include <visp3/gui/vpDisplayGTK.h>
#include <visp3/gui/vpDisplayOpenCV.h>
#include <visp3/gui/vpDisplayX.h>
#include <visp3/gui/vpPlot.h>
#include <visp3/io/vpImageIo.h>
#include <visp3/io/vpParseArgv.h>
#include <visp3/robot/vpSimulatorCamera.h>
#include <visp3/robot/vpWireFrameSimulator.h>
#include <visp3/visual_features/vpFeatureBuilder.h>
#include <visp3/vs/vpServo.h>

#define GETOPTARGS "dhp"

#ifdef VISP_HAVE_DISPLAY

/*!

  Print the program options.

  \param name : Program name.
  \param badparam : Bad parameter name.

*/
void usage(const char *name, const char *badparam)
{
  fprintf(stdout, "\n\
Demonstration of the wireframe simulator with a simple visual servoing.\n\
          \n\
The visual servoing consists in bringing the camera at a desired position\n\
from the object.\n\
          \n\
The visual features used to compute the pose of the camera and \n\
thus the control law are two lines. These features are computed thanks \n\
to the equation of a cylinder.\n\
          \n\
This demonstration explains also how to move the object around a world \n\
reference frame. Here, the movment is a rotation around the x and y axis \n\
at a given distance from the world frame. In fact the object trajectory \n\
is on a sphere whose center is the origin of the world frame.\n\
          \n\
SYNOPSIS\n\
  %s [-d] [-h]\n", name);

  fprintf(stdout, "\n\
OPTIONS:                                          \n\
  -d \n\
     Turn off the display.\n\
            \n\
  -p \n\
     Turn off the plotter.\n\
                            \n\
  -h\n\
     Print the help.\n");

  if (badparam)
    fprintf(stdout, "\nERROR: Bad parameter [%s]\n", badparam);
}

/*!

  Set the program options.

  \param argc : Command line number of parameters.
  \param argv : Array of command line parameters.
  \param display : Display activation.
  \param plot : Plotter activation.

  \return false if the program has to be stopped, true otherwise.

*/
bool getOptions(int argc, const char **argv, bool &display, bool &plot)
{
  const char *optarg_;
  int c;
  while ((c = vpParseArgv::parse(argc, argv, GETOPTARGS, &optarg_)) > 1) {

    switch (c) {
    case 'd':
      display = false;
      break;
    case 'p':
      plot = false;
      break;
    case 'h':
      usage(argv[0], NULL);
      return false;

    default:
      usage(argv[0], optarg_);
      return false;
    }
  }

  if ((c == 1) || (c == -1)) {
    // standalone param or error
    usage(argv[0], NULL);
    std::cerr << "ERROR: " << std::endl;
    std::cerr << "  Bad argument " << optarg_ << std::endl << std::endl;
    return false;
  }

  return true;
}

int main(int argc, const char **argv)
{
  try {
    bool opt_display = true;
    bool opt_plot = true;

    // Read the command line options
    if (getOptions(argc, argv, opt_display, opt_plot) == false) {
      exit(-1);
    }

    vpImage<vpRGBa> Iint(480, 640, 255);
    vpImage<vpRGBa> Iext(480, 640, 255);

#if defined VISP_HAVE_X11
    vpDisplayX display[2];
#elif defined VISP_HAVE_OPENCV
    vpDisplayOpenCV display[2];
#elif defined VISP_HAVE_GDI
    vpDisplayGDI display[2];
#elif defined VISP_HAVE_D3D9
    vpDisplayD3D display[2];
#elif defined VISP_HAVE_GTK
    vpDisplayGTK display[2];
#endif

    if (opt_display) {
      // Display size is automatically defined by the image (I) size
      display[0].init(Iint, 100, 100, "The internal view");
      display[1].init(Iext, 100, 100, "The first external view");
      vpDisplay::setWindowPosition(Iint, 0, 0);
      vpDisplay::setWindowPosition(Iext, 750, 0);
      vpDisplay::display(Iint);
      vpDisplay::flush(Iint);
      vpDisplay::display(Iext);
      vpDisplay::flush(Iext);
    }

    vpPlot *plotter = NULL;

    vpServo task;
    vpSimulatorCamera robot;
    float sampling_time = 0.020f; // Sampling period in second
    robot.setSamplingTime(sampling_time);

    // Set initial position of the object in the camera frame
    vpHomogeneousMatrix cMo(0, 0.1, 0.3, vpMath::rad(35), vpMath::rad(25), vpMath::rad(75));
    // Set desired position of the object in the camera frame
    vpHomogeneousMatrix cdMo(0.0, 0.0, 0.5, vpMath::rad(90), vpMath::rad(0), vpMath::rad(0));
    // Set initial position of the object in the world frame
    vpHomogeneousMatrix wMo(0.0, 0.0, 0, 0, 0, 0);
    // Position of the camera in the world frame
    vpHomogeneousMatrix wMc;
    wMc = wMo * cMo.inverse();

    // Create a cylinder
    vpCylinder cylinder(0, 0, 1, 0, 0, 0, 0.1);

    // Projection of the cylinder
    cylinder.track(cMo);

    // Set the current visual feature
    vpFeatureLine l[2];
    vpFeatureBuilder::create(l[0], cylinder, vpCylinder::line1);
    vpFeatureBuilder::create(l[1], cylinder, vpCylinder::line2);

    // Projection of the cylinder
    cylinder.track(cdMo);

    vpFeatureLine ld[2];
    vpFeatureBuilder::create(ld[0], cylinder, vpCylinder::line1);
    vpFeatureBuilder::create(ld[1], cylinder, vpCylinder::line2);

    task.setServo(vpServo::EYEINHAND_L_cVe_eJe);
    task.setInteractionMatrixType(vpServo::DESIRED);

    vpHomogeneousMatrix cMe;
    vpVelocityTwistMatrix cVe(cMe);
    task.set_cVe(cVe);

    vpMatrix eJe;
    robot.get_eJe(eJe);
    task.set_eJe(eJe);

    for (int i = 0; i < 2; i++)
      task.addFeature(l[i], ld[i]);

    if (opt_plot) {
      plotter = new vpPlot(2, 480, 640, 750, 550, "Real time curves plotter");
      plotter->setTitle(0, "Visual features error");
      plotter->setTitle(1, "Camera velocities");
      plotter->initGraph(0, task.getDimension());
      plotter->initGraph(1, 6);
      plotter->setLegend(0, 0, "error_feat_l1_rho");
      plotter->setLegend(0, 1, "error_feat_l1_theta");
      plotter->setLegend(0, 2, "error_feat_l2_rho");
      plotter->setLegend(0, 3, "error_feat_l2_theta");
      plotter->setLegend(1, 0, "vc_x");
      plotter->setLegend(1, 1, "vc_y");
      plotter->setLegend(1, 2, "vc_z");
      plotter->setLegend(1, 3, "wc_x");
      plotter->setLegend(1, 4, "wc_y");
      plotter->setLegend(1, 5, "wc_z");
    }

    task.setLambda(1);

    vpWireFrameSimulator sim;

    // Set the scene
    sim.initScene(vpWireFrameSimulator::CYLINDER, vpWireFrameSimulator::D_STANDARD);

    // Initialize simulator frames
    sim.set_fMo(wMo);                   // Position of the object in the world reference frame
    sim.setCameraPositionRelObj(cMo);   // Initial position of the object in the camera frame
    sim.setDesiredCameraPosition(cdMo); // Desired position of the object in the camera frame

    // Set the External camera position
    vpHomogeneousMatrix camMf(vpHomogeneousMatrix(0.0, 0, 3.5, vpMath::rad(0), vpMath::rad(30), 0));
    sim.setExternalCameraPosition(camMf);

    // Set the parameters of the cameras (internal and external)
    vpCameraParameters camera(1000, 1000, 320, 240);
    sim.setInternalCameraParameters(camera);
    sim.setExternalCameraParameters(camera);

    int max_iter = 10;

    if (opt_display) {
      max_iter = 2500;

      // Get the internal and external views
      sim.getInternalImage(Iint);
      sim.getExternalImage(Iext);

      // Display the object frame (current and desired position)
      vpDisplay::displayFrame(Iint, cMo, camera, 0.2, vpColor::none);
      vpDisplay::displayFrame(Iint, cdMo, camera, 0.2, vpColor::none);

      // Display the object frame the world reference frame and the camera
      // frame
      vpDisplay::displayFrame(Iext, camMf * sim.get_fMo() * cMo.inverse(), camera, 0.2, vpColor::none);
      vpDisplay::displayFrame(Iext, camMf * sim.get_fMo(), camera, 0.2, vpColor::none);
      vpDisplay::displayFrame(Iext, camMf, camera, 0.2, vpColor::none);

      vpDisplay::displayText(Iint, 20, 20, "Click to start visual servo", vpColor::red);

      vpDisplay::flush(Iint);
      vpDisplay::flush(Iext);

      std::cout << "Click on a display" << std::endl;
      while (!vpDisplay::getClick(Iint, false) && !vpDisplay::getClick(Iext, false)) {
      };
    }

    robot.setPosition(wMc);

    // Print the task
    task.print();

    int iter = 0;
    bool stop = false;
    vpColVector v;

    // Set the secondary task parameters
    vpColVector e1(6, 0);
    vpColVector e2(6, 0);
    vpColVector proj_e1;
    vpColVector proj_e2;
    double rapport = 0;
    double vitesse = 0.3;
    int tempo = 600;

    double t_prev, t = vpTime::measureTimeMs();

    while (iter++ < max_iter && !stop) {
      t_prev = t;
      t = vpTime::measureTimeMs();

      if (opt_display) {
        vpDisplay::display(Iint);
        vpDisplay::display(Iext);
      }

      robot.get_eJe(eJe);
      task.set_eJe(eJe);

      wMc = robot.getPosition();
      cMo = wMc.inverse() * wMo;

      cylinder.track(cMo);
      vpFeatureBuilder::create(l[0], cylinder, vpCylinder::line1);
      vpFeatureBuilder::create(l[1], cylinder, vpCylinder::line2);

      v = task.computeControlLaw();

      // Compute the velocity with the secondary task
      if (iter % tempo < 200 && iter % tempo >= 0) {
        e2 = 0;
        e1[0] = -fabs(vitesse);
        proj_e1 = task.secondaryTask(e1, true);
        rapport = -vitesse / proj_e1[0];
        proj_e1 *= rapport;
        v += proj_e1;
      }

      else if (iter % tempo < 300 && iter % tempo >= 200) {
        e1 = 0;
        e2[1] = -fabs(vitesse);
        proj_e2 = task.secondaryTask(e2, true);
        rapport = -vitesse / proj_e2[1];
        proj_e2 *= rapport;
        v += proj_e2;
      }

      else if (iter % tempo < 500 && iter % tempo >= 300) {
        e2 = 0;
        e1[0] = -fabs(vitesse);
        proj_e1 = task.secondaryTask(e1, true);
        rapport = vitesse / proj_e1[0];
        proj_e1 *= rapport;
        v += proj_e1;
      }

      else if (iter % tempo < 600 && iter % tempo >= 500) {
        e1 = 0;
        e2[1] = -fabs(vitesse);
        proj_e2 = task.secondaryTask(e2, true);
        rapport = vitesse / proj_e2[1];
        proj_e2 *= rapport;
        v += proj_e2;
      }

      robot.setVelocity(vpRobot::CAMERA_FRAME, v);

      // Update the simulator frames
      sim.set_fMo(wMo); // This line is not really requested since the object
                        // doesn't move
      sim.setCameraPositionRelObj(cMo);

      if (opt_plot) {
        plotter->plot(0, iter, task.getError());
        plotter->plot(1, iter, v);
      }

      if (opt_display) {
        // Get the internal and external views
        sim.getInternalImage(Iint);
        sim.getExternalImage(Iext);

        // Display the object frame (current and desired position)
        vpDisplay::displayFrame(Iint, cMo, camera, 0.2, vpColor::none);
        vpDisplay::displayFrame(Iint, cdMo, camera, 0.2, vpColor::none);

        // Display the object frame the world reference frame and the camera
        // frame
        vpDisplay::displayFrame(Iext, sim.getExternalCameraPosition() * sim.get_fMo() * cMo.inverse(), camera, 0.2,
                                vpColor::none);
        vpDisplay::displayFrame(Iext, sim.getExternalCameraPosition() * sim.get_fMo(), camera, 0.2, vpColor::none);
        vpDisplay::displayFrame(Iext, sim.getExternalCameraPosition(), camera, 0.2, vpColor::none);

        vpDisplay::displayText(Iint, 20, 20, "Click to stop visual servo", vpColor::red);

        std::stringstream ss;
        ss << "Loop time: " << t - t_prev << " ms";
        vpDisplay::displayText(Iint, 40, 20, ss.str(), vpColor::red);

        if (vpDisplay::getClick(Iint, false)) {
          stop = true;
        }

        vpDisplay::flush(Iext);
        vpDisplay::flush(Iint);

        vpTime::wait(t, sampling_time * 1000); // Wait ms
      }

      std::cout << "|| s - s* || = " << (task.getError()).sumSquare() << std::endl;
    }

    if (opt_plot && plotter != NULL) {
      vpDisplay::display(Iint);
      sim.getInternalImage(Iint);
      vpDisplay::displayFrame(Iint, cMo, camera, 0.2, vpColor::none);
      vpDisplay::displayFrame(Iint, cdMo, camera, 0.2, vpColor::none);
      vpDisplay::displayText(Iint, 20, 20, "Click to quit", vpColor::red);
      if (vpDisplay::getClick(Iint)) {
        stop = true;
      }
      vpDisplay::flush(Iint);

      delete plotter;
    }

    task.print();
    task.kill();

    return EXIT_SUCCESS;
  } catch (const vpException &e) {
    std::cout << "Catch an exception: " << e << std::endl;
    return EXIT_FAILURE;
  }
}
#else
int main()
{
  std::cout << "You do not have X11, or GDI (Graphical Device Interface), or GTK functionalities to display images..." << std::endl;
  std::cout << "Tip if you are on a unix-like system:" << std::endl;
  std::cout << "- Install X11, configure again ViSP using cmake and build again this example" << std::endl;
  std::cout << "Tip if you are on a windows-like system:" << std::endl;
  std::cout << "- Install GDI, configure again ViSP using cmake and build again this example" << std::endl;
  return EXIT_SUCCESS;
}

#endif
