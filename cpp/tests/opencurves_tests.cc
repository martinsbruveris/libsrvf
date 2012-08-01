#include <fstream>
#include <cmath>
#include <boost/test/unit_test.hpp>

#include <srvf/matrix.h>
#include <srvf/pointset.h>
#include <srvf/srvf.h>
#include <srvf/qmap.h>
#include <srvf/opencurves.h>
#include <srvf/rotate.h>
#include <srvf/fileio.h>
#include <srvf/render.h>
#include <srvf/plotwin.h>

#include <FL/Fl.h>


BOOST_AUTO_TEST_SUITE(opencurves_tests)

BOOST_AUTO_TEST_CASE(shooting_vector_test1)
{
  double samps1_data[] = {1.0};
  double samps2_data[] = {1.0, -1.0};
  double exp_data[] = {M_PI_2, -M_PI_2};
  double exp_params_data[] = {0.0, 0.5, 1.0};
  size_t exp_ncp = sizeof(exp_params_data) / sizeof(double);

  srvf::Pointset samps1(1, 1, samps1_data);
  srvf::Pointset samps2(1, 2, samps2_data);
  srvf::Srvf Q1(samps1);
  srvf::Srvf Q2(samps2);

  srvf::Srvf Sv = srvf::opencurves::shooting_vector(Q1, Q2);
  BOOST_CHECK_EQUAL(Sv.dim(), 1);
  BOOST_CHECK_EQUAL(Sv.ncp(), exp_ncp);
  for (size_t i=0; i<exp_ncp; ++i)
  {
    BOOST_CHECK_EQUAL(Sv.params()[i], exp_params_data[i]);
  }
  for (size_t i=0; i+1<exp_ncp; ++i)
  {
    BOOST_CHECK_EQUAL(Sv.samps()[i][0], exp_data[i]);
  }
}

BOOST_AUTO_TEST_CASE(karcher_mean_test1)
{
  std::ifstream ifs("data/rna.csv");
  BOOST_REQUIRE_EQUAL(ifs.good(), true);
  std::vector<srvf::Matrix> samps_data = srvf::io::load_csv(ifs);
  std::vector<srvf::Srvf> Qs;
  std::vector<srvf::Plf> Fs;

  for (size_t i=0; i<samps_data.size(); ++i)
  {
    srvf::Pointset samps(samps_data[i], srvf::Pointset::POINT_PER_COLUMN);
    Fs.push_back(srvf::Plf(samps));
    Fs[i].translate_to_origin();
    Fs[i].scale_to_unit_arc_length();
    Qs.push_back(srvf::plf_to_srvf(Fs[i]));
  }

  srvf::Srvf Mu = srvf::opencurves::karcher_mean(Qs, 5e-3, 5);

  srvf::plot::Plot3D plot;
  plot.insert(srvf::srvf_to_plf(Mu), srvf::plot::Color(1.0,0.0,0.0));
  for (size_t i=0; i<Qs.size(); ++i)
  {
    srvf::Matrix Ri = srvf::optimal_rotation(Mu, Qs[i]);
    Fs[i].rotate(Ri);
    plot.insert(Fs[i], srvf::plot::Color(0.0,0.0,1.0));
  }
  srvf::plot::FltkGlPlotWindow win(400, 400, 
      "Karcher Mean of Serveral RNA Molecules");
  win.add_plot(&plot);
  win.show();
  Fl::run();
}

BOOST_AUTO_TEST_SUITE_END()
