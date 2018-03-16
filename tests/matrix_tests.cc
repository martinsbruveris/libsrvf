#include <cmath>
#include <boost/test/unit_test.hpp>

#include <srvf/matrix.h>

BOOST_AUTO_TEST_SUITE(matrix_tests)

BOOST_AUTO_TEST_CASE(ctor_test1)
{
  srvf::Matrix A;
  BOOST_CHECK_EQUAL(A.rows(),0);
  BOOST_CHECK_EQUAL(A.cols(),0);
  BOOST_CHECK_EQUAL(A.size(),0);
  BOOST_CHECK_EQUAL(A.data().size(),0);
}

BOOST_AUTO_TEST_CASE(ctor_test2)
{
  srvf::Matrix A(5000,5000);
  BOOST_CHECK_EQUAL(A.rows(),5000);
  BOOST_CHECK_EQUAL(A.cols(),5000);
  BOOST_CHECK_EQUAL(A.size(),25000000);
  BOOST_CHECK_EQUAL(A.data().size(),25000000);
}

BOOST_AUTO_TEST_CASE(ctor_test3)
{
  srvf::Matrix A(3,500,0.5);
  BOOST_CHECK_EQUAL(A.rows(),3);
  BOOST_CHECK_EQUAL(A.cols(),500);
  BOOST_CHECK_EQUAL(A.size(),1500);
  BOOST_CHECK_EQUAL(A.data().size(),1500);
  // check entries using 2-argument accessor
  for (size_t i=0; i<A.rows(); ++i){
    for (size_t j=0; j<A.cols(); ++j){
      BOOST_CHECK_EQUAL(A(i,j),0.5);
    }
  }
}

BOOST_AUTO_TEST_CASE(ctor_test4)
{
  double Adata[]={
    0.0, 1.0, 2.0, 3.0,
    4.0, 5.0, 6.0, 7.0
  };
  size_t rows=2;
  size_t cols=4;
  size_t sz=rows*cols;
  srvf::Matrix A(rows, cols, Adata);
  BOOST_CHECK_EQUAL(A.rows(),rows);
  BOOST_CHECK_EQUAL(A.cols(),cols);
  for (size_t i=0; i<sz; ++i)
  {
    BOOST_CHECK_EQUAL(Adata[i],A(i));
  }
}

BOOST_AUTO_TEST_CASE(ctor_test5)
{
  srvf::Matrix A(20,20,0.125);
  srvf::Matrix B(A);
  BOOST_CHECK_EQUAL(A.rows(),B.rows());
  BOOST_CHECK_EQUAL(A.cols(),B.cols());
  for (size_t i=0; i<A.rows(); ++i)
  {
    for (size_t j=0; j<A.cols(); ++j)
    {
      BOOST_CHECK_EQUAL(A(i,j), B(i,j));
    }
  }
}

BOOST_AUTO_TEST_CASE(assignment_test1)
{
  srvf::Matrix A(20,20,0.5);
  srvf::Matrix B;
  B=A;
  BOOST_CHECK_EQUAL(A.rows(),B.rows());
  BOOST_CHECK_EQUAL(A.cols(),B.cols());
  BOOST_CHECK_EQUAL(A.size(),B.size());
  for (size_t i=0; i<A.size(); ++i)
  {
    BOOST_CHECK_EQUAL(A(i),B(i));
  }
}

BOOST_AUTO_TEST_CASE(clear_test1)
{
  srvf::Matrix A(33,147,19.35);
  A.clear();
  BOOST_CHECK_EQUAL(A.rows(),0);
  BOOST_CHECK_EQUAL(A.cols(),0);
  BOOST_CHECK_EQUAL(A.size(),0);
  BOOST_CHECK_EQUAL(A.data().size(),0);
}

BOOST_AUTO_TEST_CASE(resize_test1)
{
  srvf::Matrix A;
  A.resize(5,5);
  BOOST_CHECK_EQUAL(A.rows(),5);
  BOOST_CHECK_EQUAL(A.cols(),5);
  BOOST_CHECK_EQUAL(A.size(),25);
}

BOOST_AUTO_TEST_CASE(resize_test2)
{
  srvf::Matrix A(5,500,0.25);
  A.resize(10,250);
  BOOST_CHECK_EQUAL(A.rows(),10);
  BOOST_CHECK_EQUAL(A.cols(),250);
  BOOST_CHECK_EQUAL(A.size(),2500);
  for (size_t i=0; i<5; ++i)
  {
    for (size_t j=0; j<250; ++j)
    {
      BOOST_CHECK_EQUAL(A(i,j),0.25);
    }
  }
}

BOOST_AUTO_TEST_CASE(resize_test3)
{
  srvf::Matrix A(500,5,0.25);
  A.resize(250,10);
  BOOST_CHECK_EQUAL(A.rows(),250);
  BOOST_CHECK_EQUAL(A.cols(),10);
  BOOST_CHECK_EQUAL(A.size(),2500);
  for (size_t i=0; i<250; ++i)
  {
    for (size_t j=0; j<5; ++j)
    {
      BOOST_CHECK_EQUAL(A(i,j),0.25);
    }
  }
}

BOOST_AUTO_TEST_CASE(inplace_plus_test1)
{
  srvf::Matrix A(1,20);
  srvf::Matrix B(1,20);
  for (size_t i=0; i<20; ++i)
  {
    A(i)=(double)i;
    B(i)=-(double)i;
  }
  A += B;
  for (size_t i=0; i<20; ++i)
  {
    BOOST_CHECK_EQUAL(A(i),0.0);
  }
}

BOOST_AUTO_TEST_CASE(inplace_minus_test1)
{
  srvf::Matrix A(50,200,16.0);
  srvf::Matrix B(50,200,8.0);
  A -= B;
  for (size_t i=0; i<A.size(); ++i)
  {
    BOOST_CHECK_EQUAL(A(i),8.0);
  }
}

// Test case for A*=A
BOOST_AUTO_TEST_CASE(inplace_times_test1)
{
  srvf::Matrix A(1,20);
  for (size_t i=0; i<20; ++i)
  {
    A(i)=(double)(i);
  }
  A *= A;
  for (size_t i=0; i<20; ++i)
  {
    BOOST_CHECK_EQUAL(A(i),(double)(i*i));
  }
}

// Test case for A*=B
BOOST_AUTO_TEST_CASE(inplace_times_test2)
{
  srvf::Matrix A(100,20,2.0);
  srvf::Matrix B(100,20,0.25);
  A *= B;
  for (size_t i=0; i<A.size(); ++i)
  {
    BOOST_CHECK_EQUAL(A(i),0.5);
  }
}

// Test case for A/=A
BOOST_AUTO_TEST_CASE(inplace_div_test1)
{
  srvf::Matrix A(1,20,16.0);
  A /= A;
  for (size_t i=0; i<20; ++i)
  {
    BOOST_CHECK_EQUAL(A(i),1.0);
  }
}

// Test case for A/=B with division by zero
BOOST_AUTO_TEST_CASE(inplace_div_test2)
{
  srvf::Matrix A(10,10,4.0);
  srvf::Matrix B(10,10,2.0);
  B(5,5)=0.0;
  A /= B;
  BOOST_CHECK(std::isinf(A(5,5)));
  for (size_t i=0; i<10; ++i)
  {
    for (size_t j=0; j<10; ++j)
    {
      if (i!=5 || j!=5)
      {
        BOOST_CHECK_EQUAL(A(i,j),2.0);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(inplace_scalar_plus_test1)
{
  srvf::Matrix A(3,4,0.0);
  A+=1.0;
  for (size_t i=0; i<A.size(); ++i)
  {
    BOOST_CHECK_EQUAL(A(i),1.0);
  }
}

BOOST_AUTO_TEST_CASE(inplace_scalar_minus_test1)
{
  srvf::Matrix A(3,4,0.0);
  A-=1.0;
  for (size_t i=0; i<A.size(); ++i)
  {
    BOOST_CHECK_EQUAL(A(i),-1.0);
  }
}

BOOST_AUTO_TEST_CASE(inplace_scalar_times_test1)
{
  srvf::Matrix A(3,4,1.0);
  A*=0.25;
  for (size_t i=0; i<A.size(); ++i)
  {
    BOOST_CHECK_EQUAL(A(i),0.25);
  }
}

BOOST_AUTO_TEST_CASE(inplace_scalar_div_test1)
{
  srvf::Matrix A(3,4,1.0);
  A/=0.25;
  for (size_t i=0; i<A.size(); ++i)
  {
    BOOST_CHECK_EQUAL(A(i),4.0);
  }
}

BOOST_AUTO_TEST_CASE(scalar_plus_test1)
{
  srvf::Matrix A(20,32,0.25);
  srvf::Matrix B=A+0.75;
  for (size_t i=0; i<B.size(); ++i)
  {
    BOOST_CHECK_EQUAL(B(i),1.0);
  }
}

BOOST_AUTO_TEST_CASE(scalar_minus_test1)
{
  srvf::Matrix A(20,32,0.25);
  srvf::Matrix B=A-0.75;
  for (size_t i=0; i<B.size(); ++i)
  {
    BOOST_CHECK_EQUAL(B(i),-0.5);
  }
}

BOOST_AUTO_TEST_CASE(scalar_times_test1)
{
  srvf::Matrix A(20,32,1.0);
  srvf::Matrix B=A*0.5;
  for (size_t i=0; i<B.size(); ++i)
  {
    BOOST_CHECK_EQUAL(B(i),0.5);
  }
}

BOOST_AUTO_TEST_CASE(scalar_div_test1)
{
  srvf::Matrix A(20,32,1.0);
  srvf::Matrix B=A/0.5;
  for (size_t i=0; i<B.size(); ++i)
  {
    BOOST_CHECK_EQUAL(B(i),2.0);
  }
}

BOOST_AUTO_TEST_CASE(matmul_test1)
{
  srvf::Matrix A(2,2,1.0);
  srvf::Matrix A2 = srvf::product(A,A);
  for (size_t i=0; i<4; ++i)
  {
    BOOST_CHECK_EQUAL(A2(i),2.0);
  }
}

BOOST_AUTO_TEST_CASE(product_matvec_test1)
{
  srvf::Matrix A = srvf::Matrix::identity(2);
  srvf::Point P1(2,1.0);
  srvf::Point P2 = product(A, P1);
  srvf::Point P3 = product(P1, A);

  BOOST_REQUIRE_EQUAL((P1 == P2), true);
  BOOST_REQUIRE_EQUAL((P1 == P3), true);
}

BOOST_AUTO_TEST_CASE(transpose_test1)
{
  srvf::Matrix A(2,3);
  for (size_t i=0; i<2; ++i)
  {
    for (size_t j=0; j<3; ++j)
    {
      A(i,j)=3*i+j;
    }
  }
  srvf::Matrix At=srvf::transpose(A);
  BOOST_CHECK_EQUAL(At.rows(),A.cols());
  BOOST_CHECK_EQUAL(At.cols(),A.rows());
  for (size_t i=0; i<At.rows(); ++i)
  {
    for (size_t j=0; j<At.cols(); ++j)
    {
      BOOST_CHECK_EQUAL(At(i,j),A(j,i));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
