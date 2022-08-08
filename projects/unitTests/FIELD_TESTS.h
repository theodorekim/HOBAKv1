/*
This file is part of HOBAK.

Permission is hereby granted to use this software solely for non-commercial applications
and purposes including academic or industrial research, evaluation and not-for-profit media
production.
 
THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PIXAR OR ITS AFFILIATES, YALE, OR
THE AUTHORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
///////////////////////////////////////////////////////////////////////
// generate a K2 kernel -- this is here just so unit testing is
// self-contained
///////////////////////////////////////////////////////////////////////
FIELD_3D generateK2Kernel3D()
{
  cout << " Building kernel ..."; flush(cout);
  int filterWidth = 15;
  int filterHalf = filterWidth / 2;
  FIELD_3D k2Kernel(filterWidth, filterWidth, filterWidth);

  double dk = 0.01;
  double sigma = 1.0;
  double norm = 0;

  double kMin = 2;
  double kMax = 10;

  for (double k = kMin; k < kMax; k += dk)
    norm += k * k * exp(-sigma * k * k);

  for (int h = -filterHalf; h <= filterHalf; h++)
    for (int i = -filterHalf; i <= filterHalf; i++)
      for (int j = -filterHalf; j <= filterHalf; j++)
      {
        double r = sqrt((float)(i * i + j * j + h * h));
        double kern = 0;

        for (double k = kMin; k < kMax; k += dk)
        {
          double kr = k * r;
          kern += (sin(kr) / kr) * exp(-sigma * k * k) * k * k * k / M_PI;
        }

        double interp = ((r /filterHalf) - 0.9) / 0.1;
        REAL squared = interp * interp;
        interp = 2 * squared * interp - 3 * squared + 1;
        if (interp < 0) interp = 1;
        if (interp > 1) interp = 0;

        kern *= interp;
        k2Kernel(i + filterHalf, j + filterHalf, h + filterHalf) = kern / norm;
      }

  k2Kernel(filterHalf, filterHalf, filterHalf) = 0;

  // get the center integral
  double center = 0;
  for (double k = kMin; k < kMax; k += dk)
    center += k * k * exp(-sigma * k * k) * j0(0);
  center *= 1.0 / norm;

  // see what the current integral works out to
  FIELD_2D projection = k2Kernel.zProjection();
  double centerIntegral = projection(filterHalf, filterHalf);

  // make the center integral work out right
  double newCenter = center - centerIntegral;
  k2Kernel(filterHalf, filterHalf, filterHalf) = newCenter;
  cout << " done." << endl;

  return k2Kernel;
}

//////////////////////////////////////////////////////////////////////////////
// This tests whether the fast convolution (not FFT-based) in FIELD_3D is
// correct
//////////////////////////////////////////////////////////////////////////////
TEST_CASE("Fast FIELD_3D convolution", "[FIELD_3D convolution]" )
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING fast convolution in FIELD_3D " << endl;
  cout << "=============================================================== " << endl;

  // generate the convolution kernel
  FIELD_3D k2Kernel = generateK2Kernel3D();

  int res = 64;
  FIELD_3D source(res, res, res);
  source.setToSolidCheckboard(5,5,5);

  cout << " Performing full convolution ..."; flush(cout);
  FIELD_3D full = source.convolve(k2Kernel);
  cout << " done." << endl;

  cout << " Performing fast 15 convolution ..."; flush(cout);
  FIELD_3D fast15 = source.convolveFast15(k2Kernel);
  cout << " done." << endl;

  FIELD_3D diff = full - fast15;
  REAL diffMax = diff.absMax();

  bool success = false;

  cout << " Maximum diff found: " << diffMax << endl;
  if (diffMax < 1e-7)
  {
    cout << " Convolution unit test PASSED." << endl;
    success = true;;
  }
  else
  {
    cout << " Convolution unit test FAILED!!!!!!!!!!!!!!!!!" << endl;
  }
  REQUIRE(success);
  //FIELDVIEW3D(diff);
}
