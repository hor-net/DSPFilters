/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_RBJ_H
#define DSPFILTERS_RBJ_H

#include "DspFilters/Common.h"
#include "DspFilters/Biquad.h"
#include "DspFilters/Design.h"
#include "DspFilters/Filter.h"
#define _USE_MATH_DEFINES
#include <math.h>

namespace Dsp {

/*
 * Filter realizations based on Robert Bristol-Johnson formulae:
 *
 * http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 *
 */

namespace RBJ {


typedef int    Int32;  // Must be 32 bits

// Must be 64 bits
#ifdef _MSC_VER
typedef __int64 Int64;
#else
typedef int64_t Int64;
#endif


// Lower precision than ::fastsincos() but still decent
inline void fastersincos( double x, double *sn, double *cs )
{
  //always wrap input angle to -PI..PI
  if    (x < -M_PI) x += 2*M_PI;
  else if (x >  M_PI) x -= 2*M_PI;
  //compute sine
  if (x < 0)  *sn = 1.27323954 * x + 0.405284735 * x * x;
  else    *sn = 1.27323954 * x - 0.405284735 * x * x;
  //compute cosine: sin(x + PI/2) = cos(x)
  x += M_PI_2;
  if (x > M_PI ) x -= 2*M_PI;
  if (x < 0)  *cs = 1.27323954 * x + 0.405284735 * x * x;
  else    *cs = 1.27323954 * x - 0.405284735 * x * x;
}

// Slower than ::fastersincos() but still faster than
// sin(), and with the best accuracy of these routines.
inline void fastsincos( double x, double *sn, double *cs )
{
  double s, c;
  //always wrap input angle to -PI..PI
     if (x < -M_PI) x += 2*M_PI;
  else if (x >  M_PI) x -= 2*M_PI;
  //compute sine
  if (x < 0)
  {
    s = 1.27323954 * x + .405284735 * x * x;
    if (s < 0)  s = .225 * (s * -s - s) + s;
    else    s = .225 * (s *  s - s) + s;
  }
  else
  {
    s = 1.27323954 * x - 0.405284735 * x * x;
    if (s < 0)  s = .225 * (s * -s - s) + s;
    else    s = .225 * (s *  s - s) + s;
  }
  *sn=s;
  //compute cosine: sin(x + PI/2) = cos(x)
  x += M_PI_2;
  if (x > M_PI ) x -= 2*M_PI;
  if (x < 0)
  {
    c = 1.27323954 * x + 0.405284735 * x * x;
    if (c < 0)  c = .225 * (c * -c - c) + c;
    else    c = .225 * (c *  c - c) + c;
  }
  else
  {
    c = 1.27323954 * x - 0.405284735 * x * x;
    if (c < 0)  c = .225 * (c * -c - c) + c;
    else    c = .225 * (c *  c - c) + c;
  }
  *cs=c;
}

inline float fastsqrt1( float x )
{
  union { Int32 i; float x; } u;
  u.x = x;
  u.i = (Int32(1)<<29) + (u.i >> 1) - (Int32(1)<<22);
  return u.x;
}

inline double fastsqrt1( double x )
{
  union { Int64 i; double x; } u;
  u.x = x;
  u.i = (Int64(1)<<61) + (u.i >> 1) - (Int64(1)<<51);
  return u.x;
}

static double fastsinh (double x) noexcept
{
  auto x2 = x * x;
  auto numerator = -x * (11511339840 + x2 * (1640635920 + x2 * (52785432 + x2 * 479249)));
  auto denominator = -11511339840 + x2 * (277920720 + x2 * (-3177720 + x2 * 18361));
  return numerator / denominator;
}



//
// Raw filters
//

struct LowPass : BiquadBase
{
  void setup (double sampleRate,
              double cutoffFrequency,
              double q);
};

struct HighPass : BiquadBase
{
  void setup (double sampleRate,
              double cutoffFrequency,
              double q);
};

struct BandPass1 : BiquadBase
{
  // (constant skirt gain, peak gain = Q)
  void setup (double sampleRate,
              double centerFrequency,
              double bandWidth);
};

struct BandPass2 : BiquadBase
{
  // (constant 0 dB peak gain)
  void setup (double sampleRate,
              double centerFrequency,
              double bandWidth);
};

struct BandStop : BiquadBase
{
  void setup (double sampleRate,
              double centerFrequency,
              double bandWidth);
};

struct LowShelf : BiquadBase
{
  void setup (double sampleRate,
              double cutoffFrequency,
              double gainDb,
              double shelfSlope);
};

struct HighShelf : BiquadBase
{
  void setup (double sampleRate,
              double cutoffFrequency,
              double gainDb,
              double shelfSlope);
};
	
struct HighShelfOld : BiquadBase
{
  void setup (double sampleRate,
			  double cutoffFrequency,
			  double gainDb,
			  double shelfSlope);
};

struct BandShelf : BiquadBase
{
  void setup (double sampleRate,
              double centerFrequency,
              double gainDb,
              double bandWidth);
};

struct AllPass : BiquadBase
{
  void setup (double sampleRate,
              double phaseFrequency,
              double q);
};

//------------------------------------------------------------------------------

//
// Gui-friendly Design layer
//

namespace Design {

struct TypeIBase : DesignBase
{
  enum
  {
    NumParams = 3
  };

  static int getNumParams ()
  {
    return 3;
  }

  static const ParamInfo getParamInfo_1 ()
  {
    return ParamInfo::defaultCutoffFrequencyParam ();
  }

  static const ParamInfo getParamInfo_2 ()
  {
    return ParamInfo::defaultQParam ();
  }
};

template <class FilterClass>
struct TypeI : TypeIBase, FilterClass
{
  void setParams (const Params& params)
  {
    FilterClass::setup (params[0], params[1], params[2]);
  }
};

struct TypeIIBase : DesignBase
{
  enum
  {
    NumParams = 3
  };

  static int getNumParams ()
  {
    return 3;
  }

  static const ParamInfo getParamInfo_1 ()
  {
    return ParamInfo::defaultCenterFrequencyParam ();
  }

  static const ParamInfo getParamInfo_2 ()
  {
    return ParamInfo::defaultBandwidthParam ();
  }
};

template <class FilterClass>
struct TypeII : TypeIIBase, FilterClass
{
  void setParams (const Params& params)
  {
    FilterClass::setup (params[0], params[1], params[2]);
  }
};

struct TypeIIIBase : DesignBase
{
  enum
  {
    NumParams = 4
  };

  static int getNumParams ()
  {
    return 4;
  }

  static const ParamInfo getParamInfo_1 ()
  {
    return ParamInfo::defaultCutoffFrequencyParam ();
  }

  static const ParamInfo getParamInfo_2 ()
  {
    return ParamInfo::defaultGainParam ();
  }

  static const ParamInfo getParamInfo_3 ()
  {
    return ParamInfo::defaultSlopeParam ();
  }
};

template <class FilterClass>
struct TypeIII : TypeIIIBase, FilterClass
{
  void setParams (const Params& params)
  {
    FilterClass::setup (params[0], params[1], params[2], params[3]);
  }
};

struct TypeIVBase : DesignBase
{
  enum
  {
    NumParams = 4
  };

  static int getNumParams ()
  {
    return 4;
  }

  static const ParamInfo getParamInfo_1 ()
  {
    return ParamInfo::defaultCenterFrequencyParam ();
  }

  static const ParamInfo getParamInfo_2 ()
  {
    return ParamInfo::defaultGainParam ();
  }

  static const ParamInfo getParamInfo_3 ()
  {
    return ParamInfo::defaultBandwidthParam ();
  }
};

template <class FilterClass>
struct TypeIV : TypeIVBase , FilterClass
{
  void setParams (const Params& params)
  {
    FilterClass::setup (params[0], params[1], params[2], params[3]);
  }
};

//------------------------------------------------------------------------------

struct LowPass : TypeI <RBJ::LowPass>
{
  static Kind getKind () { return kindLowPass; }
  static const char* getName() { return "RBJ Low Pass"; }
};

struct HighPass : TypeI <RBJ::HighPass>
{
  static Kind getKind () { return kindHighPass; }
  static const char* getName() { return "RBJ High Pass"; }
};

struct BandPass1 : TypeII <RBJ::BandPass1>
{
  static Kind getKind () { return kindBandPass; }
  static const char* getName() { return "RBJ Band Pass 1"; }
};

struct BandPass2 : TypeII <RBJ::BandPass2>
{
  static Kind getKind () { return kindBandPass; }
  static const char* getName() { return "RBJ Band Pass 2"; }
};

struct BandStop : TypeII <RBJ::BandStop>
{
  static Kind getKind () { return kindBandStop; }
  static const char* getName() { return "RBJ Band Stop"; }
};

struct LowShelf : TypeIII <RBJ::LowShelf>
{
  static Kind getKind () { return kindLowShelf; }
  static const char* getName() { return "RBJ Low Shelf"; }
};

struct HighShelf : TypeIII <RBJ::HighShelf>
{
  static Kind getKind () { return kindHighShelf; }
  static const char* getName() { return "RBJ High Shelf"; }
};
	
struct HighShelfOld : TypeIII <RBJ::HighShelfOld>
{
	static Kind getKind () { return kindHighShelf; }
	static const char* getName() { return "RBJ High Shelf Old version"; }
};

struct BandShelf : TypeIV <RBJ::BandShelf>
{
  static Kind getKind () { return kindBandShelf; }
  static const char* getName() { return "RBJ Band Shelf"; }
};

struct AllPass : TypeI <RBJ::AllPass>
{
  static Kind getKind () { return kindOther; }
  static const char* getName() { return "RBJ All Pass"; }
};

}

}

}

#endif
