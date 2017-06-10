/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modifyf
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 *  2016 Jacques Desmis <jdesmis@gmail.com>
 *  2016 Ingo Weyrich <heckflosse@i-weyrich.de>

 */
#include <cmath>
#include <glib.h>
#include <glibmm.h>

#include "rtengine.h"
#include "improcfun.h"
#include "rawimagesource.h"
#include "../rtgui/thresholdselector.h"
#include "curves.h"
#include "gauss.h"
#include "iccstore.h"
#include "iccmatrices.h"
#include "color.h"
#include "rt_math.h"
#include "jaggedarray.h"
#ifdef _DEBUG
#include "mytime.h"
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include "cplx_wavelet_dec.h"

#define BENCHMARK
#include "StopWatch.h"

#define cliploc( val, minv, maxv )    (( val = (val < minv ? minv : val ) ) > maxv ? maxv : val )

#define CLIPC(a) ((a)>-42000?((a)<42000?(a):42000):-42000)  // limit a and b  to 130 probably enough ?
#define CLIPL(x) LIM(x,0.f,40000.f) // limit L to about L=120 probably enough ?
#define CLIPLOC(x) LIM(x,0.f,32767.f)
#define CLIPLIG(x) LIM(x,0.f, 99.5f)
#define CLIPCHRO(x) LIM(x,0.f, 140.f)
#define CLIPRET(x) LIM(x,-99.5f, 99.5f)

namespace
{

float calcLocalFactor (const float lox, const float loy, const float lcx, const float dx, const float lcy, const float dy, const float ach)
{
//elipse x2/a2 + y2/b2=1
//transition elipsoidal
//x==>lox y==>loy
// a==> dx  b==>dy

    float kelip = dx / dy;
    float belip = sqrt ((rtengine::SQR ((lox - lcx) / kelip) + rtengine::SQR (loy - lcy))); //determine position ellipse ==> a and b
    float aelip = belip * kelip;
    float degrad = aelip / dx;
    float ap = rtengine::RT_PI / (1.f - ach);
    float bp = rtengine::RT_PI - ap;
    return 0.5f * (1.f + xcosf (degrad * ap + bp)); //trigo cos transition

}

}
using namespace std;

namespace rtengine
{
using namespace procparams;

#define SAT(a,b,c) ((float)max(a,b,c)-(float)min(a,b,c))/(float)max(a,b,c)


extern const Settings* settings;

struct local_params {
    float yc, xc;
    float lx, ly;
    float lxL, lyT;
    float dxx, dyy;
    float iterat;
    int cir;
    float thr;
    int prox;
    int chro, cont, sens, sensv;
    float ligh;
    int qualmet;
    float threshol;
    bool exposeena;
    bool expvib;
    bool expwb;
    int blac;
    int shcomp;
    int hlcomp;
    int hlcompthr;
    double temp;
    double  green;
    double equal;
    double expcomp;
    int trans;
    float past;
    float satur;


};

void fillCurveArrayVibloc (DiagonalCurve* diagCurve, LUTf &outCurve)
{

    if (diagCurve) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = 0; i <= 0xffff; i++ ) {
            // change to [0,1] range
            // apply custom/parametric/NURBS curve, if any
            // and store result in a temporary array
            outCurve[i] = 65535.f * diagCurve->getVal ( double (i) / 65535.0 );
        }
    } else {
        outCurve.makeIdentity();
    }
}



static void calcLocalParams (int oW, int oH, const LocalrgbParams& localrgb, struct local_params& lp)
{
    int w = oW;
    int h = oH;
    int circr = localrgb.circrad;

    float thre = localrgb.thres / 100.f;
    double local_x = localrgb.locX / 2000.0;
    double local_y = localrgb.locY / 2000.0;
    double local_xL = localrgb.locXL / 2000.0;
    double local_yT = localrgb.locYT / 2000.0;
    double local_center_x = localrgb.centerX / 2000.0 + 0.5;
    double local_center_y = localrgb.centerY / 2000.0 + 0.5;
    double local_dxx = localrgb.proxi / 8000.0;//for proxi = 2==> # 1 pixel
    double local_dyy = localrgb.proxi / 8000.0;
    float iterati = (float) localrgb.proxi;
    float chromaPastel = float (localrgb.pastels)   / 100.0f;
    float chromaSatur  = float (localrgb.saturated) / 100.0f;

    if (localrgb.qualityMethod == "std") {
        lp.qualmet = 0;
    } else if (localrgb.qualityMethod == "enh") {
        lp.qualmet = 1;
    } else if (localrgb.qualityMethod == "enhden") {
        lp.qualmet = 2;
    }


    int local_chroma = localrgb.chroma;
    int local_sensi = localrgb.sensi;
    int local_contrast = localrgb.contrast;
    float local_lightness = (float) localrgb.lightness;
    int local_transit = localrgb.transit;
    int local_sensiv = localrgb.sensiv;

    lp.cir = circr;
    //  lp.actsp = acti;
    lp.xc = w * local_center_x;
    lp.yc = h * local_center_y;
    lp.lx = w * local_x;
    lp.ly = h * local_y;
    lp.lxL = w * local_xL;
    lp.lyT = h * local_yT;
    lp.chro = local_chroma;
    lp.sens = local_sensi;
    lp.cont = local_contrast;
    lp.ligh = local_lightness;

    if (lp.ligh >= -2.f && lp.ligh <= 2.f) {
        lp.ligh /= 5.f;
    }

    lp.trans = local_transit;
    lp.iterat = iterati;
    lp.dxx = w * local_dxx;
    lp.dyy = h * local_dyy;
    lp.thr = thre;

    lp.exposeena = localrgb.expexpose;
    lp.expvib = localrgb.expvibrance;
    lp.expwb = localrgb.expwb;
    lp.blac = localrgb.black;
    lp.shcomp = localrgb.shcompr;
    lp.hlcomp = localrgb.hlcompr;
    lp.hlcompthr = localrgb.hlcomprthresh;
    lp.temp = localrgb.temp;
    lp.green = localrgb.green;
    lp.equal = localrgb.equal;
    lp.expcomp = localrgb.expcomp;
    lp.sensv = local_sensiv;
    lp.past =  chromaPastel;
    lp.satur = chromaSatur;
}

static void calcTransition (const float lox, const float loy, const float ach, const local_params& lp, int &zone, float &localFactor)
{
    // returns the zone (0 = outside selection, 1 = transition zone between outside and inside selection, 2 = inside selection)
    // and a factor to calculate the transition in case zone == 1

    zone = 0;

    if (lox >= lp.xc && lox < (lp.xc + lp.lx) && loy >= lp.yc && loy < lp.yc + lp.ly) {
        float zoneVal = SQR ((lox - lp.xc) / (ach * lp.lx)) + SQR ((loy - lp.yc) / (ach * lp.ly));
        zone = zoneVal < 1.f ? 2 : 0;

        if (!zone) {
            zone = (zoneVal > 1.f && ((SQR ((lox - lp.xc) / (lp.lx)) + SQR ((loy - lp.yc) / (lp.ly))) < 1.f)) ? 1 : 0;

            if (zone) {
                localFactor = calcLocalFactor (lox, loy, lp.xc, lp.lx, lp.yc, lp.ly, ach);
            }
        }
    } else if (lox >= lp.xc && lox < lp.xc + lp.lx && loy < lp.yc && loy > lp.yc - lp.lyT) {
        float zoneVal = SQR ((lox - lp.xc) / (ach * lp.lx)) + SQR ((loy - lp.yc) / (ach * lp.lyT));
        zone = zoneVal < 1.f ? 2 : 0;

        if (!zone) {
            zone = (zoneVal > 1.f && ((SQR ((lox - lp.xc) / (lp.lx)) + SQR ((loy - lp.yc) / (lp.lyT))) < 1.f)) ? 1 : 0;

            if (zone) {
                localFactor = calcLocalFactor (lox, loy, lp.xc, lp.lx, lp.yc, lp.lyT, ach);
            }
        }
    } else if (lox < lp.xc && lox > lp.xc - lp.lxL && loy <= lp.yc && loy > lp.yc - lp.lyT) {
        float zoneVal = SQR ((lox - lp.xc) / (ach * lp.lxL)) + SQR ((loy - lp.yc) / (ach * lp.lyT));
        zone = zoneVal < 1.f ? 2 : 0;

        if (!zone) {
            zone = (zoneVal > 1.f && ((SQR ((lox - lp.xc) / (lp.lxL)) + SQR ((loy - lp.yc) / (lp.lyT))) < 1.f)) ? 1 : 0;

            if (zone) {
                localFactor = calcLocalFactor (lox, loy, lp.xc, lp.lxL, lp.yc, lp.lyT, ach);
            }
        }
    } else if (lox < lp.xc && lox > lp.xc - lp.lxL && loy > lp.yc && loy < lp.yc + lp.ly) {
        float zoneVal = SQR ((lox - lp.xc) / (ach * lp.lxL)) + SQR ((loy - lp.yc) / (ach * lp.ly));
        zone = zoneVal < 1.f ? 2 : 0;

        if (!zone) {
            zone = (zoneVal > 1.f && ((SQR ((lox - lp.xc) / (lp.lxL)) + SQR ((loy - lp.yc) / (lp.ly))) < 1.f)) ? 1 : 0;

            if (zone) {
                localFactor = calcLocalFactor (lox, loy, lp.xc, lp.lxL, lp.yc, lp.ly, ach);
            }
        }
    }
}
void ImProcFunctions::calcrgb_ref (LabImage * original, LabImage * transformed, int sk, int sx, int sy, int cx, int cy, int oW, int oH,  int fw, int fh, double & hueref, double & chromaref, double & lumaref)
{
    if (params->localrgb.enabled) {
        //always calculate hueref, chromaref, lumaref  before others operations use in normal mode for all modules exceprt denoise
        struct local_params lp;
        calcLocalParams (oW, oH, params->localrgb, lp);
        //  printf("OK calcrgb\n");
        //int sk = 1;
// double precision for large summations
        double aveA = 0.;
        double aveB = 0.;
        double aveL = 0.;
        double aveChro = 0.;
// int precision for the counters
        int nab = 0;
// single precision for the result
        float avA, avB, avL;
        int spotSize = 0.88623f * max (1,  lp.cir / sk); //18

        //O.88623 = sqrt(PI / 4) ==> sqare equal to circle
        //  printf("Spots=%i  xc=%i\n", spotSize, (int) lp.xc);
        // very small region, don't use omp here
        for (int y = max (cy, (int) (lp.yc - spotSize)); y < min (transformed->H + cy, (int) (lp.yc + spotSize + 1)); y++) {
            for (int x = max (cx, (int) (lp.xc - spotSize)); x < min (transformed->W + cx, (int) (lp.xc + spotSize + 1)); x++) {
                aveL += original->L[y - cy][x - cx];
                aveA += original->a[y - cy][x - cx];
                aveB += original->b[y - cy][x - cx];
                aveChro += sqrtf (SQR (original->b[y - cy][x - cx]) + SQR (original->a[y - cy][x - cx]));

                nab++;
            }
        }

        aveL = aveL / nab;
        aveA = aveA / nab;
        aveB = aveB / nab;
        aveChro = aveChro / nab;
        aveChro /= 327.68f;
        avA = aveA / 327.68f;
        avB = aveB / 327.68f;
        avL = aveL / 327.68f;
        hueref = xatan2f (avB, avA);   //mean hue
        chromaref = aveChro;
        lumaref = avL;
    }
}





void ImProcFunctions::WB_Local (ImageSource* imgsrc, int call, int sp, int sx, int sy, int cx, int cy, int oW, int oH,  int fw, int fh, Imagefloat* improv, Imagefloat* imagetransformed, const ColorTemp &ctemploc, int tran, Imagefloat* imageoriginal, const PreviewProps &pp, const ToneCurveParams &hrp, const ColorManagementParams &cmp, const RAWParams &raw, double &ptemp, double &pgreen)
{
    if (params->localrgb.enabled) {
        // BENCHFUN
#ifdef _DEBUG
        MyTime t1e, t2e;
        t1e.set();
// init variables to display Munsell corrections
        MunsellDebugInfo* MunsDebugInfo = new MunsellDebugInfo();
#endif



        struct local_params lp;
        calcLocalParams (oW, oH, params->localrgb, lp);

        double ave = 0.;
        int n = 0;
        float av = 0;
        constexpr float ared = (rtengine::RT_PI - 0.05f) / 100.f;

        constexpr float bred = 0.05f;

        float dhue = ared * lp.sens + bred; //delta hue lght chroma

        Imagefloat* bufimage = nullptr;

        int bfh = 0.f, bfw = 0.f;
        int del = 3; // to avoid crash

        int begy = lp.yc - lp.lyT;
        int begx = lp.xc - lp.lxL;
        int yEn = lp.yc + lp.ly;
        int xEn = lp.xc + lp.lx;
        //  printf ("by=%i bx=%i y=%i x=%i\n", begx, begy, yEn, xEn);

        if (call <= 3) { //simpleprocess, dcrop, improccoordinator
            bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
            bfw = int (lp.lx + lp.lxL) + del;

            bufimage = new Imagefloat (bfw, bfh);
#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int ir = 0; ir < bfh; ir++) //fill with 0
                for (int jr = 0; jr < bfw; jr++) {
                    bufimage->r (ir, jr) = 0.f;
                    bufimage->g (ir, jr) = 0.f;
                    bufimage->b (ir, jr) = 0.f;
                }

            if (lp.expwb && params->localrgb.wbMethod != "none") {
                /*                float hueplus = hueref + dhue;
                                float huemoins = hueref - dhue;

                                //printf("hueplus=%f huemoins=%f dhu=%f\n", hueplus, huemoins, dhue);
                                if (hueplus > rtengine::RT_PI) {
                                    hueplus = hueref + dhue - 2.f * rtengine::RT_PI;
                                }

                                if (huemoins < -rtengine::RT_PI) {
                                    huemoins = hueref - dhue + 2.f * rtengine::RT_PI;
                                }
                */
                imgsrc->getImage_local (begx, begy, yEn, xEn, cx, cy, ctemploc, tran, improv, bufimage, pp, hrp, cmp, raw);

                Whitebalance_Local (call, sp, bufimage, lp, imageoriginal, imagetransformed, cx, cy);
                /*
                                #ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16)
                #endif

                                for (int y = 0; y < imageoriginal->getHeight() ; y++) //{
                                    for (int x = 0; x < imageoriginal->getWidth(); x++) {
                                        int lox = cx + x;
                                        int loy = cy + y;

                                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                                            imageoriginal->r (y, x) =  bufimage->r (loy - begy, lox - begx);
                                            imageoriginal->g (y, x) =  bufimage->g (loy - begy, lox - begx);
                                            imageoriginal->b (y, x) =  bufimage->b (loy - begy, lox - begx);
                                        }
                                    }
                */

            }

            delete bufimage;

        }
    }
}


void ImProcFunctions::vibrancelocal (const local_params& lp, int bfw, int bfh, LabImage* lab,  LabImage* dest)
{
    if (!params->localrgb.expvibrance) {
        return;
    }

//  int skip=1; //scale==1 ? 1 : 16;
    bool skinCurveIsSet = false;
    DiagonalCurve* dcurve = nullptr;
    dcurve = new DiagonalCurve (params->localrgb.skintonescurve, CURVES_MIN_POLY_POINTS);

    if (dcurve) {
        if (!dcurve->isIdentity()) {
            skinCurveIsSet = true;
        } else {
            delete dcurve;
            dcurve = nullptr;
        }
    }

    if (!skinCurveIsSet && !params->localrgb.pastels && !params->localrgb.saturated) {
        if (dcurve) {
            delete dcurve;
            dcurve = nullptr;
        }

        return;
    }

    const int width = bfw;
    const int height = bfh;

#ifdef _DEBUG
    MyTime t1e, t2e;
    t1e.set();
    int negat = 0, moreRGB = 0, negsat = 0, moresat = 0;
#endif

    // skin hue curve
    // I use diagonal because I think it's better
    LUTf skin_curve (65536, 0);

    if (skinCurveIsSet) {
        fillCurveArrayVibloc (dcurve, skin_curve);
    }

    if (dcurve) {
        delete dcurve;
        dcurve = nullptr;
    }



    const float chromaPastel = float (params->localrgb.pastels)   / 100.0f;
    const float chromaSatur  = float (params->localrgb.saturated) / 100.0f;
    const float p00 = 0.07f;
    const float limitpastelsatur =    (static_cast<float> (params->localrgb.psthreshold.value[ThresholdSelector::TS_TOPLEFT])    / 100.0f) * (1.0f - p00) + p00;
    const float maxdp = (limitpastelsatur - p00) / 4.0f;
    const float maxds = (1.0 - limitpastelsatur) / 4.0f;
    const float p0 = p00 + maxdp;
    const float p1 = p00 + 2.0f * maxdp;
    const float p2 = p00 + 3.0f * maxdp;
    const float s0 = limitpastelsatur + maxds;
    const float s1 = limitpastelsatur + 2.0f * maxds;
    const float s2 = limitpastelsatur + 3.0f * maxds;
    const float transitionweighting = static_cast<float> (params->localrgb.psthreshold.value[ThresholdSelector::TS_BOTTOMLEFT]) / 100.0f;
    float chromamean = 0.0f;

    if (chromaPastel != chromaSatur) {
        //if sliders pastels and saturated are different: transition with a double linear interpolation: between p2 and limitpastelsatur, and between limitpastelsatur and s0
        //modify the "mean" point in function of double threshold  => differential transition
        chromamean = maxdp * (chromaSatur - chromaPastel) / (s0 - p2) + chromaPastel;

        // move chromaMean up or down depending on transitionCtrl
        if (transitionweighting > 0.0f) {
            chromamean = (chromaSatur - chromamean) * transitionweighting + chromamean;
        } else if (transitionweighting < 0.0f) {
            chromamean = (chromamean - chromaPastel)  * transitionweighting + chromamean;
        }
    }

    const float chromaPastel_a = (chromaPastel - chromamean) / (p2 - limitpastelsatur);
    const float chromaPastel_b = chromaPastel - chromaPastel_a * p2;

    const float chromaSatur_a = (chromaSatur - chromamean) / (s0 - limitpastelsatur);
    const float chromaSatur_b = chromaSatur - chromaSatur_a * s0;

    const float dhue = 0.15f; //hue transition
    const float dchr = 20.0f; //chroma transition
    const float skbeg = -0.05f; //begin hue skin
    const float skend = 1.60f; //end hue skin
    const float xx = 0.5f; //soft : between 0.3 and 1.0
    const float ask = 65535.0f / (skend - skbeg);
    const float bsk = -skbeg * ask;


    const bool highlight = params->toneCurve.hrenabled;//Get the value if "highlight reconstruction" is activated
    const bool protectskins = params->localrgb.protectskins;
    const bool avoidcolorshift = params->localrgb.avoidcolorshift;

    TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix (params->icm.working);
    //inverse matrix user select
    const double wip[3][3] = {
        {wiprof[0][0], wiprof[0][1], wiprof[0][2]},
        {wiprof[1][0], wiprof[1][1], wiprof[1][2]},
        {wiprof[2][0], wiprof[2][1], wiprof[2][2]}
    };


#ifdef _DEBUG
    MunsellDebugInfo* MunsDebugInfo = nullptr;

    if (avoidcolorshift) {
        MunsDebugInfo = new MunsellDebugInfo();
    }

    #pragma omp parallel default(shared) firstprivate(lab, dest, MunsDebugInfo) reduction(+: negat, moreRGB, negsat, moresat) if (multiThread)
#else
    #pragma omp parallel default(shared) if (multiThread)
#endif
    {

        float sathue[5], sathue2[4]; // adjust sat in function of hue

        /*
            // Fitting limitpastelsatur into the real 0.07->1.0 range
        //  limitpastelsatur = limitpastelsatur*(1.0f-p00) + p00;
            float p0,p1,p2;//adapt limit of pyramid to psThreshold
            float s0,s1,s2;
        */

#ifdef _OPENMP

        if (settings->verbose && omp_get_thread_num() == 0) {
#else

        if (settings->verbose) {
#endif
            printf ("vibrance:  p0=%1.2f  p1=%1.2f  p2=%1.2f  s0=%1.2f s1=%1.2f s2=%1.2f\n", p0, p1, p2, s0, s1, s2);
            printf ("           pastel=%f   satur=%f   limit= %1.2f   chromamean=%0.5f\n", 1.0f + chromaPastel, 1.0f + chromaSatur, limitpastelsatur, chromamean);
        }

        #pragma omp for schedule(dynamic, 16)

        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++) {
                float LL = lab->L[i][j] / 327.68f;
                float CC = sqrt (SQR (lab->a[i][j]) + SQR (lab->b[i][j])) / 327.68f;
                float HH = xatan2f (lab->b[i][j], lab->a[i][j]);

                float satredu = 1.0f; //reduct sat in function of skin

                if (protectskins) {
                    Color::SkinSat (LL, HH, CC, satredu);// for skin colors
                }

                // here we work on Chromaticity and Hue
                // variation of Chromaticity  ==> saturation via RGB
                // Munsell correction, then conversion to Lab
                float Lprov = LL;
                float Chprov = CC;
                float R, G, B;
                float2 sincosval;

                if (CC == 0.0f) {
                    sincosval.y = 1.f;
                    sincosval.x = 0.0f;
                } else {
                    sincosval.y = lab->a[i][j] / (CC * 327.68f);
                    sincosval.x = lab->b[i][j] / (CC * 327.68f);
                }

#ifdef _DEBUG
                bool neg = false;
                bool more_rgb = false;
                //gamut control : Lab values are in gamut
                Color::gamutLchonly (HH, sincosval, Lprov, Chprov, R, G, B, wip, highlight, 0.15f, 0.98f, neg, more_rgb);

                if (neg) {
                    negat++;
                }

                if (more_rgb) {
                    moreRGB++;
                }

#else
                //gamut control : Lab values are in gamut
                Color::gamutLchonly (HH, sincosval, Lprov, Chprov, R, G, B, wip, highlight, 0.15f, 0.98f);
#endif

                if (Chprov > 6.0f) {
                    const float saturation = SAT (R, G, B);

                    if (saturation > 0.0f) {
                        if (satredu != 1.0f) {
                            // for skin, no differentiation
                            sathue [0] = sathue [1] = sathue [2] = sathue [3] = sathue[4] = 1.0f;
                            sathue2[0] = sathue2[1] = sathue2[2] = sathue2[3]          = 1.0f;
                        } else {
                            //double pyramid: LL and HH
                            //I try to take into account: Munsell response (human vision) and Gamut..(less response for red): preferably using Prophoto or WideGamut
                            //blue: -1.80 -3.14  green = 2.1 3.14   green-yellow=1.4 2.1  red:0 1.4  blue-purple:-0.7  -1.4   purple: 0 -0.7
                            //these values allow a better and differential response
                            if (LL < 20.0f) { //more for blue-purple, blue and red modulate
                                if     (/*HH> -3.1415f &&*/ HH < -1.5f   ) {
                                    sathue[0] = 1.3f;    //blue
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.1f;
                                    sathue[3] = 1.05f;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.05f;
                                    sathue2[1] = 1.1f ;
                                    sathue2[2] = 1.05f;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>=-1.5f    &&*/ HH < -0.7f   ) {
                                    sathue[0] = 1.6f;    //blue purple  1.2 1.1
                                    sathue[1] = 1.4f;
                                    sathue[2] = 1.3f;
                                    sathue[3] = 1.2f ;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.2f ;
                                    sathue2[1] = 1.15f;
                                    sathue2[2] = 1.1f ;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>=-0.7f    &&*/ HH <  0.0f   ) {
                                    sathue[0] = 1.2f;    //purple
                                    sathue[1] = 1.0f;
                                    sathue[2] = 1.0f;
                                    sathue[3] = 1.0f ;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 1.0f ;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                }
                                //          else if(  HH>= 0.0f    &&   HH<= 1.4f   ) {sathue[0]=1.1f;sathue[1]=1.1f;sathue[2]=1.1f;sathue[3]=1.0f ;sathue[4]=0.4f;sathue2[0]=1.0f ;sathue2[1]=1.0f ;sathue2[2]=1.0f ;sathue2[3]=1.0f;}//red   0.8 0.7
                                else if (/*HH>= 0.0f    &&*/ HH <= 1.4f   ) {
                                    sathue[0] = 1.3f;    //red   0.8 0.7
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.1f;
                                    sathue[3] = 1.0f ;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 1.0f ;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>  1.4f    &&*/ HH <= 2.1f   ) {
                                    sathue[0] = 1.0f;    //green yellow 1.2 1.1
                                    sathue[1] = 1.0f;
                                    sathue[2] = 1.0f;
                                    sathue[3] = 1.0f ;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 1.0f ;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                } else { /*if(HH>  2.1f    && HH<= 3.1415f)*/
                                    sathue[0] = 1.4f;    //green
                                    sathue[1] = 1.3f;
                                    sathue[2] = 1.2f;
                                    sathue[3] = 1.15f;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.15f;
                                    sathue2[1] = 1.1f ;
                                    sathue2[2] = 1.05f;
                                    sathue2[3] = 1.0f;
                                }
                            } else if (LL < 50.0f) { //more for blue and green, less for red and green-yellow
                                if     (/*HH> -3.1415f &&*/ HH < -1.5f   ) {
                                    sathue[0] = 1.5f;    //blue
                                    sathue[1] = 1.4f;
                                    sathue[2] = 1.3f;
                                    sathue[3] = 1.2f ;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.2f ;
                                    sathue2[1] = 1.1f ;
                                    sathue2[2] = 1.05f;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>=-1.5f    &&*/ HH < -0.7f   ) {
                                    sathue[0] = 1.3f;    //blue purple  1.2 1.1
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.1f;
                                    sathue[3] = 1.05f;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.05f;
                                    sathue2[1] = 1.05f;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>=-0.7f    &&*/ HH <  0.0f   ) {
                                    sathue[0] = 1.2f;    //purple
                                    sathue[1] = 1.0f;
                                    sathue[2] = 1.0f;
                                    sathue[3] = 1.0f ;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 1.0f ;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                }
                                //          else if(  HH>= 0.0f    &&   HH<= 1.4f   ) {sathue[0]=0.8f;sathue[1]=0.8f;sathue[2]=0.8f;sathue[3]=0.8f ;sathue[4]=0.4f;sathue2[0]=0.8f ;sathue2[1]=0.8f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//red   0.8 0.7
                                else if (/*HH>= 0.0f    &&*/ HH <= 1.4f   ) {
                                    sathue[0] = 1.1f;    //red   0.8 0.7
                                    sathue[1] = 1.0f;
                                    sathue[2] = 0.9f;
                                    sathue[3] = 0.8f ;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 0.8f ;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                } else if (/*HH>  1.4f    &&*/ HH <= 2.1f   ) {
                                    sathue[0] = 1.1f;    //green yellow 1.2 1.1
                                    sathue[1] = 1.1f;
                                    sathue[2] = 1.1f;
                                    sathue[3] = 1.05f;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 0.9f ;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.7f ;
                                    sathue2[3] = 0.6f;
                                } else { /*if(HH>  2.1f    && HH<= 3.1415f)*/
                                    sathue[0] = 1.5f;    //green
                                    sathue[1] = 1.4f;
                                    sathue[2] = 1.3f;
                                    sathue[3] = 1.2f ;
                                    sathue[4] = 0.4f;
                                    sathue2[0] = 1.2f ;
                                    sathue2[1] = 1.1f ;
                                    sathue2[2] = 1.05f;
                                    sathue2[3] = 1.0f;
                                }

                            } else if (LL < 80.0f) { //more for green, less for red and green-yellow
                                if     (/*HH> -3.1415f &&*/ HH < -1.5f   ) {
                                    sathue[0] = 1.3f;    //blue
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.15f;
                                    sathue[3] = 1.1f ;
                                    sathue[4] = 0.3f;
                                    sathue2[0] = 1.1f ;
                                    sathue2[1] = 1.1f ;
                                    sathue2[2] = 1.05f;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>=-1.5f    &&*/ HH < -0.7f   ) {
                                    sathue[0] = 1.3f;    //blue purple  1.2 1.1
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.15f;
                                    sathue[3] = 1.1f ;
                                    sathue[4] = 0.3f;
                                    sathue2[0] = 1.1f ;
                                    sathue2[1] = 1.05f;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>=-0.7f    &&*/ HH <  0.0f   ) {
                                    sathue[0] = 1.2f;    //purple
                                    sathue[1] = 1.0f;
                                    sathue[2] = 1.0f ;
                                    sathue[3] = 1.0f ;
                                    sathue[4] = 0.3f;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 1.0f ;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                }
                                //          else if(  HH>= 0.0f    &&   HH<= 1.4f   ) {sathue[0]=0.8f;sathue[1]=0.8f;sathue[2]=0.8f ;sathue[3]=0.8f ;sathue[4]=0.3f;sathue2[0]=0.8f ;sathue2[1]=0.8f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//red   0.8 0.7
                                else if (/*HH>= 0.0f    &&*/ HH <= 1.4f   ) {
                                    sathue[0] = 1.1f;    //red   0.8 0.7
                                    sathue[1] = 1.0f;
                                    sathue[2] = 0.9f ;
                                    sathue[3] = 0.8f ;
                                    sathue[4] = 0.3f;
                                    sathue2[0] = 0.8f ;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                } else if (/*HH>  1.4f    &&*/ HH <= 2.1f   ) {
                                    sathue[0] = 1.3f;    //green yellow 1.2 1.1
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.1f ;
                                    sathue[3] = 1.05f;
                                    sathue[4] = 0.3f;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 0.9f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.7f;
                                } else { /*if(HH>  2.1f    && HH<= 3.1415f)*/
                                    sathue[0] = 1.6f;    //green - even with Prophoto green are too "little"  1.5 1.3
                                    sathue[1] = 1.4f;
                                    sathue[2] = 1.3f ;
                                    sathue[3] = 1.25f;
                                    sathue[4] = 0.3f;
                                    sathue2[0] = 1.25f;
                                    sathue2[1] = 1.2f ;
                                    sathue2[2] = 1.15f;
                                    sathue2[3] = 1.05f;
                                }
                            } else { /*if (LL>=80.0f)*/ //more for green-yellow, less for red and purple
                                if     (/*HH> -3.1415f &&*/ HH < -1.5f   ) {
                                    sathue[0] = 1.0f;    //blue
                                    sathue[1] = 1.0f;
                                    sathue[2] = 0.9f;
                                    sathue[3] = 0.8f;
                                    sathue[4] = 0.2f;
                                    sathue2[0] = 0.8f;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                } else if (/*HH>=-1.5f    &&*/ HH < -0.7f   ) {
                                    sathue[0] = 1.0f;    //blue purple  1.2 1.1
                                    sathue[1] = 1.0f;
                                    sathue[2] = 0.9f;
                                    sathue[3] = 0.8f;
                                    sathue[4] = 0.2f;
                                    sathue2[0] = 0.8f;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                } else if (/*HH>=-0.7f    &&*/ HH <  0.0f   ) {
                                    sathue[0] = 1.2f;    //purple
                                    sathue[1] = 1.0f;
                                    sathue[2] = 1.0f;
                                    sathue[3] = 0.9f;
                                    sathue[4] = 0.2f;
                                    sathue2[0] = 0.9f;
                                    sathue2[1] = 0.9f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                }
                                //          else if(  HH>= 0.0f    &&   HH<= 1.4f   ) {sathue[0]=0.8f;sathue[1]=0.8f;sathue[2]=0.8f;sathue[3]=0.8f;sathue[4]=0.2f;sathue2[0]=0.8f;sathue2[1]=0.8f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//red   0.8 0.7
                                else if (/*HH>= 0.0f    &&*/ HH <= 1.4f   ) {
                                    sathue[0] = 1.1f;    //red   0.8 0.7
                                    sathue[1] = 1.0f;
                                    sathue[2] = 0.9f;
                                    sathue[3] = 0.8f;
                                    sathue[4] = 0.2f;
                                    sathue2[0] = 0.8f;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                } else if (/*HH>  1.4f    &&*/ HH <= 2.1f   ) {
                                    sathue[0] = 1.6f;    //green yellow 1.2 1.1
                                    sathue[1] = 1.5f;
                                    sathue[2] = 1.4f;
                                    sathue[3] = 1.2f;
                                    sathue[4] = 0.2f;
                                    sathue2[0] = 1.1f;
                                    sathue2[1] = 1.05f;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                } else { /*if(HH>  2.1f    && HH<= 3.1415f)*/
                                    sathue[0] = 1.4f;    //green
                                    sathue[1] = 1.3f;
                                    sathue[2] = 1.2f;
                                    sathue[3] = 1.1f;
                                    sathue[4] = 0.2f;
                                    sathue2[0] = 1.1f;
                                    sathue2[1] = 1.05f;
                                    sathue2[2] = 1.05f;
                                    sathue2[3] = 1.0f;
                                }
                            }
                        }

                        float chmodpastel = 0.f, chmodsat = 0.f;
                        // variables to improve transitions
                        float pa, pb;// transition = pa*saturation + pb
                        float chl00 = chromaPastel * satredu * sathue[4];
                        float chl0  = chromaPastel * satredu * sathue[0];
                        float chl1  = chromaPastel * satredu * sathue[1];
                        float chl2  = chromaPastel * satredu * sathue[2];
                        float chl3  = chromaPastel * satredu * sathue[3];
                        float chs0  = chromaSatur * satredu * sathue2[0];
                        float chs1  = chromaSatur * satredu * sathue2[1];
                        float chs2  = chromaSatur * satredu * sathue2[2];
                        float chs3  = chromaSatur * satredu * sathue2[3];
                        float s3    = 1.0f;

                        // We handle only positive values here ;  improve transitions
                        if      (saturation < p00) {
                            chmodpastel = chl00 ;    //neutral tones
                        } else if (saturation < p0 )               {
                            pa = (chl00 - chl0) / (p00 - p0);
                            pb = chl00 - pa * p00;
                            chmodpastel = pa * saturation + pb;
                        } else if (saturation < p1)                {
                            pa = (chl0 - chl1) / (p0 - p1);
                            pb = chl0 - pa * p0;
                            chmodpastel = pa * saturation + pb;
                        } else if (saturation < p2)                {
                            pa = (chl1 - chl2) / (p1 - p2);
                            pb = chl1 - pa * p1;
                            chmodpastel = pa * saturation + pb;
                        } else if (saturation < limitpastelsatur)  {
                            pa = (chl2 - chl3) / (p2 - limitpastelsatur);
                            pb = chl2 - pa * p2;
                            chmodpastel = pa * saturation + pb;
                        } else if (saturation < s0)                {
                            pa = (chl3 - chs0) / (limitpastelsatur - s0) ;
                            pb = chl3 - pa * limitpastelsatur;
                            chmodsat    = pa * saturation + pb;
                        } else if (saturation < s1)                {
                            pa = (chs0 - chs1) / (s0 - s1);
                            pb = chs0 - pa * s0;
                            chmodsat    = pa * saturation + pb;
                        } else if (saturation < s2)                {
                            pa = (chs1 - chs2) / (s1 - s2);
                            pb = chs1 - pa * s1;
                            chmodsat    = pa * saturation + pb;
                        } else                                     {
                            pa = (chs2 - chs3) / (s2 - s3);
                            pb = chs2 - pa * s2;
                            chmodsat    = pa * saturation + pb;
                        }

                        if (chromaPastel != chromaSatur) {

                            // Pastels
                            if (saturation > p2 && saturation < limitpastelsatur) {
                                float newchromaPastel = chromaPastel_a * saturation + chromaPastel_b;
                                chmodpastel = newchromaPastel * satredu * sathue[3];
                            }

                            // Saturated
                            if (saturation < s0 && saturation >= limitpastelsatur) {
                                float newchromaSatur = chromaSatur_a * saturation + chromaSatur_b;
                                chmodsat = newchromaSatur * satredu * sathue2[0];
                            }
                        }// end transition

                        if (saturation <= limitpastelsatur) {
                            if (chmodpastel >  2.0f ) {
                                chmodpastel = 2.0f;    //avoid too big values
                            } else if (chmodpastel < -0.93f) {
                                chmodpastel = -0.93f;    //avoid negative values
                            }

                            Chprov *= (1.0f + chmodpastel);

                            if (Chprov < 6.0f) {
                                Chprov = 6.0f;
                            }
                        } else { //if (saturation > limitpastelsatur)
                            if (chmodsat >  1.8f ) {
                                chmodsat = 1.8f;    //saturated
                            } else if (chmodsat < -0.93f) {
                                chmodsat = -0.93f;
                            }

                            Chprov *= 1.0f + chmodsat;

                            if (Chprov < 6.0f) {
                                Chprov = 6.0f;
                            }
                        }
                    }
                }

                bool hhModified = false;

                // Vibrance's Skin curve
                if (skinCurveIsSet) {
                    if (HH > skbeg && HH < skend) {
                        if (Chprov < 60.0f) { //skin hue  : todo ==> transition
                            float HHsk = ask * HH + bsk;
                            float Hn = (skin_curve[HHsk] - bsk) / ask;
                            float Hc = (Hn * xx + HH * (1.0f - xx));
                            HH = Hc;
                            hhModified = true;
                        } else if (Chprov < (60.0f + dchr)) { //transition chroma
                            float HHsk = ask * HH + bsk;
                            float Hn = (skin_curve[HHsk] - bsk) / ask;
                            float Hc = (Hn * xx + HH * (1.0f - xx));
                            float aa = (HH - Hc) / dchr ;
                            float bb = HH - (60.0f + dchr) * aa;
                            HH = aa * Chprov + bb;
                            hhModified = true;
                        }
                    }
                    //transition hue
                    else if (HH > (skbeg - dhue) && HH <= skbeg && Chprov < (60.0f + dchr * 0.5f)) {
                        float HHsk = ask * skbeg + bsk;
                        float Hn = (skin_curve[HHsk] - bsk) / ask;
                        float Hcc = (Hn * xx + skbeg * (1.0f - xx));
                        float adh = (Hcc - (skbeg - dhue)) / (dhue);
                        float bdh = Hcc - adh * skbeg;
                        HH = adh * HH + bdh;
                        hhModified = true;
                    } else if (HH >= skend && HH < (skend + dhue) && Chprov < (60.0f + dchr * 0.5f)) {
                        float HHsk = ask * skend + bsk;
                        float Hn = (skin_curve[HHsk] - bsk) / ask;
                        float Hcc = (Hn * xx + skend * (1.0f - xx));
                        float adh = (skend + dhue - Hcc) / (dhue);
                        float bdh = Hcc - adh * skend;
                        HH = adh * HH + bdh;
                        hhModified = true;
                    }
                } // end skin hue

                //Munsell correction
//          float2 sincosval;
                if (!avoidcolorshift && hhModified) {
                    sincosval = xsincosf (HH);
                }

                float aprovn, bprovn;
                bool inGamut;

                do {
                    inGamut = true;

                    if (avoidcolorshift) {
                        float correctionHue = 0.0f;
                        float correctlum = 0.0f;

#ifdef _DEBUG
                        Color::AllMunsellLch (/*lumaMuns*/false, Lprov, Lprov, HH, Chprov, CC, correctionHue, correctlum, MunsDebugInfo);
#else
                        Color::AllMunsellLch (/*lumaMuns*/false, Lprov, Lprov, HH, Chprov, CC, correctionHue, correctlum);
#endif

                        if (correctionHue != 0.f || hhModified) {
                            sincosval = xsincosf (HH + correctionHue);
                            hhModified = false;
                        }
                    }

                    aprovn = Chprov * sincosval.y;
                    bprovn = Chprov * sincosval.x;

                    float fyy = (0.00862069f * Lprov ) + 0.137932f;
                    float fxx = (0.002f * aprovn) + fyy;
                    float fzz = fyy - (0.005f * bprovn);
                    float xx_ = 65535.f * Color::f2xyz (fxx) * Color::D50x;
                    //  float yy_ = 65535.0f * Color::f2xyz(fyy);
                    float zz_ = 65535.f * Color::f2xyz (fzz) * Color::D50z;
                    float yy_ = 65535.f * ((Lprov > Color::epskap) ? fyy * fyy*fyy : Lprov / Color::kappa);

                    Color::xyz2rgb (xx_, yy_, zz_, R, G, B, wip);

                    if (R < 0.0f || G < 0.0f || B < 0.0f) {
#ifdef _DEBUG
                        negsat++;
#endif
                        Chprov *= 0.98f;
                        inGamut = false;
                    }

                    // if "highlight reconstruction" enabled don't control Gamut for highlights
                    if ((!highlight) && (R > 65535.0f || G > 65535.0f || B > 65535.0f)) {
#ifdef _DEBUG
                        moresat++;
#endif
                        Chprov *= 0.98f;
                        inGamut = false;
                    }
                } while (!inGamut);

                //put new values in Lab
                dest->L[i][j] = Lprov * 327.68f;
                dest->a[i][j] = aprovn * 327.68f;
                dest->b[i][j] = bprovn * 327.68f;
            }
    } // end of parallelization

#ifdef _DEBUG
    t2e.set();

    if (settings->verbose) {
        printf ("Vibrance local (performed in %d usec):\n", t2e.etime (t1e));
        printf ("   Gamut: G1negat=%iiter G165535=%iiter G2negsat=%iiter G265535=%iiter\n", negat, moreRGB, negsat, moresat);

        if (MunsDebugInfo) {
            printf ("   Munsell chrominance: MaxBP=%1.2frad  MaxRY=%1.2frad  MaxGY=%1.2frad  MaxRP=%1.2frad  depass=%u\n", MunsDebugInfo->maxdhue[0], MunsDebugInfo->maxdhue[1], MunsDebugInfo->maxdhue[2], MunsDebugInfo->maxdhue[3], MunsDebugInfo->depass);
        }
    }

    if (MunsDebugInfo) {
        delete MunsDebugInfo;
    }

#endif

}



void ImProcFunctions::rgblabLocal (const local_params& lp, Imagefloat* working, int bfh, int bfw, LabImage* bufexporig, LabImage* lab, Imagefloat* orirgb,  LUTf & hltonecurve, LUTf & shtonecurve, LUTf & tonecurve,
                                   const ToneCurve & customToneCurve1, const ToneCurve & customToneCurve2, DCPProfile *dcpProf, const DCPProfile::ApplyState &asIn )
{

    const float exp_scale = pow (2.0, lp.expcomp);//lp.expcomp
    const float comp = (max (0.0, lp.expcomp) + 1.0) * lp.hlcomp / 100.0;
    const float shoulder = ((65536.0 / max (1.0f, exp_scale)) * (lp.hlcompthr / 200.0)) + 0.1;
    const float hlrange = 65536.0 - shoulder;
    const bool isProPhoto = (params->icm.working == "ProPhoto");
    // extracting datas from 'params' to avoid cache flush (to be confirmed)
    LocalrgbParams::eTCModeId curveMode = params->localrgb.curveMode;
    LocalrgbParams::eTCModeId curveMode2 = params->localrgb.curveMode2;
    //  bool highlight = params->toneCurve.hrenabled;//Get the value if "highlight reconstruction" is activated
    bool hasToneCurve1 = bool (customToneCurve1);
    bool hasToneCurve2 = bool (customToneCurve2);


#define TS 112

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        char *buffer;

        buffer = (char *) malloc (3 * sizeof (float) * TS * TS + 20 * 64 + 63);
        char *data;
        data = (char*) ( ( uintptr_t (buffer) + uintptr_t (63)) / 64 * 64);

        float *Ltemp = (float (*))data;
        float *atemp = (float (*))         ((char*)Ltemp + sizeof (float) * TS * TS + 4 * 64);
        float *btemp = (float (*))         ((char*)atemp + sizeof (float) * TS * TS + 8 * 64);
        int istart;
        int jstart;
        int tW;
        int tH;

#ifdef _OPENMP
        #pragma omp for schedule(dynamic) collapse(2)
#endif

        for (int ii = 0; ii < bfh; ii += TS)
            for (int jj = 0; jj < bfw; jj += TS) {

                istart = ii;
                jstart = jj;
                tH = min (ii + TS, bfh);
                tW = min (jj + TS, bfw);


                for (int i = istart, ti = 0; i < tH; i++, ti++) {
                    for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                        Ltemp[ti * TS + tj] = 2.f * bufexporig->L[i][j];
                        atemp[ti * TS + tj] = bufexporig->a[i][j];
                        btemp[ti * TS + tj] = bufexporig->b[i][j];;
                    }
                }



                for (int i = istart, ti = 0; i < tH; i++, ti++) {
                    for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                        float L = Ltemp[ti * TS + tj];
                        float a = atemp[ti * TS + tj];
                        float b = btemp[ti * TS + tj];
                        float tonefactor = (L < MAXVALF ? hltonecurve[L] : CurveFactory::hlcurve (exp_scale, comp, hlrange, L) );
                        Ltemp[ti * TS + tj] = L * tonefactor;
                    }
                }

                for (int i = istart, ti = 0; i < tH; i++, ti++) {
                    for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                        float L = Ltemp[ti * TS + tj];
                        float a = atemp[ti * TS + tj];
                        float b = btemp[ti * TS + tj];

                        //shadow tone curve
                        float Y = L;
                        float tonefactor = shtonecurve[Y];
                        Ltemp[ti * TS + tj] = Ltemp[ti * TS + tj] * tonefactor;
                    }
                }

                if (dcpProf) {
                    //          dcpProf->step2ApplyTile (rtemp, gtemp, btemp, tW - jstart, tH - istart, TS, asIn);
                }

                /*
                                for (int i = istart, ti = 0; i < tH; i++, ti++) {
                                    for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                        float r = rtemp[ti * TS + tj];
                                        float g = gtemp[ti * TS + tj];
                                        float b = btemp[ti * TS + tj];

                                        // clip out of gamut colors, without distorting color too bad
                                        if (r < 0) {
                                            r = 0;
                                        }

                                        if (g < 0) {
                                            g = 0;
                                        }

                                        if (b < 0) {
                                            b = 0;
                                        }

                                        if (r > 65535 || g > 65535 || b > 65535) {
                                            filmlike_clip (&r, &g, &b);
                                        }

                                        rtemp[ti * TS + tj] = r;
                                        gtemp[ti * TS + tj] = g;
                                        btemp[ti * TS + tj] = b;
                                    }
                                }
                */
                for (int i = istart, ti = 0; i < tH; i++, ti++) {
                    for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                        //brightness/contrast
                        Ltemp[ti * TS + tj] = tonecurve[Ltemp[ti * TS + tj] ];

                    }
                }


                if (hasToneCurve1) {
                    if (curveMode == LocalrgbParams::TC_MODE_STD) { // Standard
                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                const StandardToneCurveL& userToneCurve = static_cast<const StandardToneCurveL&> (customToneCurve1);
                                userToneCurve.Apply (Ltemp[ti * TS + tj]);
                            }
                        }
                    }
                }


                if (hasToneCurve2) {
                    if (curveMode2 == LocalrgbParams::TC_MODE_STD) { // Standard
                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                const StandardToneCurveL& userToneCurve = static_cast<const StandardToneCurveL&> (customToneCurve2);
                                userToneCurve.Apply (Ltemp[ti * TS + tj]);
                            }
                        }
                    }
                }


                if (lp.chro != 0) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                            float satby100 = lp.chro / 100.f;
                            float L = 2.f * Ltemp[ti * TS + tj];
                            float a = atemp[ti * TS + tj];
                            float b = btemp[ti * TS + tj];

                            atemp[ti * TS + tj] = a * (1.f + satby100);
                            btemp[ti * TS + tj] = b * (1.f + satby100);
                        }
                    }
                }

                /*
                                if (isProPhoto) { // this is a hack to avoid the blue=>black bug (Issue 2141)
                                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                            float r = rtemp[ti * TS + tj];
                                            float g = gtemp[ti * TS + tj];

                                            if (r == 0.0f || g == 0.0f) {
                                                float b = btemp[ti * TS + tj];
                                                float h, s, v;
                                                Color::rgb2hsv (r, g, b, h, s, v);
                                                s *= 0.99f;
                                                Color::hsv2rgb (h, s, v, rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj]);
                                            }
                                        }
                                    }
                                }

                */

                bool vasy = true;

                if (vasy) {
                    // ready, fill lab
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                            lab->L[i][j] = 0.5f * Ltemp[ti * TS + tj];
                            lab->a[i][j] = atemp[ti * TS + tj];
                            lab->b[i][j] = btemp[ti * TS + tj];

                            //in case off futur usage
                            /*
                            origrgb->r (i, j) = r;
                            origrgb->g (i, j) = g;
                            origrgb->b (i, j) = b;
                            */
                        }
                    }
                }
            }

        free (buffer);


    }


}

void ImProcFunctions::Rgb_Local (int call, int sp, LabImage* original, LabImage* transformed, int sx, int sy, int cx, int cy, int oW, int oH,  int fw, int fh, double &hueref, double &chromaref, double &lumaref,
                                 Imagefloat* working, LabImage* lab, Imagefloat* orirgb, LUTf & hltonecurveloc, LUTf & shtonecurveloc, LUTf & tonecurveloc,
                                 int sat, const ToneCurve & customToneCurve1, const ToneCurve & customToneCurve2,
                                 double expcomp, int hlcompr, int hlcomprthresh, DCPProfile *dcpProf, const DCPProfile::ApplyState &as)
{
    if (params->localrgb.enabled) {
        // BENCHFUN
#ifdef _DEBUG
        MyTime t1e, t2e;
        t1e.set();
// init variables to display Munsell corrections
        MunsellDebugInfo* MunsDebugInfo = new MunsellDebugInfo();
#endif

        int del = 3; // to avoid crash with [loy - begy] and [lox - begx] and bfh bfw  // with gtk2 [loy - begy-1] [lox - begx -1 ] and del = 1
        float moy = 0.f;

        struct local_params lp;
        calcLocalParams (oW, oH, params->localrgb, lp);



        double ave = 0.;
        int n = 0;
        float av = 0;
//        int levred;
//        bool noiscfactiv = false;

        // no activ actually because in Locallab ==> todo
        /*
                if (lp.qualmet == 2) { //suppress artifacts with quality enhanced
                    levred = 4;
                    noiscfactiv = true;
                }    else {
                    levred = 7;
                    noiscfactiv = false;
                }
        */
// we must here detect : general case, skin, sky,...foliages ???
// delta dhue, luminance and chroma
        constexpr float ared = (rtengine::RT_PI - 0.05f) / 100.f;

        constexpr float bred = 0.05f;

        float dhue = ared * lp.sens + bred; //delta hue lght chroma
        float dhuev = ared * lp.sensv + bred; //delta hue lght chroma


        constexpr float maxh = 3.5f; // 3.5 amplification contrast above mean

        constexpr float maxl = 2.5f; // 3 reductio contrast under mean

        float multh = (float) fabs (lp.cont) * (maxh - 1.f) / 100.f + 1.f;

        float mult = (float)fabs (lp.cont) * (maxl - 1.f) / 100.f + 1.f;


        if ((lp.chro != 0 || lp.ligh != 0.f || lp.cont != 0.f || lp.expcomp != 0.f  || customToneCurve1 || customToneCurve2) && lp.exposeena) { //interior ellipse renforced lightness and chroma  //locallutili
            //  if (lp.expvib) { //interior ellipse renforced lightness and chroma  //locallutili
            printf ("OK appel vib loc\n");
            float hueplus = hueref + dhue;
            float huemoins = hueref - dhue;

            printf ("hueplus=%f huemoins=%f dhu=%f\n", hueplus, huemoins, dhue);

            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhue - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhue + 2.f * rtengine::RT_PI;
            }

            LabImage *bufexporig = nullptr;
            LabImage *bufexpfin = nullptr;
            Imagefloat* bufworking = nullptr;
            float **buflight = nullptr;
            float **bufl_ab = nullptr;

            int bfh = 0.f, bfw = 0.f;


            if (call <= 3) { //simpleprocess, dcrop, improccoordinator
                bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
                bfw = int (lp.lx + lp.lxL) + del;


                bufexporig = new LabImage (bfw, bfh);//buffer for data in zone limit
                bufexpfin = new LabImage (bfw, bfh);//buffer for data in zone limit
                bufworking = new Imagefloat (bfw, bfh);

                buflight   = new float*[bfh];//for lightness

                for (int i = 0; i < bfh; i++) {
                    buflight[i] = new float[bfw];
                }

                bufl_ab   = new float*[bfh];//for chroma

                for (int i = 0; i < bfh; i++) {
                    bufl_ab[i] = new float[bfw];
                }

#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int ir = 0; ir < bfh; ir++) //fill with 0
                    for (int jr = 0; jr < bfw; jr++) {
                        bufexporig->L[ir][jr] = 0.f;
                        bufexporig->a[ir][jr] = 0.f;
                        bufexporig->b[ir][jr] = 0.f;
                        bufworking->r (ir, jr) = 0.f;
                        bufworking->g (ir, jr) = 0.f;
                        bufworking->b (ir, jr) = 0.f;
                        bufexpfin->L[ir][jr] = 0.f;
                        bufexpfin->a[ir][jr] = 0.f;
                        bufexpfin->b[ir][jr] = 0.f;
                        buflight[ir][jr] = 0.f;
                        bufl_ab[ir][jr] = 0.f;


                    }

                int begy = lp.yc - lp.lyT;
                int begx = lp.xc - lp.lxL;
                int yEn = lp.yc + lp.ly;
                int xEn = lp.xc + lp.lx;
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                            bufworking->r (loy - begy, lox - begx) = working->r (y, x); //fill square buffer with datas
                            bufworking->g (loy - begy, lox - begx) = working->g (y, x); //fill square buffer with datas
                            bufworking->b (loy - begy, lox - begx) = working->b (y, x); //fill square buffer with datas

                            bufexporig->L[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                            bufexporig->a[loy - begy][lox - begx] = original->a[y][x];//fill square buffer with datas
                            bufexporig->b[loy - begy][lox - begx] = original->b[y][x];//fill square buffer with datas

                        }
                    }



                //     ImProcFunctions::rgbLocal (bufworking, bufexpfin, orirgb, hltonecurveloc, shtonecurveloc, tonecurveloc, lp.chro,
                //                               customToneCurve1, customToneCurve2, lp.expcomp, lp.hlcomp, lp.hlcompthr, dcpProf, as);

                ImProcFunctions::rgblabLocal (lp, bufworking, bfh, bfw, bufexporig, bufexpfin, orirgb, hltonecurveloc, shtonecurveloc, tonecurveloc,
                                              customToneCurve1, customToneCurve2, dcpProf, as);
                //      ImProcFunctions::vibrancelocal (lp, bfw, bfh,bufexporig, bufexpfin);


                //          float maxc = -10000.f;
                //          float minc = 100000.f;
                //          float chpro = 0.f;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {

                            float lL;
                            float amplil = 140.f;
                            float lighL = bufexporig->L[loy - begy][lox - begx];
                            float lighLnew = bufexpfin->L[loy - begy][lox - begx];
                            float rL;
                            rL = CLIPRET ((bufexpfin->L[loy - begy][lox - begx] - bufexporig->L[loy - begy][lox - begx]) / 328.f);

                            buflight[loy - begy][lox - begx] = rL;


                            float chp;
                            chp = CLIPRET ((sqrt (SQR (bufexpfin->a[loy - begy][lox - begx]) + SQR (bufexpfin->b[loy - begy][lox - begx])) - sqrt (SQR (bufexporig->a[loy - begy][lox - begx]) + SQR (bufexporig->b[loy - begy][lox - begx]))) / 250.f);
                            /*
                              if (chp > maxc) {
                                   maxc = chp;
                               }

                               if (chp < minc) {
                                   minc = chp;
                               }
                            */
                            //  chpro = CLIPCHRO (amplil * ra - amplil); //ampli = 25.f arbitrary empirical coefficient between 5 and 50

                            //ra = 1.f;
                            bufl_ab[loy - begy][lox - begx] = chp;

                        }
                    }

                //   printf ("min=%2.2f max=%2.2f", minc, maxc);
                Expose_Local (1, call, buflight, bufl_ab, hueplus, huemoins, hueref, dhue, chromaref, lumaref, lp, original, transformed, bufexpfin, cx, cy);

            }

            if (call <= 3) {

                delete bufworking;
                delete bufexporig;
                delete bufexpfin;

                for (int i = 0; i < bfh; i++) {
                    delete [] buflight[i];
                }

                delete [] buflight;

                for (int i = 0; i < bfh; i++) {
                    delete [] bufl_ab[i];
                }

                delete [] bufl_ab;

            }

        }

        // vibrance
        if (lp.expvib && (lp.past != 0.f  || lp.satur != 0.f)) { //interior ellipse renforced lightness and chroma  //locallutili
            //  printf("OK appel vib loc\n");
            float hueplus = hueref + dhuev;
            float huemoins = hueref - dhuev;

          //  printf ("hueplus=%f huemoins=%f dhu=%f\n", hueplus, huemoins, dhuev);

            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhuev - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhuev + 2.f * rtengine::RT_PI;
            }

            LabImage *bufexporig = nullptr;
            LabImage *bufexpfin = nullptr;
            Imagefloat* bufworking = nullptr;
            float **buflight = nullptr;
            float **bufl_ab = nullptr;

            int bfh = 0.f, bfw = 0.f;


            if (call <= 3) { //simpleprocess, dcrop, improccoordinator
                bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
                bfw = int (lp.lx + lp.lxL) + del;


                bufexporig = new LabImage (bfw, bfh);//buffer for data in zone limit
                bufexpfin = new LabImage (bfw, bfh);//buffer for data in zone limit
                bufworking = new Imagefloat (bfw, bfh);

                buflight   = new float*[bfh];//for lightness

                for (int i = 0; i < bfh; i++) {
                    buflight[i] = new float[bfw];
                }

                bufl_ab   = new float*[bfh];//for chroma

                for (int i = 0; i < bfh; i++) {
                    bufl_ab[i] = new float[bfw];
                }

#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int ir = 0; ir < bfh; ir++) //fill with 0
                    for (int jr = 0; jr < bfw; jr++) {
                        bufexporig->L[ir][jr] = 0.f;
                        bufexporig->a[ir][jr] = 0.f;
                        bufexporig->b[ir][jr] = 0.f;
                        bufworking->r (ir, jr) = 0.f;
                        bufworking->g (ir, jr) = 0.f;
                        bufworking->b (ir, jr) = 0.f;
                        bufexpfin->L[ir][jr] = 0.f;
                        bufexpfin->a[ir][jr] = 0.f;
                        bufexpfin->b[ir][jr] = 0.f;
                        buflight[ir][jr] = 0.f;
                        bufl_ab[ir][jr] = 0.f;


                    }

                int begy = lp.yc - lp.lyT;
                int begx = lp.xc - lp.lxL;
                int yEn = lp.yc + lp.ly;
                int xEn = lp.xc + lp.lx;
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                            bufworking->r (loy - begy, lox - begx) = working->r (y, x); //fill square buffer with datas
                            bufworking->g (loy - begy, lox - begx) = working->g (y, x); //fill square buffer with datas
                            bufworking->b (loy - begy, lox - begx) = working->b (y, x); //fill square buffer with datas

                            bufexporig->L[loy - begy][lox - begx] = original->L[y][x];//fill square buffer with datas
                            bufexporig->a[loy - begy][lox - begx] = original->a[y][x];//fill square buffer with datas
                            bufexporig->b[loy - begy][lox - begx] = original->b[y][x];//fill square buffer with datas

                        }
                    }



                ImProcFunctions::vibrancelocal (lp, bfw, bfh, bufexporig, bufexpfin);


                //          float maxc = -10000.f;
                //          float minc = 100000.f;
                //          float chpro = 0.f;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {

                            float lL;
                            float amplil = 140.f;
                            float lighL = bufexporig->L[loy - begy][lox - begx];
                            float lighLnew = bufexpfin->L[loy - begy][lox - begx];
                            float rL;
                            rL = CLIPRET ((bufexpfin->L[loy - begy][lox - begx] - bufexporig->L[loy - begy][lox - begx]) / 328.f);

                            buflight[loy - begy][lox - begx] = rL;


                            float chp;
                            chp = CLIPRET ((sqrt (SQR (bufexpfin->a[loy - begy][lox - begx]) + SQR (bufexpfin->b[loy - begy][lox - begx])) - sqrt (SQR (bufexporig->a[loy - begy][lox - begx]) + SQR (bufexporig->b[loy - begy][lox - begx]))) / 250.f);
                            /*
                              if (chp > maxc) {
                                   maxc = chp;
                               }

                               if (chp < minc) {
                                   minc = chp;
                               }
                            */
                            //  chpro = CLIPCHRO (amplil * ra - amplil); //ampli = 25.f arbitrary empirical coefficient between 5 and 50

                            //ra = 1.f;
                            bufl_ab[loy - begy][lox - begx] = chp;

                        }
                    }

                //   printf ("min=%2.2f max=%2.2f", minc, maxc);
                Expose_Local (2, call, buflight, bufl_ab, hueplus, huemoins, hueref, dhuev, chromaref, lumaref, lp, original, transformed, bufexpfin, cx, cy);

            }

            if (call <= 3) {

                delete bufworking;
                delete bufexporig;
                delete bufexpfin;

                for (int i = 0; i < bfh; i++) {
                    delete [] buflight[i];
                }

                delete [] buflight;

                for (int i = 0; i < bfh; i++) {
                    delete [] bufl_ab[i];
                }

                delete [] bufl_ab;

            }

        }



    }


}

void ImProcFunctions::Whitebalance_Local (int call, int sp, Imagefloat* bufimage, const struct local_params & lp, Imagefloat* imageoriginal, Imagefloat* imagetransformed, int cx, int cy)
{
    BENCHFUN

    const float ach = (float)lp.trans / 100.f;


    /*
    constexpr float delhu = 0.1f; //between 0.05 and 0.2

    constexpr float apl = (-1.f) / delhu;
    const float bpl = - apl * hueplus;
    constexpr float amo = 1.f / delhu;
    const float bmo = - amo * huemoins;


    constexpr float pb = 4.f;
    constexpr float pa = (1.f - pb) / 40.f;

    const float ahu = 1.f / (2.8f * lp.sensbn - 280.f);
    const float bhu = 1.f - ahu * 2.8f * lp.sensbn;
    */
#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef __SSE2__
        //  float atan2Buffer[transformed->W] ALIGNED16;
        //  float sqrtBuffer[transformed->W] ALIGNED16;
        //  vfloat c327d68v = F2V (327.68f);
#endif

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int y = 0; y < imagetransformed->getHeight(); y++) {
            const int loy = cy + y;

            const bool isZone0 = loy > lp.yc + lp.ly || loy < lp.yc - lp.lyT; // whole line is zone 0 => we can skip a lot of processing

            if (isZone0) { // outside selection and outside transition zone => no effect, keep original values
                for (int x = 0; x < imagetransformed->getWidth(); x++) {
                    imagetransformed->r (y, x) = imageoriginal->r (y, x);
                    imagetransformed->g (y, x) = imageoriginal->g (y, x);
                    imagetransformed->b (y, x) = imageoriginal->b (y, x);

                }


                continue;
            }

#ifdef __SSE2__
            int i = 0;

            for (; i < imagetransformed->getWidth() - 3; i += 4) {
                //    vfloat av = LVFU (original->a[y][i]);
                //    vfloat bv = LVFU (original->b[y][i]);
                //    STVF (atan2Buffer[i], xatan2f (bv, av));
                //    STVF (sqrtBuffer[i], _mm_sqrt_ps (SQRV (bv) + SQRV (av)) / c327d68v);
            }

            for (; i < imagetransformed->getWidth(); i++) {
                //    atan2Buffer[i] = xatan2f (original->b[y][i], original->a[y][i]);
                //    sqrtBuffer[i] = sqrt (SQR (original->b[y][i]) + SQR (original->a[y][i])) / 327.68f;
            }

#endif



            for (int x = 0, lox = cx + x; x < imagetransformed->getWidth(); x++, lox++) {
                int zone = 0;
                float localFactor = 1.f;
                calcTransition (lox, loy, ach, lp, zone, localFactor);

                if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                    imagetransformed->r (y, x) = imageoriginal->r (y, x);
                    imagetransformed->g (y, x) = imageoriginal->g (y, x);
                    imagetransformed->b (y, x) = imageoriginal->b (y, x);

                    continue;
                }

#ifdef __SSE2__
                //    const float rhue = atan2Buffer[x];
                //    const float rchro = sqrtBuffer[x];
#else
                //    const float rhue = xatan2f (original->b[y][x], original->a[y][x]);
                //    const float rchro = sqrt (SQR (original->b[y][x]) + SQR (original->a[y][x])) / 327.68f;
#endif

                //prepare shape detection
                float kch = 1.f;
                float fach = 1.f;
                //  float deltachro = fabs (rchro - chromaref);
                //  float deltahue = fabs (rhue - hueref);
                /*
                                if (deltahue > rtengine::RT_PI) {
                                    deltahue = - (deltahue - 2.f * rtengine::RT_PI);
                                }

                                float deltaE = 20.f * deltahue + deltachro; //pseudo deltaE between 0 and 280

                                //kch to modulate action with chroma
                                if (deltachro < 160.f * SQR (lp.sensbn / 100.f)) {
                                    kch = 1.f;
                                } else {
                                    float ck = 160.f * SQR (lp.sensbn / 100.f);
                                    float ak = 1.f / (ck - 160.f);
                                    float bk = -160.f * ak;
                                    kch = ak * deltachro + bk;
                                }

                                if (lp.sensbn < 40.f ) {
                                    float khu = 0.f;
                                    kch = pow (kch, pa * lp.sensbn + pb);   //increase under 40


                                    // algo with detection of hue ==> artifacts for noisy images  ==> denoise before
                                    if (lp.qualmet == 1 && lp.sensbn < 20.f) { //to try...
                                        //hue detection
                                        if ((hueref + dhue) < rtengine::RT_PI && rhue < hueplus && rhue > huemoins) { //transition are good
                                            if (rhue >= hueplus - delhu )  {
                                                khu  = apl * rhue + bpl;
                                            } else if (rhue < huemoins + delhu)  {
                                                khu = amo * rhue + bmo;
                                            } else {
                                                khu = 1.f;
                                            }
                                        } else if ((hueref + dhue) >= rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                                            if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                                                khu  = apl * rhue + bpl;
                                            } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                                                khu = amo * rhue + bmo;
                                            } else {
                                                khu = 1.f;
                                            }
                                        }

                                        if ((hueref - dhue) > -rtengine::RT_PI && rhue < hueplus && rhue > huemoins ) {
                                            if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                                                khu  = apl * rhue + bpl;
                                            } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                                                khu = amo * rhue + bmo;
                                            } else {
                                                khu = 1.f;
                                            }
                                        } else if ((hueref - dhue) <= -rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                                            if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                                                khu  = apl * rhue + bpl;
                                            } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                                                khu = amo * rhue + bmo;
                                            } else {
                                                khu = 1.f;
                                            }
                                        }

                                        if (deltaE <  2.8f * lp.sensbn) {
                                            fach = khu;
                                        } else {
                                            fach = khu * (ahu * deltaE + bhu);
                                        }


                                        float kcr = 10.f;

                                        if (rchro < kcr) {
                                            fach *= (1.f / (kcr * kcr)) * rchro * rchro;
                                        }
                                    }
                                }
                */
                int begx = lp.xc - lp.lxL;
                int begy = lp.yc - lp.lyT;

                switch (zone) {

                    case 1: { // inside transition zone
                        float difr = 0.f, difg = 0.f, difb = 0.f;

                        if (call <= 3) {

                            difr = bufimage->r (loy - begy, lox - begx) - imageoriginal->r (y, x);
                            difg = bufimage->g (loy - begy, lox - begx) - imageoriginal->g (y, x);
                            difb = bufimage->b (loy - begy, lox - begx) - imageoriginal->b (y, x);
                            /*
                                                        difL = tmp1->L[loy - begy][lox - begx] - original->L[y][x];
                                                        difa = tmp1->a[loy - begy][lox - begx] - original->a[y][x];
                                                        difb = tmp1->b[loy - begy][lox - begx] - original->b[y][x];
                            */
                        }

                        difr *= localFactor;
                        difg *= localFactor;
                        difb *= localFactor;

                        imagetransformed->r (y, x) = imageoriginal->r (y, x) + difr;
                        imagetransformed->g (y, x) = imageoriginal->g (y, x) + difg;
                        imagetransformed->b (y, x) = imageoriginal->b (y, x) + difb;

                        // transformed->L[y][x] = original->L[y][x] + difL * kch * fach;


                        break;
                    }

                    case 2: { // inside selection => full effect, no transition
                        float difr = 0.f, difg = 0.f, difb = 0.f;

                        if (call <= 3) {
                            difr = bufimage->r (loy - begy, lox - begx) - imageoriginal->r (y, x);
                            difg = bufimage->g (loy - begy, lox - begx) - imageoriginal->g (y, x);
                            difb = bufimage->b (loy - begy, lox - begx) - imageoriginal->b (y, x);

                        }

                        imagetransformed->r (y, x) = imageoriginal->r (y, x) + difr;
                        imagetransformed->g (y, x) = imageoriginal->g (y, x) + difg;
                        imagetransformed->b (y, x) = imageoriginal->b (y, x) + difb;

                        //    transformed->L[y][x] = original->L[y][x] + difL * kch * fach;

                    }
                }
            }
        }
    }
}




void ImProcFunctions::Expose_Local (int sen, int call, float **buflight, float **bufchro, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, const struct local_params & lp, LabImage * original, LabImage * transformed, const LabImage * const tmp1, int cx, int cy)
{

//local exposure
    BENCHFUN {
        const float ach = (float)lp.trans / 100.f;
        float varsens =  lp.sens;

        if (sen == 1)
        {
            varsens =  lp.sens;
        }

        if (sen == 2)
        {
            varsens =  lp.sensv;
        }

        //chroma
        constexpr float amplchsens = 2.5f;
        constexpr float achsens = (amplchsens - 1.f) / (100.f - 20.f); //20. default locallab.sensih
        constexpr float bchsens = 1.f - 20.f * achsens;
        const float multchro = varsens * achsens + bchsens;

        //luma

        //skin
        constexpr float amplchsensskin = 1.6f;
        constexpr float achsensskin = (amplchsensskin - 1.f) / (100.f - 20.f); //20. default locallab.sensih
        constexpr float bchsensskin = 1.f - 20.f * achsensskin;
        const float multchroskin = varsens * achsensskin + bchsensskin;

        //transition = difficult to avoid artifact with scope on flat area (sky...)

        constexpr float delhu = 0.1f; //between 0.05 and 0.2

        const float apl = (-1.f) / delhu;
        const float bpl = - apl * hueplus;
        const float amo = 1.f / delhu;
        const float bmo = - amo * huemoins;


        const float pb = 4.f;
        const float pa = (1.f - pb) / 40.f;

        const float ahu = 1.f / (2.8f * varsens - 280.f);
        const float bhu = 1.f - ahu * 2.8f * varsens;

        const float alum = 1.f / (varsens - 100.f);
        const float blum = 1.f - alum * varsens;
        //float maxc = -100000.f;
        //float minc = 100000.f;

#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
#ifdef __SSE2__
            float atan2Buffer[transformed->W] ALIGNED16;
            float sqrtBuffer[transformed->W] ALIGNED16;
            vfloat c327d68v = F2V (327.68f);
#endif

#ifdef _OPENMP
            #pragma omp for schedule(dynamic,16)
#endif

            for (int y = 0; y < transformed->H; y++)
            {

                const int loy = cy + y;
                const bool isZone0 = loy > lp.yc + lp.ly || loy < lp.yc - lp.lyT; // whole line is zone 0 => we can skip a lot of processing

                if (isZone0) { // outside selection and outside transition zone => no effect, keep original values
                    for (int x = 0; x < transformed->W; x++) {
                        transformed->L[y][x] = original->L[y][x];
                    }

                    continue;
                }

#ifdef __SSE2__
                int i = 0;

                for (; i < transformed->W - 3; i += 4) {
                    vfloat av = LVFU (original->a[y][i]);
                    vfloat bv = LVFU (original->b[y][i]);
                    STVF (atan2Buffer[i], xatan2f (bv, av));
                    STVF (sqrtBuffer[i], _mm_sqrt_ps (SQRV (bv) + SQRV (av)) / c327d68v);
                }

                for (; i < transformed->W; i++) {
                    atan2Buffer[i] = xatan2f (original->b[y][i], original->a[y][i]);
                    sqrtBuffer[i] = sqrt (SQR (original->b[y][i]) + SQR (original->a[y][i])) / 327.68f;
                }

#endif

                for (int x = 0; x < transformed->W; x++) {
                    int lox = cx + x;
                    int begx = int (lp.xc - lp.lxL);
                    int begy = int (lp.yc - lp.lyT);

                    int zone = 0;
                    float localFactor = 1.f;
                    calcTransition (lox, loy, ach, lp, zone, localFactor);

                    if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                        transformed->L[y][x] = original->L[y][x];
                        continue;
                    }

#ifdef __SSE2__
                    float rhue = atan2Buffer[x];
                    float rchro = sqrtBuffer[x];
#else
                    float rhue = xatan2f (original->b[y][x], original->a[y][x]);
                    float rchro = sqrt (SQR (original->b[y][x]) + SQR (original->a[y][x])) / 327.68f;
#endif
                    float rL = original->L[y][x] / 327.68f;

                    float cli = 1.f;
                    float clc = 1.f;

                    //                       if (lp.curvact == true) {
                    cli = (buflight[loy - begy][lox - begx]);
                    clc = (bufchro[loy - begy][lox - begx]);
                    //printf("cl=%f",clc);

                    //                      } else {
                    //                          cli = lp.str;
                    //                           clc = params->locallab.chrrt;
                    //                       }

                    float aplus = (1.f - cli) / delhu;
                    float bplus = 1.f - aplus * hueplus;
                    float amoins = (cli - 1.f) / delhu;
                    float bmoins = 1.f - amoins * huemoins;

                    float aplusch = (1.f - clc) / delhu;
                    float bplusch = 1.f - aplusch * hueplus;
                    float amoinsch = (clc - 1.f) / delhu;
                    float bmoinsch = 1.f - amoinsch * huemoins;

                    float realstr = 1.f;
                    float realstrch = 1.f;
                    //prepare shape detection
                    float deltachro = fabs (rchro - chromaref);
                    float deltahue = fabs (rhue - hueref);

                    if (deltahue > rtengine::RT_PI) {
                        deltahue = - (deltahue - 2.f * rtengine::RT_PI);
                    }

                    float deltaE = 20.f * deltahue + deltachro; //between 0 and 280
                    float deltaL = fabs (lumaref - rL); //between 0 and 100

                    float kch = 1.f;
                    float khu = 0.f;
                    float fach = 1.f;
                    float falu = 1.f;

                    if (deltachro < 160.f * SQR (varsens / 100.f)) {
                        kch = 1.f;
                    } else {
                        float ck = 160.f * SQR (varsens / 100.f);
                        float ak = 1.f / (ck - 160.f);
                        float bk = -160.f * ak;
                        kch = ak * deltachro + bk;
                    }

                    if (varsens < 40.f ) {
                        kch = pow (kch, pa * varsens + pb);   //increase under 40
                    }

                    bool kzon = false;

                    //transition = difficult to avoid artifact with scope on flat area (sky...)
                    //hue detection
                    if ((hueref + dhue) < rtengine::RT_PI && rhue < hueplus && rhue > huemoins) { //transition are good
                        if (rhue >= hueplus - delhu)  {
                            realstr = aplus * rhue + bplus;
                            realstrch = aplusch * rhue + bplusch;
                            khu  = apl * rhue + bpl;

                        } else if (rhue < huemoins + delhu)  {
                            realstr = amoins * rhue + bmoins;
                            realstrch = amoinsch * rhue + bmoinsch;
                            khu = amo * rhue + bmo;

                        } else {
                            realstr = cli;
                            khu = 1.f;
                            realstrch = clc;

                        }

                        kzon = true;
                    } else if ((hueref + dhue) >= rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            realstr = aplus * rhue + bplus;
                            realstrch = aplusch * rhue + bplusch;
                            khu  = apl * rhue + bpl;

                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            realstr = amoins * rhue + bmoins;
                            realstrch = amoinsch * rhue + bmoinsch;
                            khu = amo * rhue + bmo;

                        } else {
                            realstr = cli;
                            khu = 1.f;
                            realstrch = clc;

                        }

                        kzon = true;
                    }

                    if ((hueref - dhue) > -rtengine::RT_PI && rhue < hueplus && rhue > huemoins) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            realstr = aplus * rhue + bplus;
                            realstrch = aplusch * rhue + bplusch;
                            khu  = apl * rhue + bpl;

                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            realstr = amoins * rhue + bmoins;
                            realstrch = amoinsch * rhue + bmoinsch;
                            khu = amo * rhue + bmo;

                        } else {
                            realstr = cli;
                            khu = 1.f;
                            realstrch = clc;

                        }

                        kzon = true;
                    } else if ((hueref - dhue) <= -rtengine::RT_PI && (rhue > huemoins  || rhue < hueplus )) {
                        if (rhue >= hueplus - delhu  && rhue < hueplus)  {
                            realstr = aplus * rhue + bplus;
                            realstrch = aplusch * rhue + bplusch;
                            khu  = apl * rhue + bpl;

                        } else if (rhue >= huemoins && rhue < huemoins + delhu)  {
                            realstr = amoins * rhue + bmoins;
                            realstrch = amoinsch * rhue + bmoinsch;
                            khu = amo * rhue + bmo;

                        } else {
                            realstr = cli;
                            khu = 1.f;
                            realstrch = clc;

                        }

                        kzon = true;
                    }

                    /*
                                        //printf("re=%f",  realstrch);
                                        if (realstrch > maxc) {
                                            maxc = realstrch;
                                        }

                                        if (realstrch < minc) {
                                            minc = realstrch;
                                        }
                    */
                    //shape detection for hue chroma and luma
                    if (varsens <= 20.f) { //to try...

                        if (deltaE <  2.8f * varsens) {
                            fach = khu;
                        } else {
                            fach = khu * (ahu * deltaE + bhu);
                        }

                        float kcr = 10.f;

                        if (rchro < kcr) {
                            fach *= (1.f / (kcr * kcr)) * rchro * rchro;
                        }

                        if (lp.qualmet >= 1) {
                        } else {
                            fach = 1.f;
                        }

                        if (deltaL <  varsens) {
                            falu = 1.f;
                        } else {
                            falu = alum * deltaL + blum;
                        }

                    }

                    //    float kdiff = 0.f;
                    // I add these functions...perhaps not good
                    if (kzon) {
                        if (varsens < 60.f) { //arbitrary value
                            if (hueref < -1.1f && hueref > -2.8f) { // detect blue sky
                                if (chromaref > 0.f && chromaref < 35.f * multchro) { // detect blue sky
                                    if ( (rhue > -2.79f && rhue < -1.11f) && (rchro < 35.f * multchro)) {
                                        realstr *= 0.9f;
                                    } else {
                                        realstr = 1.f;
                                    }
                                }
                            } else {
                                realstr = cli;
                            }

                            if (varsens < 50.f) { //&& lp.chro > 0.f
                                if (hueref > -0.1f && hueref < 1.6f) { // detect skin
                                    if (chromaref > 0.f && chromaref < 55.f * multchroskin) { // detect skin
                                        if ( (rhue > -0.09f && rhue < 1.59f) && (rchro < 55.f * multchroskin)) {
                                            realstr *= 0.7f;
                                        } else {
                                            realstr = 1.f;
                                        }
                                    }
                                } else {
                                    realstr = cli;
                                }
                            }
                        }

                    }

                    float kcr = 100.f * lp.thr;
                    float falL = 1.f;

                    if (rchro < kcr && chromaref > kcr) { // reduce artifacts in grey tones near hue spot and improve algorithm
                        falL *= pow (rchro / kcr, lp.iterat / 10.f);
                    }

                    //                     int zone;
                    //                     float localFactor;
                    //                     calcTransition (lox, loy, ach, lp, zone, localFactor);

                    if (rL > 0.1f) { //to avoid crash with very low gamut in rare cases ex : L=0.01 a=0.5 b=-0.9
                        switch (zone) {
                            case 0: { // outside selection and outside transition zone => no effect, keep original values
                                transformed->L[y][x] = original->L[y][x];
                                transformed->a[y][x] = original->a[y][x];
                                transformed->b[y][x] = original->b[y][x];

                                break;
                            }

                            case 1: { // inside transition zone
                                float factorx = localFactor;

                                float difL;
                                difL = tmp1->L[loy - begy][lox - begx] - original->L[y][x];
                                difL *= factorx * (100.f + realstr * falL) / 100.f;
                                difL *= kch * fach;

                                transformed->L[y][x] = original->L[y][x] + difL;
                                float difa, difb;

                                difa = tmp1->a[loy - begy][lox - begx] - original->a[y][x];
                                difb = tmp1->b[loy - begy][lox - begx] - original->b[y][x];
                                difa *= factorx * (100.f + realstrch * falu * falL) / 100.f;
                                difb *= factorx * (100.f + realstrch * falu * falL) / 100.f;
                                difa *= kch * fach;
                                difb *= kch * fach;
                                transformed->a[y][x] = CLIPC (original->a[y][x] + difa);
                                transformed->b[y][x] = CLIPC (original->b[y][x] + difb);


                                break;

                            }

                            case 2: { // inside selection => full effect, no transition
                                float difL;

                                difL = tmp1->L[loy - begy][lox - begx] - original->L[y][x];
                                difL *= (100.f + realstr * falL) / 100.f;
                                difL *= kch * fach;
                                transformed->L[y][x] = original->L[y][x] + difL;
                                float difa, difb;

                                difa = tmp1->a[loy - begy][lox - begx] - original->a[y][x];
                                difb = tmp1->b[loy - begy][lox - begx] - original->b[y][x];
                                difa *= (100.f + realstrch * falu * falL) / 100.f;
                                difb *= (100.f + realstrch * falu * falL) / 100.f;
                                difa *= kch * fach;
                                difb *= kch * fach;

                                transformed->a[y][x] = CLIPC (original->a[y][x] + difa);
                                transformed->b[y][x] = CLIPC (original->b[y][x] + difb);

                            }
                        }

                        //}
                    }
                }
            }

            //    printf ("minc=%f maxc=%f \n", minc, maxc);
        }

    }
}


}