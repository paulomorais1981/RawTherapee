/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
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

namespace rtengine
{
using namespace procparams;

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
    int chro, cont, sens;
    float ligh;
    int qualmet;
    float threshol;
    bool exposeena;
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


};

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
    lp.expwb = localrgb.expwb;
    lp.blac = localrgb.black;
    lp.shcomp = localrgb.shcompr;
    lp.hlcomp = localrgb.hlcompr;
    lp.hlcompthr = localrgb.hlcomprthresh;
    lp.temp = localrgb.temp;
    lp.green = localrgb.green;
    lp.equal = localrgb.equal;
    lp.expcomp = localrgb.expcomp;
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

        // const float radius = lp.rad / (sk * 1.4f); //0 to 70 ==> see skip

        double ave = 0.;
        int n = 0;
        float av = 0;
        int levred;
        bool noiscfactiv = false;

        if (lp.qualmet == 2) { //suppress artifacts with quality enhanced
            levred = 4;
            noiscfactiv = true;
        }    else {
            levred = 7;
            noiscfactiv = false;
        }

// we must here detect : general case, skin, sky,...foliages ???
// delta dhue, luminance and chroma
        constexpr float ared = (rtengine::RT_PI - 0.05f) / 100.f;

        constexpr float bred = 0.05f;

        float dhue = ared * lp.sens + bred; //delta hue lght chroma


        constexpr float maxh = 3.5f; // 3.5 amplification contrast above mean

        constexpr float maxl = 2.5f; // 3 reductio contrast under mean

        float multh = (float) fabs (lp.cont) * (maxh - 1.f) / 100.f + 1.f;

        float mult = (float)fabs (lp.cont) * (maxl - 1.f) / 100.f + 1.f;

        //     lco.dx = 1.f - 1.f / mult;

        //     lco.dy = 1.f - 1.f / multh;

        if ((lp.chro != 0 || lp.ligh != 0.f || lp.expcomp != 0.f) && lp.exposeena) { //interior ellipse renforced lightness and chroma  //locallutili

            float hueplus = hueref + dhue;
            float huemoins = hueref - dhue;

            //printf("hueplus=%f huemoins=%f dhu=%f\n", hueplus, huemoins, dhue);
            if (hueplus > rtengine::RT_PI) {
                hueplus = hueref + dhue - 2.f * rtengine::RT_PI;
            }

            if (huemoins < -rtengine::RT_PI) {
                huemoins = hueref - dhue + 2.f * rtengine::RT_PI;
            }

            LabImage *bufexporig = nullptr;
            Imagefloat* bufworking = nullptr;
            int bfh = 0.f, bfw = 0.f;


            if (call <= 3) { //simpleprocess, dcrop, improccoordinator
                bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
                bfw = int (lp.lx + lp.lxL) + del;
                bufexporig = new LabImage (bfw, bfh);//buffer for data in zone limit
                bufworking = new Imagefloat (bfw, bfh);
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
                        }
                    }
            }


            ImProcFunctions::rgbLocal (bufworking, bufexporig, orirgb, hltonecurveloc, shtonecurveloc, tonecurveloc, lp.chro,
                                       customToneCurve1, customToneCurve2, lp.expcomp, lp.hlcomp, lp.hlcompthr, dcpProf, as);


            if (call <= 3) {

                delete bufworking;
                delete bufexporig;

            }

        }


    }


}


}