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
            int GW = transformed->W;
            int GH = transformed->H;

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
            LabImage *bufexpfin = nullptr;
            Imagefloat* bufworking = nullptr;
            float **buflight = nullptr;
            float **bufl_ab = nullptr;

            int bfh = 0.f, bfw = 0.f;
            int Hd, Wd;
            Hd = GH;
            Wd = GW;


            if (call <= 3) { //simpleprocess, dcrop, improccoordinator
                bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
                bfw = int (lp.lx + lp.lxL) + del;
                Hd = bfh;
                Wd = bfw;
                float clighL = 1.f;
                float clLab = 1.f;


                bufexporig = new LabImage (bfw, bfh);//buffer for data in zone limit
                bufexpfin = new LabImage (bfw, bfh);//buffer for data in zone limit
                bufworking = new Imagefloat (bfw, bfh);

                buflight   = new float*[bfh];//for lightness reti

                for (int i = 0; i < bfh; i++) {
                    buflight[i] = new float[bfw];
                }

                bufl_ab   = new float*[bfh];//for lightness reti

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



                ImProcFunctions::rgbLocal (bufworking, bufexpfin, orirgb, hltonecurveloc, shtonecurveloc, tonecurveloc, lp.chro,
                                           customToneCurve1, customToneCurve2, lp.expcomp, lp.hlcomp, lp.hlcompthr, dcpProf, as);


                //      float maxc = -10000.f;
                //     float minc = 100000.f;

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
                            /*
                                                        if (lighL == 0.f) {
                                                            lighL = 0.001f;
                                                        }

                                                        lL = lighLnew / lighL;

                                                        if (lL <= 1.f) {//convert data near values of slider -100 + 100, to be used after to detection shape
                                                            clighL = 99.f * lL - 99.f;
                                                        } else {
                                                            clighL = CLIPLIG (amplil * lL - amplil); //ampli = 25.f arbitrary empirical coefficient between 5 and 150
                                                        }
                            */
                            /*
                                                                if (clighL > maxc) {
                                                                    maxc = clighL;
                                                                }

                                                                if (clighL < minc) {
                                                                    minc = clighL;
                                                                }
                            */
//clighL = 1.f;
                            float rL;
                            rL = CLIPRET ((bufexpfin->L[loy - begy][lox - begx] - bufexporig->L[loy - begy][lox - begx]) / 328.f);

                            buflight[loy - begy][lox - begx] = rL;


                            float ra;
                            ra = CLIPRET ((sqrt (SQR (bufexpfin->a[loy - begy][lox - begx]) + SQR (bufexpfin->b[loy - begy][lox - begx])) - sqrt (SQR (bufexporig->a[loy - begy][lox - begx]) + SQR (bufexporig->b[loy - begy][lox - begx]))) / 300.f);
                            /*
                                                        if (ra > maxc) {
                                                            maxc = ra;
                                                        }

                                                        if (ra < minc) {
                                                            minc = ra;
                                                        }
                            */
                            //ra = 1.f;
                            bufl_ab[loy - begy][lox - begx] = ra;

                        }
                    }

//               printf ("min=%2.2f max=%2.2f", minc, maxc);
                Expose_Local (call, buflight, bufl_ab, hueplus, huemoins, hueref, dhue, chromaref, lumaref, lp, original, transformed, bufexpfin, cx, cy);

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

void ImProcFunctions::Expose_Local (int call, float **buflight, float **bufchro, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, const struct local_params & lp, LabImage * original, LabImage * transformed, const LabImage * const tmp1, int cx, int cy)
{

//local exposure
    BENCHFUN {
        const float ach = (float)lp.trans / 100.f;

        //chroma
        constexpr float amplchsens = 2.5f;
        constexpr float achsens = (amplchsens - 1.f) / (100.f - 20.f); //20. default locallab.sensih
        constexpr float bchsens = 1.f - 20.f * achsens;
        const float multchro = lp.sens * achsens + bchsens;

        //luma

        //skin
        constexpr float amplchsensskin = 1.6f;
        constexpr float achsensskin = (amplchsensskin - 1.f) / (100.f - 20.f); //20. default locallab.sensih
        constexpr float bchsensskin = 1.f - 20.f * achsensskin;
        const float multchroskin = lp.sens * achsensskin + bchsensskin;

        //transition = difficult to avoid artifact with scope on flat area (sky...)

        constexpr float delhu = 0.1f; //between 0.05 and 0.2

        const float apl = (-1.f) / delhu;
        const float bpl = - apl * hueplus;
        const float amo = 1.f / delhu;
        const float bmo = - amo * huemoins;


        const float pb = 4.f;
        const float pa = (1.f - pb) / 40.f;

        const float ahu = 1.f / (2.8f * lp.sens - 280.f);
        const float bhu = 1.f - ahu * 2.8f * lp.sens;

        const float alum = 1.f / (lp.sens - 100.f);
        const float blum = 1.f - alum * lp.sens;


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

                    if (deltachro < 160.f * SQR (lp.sens / 100.f)) {
                        kch = 1.f;
                    } else {
                        float ck = 160.f * SQR (lp.sens / 100.f);
                        float ak = 1.f / (ck - 160.f);
                        float bk = -160.f * ak;
                        kch = ak * deltachro + bk;
                    }

                    if (lp.sens < 40.f ) {
                        kch = pow (kch, pa * lp.sens + pb);   //increase under 40
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

                    //shape detection for hue chroma and luma
                    if (lp.sens <= 20.f) { //to try...

                        if (deltaE <  2.8f * lp.sens) {
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

                        if (deltaL <  lp.sens) {
                            falu = 1.f;
                        } else {
                            falu = alum * deltaL + blum;
                        }

                    }


                    //    float kdiff = 0.f;
                    // I add these functions...perhaps not good
                    if (kzon) {
                        if (lp.sens < 60.f) { //arbitrary value
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

                            if (lp.sens < 50.f) { //&& lp.chro > 0.f
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
                                transformed->a[y][x] = CLIPC (original->a[y][x] + difa);
                                transformed->b[y][x] = CLIPC (original->b[y][x] + difb);

                            }
                        }

                        //}
                    }
                }
            }
        }

    }
}


}