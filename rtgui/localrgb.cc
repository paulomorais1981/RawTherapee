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
 *  2017 Jacques Desmis <jdesmis@gmail.com>
 */


#include "localrgb.h"
#include "rtimage.h"
#include <iomanip>
#include "../rtengine/rt_math.h"
#include "options.h"
#include <cmath>
#include "edit.h"
#include "guiutils.h"
#include <string>
#include <unistd.h>
#include "../rtengine/improcfun.h"

#define MINTEMP 1500   //1200
#define MAXTEMP 60000  //12000
#define CENTERTEMP 4750
#define MINGREEN 0.02
#define MAXGREEN 10.0
#define MINEQUAL 0.8
#define MAXEQUAL 1.5

using namespace rtengine;
using namespace rtengine::procparams;
extern Options options;

static double wbSlider2Temp (double sval)
{

    // slider range: 0 - 10000
    double temp;

    if (sval <= 5000) {
        // linear below center-temp
        temp = MINTEMP + (sval / 5000.0) * (CENTERTEMP - MINTEMP);
    } else {
        const double slope = (double) (CENTERTEMP - MINTEMP) / (MAXTEMP - CENTERTEMP);
        double x = (sval - 5000) / 5000; // x 0..1
        double y = x * slope + (1.0 - slope) * pow (x, 4.0);
        //double y = pow(x, 4.0);
        temp = CENTERTEMP + y * (MAXTEMP - CENTERTEMP);
    }

    if (temp < MINTEMP) {
        temp = MINTEMP;
    }

    if (temp > MAXTEMP) {
        temp = MAXTEMP;
    }

    return temp;
}

static double wbTemp2Slider (double temp)
{

    double sval;

    if (temp <= CENTERTEMP) {
        sval = ((temp - MINTEMP) / (CENTERTEMP - MINTEMP)) * 5000.0;
    } else {
        const double slope = (double) (CENTERTEMP - MINTEMP) / (MAXTEMP - CENTERTEMP);
        const double y = (temp - CENTERTEMP) / (MAXTEMP - CENTERTEMP);
        double x = pow (y, 0.25); // rough guess of x, will be a little lower
        double k = 0.1;
        bool add = true;

        // the y=f(x) function is a mess to invert, therefore we have this trial-refinement loop instead.
        // from tests, worst case is about 20 iterations, ie no problem
        for (;;) {
            double y1 = x * slope + (1.0 - slope) * pow (x, 4.0);

            if (5000 * fabs (y1 - y) < 0.1) {
                break;
            }

            if (y1 < y) {
                if (!add) {
                    k /= 2;
                }

                x += k;
                add = true;
            } else {
                if (add) {
                    k /= 2;
                }

                x -= k;
                add = false;
            }
        }

        sval = 5000.0 + x * 5000.0;
    }

    if (sval < 0) {
        sval = 0;
    }

    if (sval > 10000) {
        sval = 10000;
    }

    return sval;
}

Localrgb::Localrgb ():
    FoldableToolPanel (this, "localrgb", M ("TP_LOCALRGB_LABEL"), false, true),
    EditSubscriber (ET_OBJECTS), lastObject (-1),
    expexpose (new MyExpander (true, M ("TP_LOCALRGB_EXPO"))),
    expsettings (new MyExpander (false, M ("TP_LOCALLAB_SETTINGS"))),
    expwb (new MyExpander (true, M ("TP_LOCALRGB_WB"))),

    anbspot (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_ANBSPOT"), 0, 1, 1, 0))),
    locX (Gtk::manage (new Adjuster (M ("TP_LOCAL_WIDTH"), 0, 1500, 1, 250))),
    locXL (Gtk::manage (new Adjuster (M ("TP_LOCAL_WIDTH_L"), 0, 1500, 1, 250))),
    degree (Gtk::manage (new Adjuster (M ("TP_LOCAL_DEGREE"), -180, 180, 1, 0))),
    locY (Gtk::manage (new Adjuster (M ("TP_LOCAL_HEIGHT"), 0, 1500, 1, 250))),
    locYT (Gtk::manage (new Adjuster (M ("TP_LOCAL_HEIGHT_T"), 0, 1500, 1, 250))),
    centerX (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_CENTER_X"), -1000, 1000, 1, 0))),
    centerY (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_CENTER_Y"), -1000, 1000, 1, 0))),
    circrad (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_CIRCRADIUS"), 4, 150, 1, 18))),
    thres (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_THRES"), 1, 35, 1, 18))),
    proxi (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_PROXI"), 0, 60, 1, 20))),
    lightness (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_LIGHTNESS"), -100, 100, 1, 0))),
    contrast (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_CONTRAST"), -100, 100, 1, 0))),
    chroma (Gtk::manage (new Adjuster (M ("TP_LOCALRGB_CHROMA"), -100, 100, 1, 0))),
    sensi (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_SENSI"), 0, 100, 1, 19))),
    transit (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_TRANSIT"), 5, 95, 1, 60))),
    retrab (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_RETRAB"), 0, 10000, 1, 500))),
    expcomp (Gtk::manage (new Adjuster (M ("TP_EXPOSURE_EXPCOMP"), -4, 4, 0.05, 0))),
    hlcompr (Gtk::manage (new Adjuster (M ("TP_EXPOSURE_COMPRHIGHLIGHTS"), 0, 500, 1, 0))),
    hlcomprthresh (Gtk::manage (new Adjuster (M ("TP_EXPOSURE_COMPRHIGHLIGHTSTHRESHOLD"), 0, 100, 1, 33))),
    black (Gtk::manage (new Adjuster (M ("TP_EXPOSURE_BLACKLEVEL"), -16384, 32768, 50, 0))),
    shcompr (Gtk::manage (new Adjuster (M ("TP_EXPOSURE_COMPRSHADOWS"), 0, 100, 1, 50))),





    hueref (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_HUEREF"), -3.15, 3.15, 0.01, 0))),
    chromaref (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_CHROMAREF"), 0, 200, 0.01, 0))),
    lumaref (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_LUMAMAREF"), 0, 100, 0.01, 0))),


    Smethod (Gtk::manage (new MyComboBoxText ())),
    qualityMethod (Gtk::manage (new MyComboBoxText ())),
    shapeFrame (Gtk::manage (new Gtk::Frame (M ("TP_LOCALLAB_SHFR")))),
    artifFrame (Gtk::manage (new Gtk::Frame (M ("TP_LOCALLAB_ARTIF")))),
    superFrame (Gtk::manage (new Gtk::Frame ())),

    //  artifVBox (Gtk::manage (new Gtk::VBox ())),
    //  shapeVBox (Gtk::manage (new Gtk::VBox ())),
    //  colorVBox (Gtk::manage ( new Gtk::VBox())),

    labqual (Gtk::manage (new Gtk::Label (M ("TP_LOCALLAB_QUAL_METHOD") + ":"))),

    labmS (Gtk::manage (new Gtk::Label (M ("TP_LOCALLAB_STYPE") + ":"))),
    ctboxS (Gtk::manage (new Gtk::HBox ())),
    qualbox (Gtk::manage (new Gtk::HBox ())),
    draggedPointOldAngle (-1000.)

{
    CurveListener::setMulti (true);

    std::vector<GradientMilestone> bottomMilestones;
    bottomMilestones.push_back ( GradientMilestone (0., 0., 0., 0.) );
    bottomMilestones.push_back ( GradientMilestone (1., 1., 1., 1.) );

    ProcParams params;
    editHBox = Gtk::manage (new Gtk::HBox());
    edit = Gtk::manage (new Gtk::ToggleButton());
    edit->add (*Gtk::manage (new RTImage ("editmodehand.png")));
    edit->set_tooltip_text (M ("EDIT_OBJECT_TOOLTIP"));
    editConn = edit->signal_toggled().connect ( sigc::mem_fun (*this, &Localrgb::editToggled) );
    editHBox->pack_start (*edit, Gtk::PACK_SHRINK, 0);
    pack_start (*editHBox, Gtk::PACK_SHRINK, 0);
    int realnbspot;


    realnbspot = options.rtSettings.nspot;

    nbspot = Gtk::manage (new Adjuster (M ("TP_LOCALLAB_NBSPOT"), 1, realnbspot, 1, 1));

    if (options.rtSettings.locdelay) {

        if (nbspot->delay < 200) {
            nbspot->delay = 200;
        }
    }


    nbspot->setAdjusterListener (this);
    nbspot->set_tooltip_text (M ("TP_LOCALLAB_NBSPOT_TOOLTIP"));


    anbspot->setAdjusterListener (this);
    anbspot->set_tooltip_text (M ("TP_LOCALLAB_ANBSPOT_TOOLTIP"));
    retrab->setAdjusterListener (this);

    shapeFrame->set_label_align (0.025, 0.5);

    expsettings->signal_button_release_event().connect_notify ( sigc::bind ( sigc::mem_fun (this, &Localrgb::foldAllButMe), expsettings) );


    expexpose->signal_button_release_event().connect_notify ( sigc::bind ( sigc::mem_fun (this, &Localrgb::foldAllButMe), expexpose) );
    enableexposeConn = expexpose->signal_enabled_toggled().connect ( sigc::bind ( sigc::mem_fun (this, &Localrgb::enableToggled), expexpose) );

    expwb->signal_button_release_event().connect_notify ( sigc::bind ( sigc::mem_fun (this, &Localrgb::foldAllButMe), expwb) );
    enablewbConn = expwb->signal_enabled_toggled().connect ( sigc::bind ( sigc::mem_fun (this, &Localrgb::enableToggled), expwb) );

    ctboxS->pack_start (*labmS, Gtk::PACK_SHRINK, 4);
    ctboxS->set_tooltip_markup (M ("TP_LOCALLAB_STYPE_TOOLTIP"));

    Smethod->append (M ("TP_LOCALLAB_IND"));
    Smethod->append (M ("TP_LOCALLAB_SYM"));
    Smethod->append (M ("TP_LOCALLAB_INDSL"));
    Smethod->append (M ("TP_LOCALLAB_SYMSL"));
    Smethod->set_active (0);
    Smethodconn = Smethod->signal_changed().connect ( sigc::mem_fun (*this, &Localrgb::SmethodChanged) );

    //locX->set_tooltip_text (M("TP_LOCAL_WIDTH_TOOLTIP"));
    locX->setAdjusterListener (this);

    //locX->set_tooltip_text (M("TP_LOCAL_WIDTH_TOOLTIP"));
    locXL->setAdjusterListener (this);

    //degree->set_tooltip_text (M("TP_LOCAL_DEGREE_TOOLTIP"));
    degree->setAdjusterListener (this);

    //locY->set_tooltip_text (M("TP_LOCAL_HEIGHT_TOOLTIP"));
    locY->setAdjusterListener (this);

    //locY->set_tooltip_text (M("TP_LOCAL_HEIGHT_TOOLTIP"));
    locYT->setAdjusterListener (this);

    //centerX->set_tooltip_text (M("TP_LOCALLAB_CENTER_X_TOOLTIP"));
    centerX->setAdjusterListener (this);

    //centerY->set_tooltip_text (M("TP_LOCALLAB_CENTER_Y_TOOLTIP"));
    centerY->setAdjusterListener (this);

    circrad->setAdjusterListener (this);
    thres->setAdjusterListener (this);

    proxi->setAdjusterListener (this);


    qualityMethod->append (M ("TP_LOCALLAB_STD"));
    qualityMethod->append (M ("TP_LOCALLAB_ENH"));
    qualityMethod->append (M ("TP_LOCALLAB_ENHDEN"));
    qualityMethod->set_active (0);
    qualityMethodConn = qualityMethod->signal_changed().connect ( sigc::mem_fun (*this, &Localrgb::qualityMethodChanged) );
    qualityMethod->set_tooltip_markup (M ("TP_LOCALLAB_METHOD_TOOLTIP"));

    expcomp->setAdjusterListener (this);
    hlcomprthresh->setAdjusterListener (this);
    black->setAdjusterListener (this);
    hlcompr->setAdjusterListener (this);
    shcompr->setAdjusterListener (this);

    //lightness->set_tooltip_text (M("TP_LOCALLAB_LIGHTNESS_TOOLTIP"));
    lightness->setAdjusterListener (this);

    //contrast->set_tooltip_text (M("TP_LOCALLAB_CONTRAST_TOOLTIP"));
    contrast->setAdjusterListener (this);

    //chroma->set_tooltip_text (M("TP_LOCALLAB_CHROMA_TOOLTIP"));
    chroma->setAdjusterListener (this);

    sensi->set_tooltip_text (M ("TP_LOCALLAB_SENSI_TOOLTIP"));
    sensi->setAdjusterListener (this);

    transit->set_tooltip_text (M ("TP_LOCALLAB_TRANSIT_TOOLTIP"));
    transit->setAdjusterListener (this);


    ToolParamBlock* const shapeBox = Gtk::manage (new ToolParamBlock());

    shapeBox->pack_start (*nbspot);
    pack_start (*anbspot);

    hueref->setAdjusterListener (this);
    chromaref->setAdjusterListener (this);
    lumaref->setAdjusterListener (this);

    pack_start (*hueref);
    pack_start (*chromaref);
    pack_start (*lumaref);
    /*
        anbspot->hide();//keep anbspot  - i used it to test diffrent algo...
        hueref->hide();
        chromaref->hide();
        lumaref->hide();
    */
    ctboxS->pack_start (*Smethod);
    shapeBox->pack_start (*ctboxS);
    shapeBox->pack_start (*locX);
    shapeBox->pack_start (*locXL);
    //pack_start (*degree);
    shapeBox->pack_start (*locY);
    shapeBox->pack_start (*locYT);
    shapeBox->pack_start (*centerX);
    shapeBox->pack_start (*centerY);
    shapeBox->pack_start (*circrad);
    qualbox->pack_start (*labqual, Gtk::PACK_SHRINK, 4);
    qualbox->pack_start (*qualityMethod);
    shapeBox->pack_start (*qualbox);
    shapeBox->pack_start (*transit);

    artifFrame->set_label_align (0.025, 0.5);
    artifFrame->set_tooltip_text (M ("TP_LOCALLAB_ARTIF_TOOLTIP"));

    ToolParamBlock* const artifBox = Gtk::manage (new ToolParamBlock());

    artifBox->pack_start (*thres);
    artifBox->pack_start (*proxi);
    artifBox->pack_start (*retrab);

    artifFrame->add (*artifBox);
    shapeBox->pack_start (*artifFrame);

    expsettings->add (*shapeBox);
    expsettings->setLevel (2);
    pack_start (*expsettings);

    superFrame->set_label_align (0.025, 0.5);
//   Gtk::VBox *superVBox = Gtk::manage ( new Gtk::VBox());
//   superFrame->set_label_widget (*curvactiv);
    ToolParamBlock* const superBox = Gtk::manage (new ToolParamBlock());
    ToolParamBlock* const colorBox = Gtk::manage (new ToolParamBlock());

    colorBox->pack_start (*expcomp);
    colorBox->pack_start (*hlcompr);
    colorBox->pack_start (*hlcomprthresh);
    colorBox->pack_start (*black);
    colorBox->pack_start (*shcompr);
    colorBox->pack_start (*lightness);
    colorBox->pack_start (*contrast);
    colorBox->pack_start (*chroma);


    toneCurveMode = Gtk::manage (new MyComboBoxText ());
    toneCurveMode->append (M ("TP_EXPOSURE_TCMODE_STANDARD"));
    toneCurveMode->append (M ("TP_EXPOSURE_TCMODE_WEIGHTEDSTD"));
    toneCurveMode->append (M ("TP_EXPOSURE_TCMODE_FILMLIKE"));
    toneCurveMode->append (M ("TP_EXPOSURE_TCMODE_SATANDVALBLENDING"));
    toneCurveMode->append (M ("TP_EXPOSURE_TCMODE_LUMINANCE"));
    toneCurveMode->append (M ("TP_EXPOSURE_TCMODE_PERCEPTUAL"));
    toneCurveMode->set_active (0);
    toneCurveMode->set_tooltip_text (M ("TP_EXPOSURE_TCMODE_LABEL1"));

    curveEditorG = new CurveEditorGroup (options.lastToneCurvesDir, M ("TP_EXPOSURE_CURVEEDITOR1"));
    curveEditorG->setCurveListener (this);

    shape = static_cast<DiagonalCurveEditor*> (curveEditorG->addCurve (CT_Diagonal, "", toneCurveMode));
//   shape->setEditID(EUID_ToneCurve1, BT_IMAGEFLOAT);
    shape->setBottomBarBgGradient (bottomMilestones);
    shape->setLeftBarBgGradient (bottomMilestones);

    // This will add the reset button at the end of the curveType buttons
    curveEditorG->curveListComplete();

    colorBox->pack_start ( *curveEditorG, Gtk::PACK_SHRINK, 2);

    tcmodeconn = toneCurveMode->signal_changed().connect ( sigc::mem_fun (*this, &Localrgb::curveMode1Changed), true );

    toneCurveMode2 = Gtk::manage (new MyComboBoxText ());
    toneCurveMode2->append (M ("TP_EXPOSURE_TCMODE_STANDARD"));
    toneCurveMode2->append (M ("TP_EXPOSURE_TCMODE_WEIGHTEDSTD"));
    toneCurveMode2->append (M ("TP_EXPOSURE_TCMODE_FILMLIKE"));
    toneCurveMode2->append (M ("TP_EXPOSURE_TCMODE_SATANDVALBLENDING"));
    toneCurveMode2->append (M ("TP_EXPOSURE_TCMODE_LUMINANCE"));
    toneCurveMode2->append (M ("TP_EXPOSURE_TCMODE_PERCEPTUAL"));
    toneCurveMode2->set_active (0);
    toneCurveMode2->set_tooltip_text (M ("TP_EXPOSURE_TCMODE_LABEL2"));

    curveEditorG2 = new CurveEditorGroup (options.lastToneCurvesDir, M ("TP_EXPOSURE_CURVEEDITOR2"));
    curveEditorG2->setCurveListener (this);

    shape2 = static_cast<DiagonalCurveEditor*> (curveEditorG2->addCurve (CT_Diagonal, "", toneCurveMode2));
//   shape2->setEditID(EUID_ToneCurve2, BT_IMAGEFLOAT);
    shape2->setBottomBarBgGradient (bottomMilestones);
    shape2->setLeftBarBgGradient (bottomMilestones);

    // This will add the reset button at the end of the curveType buttons
    curveEditorG2->curveListComplete();
    curveEditorG2->setTooltip (M ("TP_EXPOSURE_CURVEEDITOR2_TOOLTIP"));

    colorBox->pack_start ( *curveEditorG2, Gtk::PACK_SHRINK, 2);

    tcmode2conn = toneCurveMode2->signal_changed().connect ( sigc::mem_fun (*this, &Localrgb::curveMode2Changed), true );

    superBox->pack_start (*sensi);
    superFrame->add (*superBox);
    colorBox->pack_start (*superFrame);

    expexpose->add (*colorBox);
    expexpose->setLevel (2);
    pack_start (*expexpose);


    ToolParamBlock* const wbBox = Gtk::manage (new ToolParamBlock());
    /*
        Gtk::HBox* spotbox = Gtk::manage (new Gtk::HBox ());
        spotbox->set_spacing(4);
        spotbox->show ();

        spotbutton = Gtk::manage (new Gtk::Button ());
        spotbutton->set_tooltip_text(M("TP_WBALANCE_SPOTWB"));
        Gtk::Image* spotimg = Gtk::manage (new RTImage ("gtk-color-picker-small.png"));
        spotimg->show ();
        spotbutton->set_image (*spotimg);
        spotbutton->show ();

        spotbox->pack_start (*spotbutton);
    */
    Gtk::Image* itempL =  Gtk::manage (new RTImage ("ajd-wb-temp1.png"));
    Gtk::Image* itempR =  Gtk::manage (new RTImage ("ajd-wb-temp2.png"));
    Gtk::Image* igreenL = Gtk::manage (new RTImage ("ajd-wb-green1.png"));
    Gtk::Image* igreenR = Gtk::manage (new RTImage ("ajd-wb-green2.png"));
    Gtk::Image* iblueredL = Gtk::manage (new RTImage ("ajd-wb-bluered1.png"));
    Gtk::Image* iblueredR = Gtk::manage (new RTImage ("ajd-wb-bluered2.png"));

    temp = Gtk::manage (new Adjuster (M ("TP_WBALANCE_TEMPERATURE"), MINTEMP, MAXTEMP, 5, CENTERTEMP, itempL, itempR, &wbSlider2Temp, &wbTemp2Slider));
    green = Gtk::manage (new Adjuster (M ("TP_WBALANCE_GREEN"), MINGREEN, MAXGREEN, 0.001, 1.0, igreenL, igreenR));
    equal = Gtk::manage (new Adjuster (M ("TP_WBALANCE_EQBLUERED"), MINEQUAL, MAXEQUAL, 0.001, 1.0, iblueredL, iblueredR));

    temp->show ();
    green->show ();
    equal->show ();
//  wbBox->pack_start (*spotbox);

    wbBox->pack_start (*temp);
    wbBox->pack_start (*green);
    wbBox->pack_start (*equal);

    temp->setAdjusterListener (this);
    green->setAdjusterListener (this);
    equal->setAdjusterListener (this);



    expwb->add (*wbBox);
    expwb->setLevel (2);
    pack_start (*expwb);

    // Instantiating the Editing geometry; positions will be initialized later
    Line *locYLine[2], *locXLine[2];
    Circle *centerCircle;

    Beziers *onebeziers[4] = {};
    Beziers *twobeziers[4] = {};
    Beziers *thrbeziers[4] = {};
    Beziers *foubeziers[4] = {};
    float innw = 0.7f;
    // Visible geometry
    locXLine[0] = new Line();
    locXLine[0]->innerLineWidth = 2;
    locXLine[1] = new Line();
    locXLine[1]->innerLineWidth = 2;
    locXLine[0]->datum  = locXLine[1]->datum = Geometry::IMAGE;

    locYLine[0] = new Line();
    locYLine[0]->innerLineWidth = 2;
    locYLine[1] = new Line();
    locYLine[1]->innerLineWidth = 2;
    locYLine[0]->datum = locYLine[1]->datum = Geometry::IMAGE;

    centerCircle = new Circle();
    centerCircle->datum = Geometry::IMAGE;
    centerCircle->radiusInImageSpace = true;
    centerCircle->radius = circrad->getValue(); //19;
    centerCircle->filled = false;

    if (options.showdelimspot) {
        onebeziers[0] = new Beziers();
        onebeziers[0]->datum = Geometry::IMAGE;
        onebeziers[0]->innerLineWidth = innw;

        onebeziers[1] = new Beziers();
        onebeziers[1]->datum = Geometry::IMAGE;
        onebeziers[1]->innerLineWidth = innw;

        onebeziers[2] = new Beziers();
        onebeziers[2]->datum = Geometry::IMAGE;
        onebeziers[2]->innerLineWidth = innw;

        onebeziers[3] = new Beziers();
        onebeziers[3]->datum = Geometry::IMAGE;
        onebeziers[3]->innerLineWidth = innw;

        twobeziers[0] = new Beziers();
        twobeziers[0]->datum = Geometry::IMAGE;
        twobeziers[0]->innerLineWidth = innw;

        twobeziers[1] = new Beziers();
        twobeziers[1]->datum = Geometry::IMAGE;
        twobeziers[1]->innerLineWidth = innw;

        twobeziers[2] = new Beziers();
        twobeziers[2]->datum = Geometry::IMAGE;
        twobeziers[2]->innerLineWidth = innw;

        twobeziers[3] = new Beziers();
        twobeziers[3]->datum = Geometry::IMAGE;
        twobeziers[3]->innerLineWidth = innw;

        thrbeziers[0] = new Beziers();
        thrbeziers[0]->datum = Geometry::IMAGE;
        thrbeziers[0]->innerLineWidth = innw;

        thrbeziers[1] = new Beziers();
        thrbeziers[1]->datum = Geometry::IMAGE;
        thrbeziers[1]->innerLineWidth = innw;

        thrbeziers[2] = new Beziers();
        thrbeziers[2]->datum = Geometry::IMAGE;
        thrbeziers[2]->innerLineWidth = innw;

        thrbeziers[3] = new Beziers();
        thrbeziers[3]->datum = Geometry::IMAGE;
        thrbeziers[3]->innerLineWidth = innw;

        foubeziers[0] = new Beziers();
        foubeziers[0]->datum = Geometry::IMAGE;
        foubeziers[0]->innerLineWidth = innw;

        foubeziers[1] = new Beziers();
        foubeziers[1]->datum = Geometry::IMAGE;
        foubeziers[1]->innerLineWidth = innw;

        foubeziers[2] = new Beziers();
        foubeziers[2]->datum = Geometry::IMAGE;
        foubeziers[2]->innerLineWidth = innw;

        foubeziers[3] = new Beziers();
        foubeziers[3]->datum = Geometry::IMAGE;
        foubeziers[3]->innerLineWidth = innw;

    }

    // oneellipse->radiusInImageSpace = true;
    // oneellipse->radius = locX->getValue();
    // oneellipse->filled = false;

    EditSubscriber::visibleGeometry.push_back ( locXLine[0] );
    EditSubscriber::visibleGeometry.push_back ( locXLine[1] );
    EditSubscriber::visibleGeometry.push_back ( locYLine[0] );
    EditSubscriber::visibleGeometry.push_back ( locYLine[1] );
    EditSubscriber::visibleGeometry.push_back ( centerCircle );

    if (options.showdelimspot) {
        EditSubscriber::visibleGeometry.push_back ( onebeziers[0] );
        EditSubscriber::visibleGeometry.push_back ( onebeziers[1] );
        EditSubscriber::visibleGeometry.push_back ( onebeziers[2] );
        EditSubscriber::visibleGeometry.push_back ( onebeziers[3] );
        EditSubscriber::visibleGeometry.push_back ( twobeziers[0] );
        EditSubscriber::visibleGeometry.push_back ( twobeziers[1] );
        EditSubscriber::visibleGeometry.push_back ( twobeziers[2] );
        EditSubscriber::visibleGeometry.push_back ( twobeziers[3] );
        EditSubscriber::visibleGeometry.push_back ( thrbeziers[0] );
        EditSubscriber::visibleGeometry.push_back ( thrbeziers[1] );
        EditSubscriber::visibleGeometry.push_back ( thrbeziers[2] );
        EditSubscriber::visibleGeometry.push_back ( thrbeziers[3] );
        EditSubscriber::visibleGeometry.push_back ( foubeziers[0] );
        EditSubscriber::visibleGeometry.push_back ( foubeziers[1] );
        EditSubscriber::visibleGeometry.push_back ( foubeziers[2] );
        EditSubscriber::visibleGeometry.push_back ( foubeziers[3] );
    }

    // MouseOver geometry
    locXLine[0] = new Line();
    locXLine[0]->innerLineWidth = 2;
    locXLine[1] = new Line();
    locXLine[1]->innerLineWidth = 2;
    locXLine[0]->datum  = locXLine[1]->datum = Geometry::IMAGE;

    locYLine[0] = new Line();
    locYLine[0]->innerLineWidth = 2;
    locYLine[1] = new Line();
    locYLine[1]->innerLineWidth = 2;
    locYLine[0]->datum = locYLine[1]->datum = Geometry::IMAGE;

    centerCircle = new Circle();
    centerCircle->datum = Geometry::IMAGE;
    centerCircle->radiusInImageSpace = true;
    centerCircle->radius = circrad->getValue();//19;
    centerCircle->filled = true;

    if (options.showdelimspot) {
        onebeziers[0]   = new Beziers();
        onebeziers[0]->datum = Geometry::IMAGE;
        onebeziers[0]->innerLineWidth = innw;

        onebeziers[1]   = new Beziers();
        onebeziers[1]->datum = Geometry::IMAGE;
        onebeziers[1]->innerLineWidth = innw;

        onebeziers[2]   = new Beziers();
        onebeziers[2]->datum = Geometry::IMAGE;
        onebeziers[2]->innerLineWidth = innw;

        onebeziers[3]   = new Beziers();
        onebeziers[3]->datum = Geometry::IMAGE;
        onebeziers[3]->innerLineWidth = innw;

        twobeziers[0] = new Beziers();
        twobeziers[0]->datum = Geometry::IMAGE;
        twobeziers[0]->innerLineWidth = innw;

        twobeziers[1] = new Beziers();
        twobeziers[1]->datum = Geometry::IMAGE;
        twobeziers[1]->innerLineWidth = innw;

        twobeziers[2] = new Beziers();
        twobeziers[2]->datum = Geometry::IMAGE;
        twobeziers[2]->innerLineWidth = innw;

        twobeziers[3] = new Beziers();
        twobeziers[3]->datum = Geometry::IMAGE;
        twobeziers[3]->innerLineWidth = innw;

        thrbeziers[0] = new Beziers();
        thrbeziers[0]->datum = Geometry::IMAGE;
        thrbeziers[0]->innerLineWidth = innw;

        thrbeziers[1] = new Beziers();
        thrbeziers[1]->datum = Geometry::IMAGE;
        thrbeziers[1]->innerLineWidth = innw;

        thrbeziers[2] = new Beziers();
        thrbeziers[2]->datum = Geometry::IMAGE;
        thrbeziers[2]->innerLineWidth = innw;

        thrbeziers[3] = new Beziers();
        thrbeziers[3]->datum = Geometry::IMAGE;
        thrbeziers[3]->innerLineWidth = innw;

        foubeziers[0] = new Beziers();
        foubeziers[0]->datum = Geometry::IMAGE;
        foubeziers[0]->innerLineWidth = innw;

        foubeziers[1] = new Beziers();
        foubeziers[1]->datum = Geometry::IMAGE;
        foubeziers[1]->innerLineWidth = innw;

        foubeziers[2] = new Beziers();
        foubeziers[2]->datum = Geometry::IMAGE;
        foubeziers[2]->innerLineWidth = innw;

        foubeziers[3] = new Beziers();
        foubeziers[3]->datum = Geometry::IMAGE;
        foubeziers[3]->innerLineWidth = innw;
    }

//   oneellipse->radiusInImageSpace = true;
//   oneellipse->radius = 10;//locX->getValue();
//    oneellipse->filled = false;

    EditSubscriber::mouseOverGeometry.push_back ( locXLine[0] );
    EditSubscriber::mouseOverGeometry.push_back ( locXLine[1] );

    EditSubscriber::mouseOverGeometry.push_back ( locYLine[0] );
    EditSubscriber::mouseOverGeometry.push_back ( locYLine[1] );

    EditSubscriber::mouseOverGeometry.push_back ( centerCircle );

    if (options.showdelimspot) {
        EditSubscriber::mouseOverGeometry.push_back ( onebeziers[0] );
        EditSubscriber::mouseOverGeometry.push_back ( onebeziers[1] );
        EditSubscriber::mouseOverGeometry.push_back ( onebeziers[2] );
        EditSubscriber::mouseOverGeometry.push_back ( onebeziers[3] );
        EditSubscriber::mouseOverGeometry.push_back ( twobeziers[0] );
        EditSubscriber::mouseOverGeometry.push_back ( twobeziers[1] );
        EditSubscriber::mouseOverGeometry.push_back ( twobeziers[2] );
        EditSubscriber::mouseOverGeometry.push_back ( twobeziers[3] );
        EditSubscriber::mouseOverGeometry.push_back ( thrbeziers[0] );
        EditSubscriber::mouseOverGeometry.push_back ( thrbeziers[1] );
        EditSubscriber::mouseOverGeometry.push_back ( thrbeziers[2] );
        EditSubscriber::mouseOverGeometry.push_back ( thrbeziers[3] );
        EditSubscriber::mouseOverGeometry.push_back ( foubeziers[0] );
        EditSubscriber::mouseOverGeometry.push_back ( foubeziers[1] );
        EditSubscriber::mouseOverGeometry.push_back ( foubeziers[2] );
        EditSubscriber::mouseOverGeometry.push_back ( foubeziers[3] );
    }

    show_all();
}

Localrgb::~Localrgb()
{
    for (std::vector<Geometry*>::const_iterator i = visibleGeometry.begin(); i != visibleGeometry.end(); ++i) {
        delete *i;
    }

    for (std::vector<Geometry*>::const_iterator i = mouseOverGeometry.begin(); i != mouseOverGeometry.end(); ++i) {
        delete *i;
    }

    delete curveEditorG;
    delete curveEditorG2;
}
void Localrgb::foldAllButMe (GdkEventButton* event, MyExpander *expander)
{
    if (event->button == 3) {
        expsettings->set_expanded (expsettings == expander);
        expexpose->set_expanded (expexpose == expander);
        expwb->set_expanded (expwb == expander);
    }
}

void Localrgb::enableToggled (MyExpander *expander)
{
    if (listener) {
        rtengine::ProcEvent event = NUMOFEVENTS;

        if (expander == expexpose) {
            event = EvLocrgbenaexpose;
        } else if (expander == expwb) {
            event = EvLocrgbenawb;

        } else {
            return;
        }

        if (expander->get_inconsistent()) {
            listener->panelChanged (event, M ("GENERAL_UNCHANGED"));
        } else if (expander->getEnabled()) {
            listener->panelChanged (event, M ("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (event, M ("GENERAL_DISABLED"));
        }

    }
}

void Localrgb::writeOptions (std::vector<int> &tpOpen)
{
    tpOpen.push_back (expsettings->get_expanded ());
    tpOpen.push_back (expexpose->get_expanded ());
    tpOpen.push_back (expwb->get_expanded ());
}

void Localrgb::updateToolState (std::vector<int> &tpOpen)
{
    if (tpOpen.size() == 3) {
        expsettings->set_expanded (tpOpen.at (0));
        expexpose->set_expanded (tpOpen.at (1));
        expwb->set_expanded (tpOpen.at (2));

    }
}

/*
void Localrgb::autoOpenCurve ()
{
    cTgainshape->openIfNonlinear();
    //  llshape->openIfNonlinear();
    //  LHshape->openIfNonlinear();

}


int localChangedUI (void* data)
{

    GThreadLock lock;
    (static_cast<Localrgb*> (data))->localComputed_ ();

    return 0;
}

int localretChangedUI (void* data)
{

    GThreadLock lock;
    (static_cast<Localrgb*> (data))->localretComputed_ ();

    return 0;
}

bool Localrgb::localretComputed_ ()
{
    disableListener ();

    //Reticurv
//update GUI and MIP specially for curve

    int *s_datc;
    s_datc = new int[70];
    int siz;
    //printf("nexts=%s\n", nextstr2.c_str());
    ImProcFunctions::strcurv_data (nextstr2, s_datc, siz);
    std::vector<double>   creti;

    for (int j = 0; j < siz; j++) {
        creti.push_back ((double) (s_datc[j]) / 1000.);
    }

    delete [] s_datc;

    cTgainshape->setCurve (creti);

    int *s_datcl;
    s_datcl = new int[70];
    int sizl;
    ImProcFunctions::strcurv_data (nextll_str2, s_datcl, sizl);
    std::vector<double>   cll;

    for (int j = 0; j < sizl; j++) {
        cll.push_back ((double) (s_datcl[j]) / 1000.);
    }

    delete [] s_datcl;

    llshape->setCurve (cll);


    int *s_datcc;
    s_datcc = new int[70];
    int sizc;
    ImProcFunctions::strcurv_data (nextcc_str2, s_datcc, sizc);
    std::vector<double>   ccc;

    for (int j = 0; j < sizc; j++) {
        ccc.push_back ((double) (s_datcc[j]) / 1000.);
    }

    delete [] s_datcc;

    ccshape->setCurve (ccc);

    int *s_datch;
    s_datch = new int[70];
    int sizh;
    ImProcFunctions::strcurv_data (nextlh_str2, s_datch, sizh);
    std::vector<double>   clh;

    for (int j = 0; j < sizh; j++) {
        clh.push_back ((double) (s_datch[j]) / 1000.);
    }

    delete [] s_datch;

    LHshape->setCurve (clh);


    enableListener ();

    //update all sliders by this strange process!
    if (anbspot->getValue() == 0) {
        anbspot->setValue (1);

        if (options.rtSettings.locdelay) {
            if (anbspot->delay < 100) {
                anbspot->delay = 100;
            }
        }

        adjusterChanged (anbspot, 1);

    } else if (anbspot->getValue() == 1) {
        anbspot->setValue (0);

        if (options.rtSettings.locdelay) {
            if (anbspot->delay < 100) {
                anbspot->delay = 100;
            }
        }

        adjusterChanged (anbspot, 0);

    }

    //update all curves
    std::vector<double>   cretirab;
    cretirab = cTgainshaperab->getCurve ();

    if (cretirab.at (5) == 0.70) {
        cretirab.at (5) = 0.9;
        cTgainshaperab->setCurve (cretirab);

        curveChanged (cTgainshaperab);
    } else if (cretirab.at (5) == 0.90) {
        cretirab.at (5) = 0.7;
        cTgainshaperab->setCurve (cretirab);
        curveChanged (cTgainshaperab);

    }


//    printf("G2 anbspot=%i\n", anbspot->getValue());

    if (listener) { //for all sliders
        listener->panelChanged (Evlocallabanbspot, "");//anbspot->getTextValue());
    }

    if (listener) {//for curve
        listener->panelChanged (EvlocallabCTgainCurverab, M (""));
    }

    if (listener) {//for curve
        listener->panelChanged (EvlocallabCTgainCurve, M (""));
    }

    if (listener) {//for curve
        listener->panelChanged (Evlocallabllshape, M (""));
    }

    if (listener) {//for curve
        listener->panelChanged (Evlocallabccshape, M (""));
    }

    if (listener) {//for curve
        listener->panelChanged (EvlocallabLHshape, M (""));
    }


}
*/

/*
bool Localrgb::localComputed_ ()
{
//update GUI and MIP
    disableListener ();

    //size spot
    circrad->setValue (nextdatasp[2]);
    //center and cursor
    locX->setValue (nextdatasp[3]);
    locY->setValue (nextdatasp[4]);
    locYT->setValue (nextdatasp[5]);
    locXL->setValue (nextdatasp[6]);
    centerX->setValue (nextdatasp[7]);
    centerY->setValue (nextdatasp[8]);

    //sliders
    lightness->setValue (nextdatasp[9]);
    contrast->setValue (nextdatasp[10]);
    chroma->setValue (nextdatasp[11]);
    sensi->setValue (nextdatasp[12]);
    transit->setValue (nextdatasp[13]);

    //inverse
    if (nextdatasp[14] == 0) {
        invers->set_active (false);
    } else {
        invers->set_active (true);
    }

    //method cursor
    if (nextdatasp[15] == 0) {
        Smethod->set_active (0);
    } else if (nextdatasp[15] == 1) {
        Smethod->set_active (1);
    } else if (nextdatasp[15] == 2) {
        Smethod->set_active (2);
    } else if (nextdatasp[15] == 3) {
        Smethod->set_active (3);
    }

    //sliders blurr
    radius->setValue (nextdatasp[17]);
    strength->setValue (nextdatasp[18]);
    sensibn->setValue (nextdatasp[19]);

    //inverse
    if (nextdatasp[20] == 0) {
        inversrad->set_active (false);
    } else {
        inversrad->set_active (true);
    }

    //sliders retinex
    str->setValue (nextdatasp[21]);
    chrrt->setValue (nextdatasp[22]);
    neigh->setValue (nextdatasp[23]);
    vart->setValue (nextdatasp[24]);
    sensih->setValue (nextdatasp[25]);

    //inverse
    if (nextdatasp[26] == 0) {
        inversret->set_active (false);
    } else {
        inversret->set_active (true);
    }

    //method retinex
    if (nextdatasp[27] == 0) {
        retinexMethod->set_active (0);
    } else if (nextdatasp[27] == 1) {
        retinexMethod->set_active (1);
    } else if (nextdatasp[27] == 2) {
        retinexMethod->set_active (2);
    }

    //sharpening
    sharradius->setValue (nextdatasp[28]);
    sharamount->setValue (nextdatasp[29]);
    shardamping->setValue (nextdatasp[30]);
    shariter->setValue (nextdatasp[31]);
    sensisha->setValue (nextdatasp[32]);

    if (nextdatasp[33] == 0) {
        inverssha->set_active (false);
    } else {
        inverssha->set_active (true);
    }

    if (nextdatasp[34] == 0) {
        qualityMethod->set_active (0);
    } else if (nextdatasp[34] == 1) {
        qualityMethod->set_active (1);
    } else if (nextdatasp[34] == 2) {
        qualityMethod->set_active (2);
    }

    thres->setValue (nextdatasp[35]);
    proxi->setValue (nextdatasp[36]);

    //denoise
    noiselumf->setValue (nextdatasp[37]);
    noiselumc->setValue (nextdatasp[38]);
    noisechrof->setValue (nextdatasp[39]);
    noisechroc->setValue (nextdatasp[40]);

    //cbdl
    multiplier[0]->setValue (nextdatasp[41]);
    multiplier[1]->setValue (nextdatasp[42]);
    multiplier[2]->setValue (nextdatasp[43]);
    multiplier[3]->setValue (nextdatasp[44]);
    multiplier[4]->setValue (nextdatasp[45]);
    threshold->setValue (nextdatasp[46]);
    sensicb->setValue (nextdatasp[47]);

    //blur luma
    if (nextdatasp[48] == 0) {
        activlum->set_active (false);
    } else {
        activlum->set_active (true);
    }

//TM
    stren->setValue (nextdatasp[49]);
    gamma->setValue (nextdatasp[50]);
    estop->setValue (nextdatasp[51]);
    scaltm->setValue (nextdatasp[52]);
    rewei->setValue (nextdatasp[53]);
    sensitm->setValue (nextdatasp[54]);
    //  usleep(10000);

    //Reticurv
    retrab->setValue (nextdatasp[55]);

    //curvactiv
    if (nextdatasp[56] == 0) {
        curvactiv->set_active (false);
    } else {
        curvactiv->set_active (true);
    }

    if (nextdatasp[57] == 0) {
        qualitycurveMethod->set_active (0);
    } else if (nextdatasp[57] == 1) {
        qualitycurveMethod->set_active (1);
    } else if (nextdatasp[57] == 2) {
        qualitycurveMethod->set_active (2);
    }

    double intermed = 0.01 * (double) nextdatasp[58];
    hueref->setValue (intermed);
    chromaref->setValue (nextdatasp[59]);
    lumaref->setValue (nextdatasp[60]);

    int *s_datc;
    s_datc = new int[70];
    int siz;
    ImProcFunctions::strcurv_data (nextstr, s_datc, siz);


    std::vector<double>   creti;

    for (int j = 0; j < siz; j++) {
        creti.push_back ((double) (s_datc[j]) / 1000.);
    }

    delete [] s_datc;

    cTgainshape->setCurve (creti);

    //LLcurv
    int *s_datcl;
    s_datcl = new int[70];
    int sizl;
    ImProcFunctions::strcurv_data (nextll_str, s_datcl, sizl);


    std::vector<double>   cll;

    for (int j = 0; j < sizl; j++) {
        cll.push_back ((double) (s_datcl[j]) / 1000.);
    }

    delete [] s_datcl;
    llshape->setCurve (cll);

    //CCcurv
    int *s_datcc;
    s_datcc = new int[70];
    int sizc;
    ImProcFunctions::strcurv_data (nextcc_str, s_datcc, sizc);


    std::vector<double>   ccc;

    for (int j = 0; j < sizc; j++) {
        ccc.push_back ((double) (s_datcc[j]) / 1000.);
    }

    delete [] s_datcc;
    ccshape->setCurve (ccc);


    //LHcurv
    int *s_datch;
    s_datch = new int[70];
    int sizh;
    ImProcFunctions::strcurv_data (nextlh_str, s_datch, sizh);


    std::vector<double>   clh;

    for (int j = 0; j < sizh; j++) {
        clh.push_back ((double) (s_datch[j]) / 1000.);
    }

    delete [] s_datch;
    LHshape->setCurve (clh);


    //  usleep(10000);


    enableListener ();

    //update all sliders by this strange process!
    if (anbspot->getValue() == 0) {
        anbspot->setValue (1);

        if (options.rtSettings.locdelay) {
            if (anbspot->delay < 100) {
                anbspot->delay = 100;
            }
        }

        adjusterChanged (anbspot, 1);

    } else if (anbspot->getValue() == 1) {
        anbspot->setValue (0);

        if (options.rtSettings.locdelay) {
            if (anbspot->delay < 100) {
                anbspot->delay = 100;
            }
        }

        adjusterChanged (anbspot, 0);

    }


    //update all curves
    std::vector<double>   cretirab;
    cretirab = cTgainshaperab->getCurve ();

    if (cretirab.at (5) == 0.70) {
        cretirab.at (5) = 0.9;
        cTgainshaperab->setCurve (cretirab);

        curveChanged (cTgainshaperab);
    } else if (cretirab.at (5) == 0.90) {
        cretirab.at (5) = 0.7;
        cTgainshaperab->setCurve (cretirab);
        curveChanged (cTgainshaperab);

    }

    //

//   printf("G1 maj anbspot=%i  cretirab=%f\n", anbspot->getValue(), cretirab.at(5));


    //add events for each cases
    if (listener) { //for all sliders
        listener->panelChanged (Evlocallabanbspot, "");//anbspot->getTextValue());
    }

    if (listener) {//for curve
        listener->panelChanged (EvlocallabCTgainCurverab, M (""));
    }

    if (listener) {//for inverse color
        listener->panelChanged (Evlocallabinvers, M ("GENERAL_ENABLED"));
    }

    if (listener) {//for curvactiv
        listener->panelChanged (Evlocallabcurvactiv, M ("GENERAL_ENABLED"));
    }

    if (listener) {//for inverse blurr
        listener->panelChanged (Evlocallabinversrad, M ("GENERAL_ENABLED"));
    }

    if (listener) {//for quality method
        listener->panelChanged (EvlocallabqualityMethod, qualityMethod->get_active_text ());
    }

    if (listener) {//for quality method
        listener->panelChanged (EvlocallabqualitycurveMethod, qualitycurveMethod->get_active_text ());
    }

    if (listener) {//for inverse retinex
        listener->panelChanged (Evlocallabinversret, M ("GENERAL_ENABLED"));
    }

    if (listener) {//for inverse sharpen
        listener->panelChanged (Evlocallabinverssha, M ("GENERAL_ENABLED"));
    }

    if (listener) {//for Smethod : position of mouse cursor
        listener->panelChanged (EvlocallabSmet, Smethod->get_active_text ());
    }

    if (listener) {//for retinex method
        listener->panelChanged (EvlocallabretinexMethod, retinexMethod->get_active_text ());
    }

    if (listener) {//for curve reti
        listener->panelChanged (EvlocallabCTgainCurve, M (""));
    }

    if (listener) {//for curve LL
        listener->panelChanged (Evlocallabllshape, M (""));
    }

    if (listener) {//for curve LH
        listener->panelChanged (EvlocallabLHshape, M (""));
    }

    if (listener) {//for curve LH
        listener->panelChanged (Evlocallabccshape, M (""));
    }

    return false;
}

void Localrgb::localChanged  (int **datasp, std::string datastr, std::string ll_str, std::string lh_str, std::string cc_str, int sp, int maxdat)
{
    for (int i = 2; i < 61; i++) {
        nextdatasp[i] = datasp[i][sp];
    }

    nextstr = datastr;
    nextll_str = ll_str;
    nextlh_str = lh_str;
    nextcc_str = cc_str;

    nextlength = maxdat;
    g_idle_add (localChangedUI, this);
}

void Localrgb::localretChanged  (int **datasp, std::string datastr, std::string ll_str, std::string lh_str, std::string cc_str, int sp, int maxdat)
{
    nextlength = maxdat;
    nextstr2 = datastr;
    nextll_str2 = ll_str;
    nextlh_str2 = lh_str;
    nextcc_str2 = cc_str;

    g_idle_add (localretChangedUI, this);
}

*/
void Localrgb::read (const ProcParams* pp, const ParamsEdited* pedited)

{

    disableListener ();
    enableexposeConn.block (true);
    tcmodeconn.block (true);
    tcmode2conn.block (true);

    if (pedited) {
        degree->setEditedState (pedited->localrgb.degree ? Edited : UnEdited);
        locY->setEditedState (pedited->localrgb.locY ? Edited : UnEdited);
        locX->setEditedState (pedited->localrgb.locX ? Edited : UnEdited);
        locYT->setEditedState (pedited->localrgb.locYT ? Edited : UnEdited);
        locXL->setEditedState (pedited->localrgb.locXL ? Edited : UnEdited);
        centerX->setEditedState (pedited->localrgb.centerX ? Edited : UnEdited);
        centerY->setEditedState (pedited->localrgb.centerY ? Edited : UnEdited);
        circrad->setEditedState (pedited->localrgb.circrad ? Edited : UnEdited);
        thres->setEditedState (pedited->localrgb.thres ? Edited : UnEdited);
        proxi->setEditedState (pedited->localrgb.proxi ? Edited : UnEdited);
        lightness->setEditedState (pedited->localrgb.lightness ? Edited : UnEdited);
        contrast->setEditedState (pedited->localrgb.contrast ? Edited : UnEdited);
        chroma->setEditedState (pedited->localrgb.chroma ? Edited : UnEdited);
        sensi->setEditedState (pedited->localrgb.sensi ? Edited : UnEdited);
        nbspot->setEditedState (pedited->localrgb.nbspot ? Edited : UnEdited);
        anbspot->setEditedState (pedited->localrgb.anbspot ? Edited : UnEdited);
        retrab->setEditedState (pedited->localrgb.retrab ? Edited : UnEdited);
        hueref->setEditedState (pedited->localrgb.hueref ? Edited : UnEdited);
        chromaref->setEditedState (pedited->localrgb.chromaref ? Edited : UnEdited);
        lumaref->setEditedState (pedited->localrgb.lumaref ? Edited : UnEdited);
        transit->setEditedState (pedited->localrgb.transit ? Edited : UnEdited);
        expexpose->set_inconsistent   (!pedited->localrgb.expexpose);
        expcomp->setEditedState (pedited->localrgb.expcomp ? Edited : UnEdited);
        hlcompr->setEditedState (pedited->localrgb.hlcompr ? Edited : UnEdited);
        hlcomprthresh->setEditedState (pedited->localrgb.hlcomprthresh ? Edited : UnEdited);
        black->setEditedState (pedited->localrgb.black ? Edited : UnEdited);
        shcompr->setEditedState (pedited->localrgb.shcompr ? Edited : UnEdited);

        shape->setUnChanged (!pedited->localrgb.curve);
        shape2->setUnChanged (!pedited->localrgb.curve2);
        expwb->set_inconsistent   (!pedited->localrgb.expwb);
        temp->setEditedState (pedited->localrgb.temp ? Edited : UnEdited);
        green->setEditedState (pedited->localrgb.green ? Edited : UnEdited);
        equal->setEditedState (pedited->localrgb.equal ? Edited : UnEdited);

        if (!pedited->localrgb.curveMode) {
            toneCurveMode->set_active (6);
        }

        if (!pedited->localrgb.curveMode2) {
            toneCurveMode2->set_active (6);
        }

        if (!pedited->localrgb.Smethod) {
            Smethod->set_active_text (M ("GENERAL_UNCHANGED"));
        }

    }

    setEnabled (pp->localrgb.enabled);

    Smethodconn.block (true);
    qualityMethodConn.block (true);

    degree->setValue (pp->localrgb.degree);
    locY->setValue (pp->localrgb.locY);
    locX->setValue (pp->localrgb.locX);
    locYT->setValue (pp->localrgb.locYT);
    locXL->setValue (pp->localrgb.locXL);
    centerX->setValue (pp->localrgb.centerX);
    centerY->setValue (pp->localrgb.centerY);
    circrad->setValue (pp->localrgb.circrad);
    thres->setValue (pp->localrgb.thres);
    proxi->setValue (pp->localrgb.proxi);
    lightness->setValue (pp->localrgb.lightness);
    contrast->setValue (pp->localrgb.contrast);
    chroma->setValue (pp->localrgb.chroma);
    transit->setValue (pp->localrgb.transit);
    nbspot->setValue (pp->localrgb.nbspot);
    anbspot->setValue (pp->localrgb.anbspot);
    retrab->setValue (pp->localrgb.retrab);
    hueref->setValue (pp->localrgb.hueref);
    chromaref->setValue (pp->localrgb.chromaref);
    lumaref->setValue (pp->localrgb.lumaref);
    sensi->setValue (pp->localrgb.sensi);
    expcomp->setValue (pp->localrgb.expcomp);
    hlcompr->setValue (pp->localrgb.hlcompr);
    hlcomprthresh->setValue (pp->localrgb.hlcomprthresh);
    black->setValue (pp->localrgb.black);
    shcompr->setValue (pp->localrgb.shcompr);
    expexpose->setEnabled (pp->localrgb.expexpose);
    shape->setCurve (pp->localrgb.curve);
    shape2->setCurve (pp->localrgb.curve2);
    toneCurveMode->set_active (pp->localrgb.curveMode);
    toneCurveMode2->set_active (pp->localrgb.curveMode2);
    expwb->setEnabled (pp->localrgb.expwb);
    temp->setValue (pp->localrgb.temp);
    green->setValue (pp->localrgb.green);
    equal->setValue (pp->localrgb.equal);

    updateGeometry (pp->localrgb.centerX, pp->localrgb.centerY, pp->localrgb.circrad, pp->localrgb.locY, pp->localrgb.degree,  pp->localrgb.locX, pp->localrgb.locYT, pp->localrgb.locXL);

    if (pp->localrgb.Smethod == "IND") {
        Smethod->set_active (0);
    } else if (pp->localrgb.Smethod == "SYM") {
        Smethod->set_active (1);
    } else if (pp->localrgb.Smethod == "INDSL") {
        Smethod->set_active (2);
    } else if (pp->localrgb.Smethod == "SYMSL") {
        Smethod->set_active (3);
    }

    SmethodChanged();
    Smethodconn.block (false);

    if (pp->localrgb.qualityMethod == "std") {
        qualityMethod->set_active (0);
    } else if (pp->localrgb.qualityMethod == "enh") {
        qualityMethod->set_active (1);
    } else if (pp->localrgb.qualityMethod == "enhden") {
        qualityMethod->set_active (2);
    }

    qualityMethodChanged ();

    qualityMethodConn.block (false);

    anbspot->hide();
    hueref->hide();
    chromaref->hide();
    lumaref->hide();
    retrab->hide();

    if (pp->localrgb.Smethod == "SYM" || pp->localrgb.Smethod == "SYMSL") {
        locXL->setValue (locX->getValue());
        locYT->setValue (locY->getValue());
    } else if (pp->localrgb.Smethod == "LOC") {
        locXL->setValue (locX->getValue());
        locYT->setValue (locX->getValue());
        locY->setValue (locX->getValue());
    } else if (pp->localrgb.Smethod == "INDSL" || pp->localrgb.Smethod == "IND") {
        locX->setValue (pp->localrgb.locX);
        locY->setValue (pp->localrgb.locY);
        locXL->setValue (pp->localrgb.locXL);
        locYT->setValue (pp->localrgb.locYT);

    }

    tcmode2conn.block (false);
    tcmodeconn.block (false);

    enableexposeConn.block (false);
    enableListener ();

}

void Localrgb::updateGeometry (const int centerX_, const int centerY_, const int circrad_, const int locY_, const double degree_, const int locX_, const int locYT_, const int locXL_, const int fullWidth, const int fullHeight)
{
    EditDataProvider* dataProvider = getEditProvider();


    if (!dataProvider) {
        return;
    }

    int imW = 0;
    int imH = 0;

    if (fullWidth != -1 && fullHeight != -1) {
        imW = fullWidth;
        imH = fullHeight;
    } else {
        dataProvider->getImageSize (imW, imH);

        if (!imW || !imH) {
            return;
        }
    }

    PolarCoord polCoord1, polCoord2, polCoord0;
    // dataProvider->getImageSize(imW, imH);
    double decayY = (locY_) * double (imH) / 2000.;
    double decayYT = (locYT_) * double (imH) / 2000.;
    double decayX = (locX_) * (double (imW)) / 2000.;
    double decayXL = (locXL_) * (double (imW)) / 2000.;
    rtengine::Coord origin (imW / 2 + centerX_ * imW / 2000.f, imH / 2 + centerY_ * imH / 2000.f);
//   printf("deX=%f dexL=%f deY=%f deyT=%f locX=%i locY=%i\n", decayX, decayXL, decayY, decayYT, locX_, locY_);

    if (Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
        decayYT = decayY;
        decayXL = decayX;
    }

//    Line *currLine;
//    Circle *currCircle;
    //  Arcellipse *currArcellipse;
//    Beziers *currBeziers;
    double decay;
    /*
    const auto updateLine = [&] (Geometry * geometry, const float radius, const float begin, const float end) {
        const auto line = static_cast<Line*> (geometry);
        line->begin = PolarCoord (radius, -degree_ + begin);
        line->begin += origin;
        line->end = PolarCoord (radius, -degree_ + end);
        line->end += origin;
    };
    */
    const auto updateLineWithDecay = [&] (Geometry * geometry, const float radius, const float decal, const float offSetAngle) {
        const auto line = static_cast<Line*> (geometry); //180
        line->begin = PolarCoord (radius, -degree_ + decal) + PolarCoord (decay, -degree_ + offSetAngle);
        line->begin += origin;//0
        line->end = PolarCoord (radius, -degree_ + (decal - 180)) + PolarCoord (decay, -degree_ + offSetAngle);
        line->end += origin;
    };

    const auto updateCircle = [&] (Geometry * geometry) {
        const auto circle = static_cast<Circle*> (geometry);
        circle->center = origin;
        circle->radius = circrad_;
    };

    const auto updateBeziers = [&] (Geometry * geometry, const double dX_, const double dI_, const double dY_,  const float begi, const float inte, const float en) {
        const auto beziers = static_cast<Beziers*> (geometry);
        beziers->begin = PolarCoord (dX_, begi);
        beziers->begin += origin;//0
        beziers->inter = PolarCoord (dI_, inte);
        beziers->inter += origin;//0
        beziers->end = PolarCoord (dY_,  en);
        beziers->end += origin;
        //  printf("dX=%f dI=%f dY=%f begx=%i begy=%i intx=%i inty=%i endx=%i endy=%i\n", dX_, dI_, dY_, beziers->begin.x, beziers->begin.y, beziers->inter.x, beziers->inter.y, beziers->end.x, beziers->end.y);
    };

    /*
        const auto updateArcellipse = [&] (Geometry * geometry, const double dX_, const double dY_, const float kbegang, const float kendang) {
            const auto arcellipse = static_cast<Arcellipse*> (geometry);
            arcellipse->center = origin;
            arcellipse->radius = dY_;
            arcellipse->radius2 = dX_;
            arcellipse->translax = (double) imW /2.; //dX_ - dY_;
            arcellipse->translay = (double) imH /2.;
            arcellipse->scalx = dX_ / dY_; // double(locX_) / double (locY_); //arcellipse->radius2 / arcellipse->radius ; // dX_ / dY_;
            arcellipse->scaly = 1.; //dX_ / dY_; //locY_/locX_;
            arcellipse->begang = kbegang * M_PI;
            arcellipse->endang = kendang * M_PI;


        };
    */
    double dimline = 100.;

    if (options.showdelimspot) {
        dimline = 500.;
    }


    decay = decayX;
    updateLineWithDecay (visibleGeometry.at (0), dimline, 90., 0.);
    updateLineWithDecay (mouseOverGeometry.at (0), dimline, 90., 0.);

    decay = decayXL;

    updateLineWithDecay (visibleGeometry.at (1), dimline, 90., 180.);
    updateLineWithDecay (mouseOverGeometry.at (1), dimline, 90., 180.);

    decay = decayYT;
    updateLineWithDecay (visibleGeometry.at (2), dimline, 180., 270.);
    updateLineWithDecay (mouseOverGeometry.at (2), dimline, 180., 270.);

    decay = decayY;

    updateLineWithDecay (visibleGeometry.at (3), dimline, 180, 90.);
    updateLineWithDecay (mouseOverGeometry.at (3), dimline, 180., 90.);


    updateCircle (visibleGeometry.at (4));
    updateCircle (mouseOverGeometry.at (4));

    if (options.showdelimspot) {
        //this decayww evaluate approximation of a point in the ellipse for an angle alpha
        //this decayww evaluate approximation of a point in the ellipse for an angle alpha
        double decay5 = 1.003819 * ((decayX * decayY) / sqrt (0.00765 * SQR (decayX) + SQR (decayY))); //0.07179 = SQR(sin(15)/cos(15))  1.0038 = 1 / cos(5)
        double decay15 = 1.03527 * ((decayX * decayY) / sqrt (0.07179 * SQR (decayX) + SQR (decayY))); //0.07179 = SQR(sin(15)/cos(15))  1.03527 = 1 / cos(15)
        double decay30 = 1.15473 * ((decayX * decayY) / sqrt (0.33335 * SQR (decayX) + SQR (decayY)));
        double decay60 = 2. * ((decayX * decayY) / sqrt (3.0 * SQR (decayX) + SQR (decayY)));
        double decay75 = 3.86398 * ((decayX * decayY) / sqrt (13.929 * SQR (decayX) + SQR (decayY)));
        double decay85 = 11.473 * ((decayX * decayY) / sqrt (130.64 * SQR (decayX) + SQR (decayY)));

        double decay5L = 1.003819 * ((decayXL * decayY) / sqrt (0.00765 * SQR (decayXL) + SQR (decayY))); //0.07179 = SQR(sin(15)/cos(15))  1.0038 = 1 / cos(5)
        double decay15L = 1.03527 * ((decayXL * decayY) / sqrt (0.07179 * SQR (decayXL) + SQR (decayY)));
        double decay30L = 1.15473 * ((decayXL * decayY) / sqrt (0.33335 * SQR (decayXL) + SQR (decayY)));
        double decay60L = 2. * ((decayXL * decayY) / sqrt (3.0 * SQR (decayXL) + SQR (decayY)));
        double decay75L = 3.86398 * ((decayXL * decayY) / sqrt (13.929 * SQR (decayXL) + SQR (decayY)));
        double decay85L = 11.473 * ((decayXL * decayY) / sqrt (130.64 * SQR (decayXL) + SQR (decayY)));

        double decay5LT = 1.003819 * ((decayXL * decayYT) / sqrt (0.00765 * SQR (decayXL) + SQR (decayYT))); //0.07179 = SQR(sin(15)/cos(15))  1.0038 = 1 / cos(5)
        double decay15LT = 1.03527 * ((decayXL * decayYT) / sqrt (0.07179 * SQR (decayXL) + SQR (decayYT)));
        double decay30LT = 1.15473 * ((decayXL * decayYT) / sqrt (0.33335 * SQR (decayXL) + SQR (decayYT)));
        double decay60LT = 2. * ((decayXL * decayYT) / sqrt (3.0 * SQR (decayXL) + SQR (decayYT)));
        double decay75LT = 3.86398 * ((decayXL * decayYT) / sqrt (13.929 * SQR (decayXL) + SQR (decayYT)));
        double decay85LT = 11.473 * ((decayXL * decayYT) / sqrt (130.64 * SQR (decayXL) + SQR (decayYT)));

        double decay5T = 1.003819 * ((decayX * decayYT) / sqrt (0.00765 * SQR (decayX) + SQR (decayYT))); //0.07179 = SQR(sin(15)/cos(15))  1.0038 = 1 / cos(5)
        double decay15T = 1.03527 * ((decayX * decayYT) / sqrt (0.07179 * SQR (decayX) + SQR (decayYT)));
        double decay30T = 1.15473 * ((decayX * decayYT) / sqrt (0.33335 * SQR (decayX) + SQR (decayYT)));
        double decay60T = 2. * ((decayX * decayYT) / sqrt (3.0 * SQR (decayX) + SQR (decayYT)));
        double decay75T = 3.86398 * ((decayX * decayYT) / sqrt (13.929 * SQR (decayX) + SQR (decayYT)));
        double decay85T = 11.473 * ((decayX * decayYT) / sqrt (130.64 * SQR (decayX) + SQR (decayYT)));

        double decay45 = (1.414 * decayX * decayY) / sqrt (SQR (decayX) + SQR (decayY));
        double decay45L = (1.414 * decayXL * decayY) / sqrt (SQR (decayXL) + SQR (decayY));
        double decay45LT = (1.414 * decayXL * decayYT) / sqrt (SQR (decayXL) + SQR (decayYT));
        double decay45T = (1.414 * decayX * decayYT) / sqrt (SQR (decayX) + SQR (decayYT));

        //printf("decayX=%f decayY=%f decay10=%f decay45=%f oriX=%i origY=%i\n", decayX, decayY, decay10, decay45, origin.x, origin.y);
        updateBeziers (visibleGeometry.at (5), decayX, decay5  , decay15, 0., 5., 15.);
        updateBeziers (mouseOverGeometry.at (5), decayX, decay5 , decay15, 0., 5., 15.);

        updateBeziers (visibleGeometry.at (6), decay15, decay30 , decay45, 15., 30., 45.);
        updateBeziers (mouseOverGeometry.at (6), decay15, decay30 , decay45, 15., 30., 45.);

        updateBeziers (visibleGeometry.at (7), decay45, decay60 , decay75, 45., 60., 75.);
        updateBeziers (mouseOverGeometry.at (7), decay45, decay60 , decay75, 45., 60., 75.);

        updateBeziers (visibleGeometry.at (8), decay75, decay85 , decayY, 75., 85., 90.);
        updateBeziers (mouseOverGeometry.at (8), decay75, decay85 , decayY, 75., 85., 90.);

        updateBeziers (visibleGeometry.at (9), decayY, decay85L  , decay75L, 90., 95., 105.);
        updateBeziers (mouseOverGeometry.at (9), decayY, decay85L , decay75L, 90., 95., 105.);

        updateBeziers (visibleGeometry.at (10), decay75L, decay60L  , decay45L, 105., 120., 135.);
        updateBeziers (mouseOverGeometry.at (10), decay75L, decay60L , decay45L, 105., 120., 135.);

        updateBeziers (visibleGeometry.at (11), decay45L, decay30L  , decay15L, 135., 150., 165.);
        updateBeziers (mouseOverGeometry.at (11), decay45L, decay30L , decay15L, 135., 150., 165.);

        updateBeziers (visibleGeometry.at (12), decay15L, decay5L  , decayXL, 165., 175., 180.);
        updateBeziers (mouseOverGeometry.at (12), decay15L, decay5L , decayXL, 165., 175., 180.);


        updateBeziers (visibleGeometry.at (13), decayXL, decay5LT  , decay15LT, 180., 185., 195.);
        updateBeziers (mouseOverGeometry.at (13), decayXL, decay5LT , decay15LT, 180., 185., 195.);

        updateBeziers (visibleGeometry.at (14), decay15LT, decay30LT  , decay45LT, 195., 210., 225.);
        updateBeziers (mouseOverGeometry.at (14), decay15LT, decay30LT , decay45LT, 195., 210., 225.);

        updateBeziers (visibleGeometry.at (15), decay45LT, decay60LT  , decay75LT, 225., 240., 255.);
        updateBeziers (mouseOverGeometry.at (15), decay45LT, decay60LT , decay75LT, 225., 240., 255.);

        updateBeziers (visibleGeometry.at (16), decay75LT, decay85LT  , decayYT, 255., 265., 270.);
        updateBeziers (mouseOverGeometry.at (16), decay75LT, decay85LT , decayYT, 255., 265., 270.);

        updateBeziers (visibleGeometry.at (17), decayYT, decay85T  , decay75T, 270., 275., 285.);
        updateBeziers (mouseOverGeometry.at (17), decayYT, decay85T , decay75T, 270., 275., 285.);

        updateBeziers (visibleGeometry.at (18), decay75T, decay60T  , decay45T, 285., 300., 315.);
        updateBeziers (mouseOverGeometry.at (18), decay75T, decay60T , decay45T, 285., 300., 315.);

        updateBeziers (visibleGeometry.at (19), decay45T, decay30T  , decay15T, 315., 330., 345.);
        updateBeziers (mouseOverGeometry.at (19), decay45T, decay30T , decay15T, 315., 330., 345.);

        updateBeziers (visibleGeometry.at (20), decay15T, decay5T  , decayX, 345., 355., 360.);
        updateBeziers (mouseOverGeometry.at (20), decay15T, decay5T , decayX, 345., 355., 360.);

    }

    //  updateArcellipse (visibleGeometry.at (5), decayX, decayY, 0., 0.5);
    //  updateArcellipse (mouseOverGeometry.at (5), decayX, decayY, 0., 0.5);

}

void Localrgb::write (ProcParams* pp, ParamsEdited* pedited)
{
    pp->localrgb.degree = degree->getValue ();
    pp->localrgb.locY = locY->getIntValue ();
    pp->localrgb.locX = locX->getValue ();
    pp->localrgb.locYT = locYT->getIntValue ();
    pp->localrgb.locXL = locXL->getValue ();
    pp->localrgb.centerX = centerX->getIntValue ();
    pp->localrgb.centerY = centerY->getIntValue ();
    pp->localrgb.circrad = circrad->getIntValue ();
    pp->localrgb.proxi = proxi->getIntValue ();
    pp->localrgb.thres = thres->getIntValue ();
    pp->localrgb.lightness = lightness->getIntValue ();
    pp->localrgb.contrast = contrast->getIntValue ();
    pp->localrgb.chroma = chroma->getIntValue ();
    pp->localrgb.sensi = sensi->getIntValue ();
    pp->localrgb.transit = transit->getIntValue ();
    pp->localrgb.nbspot = nbspot->getIntValue ();
    pp->localrgb.anbspot = anbspot->getIntValue ();
    pp->localrgb.retrab = retrab->getIntValue ();
    pp->localrgb.hueref = hueref->getValue ();
    pp->localrgb.chromaref = chromaref->getValue ();
    pp->localrgb.lumaref = lumaref->getValue ();
    pp->localrgb.expcomp = expcomp->getValue ();
    pp->localrgb.black = (int)black->getValue ();
    pp->localrgb.hlcompr = (int)hlcompr->getValue ();
    pp->localrgb.hlcomprthresh = (int)hlcomprthresh->getValue ();
    pp->localrgb.shcompr = (int)shcompr->getValue ();
    pp->localrgb.enabled = getEnabled();
    pp->localrgb.curve = shape->getCurve ();
    pp->localrgb.curve2 = shape2->getCurve ();
    pp->localrgb.temp = temp->getValue ();
    pp->localrgb.green = green->getValue ();
    pp->localrgb.equal = equal->getValue ();

    int tcMode = toneCurveMode->get_active_row_number();

    if      (tcMode == 0) {
        pp->localrgb.curveMode = LocalrgbParams::TC_MODE_STD;
    } else if (tcMode == 1) {
        pp->localrgb.curveMode = LocalrgbParams::TC_MODE_WEIGHTEDSTD;
    } else if (tcMode == 2) {
        pp->localrgb.curveMode = LocalrgbParams::TC_MODE_FILMLIKE;
    } else if (tcMode == 3) {
        pp->localrgb.curveMode = LocalrgbParams::TC_MODE_SATANDVALBLENDING;
    } else if (tcMode == 4) {
        pp->localrgb.curveMode = LocalrgbParams::TC_MODE_LUMINANCE;
    } else if (tcMode == 5) {
        pp->localrgb.curveMode = LocalrgbParams::TC_MODE_PERCEPTUAL;
    }

    tcMode = toneCurveMode2->get_active_row_number();

    if      (tcMode == 0) {
        pp->localrgb.curveMode2 = LocalrgbParams::TC_MODE_STD;
    } else if (tcMode == 1) {
        pp->localrgb.curveMode2 = LocalrgbParams::TC_MODE_WEIGHTEDSTD;
    } else if (tcMode == 2) {
        pp->localrgb.curveMode2 = LocalrgbParams::TC_MODE_FILMLIKE;
    } else if (tcMode == 3) {
        pp->localrgb.curveMode2 = LocalrgbParams::TC_MODE_SATANDVALBLENDING;
    } else if (tcMode == 4) {
        pp->localrgb.curveMode2 = LocalrgbParams::TC_MODE_LUMINANCE;
    } else if (tcMode == 5) {
        pp->localrgb.curveMode2 = LocalrgbParams::TC_MODE_PERCEPTUAL;
    }

    pp->localrgb.expexpose      = expexpose->getEnabled();
    pp->localrgb.expwb      = expwb->getEnabled();

    if (pedited) {
        pedited->localrgb.degree = degree->getEditedState ();
        pedited->localrgb.Smethod  = Smethod->get_active_text() != M ("GENERAL_UNCHANGED");
        pedited->localrgb.qualityMethod    = qualityMethod->get_active_text() != M ("GENERAL_UNCHANGED");
        pedited->localrgb.locY = locY->getEditedState ();
        pedited->localrgb.locX = locX->getEditedState ();
        pedited->localrgb.locYT = locYT->getEditedState ();
        pedited->localrgb.locXL = locXL->getEditedState ();
        pedited->localrgb.centerX = centerX->getEditedState ();
        pedited->localrgb.centerY = centerY->getEditedState ();
        pedited->localrgb.circrad = circrad->getEditedState ();
        pedited->localrgb.proxi = proxi->getEditedState ();
        pedited->localrgb.thres = thres->getEditedState ();
        pedited->localrgb.lightness = lightness->getEditedState ();
        pedited->localrgb.contrast = contrast->getEditedState ();
        pedited->localrgb.chroma = chroma->getEditedState ();
        pedited->localrgb.sensi = sensi->getEditedState ();
        pedited->localrgb.transit = transit->getEditedState ();
        pedited->localrgb.nbspot = nbspot->getEditedState ();
        pedited->localrgb.anbspot = anbspot->getEditedState ();
        pedited->localrgb.retrab = retrab->getEditedState ();
        pedited->localrgb.hueref = hueref->getEditedState ();
        pedited->localrgb.chromaref = chromaref->getEditedState ();
        pedited->localrgb.lumaref = lumaref->getEditedState ();
        pedited->localrgb.expcomp    = expcomp->getEditedState ();
        pedited->localrgb.black      = black->getEditedState ();
        pedited->localrgb.hlcompr    = hlcompr->getEditedState ();
        pedited->localrgb.hlcomprthresh = hlcomprthresh->getEditedState ();
        pedited->localrgb.shcompr    = shcompr->getEditedState ();
        pedited->localrgb.curve      = !shape->isUnChanged ();
        pedited->localrgb.curve2     = !shape2->isUnChanged ();
        pedited->localrgb.curveMode  = toneCurveMode->get_active_row_number() != 6;
        pedited->localrgb.curveMode2 = toneCurveMode2->get_active_row_number() != 6;

        pedited->localrgb.enabled = !get_inconsistent();

        pedited->localrgb.expexpose     = !expexpose->get_inconsistent();
        pedited->localrgb.expwb     = !expwb->get_inconsistent();
        pedited->localrgb.temp    = temp->getEditedState ();
        pedited->localrgb.green    = green->getEditedState ();
        pedited->localrgb.equal    = equal->getEditedState ();

    }


    if (qualityMethod->get_active_row_number() == 0) {
        pp->localrgb.qualityMethod = "std";
    } else if (qualityMethod->get_active_row_number() == 1) {
        pp->localrgb.qualityMethod = "enh";
    } else if (qualityMethod->get_active_row_number() == 2) {
        pp->localrgb.qualityMethod = "enhden";
    }


    if (Smethod->get_active_row_number() == 0) {
        pp->localrgb.Smethod = "IND";
    } else if (Smethod->get_active_row_number() == 1) {
        pp->localrgb.Smethod = "SYM";
    } else if (Smethod->get_active_row_number() == 2) {
        pp->localrgb.Smethod = "INDSL";
    } else if (Smethod->get_active_row_number() == 3) {
        pp->localrgb.Smethod = "SYMSL";
    }

    if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
        //   if(Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 1) {
        pp->localrgb.locX = locX->getValue();
        pp->localrgb.locY = locY->getValue();

        pp->localrgb.locXL = pp->localrgb.locX;
        pp->localrgb.locYT = pp->localrgb.locY;
    } else {
        pp->localrgb.locXL = locXL->getValue();
        pp->localrgb.locX = locX->getValue();
        pp->localrgb.locY = locY->getValue();
        pp->localrgb.locYT = locYT->getValue();

    }
}

void Localrgb::qualityMethodChanged()
{
    if (!batchMode) {
    }

    if (listener) {
        listener->panelChanged (EvlocalrgbqualityMethod, qualityMethod->get_active_text ());
    }
}



void Localrgb::curveChanged (CurveEditor* ce)
{
    if (listener) {
        if (ce == shape) {
            listener->panelChanged (EvlocalrgbCurve1, M ("HISTORY_CUSTOMCURVE"));
        } else if (ce == shape2) {
            listener->panelChanged (EvlocalrgbCurve2, M ("HISTORY_CUSTOMCURVE"));
        }
    }

    /*
        if (listener && getEnabled()) {
            if (ce == cTgainshape) {
                listener->panelChanged (EvlocallabCTgainCurve, M ("HISTORY_CUSTOMCURVE"));//HISTORY_CUSTOMCURVE
                int strval = retrab->getValue();
                //update MIP
                retrab->setValue (strval + 1);
                adjusterChanged (retrab, strval + 1);
                usleep (10000); //to test
                retrab->setValue (strval);

                adjusterChanged (retrab, strval);
            }

            else if (ce == cTgainshaperab) {
                listener->panelChanged (EvlocallabCTgainCurverab, M (""));
            } else if (ce == LHshape) {
                listener->panelChanged (EvlocallabLHshape, M (""));
                int strval = retrab->getValue();
                //update MIP
                retrab->setValue (strval + 1);
                adjusterChanged (retrab, strval + 1);
                usleep (10000); //to test
                retrab->setValue (strval);

                adjusterChanged (retrab, strval);


            } else if (ce == llshape) {
                listener->panelChanged (Evlocallabllshape, M ("HISTORY_CUSTOMCURVE"));
                int strval = retrab->getValue();
                //update MIP
                retrab->setValue (strval + 1);
                adjusterChanged (retrab, strval + 1);
                usleep (10000); //to test
                retrab->setValue (strval);

                adjusterChanged (retrab, strval);

            } else if (ce == ccshape) {
                listener->panelChanged (Evlocallabccshape, M ("HISTORY_CUSTOMCURVE"));
                int strval = retrab->getValue();
                //update MIP
                retrab->setValue (strval + 1);
                adjusterChanged (retrab, strval + 1);
                usleep (10000); //to test
                retrab->setValue (strval);

                adjusterChanged (retrab, strval);

            }

        }
        */
}

void Localrgb::curveMode1Changed ()
{
    if (listener) {
        Glib::signal_idle().connect (sigc::mem_fun (*this, &Localrgb::curveMode1Changed_));
    }
}

bool Localrgb::curveMode1Changed_ ()
{
    if (listener) {
        listener->panelChanged (EvlocalrgbCurveMode1, toneCurveMode->get_active_text());
    }

    return false;
}

void Localrgb::curveMode2Changed ()
{
    //if (listener)  listener->panelChanged (EvToneCurveMode, toneCurveMode->get_active_text());
    if (listener) {
        Glib::signal_idle().connect (sigc::mem_fun (*this, &Localrgb::curveMode2Changed_));
    }
}

bool Localrgb::curveMode2Changed_ ()
{
    if (listener) {
        listener->panelChanged (EvlocalrgbCurveMode2, toneCurveMode2->get_active_text());
    }

    return false;
}


void Localrgb::SmethodChanged ()
{
    if (!batchMode) {
        if (Smethod->get_active_row_number() == 0) { //IND 0
            locX->hide();
            locXL->hide();
            locY->hide();
            locYT->hide();
            centerX->hide();
            centerY->hide();
        } else if (Smethod->get_active_row_number() == 1) {         // 1 SYM
            locX->hide();
            locXL->hide();
            locY->hide();
            locYT->hide();
            centerX->hide();
            centerY->hide();

        } else if (Smethod->get_active_row_number() == 2) {         //2 SYM
            locX->show();
            locXL->show();
            locY->show();
            locYT->show();
            centerX->show();
            centerY->show();

        } else if (Smethod->get_active_row_number() == 3) {         // 3 SYM
            locX->show();
            locXL->hide();
            locY->show();
            locYT->hide();
            centerX->show();
            centerY->show();

        }

    }

    if (listener && getEnabled()) {
        if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
            listener->panelChanged (EvlocalrgbSmet, Smethod->get_active_text ());
            locXL->setValue (locX->getValue());
            locYT->setValue (locY->getValue());
        }
        //   else if(Smethod->get_active_row_number()==2) {
        //          listener->panelChanged (EvlocallabSmet, Smethod->get_active_text ());
        //           locXL->setValue (locX->getValue());
        //           locYT->setValue (locX->getValue());
        //          locY->setValue (locX->getValue());
        //     }
        else

        {
            listener->panelChanged (EvlocalrgbSmet, Smethod->get_active_text ());

        }
    }
}


void Localrgb::setDefaults (const ProcParams * defParams, const ParamsEdited * pedited)
{
    degree->setDefault (defParams->localrgb.degree);
    locY->setDefault (defParams->localrgb.locY);
    locX->setDefault (defParams->localrgb.locX);
    locYT->setDefault (defParams->localrgb.locYT);
    locXL->setDefault (defParams->localrgb.locXL);
    centerX->setDefault (defParams->localrgb.centerX);
    centerY->setDefault (defParams->localrgb.centerY);
    circrad->setDefault (defParams->localrgb.circrad);
    thres->setDefault (defParams->localrgb.thres);
    proxi->setDefault (defParams->localrgb.proxi);
    lightness->setDefault (defParams->localrgb.lightness);
    contrast->setDefault (defParams->localrgb.contrast);
    chroma->setDefault (defParams->localrgb.chroma);

    expcomp->setDefault (defParams->localrgb.expcomp);
    black->setDefault (defParams->localrgb.black);
    hlcompr->setDefault (defParams->localrgb.hlcompr);
    hlcomprthresh->setDefault (defParams->localrgb.hlcomprthresh);
    shcompr->setDefault (defParams->localrgb.shcompr);


    sensi->setDefault (defParams->localrgb.sensi);
    transit->setDefault (defParams->localrgb.transit);
    nbspot->setDefault (defParams->localrgb.nbspot);
    anbspot->setDefault (defParams->localrgb.anbspot);
    retrab->setDefault (defParams->localrgb.retrab);
    hueref->setDefault (defParams->localrgb.hueref);
    chromaref->setDefault (defParams->localrgb.chromaref);
    lumaref->setDefault (defParams->localrgb.lumaref);
    temp->setDefault (defParams->localrgb.temp);
    green->setDefault (defParams->localrgb.green);
    equal->setDefault (defParams->localrgb.equal);


    if (pedited) {
        degree->setDefaultEditedState (pedited->localrgb.degree ? Edited : UnEdited);
        locY->setDefaultEditedState (pedited->localrgb.locY ? Edited : UnEdited);
        locX->setDefaultEditedState (pedited->localrgb.locX ? Edited : UnEdited);
        locYT->setDefaultEditedState (pedited->localrgb.locYT ? Edited : UnEdited);
        locXL->setDefaultEditedState (pedited->localrgb.locXL ? Edited : UnEdited);
        centerX->setDefaultEditedState (pedited->localrgb.centerX ? Edited : UnEdited);
        centerY->setDefaultEditedState (pedited->localrgb.centerY ? Edited : UnEdited);
        circrad->setDefaultEditedState (pedited->localrgb.circrad ? Edited : UnEdited);
        thres->setDefaultEditedState (pedited->localrgb.thres ? Edited : UnEdited);
        proxi->setDefaultEditedState (pedited->localrgb.proxi ? Edited : UnEdited);
        lightness->setDefaultEditedState (pedited->localrgb.lightness ? Edited : UnEdited);
        contrast->setDefaultEditedState (pedited->localrgb.contrast ? Edited : UnEdited);
        chroma->setDefaultEditedState (pedited->localrgb.chroma ? Edited : UnEdited);
        sensi->setDefaultEditedState (pedited->localrgb.sensi ? Edited : UnEdited);
        transit->setDefaultEditedState (pedited->localrgb.transit ? Edited : UnEdited);
        nbspot->setDefaultEditedState (pedited->localrgb.nbspot ? Edited : UnEdited);
        anbspot->setDefaultEditedState (pedited->localrgb.anbspot ? Edited : UnEdited);
        retrab->setDefaultEditedState (pedited->localrgb.retrab ? Edited : UnEdited);
        hueref->setDefaultEditedState (pedited->localrgb.hueref ? Edited : UnEdited);
        chromaref->setDefaultEditedState (pedited->localrgb.chromaref ? Edited : UnEdited);
        lumaref->setDefaultEditedState (pedited->localrgb.lumaref ? Edited : UnEdited);
        expcomp->setDefaultEditedState (pedited->localrgb.expcomp ? Edited : UnEdited);
        black->setDefaultEditedState (pedited->localrgb.black ? Edited : UnEdited);
        hlcompr->setDefaultEditedState (pedited->localrgb.hlcompr ? Edited : UnEdited);
        hlcomprthresh->setDefaultEditedState (pedited->localrgb.hlcomprthresh ? Edited : UnEdited);
        shcompr->setDefaultEditedState (pedited->localrgb.shcompr ? Edited : UnEdited);
        temp->setDefaultEditedState (pedited->localrgb.temp ? Edited : UnEdited);
        green->setDefaultEditedState (pedited->localrgb.green ? Edited : UnEdited);
        equal->setDefaultEditedState (pedited->localrgb.equal ? Edited : UnEdited);

    } else {
        degree->setDefaultEditedState (Irrelevant);
        locY->setDefaultEditedState (Irrelevant);
        locX->setDefaultEditedState (Irrelevant);
        locYT->setDefaultEditedState (Irrelevant);
        locXL->setDefaultEditedState (Irrelevant);
        centerX->setDefaultEditedState (Irrelevant);
        centerY->setDefaultEditedState (Irrelevant);
        circrad->setDefaultEditedState (Irrelevant);
        thres->setDefaultEditedState (Irrelevant);
        proxi->setDefaultEditedState (Irrelevant);
        lightness->setDefaultEditedState (Irrelevant);
        contrast->setDefaultEditedState (Irrelevant);
        chroma->setDefaultEditedState (Irrelevant);
        sensi->setDefaultEditedState (Irrelevant);
        transit->setDefaultEditedState (Irrelevant);
        nbspot->setDefaultEditedState (Irrelevant);
        anbspot->setDefaultEditedState (Irrelevant);
        retrab->setDefaultEditedState (Irrelevant);
        hueref->setDefaultEditedState (Irrelevant);
        chromaref->setDefaultEditedState (Irrelevant);
        lumaref->setDefaultEditedState (Irrelevant);
        expcomp->setDefaultEditedState (Irrelevant);
        black->setDefaultEditedState (Irrelevant);
        hlcompr->setDefaultEditedState (Irrelevant);
        hlcomprthresh->setDefaultEditedState (Irrelevant);
        shcompr->setDefaultEditedState (Irrelevant);
        temp->setDefaultEditedState (Irrelevant);
        green->setDefaultEditedState (Irrelevant);
        equal->setDefaultEditedState (Irrelevant);


    }

}

void Localrgb::adjusterChanged (Adjuster * a, double newval)
{

    updateGeometry (int (centerX->getValue()), int (centerY->getValue()), int (circrad->getValue()), (int)locY->getValue(), degree->getValue(), (int)locX->getValue(), (int)locYT->getValue(), (int)locXL->getValue());
    anbspot->hide();
    retrab->hide();
    hueref->hide();
    chromaref->hide();
    lumaref->hide();

    if (listener && getEnabled()) {
        if (a == degree) {
            listener->panelChanged (EvlocalrgbDegree, degree->getTextValue());
        } else if (a == locY) {
            if (Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 2) { // 0 2
                listener->panelChanged (EvlocalrgblocY, locY->getTextValue());
            } else if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
                listener->panelChanged (EvlocalrgblocY, locY->getTextValue());
                locYT->setValue (locY->getValue());
            }
        } else if (a == locX) {
            //listener->panelChanged (EvlocallablocX, locX->getTextValue());
            if (Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 2) {
                listener->panelChanged (EvlocalrgblocX, locX->getTextValue());
            } else if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
                listener->panelChanged (EvlocalrgblocX, locX->getTextValue());
                locXL->setValue (locX->getValue());
            }
        } else if (a == locYT) {
            if (Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 2) {
                listener->panelChanged (EvlocalrgblocYT, locYT->getTextValue());
            } else if (Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
                listener->panelChanged (EvlocalrgblocYT, locYT->getTextValue());
                locYT->setValue (locY->getValue());
            }
        } else if (a == locXL) {
            if (Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 2) {
                listener->panelChanged (EvlocalrgblocXL, locXL->getTextValue());
            } else if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
                listener->panelChanged (EvlocalrgblocXL, locXL->getTextValue());
                locXL->setValue (locX->getValue());
            }
        } else if (a == lightness) {
            listener->panelChanged (Evlocalrgblightness, lightness->getTextValue());
        } else if (a == contrast) {
            listener->panelChanged (Evlocalrgbcontrast, contrast->getTextValue());
        } else if (a == chroma) {
            listener->panelChanged (Evlocalrgbchroma, chroma->getTextValue());
        } else if (a == sensi) {
            listener->panelChanged (Evlocalrgbsensi, sensi->getTextValue());
        } else if (a == expcomp) {
            listener->panelChanged (Evlocalrgbexpcomp, expcomp->getTextValue());
        } else if (a == hlcompr) {
            listener->panelChanged (Evlocalrgbhlcompr, hlcompr->getTextValue());
        } else if (a == hlcomprthresh) {
            listener->panelChanged (Evlocalrgbhlcomprthresh, hlcomprthresh->getTextValue());
        } else if (a == black) {
            listener->panelChanged (Evlocalrgbblack, black->getTextValue());
        } else if (a == shcompr) {
            listener->panelChanged (Evlocalrgbshcompr, shcompr->getTextValue());
        } else if (a == transit) {
            listener->panelChanged (Evlocalrgbtransit, transit->getTextValue());
        } else if (a == nbspot) {
            listener->panelChanged (Evlocalrgbnbspot, nbspot->getTextValue());
        } else if (a == retrab) {
            listener->panelChanged (Evlocalrgbanbspot, "");//anbspot->getTextValue());
        } else if (a == anbspot) {
            listener->panelChanged (Evlocalrgbretrab, "");//anbspot->getTextValue());
        } else if (a == hueref) {
            listener->panelChanged (Evlocalrgbhueref, "");//anbspot->getTextValue());
        } else if (a == chromaref) {
            listener->panelChanged (Evlocalrgbchromaref, "");//anbspot->getTextValue());
        } else if (a == temp) {
            listener->panelChanged (Evlocalrgbtemp, temp->getTextValue());
        } else if (a == green) {
            listener->panelChanged (Evlocalrgbgreen, green->getTextValue());
        } else if (a == equal) {
            listener->panelChanged (Evlocalrgbequal, equal->getTextValue());
        } else if (a == lumaref) {
            listener->panelChanged (Evlocalrgblumaref, "");//anbspot->getTextValue());
        } else if (a == circrad) {
            listener->panelChanged (Evlocalrgbcircrad, circrad->getTextValue());
        } else if (a == thres) {
            listener->panelChanged (Evlocalrgbthres, thres->getTextValue());
        } else if (a == proxi) {
            listener->panelChanged (Evlocalrgbproxi, proxi->getTextValue());
        } else if (a == centerX || a == centerY) {
            listener->panelChanged (EvlocalrgbCenter, Glib::ustring::compose ("X=%1\nY=%2", centerX->getTextValue(), centerY->getTextValue()));
        } else {
            //
        }
    }

}

void Localrgb::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvlocalrgbEnabled, M ("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvlocalrgbEnabled, M ("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvlocalrgbEnabled, M ("GENERAL_DISABLED"));
        }
    }
}
void Localrgb::setAdjusterBehavior (bool degreeadd, bool locYadd, bool locXadd, bool locYTadd, bool locXLadd, bool centeradd, bool lightnessadd, bool contrastadd, bool chromaadd, bool sensiadd, bool transitadd, bool radiusadd,  bool strengthadd)
{
    /*
    degree->setAddMode (degreeadd);
    locY->setAddMode (locYadd);
    locX->setAddMode (locXadd);
    locYT->setAddMode (locYTadd);
    locXL->setAddMode (locXLadd);
    centerX->setAddMode (centeradd);
    centerY->setAddMode (centeradd);
    lightness->setAddMode (lightnessadd);
    contrast->setAddMode (contrastadd);
    chroma->setAddMode (chromaadd);
    sensi->setAddMode (sensiadd);
    transit->setAddMode (transitadd);
    radius->setAddMode (radiusadd);
    strength->setAddMode (strengthadd);
    */
}

void Localrgb::trimValues (rtengine::procparams::ProcParams * pp)
{

    degree->trimValue (pp->localrgb.degree);
    locY->trimValue (pp->localrgb.locY);
    locX->trimValue (pp->localrgb.locX);
    locYT->trimValue (pp->localrgb.locYT);
    locXL->trimValue (pp->localrgb.locXL);
    centerX->trimValue (pp->localrgb.centerX);
    centerY->trimValue (pp->localrgb.centerY);
    circrad->trimValue (pp->localrgb.circrad);
    thres->trimValue (pp->localrgb.thres);
    proxi->trimValue (pp->localrgb.proxi);
    lightness->trimValue (pp->localrgb.lightness);
    expcomp->trimValue (pp->localrgb.expcomp);
    hlcompr->trimValue (pp->localrgb.hlcompr);
    hlcomprthresh->trimValue (pp->localrgb.hlcomprthresh);
    black->trimValue (pp->localrgb.black);
    shcompr->trimValue (pp->localrgb.shcompr);

    contrast->trimValue (pp->localrgb.contrast);
    chroma->trimValue (pp->localrgb.chroma);
    sensi->trimValue (pp->localrgb.sensi);
    transit->trimValue (pp->localrgb.transit);
    nbspot->trimValue (pp->localrgb.nbspot);
    anbspot->trimValue (pp->localrgb.anbspot);
    retrab->trimValue (pp->localrgb.retrab);
    hueref->trimValue (pp->localrgb.hueref);
    chromaref->trimValue (pp->localrgb.chromaref);
    lumaref->trimValue (pp->localrgb.lumaref);
    temp->trimValue (pp->localrgb.temp);
    green->trimValue (pp->localrgb.green);
    equal->trimValue (pp->localrgb.equal);

}

void Localrgb::setBatchMode (bool batchMode)
{
    removeIfThere (this, edit, false);
    ToolPanel::setBatchMode (batchMode);
    degree->showEditedCB ();
    locY->showEditedCB ();
    locX->showEditedCB ();
    locYT->showEditedCB ();
    locXL->showEditedCB ();
    centerX->showEditedCB ();
    centerY->showEditedCB ();
    circrad->showEditedCB ();
    thres->showEditedCB ();
    proxi->showEditedCB ();
    lightness->showEditedCB ();
    contrast->showEditedCB ();
    chroma->showEditedCB ();
    expcomp->showEditedCB ();
    black->showEditedCB ();
    hlcompr->showEditedCB ();
    hlcomprthresh->showEditedCB ();
    shcompr->showEditedCB ();

    sensi->showEditedCB ();
    transit->showEditedCB ();
    Smethod->append (M ("GENERAL_UNCHANGED"));
    nbspot->showEditedCB ();
    anbspot->showEditedCB ();
    retrab->showEditedCB ();
    hueref->showEditedCB ();
    chromaref->showEditedCB ();
    lumaref->showEditedCB ();
    temp->showEditedCB ();
    green->showEditedCB ();
    equal->showEditedCB ();

}

void Localrgb::setEditProvider (EditDataProvider * provider)
{
    EditSubscriber::setEditProvider (provider);
    shape->setEditProvider (provider);
    shape2->setEditProvider (provider);

}

void Localrgb::editToggled ()
{
    if (edit->get_active()) {
        subscribe();
    } else {
        unsubscribe();
    }
}

void Localrgb::colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller)
{

    float R, G, B;

    if (elemType == ColorCaller::CCET_VERTICAL_BAR) {
        valY = 0.5;
    }

    if (callerId == 1) {         // ch - main curve

        Color::hsv2rgb01 (float (valX), float (valY), 0.5f, R, G, B);
    } else if (callerId == 2) {  // cc - bottom bar

        float value = (1.f - 0.7f) * float (valX) + 0.7f;
        // whole hue range
        // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
        Color::hsv2rgb01 (float (valY), float (valX), value, R, G, B);
    } else if (callerId == 3) {  // lc - bottom bar

        float value = (1.f - 0.7f) * float (valX) + 0.7f;
        // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
        Color::hsv2rgb01 (float (valY), float (valX), value, R, G, B);
    } else if (callerId == 4) {  // LH - bottom bar
        Color::hsv2rgb01 (float (valX), 0.5f, float (valY), R, G, B);
    } else if (callerId == 5) {  // HH - bottom bar
        float h = float ((valY - 0.5) * 0.3 + valX);

        if (h > 1.0f) {
            h -= 1.0f;
        } else if (h < 0.0f) {
            h += 1.0f;
        }

        Color::hsv2rgb01 (h, 0.5f, 0.5f, R, G, B);
    }

    caller->ccRed = double (R);
    caller->ccGreen = double (G);
    caller->ccBlue = double (B);
}



CursorShape Localrgb::getCursor (int objectID)
{
    switch (objectID) {
        case (2): {
            int angle = degree->getIntValue();

            if (angle < -135 || (angle >= -45 && angle <= 45) || angle > 135) {
                return CSMove1DV;
            }

            return CSMove1DH;
        }

        case (3): {
            int angle = degree->getIntValue();

            if (angle < -135 || (angle >= -45 && angle <= 45) || angle > 135) {
                return CSMove1DV;
            }

            return CSMove1DH;
        }

        case (0): {
            int angle = degree->getIntValue();

            if (angle < -135 || (angle >= -45 && angle <= 45) || angle > 135) {
                return CSMove1DH;
            }

            return CSMove1DV;
        }

        case (1): {
            int angle = degree->getIntValue();

            if (angle < -135 || (angle >= -45 && angle <= 45) || angle > 135) {
                return CSMove1DH;
            }

            return CSMove1DV;
        }

        case (4):
            return CSMove2D;

        default:
            return CSOpenHand;
    }
}

bool Localrgb::mouseOver (int modifierKey)
{
    EditDataProvider* editProvider = getEditProvider();

    if (editProvider && editProvider->object != lastObject) {
        if (lastObject > -1) {
            if (lastObject == 2 || lastObject == 3) {
                EditSubscriber::visibleGeometry.at (2)->state = Geometry::NORMAL;
                EditSubscriber::visibleGeometry.at (3)->state = Geometry::NORMAL;

            } else if (lastObject == 0 || lastObject == 1) {
                EditSubscriber::visibleGeometry.at (0)->state = Geometry::NORMAL;
                EditSubscriber::visibleGeometry.at (1)->state = Geometry::NORMAL;

            }

            else {
                EditSubscriber::visibleGeometry.at (4)->state = Geometry::NORMAL;
//               EditSubscriber::visibleGeometry.at (lastObject)->state = Geometry::NORMAL;
            }
        }

        if (editProvider->object > -1) {
            if (editProvider->object == 2 || editProvider->object == 3) {
                EditSubscriber::visibleGeometry.at (2)->state = Geometry::PRELIGHT;
                EditSubscriber::visibleGeometry.at (3)->state = Geometry::PRELIGHT;

            } else if (editProvider->object == 0 || editProvider->object == 1) {
                EditSubscriber::visibleGeometry.at (0)->state = Geometry::PRELIGHT;
                EditSubscriber::visibleGeometry.at (1)->state = Geometry::PRELIGHT;

            }

            else {
                EditSubscriber::visibleGeometry.at (4)->state = Geometry::PRELIGHT;
                //              EditSubscriber::visibleGeometry.at (editProvider->object)->state = Geometry::PRELIGHT;
            }
        }

        lastObject = editProvider->object;
        return true;
    }

    return false;
}

bool Localrgb::button1Pressed (int modifierKey)
{
    if (lastObject < 0) {
        return false;
    }

    EditDataProvider *provider = getEditProvider();

    if (! (modifierKey & GDK_CONTROL_MASK)) {
        // button press is valid (no modifier key)
        PolarCoord pCoord;
        //  EditDataProvider *provider = getEditProvider();
        int imW, imH;
        provider->getImageSize (imW, imH);
        double halfSizeW = imW / 2.;
        double halfSizeH = imH / 2.;
        draggedCenter.set (int (halfSizeW + halfSizeW * (centerX->getValue() / 1000.)), int (halfSizeH + halfSizeH * (centerY->getValue() / 1000.)));

        // trick to get the correct angle (clockwise/counter-clockwise)
        rtengine::Coord p1 = draggedCenter;
        rtengine::Coord p2 = provider->posImage;
        int p = p1.y;
        p1.y = p2.y;
        p2.y = p;
        pCoord = p2 - p1;
        draggedPointOldAngle = pCoord.angle;
        draggedPointAdjusterAngle = degree->getValue();

        if (Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 2) {
            if (lastObject == 2) {
                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;
                double verti = double (imH);
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;
                // compute the projected value of the dragged point
                draggedlocYOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

                if (lastObject == 2) {
                    //draggedlocYOffset = -draggedlocYOffset;
                    draggedlocYOffset -= (locYT->getValue() / 2000. * verti);

                }
            } else if (lastObject == 3) {
                // Dragging a line to change the angle
                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;

                double verti = double (imH);

                // trick to get the correct angle (clockwise/counter-clockwise)
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;

                // draggedPoint.setFromCartesian(centerPos, currPos);
                // compute the projected value of the dragged point
                draggedlocYOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

                if (lastObject == 3) {
                    draggedlocYOffset = -draggedlocYOffset;
                    draggedlocYOffset -= (locY->getValue() / 2000. * verti);

                }

            }

        } else if (Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
            if (lastObject == 2 || lastObject == 3) {
                // Dragging a line to change the angle
                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;
                double verti = double (imH);
                // trick to get the correct angle (clockwise/counter-clockwise)
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;

                //    draggedPoint.setFromCartesian(centerPos, currPos);
                // compute the projected value of the dragged point
                draggedlocYOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

                if (lastObject == 3) {
                    draggedlocYOffset = -draggedlocYOffset;
                }

                draggedlocYOffset -= (locY->getValue() / 2000. * verti);
            }
        }

        if (Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 2) {
            if (lastObject == 0) {
                // Dragging a line to change the angle

                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;

                double horiz = double (imW);

                // trick to get the correct angle (clockwise/counter-clockwise)
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;

                //     draggedPoint.setFromCartesian(centerPos, currPos);
                // compute the projected value of the dragged point
                //printf ("rad=%f ang=%f\n", draggedPoint.radius, draggedPoint.angle - degree->getValue());
                draggedlocXOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);
                //  if (lastObject==1)
                //      draggedlocXOffset = -draggedlocXOffset;//-
                draggedlocXOffset -= (locX->getValue() / 2000. * horiz);
            } else if (lastObject == 1) {

                // Dragging a line to change the angle
                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;
                double horiz = double (imW);
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;

                //     draggedPoint.setFromCartesian(centerPos, currPos);
                // printf ("rad=%f ang=%f\n", draggedPoint.radius, draggedPoint.angle - degree->getValue());
                draggedlocXOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);

                if (lastObject == 1) {
                    draggedlocXOffset = -draggedlocXOffset;    //-
                }

                draggedlocXOffset -= (locXL->getValue() / 2000. * horiz);
            }

        } else if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {

            if (lastObject == 0 || lastObject == 1) {
                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;
                double horiz = double (imW);
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;

                //    draggedPoint.setFromCartesian(centerPos, currPos);
                //printf ("rad=%f ang=%f\n", draggedPoint.radius, draggedPoint.angle - degree->getValue());
                draggedlocXOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);

                if (lastObject == 1) {
                    draggedlocXOffset = -draggedlocXOffset;    //-
                }

                draggedlocXOffset -= (locX->getValue() / 2000. * horiz);
            }
        }

        /*  else if(Smethod->get_active_row_number()==2) {
                if (lastObject==0 || lastObject==1 || lastObject==2 || lastObject==3) {
                if (lastObject==2 || lastObject==3) {
                    // Dragging a line to change the angle
                    PolarCoord draggedPoint;
                    Coord currPos;
                    currPos = provider->posImage;
                    Coord centerPos = draggedCenter;
                    double verti = double(imH);
                    // trick to get the correct angle (clockwise/counter-clockwise)
                    int p = centerPos.y;
                    centerPos.y = currPos.y;
                    currPos.y = p;

                    draggedPoint.setFromCartesian(centerPos, currPos);
                    // compute the projected value of the dragged point
                    draggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle-degree->getValue())/180.*M_PI);
                    if (lastObject==3)
                        draggedlocYOffset = -draggedlocYOffset;
                    draggedlocYOffset -= (locY->getValue() / 200. * verti);
                }


                if (lastObject==0 || lastObject==1) {
                    PolarCoord draggedPoint;
                    Coord currPos;
                    currPos = provider->posImage;
                    Coord centerPos = draggedCenter;
                    double horiz = double(imW);
                    int p = centerPos.y;
                    centerPos.y = currPos.y;
                    currPos.y = p;
                    draggedPoint.setFromCartesian(centerPos, currPos);
                    printf("rad=%f ang=%f\n",draggedPoint.radius,draggedPoint.angle-degree->getValue());
                    draggedlocXOffset = draggedPoint.radius * sin((draggedPoint.angle-degree->getValue()+90.)/180.*M_PI);
                    if (lastObject==1)
                        draggedlocXOffset = -draggedlocXOffset;//-
                    draggedlocXOffset -= (locX->getValue() / 200. * horiz);
                }

                }
            }
            */
        //    EditSubscriber::dragging = true;
        EditSubscriber::action = ES_ACTION_DRAGGING;
        return false;
    } else {
        // this will let this class ignore further drag events
        if (lastObject > -1) { // should theoretically always be true
            if (lastObject == 2 || lastObject == 3) {
                EditSubscriber::visibleGeometry.at (2)->state = Geometry::NORMAL;
                EditSubscriber::visibleGeometry.at (3)->state = Geometry::NORMAL;
            }

            if (lastObject == 0 || lastObject == 1) {
                EditSubscriber::visibleGeometry.at (0)->state = Geometry::NORMAL;
                EditSubscriber::visibleGeometry.at (1)->state = Geometry::NORMAL;

            } else {
                EditSubscriber::visibleGeometry.at (4)->state = Geometry::NORMAL;
//               EditSubscriber::visibleGeometry.at (lastObject)->state = Geometry::NORMAL;
            }
        }

        lastObject = -1;
        return true;
    }
}

bool Localrgb::button1Released()
{
    draggedPointOldAngle = -1000.;
    EditSubscriber::action = ES_ACTION_NONE;

    return true;
}

bool Localrgb::drag1 (int modifierKey)
{
    // compute the polar coordinate of the mouse position
    EditDataProvider *provider = getEditProvider();
    int imW, imH;
    provider->getImageSize (imW, imH);
    double halfSizeW = imW / 2.;
    double halfSizeH = imH / 2.;

    if (Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 2) {
        if (lastObject == 2) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double verti = double (imH);
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            //  draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedlocYOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

            if (lastObject == 2) {
                currDraggedlocYOffset -= draggedlocYOffset;
            }

            //else if (lastObject==3)
            // Dragging the lower locY bar
            //  currDraggedlocYOffset = -currDraggedlocYOffset + draggedlocYOffset;
            currDraggedlocYOffset = currDraggedlocYOffset * 2000. / verti;

            if (int (currDraggedlocYOffset) != locYT->getIntValue()) {
                locYT->setValue ((int (currDraggedlocYOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();
                updateGeometry (centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue() );

                if (listener) {
                    listener->panelChanged (EvlocalrgblocY, locYT->getTextValue());
                }

                return true;
            }
        } else if (lastObject == 3) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double verti = double (imH);
            // trick to get the correct angle (clockwise/counter-clockwise)
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            //  draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedlocYOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

            //  if (lastObject==2)
            // Dragging the upper locY bar
            //      currDraggedlocYOffset -= draggedlocYOffset;
            //  else
            if (lastObject == 3)
                // Dragging the lower locY bar
            {
                currDraggedlocYOffset = -currDraggedlocYOffset + draggedlocYOffset;
            }

            currDraggedlocYOffset = currDraggedlocYOffset * 2000. / verti;

            if (int (currDraggedlocYOffset) != locY->getIntValue()) {

                locY->setValue ((int (currDraggedlocYOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();

                updateGeometry (centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged (EvlocalrgblocY, locY->getTextValue());
                }

                return true;
            }
        }

    } else if (Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
        if (lastObject == 2 || lastObject == 3) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double verti = double (imH);
            // trick to get the correct angle (clockwise/counter-clockwise)
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            //   draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedlocYOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue()) / 180.*rtengine::RT_PI);

            if (lastObject == 2)
                // Dragging the upper locY bar
            {
                currDraggedlocYOffset -= draggedlocYOffset;
            } else if (lastObject == 3)
                // Dragging the lower locY bar
            {
                currDraggedlocYOffset = -currDraggedlocYOffset + draggedlocYOffset;
            }

            currDraggedlocYOffset = currDraggedlocYOffset * 2000. / verti;

            if (int (currDraggedlocYOffset) != locY->getIntValue()) {
                locY->setValue ((int (currDraggedlocYOffset)));
                //Smethod->get_active_row_number()==2
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();

                updateGeometry (centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(),  locYT->getValue(), locXL->getValue());

                if (listener) {
                    if (Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
                        listener->panelChanged (EvlocalrgblocY, locY->getTextValue());
                    }

                    //  else listener->panelChanged (EvlocallablocY, locX->getTextValue());

                }

                return true;
            }
        }

    }

    if (Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 2) {
        //else if (lastObject==0) {
        if (lastObject == 0) {// >=4
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double horiz = double (imW);
            // trick to get the correct angle (clockwise/counter-clockwise)
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            //    draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedStrOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);

            if (lastObject == 0) //>=4
                // Dragging the upper locY bar
            {
                currDraggedStrOffset -= draggedlocXOffset;
            } else if (lastObject == 1)
                // Dragging the lower locY bar
            {
                currDraggedStrOffset = - currDraggedStrOffset - draggedlocXOffset;    //-
            }

            currDraggedStrOffset = currDraggedStrOffset * 2000. / horiz;

            if (int (currDraggedStrOffset) != locX->getIntValue()) {
                locX->setValue ((int (currDraggedStrOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();
                updateGeometry (centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged (EvlocalrgblocX, locX->getTextValue());
                }

                return true;
            }
        } else if (lastObject == 1) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double horiz = double (imW);
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            //draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedStrOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);

            if (lastObject == 0)
                // Dragging the upper locY bar
            {
                currDraggedStrOffset -= draggedlocXOffset;
            } else if (lastObject == 1)
                // Dragging the lower locY bar
            {
                currDraggedStrOffset = - currDraggedStrOffset - draggedlocXOffset;    //-
            }

            currDraggedStrOffset = currDraggedStrOffset * 2000. / horiz;

            if (int (currDraggedStrOffset) != locXL->getIntValue()) {
                locXL->setValue ((int (currDraggedStrOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();
                updateGeometry (centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged (EvlocalrgblocX, locX->getTextValue());
                }

                return true;
            }
        }

    } else if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
        if (lastObject == 0 || lastObject == 1) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double horiz = double (imW);
            // trick to get the correct angle (clockwise/counter-clockwise)
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            // draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedStrOffset = draggedPoint.radius * sin ((draggedPoint.angle - degree->getValue() + 90.) / 180.*rtengine::RT_PI);

            if (lastObject == 0)
                // Dragging the upper locY bar
            {
                currDraggedStrOffset -= draggedlocXOffset;
            } else if (lastObject == 1)
                // Dragging the lower locY bar
            {
                currDraggedStrOffset = - currDraggedStrOffset - draggedlocXOffset;    //-
            }

            currDraggedStrOffset = currDraggedStrOffset * 2000. / horiz;

            if (int (currDraggedStrOffset) != locX->getIntValue()) {
                locX->setValue ((int (currDraggedStrOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();
                updateGeometry (centX, centY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged (EvlocalrgblocX, locX->getTextValue());
                }

                return true;
            }
        }
    }

    /*  else if(Smethod->get_active_row_number()==2) {
            if (lastObject==0 || lastObject==1 || lastObject==2 || lastObject==3) {
        if (lastObject==2 || lastObject==3) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            Coord currPos;
            currPos = provider->posImage+provider->deltaImage;
            Coord centerPos = draggedCenter;
            double verti = double(imH);
            // trick to get the correct angle (clockwise/counter-clockwise)
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle-degree->getValue())/180.*M_PI);
            double currDraggedStrOffset = draggedPoint.radius * sin((draggedPoint.angle-degree->getValue() +90.)/180.*M_PI);

            if (lastObject==2)
                currDraggedlocYOffset -= draggedlocYOffset;
            else if (lastObject==3)
                currDraggedlocYOffset = -currDraggedlocYOffset + draggedlocYOffset;
            currDraggedlocYOffset = currDraggedlocYOffset * 200. / verti;
        //  if (int(currDraggedlocYOffset) != locY->getIntValue()) {
        //      locY->setValue((int(currDraggedlocYOffset)));
            if (int(currDraggedlocYOffset) != locX->getIntValue()) {//locX
        //  if (int(currDraggedStrOffset) != locX->getIntValue()) {//locX
                locX->setValue((int(currDraggedlocYOffset)));
                double centX,centY;
                centX=centerX->getValue();
                centY=centerY->getValue();

            //  updateGeometry (centX, centY, locY->getValue(), degree->getValue(), locX->getValue(),  locYT->getValue(), locXL->getValue());
                updateGeometry (centX, centY, locX->getValue(), degree->getValue(), locX->getValue(),  locX->getValue(), locX->getValue());
                if (listener) {
                    if(Smethod->get_active_row_number()==1) listener->panelChanged (EvlocallablocY, locY->getTextValue());

                    }
                return true;
            }
        }
            if (lastObject==0 || lastObject==1) {
                // Dragging the upper or lower locY bar
                PolarCoord draggedPoint;
                Coord currPos;
                currPos = provider->posImage+provider->deltaImage;
                Coord centerPos = draggedCenter;
                double horiz = double(imW);
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint.setFromCartesian(centerPos, currPos);
                double currDraggedStrOffset = draggedPoint.radius * sin((draggedPoint.angle-degree->getValue() +90.)/180.*M_PI);
                if (lastObject==0)
                    currDraggedStrOffset -= draggedlocXOffset;
                else if (lastObject==1)
                    currDraggedStrOffset = - currDraggedStrOffset - draggedlocXOffset;//-
                    currDraggedStrOffset = currDraggedStrOffset * 200. / horiz;

                if (int(currDraggedStrOffset) != locX->getIntValue()) {
                    locX->setValue((int(currDraggedStrOffset)));
                    double centX,centY;
                    centX=centerX->getValue();
                    centY=centerY->getValue();
                    updateGeometry (centX, centY, locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(),locXL->getValue());
                    if (listener)
                        listener->panelChanged (EvlocallablocX, locX->getTextValue());
                    return true;
                }
            }


            }
        }
        */
    //else if (lastObject==4) {
    if (lastObject == 4) {

        // Dragging the circle to change the center
        rtengine::Coord currPos;
        draggedCenter += provider->deltaPrevImage;
        currPos = draggedCenter;
        currPos.clip (imW, imH);
        int newCenterX = int ((double (currPos.x) - halfSizeW) / halfSizeW * 1000.);
        int newCenterY = int ((double (currPos.y) - halfSizeH) / halfSizeH * 1000.);

        if (newCenterX != centerX->getIntValue() || newCenterY != centerY->getIntValue()) {
            centerX->setValue (newCenterX);
            centerY->setValue (newCenterY);
            updateGeometry (newCenterX, newCenterY, circrad->getValue(), locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

            if (listener) {
                listener->panelChanged (EvlocalrgbCenter, Glib::ustring::compose ("X=%1\nY=%2", centerX->getTextValue(), centerY->getTextValue()));
            }

            return true;
        }
    }

    return false;
}

void Localrgb::switchOffEditMode ()
{
    if (edit->get_active()) {
        // switching off the toggle button
        bool wasBlocked = editConn.block (true);
        edit->set_active (false);

        if (!wasBlocked) {
            editConn.block (false);
        }
    }

    EditSubscriber::switchOffEditMode();  // disconnect
}

