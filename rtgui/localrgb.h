/*
 *  This file is part of RawTherapee.
 */

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "edit.h"
#include "guiutils.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "toolpanel.h"
#include "../rtengine/imagedata.h"
#include <memory>
#include "options.h"
#include <string>
#include "../rtengine/improcfun.h"


class Localrgb :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel,
    public rtengine::localListener,
    public CurveListener,
    public EditSubscriber,
    public ColorProvider,
    public rtengine::localrgbListener

{
private:
    int lastObject;
    void foldAllButMe (GdkEventButton* event, MyExpander *expander);
    void enableToggled (MyExpander *expander);

    IdleRegister idle_register;
    IdleRegister idle_register2;

//protected:
    Gtk::HBox *editHBox;
    Gtk::ToggleButton* edit;

    MyExpander* const expexpose;
    MyExpander* const expsettings;
    MyExpander* const expwb;

    Adjuster* nbspot;

    Adjuster* const anbspot;
    Adjuster* const locX;
    Adjuster* const locXL;
    Adjuster* const degree;
    Adjuster* const locY;
    Adjuster* const locYT;
    Adjuster* const centerX;
    Adjuster* const centerY;
    Adjuster* const circrad;
    Adjuster* const thres;
    Adjuster* const proxi;
    Adjuster* const lightness;
    Adjuster* const contrast;
    Adjuster* const chroma;
    Adjuster* const sensi;
    Adjuster* const transit;
    Adjuster* const retrab;
    Adjuster* const expcomp;
    Adjuster* const hlcompr;
    Adjuster* const hlcomprthresh;
    Adjuster* const black;
    Adjuster* const shcompr;
    Adjuster* const hueref;
    Adjuster* const chromaref;
    Adjuster* const lumaref;

    MyComboBoxText*   const Smethod;
    MyComboBoxText*   const qualityMethod;
    MyComboBoxText*   const wbMethod;


    Gtk::Frame* const shapeFrame;
    Gtk::Frame* const artifFrame;
    Gtk::Frame* const superFrame;

    Gtk::Label* const labqual;
    Gtk::Label* const labmS;

    Gtk::HBox* const ctboxS;
    Gtk::HBox* const qualbox;

    Adjuster* temp;
    Adjuster* green;
    Adjuster* equal;
    Gtk::Label* ttLabels;

    MyComboBoxText* toneCurveMode;
    MyComboBoxText* toneCurveMode2;
    sigc::connection tcmodeconn, tcmode2conn;
    CurveEditorGroup* curveEditorG;
    CurveEditorGroup* curveEditorG2;
    DiagonalCurveEditor* shape;
    DiagonalCurveEditor* shape2;


    Gtk::Button* spotbutton;

    sigc::connection enableexposeConn;
    sigc::connection  editConn;
    sigc::connection enablewbConn;

//   Gtk::HBox* const qualcurvbox;


//    Gtk::VBox* const artifVBox;
//    Gtk::VBox* const shapeVBox;
//    Gtk::VBox* const colorVBox;
    sigc::connection  Smethodconn;
    sigc::connection qualityMethodConn;
    sigc::connection wbMethodConn;


    int nextdatasp[61];
    int nextlength;
    std::string nextstr;
    std::string nextstr2;
    std::string nextll_str;
    std::string nextll_str2;
    std::string nextlh_str;
    std::string nextlh_str2;
    std::string nextcc_str;
    std::string nextcc_str2;
    double nexttemp;
    double nexttint;
    double nextequal;
    double next_temp;
    double next_green;
    int next_wbauto;
	int nextmeth;
    double draggedPointOldAngle;
    double draggedPointAdjusterAngle;
    double draggedFeatherOffset;
    double draggedlocYOffset;
    double draggedlocXOffset;
    double draggedlocYTOffset;
    double draggedlocXLOffset;
    rtengine::Coord draggedCenter;
    bool lastavoid, lastinvers, lastinversrad, lastinversret, lastactivlum, lastinverssha, lastcurvactiv;
    int lastanbspot;

    void editToggled ();

public:

    Localrgb ();
    ~Localrgb ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);

    void setBatchMode   (bool batchMode);

    void updateGeometry (const int centerX_, const int centerY_, const int circrad_, const int locY_, const double degree_, const int locX_, const int locYT_, const int locXL_, const int fullWidth = -1, const int fullHeight = -1);
    void SmethodChanged      ();
    void writeOptions (std::vector<int> &tpOpen);
    void updateToolState (std::vector<int> &tpOpen);

    void adjusterChanged (Adjuster* a, double newval);
    void enabledChanged ();
    void setAdjusterBehavior (bool degreeadd, bool locYadd, bool locXadd, bool locYTadd, bool locXLadd, bool centeradd, bool lightnessadd, bool contrastadd, bool chromaadd, bool sensiadd, bool transitadd, bool radiusadd, bool strengthadd);
    void trimValues          (rtengine::procparams::ProcParams* pp);
    void avoidChanged ();
    void activlumChanged ();
    void inversChanged ();
    void curvactivChanged ();
    void inversradChanged ();
    void inversretChanged ();
    void inversshaChanged ();
    void curveChanged (CurveEditor* ce);
    void curveMode1Changed ();
    bool curveMode1Changed_ ();
    void curveMode2Changed ();
    bool curveMode2Changed_ ();

//   void autoOpenCurve ();
//   void localChanged           (int **datasp, std::string datastr, std::string ll_str, std::string lh_str, std::string cc_str, int sp, int maxdat);
//   void localretChanged           (int **datasp, std::string datastr, std::string ll_str, std::string lh_str, std::string cc_str, int sp, int maxdat);
    bool localComputed_         ();
    bool localretComputed_         ();
    bool localwbComputed_ ();
    void setEditProvider (EditDataProvider* provider);
    void retinexMethodChanged();
    void qualityMethodChanged();
    void wbMethodChanged();
    void qualitycurveMethodChanged();
    void lumaneutralPressed ();
    void lumacontrastPlusPressed ();
    void lumacontrastMinusPressed ();
    void neutral_pressed       ();
    virtual void colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller);
    void temptintChanged (double ctemp, double ctint, double cequal, int meth);
    bool temptintComputed_ ();
    void updateLabel      ();
    void WBChanged           (double temp, double green, int wbauto);

    // EditSubscriber interface
    CursorShape getCursor (int objectID);
    bool mouseOver (int modifierKey);
    bool button1Pressed (int modifierKey);
    bool button1Released();
    bool drag1 (int modifierKey);
    void switchOffEditMode ();
};

