#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Resize Controls>
#include <Rewrite Control Positions>


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 



// creates selection in analysis pull down menu
// /q stops execution of menu opening from being stored in command window

menu "Optical Constants"
	"Load & Calculate Panel",/q, opencalculatepanel()
end
 
// creates panel, defines operations
// this function executes when the menu item optical constants/load & calculate panel is chosen 
 
function opencalculatepanel()

// there can be conflicts when waves that are used in graphs or tables
// are attempted to be killed; begin by closing all graphs and tables with this
// homemade command

setdatafolder root:
killallgraphsandtablesNEW()

// begin by closing "calculate panel" if it is open, then recreating it
// make space for the panel to live, name it

dowindow/k calculatepanel
newpanel /w=(1063,0,1438,1000)/n=calculatepanel
drawtext 53,24,"Enter Experimental Parameters and Load Waves"
drawtext 225,690,"Funded by Center" 
drawtext 225, 705, "For Emergent Materials" 
drawtext 225, 720, "NSF grant DMR-1420451"
	
// to test what has been done here, you can place the line "return -1" and 
// see what the above commands do; this is generally true
// or click on the left bar of this window, placing a red dot and make sure "enable debugger"
// is chosen in the procedure menu item; then you can iterate through the code by using the yellow arrow,
// beginning at the red dot

	
// create a button, with title "Load Waves" that execeutes the procedure loadwavesbuttonproc when pressed		

	Button loadwavesbutton,pos={147,135},size={90,40},proc=loadwavesbuttonproc,title="Load Waves"
	SetVariable ltempctrl title="custom lower temp ",pos={245,138},size={120,15},proc=SetVarProc, value=_STR:"", disable=2
	SetVariable utempctrl title="custom upper temp ",pos={245,158},size={120,15},proc=SetVarProc, value=_STR:"",disable=2
	
// if using panel/show tools to set the positions by hand, the top left height width goes like pos{left,top}, size{height, width}	

// the first name is the name of the control
// these are all input parameters, so the schematic is to use the general name that is used in the program
// sam for sample name, for example, with the suffix ctrl
// the limits are optional, but used because using the limits and setting the increment to zero eliminates the arrows
 	
	SetVariable samctrl title="sample ",pos={8,38},size={125,15},proc=SetVarProc, value=_STR:""
	SetVariable refctrl title="reference ",pos={147,38},size={136,15},proc=SetVarProc, value=_STR:""
	SetVariable stepctrl title="step ",pos={303,38},size={63,15},proc=SetVarProc, value=_num:0, limits={-inf,inf,0}
	SetVariable sspctrl title="sample start position ",pos={7,65},size={173,15},proc=SetVarProc, value=_num:0, limits={-inf,inf,0}
	SetVariable sepctrl title="sample end position ",pos={199,67},size={167,15},proc=SetVarProc, value=_num:0, limits={-inf,inf,0}
	SetVariable rspctrl title="ref start position ",pos={24,87},size={156,15},proc=SetVarProc, value=_num:0, limits={-inf,inf,0}
	SetVariable repctrl title="ref end position ",pos={218,86},size={148,15},proc=SetVarProc, value=_num:0, limits={-inf,inf,0}
	SetVariable tlctrl title="temperature list ",pos={8,110},size={358,15},proc=SetVarProc, value=_STR:"temperature list;separated by ';' ex. 10;20 "

	DrawText 125,195,"Calculational Parameters"

	SetVariable decctrl title="windowing decimal",pos={8,205},size={115,15},proc=SetVarProc, value=_num:0, limits={-inf,inf,0}
	SetVariable dctrl title="thickness (um) ",pos={130,205},size={115,15},proc=SetVarProc, value=_num:0, limits={-inf,inf,0}
	SetVariable freqsctrl title="# freqs to display",pos={253,205},size={115,15},proc=SetVarProc, value=_num:0, limits={-inf,inf,0}
	SetVariable nspctrl title="pad start position ",pos={14,230},size={167,15},proc=SetVarProc, value=_num:0, limits={-inf,inf,0}
	SetVariable nepctrl title="pad end position ",pos={200,230},size={156,15},proc=SetVarProc, value=_num:0, limits={-inf,inf,0}
	SetVariable lfreqctrl title="lower freq",pos={23,255},size={85,15},proc=SetVarProc, value=_num:0, limits={-inf,inf,0}
	SetVariable ufreqctrl title="upper freq",pos={125,255},size={85,15},proc=SetVarProc, value=_num:0, limits={-inf,inf,0}
	SetVariable nctrl title="index guess (real only)",pos={225,255},size={125,15},proc=SetVarProc, value=_num:0, limits={-inf,inf,0}
	
	Button bringgraphstofrontbutton,pos={22,280},size={90,40},fsize=12,proc=bringgraphstofrontbuttonproc,title="Bring Graphs \rTo Front"
	Button calculatebutton,pos={122,280},size={141,40},proc=slabcalcbuttonproc,title="Calculate"
	Button printbutton,pos={274,280},size={81,40},proc=reportproc,title="Print \r Graphs",disable=2
	
	drawtext 22,340,"Graphs"
	
// cbstr is used here so that the default graphs (selected by setting value=1 below) are opened when the calculate button is hit
// the command is also so written that it keeps graphs for all checked boxes
// this is unlike other reopening of graphs, like when a log y is chosen, that typically already have graphs open when they execute
// since the calculate button can start from no graphs present, but a number of checked boxes, it is special	
	
	string/g cbstr = "sigma2freqctrl,eps2freqctrl,eps2tempctrl,nfreqctrl,kfreqctrl,ntempctrl,ktempctrl,sigma1freqctrl,eps1freqctrl,sigma1tempctrl,eps1tempctrl,trfreqctrl,trtempctrl,alphafreqctrl,alphatempctrl,eps1imagectrl,eps2imagectrl,nimagectrl,kimagectrl,sigma1imagectrl,sigma2imagectrl,alphaimagectrl"
	checkbox eps1freqctrl title="\Z13\F'symbol'e\F'geneva'\B1\M\Z11 vs f", pos={14,345}, size={115,45},proc=checkboxproc, value=1, userdata = "eps1freqgraph"
	checkbox eps2freqctrl title="\Z13\F'symbol'e\F'geneva'\B2\M \Z11vs f", pos={14,365}, size={115,45},proc=checkboxproc, value=0, userdata = "eps2freqgraph"
	checkbox eps1tempctrl title="\Z13\F'symbol'e\F'geneva'\B1\M \Z11vs T", pos={14,385}, size={115,45},proc=checkboxproc, value=1, userdata = "eps1tempgraph"
	checkbox eps2tempctrl title="\Z13\F'symbol'e\F'geneva'\B2\M \Z11vs T", pos={14,405}, size={115,45},proc=checkboxproc, value=0, userdata = "eps2tempgraph"
	checkbox nfreqctrl title="\Z11n vs f", pos={14,425}, size={115,45},proc=checkboxproc, value=0, userdata = "nfreqgraph"
	checkbox kfreqctrl title="\Z11k vs f", pos={14,445}, size={115,45},proc=checkboxproc, value=0, userdata = "kfreqgraph"
	checkbox ntempctrl title="\Z11n vs T", pos={14,465}, size={115,45},proc=checkboxproc, value=0, userdata = "ntempgraph"
	checkbox ktempctrl title="\Z11k vs T", pos={14,485}, size={115,45},proc=checkboxproc, value=0, userdata = "ktempgraph"
	checkbox sigma1freqctrl title="\Z13\F'symbol's\F'geneva'\B1\M\Z11 vs f", pos={14,505}, size={115,45},proc=checkboxproc, value=1, userdata = "sigma1freqgraph"
	checkbox sigma2freqctrl title="\Z13\F'symbol's\F'geneva'\B2\M \Z11vs f", pos={14,525}, size={115,45},proc=checkboxproc, value=0, userdata = "sigma2freqgraph"
	checkbox sigma1tempctrl title="\Z13\F'symbol's\F'geneva'\B1\M \Z11vs T", pos={14,545}, size={115,45},proc=checkboxproc, value=1, userdata = "sigma1tempgraph"
	checkbox sigma2tempctrl title="\Z13\F'symbol's\F'geneva'\B2\M \Z11vs T", pos={14,565}, size={115,45},proc=checkboxproc, value=0, userdata = "sigma2tempgraph"
	checkbox trfreqctrl title="\Z11tr vs f", pos={14,585}, size={115,45},proc=checkboxproc, value=1, userdata = "trfreqgraph"
	checkbox trtempctrl title="\Z11tr vs T", pos={14,605}, size={115,45},proc=checkboxproc, value=1, userdata = "trtempgraph"	
	checkbox alphafreqctrl title="\Z13\F'symbol'a\F'geneva'\B\M vs f", pos={14,625}, size={115,45},proc=checkboxproc, value=0, userdata = "alphafreqgraph"
	checkbox alphatempctrl title="\Z13\F'symbol'a\F'geneva'\B\M  vs T", pos={14,645}, size={115,45},proc=checkboxproc, value=0, userdata = "alphatempgraph"	
	drawtext 16,680,"VRH dim:"
	checkbox vrhdim1ctrl title="1:", pos={-12,665}, size={115,45},proc=vrhproc, value=0, side=1,fsize=10
	checkbox vrhdim2ctrl title="2:", pos={18,665}, size={115,45},proc=vrhproc, value=0, side=1,fsize=10
	checkbox vrhdim3ctrl title="3:", pos={48,665}, size={115,45},proc=vrhproc, value=0, side=1,fsize=10
	checkbox ratiotolowtempctrl title="ratio to low temp", pos={14,685}, size={115,45},proc=ratiotolowtempproc, value=0, fsize=10
	checkbox losstanctrl title="loss tan", pos={14,705}, size={115,45},proc=losstanproc, value=0, fsize=10
	checkbox zsctrl title="surface impedance", pos={14,725}, size={115,45},proc=zsproc, value=0, fsize=10
	checkbox rctrl title="reflectivity", pos={14,745}, size={115,45},proc=rproc, value=0, fsize=10
	checkbox energylossctrl title="energy loss", pos={135,685}, size={115,45},proc=energylossproc, value=0, fsize=10	
	checkbox ebphasenctrl title="EB Phase", pos={135,705}, size={115,45},proc=ebphaseproc, value=0, fsize=10
	
// when loading the panel, set all logy variables to zero
	
	drawtext 81,340,"Log y"
	
// a special folder is made for the logy (and logx, later) global variables	
// a function is run to create them all and initialize them all to zero
// overwritefolder is a user created function that kills a folder if it exists, and then recreates it
// this is useful for resetting things
	
	string/g logypath = "root:globals:logy:"
       overwritefolder(logypath)
       setlogyvartozero()
       
// the log checkboxes call specific procedures when clicked
       
       checkbox eps1freqlogyctrl title="", pos={90,345}, size={115,45},proc=logycheckboxproc, value=0
       checkbox eps2freqlogyctrl title="", pos={90,365}, size={115,45},proc=logycheckboxproc, value=0
       checkbox eps1templogyctrl title="", pos={90,385}, size={115,45},proc=logycheckboxproc, value=0
       checkbox eps2templogyctrl title="", pos={90,405}, size={115,45},proc=logycheckboxproc, value=0
       checkbox nfreqlogyctrl title="", pos={90,425}, size={115,45},proc=logycheckboxproc, value=0
       checkbox kfreqlogyctrl title="", pos={90,445}, size={115,45},proc=logycheckboxproc, value=0
	checkbox ntemplogyctrl title="", pos={90,465}, size={115,45},proc=logycheckboxproc, value=0
	checkbox ktemplogyctrl title="", pos={90,485}, size={115,45},proc=logycheckboxproc, value=0
	checkbox sigma1freqlogyctrl title="", pos={90,505}, size={115,45},proc=logycheckboxproc, value=0
	checkbox sigma2freqlogyctrl title="", pos={90,525}, size={115,45},proc=logycheckboxproc, value=0
	checkbox sigma1templogyctrl title="", pos={90,545}, size={115,45},proc=logycheckboxproc, value=0
	checkbox sigma2templogyctrl title="", pos={90,565}, size={115,45},proc=logycheckboxproc, value=0
	checkbox trfreqlogyctrl title="", pos={90,585}, size={115,45},proc=logycheckboxproc, value=0
	checkbox trtemplogyctrl title="", pos={90,605}, size={115,45},proc=logycheckboxproc, value=0
	checkbox alphafreqlogyctrl title="", pos={90,625}, size={115,45},proc=logycheckboxproc, value=0
	checkbox alphatemplogyctrl title="", pos={90,645}, size={115,45},proc=logycheckboxproc, value=0


	drawtext 125,340,"Log x"
	
	string/g logxpath = "root:globals:logx:"
       overwritefolder(logxpath)
       setlogxvartozero()
       
       checkbox eps1freqlogxctrl title="", pos={132,345}, size={115,45},proc=logxcheckboxproc, value=0
       checkbox eps2freqlogxctrl title="", pos={132,365}, size={115,45},proc=logxcheckboxproc, value=0
       checkbox eps1templogxctrl title="", pos={132,385}, size={115,45},proc=logxcheckboxproc, value=0
       checkbox eps2templogxctrl title="", pos={132,405}, size={115,45},proc=logxcheckboxproc, value=0
       checkbox nfreqlogxctrl title="", pos={132,425}, size={115,45},proc=logxcheckboxproc, value=0
       checkbox kfreqlogxctrl title="", pos={132,445}, size={115,45},proc=logxcheckboxproc, value=0
	checkbox ntemplogxctrl title="", pos={132,465}, size={115,45},proc=logxcheckboxproc, value=0
	checkbox ktemplogxctrl title="", pos={132,485}, size={115,45},proc=logxcheckboxproc, value=0
	checkbox sigma1freqlogxctrl title="", pos={132,505}, size={115,45},proc=logxcheckboxproc, value=0
	checkbox sigma2freqlogxctrl title="", pos={132,525}, size={115,45},proc=logxcheckboxproc, value=0
	checkbox sigma1templogxctrl title="", pos={132,545}, size={115,45},proc=logxcheckboxproc, value=0
	checkbox sigma2templogxctrl title="", pos={132,565}, size={115,45},proc=logxcheckboxproc, value=0
	checkbox trfreqlogxctrl title="", pos={132,585}, size={115,45},proc=logxcheckboxproc, value=0
	checkbox trtemplogxctrl title="", pos={132,605}, size={115,45},proc=logxcheckboxproc, value=0
	checkbox alphafreqlogxctrl title="", pos={132,625}, size={115,45},proc=logxcheckboxproc, value=0
	checkbox alphatemplogxctrl title="", pos={132,645}, size={115,45},proc=logxcheckboxproc, value=0


	drawtext 213,340,"Image Plots"
	
	// the mode is the default color chosen; see popupmenu documentation for syntax and indexing of colors
	// the lines after extract the default mode and initizialize imagecolor to this name (again, see popupmenu)
	// this is used at first in the calculate button to set the image plots' colortable
	// colortablepopupproc will also update graphs with newly selected color tables
	
	popupmenu imageplotpopup pos={198,345},bodywidth=80,proc=colortablepopupproc,value="*COLORTABLEPOPNONAMES*",mode=7
	popupmenu contourpopup pos={283,345},bodywidth=80,proc=colortablepopupproc,value="*COLORTABLEPOPNONAMES*",mode=17
	controlinfo imageplotpopup
	string/g imagecolor = stringfromlist(v_value-1,ctablist())
//	variable/g imagenum = v_value
	controlinfo contourpopup
	string/g contourcolor = stringfromlist(v_value-1,ctablist())
	
	drawtext 240,387,"Log y"
	drawtext 280,387,"Log x"
	
	checkbox eps1imagectrl title="\Z15\F'symbol'e\F'geneva'\B1\M", pos={190,392}, size={115,45},proc=checkboxproc, value=0,userdata="eps1imagegraph"
	checkbox eps2imagectrl title="\Z15\F'symbol'e\F'geneva'\B2\M", pos={190,417}, size={115,45},proc=checkboxproc, value=0,userdata="eps2imagegraph"
	checkbox nimagectrl title="\Z15n", pos={190,442}, size={115,45},proc=checkboxproc, value=0,userdata="nimagegraph"
	checkbox kimagectrl title="\Z15k", pos={190,467}, size={115,45},proc=checkboxproc, value=0,userdata="kimagegraph"
	checkbox sigma1imagectrl title="\Z15\F'symbol's\F'geneva'\B1\M", pos={190,492}, size={115,45},proc=checkboxproc, value=0,userdata="sigma1imagegraph"
       checkbox sigma2imagectrl title="\Z15\F'symbol's\F'geneva'\B2\M", pos={190,517}, size={115,45},proc=checkboxproc, value=0,userdata="sigma2imagegraph"
	checkbox alphaimagectrl title="\Z15\F'symbol'a\F'geneva'\B\M", pos={190,542}, size={115,45},proc=checkboxproc, value=0,userdata="alphaimagegraph"
	
	checkbox eps1imagelogyctrl title="", pos={245,392}, size={115,45},proc=logycheckboxproc, value=0
	checkbox eps2imagelogyctrl title="", pos={245,417}, size={115,45},proc=logycheckboxproc, value=0
	checkbox nimagelogyctrl title="", pos={245,442}, size={115,45},proc=logycheckboxproc, value=0
	checkbox kimagelogyctrl title="", pos={245,467}, size={115,45},proc=logycheckboxproc, value=0
	checkbox sigma1imagelogyctrl title="", pos={245,492}, size={115,45},proc=logycheckboxproc, value=0
	checkbox sigma2imagelogyctrl title="", pos={245,517}, size={115,45},proc=logycheckboxproc, value=0
	checkbox alphaimagelogyctrl title="", pos={245,542}, size={115,45},proc=logycheckboxproc, value=0
	
	
	checkbox eps1imagelogxctrl title="", pos={285,392}, size={115,45},proc=logxcheckboxproc, value=0
	checkbox eps2imagelogxctrl title="", pos={285,417}, size={115,45},proc=logxcheckboxproc, value=0
	checkbox nimagelogxctrl title="", pos={285,442}, size={115,45},proc=logxcheckboxproc, value=0
	checkbox kimagelogxctrl title="", pos={285,467}, size={115,45},proc=logxcheckboxproc, value=0
	checkbox sigma1imagelogxctrl title="", pos={285,492}, size={115,45},proc=logxcheckboxproc, value=0
	checkbox sigma2imagelogxctrl title="", pos={285,517}, size={115,45},proc=logxcheckboxproc, value=0
	checkbox alphaimagelogxctrl title="", pos={285,542}, size={115,45},proc=logxcheckboxproc, value=0
	
	// tempticksctrl is disabled until temperatures are entered in the "temperature list" input
	// this make sense in normal use, but when modifying the panel, it is possible to have temperatures input without having 
	// entered them in *that* version of the panel. i.e. a tl string is not deleted when the panel is recompiled.
	
	setvariable tempticksctrl title="temperature ticks",pos={60,760},size={275,15},proc=tempticksproc, value=_STR:"disabled until temp list entered", disable=2
	
	// notable temps help identify temperatures of interest on image plots
	
	string/g notabletemps = ""
	setvariable notabletempsctrl title="notable temps",pos={170,570},size={150,15},proc=setvarproc, value=_STR:""
	
	// thickness of line at notable temperature
	// see popupmenu documentation for info on *LINESTYLEPOP*, mode, etc.
	
	variable/g notabledash = 12
	popupmenu notablelinepopup pos={200,595},bodywidth=80,proc=linestylepopupproc,value="*LINESTYLEPOP*",mode=notabledash
	
	// color of a line at the notable temperature
	
	popupmenu colorpopup pos={290,595},bodywidth=80,proc=colorpopupproc,value="*COLORPOP*",popcolor=(0,0,0)
	controlinfo colorpopup
	// there are other ways to extract this information, but s_value returns as a string of three numbers within paretheses
	// when pulling strings from this list, stringfromlist command will not eliminate the parentheses
	// the following lines take that string, eliminate the first character "(" by only taking characters [1,inf]
	// then it removes the ending as well
	string/g notablecolor = s_value
	notablecolor = notablecolor[1,inf]
	notablecolor = removeending(notablecolor)
	
	// initialize the thickness of the notable line
	variable/g notablethickness = .5
	setvariable notablethicknessctrl title="line thickness",pos={199,625},size={90,15},proc=setvarproc, value=_num:notablethickness, limits={-inf,inf,0}
	
	
	button derivbutton,pos={160,625},size={160,20},proc=derivproc,title="Differentiate Something"
	button derivhelpbutton,pos={325,625},size={18,20},proc=derivhelpproc,title="?"
	
// i, j, and k are created as globals when the panel is opened (i.e. here)
// these variables are intended to be used as looping variables and are created globally so that they don't 
// have to be defined within other functions	
// also define npnts which is useful in temperature loops

variable/g i, j, k
               	
// also, create constants: (units are all kg, um, ps)

variable/g eps0 = 8.854187817*10^-12 
variable/g ii = sqrt(-1)
variable/g c = 299.792458
variable/g z0 = 376.73031	

// return 0 just ensures exit of this function
// doesn't seem to make much difference, but people on igorexchange tend to use it
	 
return 0
end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// procedures used in panel:

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function setvarproc(sva) : setvariablecontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// the addition of : setvariablecontrol in the function command tells igor that input will be interpreted
// as per the documentation for the set variable control
// i.e. in conjunction with the struct usage, many parameters are extractable
// from sva (sva.eventcode, sva.sval, etc.)

// the ampersand below establishes the struct environment and 
// the type of object sva is (this allows extraction of parameters from sva
// via the struct syntax, i.e. sva.eventcode yields the eventcode)

        struct wmsetvariableaction &sva
 
// the use of switch will evaluate different cases based on the manipulation
// of the control from the panel 
        
// case 2 is enter, case 3 is live update
// case 3 has been working fine so far
 
        switch(sva.eventcode)
        
               case 1:
               case 2:
               case 3:
               
// this is one giant if statement, comparing the control name of the control altered (sva.ctrlname) to  
// various control names. cmpstr returns 0 when the two strings being compared are equal
// in that case, the code until the next part of the if-elseif statement runs, usually just creating and defining 
// a global string or variable
// note that to extract numbers, variable/g is used with sva.dval, not sva.sval

                  // sam: sample name (string)
                  if(cmpstr(sva.ctrlname,"samctrl")==0)
                  string/g sam = sva.sval 
                  
                  // ref: reference name (string)               
                  elseif(cmpstr(sva.ctrlname,"refctrl")==0)
                  string/g ref = sva.sval
                  
                  // step: delay stage increments (don't multiply by two for incoming and outgoing shift) (variable)
                  elseif(cmpstr(sva.ctrlname,"stepctrl")==0)
                  variable/g step = sva.dval
                  
                  // ssp (etc.) = delay stage readouts for system start point (variable) (rep = reference end point; usually the end point is 5 less than what you set it to be
                  // so that there is an even number of points, so igor can do fft; if in doubt, check the very last point of the raw data - not what the summary at the top says,
                  // but the actual very last point) 
                  elseif(cmpstr(sva.ctrlname,"sspctrl")==0)
                  variable/g ssp = sva.dval
                  
                  elseif(cmpstr(sva.ctrlname,"sepctrl")==0)
                  variable/g sep = sva.dval
                  
                  elseif(cmpstr(sva.ctrlname,"rspctrl")==0)
                  variable/g rsp = sva.dval
                  
                  elseif(cmpstr(sva.ctrlname,"repctrl")==0)
                  variable/g rep = sva.dval
                 
                 // tl: temperature list (string, separated by ";") 
                  // the tlctrl if statement is different to update the tempticks setvariable box with the min and max
                  // values of the temperature list; this way without any user effort your graphs at least have the endpoints on
                  // the temperature color bar                  
                  elseif(cmpstr(sva.ctrlname,"tlctrl")==0)
                  string/g tl = sva.sval
                  variable/g npnts = itemsinlist(tl)
                  setdatafolder root:
		     string/g tempticks = ""
		     tempticks = addlistitem(stringfromlist(itemsinlist(tl)-1,tl),tempticks)
		     tempticks = addlistitem(stringfromlist(0,tl),tempticks)

                  // the tempticksctrl was previously disabled (see opencalculatepanel to see its initial disabled state)
                  // this command makes it usable for input
		     
		     setvariable tempticksctrl title="temperature ticks",pos={60,760},size={275,15},proc=tempticksproc, value=_STR:tempticks,disable=0
                  
                  // make a global temperature wave
                  // use this program; it has an if statement to convert
                  // non-decimal temperatures, like 7p6 to 7.6
                   temperaturewaveNEW()
                  
                  // decimal: percent from edges that exponentially decays (variable)
                  // see trwindower function for use
                  elseif(cmpstr(sva.ctrlname,"decctrl")==0)
                  variable/g dec = sva.dval
                  
                  // d: slab thickness (variable)
                  elseif(cmpstr(sva.ctrlname,"dctrl")==0)
                  variable/g d = sva.dval
                  
                  // freqs: the number of frequencies between high and low frequency points to calculate
                  elseif(cmpstr(sva.ctrlname,"freqsctrl")==0)
                  variable/g freqs = sva.dval
                  
                  // newstartpt (etc.) = endpoints for padded wave (variable)
                  elseif(cmpstr(sva.ctrlname,"nepctrl")==0)
                  variable/g nep = sva.dval
                  
                  elseif(cmpstr(sva.ctrlname,"nspctrl")==0)
                  variable/g nsp = sva.dval
                  
                  // lfreq (etc.): useful if divergent solutions at very low frequency impede solution at more stable frequencies (variable)
                  // lfreq: lower freq; ufreq: upper freq
                  // graphs etc. will reflect these bounds
                  elseif(cmpstr(sva.ctrlname,"lfreqctrl")==0)
                  variable/g lfreq = sva.dval
                  
                  elseif(cmpstr(sva.ctrlname,"ufreqctrl")==0)
                  variable/g ufreq = sva.dval
                  
                  // n: guess for index (variable)
                  elseif(cmpstr(sva.ctrlname,"nctrl")==0)
                  variable/g n = sva.dval
                  
                  // notable temps: temperatures of interest
                  // lines can be shown on image plots to guide the eye to these temperatures
                  // since graphs need to be recreated when new notable temperatures are added,
                  // recalculateandregraph is run to update all the graphs
                  // i'm sure there is a more efficient way to do this, but it only takes ~1 sec to recalculate and regraph everything
                  // and it seems there would be a web of code to introduce to single out which calculations need to be done when this input
                  // is updated and which do not; for now, the calculation time is tolerable
                  elseif(cmpstr(sva.ctrlname,"notabletempsctrl")==0)
                  string/g notabletemps = sva.sval
                  recalculateandregraphNEW()
                  
                  // notable temps line thickness control; see above for recalculation philosophy
                  elseif(cmpstr(sva.ctrlname,"notablethicknessctrl")==0)
                  variable/g notablethickness = sva.dval
                  recalculateandregraphNEW()
		  
		  
	  endif
            
// break is not essential in switches, but recommended        
 
break
                     
endswitch

        return 0
        
end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
  
function loadwavesbuttonproc(ba) : buttoncontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// this function governs what happens when the load waves button is pressed
// it uses "struct" which creates a packet of information 
// here, we use the eventcode associated with the button
// when the button is clicked, the loadwavesbuttonproc runs 
// and based on the eventcode generated, in this example, mouse down and mouse up
// (i.e. a click), different things will happen
// eventcode 2 is mouse up - i.e. the release of a click

struct wmbuttonaction &ba

switch(ba.eventCode)
        
      case 2:
      
// begin by ensuring we are in the root folder
// the load button is intended to be pressed after input such as sample name, step size, etc.
// are input; so when we go to the root folder, we expect these global variables to be present
// we load global variables and strings with the nvar and svar commands
      
setdatafolder root:

// load previously defined global variables and strings with nvar and svar
// this is distinct to creating globals or creating local variables and strings, both
// of which use variable and string, possibly with /g
nvar step, c, i
svar  sam, ref, tl
               	
// the timestep variable is occasionally useful
// and the number of temperatures, defined as npnts is used in most temperature loops
variable/g timestep = 2*step / c
variable/g npnts = itemsinlist(tl)

// kill any open graphs or tables via a user defined function
// this is necessary when killing waves, overwriting folders, etc.
// because an error will pop up if a graph is open while you try to kill a wave of that graph
killallgraphsandtablesNEW()

// and load a function called loadfiles in a loop
// the use of samref is to run the load waves function twice, taking the first string (sam)
// from the list (sam;ref;) the first time, then ref the second time
string/g samref = sam + ";" + ref
for(i=0; i<2; i+=1)
string samrefstr = stringfromlist(i,samref)
loadfilesNEW(samrefstr)
endfor

// pulserefresh loads the raw data into one folder
// and copies the raw data to another folder, pulses
// this allows the user to try different paddings, windowings, etc. from the raw data
// without having to load the raw data from the hard drive each time
// rather, the raw data is loaded into igor and then copied to the pulses folder each time different parameters are chosen
// and calculations proceed from there
pulserefreshNEW()

// i, j, and k are created as globals when the panel is opened (i.e. here)
// these variables are intended to be used as looping variables and are created globally so that they don't 
// have to be defined within other functions	
// also define npnts which is useful in temperature loops

variable/g i, j, k
               	
// also, create constants: (units are all kg, um, ps)

variable/g eps0 = 8.854187817*10^-12 
variable/g ii = sqrt(-1)
variable/g c = 299.792458
variable/g z0 = 376.73031

                     break
        endswitch
        return 0
end


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

Function slabcalcbuttonproc(ba) : ButtonControl

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


        STRUCT WMButtonAction &ba 
        switch(ba.eventCode)
        
 // case 2 tells you that the mouse click has been released (click, then let go)
        
               case 2:
               
// ensure we're in root, then load all the globals necessary to calculate
// optical constants, make graphs, etc.
               
                     setdatafolder root:
                     svar sam, ref, tl, sampath, refpath, tempticks, cbstr, notabletemps, notablecolor
                     nvar step, dec, ssp, sep, rsp,rep,nsp, nep,lfreq,ufreq, c, notablethickness, notabledash

// clear all graphs and tables
killallgraphsandtablesNEW()

// load raw data                    
pulserefreshNEW()

svar sam, ref
nvar i

// set the temp ticks with the first and last temperatures; this suffices as
// temporary color bar scaling until more specific temperatures are input
tempticksNEW(tempticks)

// calculate the transmission
trNEW()

// calculate the index for slab geometry and from it other optical constants
// slab index goes before trtemperatureNEW so that the left and right endpoints used for the
// index wave are imported to the trtemperature function
slabindexNEW()

// make waves of the transmission versus temperature
trtemperatureNEW()


slabepsNEW()
slabsigmaNEW()
powerabsNEW()
// calculate the optical constants vs temperature
sigmavstempNEW()
epsvstempNEW()
indexvstempNEW()
alphavstemp()
// create image plots
// a new technique was used here to make image plots for epsilon, index, and sigma in one function, more concisely
ratiotolowtempNEW()
ratiotolowtemptrNEW()
imagegraphmaker()


// visible graphs are now hidden
hideallgraphsandtablesNEW()

// the string containing the visible graphs must be set manually to contain nothing by using ""
// otherwise, igor thinks it is null and has trouble working with it
// this loop sees which boxes for graphs are checked and makes those graphs visible


string visiblegraphs = ""
variable numctrls = itemsinlist(cbstr,",")

for(i=0;i<numctrls;i+=1)

	// get info about each control
	controlinfo/w=calculatepanel $stringfromlist(i,cbstr,",")

	if(v_value==1) // (if it's checked)

	// create the string name of the graph by changing "ctrl" to "graph"
	// the naming conventions were chosen to allow this
	string bringtofrontstr = replacestring("ctrl",stringfromlist(i,cbstr,","),"graph")
	
	// and also make the graph visible; this way tilewindows can see it and tile it
	// the $ is used to access the specific string that is represented by the reference "bringtofrontstr"
	dowindow/hide=0/f $bringtofrontstr
	
	// add it to the list of visible graphs
	visiblegraphs=addlistitem(bringtofrontstr,visiblegraphs,",")

	endif
	
endfor

// alphabetize graph list so they graphs are displayed in a consistent order
visiblegraphs = sortlist(visiblegraphs,",",16)

// this if statement prevents igor from tiling procedure and help windows when nothing is checked (i.e. strlen(visiblegraphs)==0).
// the tilewindows command is only available via the command line, so we can make a string and then execute that string in the command line
if(strlen(visiblegraphs)>0)
string tilecmd = "tilewindows/w=(0,0,1063,800) "+visiblegraphs
execute tilecmd
endif

// break is not strictly necessary but recommended and usually necessary
                     break

endswitch
return 0
end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

Function checkboxproc(cba) : checkboxcontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// when checking and unchecking boxes that determine is graphs are shown or not,
// there is no recalculation necessary
// we simply add or remove and then retile

        STRUCT wmcheckboxaction &cba
 
        switch(cba.eventCode)
               case 1:
               case 2:
               
               
               if(cba.checked==1)
               
               // winlist grabs the visible graphs, and returns their names in a string, separated by comma
               // the comma separation is used for compatibility with the tilewindows command
               string visiblegraphs = winlist("*",",","win:1,visible:1")
               
               // we then add the user data from the checked box, which is defined to be the name of the graph to the current list
               // of visible graphs
               visiblegraphs = addlistitem(cba.userdata,visiblegraphs,",")
               
               // alphabetize the list
               // 16 determines the type of ordering; see documentation
               visiblegraphs = sortlist(visiblegraphs,",",16)
               
               // make the graph of the newly checked box visible
               dowindow/hide=0/f $cba.userdata
               
               // retile these graphs
               string tilecmd = "tilewindows/w=(0,0,1063,800) "+visiblegraphs
               execute tilecmd             


               elseif(cba.checked==0)
               
               // we remove the graph of the unchecked box
               visiblegraphs = removefromlist(cba.userdata,winlist("*",",","win:1,visible:1"),",")
               
               if(strlen(visiblegraphs)==0)
               // if it was the last visible graph, we only hide it; we do not retile anything as this can cause procedure, panel and help
               // windows to be retiled and is annoying
               dowindow/hide=1/f $cba.userdata
               
               else
               // else we hide and then retile
               visiblegraphs = removefromlist(cba.userdata,winlist("*",",","win:1,visible:1"),",")
               visiblegraphs = sortlist(visiblegraphs,",",16)
               dowindow/hide=1 $cba.userdata
               tilecmd = "tilewindows/w=(0,0,1063,800) "+visiblegraphs
               execute tilecmd
               endif
               endif
                     break
        endswitch
        return 0
end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function logycheckboxproc(cba) : checkboxcontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// the logy systemization is as follows:
// when the panel is opened, all logy variables are set to 0 via setlogyvartozero()
// when a logy box is checked, the global in the globals:logy folder is updated to reflect the current
// state of the checkbox (1 for checked, 0 for unchecked)
// the recalculateandregraph function is then run after any checking or unchecking of a logy checkbox
// before each graph is created in the recalculateandregraph function, the globals:logy folder is accessed and the value of that graphs logy variable 
// is loaded and fed to the graphing function

// i think this procedure can be shortened and the runtime reduced using a loop, $ taking strings to variables, and maybe only rerunning of the 
// relevant optical constant calculation and regraph; but it works for now.

        struct wmcheckboxaction &cba
 
        switch(cba.eventCode)
               case 1:
               case 2:
               setdatafolder root:
               svar logypath
		 
               nvar i
               setdatafolder $logypath
               
                            
                  if(cmpstr(cba.ctrlname,"sigma2freqlogyctrl")==0)
                  variable/g sigma2freqlogy = cba.checked
                
                  elseif(cmpstr(cba.ctrlname,"eps2freqlogyctrl")==0)
                  variable/g eps2freqlogy= cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"eps2templogyctrl")==0)
                  variable/g eps2templogy = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"nfreqlogyctrl")==0)
                  variable/g nfreqlogy = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"kfreqlogyctrl")==0)
                  variable/g kfreqlogy = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"ntemplogyctrl")==0)
                  variable/g ntemplogy = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"ktemplogyctrl")==0)
                  variable/g ktemplogy = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"sigma1freqlogyctrl")==0)
                  variable/g sigma1freqlogy = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"eps1freqlogyctrl")==0)
                  variable/g eps1freqlogy = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"sigma1templogyctrl")==0)
                  variable/g sigma1templogy = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"eps1templogyctrl")==0)
                  variable/g eps1templogy = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"trfreqlogyctrl")==0)
                  variable/g trfreqlogy = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"trtemplogyctrl")==0)
                  variable/g trtemplogy = cba.checked
                  
                   elseif(cmpstr(cba.ctrlname,"sigma2templogyctrl")==0)
                  variable/g sigma2templogy = cba.checked     
                  
                  elseif(cmpstr(cba.ctrlname,"alphafreqlogyctrl")==0)
                  variable/g alphafreqlogy = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"alphatemplogyctrl")==0)
                  variable/g alphatemplogy = cba.checked             
                  
                  elseif(cmpstr(cba.ctrlname,"eps1imagelogyctrl")==0)
                  variable/g eps1imagelogy = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"eps2imagelogyctrl")==0)
                  variable/g eps2imagelogy = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"nimagelogyctrl")==0)
                  variable/g nimagelogy = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"kimagelogyctrl")==0)
                  variable/g kimagelogy = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"sigma1imagelogyctrl")==0)
                  variable/g sigma1imagelogy = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"sigma2imagelogyctrl")==0)
                  variable/g sigma2imagelogy = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"alphaimagelogyctrl")==0)
                  variable/g alphaimagelogy = cba.checked
                  
                  endif

setdatafolder root:
recalculateandregraphNEW()

                     break
        endswitch
        return 0
end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function logxcheckboxproc(cba) : checkboxcontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// see logycheckboxproc comments

        struct wmcheckboxaction &cba
 
        switch(cba.eventcode)
               case 1:
               case 2:
               setdatafolder root:
               svar logxpath//, logxstr
		 
               nvar i
               setdatafolder $logxpath
               
                            
                  if(cmpstr(cba.ctrlname,"sigma2freqlogxctrl")==0)
                  variable/g sigma2freqlogx = cba.checked
                
                  elseif(cmpstr(cba.ctrlname,"eps2freqlogxctrl")==0)
                  variable/g eps2freqlogx= cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"eps2templogxctrl")==0)
                  variable/g eps2templogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"nfreqlogxctrl")==0)
                  variable/g nfreqlogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"kfreqlogxctrl")==0)
                  variable/g kfreqlogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"ntemplogxctrl")==0)
                  variable/g ntemplogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"ktemplogxctrl")==0)
                  variable/g ktemplogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"sigma1freqlogxctrl")==0)
                  variable/g sigma1freqlogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"eps1freqlogxctrl")==0)
                  variable/g eps1freqlogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"sigma1templogxctrl")==0)
                  variable/g sigma1templogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"eps1templogxctrl")==0)
                  variable/g eps1templogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"trfreqlogxctrl")==0)
                  variable/g trfreqlogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"trtemplogxctrl")==0)
                  variable/g trtemplogx = cba.checked
                  
                   elseif(cmpstr(cba.ctrlname,"sigma2templogxctrl")==0)
                  variable/g sigma2templogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"alphafreqlogxctrl")==0)
                  variable/g alphafreqlogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"alphatemplogxctrl")==0)
                  variable/g alphatemplogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"eps1imagelogxctrl")==0)
                  variable/g eps1imagelogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"eps2imagelogxctrl")==0)
                  variable/g eps2imagelogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"nimagelogxctrl")==0)
                  variable/g nimagelogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"kimagelogxctrl")==0)
                  variable/g kimagelogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"sigma1imagelogxctrl")==0)
                  variable/g sigma1imagelogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"sigma2imagelogxctrl")==0)
                  variable/g sigma2imagelogx = cba.checked
                  
                  elseif(cmpstr(cba.ctrlname,"alphaimagelogxctrl")==0)
                  variable/g alphaimagelogx = cba.checked
                  endif

setdatafolder root:
recalculateandregraphNEW()

            
              // endif
                     break
        endswitch
        return 0
end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function colortablepopupproc(pua) : popupmenucontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// colortablepopupproc updates the imagecolor string 
// which sets the color table of the image plots

struct wmpopupaction &pua

switch(pua.eventcode)

case 2:
setdatafolder root:

if(cmpstr(pua.ctrlname,"imageplotpopup")==0)
string/g imagecolor = pua.popstr
recalculateandregraphNEW()
elseif(cmpstr(pua.ctrlname,"contourpopup")==0)
string/g contourcolor = pua.popstr
endif
recalculateandregraphNEW()
endswitch
end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function linestylepopupproc(pua) : popupmenucontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// this updates the linestyle of the line at the notable temperatures

struct wmpopupaction &pua

switch(pua.eventcode)

case 2:
setdatafolder root:
variable/g notabledash = pua.popnum
recalculateandregraphNEW()
endswitch
end


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


function colorpopupproc(pua) : popupmenucontrol

// this updates the color of the line at the notable temperatures

struct wmpopupaction &pua

switch(pua.eventcode)

case 2:
setdatafolder root:
string/g notablecolor = pua.popstr
notablecolor = notablecolor[1,inf]
notablecolor = removeending(notablecolor)
recalculateandregraphNEW()
endswitch
end


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function tempticksproc(sva) : setvariablecontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

        struct wmsetvariableaction &sva
 
        switch(sva.eventcode)
               case 1:
               case 2:
               case 3:  

// set the value of the global string tempticks to be 
// the user-input one from the panel

setdatafolder root:
string/g tempticks = sva.sval

// make the label and position waves to be used for temperature ticks

tempticksNEW(tempticks)
recalculateandregraphNEW()
		
                     break
        endswitch
        return 0
end
 

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function setlogyvartozero()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// this and the logx version are run at the beginning to initialize all log variables to 0

svar logypath
setdatafolder $logypath

variable/g sigma2freqlogy = 0
variable/g  eps2freqlogy= 0
variable/g eps2templogy = 0
variable/g nfreqlogy = 0
variable/g kfreqlogy = 0
variable/g ntemplogy = 0
variable/g ktemplogy = 0
variable/g sigma1freqlogy = 0
variable/g eps1freqlogy = 0
variable/g sigma1templogy = 0
variable/g eps1templogy = 0
variable/g trfreqlogy = 0
variable/g trtemplogy = 0
variable/g sigma2templogy = 0
variable/g alphafreqlogy = 0
variable/g alphatemplogy = 0

variable/g eps1imagelogy = 0
variable/g eps2imagelogy = 0
variable/g nimagelogy = 0
variable/g kimagelogy = 0
variable/g sigma1imagelogy = 0
variable/g sigma2imagelogy = 0
variable/g alphaimagelogy = 0

variable/g eps1tempratioimagelogy=0
variable/g sigma1tempratioimagelogy=0
variable/g trmagtempratioimagelogy=0

variable/g losstanvsfreqimagelogy = 0

setdatafolder root:

end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function setlogxvartozero()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

svar logxpath
setdatafolder $logxpath

variable/g sigma2freqlogx = 0
variable/g  eps2freqlogx= 0
variable/g eps2templogx = 0
variable/g nfreqlogx = 0
variable/g kfreqlogx = 0
variable/g ntemplogx = 0
variable/g ktemplogx = 0
variable/g sigma1freqlogx = 0
variable/g eps1freqlogx = 0
variable/g sigma1templogx = 0
variable/g eps1templogx = 0
variable/g trfreqlogx = 0
variable/g trtemplogx = 0
variable/g sigma2templogx = 0
variable/g alphaimagelogx = 0
variable/g alphafreqlogx = 0
variable/g alphatemplogx = 0

variable/g eps1imagelogx = 0
variable/g eps2imagelogx = 0
variable/g nimagelogx = 0
variable/g kimagelogx = 0
variable/g sigma1imagelogx = 0
variable/g sigma2imagelogx = 0

variable/g eps1tempratioimagelogx=0
variable/g sigma1tempratioimagelogx=0
variable/g trmagtempratioimagelogx=0

variable/g losstanvsfreqimagelogx = 0


setdatafolder root:

end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function vrhproc(cba) : checkboxcontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// see logycheckboxproc comments

        struct wmcheckboxaction &cba
 
        switch(cba.eventcode)
               case 1:
               case 2:
               setdatafolder root:
               svar sigma1temppath, sigma1path, sam, tl, notabletemps, notablecolor
               nvar npnts, i, freqs, alfreq,aufreq, notabledash, notablethickness, k
               
string sigma1vrhgraphstr = "sigma1vrhgraph"
windowkill("sigma1vrhgraph")
               
if(cba.checked==1)               
               
               if(cmpstr(cba.ctrlname,"vrhdim1ctrl")==0)
               variable exponent = 1/2
               string exponentstr = "1/2"
               string dimstr = "1D"
	        
	        elseif(cmpstr(cba.ctrlname,"vrhdim2ctrl")==0)
               exponent = 1/3
               exponentstr = "1/3"
               dimstr = "2D"
               
                elseif(cmpstr(cba.ctrlname,"vrhdim3ctrl")==0)
               exponent = 1/4
               exponentstr = "1/4"
               dimstr = "3D"
                             
                 endif 
                 
               vrhNEW(exponent)
               wave tempvrh
               
// this is kind of inconvenient, but since the suffixes of the waves are these random
// frequencies, i'm just recreating the waves                
               
setdatafolder $sigma1path

// these definitions have zeros because we're choosing a basis wave
// to extract other numbers from 
string sigma1str0 = sam+"_sigma1_"+stringfromlist(0,tl)
wave sigma1wave0 = $sigma1str0
variable freqstep = (aufreq-alfreq) / freqs

display/n=$sigma1vrhgraphstr/hide=0

setdatafolder $sigma1temppath
for(i=0;i<freqs;i+=1)

string sigma1tempstr = sam+"_sigma1temp_"+num2str((leftx(sigma1wave0)+i*freqstep)*1000)
wave sigma1tempwave = $sigma1tempstr
appendtograph/w=sigma1vrhgraph sigma1tempwave vs tempvrh
endfor

string bottomaxisstr = "\Z14Temperature"+"\S-"+exponentstr+"\M [K]"

label left "\Z18\F'symbol's\F'geneva'\B1\M [(\F'symbol'W\F'geneva' cm)\S-1\M]"
Label bottom bottomaxisstr
ModifyGraph mirror=2
ModifyGraph fSize(left)=14
ModifyGraph fSize(bottom)=14
ModifyGraph lsize=2
SetAxis/A=2 left
ApplyColorTableToTopGraphNEW("Rainbow")
modifygraph grid(left)=1,fstyle=1,axthick=1.5
ColorScale/C/N=text1/F=0/A=LT vert=0, ctab={0,1,Rainbow,0};DelayUpdate
ColorScale/C/N=text1 "Frequency [THz]"
ColorScale/C/N=text1 widthPct=35,heightPct=3
ColorScale/C/N=text1/B=1
TextBox/C/N=text0/F=0/A=RT "\\Z22"+dimstr+" VRH"
TextBox/C/N=text0/B=1
//TextBox/C/N=text0/A=LC/X=3.32/Y=1.39
ColorScale/C/N=text1  ctab={alfreq,aufreq,Rainbow,0}
ModifyGraph log(left)=1
ModifyGraph mode=4,marker=19,mrkThick=1

variable rvar, gvar, bvar
rvar = str2num(stringfromlist(0, notablecolor,","))
gvar = str2num(stringfromlist(1, notablecolor,","))
bvar = str2num(stringfromlist(2, notablecolor,","))

doupdate
getaxis left
if(strlen(notabletemps)>0)
for(k=0;k<itemsinlist(notabletemps);k+=1)
ShowTools/A arrow
SetDrawEnv ycoord= left;DelayUpdate
SetDrawEnv xcoord= bottom
SetDrawEnv linefgc=(rvar,gvar,bvar),dash= notabledash-1,linethick= notablethickness
DrawLine (str2num(stringfromlist(k,notabletemps))^-exponent), v_min, (str2num(stringfromlist(k,notabletemps))^-exponent), v_max
hidetools/a
endfor
endif



setdatafolder root:

//else
//dowindow/f sigma1vrhgraph  
//if (V_flag == 1)
//display/n=$sigma1vrhgraphstr/hide=1  
//endif
endif

                     break
        endswitch
        return 0
end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function ratiotolowtempproc(cba) : checkboxcontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// see logycheckboxproc comments

        struct wmcheckboxaction &cba
 
        switch(cba.eventcode)
               case 1:
               case 2:
               
if(cba.checked==1)               
ratiotolowtempNEW()
ratiotolowtemptrNEW()
imagegraphmaker()
dowindow/hide=0 eps1tempratioimagegraph
dowindow/hide=0 sigma1tempratioimagegraph
dowindow/hide=0 trmagtempratioimagegraph
		 string tilecmd = "tilewindows/w=(0,0,1063,800) "+"eps1tempratiograph,sigma1tempratiograph, trmagtempratiograph,eps1tempratioimagegraph,sigma1tempratioimagegraph,trmagtempratioimagegraph"
		 execute tilecmd
else

dowindow/hide=1 eps1tempratiograph
dowindow/hide=1 sigma1tempratiograph
dowindow/hide=1 trmagtempratiograph
dowindow/hide=1 eps1tempratioimagegraph
dowindow/hide=1 sigma1tempratioimagegraph
dowindow/hide=1 trmagtempratioimagegraph

endif

                     break
        endswitch
        return 0
end


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

Function losstanproc(cba) : checkboxcontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// when checking and unchecking boxes that determine is graphs are shown or not,
// there is no recalculation necessary
// we simply add or remove and then retile

        STRUCT wmcheckboxaction &cba
 
        switch(cba.eventCode)
               case 1:
               case 2:
               
               if(cba.checked==1)
               
auxconstantsNEW()
auximagemaker()
dowindow/hide=0 losstanfreqgraph
dowindow/hide=0 losstanvsfreqimagegraph
 string tilecmd = "tilewindows/w=(0,0,1063,800) "+"losstanfreqgraph,losstanvsfreqimagegraph"
 execute tilecmd
 elseif(cba.checked==0)
dowindow/hide=1 losstanfreqgraph
dowindow/hide=1losstanvsfreqimagegraph
endif
                     break
        endswitch
        return 0
end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

Function zsproc(cba) : checkboxcontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// when checking and unchecking boxes that determine is graphs are shown or not,
// there is no recalculation necessary
// we simply add or remove and then retile

        STRUCT wmcheckboxaction &cba
 
        switch(cba.eventCode)
               case 1:
               case 2:
               
               if(cba.checked==1)
               
auxconstantsNEW()
auximagemaker()
dowindow/hide=0 xsvsfreqgraph
dowindow/hide=0 rsvsfreqgraph
dowindow/hide=0 xsvsfreqimagegraph
dowindow/hide=0 rsvsfreqimagegraph
 string tilecmd = "tilewindows/w=(0,0,1063,800) "+"xsvsfreqgraph,rsvsfreqgraph, xsvsfreqimagegraph, rsvsfreqimagegraph"
 execute tilecmd
 elseif(cba.checked==0)
dowindow/hide=1 xsvsfreqgraph
dowindow/hide=1 rsvsfreqgraph
dowindow/hide=1 xsvsfreqimagegraph
dowindow/hide=1 rsvsfreqimagegraph
endif
                     break
        endswitch
        return 0
end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

Function rproc(cba) : checkboxcontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// when checking and unchecking boxes that determine is graphs are shown or not,
// there is no recalculation necessary
// we simply add or remove and then retile

        STRUCT wmcheckboxaction &cba
 
        switch(cba.eventCode)
               case 1:
               case 2:
               
               if(cba.checked==1)
               
auxconstantsNEW()
auximagemaker()
dowindow/hide=0 rmagvsfreqgraph
dowindow/hide=0 rphasevsfreqgraph
dowindow/hide=0 rmagvsfreqimagegraph
dowindow/hide=0 rphasevsfreqimagegraph
 string tilecmd = "tilewindows/w=(0,0,1063,800) "+"rmagvsfreqgraph,rphasevsfreqgraph, rmagvsfreqimagegraph, rphasevsfreqimagegraph"
 execute tilecmd
 elseif(cba.checked==0)
dowindow/hide=1 rmagvsfreqgraph
dowindow/hide=1 rphasevsfreqgraph
dowindow/hide=1 rmagvsfreqimagegraph
dowindow/hide=1 rphasevsfreqimagegraph
endif
                     break
        endswitch
        return 0
end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

Function energylossproc(cba) : checkboxcontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// when checking and unchecking boxes that determine is graphs are shown or not,
// there is no recalculation necessary
// we simply add or remove and then retile

        STRUCT wmcheckboxaction &cba
 
        switch(cba.eventCode)
               case 1:
               case 2:
               
               if(cba.checked==1)
               
auxconstantsNEW()
auximagemaker()
dowindow/hide=0 energylossrevsfreqgraph
dowindow/hide=0 energylossimvsfreqgraph
dowindow/hide=0 energylossrevsfreqimagegraph
dowindow/hide=0 energylossimvsfreqimagegraph
 string tilecmd = "tilewindows/w=(0,0,1063,800) "+"energylossrevsfreqgraph,energylossimvsfreqgraph,energylossrevsfreqimagegraph, energylossimvsfreqimagegraph"
 execute tilecmd
 elseif(cba.checked==0)
dowindow/hide=1 energylossrevsfreqgraph
dowindow/hide=1 energylossimvsfreqgraph
dowindow/hide=1 energylossrevsfreqimagegraph
dowindow/hide=1 energylossimvsfreqimagegraph
endif
                     break
        endswitch
        return 0
end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

Function ebphaseproc(cba) : checkboxcontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// when checking and unchecking boxes that determine is graphs are shown or not,
// there is no recalculation necessary
// we simply add or remove and then retile

        STRUCT wmcheckboxaction &cba
 
        switch(cba.eventCode)
               case 1:
               case 2:
               
               if(cba.checked==1)
               
auxconstantsNEW()
auximagemaker()
dowindow/hide=0 ebphasevsfreqgraph
dowindow/hide=0 ebphasevsfreqimagegraph
 string tilecmd = "tilewindows/w=(0,0,1063,800) "+"ebphasevsfreqgraph, ebphasevsfreqimagegraph"
 execute tilecmd
 elseif(cba.checked==0)
dowindow/hide=1 ebphasevsfreqgraph
dowindow/hide=1 ebphasevsfreqimagegraph

endif
                     break
        endswitch
        return 0
end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 



function bringgraphstofrontbuttonproc(ba) : buttoncontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// a handle button if the graphs begin to become invisible

        STRUCT WMButtonAction &ba 
        switch(ba.eventCode)
        
 // case 2 tells you that the mouse click has been released (click, then let go)
        
               case 2:
                     setdatafolder root:
			bringgraphstofrontNEW()

                     break
        endswitch
        return 0
end


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
  
function printtopdfproc(ba) : buttoncontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

struct wmbuttonaction &ba

switch(ba.eventCode)
        
      case 2:
      
setdatafolder root:

reportinput()
//	string notebookname = "THz Summary"
//	variable visorhid = 1, graphsperpara = 1, graphmode = 0, scaling = 100
//AddGraphsToNotebook(notebookname, visorhid, graphsperpara, graphmode, scaling)

                     break
        endswitch
        return 0
end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function derivproc(ba) : buttoncontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// a handle button if the graphs begin to become invisible

        STRUCT WMButtonAction &ba 
        switch(ba.eventCode)
        
 // case 2 tells you that the mouse click has been released (click, then let go)
        
               case 2:
			deriv()
return 0
                     break
        endswitch
        return 0
end





/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function derivhelpproc(ba) : buttoncontrol

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

// a handle button if the graphs begin to become invisible

        STRUCT WMButtonAction &ba 
        switch(ba.eventCode)
        
 // case 2 tells you that the mouse click has been released (click, then let go)
        
               case 2:
 
               
	String nb = "DifferentiateSomethingHelp"
	windowkill("DifferentiateSomethingHelp")
	NewNotebook/N=$nb/F=1/V=1/K=1/W=(355,234,855,529)
	notebook $nb, writeprotect=1
	Notebook $nb defaultTab=36, statusWidth=252
	Notebook $nb showRuler=1, rulerUnits=1, updating={1, 1}
	Notebook $nb newRuler=Normal, justification=0, margins={0,0,468}, spacing={0,0,0}, tabs={}, rulerDefaults={"Geneva",10,0,(0,0,0)}
	Notebook $nb ruler=Normal, text="How to use the \"Differentiate Something\" button:\r"
	Notebook $nb text="\r"
	Notebook $nb text="1. Within the data browser, drag the red arrow to the folder which contains the waves you want to differ"
	Notebook $nb text="entiate. \r"
	Notebook $nb text="\r"
	Notebook $nb text="2. THEN click the \"Differentiate Something\" button. \r"
	Notebook $nb text="\r"
	Notebook $nb text="3. Enter the base name when prompted. i.e., if a typical wave is greatwave_n_15, the base name would be "
	Notebook $nb text="greatwave_n (do not include the underscore preceding the temperature).\r"
	Notebook $nb text="\r"
	Notebook $nb text="4. Derivatives will be shown in a plot, and will be contained in the derivs folder. \r"
	Notebook $nb text="\r"
	Notebook $nb text="5. For multiple derivatives, copy the waves you just differentiated to a new folder, and repeat the proc"
	Notebook $nb text="ess, including the preceding \"d_\" in the new base name.\r"
return 0
                     break
        endswitch
        return 0
end



