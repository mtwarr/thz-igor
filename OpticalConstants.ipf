#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <WMBatchCurveFitIM>
#include <Global Fit 2>
#include <KBColorizeTraces>

// this software was funded by the Center for Emergent Materials at The Ohio State University
// NSF grant DMR-1420451 and posted with permission

// conventional practices:
// add a setdatafolder root: line to beginning and end of all functions (once input parameters have been declared)
// include colon at the end of strings specifying folder locations
// color bar scales go from 0 to 1 in terms of position

// the slab function is one of two start points, the other being film. these are for the
// two geometries we use most. 



/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function loadfilesNEWNEW(sam0, ref0, tl0,step0,ssp0, sep0, rsp0,rep0)

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



// loadfilesNEWNEW loads the raw data, and compiles associated information
// the data entry here is decoupled from the rest of the analysis since this is a 
// once-and-for-all load (for a given set of temperatures)

// sam: sample name (string)
// ref: ref name (string)
// tl: temperature list (string)
// stepsize: delay stage increments (don't multiply by two for incoming and outgoing shift) (variable)
// sysstartpt (etc.) = delay stage readouts for these values (variable)

string sam0, ref0, tl0
variable step0, ssp0, sep0,rsp0,rep0

string/g sam=sam0, ref=ref0, tl=tl0
variable/g step=step0, ssp=ssp0, sep=sep0, rsp=rsp0,rep=rep0, i, j, k
variable/g timestep = 2*step / 299.792  

// npnts is the number of temps scanned
variable/g npnts = itemsinlist(tl)

string/g samref = sam + ";" + ref

// kill any open graphs or tables
killallgraphsandtablesNEW()

// and load a function called loadfiles in a loop

for(i=0; i<2; i+=1)
string samrefstr = stringfromlist(i,samref)
loadfilesNEW(samrefstr)
endfor

// make a global temperature wave
// use this program; it has an if statement to convert
// non-decimal temperatures, like 7p6 to 7.6
temperaturewaveNEW()
 
pulserefreshNEW()
 
// now that the files are loaded, we can get to the physics //

end




/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function slabNEW(dec0,nsp0, nep0,d0, n0,lfreq0,ufreq0,freqs0,tempticks0)

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\





// the 0s are just naming convention, to free up the actual names for global use in the code

// decimal: percent from edges that exponentially decays (variable)
// newstartpt (etc.) = endpoints for padded wave (variable)
// d: slab thickness (variable)
// n: guess for index (variable)
// lowerusablefreq (etc.): useful if divergent solutions at very low frequency impede solution at more stable frequencies (variable)
// freqs: the number of frequencies between high and low frequency points to calculate

// define variable types

variable dec0, nsp0, nep0, d0, n0, lfreq0, ufreq0, freqs0
string tempticks0

// pulserefreshNEW() reloads raw pulses into the folder root:pulses:*:
// so that fresh pulses are used for each computation with adjusted parameters

pulserefreshNEW()

// the input is declared global by a separate function, f1
// (see functions after this main body; they're arranged alphabetically)
// this allows the functions that call each other to not need to load variables and strings

slabglobalNEW(dec0, nsp0, nep0, d0, n0, lfreq0, ufreq0, freqs0,tempticks0)




// you can return -1, which exits the function, here to see that these strings and variables are created in the data browser
// you might have to check the boxes from strings, variables, etc.
// erase the backspaces below to get to the command
// return -1
 
// here is a trick to load the files from a given folder:
// define a string called samref and a variable i

// need to load these global variables first

// kill any open graphs or tables
killallgraphsandtablesNEW()

svar sam, ref, tempticks

// calculate the transmission

tempticksNEW(tempticks)

trNEW()

// make waves of the transmission versus temperature
// parameter is how many frequencies to plot versus temperature
trtemperatureNEW()

// calculate the index for slab geometry
slabindexNEW()
slabepsNEW()
slabsigmaNEW()
sigmavstempNEW()
epsvstempNEW()
indexvstempNEW()

//string visiblegraphlist = "sigma1graph,eps1graph, sigma1tempgraph,eps1tempgraph, trgraph"
//string dowindowcmd = "dowindow/f sigma1graph"
//string tilecmd = "tilewindows "+visiblegraphlist
//execute dowindowcmd
//execute tilecmd

// hideallgraphsandtablesNEW()



end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function pulserefreshNEW()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


// pulserefresh reloads pulses from the raw data
// before making calculations with them
// the raw data should not be adjusted after the first loading
// custom parameter adjustment (padding for instance) will refresh the 
// pulses before operating

setdatafolder root: 
svar sam, ref, tl, samref
nvar npnts, i, j
	            
// kill pulsegraph so that we may overwrites folders and waves	            
	            
windowkill("pulsegraph")

// since we are copying from the raw data, we need to create the
// folder where the calculational procedures will draw the pulses from.
// this is important so that whenever we calculate, we are starting from raw pulses,
// not from a previously adjusted pulse.
        
string/g sampath = "root:pulses:"+sam+":"
string/g refpath = "root:pulses:"+ref+":"                                                                                                                                        
overwritefolder(sampath)
overwritefolder(refpath)

// use of samref loop:
// samref is a list that has the sam string as the first entry and
// the ref string as the second; this general way of writing loops
// allows the code to be more concise and makes it so only one interior loop
// needs to be written for both sam and ref cases
// note that samref is only used to take the list entry and use it to build up names                                                                

for(i=0;i<itemsinlist(samref);i+=1)                     
for(j=0;j<npnts;j+=1)

setdatafolder $"root:raw_pulses:"+stringfromlist(i,samref)+":"
string pulsestr = stringfromlist(i,samref)+"_"+stringfromlist(j,tl)
wave pulsewave =$pulsestr
// setscale/P x 0,timestep,"", pulsewave
duplicate/o pulsewave $pulsestr+"_copy"
movewave $pulsestr+"_copy", $"root:pulses:"+stringfromlist(i,samref)+":"
setdatafolder $"root:pulses:"+stringfromlist(i,samref)+":"
rename $pulsestr+"_copy", $pulsestr

endfor
endfor

setdatafolder root:

end




/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function slabglobalNEW(dec0,nsp0, nep0, d0, n0,  lfreq0, ufreq0, freqs0,tempticks0)

// defines global parameters from user input and physical constants

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 



// identify them to the function

variable dec0, nsp0, nep0, d0, n0, lfreq0,ufreq0, freqs0
string tempticks0

// declare strings, variables, and constants globally  //

// convert user input to global strings
variable/g lfreq=lfreq0,ufreq=ufreq0,dec=dec0,nsp=nsp0,nep=nep0,d=d0,n=n0, freqs = freqs0
string/g tempticks=tempticks0

// constants; units are all kg, um, ps
variable/g eps0 = 8.854187817*10^-12 
variable/g ii = sqrt(-1)
variable/g c = 299.792458
variable/g z0 = 376.73031

end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function temperaturewaveNEW()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


setdatafolder root:
svar tl
nvar i, npnts

make/o/n=(npnts) temperature

for(i=0;i<npnts;i+=1)
temperature[i] = str2num(replacestring("p",stringfromlist(i,tl),"."))
endfor

setdatafolder root:
end



/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function trNEW()

// calculates transmission from the pulses
// includes windowing and padding

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 



setdatafolder root:

svar sam, ref, tl, sampath, refpath,logypath, logxpath
nvar step, dec, ssp, sep, rsp,rep,nsp, nep,lfreq,ufreq, c

wave/t temptickslabel
wave temptickspos

// create folders for data organization
// overwrites if a folder already exists
// creates if no folder of the name exists

     	     string/g trpath = "root:tr:"
     	     overwritefolder(trpath)   
     	     
     	     string/g trmagpath = "root:tr:mag:"
     	     overwritefolder(trmagpath)             
     	     
            string/g samfftpath = sampath+"fft:"
	     overwritefolder(samfftpath)
	     
            string/g reffftpath = refpath+"fft:"
	     overwritefolder(reffftpath)

            
// create variable to loop over (i) and number of temperatures (npnts)
// for use in loops

variable i, npnts
npnts = itemsinlist(tl)

// loop through waves, pad, fft

for(i=0;i<npnts;i+=1)

	setdatafolder $sampath
	
// this is the basic structure that loads waves
// create a string that matches the name of a wave in a given data folder
// create a local, temporary reference for that wave
	
	string samstr = sam +"_"+ stringfromlist(i,tl)
	wave samwave = $samstr
	
// window the pulse

	trwindowerNEW(samwave,dec)
	
// do the same for the reference pulse	
	
	setdatafolder $refpath
	
	string refstr = ref +"_"+ stringfromlist(i,tl)
	wave refwave = $refstr
	
	trwindowerNEW(refwave,dec)
	
// custom pad the waves so that start and end points 	
	
	wavepadderNEW(samwave, step, ssp, sep,nsp, nep)
	wavepadderNEW(refwave, step, rsp, rep,nsp, nep)
	
// scale these waves based on stepsize	
	
	setscale/p x, 0, 2*step/ 299.792, samwave
	setscale/p x, 0, 2*step/ 299.792, refwave

	setdatafolder $trpath

// create local references for fft's of sam and ref waves
// though the same names would be created anyways, this allows easy
// division to create the transmission wave
	
	string samfftstr = nameofwave(samwave)+"_fft"
	make/c/o $samfftstr/wave= samfftwave
	
	string reffftstr = nameofwave(refwave)+"_fft"
	make/c/o $reffftstr/wave=reffftwave
	
	fft/out=1/dest=samfftwave samwave
	fft/out=1/dest=reffftwave refwave
	
	string trstr = sam+"_tr_" +stringfromlist(i,tl)
	make/c/o/n=(numpnts(samfftwave)) $trstr/wave=trwave

	trwave = samfftwave / reffftwave
	
	setscale/P x 0,dimdelta(samfftwave,0),"",trwave

setdatafolder $trmagpath
string trmagstr = sam+"_trmag_"+stringfromlist(i,tl)
make/o/n=(numpnts(trwave)) $trmagstr/wave=trmagwave
trmagwave = cabs(trwave)
setscale/P x 0,dimdelta(samfftwave,0),"",trmagwave
setdatafolder $trpath
	
	
// magnetic susceptibility test
// -ln (T(omega)) is supposed to be proportional to 
// the magnetic susceptibility 

// string magsuscstr0 = "magsusc_"+stringfromlist(i,tl)
// make/c/o/n=(numpnts(trwave)) $magsuscstr0/wave= magsuscwave0

// setscale/P x 0,dimdelta(trwave,0),"",magsuscwave0
// magsuscwave0 = -ln(trwave) / (2*pi*x)
	
	
		
// the following command utilizes the dfref system to
// move waves

dfref samfftdfref = $samfftpath
dfref reffftdfref = $reffftpath
dfref trdfref = $trpath

movewave samfftwave, samfftdfref
movewave reffftwave, reffftdfref
movewave trwave, trdfref	
	
endfor

string trbefore = sam+"_tr_"

setdatafolder $logypath
nvar trfreqlogy
setdatafolder $logxpath
nvar trfreqlogx
setdatafolder root:

graphervsfreq("trfreqgraph",trpath, trbefore,"\Z14Transmission Magnitude","\\Z30 T","ColdWarm",3,trfreqlogy,trfreqlogx)

setdatafolder $"root:"

end



/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function trwindowerNEW(wwave, dec)

// windows the raw data

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


// only have to load one wave
// decimal is fraction of wave to window over at the edges

wave wwave
variable dec

variable length = numpnts(wwave)

	variable edge = dec * length

	make/o/n=(length) windowfn	

	
variable j
	
		for(j=0;j<numpnts(windowfn);j+=1)

			if(j<edge)
			windowfn[j] = .5*(1-cos(2 * pi * j / (edge*2) ))
// this uses the n = 0, 1, 2,...,N-1 version			
		
		
			elseif(edge<=j&&j<=(length - edge))
			windowfn[j] = 1
			
			else
// this uses a modified version of the if case, to allow counting from the end of the wave
// length - 1 is used because of zero indexing			
			windowfn[j] = .5*(1-cos(2 *pi * (length -j) / (edge*2)))

			endif
			
endfor

wwave = wwave * windowfn


end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function wavepadderNEW(systemwave, stepsize, systemstartpoint, systemendpoint, newstartpoint, newendpoint)

// pads windowed data

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


wave systemwave
variable stepsize, systemstartpoint, systemendpoint, newstartpoint, newendpoint

insertpoints 0, (systemstartpoint - newstartpoint)/stepsize,  systemwave
insertpoints (systemendpoint-systemstartpoint)/stepsize + (systemstartpoint-newstartpoint)/stepsize + 1, (newendpoint - systemendpoint)/stepsize, systemwave

end




 /////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function slabindexNEW()

// solve index of refraction from raw transmission data in slab geoemetry
// create plots

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 







// root will usually be set at the beginning of functions to ensure
// global strings are loaded correctly

setdatafolder root:

	svar sam, tl, trpath, logypath, logxpath
	nvar d, n, lfreq, ufreq
	
variable i,npnts =itemsinlist(tl)

string/g npath = "root:index:n:"
overwritefolder(npath)

string/g kpath = "root:index:k:"
overwritefolder(kpath)


setdatafolder trpath
	 
for(i=0;i<npnts;i+=1)
	
	string trstr = sam +"_tr_"+ stringfromlist(i,tl)
	wave/c trwave = $trstr

	variable fstep
	fstep = dimdelta(trwave,0)
	
// set up upper and lower limits of frequency to set usable limit
// this is because the equations might not be solvable in noisy region

// lindex: lower index
// uindex: upper index

	variable lindex, uindex

// ceiling and floor are used to choose
// frequencies above and below this frequency

	lindex = ceil(lfreq / fstep)
	uindex = floor(ufreq / fstep)

// alfreqinloop: actual lower freq
// necessary to set starting point on the waves that will be created
// it is called in loop because i will make a alfreq later that is in root
// i just didn't want to have to change the folder within the loop

// don't need any +1 for zero indexing
// if the lower index is 0, the actuallowerfreq will be 0, etc.

	variable alfreqinloop
	alfreqinloop =fstep*lindex

// solved index will only exist frequency window established above
// since equations only solved for these frequencies

// generally have wave names take form sam+"_"+stringfromlist(i,tl)+"additional thing"
// for parsing wave names later
	
	string nstr = sam + "_n_" + stringfromlist(i,tl)
	make/o/n=(uindex-lindex+1) $nstr/wave=nwave
	
	string kstr = sam + "_k_"+ stringfromlist(i,tl)
	make/o/n=(uindex-lindex+1) $kstr/wave=kwave
	
// parameters are: re and im of transmission
// d (thickness) and n (guess for index)
// define starting guesses for solution wave
	
	make/o/n=4 known
	make/o/n=2 soln
	soln[0]=n; soln[1]=0

// this freq will need a factor of 2*pi when used as omega
// p used below will match the position number in the wave (0,1,2,...)

	make/o/n=(numpnts(trwave)) freq
	freq=fstep*p

// ceiling and floor were used to define the limits
// the bounds are inclusive

variable j

// this for loop that solves frequency by frequency 
// is inside the loop that loaded the transmission for a 
// particular frequency
 
	For(j=lindex;j<=uindex; j+=1)
	
// define known parameters	
	
		known[0]=real(trwave[j])
		known[1]=Imag(trwave[j])
		known[2]=freq[j]
		known[3]=d
		
		findroots/q/x=soln/f=1/I=200/t=1e-10 indexreeq, known, indeximeq, known

// since the real and imaginary parts are shorter to account for only the 
// frequencies we are interested in, lowerindex is subtracted away so that the first 
// iteration of the loop assigns the zeroth value of real part
// the scaled start frequency will insure this plots right
// (see set scale below)
		
		nwave[j-lindex]=soln[0]
		kwave[j-lindex]=soln[1]
		
	endfor

	setscale/p x alfreqinloop, fstep,"", nwave, kwave
	
// not sure why movewave is used here without dfref like in
// tr function; but i'm going to leave it since it's been working
	
	movewave nwave, $npath
	movewave kwave, $kpath
	
	killwaves/z soln, known, w_yatroot, freq
	
	endfor


variable k

setdatafolder $logypath
nvar nfreqlogy, kfreqlogy
setdatafolder $logxpath
nvar nfreqlogx, kfreqlogx
setdatafolder root:

string nbefore = sam+"_n_"
graphervsfreq("nfreqgraph",npath, nbefore,"\Z14n","\\Z30 n","ColdWarm",3,nfreqlogy,nfreqlogx)
string kbefore = sam+"_k_"
graphervsfreq("kfreqgraph",kpath, kbefore,"\Z14k","\\Z30k","ColdWarm",3,kfreqlogy,kfreqlogx)

// we will take the left and right x points from an arbitrary
// index wave so that we may set the actual upper and lower frequencies
// by instead utilizing lfreq and ufreq elsewhere in the code, we may try to generate plots
// for frequencies higher than ufreq, which creates indexing issues

setdatafolder $npath
string nstr0 = sam+"_n_"+stringfromlist(0,tl)
wave nwave0 = $nstr0
setdatafolder root:
variable/g alfreq = leftx(nwave0)
variable/g aufreq = rightx(nwave0)-deltax(nwave0)


end



/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function indexreeq(w,n,k)

// this function is called by the function that solves for the index of refraction
// slabindex()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 




// real part of equation of fresnel transmission through slab

// where is the freq. dep in these equations? o
// put in in params in the loop (in indexcalc function)
// A: Re[T], B: Im[T], o: freq[i], d: d
// o is for omega

	wave w
	variable n, k
	variable A,B,o,d
	variable c= 299.792458
	A=w[0]; B=w[1]; o=w[2]; d=w[3]

return A+2*A*n+A*n^2-2*B*k-2*B*n*k-A*k^2-4*exp(-2*pi*o*d*k/c)*(n*cos(2*pi*o*d*(n-1)/c) - k*sin(2*pi*o*d*(n-1)/c))



end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function indeximeq(w,n,k)

// this function is called by the function that solves for the index of refraction
// slabindex()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

	wave w
	variable n, k
	variable A,B,o,d
	variable c= 299.792458
	A=w[0]; B=w[1]; o=w[2]; d=w[3]
	
return B+2*B*n+B*n^2+2*A*k+2*A*n*k-B*k^2-4*exp(-2*pi*o*d*k/c)*(k*cos(2*pi*o*d*(n-1)/c)+n*sin(2*pi*o*d*(n-1)/c))

end



/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function slabepsNEW()

// calculates epsilon from the index of refraction

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


// this tells you what global variables will be loaded	  
	  
	  setdatafolder root:
	  svar sam, tl,npath,kpath,logypath, logxpath
	  nvar lfreq, ufreq, npnts, i
	  
// delete folders if they exist, create blank ones in their place	  
	      
string/g eps1path = "root:eps:one:"
overwritefolder(eps1path)

string/g eps2path = "root:eps:two:"
overwritefolder(eps2path)	                                      
            
// eps calculation from index of refraction calculation            
            
for(i=0;i<npnts;i+=1)

// load index waves

setdatafolder $npath

string nstr = sam +"_n_"+stringfromlist(i,tl)
wave nwave = $nstr	

setdatafolder $kpath
	
string kstr =sam +"_k_"+stringfromlist(i,tl)
wave kwave = $kstr


// frequency waves need to be defined in the determination
// of sigma since the relation between epsilon and sigma
// depends on frequency

variable df = dimdelta(nwave,0)
make/o/n=(numpnts(nwave)) freq
freq = leftx(nwave) + df*p

// set sigma folder
// create wave

setdatafolder $eps1path	
	
string eps1str = sam +"_eps1_"+stringfromlist(i,tl)
make/o/n=(numpnts(nwave)) $eps1str/wave= eps1wave
	
// epsilon_0 * THz gives 8.854*10^-12 * 10^12	= 8.854
// converting to cgs goes like: sigma_cgs --> sigma_cgs / 4 pi e0 --> n k f * 4 pi * e0 = n k f 4 pi 8.854 where f ~ 1
// divide by 100 to get inverse ohm cm
	
eps1wave = nwave^2-kwave^2
setscale/P x leftx(nwave),dimdelta(nwave,0),"", eps1wave
	
setdatafolder $eps2path	
	string eps2str = sam +"_eps2_"+stringfromlist(i,tl)
	make/o/n=(numpnts(nwave)) $eps2str/wave= eps2wave
	eps2wave = 2*nwave*kwave
	setscale/P x leftx(nwave),dimdelta(nwave,0),"", eps2wave

endfor

setdatafolder $npath
killwaves freq
setdatafolder $"root:"

setdatafolder $logypath
nvar eps1freqlogy, eps2freqlogy
setdatafolder $logxpath
nvar eps1freqlogx, eps2freqlogx
setdatafolder root:

string eps1before = sam+"_eps1_"
graphervsfreq("eps1freqgraph",eps1path,eps1before,"\Z18\F'symbol'e\F'geneva'\B1\M ","\Z30 \F'Symbol'e\Z30\B1","ColdWarm",3,eps1freqlogy,eps1freqlogx)
string eps2before = sam+"_eps2_"
graphervsfreq("eps2freqgraph",eps2path,eps2before,"\Z18\F'symbol'e\F'geneva'\B2\M ","\Z30 \F'Symbol'e\Z30\B2","ColdWarm",3,eps2freqlogy,eps2freqlogx)

setdatafolder $"root:"


end



/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function slabsigmaNEW()

// this calculates sigma based on an index of refraction

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


// this tells you what global variables will be loaded	  
	  
	  setdatafolder root:
	  svar sam, tl,npath,kpath,logypath, logxpath
	  nvar lfreq, ufreq, npnts, i
	  
// delete folders if they exist, create blank ones in their place	  
	      
string/g sigma1path = "root:sigma:one:"
overwritefolder(sigma1path)

string/g sigma2path = "root:sigma:two:"
overwritefolder(sigma2path)	                                      
            
// sigma calculation from index of refraction calculation            
            
for(i=0;i<npnts;i+=1)

// load index waves

setdatafolder $npath

string nstr = sam +"_n_"+stringfromlist(i,tl)
wave nwave = $nstr	

setdatafolder $kpath
	
string kstr =sam +"_k_"+stringfromlist(i,tl)
wave kwave = $kstr


// frequency waves need to be defined in the determination
// of sigma since the relation between epsilon and sigma
// depends on frequency

// set data folder to root: so that freq is created there
// this is done in lieu of movewaves since movewaves does not overwrite

variable df = dimdelta(nwave,0)
setdatafolder root:
make/o/n=(numpnts(nwave)) freq
freq = leftx(nwave) + df*p

// set sigma folder
// create wave

setdatafolder $sigma1path	
	
string sigma1str = sam +"_sigma1_"+stringfromlist(i,tl)
make/o/n=(numpnts(nwave)) $sigma1str/wave= sigma1wave
	
// epsilon_0 * THz gives 8.854*10^-12 * 10^12	= 8.854
// converting to cgs goes like: sigma_cgs --> sigma_cgs / 4 pi e0 --> n k f * 4 pi * e0 = n k f 4 pi 8.854 where f ~ 1
// divide by 100 to get inverse ohm cm
	
sigma1wave = nwave*kwave*(freq)*8.854*4*pi / 100
setscale/P x leftx(nwave),dimdelta(nwave,0),"", sigma1wave
	
setdatafolder $sigma2path	
	string sigma2str = sam +"_sigma2_"+stringfromlist(i,tl)
	make/o/n=(numpnts(nwave)) $sigma2str/wave= sigma2wave
	sigma2wave = (1 - (nwave^2 - kwave^2))*(freq)*2*pi*8.854 / 100
	setscale/P x leftx(nwave),dimdelta(nwave,0),"", sigma2wave

endfor

setdatafolder $npath
movewave freq, root:
setdatafolder $"root:"

setdatafolder $logypath
nvar sigma1freqlogy, sigma2freqlogy
setdatafolder $logxpath
nvar sigma1freqlogx,sigma2freqlogx
setdatafolder root:

// the "before" strings are local, not global

string sigma1before = sam+"_sigma1_"
graphervsfreq("sigma1freqgraph",sigma1path,sigma1before,"\Z18\F'symbol's\F'geneva'\B1\M [(\F'symbol'W\F'geneva' cm)\S-1\M]","\Z30 \F'Symbol's\Z30\B1","ColdWarm",3,sigma1freqlogy,sigma1freqlogx)
string sigma2before = sam+"_sigma2_"
graphervsfreq("sigma2freqgraph",sigma2path,sigma2before,"\Z18\F'symbol's\F'geneva'\B2\M [(\F'symbol'W\F'geneva' cm)\S-1\M]","\Z30 \F'Symbol's\Z30\B2","ColdWarm",3,sigma2freqlogy,sigma2freqlogx)

setdatafolder $"root:"

end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function trtemperatureNEW()

// calculates the transmission for a given frequency as a function of temperature

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


// enter freq in THz; output will be GHz

// make freqs0 a global variable
setdatafolder root:
svar sam, trpath, tl,logypath,logxpath, npath, trmagpath
nvar aufreq, alfreq, freqs, i, j, k
wave temperature

variable freqstep = (aufreq-alfreq) / freqs

// i loop: for each frequency, create a wave (fntempwave)
// j loop: loop over each temperature, find the transmission wave for that temperature
// and assign the jth point of fntempwave to be the transmission of the ith frequency at 
// temperature j

// temp here refers to "temperature", not "temporary"

windowkill("trtempgraph")

string/g trtemppath = trpath + "vstemp:"
overwritefolder(trtemppath)

// these definitions have zeros because we're choosing a basis wave
// to extract other numbers from 
// need to run after index is calculated since the leftx of the index
// is only calculated after index; trtemperature should be run after for consistency
// in upper and lower freqs in the plots
setdatafolder $npath
string nstr0 = sam+"_n_"+stringfromlist(0,tl)
wave nwave0 = $nstr0


/////////////

//for(i=0;i<freqs;i+=1)
//
//setdatafolder $sigma1temppath
//
//string sigma1tempstr = sam+"_sigma1temp_"+num2str((leftx(sigma1wave0)+i*freqstep)*1000)
//make/o/n=(itemsinlist(tl)) $sigma1tempstr/wave=sigma1tempwave
//
//setdatafolder $sigma1path	
//
//	for(j=0;j<itemsinlist(tl);j+=1)
//
//   		    string sigma1str = sam+"_sigma1_"+stringfromlist(j,tl)
//		    wave sigma1wave = $sigma1str
//		    sigma1tempwave[j]=sigma1wave(freqstep * i + leftx(sigma1wave0))
//       endfor
//
//endfor
//
//setdatafolder $sigma1temppath
//
//for(k=0;k<freqs;k+=1)
//	string sigma1tempstr1 =sam+"_sigma1temp_"+num2str((leftx(sigma1wave0)+k*freqstep)*1000)
//	make/o/n=(itemsinlist(tl)) $sigma1tempstr1/wave=sigma1tempwave1
//	 appendtograph/w=sigma1tempgraph sigma1tempwave1 vs temperature
//endfor


/////////////


//for(i=0;i<freqs;i+=1)
//
//	setdatafolder $trtemppath
//	string trvstempstr = sam+"_trtemp_"+num2str((leftx(nwave0)+i*freqstep)*1000)
//	make/c/o/n=(itemsinlist(tl)) $trvstempstr/wave=trvstempwave
//	setdatafolder $trpath 
//		    	
//		for(j=0;j<itemsinlist(tl);j+=1)
//
//   		    string trstr = sam+"_tr_"+stringfromlist(j,tl)
//		    wave trwave = $trstr
//	
//		    trvstempwave[j]=trwave(freqstep * i + leftx(nwave0))
//	
//             endfor
//
//endfor
//
//setdatafolder root:
//
//string trtempgraphstr = "trtempgraph"
//display/n=$trtempgraphstr/hide=1
//
//setdatafolder $trtemppath
//variable npnts = itemsinlist(tl)
//
//for(k=0;k<freqs;k+=1)
//	 string fntempstring2 = sam+"_trtemp_"+num2str((leftx(nwave0)+k*freqstep)*1000)
//	wave fntempwave2=$fntempstring2
//	 appendtograph/w=trtempgraph fntempwave2 vs temperature
//endfor

for(i=0;i<freqs;i+=1)

	setdatafolder $trtemppath
	string trvstempstr = sam+"_trtemp_"+num2str((leftx(nwave0)+i*freqstep)*1000)
	make/c/o/n=(itemsinlist(tl)) $trvstempstr/wave=trvstempwave
	setdatafolder $trmagpath 
		    	
		for(j=0;j<itemsinlist(tl);j+=1)

   		    string trmagstr = sam+"_trmag_"+stringfromlist(j,tl)
		    wave trmagwave = $trmagstr
	
		    trvstempwave[j]=trmagwave(freqstep * i + leftx(nwave0))
	
             endfor

endfor

setdatafolder root:

string trtempgraphstr = "trtempgraph"
display/n=$trtempgraphstr/hide=1

setdatafolder $trtemppath
variable npnts = itemsinlist(tl)

for(k=0;k<freqs;k+=1)
	 string fntempstring2 = sam+"_trtemp_"+num2str((leftx(nwave0)+k*freqstep)*1000)
	wave fntempwave2=$fntempstring2
	 appendtograph/w=trtempgraph fntempwave2 vs temperature
endfor





setdatafolder $logypath
nvar trtemplogy
setdatafolder $logxpath
nvar trtemplogx
setdatafolder root:

graphervstemp("trtempgraph","Transmission Magnitude","\\Z30 T","Rainbow",3,trtemplogy,trtemplogx)


end


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function sigmavstempNEW()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


setdatafolder root:
svar sigma1path, sigma2path, tl, sam,logypath,logxpath
nvar alfreq, aufreq, freqs, i, j, k
wave temperature

windowkill("sigma1tempgraph")
windowkill("sigma2tempgraph")
windowkill("sigmamagtempgraph")

string/g sigma1temppath = sigma1path + "vstemp:"
overwritefolder(sigma1temppath)

string/g sigma2temppath = sigma2path+ "vstemp:"
overwritefolder(sigma2temppath)

sigmamagNEW()
svar sigmamagpath
string/g sigmamagtemppath = sigmamagpath + "vstemp:"
overwritefolder(sigmamagtemppath)
        
setdatafolder $sigma1path

// these definitions have zeros because we're choosing a basis wave
// to extract other numbers from 
string sigma1str0 = sam+"_sigma1_"+stringfromlist(0,tl)
wave sigma1wave0 = $sigma1str0


// POTENTIAL DEBUGGING
// i'm trying something because there is an error where there is a frequency cut generated under particular conditions
// at a frequency greater than upper freq.
variable freqstep = (aufreq-alfreq) / freqs

/////// sigma1 loop

string sigma1tempgraphstr = "sigma1tempgraph"
display/n=$sigma1tempgraphstr/hide=1

for(i=0;i<freqs;i+=1)

setdatafolder $sigma1temppath

string sigma1tempstr = sam+"_sigma1temp_"+num2str((leftx(sigma1wave0)+i*freqstep)*1000)
make/o/n=(itemsinlist(tl)) $sigma1tempstr/wave=sigma1tempwave

setdatafolder $sigma1path	

	for(j=0;j<itemsinlist(tl);j+=1)

   		    string sigma1str = sam+"_sigma1_"+stringfromlist(j,tl)
		    wave sigma1wave = $sigma1str
		    sigma1tempwave[j]=sigma1wave(freqstep * i + leftx(sigma1wave0))
       endfor

endfor

setdatafolder $sigma1temppath

for(k=0;k<freqs;k+=1)
	string sigma1tempstr1 =sam+"_sigma1temp_"+num2str((leftx(sigma1wave0)+k*freqstep)*1000)
	make/o/n=(itemsinlist(tl)) $sigma1tempstr1/wave=sigma1tempwave1
	 appendtograph/w=sigma1tempgraph sigma1tempwave1 vs temperature
endfor

setdatafolder $logypath
nvar sigma1templogy
setdatafolder $logxpath
nvar sigma1templogx
setdatafolder root:

graphervstemp("sigma1tempgraph","\Z18\F'symbol's\F'geneva'\B1\M [(\F'symbol'W\F'geneva' cm)\S-1\M]","\Z30 \F'Symbol's\Z30\B1","Rainbow",3,sigma1templogy,sigma1templogx)

setdatafolder $sigma2path


// these definitions have zeros because we're choosing a basis wave
// to extract other numbers from 

/////// sigma2 loop

string sigma2tempgraphstr = "sigma2tempgraph"
display/n=$sigma2tempgraphstr/hide=1


for(i=0;i<freqs;i+=1)

setdatafolder $sigma2temppath

string sigma2tempstr = sam+"_sigma2temp_"+num2str((leftx(sigma1wave0)+i*freqstep)*1000)
make/o/n=(itemsinlist(tl)) $sigma2tempstr/wave=sigma2tempwave

setdatafolder $sigma2path	

	for(j=0;j<itemsinlist(tl);j+=1)

   		    string sigma2str = sam+"_sigma2_"+stringfromlist(j,tl)
		    wave sigma2wave = $sigma2str
		    sigma2tempwave[j]=sigma2wave(freqstep * i + leftx(sigma1wave0))
       endfor

endfor

setdatafolder $sigma2temppath

for(k=0;k<freqs;k+=1)
	string sigma2tempstr1 =sam+"_sigma2temp_"+num2str((leftx(sigma1wave0)+k*freqstep)*1000)
	make/o/n=(itemsinlist(tl)) $sigma2tempstr1/wave=sigma2tempwave1
	 appendtograph/w=sigma2tempgraph sigma2tempwave1 vs temperature
endfor

setdatafolder $logypath
nvar sigma2templogy
setdatafolder $logxpath
nvar sigma2templogx
setdatafolder root:

graphervstemp("sigma2tempgraph","\Z18\F'symbol's\F'geneva'\B2\M [(\F'symbol'W\F'geneva' cm)\S-1\M]","\Z30 \F'Symbol's\Z30\B2","Rainbow",3,sigma2templogy,sigma2templogx)

end



/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function epsvstempNEW()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


setdatafolder root:
svar eps1path, eps2path, tl, sam,logypath,logxpath
nvar alfreq, aufreq, freqs, i, j, k
wave temperature

windowkill("eps1tempgraph")
windowkill("eps2tempgraph")
windowkill("epsmagtempgraph")

string/g eps1temppath = eps1path + "vstemp:"
overwritefolder(eps1temppath)

string/g eps2temppath = eps2path+ "vstemp:"
overwritefolder(eps2temppath)

        
setdatafolder $eps1path

// these definitions have zeros because we're choosing a basis wave
// to extract other numbers from 
string eps1str0 = sam+"_eps1_"+stringfromlist(0,tl)
wave eps1wave0 = $eps1str0
variable freqstep = (aufreq-alfreq) / freqs

/////// eps1 loop

string eps1tempgraphstr = "eps1tempgraph"
display/n=$eps1tempgraphstr/hide=1

for(i=0;i<freqs;i+=1)

setdatafolder $eps1temppath

string eps1tempstr = sam+"_eps1temp_"+num2str((leftx(eps1wave0)+i*freqstep)*1000)
make/o/n=(itemsinlist(tl)) $eps1tempstr/wave=eps1tempwave

setdatafolder $eps1path	

	for(j=0;j<itemsinlist(tl);j+=1)

   		    string eps1str = sam+"_eps1_"+stringfromlist(j,tl)
		    wave eps1wave = $eps1str
		    eps1tempwave[j]=eps1wave(freqstep * i + leftx(eps1wave0))
       endfor

endfor

setdatafolder $eps1temppath

for(k=0;k<freqs;k+=1)
	string eps1tempstr1 =sam+"_eps1temp_"+num2str((leftx(eps1wave0)+k*freqstep)*1000)
	make/o/n=(itemsinlist(tl)) $eps1tempstr1/wave=eps1tempwave1
	 appendtograph/w=eps1tempgraph eps1tempwave1 vs temperature
endfor

setdatafolder $logypath
nvar eps1templogy
setdatafolder $logxpath
nvar eps1templogx
setdatafolder root:

graphervstemp("eps1tempgraph","\Z18\F'symbol'e\F'geneva'\B1\M","\Z30 \F'Symbol'e\Z30\B1","Rainbow",3,eps1templogy,eps1templogx)

setdatafolder $eps2path


// these definitions have zeros because we're choosing a basis wave
// to extract other numbers from 

/////// eps2 loop

string eps2tempgraphstr = "eps2tempgraph"
display/n=$eps2tempgraphstr/hide=1


for(i=0;i<freqs;i+=1)

setdatafolder $eps2temppath

string eps2tempstr = sam+"_eps2temp_"+num2str((leftx(eps1wave0)+i*freqstep)*1000)
make/o/n=(itemsinlist(tl)) $eps2tempstr/wave=eps2tempwave

setdatafolder $eps2path	

	for(j=0;j<itemsinlist(tl);j+=1)

   		    string eps2str = sam+"_eps2_"+stringfromlist(j,tl)
		    wave eps2wave = $eps2str
		    eps2tempwave[j]=eps2wave(freqstep * i + leftx(eps1wave0))
       endfor

endfor

setdatafolder $eps2temppath

for(k=0;k<freqs;k+=1)
	string eps2tempstr1 =sam+"_eps2temp_"+num2str((leftx(eps1wave0)+k*freqstep)*1000)
	make/o/n=(itemsinlist(tl)) $eps2tempstr1/wave=eps2tempwave1
	 appendtograph/w=eps2tempgraph eps2tempwave1 vs temperature
endfor

setdatafolder $logypath
nvar eps2templogy
setdatafolder $logxpath
nvar eps2templogx
setdatafolder root:

graphervstemp("eps2tempgraph","\Z18\F'symbol'e\F'geneva'\B2\M","\Z30 \F'Symbol'e\Z30\B2","Rainbow",3,eps2templogy,eps2templogx)

end


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function indexvstempNEW()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


setdatafolder root:
svar npath, kpath, tl, sam,logypath, logxpath
nvar alfreq, aufreq, freqs, i, j, k
wave temperature

windowkill("ntempgraph")
windowkill("ktempgraph")
windowkill("epsmagtempgraph")

string/g ntemppath = npath + "vstemp:"
overwritefolder(ntemppath)

string/g ktemppath = kpath+ "vstemp:"
overwritefolder(ktemppath)
        
setdatafolder $npath

// these definitions have zeros because we're choosing a basis wave
// to extract other numbers from 
string nstr0 = sam+"_n_"+stringfromlist(0,tl)
wave nwave0 = $nstr0
variable freqstep = (aufreq-alfreq) / freqs

/////// n loop

string ntempgraphstr = "ntempgraph"
display/n=$ntempgraphstr/hide=1

for(i=0;i<freqs;i+=1)

setdatafolder $ntemppath

string ntempstr = sam+"_ntemp_"+num2str((leftx(nwave0)+i*freqstep)*1000)
make/o/n=(itemsinlist(tl)) $ntempstr/wave=ntempwave

setdatafolder $npath	

	for(j=0;j<itemsinlist(tl);j+=1)

   		    string nstr = sam+"_n_"+stringfromlist(j,tl)
		    wave nwave = $nstr
		    ntempwave[j]=nwave(freqstep * i + leftx(nwave0))
       endfor

endfor

setdatafolder $ntemppath

for(k=0;k<freqs;k+=1)
	string ntempstr1 =sam+"_ntemp_"+num2str((leftx(nwave0)+k*freqstep)*1000)
	make/o/n=(itemsinlist(tl)) $ntempstr1/wave=ntempwave1
	 appendtograph/w=ntempgraph ntempwave1 vs temperature
endfor

setdatafolder $logypath
nvar ntemplogy
setdatafolder $logxpath
nvar ntemplogx
setdatafolder root:

graphervstemp("ntempgraph","\Z18n","\Z30n ","Rainbow",3,ntemplogy,ntemplogx)

setdatafolder $kpath


// these definitions have zeros because we're choosing a basis wave
// to extract other numbers from 

/////// k loop

string ktempgraphstr = "ktempgraph"
display/n=$ktempgraphstr/hide=1


for(i=0;i<freqs;i+=1)

setdatafolder $ktemppath

string ktempstr = sam+"_ktemp_"+num2str((leftx(nwave0)+i*freqstep)*1000)
make/o/n=(itemsinlist(tl)) $ktempstr/wave=ktempwave

setdatafolder $kpath	

	for(j=0;j<itemsinlist(tl);j+=1)

   		    string kstr = sam+"_k_"+stringfromlist(j,tl)
		    wave kwave = $kstr
		    ktempwave[j]=kwave(freqstep * i + leftx(nwave0))
       endfor

endfor

setdatafolder $ktemppath

for(k=0;k<freqs;k+=1)
	string ktempstr1 =sam+"_ktemp_"+num2str((leftx(nwave0)+k*freqstep)*1000)
	make/o/n=(itemsinlist(tl)) $ktempstr1/wave=ktempwave1
	 appendtograph/w=ktempgraph ktempwave1 vs temperature
endfor

setdatafolder $logypath
nvar ktemplogy
setdatafolder $logxpath
nvar ktemplogx
setdatafolder root:

graphervstemp("ktempgraph","\Z18k","\Z30k","Rainbow",3,ktemplogy,ktemplogx)

end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  

function imagegraphmaker()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


// this function handles all three optical constants (eps, index, sigma) in one function,

setdatafolder root:
nvar i, j
svar eps1path, eps2path, npath, kpath, sigma1path, sigma2path,tl,sam, logypath, logxpath, imagecolor, contourcolor
wave freq, temperature

reverse temperature/d=revtemp

imagefreqtempwaves()
wave imagefreq, imagetemp

string ocpath = "eps:one:;eps:two:;index:n:;index:k:;sigma:one:;sigma:two:;alpha:;eps:one:ratio:;sigma:one:ratio:;tr:mag:ratio:;"
string ocname ="eps1;eps2;n;k;sigma1;sigma2;alpha;eps1tempratio;sigma1tempratio;trmagtempratio;"
string graphtitle = "\Z16\F'Symbol'e\Z16\B1;\Z16\F'Symbol'e\Z16\B2;\Z16n;\Z16k;\Z16\F'Symbol's\Z16\B1;\Z16\F'Symbol's\Z16\B2;\Z16\F'Symbol'a\Z16\B;Ratio: \Z16\F'Symbol'e\Z16\B1;Ratio: \Z16\F'Symbol's\Z16\B1;Ratio: Tr;"

for(i=0;i<itemsinlist(ocpath);i+=1)

setdatafolder $logypath
string logystr = stringfromlist(i,ocname)+"imagelogy"
nvar logyvar = $logystr

setdatafolder $logxpath
string logxstr = stringfromlist(i,ocname)+"imagelogx"
nvar logxvar = $logxstr

string pathstr = "root:"+stringfromlist(i,ocpath)
setdatafolder $pathstr
string imagestr = ""
string imagegraph = stringfromlist(i,ocname)+"imagegraph"
windowkill(imagegraph)


for(j=0;j<itemsinlist(tl);j+=1)
string appendstr = sam +"_"+stringfromlist(i,ocname)+"_"+stringfromlist(j,tl)
imagestr = addlistitem(appendstr,imagestr)
endfor
concatenate/o imagestr, imagematrix
display/hide=1/n=$imagegraph;appendimage imagematrix vs {imagefreq,imagetemp}
appendmatrixcontour/f="%.17s=%g"/w=$imagegraph imagematrix vs {freq,revtemp}
ModifyContour imagematrix ctabLines={*,*,$contourcolor,0}
ModifyImage imagematrix ctab= {*,*,$imagecolor,0}
imagegrapher(stringfromlist(i,graphtitle), logyvar, logxvar)
endfor
end


end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function graphervsfreq(name,wavepath,startstring,leftlabel,txtbox,colortable, cmode,logy,logx)

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

// name is input as a string (in quotes)
// it is used in the display and appendtograph commands with $
// wavepath is the path to where the waves live
// startstring is the name of the string preceding the temperature (should include final underscore)
// leftlabel labels the vertical axis
// txtbox is a symbol that is placed within the graph for ease of identification
// colortable defines which color scheme will color the waves (see graph > packages > make traces different
// cmode is which complex mode (real, imag, mag, etc.) that is plotted

string name, wavepath, startstring, leftlabel,txtbox, colortable
variable cmode, logy,logx
setdatafolder root:
svar tl
nvar lfreq, ufreq, npnts,k
wave/t temptickslabel
wave temptickspos

display/n=$name/hide=1
setdatafolder $wavepath

for(k=0;k<npnts;k+=1)

string graphstr = startstring + stringfromlist(k,tl) 
wave graphwave = $graphstr
setdatafolder $wavepath
appendtograph/w=$name graphwave
endfor

label left leftlabel
Label bottom "\Z14Frequency [THz]"
ModifyGraph mirror=2
TextBox/C/N=text0/F=0/A=LT txtbox
ModifyGraph fSize(left)=14
ModifyGraph fSize(bottom)=14
ModifyGraph lsize=2
ModifyGraph cmplxmode= cmode
setaxis bottom lfreq,ufreq
//modifygraph log(left)=1
SetAxis/A=2 left
ApplyColorTableToTopGraphNEW(colortable)
modifygraph grid(left)=1,fstyle=1,axthick=1.5
ColorScale/C/N=text1/F=0/B=1 vert=0, ctab={0,1,ColdWarm,0};DelayUpdate
ColorScale/C/N=text1 userTicks={temptickspos,temptickslabel}
ColorScale/C/N=text1 "Temperature [K]"
ColorScale/C/N=text1 heightPct=3
ColorScale/C/N=text1/A=RT/X=1.31/Y=3.47
ModifyGraph log(left)=logy
modifygraph log(bottom)=logx
end



/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function graphervstemp(name,leftlabel,txtbox,colortable, cmode,logy,logx)

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



// name is input as a string (in quotes)
// it is used in the display and appendtograph commands with $
// wavepath is the path to where the waves live
// startstring is the name of the string preceding the temperature (should include final underscore)
// leftlabel labels the vertical axis
// txtbox is a symbol that is placed within the graph for ease of identification
// colortable defines which color scheme will color the waves (see graph > packages > make traces different
// cmode is which complex mode (real, imag, mag, etc.) that is plotted

string name, leftlabel,txtbox, colortable
variable cmode,logy,logx
setdatafolder root:
svar tl,notabletemps, notablecolor
nvar lfreq, ufreq, i, npnts, notabledash, notablethickness, k
wave temperature


dowindow/hide=1 $name


label left leftlabel
Label bottom "\Z14Temperature [K]"
ModifyGraph mirror=2
TextBox/C/N=text0/F=0/A=LT txtbox
ModifyGraph fSize(left)=14
ModifyGraph fSize(bottom)=14
ModifyGraph lsize=2
ModifyGraph cmplxmode= cmode
//modifygraph log(left)=1
SetAxis/A=2 left
ApplyColorTableToTopGraphNEW(colortable)
modifygraph grid(left)=1,fstyle=1,axthick=1.5
ColorScale/C/N=text1/F=0/A=LT vert=0, ctab={0,1,Rainbow,0};DelayUpdate
ColorScale/C/N=text1 "Frequency [THz]"
ColorScale/C/N=text1 widthPct=35,heightPct=3
ColorScale/C/N=text1/B=1
TextBox/C/N=text0/B=1
TextBox/C/N=text0/A=LC/X=3.32/Y=1.39
ColorScale/C/N=text1  ctab={lfreq,ufreq,Rainbow,0}
ModifyGraph log(left)=logy
modifygraph log(bottom)=logx
doupdate
getaxis/q left

variable rvar, gvar, bvar
rvar = str2num(stringfromlist(0, notablecolor,","))
gvar = str2num(stringfromlist(1, notablecolor,","))
bvar = str2num(stringfromlist(2, notablecolor,","))

if(strlen(notabletemps)>0)
for(k=0;k<itemsinlist(notabletemps);k+=1)
ShowTools/A arrow
SetDrawEnv ycoord= left;DelayUpdate
SetDrawEnv xcoord= bottom
SetDrawEnv linefgc=(rvar,gvar,bvar),dash= notabledash-1,linethick= notablethickness
DrawLine str2num(stringfromlist(k,notabletemps)), v_min, str2num(stringfromlist(k,notabletemps)), v_max
hidetools/a
endfor
endif

end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function imagegrapher(graphtitle,logyvar,logxvar)

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


string graphtitle
variable logyvar, logxvar
setdatafolder root:
svar notabletemps, notablecolor
nvar k, lfreq,ufreq, notabledash, notablethickness
Label left "Temperature [K]"
Label bottom "Frequency [THz]"
ModifyGraph fStyle=1,fSize=14,axThick=2
TextBox/C/N=text0/S=1/D={1,3}/A=MT graphtitle
modifygraph log(left) = logyvar
modifygraph log(bottom) = logxvar
doupdate

variable rvar, gvar, bvar
rvar = str2num(stringfromlist(0, notablecolor,","))
gvar = str2num(stringfromlist(1, notablecolor,","))
bvar = str2num(stringfromlist(2, notablecolor,","))

if(strlen(notabletemps)>0)
for(k=0;k<itemsinlist(notabletemps);k+=1)
ShowTools/A arrow
SetDrawEnv ycoord= left;DelayUpdate
SetDrawEnv xcoord= bottom
SetDrawEnv linefgc=(rvar,gvar,bvar),dash= notabledash-1,linethick= notablethickness
DrawLine lfreq,str2num(stringfromlist(k,notabletemps)),ufreq,str2num(stringfromlist(k,notabletemps))
hidetools/a
endfor
endif


end


  
//from igor exchange: http://www.igorexchange.com/node/2555

// the WIN:3 command kills graphs and tables (see winlist documentation)

Function killallgraphsandtablesNEW()
	string fulllist = winlist("*", ";","win:3"), name
	variable i
 
for(i=0; i<itemsinlist(fulllist); i +=1)
name= stringfromlist(i, fulllist)
// dowindow/k kills the open window
dowindow/k $name	
endfor

end


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function pulsegraphNEW()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\




setdatafolder root: 
svar sam, ref, tl, samref
nvar step, npnts, i, timestep, j
	            
windowkill("pulsegraph")
        
 string/g sampath = "root:pulses:"+sam+":"
 string/g refpath = "root:pulses:"+ref+":"
                                                                                                                                        
overwritefolder(sampath)
overwritefolder(refpath)
                                                                    
display/n=pulses/hide=1

for(i=0;i<itemsinlist(samref);i+=1)                     
for(j=0;j<npnts;j+=1)

setdatafolder $"root:raw_pulses:"+stringfromlist(i,samref)+":"
string pulsestr = stringfromlist(i,samref)+"_"+stringfromlist(j,tl)
wave pulsewave =$pulsestr
setscale/P x 0,timestep,"", pulsewave
duplicate/o pulsewave $stringfromlist(i,samref)+"_"+stringfromlist(j,tl)+"_copy"
movewave $stringfromlist(i,samref)+"_"+stringfromlist(j,tl)+"_copy", $"root:pulses:"+stringfromlist(i,samref)+":"


endfor
endfor

 end






 
/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function tempticksNEW(tempticks0)

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\




string tempticks0

setdatafolder root:
string/g tempticks=tempticks0

svar tl
nvar i,npnts
wave temperature

// the colors of the waves are assigned with a linear spread 
// i.e. 10 K and 15 K are spaced as far apart in "color space"
// as 15 K and 50 K if only 10, 15 and 50 K 

// the trick to get a good temperature bar that reflects the true spacing
// is to use the position of a temperature in the list of temperatures

// but the waves must have the true number of ticks

variable numtempticks = itemsinlist(tempticks)
make/t/o/n=(numtempticks) temptickslabel
make/o/n=(numtempticks) temptickspos

for(i=0;i<numtempticks;i+=1)
temptickslabel[i]=stringfromlist(i,tempticks)
temptickspos[i] = whichlistitem(stringfromlist(i,tempticks),tl) / (npnts-1)
endfor




end


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function overwritefolder(foldername)

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


string foldername

	     if(datafolderexists(foldername)==1)
            killdatafolder $foldername
	     newdatafolderpathNEW(foldername)
	     else
	     newdatafolderpathNEW(foldername)
	     endif
end


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function windowkill(graphname)

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

string graphname

	  dowindow $graphname
        if  (V_Flag == 1 || v_flag==2)
        killwindow $graphname
        endif

end
        


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  
 
function recalculateandregraphNEW()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


setdatafolder root: 
nvar i

// this line captures all of the graphs that were open
// they will subsequently be killed and then recreated at the end of the function
dowindow/hide=1 sigma1vrhgraph
dowindow/hide=1 eps1tempratiograph
dowindow/hide=1 eps2tempratiograph
dowindow/hide=1 ntempratiograph 
dowindow/hide=1 ktempratiograph 
dowindow/hide=1 sigma1tempratiograph
dowindow/hide=1 sigma2tempratiograph
dowindow/hide=1 trtempratiograph
string visiblegraphs = winlist("*",",","win:1,visible:1")

// alphabetize graph list so that graphs are displayed in the same order

visiblegraphs = sortlist(visiblegraphs,",",16)

             
killallgraphsandtablesNEW()             

// start from scratch; i'm sure there is a better way to do this, but for now this is sufficient

// trtemperature is run after slabindex since the lower frequency is determined in slabindex (it is 
// the first freq. data point above the user input lower freq; see the bounds set in slabindex)
pulserefreshNEW()
trNEW()
slabepsNEW()
slabindexNEW()
slabsigmaNEW()
powerabsNEW()
epsvstempNEW()
trtemperatureNEW()
indexvstempNEW()
sigmavstempNEW()
alphavstemp()
ratiotolowtempNEW()
ratiotolowtemptrNEW()
imagegraphmaker()


// the graphs that were created were created with the visible setting
// the next function hides them all

hideallgraphsandtablesNEW()

// the loop and dowindow/hide=0 command reopens the windows that were open at the beginning of this function

		 for(i=0;i<itemsinlist(visiblegraphs,",");i+=1)
		 dowindow/hide=0 $stringfromlist(i,visiblegraphs,",")
		 endfor
		 
		 if(strlen(visiblegraphs)>0)		 
		 string tilecmd = "tilewindows/w=(0,0,1063,800) "+visiblegraphs
		 execute tilecmd
		 endif

end



/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

Function ApplyColorTableToTopGraphNEW(ctabname)

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


	String ctabname
 
	String graphName = WinName(0, 1)
	if (strlen(graphName) == 0)
		return -1
	endif
 
	Variable numTraces = ItemsInList(TraceNameList(graphName,";",3))
 
	if (numTraces <= 0)
		return -1
	endif
 
	Variable denominator= numTraces-1
	if( denominator < 1 )
		denominator= 1    // avoid divide by zero, use just the first color for 1 trace
	endif
 
	ColorTab2Wave $ctabname	// creates M_colors
	Wave M_colors
	Variable numRows= DimSize(M_colors,0)
	Variable red, green, blue
	Variable i, index
	for(i=0; i<numTraces; i+=1)
		index = round(i/denominator * (numRows-1))	// spread entire color range over all traces.
		ModifyGraph/W=$graphName rgb[i]=(M_colors[index][0], M_colors[index][1], M_colors[index][2])
	endfor
	return 0
End



/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

Function/S FilePathToWaveNameNEW(path)

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



	String path		// Assumed to use colon separators
 
	// Get last element of path without extension
	String name = ParseFilePath(3, path, ":", 0, 0)	// e.g., MyFile.txt
 
	name = CleanupName(name, 0)
 
	return name
End


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function hideallgraphsandtablesNEW()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


	string fulllist = winlist("*", ";","win:3"), name
	variable i
 
for(i=0; i<itemsinlist(fulllist); i +=1)
name= stringfromlist(i, fulllist)
// dowindow/k kills the open window


dowindow/hide=1 $name
//if(v_flag == 2)	
//dowindow/hide=1 $name
//endif
endfor

end



/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

Function/S loadfilesNEW(samref)

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



// identify samref as a string and import global variables

string samref

setdatafolder root:
svar sam, ref
variable npnts, i, tl
	
// this next function will kill a folder root:sam (for example) if it exists
// this allows quick change of parameters and reloading	

if(datafolderexists("root:raw_pulses:"+samref+":")==1)
	killdatafolder $"root:raw_pulses:"+samref+":"
endif 	

// cmpstr is the compare string function
// 0 is defined as equal
// the 1 here is necessary to set the folder being loaded into
	
if(cmpstr(sam,samref)==0)
	string/g rawsampath = "root:raw_pulses:"+samref+":"
	newdatafolderpathNEW(rawsampath,set=1)	
else
	string/g rawrefpath = "root:raw_pulses:"+samref+":"
	newdatafolderpathNEW(rawrefpath,set=1)
endif	

	
	Variable refNum
	String message = "Select one or more files"
	String outputPaths,tempname
	String fileFilters = "Data Files (*.txt,*.dat,*.csv,*.tsv):.txt,.dat,.csv,.tsv;"
	fileFilters += "All Files:.*;"
 
	Open /D /R /MULT=1 /F=fileFilters /M=message refNum
	outputPaths = S_fileName
 
	if (strlen(outputPaths) == 0)
		Print "Cancelled"
	else
		Variable numFilesSelected = ItemsInList(outputPaths, "\r")
		for(i=0; i<numFilesSelected; i+=1)
			String path = StringFromList(i, outputPaths, "\r")
			// Printf "%d: %s\r", i, path
			// Add commands here to load the actual waves.  An example command
			// is included below but you will need to modify it depending on how
			// the data you are loading is organized.
			tempname=FilePathToWaveNameNEW(path)
			LoadWave/q/G/K=0/L={0,0,0,1,1}/A=$tempname path
			Duplicate/o $tempname+"0",$tempname
			killwaves $tempname+"0"
		endfor
	endif
 
setdatafolder("root:")
	
 
	return outputPaths		// Will be empty if user canceled
	
	
	end
	

	
/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\	
	
Function newdatafolderpathNEW(path[,set])	

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\	

//from IgorExchange
// http://www.igorexchange.com/node/402
// use optional set by inputting set=1

	string path
	variable set
 
	variable depth=itemsinlist(path,":"),i
	string partial=stringfromlist(0,path,":")
	if(strlen(partial)==0)	//path is relative, beginning with a :
		partial="";i=1
	elseif(cmpstr("root",partial)==0) //path is full from root
		partial="root";i=1
	else						//path is relative, with no initial :
		partial="";i=0
	endif
	for(i=i;i<depth;i+=1)
		partial+=":"+possiblyquotename(cleanupname(StringFromList(i,path,":"),1))
		newdatafolder/o $partial
	endfor
	if(set)
		setdatafolder $partial
	endif
end




 
 /////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
 
 function vrhNEW(exponent)
 
 /////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


variable exponent
setdatafolder root:
wave temperature
nvar npnts

duplicate/o temperature, tempvrh
tempvrh = temperature^(-exponent)

end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function sigmamagNEW()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 



setdatafolder root:
svar sigma1path, sigma2path, tl, sam

string/g sigmamagpath = "root:sigma:mag:"
overwritefolder(sigmamagpath)

variable i
variable npnts = itemsinlist(tl)

for(i=0;i<npnts;i+=1)

	setdatafolder $sigma1path
	string sigma1str = sam+"_sigma1_"+stringfromlist(i,tl)
	wave sigma1wave = $sigma1str

	setdatafolder $sigma2path
	string sigma2str = sam+"_sigma2_"+stringfromlist(i,tl)
	wave sigma2wave = $sigma2str
	
	setdatafolder $sigmamagpath

	string sigmamagstr= sam+"_sigmamag_"+stringfromlist(i,tl)
	make/o/n=(numpnts(sigma1wave)) $sigmamagstr/wave=sigmamagwave
	
	sigmamagwave = (sigma1wave^2+sigma2wave^2)^.5
	setscale/P x leftx(sigma1wave),dimdelta(sigma1wave,0),"", sigmamagwave
endfor

setdatafolder root:
end


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  

function ratiotolowtempNEW()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 



// this function handles all three optical constants (eps, index, sigma) in one function,
// utilizing 

setdatafolder root:
nvar i, j
svar eps1path, eps2path, npath, kpath, sigma1path, sigma2path,tl,sam

windowkill("eps1tempratiograph")
windowkill("eps2tempratiograph")
windowkill("ntempratiograph")
windowkill("ktempratiograph")
windowkill("sigma1tempratiograph")
windowkill("sigma2tempratiograph")
windowkill("eps1tempratioimagegraph")
windowkill("sigma1tempratioimagegraph")
windowkill("trmagtempratioimagegraph")


string oc = "eps;index;sigma"
for(i=0;i<itemsinlist(oc);i+=1)

string pathstr = "root:"+stringfromlist(i,oc)+":"

// if statement since "one" and "two" are used for eps and sigma but n and k for index
// make ratio folders for real and imaginary parts

if(cmpstr(stringfromlist(i,oc),"eps")==0 || cmpstr(stringfromlist(i,oc),"sigma")==0)
string realratiopath = pathstr+"one:ratio:"
string imagratiopath = pathstr+"two:ratio:"
else
realratiopath = pathstr+"n:ratio:"
imagratiopath = pathstr+"k:ratio:"
endif
overwritefolder(realratiopath)
overwritefolder(imagratiopath)


// get low temperature wave
// only do this once, so it doesn't need to be in the temperature loop


if(cmpstr(stringfromlist(i,oc),"eps")==0 || cmpstr(stringfromlist(i,oc),"sigma")==0)
string lowtemprealstr = sam+"_"+stringfromlist(i,oc)+"1_"+stringfromlist(0,tl)
string lowtempimagstr = sam+"_"+stringfromlist(i,oc)+"2_"+stringfromlist(0,tl)
else
lowtemprealstr = sam+"_n_"+stringfromlist(0,tl)
lowtempimagstr = sam+"_k_"+stringfromlist(0,tl)
endif



// temperature loop for pulling high temp waves and ratioing them to the 
// lowest; loop starts at 0 for integration with graphervsfreq, but really
//  loop should start at 1 since we don't need to ratio the lowest to itself

for(j=0;j<itemsinlist(tl);j+=1)
if(cmpstr(stringfromlist(i,oc),"eps")==0 || cmpstr(stringfromlist(i,oc),"sigma")==0)

string hightemprealstr = sam+"_"+stringfromlist(i,oc)+"1_"+stringfromlist(j,tl)
string hightempimagstr = sam+"_"+stringfromlist(i,oc)+"2_"+stringfromlist(j,tl)
string hightemprealratiostr =  sam+"_"+stringfromlist(i,oc)+"1tempratio_"+stringfromlist(j,tl)
string hightempimagratiostr =  sam+"_"+stringfromlist(i,oc)+"2tempratio_"+stringfromlist(j,tl)



setdatafolder $pathstr+"one:"
wave lowtemprealwave = $lowtemprealstr
wave hightemprealwave = $hightemprealstr

setdatafolder $pathstr+"two:"
wave lowtempimagwave = $lowtempimagstr
wave hightempimagwave = $hightempimagstr


else

hightemprealstr = sam+"_n_"+stringfromlist(j,tl)
hightempimagstr = sam+"_k_"+stringfromlist(j,tl)

setdatafolder $pathstr+"n:"
wave lowtemprealwave = $lowtemprealstr
wave hightemprealwave = $hightemprealstr
setdatafolder $pathstr+"k:"
wave lowtempimagwave = $lowtempimagstr
wave hightempimagwave = $hightempimagstr
hightemprealratiostr =  sam+"_n_"+"tempratio_"+stringfromlist(j,tl)
hightempimagratiostr =  sam+"_k_"+"tempratio_"+stringfromlist(j,tl)


endif



setdatafolder $realratiopath
make/n=(numpnts(lowtemprealwave))/o $hightemprealratiostr/wave=lowtemprealratiowave
lowtemprealratiowave = hightemprealwave / lowtemprealwave
setdatafolder $imagratiopath
make/n=(numpnts(lowtemprealwave))/o $hightempimagratiostr/wave=lowtempimagratiowave
lowtempimagratiowave = hightempimagwave / lowtempimagwave
setscale/P x leftx(lowtemprealwave),dimdelta(lowtemprealwave,0),"", lowtemprealratiowave, lowtempimagratiowave


endfor
endfor

// the graphing code will be standalone since it is optical constant specific

// make global path strings
setdatafolder root:
string/g eps1tempratiopath = "root:eps:one:ratio:"
string/g eps2tempratiopath = "root:eps:two:ratio:"
string/g ntempratiopath = "root:index:n:ratio:"
string/g ktempratiopath = "root:index:k:ratio:"
string/g sigma1tempratiopath = "root:sigma:one:ratio:"
string/g sigma2tempratiopath = "root:sigma:two:ratio:"

// default graphs will just have log settings set to zero

// graph only eps1 and sigma 1; the rest are similar to these two based on ratioing

graphervsfreq("eps1tempratiograph",eps1tempratiopath,sam+"_eps1tempratio_","\Z18Ratio: \F'symbol'e\F'geneva'\B1\M ","\Z30 \F'Symbol'e\Z30\B1","ColdWarm",3,0,0)
//graphervsfreq("eps2tempratiograph",eps2tempratiopath,sam+"_eps2tempratio_","\Z18Ratio:\F'symbol'e\F'geneva'\B2\M ","\Z30 \F'Symbol'e\Z30\B2","ColdWarm",3,0,0)
//graphervsfreq("ntempratiograph",ntempratiopath,sam+"_ntempratio_","\Z18Ratio: n","\Z30 n","ColdWarm",3,0,0)
//graphervsfreq("ktempratiograph",ktempratiopath,sam+"_ktempratio_","\Z18Ratio: k","\Z30 k","ColdWarm",3,0,0)
graphervsfreq("sigma1tempratiograph",sigma1tempratiopath,sam+"_sigma1tempratio_","\Z18Ratio:\F'symbol's\F'geneva'\B1\M [(\F'symbol'W\F'geneva' cm)\S-1\M]","\Z30 \F'Symbol's\Z30\B1","ColdWarm",3,0,0)
//graphervsfreq("sigma2tempratiograph",sigma2tempratiopath,sam+"_sigma2tempratio_","\Z18Ratio:\F'symbol's\F'geneva'\B2\M [(\F'symbol'W\F'geneva' cm)\S-1\M]","\Z30 \F'Symbol's\Z30\B2","ColdWarm",3,0,0)



dowindow/hide=0 eps1tempratiograph
dowindow/hide=0 eps2tempratiograph
dowindow/hide=0 ntempratiograph
dowindow/hide=0 ktempratiograph
dowindow/hide=0 sigma1tempratiograph
dowindow/hide=0 sigma2tempratiograph


setdatafolder root:

end


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  

function ratiotolowtemptrNEW()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


// this function handles all three optical constants (eps, index, sigma) in one function,
// utiliziing 

setdatafolder root:
nvar i, j, lfreq, ufreq
svar trpath,tl,sam

windowkill("trmagtempratiograph")
windowkill("trmagtempratioimagegraph")

string oc = "trmag;"
string pathstr = "root:tr:mag:"
for(i=0;i<itemsinlist(oc);i+=1)



// if statement since "one" and "two" are used for eps and sigma but n and k for index
// make ratio folders for real and imaginary parts

string ratiopath = stringfromlist(i,pathstr)+"ratio:"
overwritefolder(ratiopath)

// get low temperature wave
// only do this once, so it doesn't need to be in the temperature loop


string lowtempstr = sam+"_"+stringfromlist(i,oc)+"_"+stringfromlist(0,tl)



// temperature loop for pulling high temp waves and ratioing them to the 
// lowest; loop starts at 0 for integration with graphervsfreq, but really
//  loop should start at 1 since we don't need to ratio the lowest to itself

for(j=0;j<itemsinlist(tl);j+=1)
setdatafolder $pathstr

string hightempstr = sam+"_"+stringfromlist(i,oc)+"_"+stringfromlist(j,tl)
string hightempratiostr =  sam+"_"+stringfromlist(i,oc)+"tempratio_"+stringfromlist(j,tl)

wave lowtempwave = $lowtempstr
wave hightempwave = $hightempstr


setdatafolder $ratiopath
make/n=(numpnts(lowtempwave))/o $hightempratiostr/wave=lowtempratiowave
lowtempratiowave = hightempwave / lowtempwave
setscale/P x leftx(lowtempwave),dimdelta(lowtempwave,0),"", lowtempratiowave


if(cmpstr(stringfromlist(i,oc),"trmag")==0)
variable lindex, uindex
variable fstep = dimdelta(lowtempratiowave,0)
// ceiling and floor are used to choose
// frequencies above and below this frequency

lindex = ceil(lfreq / fstep)
uindex = floor(ufreq / fstep)

deletepoints uindex, inf, lowtempratiowave
deletepoints 0, lindex-1, lowtempratiowave

endif


endfor
endfor

// the graphing code will be standalone since it is optical constant specific



// make global path strings
setdatafolder root:
string/g trmagtempratiopath = "root:tr:mag:ratio:"

graphervsfreq("trmagtempratiograph",ratiopath,sam+"_trmagtempratio_","\Z18Ratio:Transmission","\Z30 T","ColdWarm",3,0,0)
dowindow/hide=0 trmagtempratiograph
setdatafolder root:

end


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  

function powerabsNEW()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


// units of cm^-1

setdatafolder root:
svar sam,tl,kpath, logypath, logxpath
nvar lfreq,ufreq, npnts, i, j
wave freq


string/g alphapath = "root:alpha:"
windowkill("alphafreqgraph")
overwritefolder(alphapath)                                              

for(i=0;i<npnts;i+=1)

setdatafolder $kpath
string kstr = sam+"_k_"+stringfromlist(i,tl)
wave kwave = $kstr

setdatafolder $alphapath
string alphastr = sam+"_alpha_"+stringfromlist(i,tl)
make/o/n=(numpnts(kwave)) $alphastr/wave=alphawave

// .02998 = c [cm / ps]
//for(k=0;k<numpnts(tempkwave);k+=1)
alphawave = 4 * pi * freq * kwave / .02998 
//endfor
setscale/p x leftx(kwave),dimdelta(kwave,0),"", alphawave

endfor

setdatafolder $logypath
nvar alphafreqlogy
setdatafolder $logxpath
nvar alphafreqlogx
setdatafolder root:

graphervsfreq("alphafreqgraph",alphapath,sam+"_alpha_","\Z18\F'symbol'a\F'geneva'\B\M ","\Z30 \F'Symbol'a\Z30\B","ColdWarm",3,alphafreqlogy,alphafreqlogx)

setdatafolder root:

end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 

function alphavstemp()

// calculates the absorption for a given frequency as a function of temperature

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


// enter freq in THz; output will be GHz

setdatafolder root:
svar sam, alphapath, tl,logypath,logxpath
nvar aufreq, alfreq, freqs, i, j, k
wave temperature

variable freqstep = (aufreq-alfreq) / freqs

// i loop: for each frequency, create a wave (fntempwave)
// j loop: loop over each temperature, find the absorption wave for that temperature
// and assign the jth point of fntempwave to be the absorption of the ith frequency at 
// temperature j

// temp here refers to "temperature", not "temporary"

windowkill("alphatempgraph")


string/g alphatemppath = alphapath + "vstemp:"
overwritefolder(alphatemppath)

// these definitions have zeros because we're choosing a basis wave
// to extract other numbers from 
setdatafolder $alphapath 
string alphastr0 = sam+"_alpha_"+stringfromlist(0,tl)
wave alphawave0 = $alphastr0

for(i=0;i<freqs;i+=1)

	setdatafolder $alphatemppath
	string alphavstempstr = sam+"_alphatemp_"+num2str((leftx(alphawave0)+i*freqstep)*1000)
	make/o/n=(itemsinlist(tl)) $alphavstempstr/wave=alphavstempwave
	
		for(j=0;j<itemsinlist(tl);j+=1)
		    setdatafolder $alphapath 
   		    string alphastr = sam+"_alpha_"+stringfromlist(j,tl)
		    wave alphawave = $alphastr
	
		    alphavstempwave[j]=alphawave(freqstep * i + leftx(alphawave0))
	
             endfor

endfor

setdatafolder root:

string alphatempgraphstr = "alphatempgraph"
display/n=$alphatempgraphstr/hide=1

setdatafolder $alphatemppath
variable npnts = itemsinlist(tl)

for(k=0;k<freqs;k+=1)
	 string fntempstring2 = sam+"_alphatemp_"+num2str((leftx(alphawave0)+k*freqstep)*1000)
	wave fntempwave2=$fntempstring2
	 appendtograph/w=alphatempgraph fntempwave2 vs temperature
endfor

setdatafolder $logypath
nvar alphatemplogy
setdatafolder $logxpath
nvar alphatemplogx
setdatafolder root:

graphervstemp("alphatempgraph","\Z23\F'symbol'a\F'geneva'\B\M\Z17 [cm\S-1\M\Z17]","\Z30\F'symbol'a\F'geneva'\B\M","Rainbow",3,alphatemplogy,alphatemplogx)


end



/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  

function losstanNEW()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 



setdatafolder root:
svar sam,tl, eps1path, eps2path
nvar i, npnts

windowkill("losstanfreqgraph")
windowkill("losstanvsfreqimagegraph")

string/g losstanvsfreqpath = "root:eps:losstanvsfreq:"
overwritefolder(losstanvsfreqpath)

for(i=0;i<npnts;i+=1)

setdatafolder $eps1path
string eps1str = sam+"_eps1_"+stringfromlist(i,tl)
wave eps1wave = $ eps1str
setdatafolder $eps2path
string eps2str = sam+"_eps2_"+stringfromlist(i,tl)
wave eps2wave = $ eps2str

setdatafolder $losstanvsfreqpath
string losstanvsfreqstr = sam+"_losstanvsfreq_"+stringfromlist(i,tl)
make/n=(numpnts(eps1wave))/o  $losstanvsfreqstr/wave=losstanvsfreqwave
losstanvsfreqwave = eps2wave/eps1wave

setscale/p x leftx(eps1wave),dimdelta(eps1wave,0),"",losstanvsfreqwave

endfor

graphervsfreq("losstanfreqgraph",losstanvsfreqpath,sam+"_losstanvsfreq_","\Z18Loss Tangent","\Z18Loss Tan","ColdWarm",3,0,0)


end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  

function energylossNEW()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 



setdatafolder root:
svar sam,tl, eps1path, eps2path
nvar i, npnts, ii

windowkill("energylossvsfreqgraph")
windowkill("energylossvsfreqimagegraph")

string/g energylossvsfreqpath = "root:eps:energylossvsfreq:"
overwritefolder(energylossvsfreqpath)

for(i=0;i<npnts;i+=1)

setdatafolder $eps1path
string eps1str = sam+"_eps1_"+stringfromlist(i,tl)
wave eps1wave = $ eps1str
setdatafolder $eps2path
string eps2str = sam+"_eps2_"+stringfromlist(i,tl)
wave eps2wave = $ eps2str

setdatafolder $energylossvsfreqpath
string energylossvsfreqstr = sam+"_energylossvsfreq_"+stringfromlist(i,tl)
make/n=(numpnts(eps1wave))/o/c  $energylossvsfreqstr/wave=energylossvsfreqwave 
energylossvsfreqwave = (eps1wave+sqrt(-1)*eps2wave)^-1
setscale/p x leftx(eps1wave),dimdelta(eps1wave,0),"",energylossvsfreqwave



endfor

graphervsfreq("energylossvsfreqgraph",energylossvsfreqpath,sam+"_energylossvsfreq_","\Z18Energy Loss","\Z18Energy Loss","ColdWarm",3,0,0)


end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  

function auxconstantsNEW()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 



setdatafolder root:
svar sam,tl, npath, kpath, eps1path, eps2path
nvar i, npnts, ii, z0

windowkill("xsvsfreqgraph")
windowkill("rsvsfreqgraph")
windowkill("rmagvsfreqgraph")
windowkill("rphasevsfreqgraph")
windowkill("energylossrevsfreqgraph")
windowkill("energylossimvsfreqgraph")
windowkill("ebphasevsfreqgraph")
windowkill("xsvsfreqimagegraph")
windowkill("rsvsfreqimagegraph")
windowkill("rmagvsfreqimagegraph")
windowkill("rphasevsfreqimagegraph")
windowkill("energylossrevsfreqimagegraph")
windowkill("energylossimvsfreqimagegraph")
windowkill("ebphasevsfreqimagegraph")
windowkill("losstanfreqgraph")
windowkill("losstanvsfreqimagegraph")

string/g zsvsfreqpath = "root:zs:vsfreq:"
overwritefolder(zsvsfreqpath)
string/g xsvsfreqpath = "root:zs:vsfreq:xs:"
overwritefolder(xsvsfreqpath)
string/g rsvsfreqpath = "root:zs:vsfreq:rs:"
overwritefolder(rsvsfreqpath)

string/g rvsfreqpath = "root:r:vsfreq:"
overwritefolder(rvsfreqpath)
string/g rmagvsfreqpath = "root:r:vsfreq:mag:"
overwritefolder(rmagvsfreqpath)
string/g rphasevsfreqpath = "root:r:vsfreq:phase:"
overwritefolder(rphasevsfreqpath)

string/g energylossvsfreqpath = "root:energyloss:vsfreq:"
overwritefolder(energylossvsfreqpath)
string/g energylossrevsfreqpath = "root:energyloss:vsfreq:re:"
overwritefolder(energylossrevsfreqpath)
string/g energylossimvsfreqpath = "root:energyloss:vsfreq:im:"
overwritefolder(energylossimvsfreqpath)

string/g ebphasevsfreqpath = "root:ebphase:"
overwritefolder(ebphasevsfreqpath)

string/g losstanvsfreqpath = "root:eps:losstanvsfreq:"
overwritefolder(losstanvsfreqpath)



for(i=0;i<npnts;i+=1)

setdatafolder $npath
string nstr = sam+"_n_"+stringfromlist(i,tl)
wave nwave = $ nstr
setdatafolder $kpath
string kstr = sam+"_k_"+stringfromlist(i,tl)
wave kwave = $ kstr

setdatafolder $zsvsfreqpath
string zsvsfreqstr = sam+"_zsvsfreq_"+stringfromlist(i,tl)
make/n=(numpnts(nwave))/o/c  $zsvsfreqstr/wave=zsvsfreqwave 
// this is probably not si units...maybe just cgs
zsvsfreqwave = z0*(nwave+sqrt(-1)*kwave)^-1
setscale/p x leftx(nwave),dimdelta(nwave,0),"",zsvsfreqwave

setdatafolder $rsvsfreqpath
string rsvsfreqstr = sam+"_rsvsfreq_"+stringfromlist(i,tl)
make/n=(numpnts(nwave))/o  $rsvsfreqstr/wave=rsvsfreqwave 
rsvsfreqwave = real(zsvsfreqwave)
setscale/p x leftx(nwave),dimdelta(nwave,0),"",rsvsfreqwave

setdatafolder $xsvsfreqpath
string xsvsfreqstr = sam+"_xsvsfreq_"+stringfromlist(i,tl)
make/n=(numpnts(nwave))/o  $xsvsfreqstr/wave=xsvsfreqwave 
xsvsfreqwave = imag(zsvsfreqwave)
setscale/p x leftx(nwave),dimdelta(nwave,0),"",xsvsfreqwave

setdatafolder $ebphasevsfreqpath
string ebphasevsfreqstr = sam+"_ebphasevsfreq_"+stringfromlist(i,tl)
make/n=(numpnts(nwave))/o  $ebphasevsfreqstr/wave=ebphasevsfreqwave 
ebphasevsfreqwave = atan(kwave/nwave)
setscale/p x leftx(nwave),dimdelta(nwave,0),"",ebphasevsfreqwave


setdatafolder $rvsfreqpath
string rvsfreqstr = sam+"_rvsfreq_"+stringfromlist(i,tl)
make/n=(numpnts(nwave))/o/c  $rvsfreqstr/wave=rvsfreqwave 
rvsfreqwave = (1-(nwave+sqrt(-1)*kwave))/(1+(nwave+sqrt(-1)*kwave))
setscale/p x leftx(nwave),dimdelta(nwave,0),"",rvsfreqwave

setdatafolder $rmagvsfreqpath
string rmagvsfreqstr = sam+"_rmagvsfreq_"+stringfromlist(i,tl)
make/n=(numpnts(nwave))/o  $rmagvsfreqstr/wave=rmagvsfreqwave 
rmagvsfreqwave = real(r2polar(rvsfreqwave))
setscale/p x leftx(nwave),dimdelta(nwave,0),"",rmagvsfreqwave

setdatafolder $rphasevsfreqpath
string rphasevsfreqstr = sam+"_rphasevsfreq_"+stringfromlist(i,tl)
make/n=(numpnts(nwave))/o  $rphasevsfreqstr/wave=rphasevsfreqwave 
rphasevsfreqwave = imag(r2polar(rvsfreqwave))
setscale/p x leftx(nwave),dimdelta(nwave,0),"",rphasevsfreqwave



setdatafolder $eps1path
string eps1str = sam+"_eps1_"+stringfromlist(i,tl)
wave eps1wave = $ eps1str
setdatafolder $eps2path
string eps2str = sam+"_eps2_"+stringfromlist(i,tl)
wave eps2wave = $ eps2str

setdatafolder $energylossvsfreqpath
string energylossvsfreqstr = sam+"_energylossvsfreq_"+stringfromlist(i,tl)
make/n=(numpnts(eps1wave))/o/c  $energylossvsfreqstr/wave=energylossvsfreqwave 
energylossvsfreqwave = (eps1wave+sqrt(-1)*eps2wave)^-1
setscale/p x leftx(eps1wave),dimdelta(eps1wave,0),"",energylossvsfreqwave

setdatafolder $energylossrevsfreqpath
string energylossrevsfreqstr = sam+"_energylossrevsfreq_"+stringfromlist(i,tl)
make/n=(numpnts(nwave))/o  $energylossrevsfreqstr/wave=energylossrevsfreqwave 
energylossrevsfreqwave = real(energylossvsfreqwave)
setscale/p x leftx(eps1wave),dimdelta(eps1wave,0),"",energylossrevsfreqwave

setdatafolder $energylossimvsfreqpath
string energylossimvsfreqstr = sam+"_energylossimvsfreq_"+stringfromlist(i,tl)
make/n=(numpnts(nwave))/o  $energylossimvsfreqstr/wave=energylossimvsfreqwave 
energylossimvsfreqwave = imag(energylossvsfreqwave)
setscale/p x leftx(eps1wave),dimdelta(eps1wave,0),"",energylossimvsfreqwave

setdatafolder $losstanvsfreqpath
string losstanvsfreqstr = sam+"_losstanvsfreq_"+stringfromlist(i,tl)
make/n=(numpnts(eps1wave))/o  $losstanvsfreqstr/wave=losstanvsfreqwave
losstanvsfreqwave = eps2wave/eps1wave

setscale/p x leftx(eps1wave),dimdelta(eps1wave,0),"",losstanvsfreqwave




endfor

setdatafolder $energylossrevsfreqpath
graphervsfreq("energylossrevsfreqgraph",energylossrevsfreqpath,sam+"_energylossrevsfreq_","\Z18Re[Energy Loss]","\Z18Re[Energy Loss]","ColdWarm",3,0,0)
setdatafolder $energylossimvsfreqpath
graphervsfreq("energylossimvsfreqgraph",energylossimvsfreqpath,sam+"_energylossimvsfreq_","\Z18Im[Energy Loss]","\Z18Im[Energy Loss]","ColdWarm",3,0,0)
setdatafolder $rsvsfreqpath
graphervsfreq("rsvsfreqgraph",rsvsfreqpath,sam+"_rsvsfreq_","\Z18R\BS","\Z18R\BS","ColdWarm",3,0,0)
setdatafolder $xsvsfreqpath
graphervsfreq("xsvsfreqgraph",xsvsfreqpath,sam+"_xsvsfreq_","\Z18X\BS","\Z18X\BS","ColdWarm",3,0,0)
setdatafolder $ebphasevsfreqpath
graphervsfreq("ebphasevsfreqgraph",ebphasevsfreqpath,sam+"_ebphasevsfreq_","\Z18EB Phase","\Z18EB Phase","ColdWarm",3,0,0)
setdatafolder $rmagvsfreqpath
graphervsfreq("rmagvsfreqgraph",rmagvsfreqpath,sam+"_rmagvsfreq_","\Z18Reflectivity Magnitude","\Z18Reflectivity \r Magnitude","ColdWarm",3,0,0)
setdatafolder $rphasevsfreqpath
graphervsfreq("rphasevsfreqgraph",rphasevsfreqpath,sam+"_rphasevsfreq_","\Z18Reflectivity Phase","\Z18Reflectivity \r Phase","ColdWarm",3,0,0)
setdatafolder $losstanvsfreqpath
graphervsfreq("losstanfreqgraph",losstanvsfreqpath,sam+"_losstanvsfreq_","\Z18Loss Tangent","\Z18Loss Tan","ColdWarm",3,0,0)


setdatafolder root:

end


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  

function bringgraphstofrontNEW()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


setdatafolder root:
nvar i

string visiblegraphs = winlist("*",",","win:1,visible:1")
visiblegraphs = sortlist(visiblegraphs,",",16)

		 for(i=0;i<itemsinlist(visiblegraphs,",");i+=1)
		 dowindow/f $stringfromlist(i,visiblegraphs,",")
		 endfor
		 
end

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  

function alphaimage()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


setdatafolder root:

svar tl, sam, imagecolor, logypath, logxpath
nvar i,freqs, npnts, lfreq,ufreq
wave temperature, freq

imagefreqtempwaves()
setdatafolder root:
wave imagetemp, imagefreq


//reverse imagetemp

string/g alphapath = "root:alpha:"
windowkill("alphaimagegraph")
string alphaliststr = ""

setdatafolder $alphapath

for(i=0;i<npnts;i+=1)
string alphastr = sam + "_alpha_"+stringfromlist(i,tl)
alphaliststr = addlistitem(alphastr,alphaliststr,";")
endfor
concatenate/o alphaliststr, alphaimagematrix
display/hide=0/n=alphaimagegraph;appendimage alphaimagematrix vs {::imagefreq,::imagetemp}
modifyimage alphaimagematrix ctab= {*,*,$imagecolor,0}
setdatafolder root:
reverse/p temperature/d=revtemp 
reverse/p freq/d=revfreq

appendmatrixcontour/f="%.17s=%g"/w=alphaimagegraph alphaimagematrix vs {freq,revtemp}

setdatafolder $logypath
nvar alphaimagelogy
setdatafolder $logxpath
nvar alphaimagelogx
setdatafolder root:
imagegrapher("\Z16\F'Symbol'a\Z16\B",alphaimagelogy,alphaimagelogx)
//Label left "\Z14Temperature [K]"
//Label bottom "\Z14Frequency [THz]"
//TextBox/C/N=text0/F=0/A=MT "Absorption [cm^-1]"
//TextBox/C/N=text0/F=2
//ModifyGraph fSize(left)=14
//ModifyGraph fSize(bottom)=14
//ModifyGraph lsize=2

end




/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  

function imagefreqtempwaves()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


setdatafolder root:
nvar i
wave temperature, freq

make/n=(numpnts(temperature)+1)/o imagetemp
make/n=(numpnts(freq)+1)/o imagefreq

// make the image temp wave as the average of nearby temperature points
// put in 0 and an artificial upper point that a distance of half length between the two
// last points above the last point

for(i=0;i<numpnts(imagetemp);i+=1)

// average with zero by hand
if(i==0)
imagetemp[i]=temperature[0]-(temperature[1]-temperature[0])/2
if(imagetemp[0]<0)
imagetemp[0]=0
endif

// standard temps in middle of range
elseif(i<(numpnts(imagetemp)-1))
imagetemp[i] = (temperature[i-1]+temperature[i])/2

elseif(i==(numpnts(imagetemp)-1))
imagetemp[i] = temperature[i-1]+(temperature[i-1]-temperature[i-2])/2
endif
endfor

reverse imagetemp

variable freqstep = freq[1]-freq[0]
imagefreq[0] = freq[0] -freqstep/2

for(i=1;i<numpnts(imagefreq);i+=1)
imagefreq[i]=imagefreq[0]+freqstep*i
endfor

end


function printtl()

svar tl
nvar i, npnts

print itemsinlist(tl)

for(i=0;i<itemsinlist(tl);i+=1)
print stringfromlist(i,tl)
endfor 

print stringfromlist(npnts-1,tl)

end



/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  

function auximagemaker()

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


// this function handles all three optical constants (eps, index, sigma) in one function,

setdatafolder root:
nvar i, j
svar rsvsfreqpath, xsvsfreqpath, rmagvsfreqpath, rphasevsfreqpath, energylossrevsfreqpath,energylossimvsfreqpath, ebphasevsfreqpath, tl,sam, imagecolor, contourcolor
wave freq, temperature

reverse temperature/d=revtemp

imagefreqtempwaves()
wave imagefreq, imagetemp

string ocpath = "zs:vsfreq:rs:;zs:vsfreq:xs:;r:vsfreq:mag:;r:vsfreq:phase:;energyloss:vsfreq:re:;energyloss:vsfreq:im:;ebphase:;eps:losstanvsfreq:;"
string ocname ="rsvsfreq;xsvsfreq;rmagvsfreq;rphasevsfreq;energylossrevsfreq;energylossimvsfreq;ebphasevsfreq;losstanvsfreq;"
string graphtitle = "R\BS;X\S;Reflectivity\rMagnitude;Reflectivity\rPhase;Re[Energy Loss];Im[Energy Loss];EB Phase;Loss Tan;"








for(i=0;i<itemsinlist(ocpath);i+=1)

string pathstr = "root:"+stringfromlist(i,ocpath)
setdatafolder $pathstr
string imagestr = ""
string imagegraph = stringfromlist(i,ocname)+"imagegraph"
windowkill(imagegraph)


for(j=0;j<itemsinlist(tl);j+=1)
string appendstr = sam +"_"+stringfromlist(i,ocname)+"_"+stringfromlist(j,tl)
imagestr = addlistitem(appendstr,imagestr)
endfor
concatenate/o imagestr, imagematrix
display/hide=1/n=$imagegraph;appendimage imagematrix vs {imagefreq,imagetemp}
appendmatrixcontour/f="%.17s=%g"/w=$imagegraph imagematrix vs {freq,revtemp}
ModifyContour imagematrix ctabLines={*,*,$contourcolor,0}
ModifyImage imagematrix ctab= {*,*,$imagecolor,0}
// set log settings to 0
imagegrapher(stringfromlist(i,graphtitle), 0, 0)
endfor
end


end



/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  

Function AddGraphsToNotebook(nb, graphList, desiredNumGraphsInParagraph, graphicsMode, scaling)

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  

//	Adds graphs from list to notebook.
//	Example:
//		String graphList = "Graph0;Graph1;Graph2;"
//		AddGraphsToNotebook("MyNotebook", graphList, 2, 0, 50)

	String nb								// Name of notebook, e.g., "MyNotebook"
	string graphList							// e.g., "Graph0;Graph1;" or WinList("*", ";", "") for all graphs
	Variable desiredNumGraphsInParagraph	// e.g., 2
	Variable graphicsMode					// e.g., 0; See documentation for Notebook picture= command.
	Variable scaling							// e.g., 50 for 50 percent
 
	DoWindow/F $nb
	if (V_flag == 0)
		NewNotebook /F=1 /N=$nb /W=(10,50,700,500)	
	endif
 
	Variable index
	Variable numGraphsInParagraph
 
	numGraphsInParagraph = 0
	do
		String graphName = StringFromList(index, graphList)
		if (strlen(graphName) == 0)
			break					// All done
		endif
 
		// Add graph
		Notebook $nb, scaling={scaling,scaling}, picture={$graphName, graphicsMode, 1}
 
		if (numGraphsInParagraph < (desiredNumGraphsInParagraph-1))
			Notebook $nb, text="\t"			// Add tab between graphs
		endif
 
		numGraphsInParagraph += 1
		if (numGraphsInParagraph == desiredNumGraphsInParagraph)
			Notebook $nb, text="\r"			// Go to the next paragraph
			Notebook $nb, text="\r"			// Add a blank line between graph rows
			numGraphsInParagraph = 0
		endif
 
		index += 1
	while(1)
End

 
function reportinput()
	string notebookname = "THz Summary"
	variable visorhid = 1, graphsperpara = 1, graphmode = 0, scaling = 100
	Prompt notebookname, "notebook name:"
	Prompt visorhid, "include graphs - vislble (1) / all (0):"
	prompt graphsperpara, "graphs per paragraph:"
	Prompt graphmode, "graph mode : "
	Prompt scaling, "scaling : "
	DOPrompt "report settings", notebookname, visorhid,graphsperpara,graphmode, scaling
end


/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

//	SaveAllGraphsAsGraphicsFiles(pathName, extension, visibleOnly, format, bitmapRes)
//	NOTE: This overwrites existing files with the same name.
//	Example:
//		SaveAllGraphsAsGraphicsFiles("home", ".png", 1, -5, 288)
Function SaveAllGraphsAsGraphicsFiles(extension, visibleOnly, format, bitmapRes)

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	String extension				// Extension for file name, e.g., ".tif", ".pdf", ".png" . . .
	Variable visibleOnly			// 1 to save only visible graphs, 0 to save visible and hidden graphs.
	Variable format				// For SavePICT /E flag. See SavePICT documentation.
	Variable bitmapRes			// For SavePICT /RES flag (e.g., 288 for 288 dpi).
 
	Variable index = 0
		
		string pathname
		newpath/o pathname
	
	do
		String graphName = WinName(index, 1, visibleOnly)
		if (strlen(graphName) == 0)
			break					// All done
		endif
	
	
		
		string layoutname = graphname
		
		  dowindow $layoutname
               if  (V_Flag == 1)
               killwindow $layoutname
               endif
		
	      newlayout/hide=1/p=landscape/n=$layoutname
	      appendlayoutobject/f=0/w=$layoutname graph $graphname 
	      dowindow/f $layoutname
 
		String fileName = graphName + extension
		SavePICT /N=$layoutName /O /P=$pathName /E=(format)/w=(0,0,0,0) /RES=(bitmapRes)
 killwindow $layoutname
		index += 1
	while(1)
End

/////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function deriv()
            
            
            
string wavenamewotemp            
prompt wavenamewotemp, "Wave Name without Temperature:"
doprompt "Wave Base Name", wavenamewotemp
        
            
           string origdf = getdatafolder(1)
           setdatafolder root:
	     svar tl
	     nvar npnts, i,j, alfreq,aufreq
            
            string/g derivpath = "root:deriv" 
            newdatafolderpathNEW(derivpath)
            
string derivs = "derivs" 
            
windowkill(derivs)
display/n=$derivs
            
for(i=0;i<npnts;i+=1)
	
      setdatafolder $origdf
      string origstr = wavenamewotemp +"_"+stringfromlist(i,tl)
      wave origwave = $origstr
      setdatafolder $derivpath
	string derivstr = "d_"+nameofwave(origwave)
	make/o/n=(numpnts(origwave)) $derivstr/wave=derivwave
	
differentiate origwave/d=derivwave	
	
//	if(order == 1)
//	differentiate origwave /d=derivwave
//	else
//	duplicate/o origwave multideriv
//	for(j=1;j<=order;j+=1)
//	differentiate multideriv
//	print j
//	endfor 
//	derivwave = multideriv
//	endif

	appendtograph derivwave


	
endfor
 
string leftaxis =  "\Z14Derivative of " + wavenamewotemp 
Label left leftaxis
Label bottom "\Z14Frequency [THz]"
ModifyGraph mirror=2
TextBox/C/N=text0/F=0/A=LT "\\Z30Deriv"
ModifyGraph fSize(left)=14
ModifyGraph fSize(bottom)=14
ModifyGraph lsize=2
ApplyColorTableToTopGraphNEW("ColdWarm")
setaxis bottom alfreq,aufreq

setdatafolder root:

end