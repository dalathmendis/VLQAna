from ROOT import TH1D,TCanvas,TH1

from array import array
from ROOT import TH1D,TFile,TCanvas, THStack,TF1, TH1,TLegend,kRed,kBlue,TPad,gPad,TLine,kBlack,TMath ,TGraph,TMultiGraph,TLatex,kGreen
from ROOT import gROOT,gStyle,gPad,gStyle,kTRUE

gStyle.SetOptStat(0)

path1 = '../'
path2 = '/uscms_data/d3/dmendis/80x/CMSSW_8_0_15/src/Analysis/VLQAna/test/Condor/Macro/mutest2/CR_Zmumu/'


leg = TLegend(0.75,0.75,0.95,0.92)
leg.SetBorderSize(1)
leg.SetFillColor(10)

f1 =TFile(path1+'cutflow4.root')

def binomialUnc(eff, Ngen):
    unc = TMath.Sqrt(eff*(1-eff)/Ngen)
    return unc

def error(B,S,dB,dS,SOB):
    p =TMath.Sqrt(S+B)
    e1=TMath.Power((1/p - (S*p/2)),2)
    e2=TMath.Power( (S*p/2),2)
    e =TMath.Sqrt(e1*dS*dS+e2*dB*dB)
    return e*SOB

def setCosmetics(hist, legname, hname, color):
    #legname.split('_',0)[0]
    hist.SetLineColor(color)
    hist.GetYaxis().SetTitle('S/sqrt(S+B)')
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetTitleFont(42)
    hist.GetXaxis().SetTitleOffset(1.2)
    hist.GetXaxis().SetTitleFont(42)
    hist.SetTitle('')
    hist.SetName(hname)
    hist.SetLineWidth(2)
    leg.AddEntry(hist, legname, 'l')



hSig = TH1D('', 'Significance; Categories; S/sqrt(S+B) ', 16, 0.5, 16.5)
hSig.GetXaxis().SetBinLabel(1, "T1Z1H1b1 ") 
hSig.GetXaxis().SetBinLabel(2, "T1Z1H1b2") 
hSig.GetXaxis().SetBinLabel(3, "T1Z1H0b1 ") 
hSig.GetXaxis().SetBinLabel(4, "T1Z1H0b2  ") 
hSig.GetXaxis().SetBinLabel(5, "T0Z1H1b1 ") 
hSig.GetXaxis().SetBinLabel(6, "T0Z1H1b2 ") 
hSig.GetXaxis().SetBinLabel(7, "T0Z1H0b1") 
hSig.GetXaxis().SetBinLabel(8, "T0Z1H0b2 " ) 
hSig.GetXaxis().SetBinLabel(9, "T1Z0H1b1 " ) 
hSig.GetXaxis().SetBinLabel(10, "T1Z0H1b2 " ) 
hSig.GetXaxis().SetBinLabel(11, "T1Z0H0b1 " ) 
hSig.GetXaxis().SetBinLabel(12, "T1Z0H0b2 " ) 
hSig.GetXaxis().SetBinLabel(13, "T0Z0H1b1 " ) 
hSig.GetXaxis().SetBinLabel(14, "T0Z0H1b2 " ) 
hSig.GetXaxis().SetBinLabel(15, "T0Z0H0b1 " ) 
hSig.GetXaxis().SetBinLabel(16, "T0Z0H0b2 " ) 




h1= TH1D ()
h2 = TH1D()
h1= f1.Get('mumuT0Z0H0b2__DY')
h2= f1.Get('mumuT0Z0H0b2__Top')
h1.Add(h2,1)

mass = [800]
channel = ['TT_tZtZ', 'TT_tZtH', 'TT_tZbW']
templates = []

Nbins = h1.GetNbinsX()

print "n bins ", Nbins 
icol=0
for ch in channel:
    h_eff = hSig.Clone()
    h_eff.SetDirectory(0)
#    icol=0
    for m in mass:
        h_eff.Reset()
        # open the right file 
        h_sig = f1.Get('mumuT0Z0H0b2__'+ch+'_M'+str(m))
        optsOverb = 0
        optVal = 0
        optBin = 0
        optSig = 0
        optBkg = 0

        # get the histogram and normalization    
        for ibin in range(1,Nbins+1):
            #bkg = h1.Integral(ibin, ibin+1)
            #sig = h_sig.Integral(ibin, ibin+1)
            bkg = h1.GetBinContent(ibin)                                                                                                                                 
            sig = h_sig.GetBinContent(ibin)
            xAxis = h1.GetBinCenter(ibin)
            if bkg+sig != 0. :
                sOverb = sig/TMath.Sqrt(bkg+sig)#+0.5*bkg)                                                                                                                     
                #unc = error()
            
            else: sOverb = 0.
         
            optsOverb = sOverb; optSig = sig; optBkg = bkg; optVal = xAxis; optBin = ibin
#        print "channel",  " bin "   , " S/sqrt(S+B " , " Signal Events ", "Bkg Events " 
           # print '{0:<5}  &{1:<5.1f} & {2:<5.1f} & {3:<5.1f} & {4:<5.1f}  & {5:<5.1f}  \\\\\ '.format(ch+str(m),optVal, optsOverb, optSig, optBkg,h_eff.GetBinError(ibin))
        
 
            h_eff.SetBinContent(optBin, optsOverb)
            print '{0:<5}  &{1:<5.1f} & {2:<5.1f} & {3:<5.1f} & {4:<5.1f}  & {5:<5.1f}  \\\\\ '.format(ch+str(m),optVal, optsOverb, optSig, optBkg,h_eff.GetBinError(ibin))
            
            unc =error(optBkg,optSig,h1.GetBinError(optBin),h_sig.GetBinError(optBin),optsOverb)
            h_eff.SetBinError(optBin, unc)
            print " error ", unc
    # make cosmetic changes to the histogram and iterate to next one
            #icol = icol+1
    icol = icol+1
    setCosmetics(h_eff, ch+"_"+str(m),ch+"_"+str(m), kBlack+icol)
    #icol = icol+1

    templates.append(h_eff)


c1 = TCanvas('c1', 'c1', 800, 600)
for h in templates :
    print h.GetName()
    h.SetMaximum(2.5);
    h.SetMinimum(0.0);
    #h.Draw("L same")
    #h.Draw("E same") 
    h.Draw("hist same")
leg.Draw()

prel = TLatex()
prel.SetNDC(kTRUE)
prel.SetTextFont(52)
prel.SetTextSize(0.05)
prel.DrawLatex(0.34,0.93,"Siginificance of Categories")

cms = TLatex()
cms.SetNDC(kTRUE)
cms.SetTextFont(61)
cms.SetTextSize(0.05)
cms.DrawLatex(0.10,0.93,"CMS(Preliminary)")
leg.Draw()

ll = TLatex()
ll.SetNDC(kTRUE)
ll.SetTextSize(0.05)
ll.DrawLatex(0.64,0.93, "12.9 fb^{-1}(13 TeV)");
raw_input("hold on")




