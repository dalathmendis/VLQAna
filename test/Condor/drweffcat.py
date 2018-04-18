#! /usr/bin/env python
import sys
from ROOT import TH1D,TFile,TMath,TCanvas,TLegend,TLatex,TLine
from ROOT import gROOT,gStyle,gPad,gStyle
from ROOT import Double,kBlue,kRed,kOrange,kMagenta,kYellow,kCyan,kGreen,kGray,kTRUE
 
gROOT.Macro("~/rootlogon.C")
gStyle.SetOptStat(0)
 
def overUnderFlow(hist):
    xbins = hist.GetNbinsX()
    hist.SetBinContent(xbins, hist.GetBinContent(xbins)+hist.GetBinContent(xbins+1))
    hist.SetBinContent(1, hist.GetBinContent(0)+hist.GetBinContent(1))
    hist.SetBinError(xbins, TMath.Sqrt(TMath.Power(hist.GetBinError(xbins),2)+TMath.Power(hist.GetBinError(xbins+1),2)))
    hist.SetBinError(1, TMath.Sqrt(TMath.Power(hist.GetBinError(0),2)+TMath.Power(hist.GetBinError(1),2)))
    hist.SetBinContent(xbins+1, 0.)
    hist.SetBinContent(0, 0.)
    hist.SetBinError(xbins+1, 0.)
    hist.SetBinError(0, 0.)
 
def binomialUnc(eff, Ngen):
    unc = TMath.Sqrt(eff*(1-eff)/Ngen)
    return unc
 
# Legend
leg = TLegend(0.72,0.91,0.92,0.70)
leg.SetBorderSize(0)
leg.SetFillColor(10)
leg.SetLineColor(10)
leg.SetLineWidth(0)
 
def setCosmetics(hist, legname, hname, color, var):
    if var == 'st': Var = 'S_{T}'
    elif var == 'chi_mass': Var = 'M_{#Chi^{2}}'
    legname.split('_',0)[0]
    hist.SetLineColor(color)
   # hist.GetYaxis().SetTitle('Efficieny ( '+Var+' ) %')
    hist.GetYaxis().SetTitle('Efficieny ('+Category+'/ Tot Sig Events) %')
   
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetTitleFont(42)
    hist.GetXaxis().SetTitleOffset(1.2)
    hist.GetXaxis().SetTitleFont(42)
    hist.SetTitle('')
    hist.SetName(hname)
    hist.SetLineWidth(2)
   # hist.GetYaxis().SetRangeUser(0,60)
    leg.AddEntry(hist, legname, 'l')
     
# ===============
# options
# ===============
from optparse import OptionParser
parser = OptionParser()
 
parser.add_option('--var', metavar='V', type='string', action='store',
                  default='st',#'chi_mass''st'
                  dest='var',
                  help='variable to plot')
parser.add_option('--lep', metavar='L', type='string', action='store',
                  default='mu',
                  dest='lep',
                  help='ele or mu'),
parser.add_option('--plotDir', metavar='P', type='string', action='store',
                  default='efficiencynew',
                  dest='plotDir',
                  help='output directory of plots')

parser.add_option('--Cat', metavar='P', type='string', action='store',
                  default='CatB',
                  dest='category',
                  help='name of category')



(options,args) = parser.parse_args()
# ==========end: options =============
var = options.var
lep = options.lep
plotDir = options.plotDir
Category = options.category

title = ''
if lep == 'ele':
    title = 'e^{#pm}e^{#mp}+jets'
elif lep == 'mu':
    title = '#mu^{#pm}#mu^{#mp}+jets'
else:    
    print 'exiting: choose the lep string either as ele or mu'
    exit()
 
#dirIn = '/uscms_data/d3/dmendis/80x/CMSSW_8_0_15/src/Analysis/VLQAna/test/Condor/Histo12MU4/'
#d_BB = dirIn+'BB/'+lep+'/'
#d_TT = dirIn+'TT/'+lep+'/'
#d_BB ='/uscms_data/d3/dmendis/80x/CMSSW_8_0_15/src/Analysis/VLQAna/test/Condor/signal/withcut/'
if lep=='mu' : name='Histo3MUb/'
else: name='Histo3ELb/'
#d_TT ='/uscms_data/d3/dmendis/80xB2Gana6/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/'+name
d_TT ='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/'+name  
print "path =", d_TT
mass = [800 ,900, 1000,1100, 1200,1300,1400,1500,1600,1700,1800]
#channel = ['BpBp_bZbH', 'BpBp_bZbZ', 'TpTp_tZtH', 'TpTp_tZtZ']
channel = ['TpTp_tZtH', 'TpTp_tZtZ','TpTp_tZbW'] 
lepdir = ['ele', 'mu']
templates = []
 
hEff = TH1D('', 'Signal Efficiency; T Mass [GeV]; efficiency (%)', 11, 750, 1850)
 
ibin = 0; icol = 0
for ch in channel:
    h_eff = hEff.Clone()
    h_eff.SetDirectory(0)
    for m in mass:
         
        # open the right file 
        if 'BpBp' in ch:
            f = TFile(d_BB+ch+'_M-'+str(m)+'.root')
        elif 'TpTp' in ch:
            f = TFile(d_TT+ch+'_M-'+str(m)+'.root')
        h2=TH1D()             
        h3=TH1D()
        h4=TH1D()
        h5=TH1D()
        h6=TH1D()

        # get the histogram and normalization    
        h1 = f.Get('ana/sig/'+var)
        h2 = f.Get('ana/cat/ST_sigT1Z1H1b1')
        h3 = f.Get('ana/cat/ST_sigT1Z1H1b2')
    
        h4 = f.Get('ana/cat/ST_sigT1Z1H0b1')
        h5 = f.Get('ana/cat/ST_sigT1Z1H0b2')
        h6 = f.Get('ana/cat/ST_sigT0Z1H1b2')
        #h2.Add(h3,1)
        h4.Add(h5,1) 
        h4.Add(h6,1)
        #h1 = f.Get('ana/pre/'+var+'_pre')
        h1_norm = f.Get('ana/signalEvts')
        nGen = h1_norm.GetBinContent(1)
        sig = h1.GetEntries()
        sigcat1=h2.GetEntries()
        sigcat2=h3.GetEntries()
        sigcat3=h4.GetEntries()    
  #print "h3 entries", h3.GetEntries()       
        # print " entries", h2.GetEntries()
 #h1.Draw()
        #raw_input('hold')
        # compute eff and uncertainy per mass point
       # eff = sig/nGen
       # print sigcat1 , "/",sig
        if Category == "CatA": eff = sigcat1/sig
        elif Category == "CatB": eff = sigcat2/sig
        else :eff = sigcat3/sig
       # unc = binomialUnc(eff, nGen)
        unc = binomialUnc(eff, sig)
 
      #  print 'channel: {0:<5}, Mass: {1:<5.0f}, Signal: {2:<5.0f}, NGen = {3:<5.0f}, Eff(%) = {4:<5.2f}+/-{5:<5.2f}'.format(ch, m, sig, nGen, eff*100, unc*100)
        print 'channel: {0:<5}, Mass: {1:<5.0f}, Signal: {2:<5.0f}, NGen = {3:<5.0f}, Eff(%) = {4:<5.2f}+/-{5:<5.2f}'.format(ch, m, sigcat3, sig, eff*100, unc*100)
        ibin += m%2+1
 
        # store the efficiencies to new historgram
        h_eff.SetBinContent(ibin, eff*100)
        h_eff.SetBinError(ibin, unc*100)
 
    # make cosmetic changes to the histogram and iterate to next one
    setCosmetics(h_eff, ch, var+'_Eff_'+ch, kGreen+icol*4, var)
    templates.append(h_eff)
    ibin = 0
    icol = icol+1
      
#print templates        
 
# ==============================
# plot
# ==============================
c1 = TCanvas('c1', 'c1', 800, 600)
 
for h in templates :
    print h.GetName()
    #h.SetMaximum(100.0);
    #h.SetMinimum(25.0);
    h.Draw("L same")
    h.GetYaxis().SetRangeUser(0,60)
leg.Draw()
     
prel = TLatex()
prel.SetNDC(kTRUE)
prel.SetTextFont(52)
prel.SetTextSize(0.05)
prel.DrawLatex(0.19,0.93,"Simulation Efficiency("+Category+")")
 
cms = TLatex()
cms.SetNDC(kTRUE)
cms.SetTextFont(61)
cms.SetTextSize(0.05)
cms.DrawLatex(0.10,0.93,"CMS")
leg.Draw()
 
ll = TLatex()
ll.SetNDC(kTRUE)
ll.SetTextSize(0.05)
ll.DrawLatex(0.60,0.93, "35.9 fb^{-1}(13 TeV) (2016)");
 
chan = TLatex()
chan.SetNDC(kTRUE)
chan.SetTextSize(0.05)
chan.DrawLatex(0.20, 0.84, title)
 
#c1.SaveAs(plotDir+'/'+Category+lep+'.pdf')
raw_input('---')
