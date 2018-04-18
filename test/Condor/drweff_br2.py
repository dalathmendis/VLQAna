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
#leg = TLegend(0.72,0.91,0.92,0.70)
leg = TLegend(0.57,0.89,0.89,0.65)   
leg.SetBorderSize(0)
leg.SetFillColor(10)
leg.SetLineColor(10)
leg.SetLineWidth(0)
 
def setCosmetics(hist, legname, hname, color, var):
    if var == 'st': Var = 'S_{T}'
    elif var == 'chi_mass': Var = 'M_{#Chi^{2}}'
    legname.split('_',0)[0]
    hist.SetLineColor(color)
    #hist.GetYaxis().SetTitle('Efficieny ( '+Var+' ) %')
    hist.GetYaxis().SetTitle('Efficieny (%)') 
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetTitleFont(42)
    hist.GetXaxis().SetTitleOffset(1.2)
    hist.GetXaxis().SetTitleFont(42)
    hist.SetTitle('')
    hist.SetName(hname)
    hist.SetLineWidth(2)
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
                  default='sigeffnewel/14.01.18',
                  dest='plotDir',
                  help='output directory of plots')
(options,args) = parser.parse_args()
# ==========end: options =============
var = options.var
lep = options.lep
plotDir = options.plotDir
title = ''
if lep == 'ele':
    title = 'e^{#pm}e^{#mp}+jets'
elif lep == 'mu':
    title = '#mu^{#pm}#mu^{#mp}+jets'
else:    
    print 'exiting: choose the lep string either as ele or mu'
    exit()
 
#dirIn = '/uscms_data/d3/dmendis/80x/CMSSW_8_0_15/src/Analysis/VLQAna/test/Condor/Histo12MU4/'
#d_BB = '/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/bprime/el/'
d_BB = '/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/bprime/mu/'
#d_TT = dirIn+'TT/'+lep+'/'
#d_BB ='/uscms_data/d3/dmendis/80x/CMSSW_8_0_15/src/Analysis/VLQAna/test/Condor/signal/withcut/'
#d_TT ='/uscms_data/d3/dmendis/80xB2Gana3/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Histo5p7EL/'
#d_TT = '/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Histo3MUb/'
#d_TT='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/runonskims/batch2/newera/Histo3ELsysnody/'
d_TT='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/runonskims/batch2/newera/Histo2MUsysnodypdf/'
#d_TT='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/runonskims/batch2/newera/Ele32_1/'

#mass = [800 , 1000 , 1200 ,1400 ,1600 ,1800]
mass = [900 , 1100 , 1300 ,1500 ,1700 ] 
#channel = ['BpBp_bZbH', 'BpBp_bZbZ', 'TpTp_tZtH', 'TpTp_tZtZ']
channel = ['TpTp_tZtH', 'TpTp_tZtZ','TpTp_tZbW','BpBp_bZbH', 'BpBp_bZbZ']#, 'BpBp_bZtW'] 
channel_m = ['TpTp_tZtH', 'TpTp_tZtZ','TpTp_tZbW', 'TpTp_tHtH','TpTp_bWbW']
channel_m2 = ['BpBp_bZbH', 'BpBp_bZbZ','BpBp_bZtW', 'BpBp_bHbH','BpBp_tWtW']
lepdir = ['ele', 'mu']
templates = []
 
#hEff = TH1D('', 'Signal Efficiency; VLQ Mass [GeV]; efficiency (%)', 13, 650, 1950)
hEff = TH1D('', 'Signal Efficiency; VLQ Mass [GeV]; efficiency (%)', 11, 750, 1850) 
ibin = 0; icol = 0
for ch in channel:
    h_eff = hEff.Clone()
    h_eff.SetDirectory(0)
    for m in mass:
      #  if 'BpBp' in ch and m == '1600': continue
        # open the right file 
        if 'BpBp' in ch:
            #f = TFile(d_BB+ch+'_M-'+str(m)+'.root')
            if ch == 'BpBp_bZbZ':
                f = TFile(d_BB+'bprime'+str(m)+'_'+ch[5:]+'.root')

            elif ch == 'BpBp_bZbH':
                f  = TFile(d_BB+'bprime'+str(m)+'_'+channel_m2[0][5:]+'.root')
                f2 = TFile(d_BB+'bprime'+str(m)+'_'+channel_m2[1][5:]+'.root')
                f3 = TFile(d_BB+'bprime'+str(m)+'_'+channel_m2[3][5:]+'.root')

            elif ch == 'BpBp_bZtW':
                f  = TFile(d_BB+'bprime'+str(m)+'_'+channel_m2[2][5:]+'.root')
                f2 = TFile(d_BB+'bprime'+str(m)+'_'+channel_m2[1][5:]+'.root')
                f3 = TFile(d_BB+'bprime'+str(m)+'_'+channel_m2[4][5:]+'.root')

        elif 'TpTp' in ch:
            if ch == 'TpTp_tZtZ':
                f = TFile(d_TT+ch+'_M-'+str(m)+'.root')

            elif ch == 'TpTp_tZtH':
                f  = TFile(d_TT+channel_m[0]+'_M-'+str(m)+'.root')
                f2 = TFile(d_TT+channel_m[1]+'_M-'+str(m)+'.root')
                f3 = TFile(d_TT+channel_m[3]+'_M-'+str(m)+'.root')

            elif ch == 'TpTp_tZbW':
                f  = TFile(d_TT+channel_m[2]+'_M-'+str(m)+'.root')
                f2 = TFile(d_TT+channel_m[1]+'_M-'+str(m)+'.root')
                f3 = TFile(d_TT+channel_m[4]+'_M-'+str(m)+'.root')
        #print f ,",",f2, ",",f3
        # get the histogram and normalization
        if ch == 'TpTp_tZtZ' or 'BpBp_bZbZ':    
            #h1 = f.Get('ana/sig/'+var)
            h1 = f.Get('ana/cutflow')
            h1_norm = f.Get('ana/signalEvts')
        else:
           # h1 = f.Get('ana/sig/'+var)
           # h2 = f2.Get('ana/sig/'+var)
           # h3 = f3.Get('ana/sig/'+var)
            h1 = f.Get('ana/cutflow')
            h2 = f2.Get('ana/cutflow')
            h3 = f3.Get('ana/cutflow')

        
            h1.Scale(0.5)
            h2.Scale(0.25)
            h3.Scale(0.25)
            h1.Add(h2,1)
            h1.Add(h3,1)
        #h1 = f.Get('ana/pre/'+var+'_pre')
            h1_norm  = f.Get('ana/signalEvts')
            h1_norm2 = f2.Get('ana/signalEvts')
            h1_norm3 = f3.Get('ana/signalEvts')
            h1_norm.Scale(0.5)
            h1_norm2.Scale(0.25)
            h1_norm3.Scale(0.25)
            h1_norm.Add(h1_norm2,1)
            h1_norm.Add(h1_norm3,1)

       # print h1.GetName() ,",",h2, ",",h3
        nGen = h1_norm.GetBinContent(1)
#        sig = h1.GetEntries()
      #  sig = h1.Integral()
      #  print "bin label = ", h1.GetXaxis().GetBinLabel(2)
        sig = h1.GetBinContent(2)
        #h1.Draw()
        #raw_input('hold')
        # compute eff and uncertainy per mass point
        eff = sig/nGen
        unc = binomialUnc(eff, nGen)
         
       # print 'channel: {0:<5}, Mass: {1:<5.0f}, Signal: {2:<5.0f}, NGen = {3:<5.0f}, Eff(%) = {4:<5.2f}+/-{5:<5.2f}'.format(ch, m, sig, nGen, eff*100, unc*100)
        print ' {0:<5}-M{1:<5.0f}& {2:<5.0f}& {3:<5.0f} & {4:<5.2f} $\pm$ {5:<5.2f} \\\\'.format(ch[5:], m, sig, nGen, eff*100, unc*100)
        # ibin += m%2 + 1
        ibin +=  2  
        # store the efficiencies to new historgram
        h_eff.SetBinContent(ibin, eff*100)
        h_eff.SetBinError(ibin, unc*100)
 
    # make cosmetic changes to the histogram and iterate to next one
    if ch == 'TpTp_tZtZ':
        leg_name= "#bf{#it{#Beta}(T #rightarrow tZ)} = 1.0"
    elif ch == 'TpTp_tZtH':
        leg_name = "#bf{#it{#Beta}(T #rightarrow tZ)}  = #bf{#it{#Beta}(T #rightarrow tH)}  = 0.5"
    elif ch == 'TpTp_tZbW':
        leg_name ="#bf{#it{#Beta}(T #rightarrow tZ)}  = #bf{#it{#Beta}(T #rightarrow bW)}  = 0.5"

    elif ch == 'BpBp_bZbZ':
        leg_name= "#bf{#it{#Beta}(B #rightarrow bZ)} = 1.0"

    elif ch == 'BpBp_bZbH':
        leg_name = "#bf{#it{#Beta}(B #rightarrow bZ)}  = #bf{#it{#Beta}(B #rightarrow bH)}  = 0.5"
    elif ch == 'BpBp_bZtW':
        leg_name ="#bf{#it{#Beta}(T #rightarrow bZ)}  = #bf{#it{#Beta}(T #rightarrow tW)}  = 0.5"

    if 'BpBp' in ch:
        setCosmetics(h_eff, leg_name, var+'_Eff_'+ch, kOrange+icol*8, var)
    elif 'TpTp' in ch:
        setCosmetics(h_eff, leg_name, var+'_Eff_'+ch, kBlue+icol*4, var)
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
    h.SetMaximum(40.0);
    h.Draw("L same")
leg.Draw()
     
prel = TLatex()
prel.SetNDC(kTRUE)
prel.SetTextFont(52)
prel.SetTextSize(0.04)
#prel.DrawLatex(0.20,0.93,"Simulation HLT_Ele23_Ele12 (dR>0.3)")
if lep == 'mu':
    #title1= "HLT_IsoMu24 or HLT_IsoTkMu24"
    title1= "HLT_IsoMu24 or HLT_IsoTkMu24"  
elif  lep == 'ele':

    title1= "HLT_Ele115_CaloIdVT_GsfTrkIdT"
   # title1= "HLT_Ele32_eta2p1_WPTight_Gsf"
prel.DrawLatex(0.11,0.85,title1)

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
chan.DrawLatex(0.45, 0.80, title)
 
c1.SaveAs(plotDir+'/'+var+'mu_eff2_'+lep+'.pdf')
raw_input('---')
