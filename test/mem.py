import ROOT, json, sys
ROOT.gSystem.Load("libTTHCommonClassifier")
ROOT.gSystem.Load("libTTHMEIntegratorStandalone")
CvectorLorentz = getattr(ROOT, "std::vector<TLorentzVector>")
Cvectordouble = getattr(ROOT, "std::vector<double>")
CvectorJetType = getattr(ROOT, "std::vector<int>")

ncalls = [1000, 2000, 4000, 8000, 16000, 32000]
f = ROOT.MEMClassifier()

if len(sys.argv)==2:
    inf = open(sys.argv[1])
else:
    inf = sys.stdin

def make_p4(pt, eta, phi, m):
    v = ROOT.TLorentzVector()
    v.SetPtEtaPhiM(pt, eta, phi, m)
    return v

for ev in inf.readlines():
    print "----"
    print ev
    jsev = json.loads(ev)
    jets_p4 = jsev["input"]["selectedJetsP4"]
    jets_csv = jsev["input"]["selectedJetsCSV"]
    jets_btag = jsev["input"]["selectedJetsBTag"]

    c_jets_p4 = CvectorLorentz()
    c_jets_csv = Cvectordouble()
    c_jets_jettype = CvectorJetType()
    for ij, j in enumerate(jets_p4):
        c_jets_p4.push_back(make_p4(*j))
        c_jets_csv.push_back(jets_csv[ij])
        c_jets_jettype.push_back(f.RESOLVED)

    c_leps_p4 = CvectorLorentz()
    c_leps_charge = Cvectordouble()
    leps_p4 = jsev["input"]["selectedLeptonsP4"]
    leps_charge = jsev["input"]["selectedLeptonsCharge"]
    for il, l in enumerate(leps_p4):
        c_leps_p4.push_back(make_p4(*l))
        c_leps_charge.push_back(leps_charge[il])

    met_p4 = make_p4(jsev["input"]["metP4"][0], 0.0, jsev["input"]["metP4"][1], 0.0)

    c_loosejets_p4 = CvectorLorentz()
    c_loosejets_csv = Cvectordouble()
    out = {}
    print "tthbb13 code cat={0}".format(jsev["event"]["cat"])
    for _ncalls in ncalls:
        if jsev["output"]["mem_cfg"] == "SL_2w2h2t":
            ret = f.GetOutput(f.SL_2W2H2T, c_leps_p4, c_leps_charge, c_jets_p4, c_jets_csv, c_jets_jettype, c_loosejets_p4, c_loosejets_csv, met_p4, _ncalls)
        elif jsev["output"]["mem_cfg"] == "SL_1w2h2t":
            ret = f.GetOutput(f.SL_1W2H2T, c_leps_p4, c_leps_charge, c_jets_p4, c_jets_csv, c_jets_jettype, c_loosejets_p4, c_loosejets_csv, met_p4, _ncalls)
        elif jsev["output"]["mem_cfg"] == "SL_0w2h2t":
            ret = f.GetOutput(f.SL_0W2H2T, c_leps_p4, c_leps_charge, c_jets_p4, c_jets_csv, c_jets_jettype, c_loosejets_p4, c_loosejets_csv, met_p4, _ncalls)
        out["output__{0}".format(_ncalls)] = {
            "p_sig": ret.p_sig,
            "p_err_sig": ret.p_err_sig,
            "p_bkg": ret.p_bkg,
            "p_err_bkg": ret.p_err_bkg
        }
    print "out:", json.dumps(out)

