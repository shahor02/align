{
//=========Macro generated from canvas: c1/c1
//=========  (Wed Jun 17 15:14:11 2015) by ROOT version5.34/08
   TCanvas *c1 = new TCanvas("c1", "c1",49,52,879,905);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c1->ToggleEventStatus();
   c1->SetHighLightColor(2);
   c1->Range(-617.6471,-650,558.8235,600);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetGridx();
   c1->SetGridy();
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.08);
   c1->SetBottomMargin(0.12);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   TH2F *hh = new TH2F("hh","",100,-500,500,100,-500,500);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   hh->SetLineColor(ci);
   hh->SetMarkerColor(4);
   hh->SetMarkerStyle(21);
   hh->GetXaxis()->SetTitle("X,cm");
   hh->GetXaxis()->SetLabelFont(42);
   hh->GetXaxis()->SetLabelOffset(0.03);
   hh->GetXaxis()->SetLabelSize(0.03);
   hh->GetXaxis()->SetTitleSize(0.035);
   hh->GetXaxis()->SetTickLength(-0.03);
   hh->GetXaxis()->SetTitleOffset(1.35);
   hh->GetXaxis()->SetTitleFont(42);
   hh->GetYaxis()->SetTitle("Y,cm");
   hh->GetYaxis()->CenterTitle(true);
   hh->GetYaxis()->SetLabelFont(42);
   hh->GetYaxis()->SetLabelOffset(0.02);
   hh->GetYaxis()->SetLabelSize(0.03);
   hh->GetYaxis()->SetTitleSize(0.035);
   hh->GetYaxis()->SetTickLength(-0.02);
   hh->GetYaxis()->SetTitleOffset(1.4);
   hh->GetYaxis()->SetTitleFont(42);
   hh->GetZaxis()->SetLabelFont(42);
   hh->GetZaxis()->SetLabelOffset(0.02);
   hh->GetZaxis()->SetLabelSize(0.035);
   hh->GetZaxis()->SetTitleSize(0.035);
   hh->GetZaxis()->SetTickLength(-0.02);
   hh->GetZaxis()->SetTitleFont(42);
   hh->Draw("");
   
   TGraph *graph = new TGraph(15);
   graph->SetName("Graph");
   graph->SetTitle("Graph");

   ci = TColor::GetColor("#000099");
   graph->SetLineColor(ci);
   graph->SetMarkerColor(2);
   graph->SetMarkerStyle(20);
   graph->SetPoint(0,419.6500072,28.7000041);
   graph->SetPoint(1,346.0999935,155.1499982);
   graph->SetPoint(2,253.5999998,263.3999938);
   graph->SetPoint(3,149.0500065,319.9999926);
   graph->SetPoint(4,29.29997136,327.8999963);
   graph->SetPoint(5,-82.75000477,332.9999936);
   graph->SetPoint(6,-195.350001,258.8999951);
   graph->SetPoint(7,-277.8999953,145.6999977);
   graph->SetPoint(8,-395.3500057,30.80000389);
   graph->SetPoint(9,-399.4000057,-100.700004);
   graph->SetPoint(10,-331.8199939,-220.599996);
   graph->SetPoint(11,-230.8499999,-322.7999949);
   graph->SetPoint(12,-128.6700061,-405.5999932);
   graph->SetPoint(13,355.7499944,-203.0999963);
   graph->SetPoint(14,411.5000049,-93.45000196);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,-481.305,501.555);
   Graph_Graph1->SetMinimum(-479.46);
   Graph_Graph1->SetMaximum(406.86);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1->SetLineColor(ci);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelOffset(0.03);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetXaxis()->SetTickLength(-0.03);
   Graph_Graph1->GetXaxis()->SetTitleOffset(1.35);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelOffset(0.02);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetYaxis()->SetTickLength(-0.02);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelOffset(0.02);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetZaxis()->SetTickLength(-0.02);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph1);
   
   graph->Draw("p ");
   
   graph = new TGraph(18);
   graph->SetName("Graph");
   graph->SetTitle("Graph");

   ci = TColor::GetColor("#000099");
   graph->SetLineColor(ci);
   graph->SetMarkerColor(4);
   graph->SetMarkerStyle(21);
   graph->SetPoint(0,444.6000085,29.29999793);
   graph->SetPoint(1,332.1750122,171.1499944);
   graph->SetPoint(2,243.735003,263.9999931);
   graph->SetPoint(3,154.7999982,344.5999998);
   graph->SetPoint(4,38.89996996,332.2499871);
   graph->SetPoint(5,-75.1999948,346.7499997);
   graph->SetPoint(6,-178.300005,270.9999939);
   graph->SetPoint(7,-253.7000108,165.399994);
   graph->SetPoint(8,-374.165006,63.3899985);
   graph->SetPoint(9,-351.1000068,-47.29999864);
   graph->SetPoint(10,-305.6500125,-158.6499929);
   graph->SetPoint(11,-248.9000031,-293.9049938);
   graph->SetPoint(12,-111.0999975,-366.5000002);
   graph->SetPoint(13,-2.924999967,-353.1999871);
   graph->SetPoint(14,144.8999975,-355.66);
   graph->SetPoint(15,253.9500033,-270.7999935);
   graph->SetPoint(16,347.2500119,-201.5999943);
   graph->SetPoint(17,418.0500076,-113.5999974);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,-456.0415,526.4765);
   Graph_Graph2->SetMinimum(-437.825);
   Graph_Graph2->SetMaximum(418.075);
   Graph_Graph2->SetDirectory(0);
   Graph_Graph2->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph2->SetLineColor(ci);
   Graph_Graph2->GetXaxis()->SetLabelFont(42);
   Graph_Graph2->GetXaxis()->SetLabelOffset(0.03);
   Graph_Graph2->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetXaxis()->SetTickLength(-0.03);
   Graph_Graph2->GetXaxis()->SetTitleOffset(1.35);
   Graph_Graph2->GetXaxis()->SetTitleFont(42);
   Graph_Graph2->GetYaxis()->SetLabelFont(42);
   Graph_Graph2->GetYaxis()->SetLabelOffset(0.02);
   Graph_Graph2->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetYaxis()->SetTickLength(-0.02);
   Graph_Graph2->GetYaxis()->SetTitleFont(42);
   Graph_Graph2->GetZaxis()->SetLabelFont(42);
   Graph_Graph2->GetZaxis()->SetLabelOffset(0.02);
   Graph_Graph2->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetZaxis()->SetTickLength(-0.02);
   Graph_Graph2->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph2);
   
   graph->Draw("p ");
   
   graph = new TGraph(18);
   graph->SetName("Graph");
   graph->SetTitle("Graph");

   ci = TColor::GetColor("#000099");
   graph->SetLineColor(ci);

   ci = TColor::GetColor("#006600");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(22);
   graph->SetPoint(0,394.9500055,24.99999964);
   graph->SetPoint(1,354.1500109,138.3999944);
   graph->SetPoint(2,258.4000034,259.449995);
   graph->SetPoint(3,144.4499973,294.9999976);
   graph->SetPoint(4,19.19997027,319.7999895);
   graph->SetPoint(5,-86.09999752,317.7999991);
   graph->SetPoint(6,-205.5000035,245.7999949);
   graph->SetPoint(7,-293.2500112,127.3999929);
   graph->SetPoint(8,-408.5500047,1.100000381);
   graph->SetPoint(9,-434.1000066,-144.1000011);
   graph->SetPoint(10,-349.9500129,-267.899996);
   graph->SetPoint(11,-212.8000019,-342.8999965);
   graph->SetPoint(12,-141.1999967,-433.5000026);
   graph->SetPoint(13,0.2749999985,-433.799988);
   graph->SetPoint(14,127.1049969,-354.5555);
   graph->SetPoint(15,273.8000033,-295.1999937);
   graph->SetPoint(16,359.1500121,-201.1999942);
   graph->SetPoint(17,401.0500059,-79.44999897);
   
   TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","Graph",100,-517.615,484.565);
   Graph_Graph3->SetMinimum(-509.16);
   Graph_Graph3->SetMaximum(395.16);
   Graph_Graph3->SetDirectory(0);
   Graph_Graph3->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph3->SetLineColor(ci);
   Graph_Graph3->GetXaxis()->SetLabelFont(42);
   Graph_Graph3->GetXaxis()->SetLabelOffset(0.03);
   Graph_Graph3->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph3->GetXaxis()->SetTickLength(-0.03);
   Graph_Graph3->GetXaxis()->SetTitleOffset(1.35);
   Graph_Graph3->GetXaxis()->SetTitleFont(42);
   Graph_Graph3->GetYaxis()->SetLabelFont(42);
   Graph_Graph3->GetYaxis()->SetLabelOffset(0.02);
   Graph_Graph3->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph3->GetYaxis()->SetTickLength(-0.02);
   Graph_Graph3->GetYaxis()->SetTitleFont(42);
   Graph_Graph3->GetZaxis()->SetLabelFont(42);
   Graph_Graph3->GetZaxis()->SetLabelOffset(0.02);
   Graph_Graph3->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph3->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph3->GetZaxis()->SetTickLength(-0.02);
   Graph_Graph3->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph3);
   
   graph->Draw("p ");
   
   graph = new TGraph(18);
   graph->SetName("Graph");
   graph->SetTitle("Graph");

   ci = TColor::GetColor("#000099");
   graph->SetLineColor(ci);

   ci = TColor::GetColor("#006600");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(26);
   graph->SetPoint(0,372.1000061,65.59999847);
   graph->SetPoint(1,327.2000122,188.8999939);
   graph->SetPoint(2,242.8000031,289.3999939);
   graph->SetPoint(3,129.1999969,355);
   graph->SetPoint(4,-2.99000003e-05,377.7999878);
   graph->SetPoint(5,-129.1999969,355);
   graph->SetPoint(6,-242.8000031,289.3999939);
   graph->SetPoint(7,-327.2000122,188.8999939);
   graph->SetPoint(8,-372.1000061,65.59999847);
   graph->SetPoint(9,-372.1000061,-65.59999847);
   graph->SetPoint(10,-327.2000122,-188.8999939);
   graph->SetPoint(11,-242.8000031,-289.3999939);
   graph->SetPoint(12,-129.1999969,-355);
   graph->SetPoint(13,0,-377.7999878);
   graph->SetPoint(14,129.1999969,-355);
   graph->SetPoint(15,242.8000031,-289.3999939);
   graph->SetPoint(16,327.2000122,-188.8999939);
   graph->SetPoint(17,372.1000061,-65.59999847);
   
   TH1F *Graph_Graph4 = new TH1F("Graph_Graph4","Graph",100,-446.52,446.52);
   Graph_Graph4->SetMinimum(-453.36);
   Graph_Graph4->SetMaximum(453.36);
   Graph_Graph4->SetDirectory(0);
   Graph_Graph4->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph4->SetLineColor(ci);
   Graph_Graph4->GetXaxis()->SetLabelFont(42);
   Graph_Graph4->GetXaxis()->SetLabelOffset(0.03);
   Graph_Graph4->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph4->GetXaxis()->SetTickLength(-0.03);
   Graph_Graph4->GetXaxis()->SetTitleOffset(1.35);
   Graph_Graph4->GetXaxis()->SetTitleFont(42);
   Graph_Graph4->GetYaxis()->SetLabelFont(42);
   Graph_Graph4->GetYaxis()->SetLabelOffset(0.02);
   Graph_Graph4->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph4->GetYaxis()->SetTickLength(-0.02);
   Graph_Graph4->GetYaxis()->SetTitleFont(42);
   Graph_Graph4->GetZaxis()->SetLabelFont(42);
   Graph_Graph4->GetZaxis()->SetLabelOffset(0.02);
   Graph_Graph4->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph4->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph4->GetZaxis()->SetTickLength(-0.02);
   Graph_Graph4->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph4);
   
   graph->Draw("p ");
   
   graph = new TGraph(18);
   graph->SetName("Graph");
   graph->SetTitle("Graph");

   ci = TColor::GetColor("#000099");
   graph->SetLineColor(ci);
   graph->SetMarkerColor(2);
   graph->SetMarkerStyle(4);
   graph->SetPoint(0,372.1000061,65.59999847);
   graph->SetPoint(1,327.2000122,188.8999939);
   graph->SetPoint(2,242.8000031,289.3999939);
   graph->SetPoint(3,129.1999969,355);
   graph->SetPoint(4,-2.99000003e-05,377.7999878);
   graph->SetPoint(5,-129.1999969,355);
   graph->SetPoint(6,-242.8000031,289.3999939);
   graph->SetPoint(7,-327.2000122,188.8999939);
   graph->SetPoint(8,-372.1000061,65.59999847);
   graph->SetPoint(9,-372.1000061,-65.59999847);
   graph->SetPoint(10,-327.2000122,-188.8999939);
   graph->SetPoint(11,-242.8000031,-289.3999939);
   graph->SetPoint(12,-129.1999969,-355);
   graph->SetPoint(13,0,-377.7999878);
   graph->SetPoint(14,129.1999969,-355);
   graph->SetPoint(15,242.8000031,-289.3999939);
   graph->SetPoint(16,327.2000122,-188.8999939);
   graph->SetPoint(17,372.1000061,-65.59999847);
   
   TH1F *Graph_Graph5 = new TH1F("Graph_Graph5","Graph",100,-446.52,446.52);
   Graph_Graph5->SetMinimum(-453.36);
   Graph_Graph5->SetMaximum(453.36);
   Graph_Graph5->SetDirectory(0);
   Graph_Graph5->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph5->SetLineColor(ci);
   Graph_Graph5->GetXaxis()->SetLabelFont(42);
   Graph_Graph5->GetXaxis()->SetLabelOffset(0.03);
   Graph_Graph5->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph5->GetXaxis()->SetTickLength(-0.03);
   Graph_Graph5->GetXaxis()->SetTitleOffset(1.35);
   Graph_Graph5->GetXaxis()->SetTitleFont(42);
   Graph_Graph5->GetYaxis()->SetLabelFont(42);
   Graph_Graph5->GetYaxis()->SetLabelOffset(0.02);
   Graph_Graph5->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph5->GetYaxis()->SetTickLength(-0.02);
   Graph_Graph5->GetYaxis()->SetTitleFont(42);
   Graph_Graph5->GetZaxis()->SetLabelFont(42);
   Graph_Graph5->GetZaxis()->SetLabelOffset(0.02);
   Graph_Graph5->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph5->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph5->GetZaxis()->SetTickLength(-0.02);
   Graph_Graph5->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph5);
   
   graph->Draw("p ");
   
   graph = new TGraph(18);
   graph->SetName("Graph");
   graph->SetTitle("Graph");

   ci = TColor::GetColor("#000099");
   graph->SetLineColor(ci);
   graph->SetMarkerColor(4);
   graph->SetMarkerStyle(25);
   graph->SetPoint(0,371.3999939,65.48999786);
   graph->SetPoint(1,326.6000061,188.6000061);
   graph->SetPoint(2,242.3999939,288.8999939);
   graph->SetPoint(3,129,354.3999939);
   graph->SetPoint(4,-2.99000003e-05,377.1000061);
   graph->SetPoint(5,-129,354.3999939);
   graph->SetPoint(6,-242.3999939,288.8999939);
   graph->SetPoint(7,-326.6000061,188.6000061);
   graph->SetPoint(8,-371.3999939,65.48999786);
   graph->SetPoint(9,-371.3999939,-65.48999786);
   graph->SetPoint(10,-326.6000061,-188.6000061);
   graph->SetPoint(11,-242.3999939,-288.8999939);
   graph->SetPoint(12,-129,-354.3999939);
   graph->SetPoint(13,0,-377.1000061);
   graph->SetPoint(14,129,-354.3999939);
   graph->SetPoint(15,242.3999939,-288.8999939);
   graph->SetPoint(16,326.6000061,-188.6000061);
   graph->SetPoint(17,371.3999939,-65.48999786);
   
   TH1F *Graph_Graph6 = new TH1F("Graph_Graph6","Graph",100,-445.68,445.68);
   Graph_Graph6->SetMinimum(-452.52);
   Graph_Graph6->SetMaximum(452.52);
   Graph_Graph6->SetDirectory(0);
   Graph_Graph6->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph6->SetLineColor(ci);
   Graph_Graph6->GetXaxis()->SetLabelFont(42);
   Graph_Graph6->GetXaxis()->SetLabelOffset(0.03);
   Graph_Graph6->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph6->GetXaxis()->SetTickLength(-0.03);
   Graph_Graph6->GetXaxis()->SetTitleOffset(1.35);
   Graph_Graph6->GetXaxis()->SetTitleFont(42);
   Graph_Graph6->GetYaxis()->SetLabelFont(42);
   Graph_Graph6->GetYaxis()->SetLabelOffset(0.02);
   Graph_Graph6->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph6->GetYaxis()->SetTickLength(-0.02);
   Graph_Graph6->GetYaxis()->SetTitleFont(42);
   Graph_Graph6->GetZaxis()->SetLabelFont(42);
   Graph_Graph6->GetZaxis()->SetLabelOffset(0.02);
   Graph_Graph6->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph6->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph6->GetZaxis()->SetTickLength(-0.02);
   Graph_Graph6->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph6);
   
   graph->Draw("p ");
   
   TArc *arc = new TArc(0,0,375.1,0,360);
   arc->SetFillStyle(0);
   arc->Draw();
   TLatex *   tex = new TLatex(-46.21849,-21.35356,"#eta~0");
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(-48.90756,70.53676,"#eta~+1");

   ci = TColor::GetColor("#006600");
   tex->SetTextColor(ci);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(-47.56303,-104.4924,"#eta~-1");
   tex->SetTextColor(4);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(-344.7059,448.308,"TOF shifts in current alignment (#times 50)");
   tex->SetTextSize(0.03733956);
   tex->SetLineWidth(2);
   tex->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
   c1->ToggleToolBar();
}
