Double_t DXs(Double_t E0, Double_t y){
  Double_t Z = 18;
  Double_t E = E0 *(1-y);
  Double_t k = E0 *y;
  Double_t r0 = 2.818e-13; //cm
  Double_t m0 = 0.511;

  Double_t sigma;
  sigma = 2.* Z*Z*r0*r0/137./y *1.4/39.948/1.6605e-24;
  Double_t b = 2*E0*E*pow(Z,1/3.)/111./k/m0;
  Double_t M = 1./(pow(k*m0/2./E0/E,2)+pow(pow(Z,1/3.)/111.,2));

  sigma = sigma * ((1+pow(E/E0,2)-2./3.*E/E0) *(log(M)+1-2/b*TMath::ATan(b)) 
		   + E/E0 *(2/pow(b,2)*log(1+b*b)+4*(2.-b*b)/(3.*b*b)*TMath::ATan(b)-8/3./b/b+2/9.));
  return sigma;
}

void plot(){

  TFile *file = new TFile("MFP.root");
  TGraph *gpair = (TGraph*)file->Get("PairProduction");

  Double_t x[10000],y[10][10000];
  Double_t sum[5]={0,0,0,0,0};
  Double_t sum1[5]={0,0,0,0,0};
  for (Int_t i=0;i!=10000;i++){
    x[i]=(i+0.5)/10000.;
    
    y[0][i] = DXs(100,x[i]);
    y[1][i] = DXs(600,x[i]);
    y[2][i] = DXs(1700,x[i]);
    y[3][i] = DXs(3000,x[i]);
    y[4][i] = DXs(7500,x[i]);

    

  }
  TGraph *g1 = new TGraph(10000,x,y[0]);
  TGraph *g2 = new TGraph(10000,x,y[1]);
  TGraph *g3 = new TGraph(10000,x,y[2]);
  TGraph *g4 = new TGraph(10000,x,y[3]);
  TGraph *g5 = new TGraph(10000,x,y[4]);
  g1->Draw("AL");
  g2->Draw("Lsame");
  g3->Draw("Lsame");
  g4->Draw("Lsame");
  g5->Draw("Lsame");
  
  g1->SetLineColor(1);
  g2->SetLineColor(2);
  g3->SetLineColor(4);
  g4->SetLineColor(6);
  g5->SetLineColor(8);

  g1->SetLineWidth(2.0);
  g2->SetLineWidth(2.0);
  g3->SetLineWidth(2.0);
  g4->SetLineWidth(2.0);
  g5->SetLineWidth(2.0);

  g1->SetTitle("Bremsstrahlung");
  g1->GetXaxis()->SetTitle("y=k/E_{e}");
  g1->GetYaxis()->SetTitle("d#sigma/dy (cm^{-1})");

  for (Int_t i = 0; i!=1000;i++){
    Double_t y_min = 3/100.;
    Double_t x1 = y_min + (1-y_min)/1000.*(i+0.5);
    sum[0] += g1->Eval(x1) * (1-y_min)/1000.;
    sum1[0]+= g1->Eval(x1)*(1-y_min)/1000. /(1./gpair->Eval(x1*100) + 1./g1->Eval(x1)/(1-y_min));
    

    y_min = 3./600.;
    x1 = y_min + (1-y_min)/1000.*(i+0.5);
    sum[1] += g2->Eval(x1) * (1-y_min)/1000.;
    sum1[1]+= g2->Eval(x1)*(1-y_min)/1000. /(1./gpair->Eval(x1*600) + 1./g2->Eval(x1)/(1-y_min));

    y_min = 3./1700.;
    x1 = y_min + (1-y_min)/1000.*(i+0.5);
    sum[2] += g3->Eval(x1) * (1-y_min)/1000.;
    sum1[2]+= g3->Eval(x1)*(1-y_min)/1000. /(1./gpair->Eval(x1*1700) + 1./g3->Eval(x1)/(1-y_min));

    y_min = 3./3000.;
    x1 = y_min + (1-y_min)/1000.*(i+0.5);
    sum[3] += g4->Eval(x1) * (1-y_min)/1000.;
    sum1[3]+= g4->Eval(x1)*(1-y_min)/1000. /(1./gpair->Eval(x1*3000) + 1./g4->Eval(x1)/(1-y_min));

    y_min = 3./7500.;
    x1 = y_min + (1-y_min)/1000.*(i+0.5);
    sum[4] += g5->Eval(x1) * (1-y_min)/1000.;
    sum1[4]+= g5->Eval(x1)*(1-y_min)/1000. /(1./gpair->Eval(x1*7500) + 1./g5->Eval(x1)/(1-y_min));
  }
  //  cout << sum[0] << " " << sum[1] << " " << sum[2] << " " << sum[3] << " " << sum[4] << endl;
  cout << 1./sum[0] << " " << 1./sum[1] << " " << 1./sum[2] << " " << 1./sum[3] << " " << 1./sum[4] << endl;
  cout << 1./(sum1[0]/sum[0]) << " " << 1./(sum1[1]/sum[1]) << " " << 1./(sum1[2]/sum[2]) << " " << 1./(sum1[3]/sum[3]) << " " << 1./(sum1[4]/sum[4]) << endl;

  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->SetFillColor(10);
  le1->SetBorderSize(0);
  le1->AddEntry(g1,"100 MeV","l");
  le1->AddEntry(g2,"600 MeV","l");
  le1->AddEntry(g3,"1700 MeV","l");
  le1->AddEntry(g4,"3000 MeV","l");
  le1->AddEntry(g5,"7500 MeV","l");
  le1->Draw();

  //TF1 *f1 = new TF1("f1","4./3.-4./3.*x+x*x");
  //f1->Draw();
}
