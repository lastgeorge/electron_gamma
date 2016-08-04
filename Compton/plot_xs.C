void plot_xs(){
  Double_t alpha = 1/137.;
  Double_t mass_e = 0.511; // MeV
  
  Double_t conv_MeV_fm = 197.3; // fm * MeV

  Double_t x[100],y[100];
  
  for (Int_t i=0;i!=100;i++){
   
    Double_t E = 1 + (1000-1)/100.*i; // MeV
    x[i] = E;

    Double_t Xs;
    Xs = 2*3.1415926 * alpha * alpha / mass_e/mass_e / 2. *197.3*197.3;
    
    Double_t factor;
    Double_t m = mass_e;
    factor = (m *(2 *E *(pow(E,3) + 9 *pow(E,2) *m + 8 *E *pow(m,2) + 2 *pow(m,3)) 
		  + pow(2 *E + m,2) *  (pow(E,2) - 2 *E *m - 2 *pow(m,2)) *(log((2 *E + m)/m))))/(pow(E,3)* pow(2 *E + m,2));
    Xs = Xs * factor;
    
    // 1.4 g/cm^3, 39.948 u, 1u = 1.66054e-24 g 
    // electrons 18
    
    Double_t path;
    path = 1./(1.4 / 39.948 / 1.6605e-24 * 18. * Xs * 1e-26); // cm 
    y[i] = path;
  }
  TGraph *g1 = new TGraph(100,x,y);
  g1->Draw("AL");
  g1->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  g1->GetYaxis()->SetTitle("M.F.P (cm)");
  g1->SetTitle("Compton Scattering");
  

  //cout << Xs * 82 << endl;

}
