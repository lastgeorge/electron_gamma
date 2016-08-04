Double_t Phi_1(Double_t delta){
  Double_t result;
  if (delta <=1){
    result = 20.867 - 3.242 * delta + 0.625 * delta * delta;
  }else{
    result = 21.12-4.184*log(delta + 0.952);
  }
  return result;
}

Double_t Phi_2(Double_t delta){
  Double_t result;
  if (delta <=1){
    result = 20.209 - 1.930 * delta -0.086 * delta * delta;
  }else{
    result = 21.12-4.184*log(delta + 0.952);
  }
  return result;
}

Double_t fc(Double_t Z){
  Double_t alpha = 1./137.;
  Double_t result = pow(alpha*Z,2) * (1/(1.0+pow(alpha*Z,2)) + 0.20206 - 0.0369 * pow(alpha * Z,2) + 0.0083 * pow(alpha*Z,4) -0.0020 * pow(alpha*Z,6));
  return result;
}

Double_t F(Double_t Z, Double_t Egamma){
  Double_t result;
  if (Egamma >=50){
    result = 8/3.*log(Z) + 8 * fc(Z);
  }else{
    result = 8/3.*log(Z);
  }
  return result;
}

Double_t delta_min(Double_t Z, Double_t Egamma){
  Double_t result;
  result = 136./pow(Z,1./3.) * 4 * 0.511/Egamma;
  return result;
}

Double_t eps_min(Double_t Z, Double_t Egamma){
  Double_t eps_0 = 0.511/Egamma;
  Double_t a = 136/pow(Z,1/3.)*eps_0 / (exp((42.24-F(Z,Egamma))/8.368)-0.952);
  // cout << a << endl;
  if (1-4.*a>=0){
    Double_t eps_1 = (1-sqrt(1-4.*a))/2.;
    if (eps_0 < eps_1){
      return eps_0;
    }else{
      return eps_1;
    }
  }else{
    return eps_0;
  }
}

Double_t csi(Double_t Z){
  Double_t result = log(1440./pow(Z,2./3.))/(log(183./pow(Z,1./3.))-fc(Z));
  return result;
}

Double_t F_1(Double_t Z, Double_t delta,Double_t Egamma){
  Double_t result;
  result = 3*Phi_1(delta) - Phi_2(delta) - F(Z,Egamma);
  return result;
}

Double_t F_2(Double_t Z, Double_t delta, Double_t Egamma){
  Double_t result;
  result = 3./2.*Phi_1(delta)-1./2.*Phi_2(delta)-F(Z,Egamma);
  return result;
}

Double_t dxs_deps(Double_t Z, Double_t eps,Double_t Egamma){
  Double_t factor;

  // Double_t epsmin = eps_min(Z,Egamma);
  // Double_t deltamin = delta_min(Z,Egamma);
  // Double_t delta = 136./pow(Z,1./3.)*0.511/Egamma/eps/(1.0-eps);

  // Double_t F10 = F_1(Z,deltamin,Egamma);
  // Double_t F20 = F_2(Z,deltamin,Egamma);
  
  // Double_t N1 = pow(0.5-epsmin,2)*F10;
  // Double_t N2 = 1.5*F20;

  // Double_t f1 = 3./pow(0.5-epsmin,3)*pow(0.5-eps,2);
  // Double_t f2 = 1./(0.5-epsmin);

  // Double_t g1 = F_1(Z,delta,Egamma)/F10;
  // Double_t g2 = F_2(Z,delta,Egamma)/F20;

  // factor = 1./137.*pow(2.818,2) * Z * (Z + csi(Z)) * 2./9.*(0.5-epsmin) * (N1*f1*g1+N2*f2*g2); // *10e-26 cm^2
  // // 1.4 /39.948 / 1.6605e-24
  // factor = factor * 2.1105e-4; // cm^-1
  
  Double_t delta = 136./pow(Z,1./3.)*0.511/Egamma/eps/(1.0-eps);
  factor = 1./137.*pow(2.818,2) * Z * (Z+csi(Z));
  factor = factor * ((eps*eps+pow(1-eps,2))*(Phi_1(delta)-F(Z,Egamma)/2.) + 2./3.*eps*(1-eps)*(Phi_2(delta) - F(Z,Egamma)/2.) );
  factor = factor * 2.1105e-4;
  

  return factor;
}


Double_t Compton_xs(Double_t E){
  Double_t alpha = 1/137.;
  Double_t mass_e = 0.511; // MeV
  
  Double_t conv_MeV_fm = 197.3; // fm * MeV
  
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
  
  Double_t result = 1./path;
  return result;
  
  
}


void plot(){
  // Integral needs to be between 0.511/Egamma to 2.511/Egamma ...

  Double_t x[1000],y[1000],y1[1000];
  for (Int_t i=0;i!=1000;i++){
    x[i] = 5 + (10000-5)/1000.*i;//0.+0.5/100.*(i+0.5);
    y[i] = 0;
    y1[i] = 0.;
    for (Int_t j=0;j!=100;j++){
      Double_t eps = (0.511 + (2.511-0.511)/100.*(j+0.5))/x[i];
      y[i] += dxs_deps(18,eps,x[i]);
      eps = (j+0.5)/100.;
      y1[i] += dxs_deps(18,eps,x[i]);
    }
    y[i] *= (2.511-0.511)/x[i]/100.*2;
    
    y1[i] *= 1./100.;
    //y[i] = y[i] / y1[i];
  
    y[i] = (Compton_xs(x[i]) + y[i])/(y1[i]+Compton_xs(x[i]));

    y[i] = Compton_xs(x[i]);
    
    //y[i] = (1-exp(-1. *Compton_xs(x[i]))) + (1-exp(-1*y[i]));   

    //    y[i] = dxs_deps(18,x[i],x[i]);//eps_min(18,x[i]);//Phi_2(x[i]);
  }
  TGraph *g1 = new TGraph(1000,x,y);
  g1->Draw("AL");
  g1->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
  g1->SetTitle("1/L_{MFP} (cm^{-1})");
  // TLine *l1 = new TLine(0.,1./18.,10000.,1./18.);
  // l1->Draw();
  // l1->SetLineColor(2);

  g1->SetTitle("Percentage of MIP-like");

  TGraph *g2 = new TGraph(1000,x,y);
  TGraph *g3 = new TGraph(1000,x,y1);
  TFile *file = new TFile("MFP.root","RECREATE");
  g2->Write("Compton");
  g3->Write("PairProduction");
  file->Write();
  file->Close();

  // TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  // le1->SetFillColor(10);
  // le1->SetBorderSize(0);
  // le1->AddEntry(l1,"14 cm radiation length","l");
  // le1->AddEntry(g1,"Bethe-Heitler","l");
  // le1->Draw();
}

