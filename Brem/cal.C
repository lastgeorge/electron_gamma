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

void cal(Double_t Emax = 1000., Double_t Emin = 3.0){
  Double_t x[10000],y[10000];
  for (Int_t i=0;i!=10000;i++){
    x[i]=(i+0.5)/10000.;
    y[i] = DXs(Emax,x[i]);
  }
  TGraph *g1 = new TGraph(10000,x,y);

  Double_t sum = 0;
  for (Int_t i = 0; i!=1000;i++){
    Double_t y_min = Emin/Emax;
    Double_t x1 = y_min + (1-y_min)/1000.*(i+0.5);
    sum += g1->Eval(x1) * (1-y_min)/1000.;
  }

  cout << "1/MFP and MFP (cm)" << endl;
  cout << sum << " " << 1./ sum << endl;
  

}
