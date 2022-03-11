#define USE_TWISS
//hoge
const Int_t       Npar  = 1e3;
const Double_t    EmitX = 2.99e-7;  // m*rad
const Double_t    EmitY = 2.52e-7;  // m*rad

const Double_t    I[3] = { .259, -.719, .009 };  // A

Double_t          BetaXinit  = 1.    , BetaYinit  = 1.    ,
                  AlphaXinit = 0.    , AlphaYinit = 0.    ,
                  GammaXinit = 1.    , GammaYinit = 1.    ;
Double_t          BetaXini2  =  1.172, BetaYini2  =  1.760,
                  AlphaXini2 = -2.355, AlphaYini2 = -3.718;

TMatrixD          MD1(4,4), MD2(4,4), MD3(4,4), MD4(4,4),
                  MQ1(4,4), MQ2(4,4), MQ3(4,4), MM (4,4);

TMatrixD          MSigmaInit(4,4)   , MSigmaEnd(4,4);

TVectorD          Vinit[Npar],  // x, x', y, y'
                  Vend [Npar];  // x, x', y, y'

vector<Double_t>  Xinit(Npar), XPinit(Npar),
                  Yinit(Npar), YPinit(Npar),
                  Xend (Npar), XPend (Npar),
                  Yend (Npar), YPend (Npar);

TCanvas           *canv;
TH2D              *frame[4];
TGraph            *GraphInit[4];  // xx', yy', xy, x'y'
TGraph            *GraphEnd [4];  // xx', yy', xy, x'y'


void init ();
void draw_dist     (Int_t flag /* 0:init, 1:end */);
void calc_sigmamat (Int_t flag /* 0:init, 1:end */);
void def_Mdrift    (TMatrixD &MD, Double_t L /* [m] */);
void def_Mquad     (TMatrixD &MQ, Double_t L /* [m] */, Double_t I /* [A] */);


void hw_hirayamasemi_1 ()
{
  init();

#ifdef USE_TWISS
  TMatrixD  CONV_TWISS(4,4);
#if 0
  CONV_TWISS(0,0) =  .20738280;
  CONV_TWISS(0,1) = 1.06254053;
  CONV_TWISS(1,0) = 1.32331656;
  CONV_TWISS(1,1) = 1.95810593;
  CONV_TWISS(2,2) =  .65479186;
  CONV_TWISS(2,3) = 1.15379705;
  CONV_TWISS(3,2) = 2.03881430;
  CONV_TWISS(3,3) = 2.06535543;
#else
  TMatrixD  Conv_TwissBeta1(4,4),
            Conv_TwissAlph1(4,4),
            Conv_TwissBeta2(4,4),
            Conv_TwissAlph2(4,4);
  Conv_TwissBeta1(0,0) = 1./sqrt(BetaXinit);
  Conv_TwissBeta1(1,1) =    sqrt(BetaXinit);
  Conv_TwissBeta1(2,2) = 1./sqrt(BetaYinit);
  Conv_TwissBeta1(3,3) =    sqrt(BetaYinit);
  Conv_TwissAlph1(0,0) = 1.;
  Conv_TwissAlph1(1,1) = 1.;
  Conv_TwissAlph1(2,2) = 1.;
  Conv_TwissAlph1(3,3) = 1.;
  Conv_TwissAlph1(1,0) = AlphaXinit;
  Conv_TwissAlph1(3,2) = AlphaYinit;
  Conv_TwissBeta2(0,0) = 1./sqrt(BetaXini2);
  Conv_TwissBeta2(1,1) =    sqrt(BetaXini2);
  Conv_TwissBeta2(2,2) = 1./sqrt(BetaYini2);
  Conv_TwissBeta2(3,3) =    sqrt(BetaYini2);
  Conv_TwissAlph2(0,0) = 1.;
  Conv_TwissAlph2(1,1) = 1.;
  Conv_TwissAlph2(2,2) = 1.;
  Conv_TwissAlph2(3,3) = 1.;
  Conv_TwissAlph2(1,0) = AlphaXini2;
  Conv_TwissAlph2(3,2) = AlphaYini2;

  Conv_TwissAlph2.Invert();
  Conv_TwissBeta2.Invert();
  CONV_TWISS = Conv_TwissBeta2*Conv_TwissAlph2*Conv_TwissAlph1*Conv_TwissBeta1;
#endif  // #if 0or1
  cout << "\n<<< conversion matrix for Twiss >>>";
  CONV_TWISS.Print();
#endif

  // initial beam dist
  for (int ipar=0; ipar<Npar; ipar++)
  {
    Vinit[ipar].ResizeTo(4);
    Vinit[ipar][0] = gRandom->Gaus(0, sqrt(BetaXinit *EmitX));
    Vinit[ipar][1] = gRandom->Gaus(0, sqrt(GammaXinit*EmitX));
    Vinit[ipar][2] = gRandom->Gaus(0, sqrt(BetaYinit *EmitY));
    Vinit[ipar][3] = gRandom->Gaus(0, sqrt(GammaYinit*EmitY));
  }

#ifdef USE_TWISS
  for (int ipar=0; ipar<Npar; ipar++)
  {
    Vinit[ipar] = CONV_TWISS*Vinit[ipar];
  }
#endif

  // draw phase space dist.
  draw_dist(0);

  // sigma matrix calc.
  calc_sigmamat(0);

  // define TMatrix
  def_Mdrift(MD1,   2.35e-3);
  def_Mdrift(MD2,  12.25e-3);
  def_Mdrift(MD3,  22.45e-3);
  def_Mdrift(MD4, 314.45e-3);
  def_Mquad (MQ1,  56.00e-3, I[0]);
  def_Mquad (MQ2,  56.00e-3, I[1]);
  def_Mquad (MQ3,  56.00e-3, I[2]);
  //cout << "\n<<< MD1 >>>"; MD1.Print();
  //cout << "\n<<< MQ1 >>> : " << I[0] << " A"; MQ1.Print();
  //cout << "\n<<< MD2 >>>"; MD2.Print();
  //cout << "\n<<< MQ2 >>> : " << I[1] << " A"; MQ2.Print();
  //cout << "\n<<< MD3 >>>"; MD3.Print();
  //cout << "\n<<< MQ3 >>> : " << I[2] << " A"; MQ3.Print();
  //cout << "\n<<< MD4 >>>"; MD4.Print();
  MM = MD4*MQ3*MD3*MQ2*MD2*MQ1*MD1;

  // calc. transport
  for (int ipar=0; ipar<Npar; ipar++)
  {
    Vend[ipar].ResizeTo(4);
    Vend[ipar] = MM*Vinit[ipar];
  }

  // draw phase space dist.
  draw_dist(1);

  // sigma matrix calc.
  calc_sigmamat(1);
}


void init ()
{
  gROOT->SetStyle("ATLAS");
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetMarkerSize(.6);
  gStyle->SetNdivisions(507, "XY");
}


void draw_dist (Int_t flag)
{
  if (!canv)
  {
    canv = new TCanvas("canv", "", 700*1.5, 500*1.5);
    canv->Divide(2,2);
  }

  if (!frame[0])
  {
    frame[0] = new TH2D("frame0", ";x [m];x' [rad]"   , 100, -.004, .004, 100, -.04 , .04 );
    frame[1] = new TH2D("frame1", ";y [m];y' [rad]"   , 100, -.004, .004, 100, -.04 , .04 );
    //frame[2] = new TH2D("frame2", ";x [m];y [m]"      , 100, -.004, .004, 100, -.004, .004);
    frame[2] = new TH2D("frame2", ";x [m];y [m]"      , 100, -.02, .02, 100, -.02, .02);
    frame[3] = new TH2D("frame3", ";x' [rad];y' [rad]", 100, -.01 , .01 , 100, -.01 , .01 );
  }

  if (flag==0)
  {
    for (int ipar=0; ipar<Npar; ipar++)
    {
      Xinit [ipar] = Vinit[ipar][0];
      XPinit[ipar] = Vinit[ipar][1];
      Yinit [ipar] = Vinit[ipar][2];
      YPinit[ipar] = Vinit[ipar][3];
    }
    GraphInit[0] = new TGraph(Npar, &Xinit [0], &XPinit[0]);
    GraphInit[1] = new TGraph(Npar, &Yinit [0], &YPinit[0]);
    GraphInit[2] = new TGraph(Npar, &Xinit [0], &Yinit [0]);
    GraphInit[3] = new TGraph(Npar, &XPinit[0], &YPinit[0]);
  }
  else
  {
    for (int ipar=0; ipar<Npar; ipar++)
    {
      Xend  [ipar] = Vend [ipar][0];
      XPend [ipar] = Vend [ipar][1];
      Yend  [ipar] = Vend [ipar][2];
      YPend [ipar] = Vend [ipar][3];
    }
    GraphEnd [0] = new TGraph(Npar, &Xend  [0], &XPend [0]);
    GraphEnd [1] = new TGraph(Npar, &Yend  [0], &YPend [0]);
    GraphEnd [2] = new TGraph(Npar, &Xend  [0], &Yend  [0]);
    GraphEnd [3] = new TGraph(Npar, &XPend [0], &YPend [0]);
  }

  for (int i=0; i<4; i++)
  {
    canv->cd(i+1);
    frame[i]->Draw("axis");
    if (flag==0)
      GraphInit[i]->Draw("same p");
    else
    {
      GraphEnd [i]->SetMarkerColor(kRed);
      GraphEnd [i]->Draw("same p");
    }
  }
}


void calc_sigmamat (Int_t flag)
{
  TVectorD *Vtmp;
  if (flag==0) Vtmp = Vinit;
  else         Vtmp = Vend ;

  TMatrixD MSigmatmp(4,4);

  for (int ipar=0; ipar<Npar; ipar++)
  {
    MSigmatmp(0,0) += Vtmp[ipar][0]*Vtmp[ipar][0]/Npar;  // xx
    MSigmatmp(0,1) += Vtmp[ipar][0]*Vtmp[ipar][1]/Npar;  // xx'
    MSigmatmp(0,2) += Vtmp[ipar][0]*Vtmp[ipar][2]/Npar;  // xy
    MSigmatmp(0,3) += Vtmp[ipar][0]*Vtmp[ipar][3]/Npar;  // xy'
    MSigmatmp(1,1) += Vtmp[ipar][1]*Vtmp[ipar][1]/Npar;  // x'x'
    MSigmatmp(1,2) += Vtmp[ipar][1]*Vtmp[ipar][2]/Npar;  // x'y
    MSigmatmp(1,3) += Vtmp[ipar][1]*Vtmp[ipar][3]/Npar;  // x'y'
    MSigmatmp(2,2) += Vtmp[ipar][2]*Vtmp[ipar][2]/Npar;  // yy
    MSigmatmp(2,3) += Vtmp[ipar][2]*Vtmp[ipar][3]/Npar;  // yy'
    MSigmatmp(3,3) += Vtmp[ipar][3]*Vtmp[ipar][3]/Npar;  // y'y'
  }
  MSigmatmp(1,0) = MSigmatmp(0,1);
  MSigmatmp(2,0) = MSigmatmp(0,2);
  MSigmatmp(2,1) = MSigmatmp(1,2);
  MSigmatmp(3,0) = MSigmatmp(0,3);
  MSigmatmp(3,1) = MSigmatmp(1,3);
  MSigmatmp(3,2) = MSigmatmp(2,3);

  cout << "\n<<< sigma matrix >>> : " << flag;
  MSigmatmp.Print();

  if (flag==0) MSigmaInit = MSigmatmp;
  else         MSigmaEnd  = MSigmatmp;

  Double_t EmitXtmp = sqrt(MSigmatmp(0,0)*MSigmatmp(1,1)-MSigmatmp(0,1)*MSigmatmp(1,0));
  Double_t EmitYtmp = sqrt(MSigmatmp(2,2)*MSigmatmp(3,3)-MSigmatmp(2,3)*MSigmatmp(3,2));

  cout << "<<< Twiss Parameters >>>" << endl;
  cout << "  Emit _X = " <<  EmitXtmp                << endl;
  cout << "  Alpha_X = " << -MSigmatmp(0,1)/EmitXtmp << endl;
  cout << "  Beta _X = " <<  MSigmatmp(0,0)/EmitXtmp << endl;
  cout << "  Gamma_X = " <<  MSigmatmp(1,1)/EmitXtmp << endl;
  cout << "  Emit _Y = " <<  EmitYtmp                << endl;
  cout << "  Alpha_Y = " << -MSigmatmp(2,3)/EmitYtmp << endl;
  cout << "  Beta _Y = " <<  MSigmatmp(2,2)/EmitYtmp << endl;
  cout << "  Gamma_Y = " <<  MSigmatmp(3,3)/EmitYtmp << endl;
  cout << endl;
}


void def_Mdrift (TMatrixD &MD, Double_t L)
{
  MD(0,0) = 1.; MD(1,1) = 1.;
  MD(2,2) = 1.; MD(3,3) = 1.;
  MD(0,1) = L ; MD(2,3) = L ;
}


void def_Mquad (TMatrixD &MQ, Double_t L, Double_t I)
{
  if (I==0)
  {
    cout << "current error" << endl;
    return;
  }

  Double_t sqrtK = sqrt(abs(I)*.3*.095/297.e-6);
  cout << "K = " << (I/abs(I))*sqrtK*sqrtK << endl;

  if (I>0)  // horizontal focusing
  {
    MQ(0,0) = cos (sqrtK*L);
    MQ(2,2) = cosh(sqrtK*L);
    MQ(0,1) = sin (sqrtK*L)/sqrtK;
    MQ(2,3) = sinh(sqrtK*L)/sqrtK;
    MQ(1,1) = MQ(0,0);
    MQ(3,3) = MQ(2,2);
    MQ(1,0) = -sqrtK*sqrtK*MQ(0,1);
    MQ(3,2) =  sqrtK*sqrtK*MQ(2,3);;
  }
  else      // vertical focusing
  {
    MQ(2,2) = cos (sqrtK*L);
    MQ(0,0) = cosh(sqrtK*L);
    MQ(2,3) = sin (sqrtK*L)/sqrtK;
    MQ(0,1) = sinh(sqrtK*L)/sqrtK;
    MQ(1,1) = MQ(0,0);
    MQ(3,3) = MQ(2,2);
    MQ(1,0) =  sqrtK*sqrtK*MQ(0,1);
    MQ(3,2) = -sqrtK*sqrtK*MQ(2,3);;
  }
}
