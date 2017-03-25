#include "ITensor/all.h"
#include "fstream"
#include "string"


using namespace itensor;
using namespace std;
using std::vector;

inline int Factorial(int x) {
  if (x==0) return 1;
  return (x == 1 ? x : x * Factorial(x - 1));
}
inline int coeff(int d, int y){
  return Factorial(d-1)/(Factorial(y-1)*Factorial(d-y));
}

int main()
  {

  //Real pa1[784][2];
  Real pa1[300][5],pa2[300][5];
  Real sam1x[300][3],sam1y[300][3],sam2x[300][3],sam2y[300][3];
  //Real set1[300][3],set2[300][5];
  int i = 0,train_set = 200, set_N = 300;
  int d = 3;
  Real pi = 3.1415926;
  ifstream file1,file2;
  //char filename[512] = {'s'};
  //string filename ＝ 's.txt';

  //string filename1 = "gaussion1.txt";
  //string filename2 = "gaussion2.txt";
  string filename1 = "circle1.txt";
  string filename2 = "circle2.txt";
  file1.open(filename1,ios::in);
  file2.open(filename2,ios::in);

  ofstream outfile_Bl1,outfile_Bl2,outfile_dots;
  outfile_Bl1.open("Bl1.txt"),outfile_Bl2.open("Bl2.txt"),outfile_dots.open("classification.txt");

    i = 0;
    if(file1.fail())
      {
        cout<<"file not found."<<endl;
        file1.close();
        cin.get();
        cin.get();
      }
    else
    {
      while(!file1.eof()&&i<set_N)   //读取数据到数组,file.eof()判断文件是否为空
        {
          //if (i%2 == 0) file>>pa1[i][0];
          //else file>>pa1[i][1];
          file1>>pa1[i][0];
          file1>>pa1[i][1];
          file1>>pa1[i][2];
          //file1>>pa1[i][5];
          //file1>>pa1[i][6];
          i++;
          //cout<<i<<" "<<pa1[i]<<endl;
        }
    }
      file1.close(); //close file/

      i = 0;
      if(file2.fail())
        {
          cout<<"file not found."<<endl;
          file2.close();
          cin.get();
          cin.get();
        }
      else
      {
        while(!file2.eof()&&i<set_N)   //读取数据到数组,file.eof()判断文件是否为空
          {
            //if (i%2 == 0) file>>pa1[i][0];
            //else file>>pa1[i][1];
            file2>>pa2[i][0];
            file2>>pa2[i][1];
            file2>>pa2[i][2];
            //file2>>pa2[i][5];
            //file2>>pa2[i][6];
            i++;
            //cout<<i<<" "<<pa1[i]<<endl;
          }
      }
        file2.close(); //close file/

cout<<"what0"<<endl;
cout<<coeff(4,2)<<endl;
  for(int j=0;j<300;j++)
    {
      for(int sj = 1; sj<=d; sj++) {
        sam1x[j][sj] = pow(pa1[j][0],sj-1);
        sam1y[j][sj] = pow(pa1[j][1],sj-1);

        sam2x[j][sj] = pow(pa2[j][0],sj-1);
        sam2y[j][sj] = pow(pa2[j][1],sj-1);
      }
    }
cout<<"what1"<<endl;

  for(int j=0;j<100;j++)
  {
    cout<<j<<" "<<pa1[j][0]<<" "<<pa1[j][1]<<" "<<pa1[j][2]<<endl;
    cout<<j<<" "<<pa2[j][0]<<" "<<pa2[j][1]<<" "<<pa2[j][2]<<endl;
  }

  int N = 2;
  Real cutoff = 1E-12;
  auto cost = 0.0;
  auto sites = SpinHalf(N);
  auto state = InitState(sites);          //??
  auto psi = MPS(state);

  auto l =  Index("classification",1);
  auto l1 = Index("classification1",1);
  auto l2 = Index("classification2",1);
  auto s1 = Index("x1",d);
  auto s2= Index("x1",d);
  auto sl = Index("sl",d);
  auto sr = Index("sr",d);

  //auto Bl = ITensor(s1,s2,l);
  auto Bl1 = ITensor(s1,s2,l1);
  auto Bl2 = ITensor(s1,s2,l2);
  auto ls = ITensor(s1,sl);
  auto rs = ITensor(s2,sr);
  //auto core1 = ITensor(ls,rs,l);
  randomize(Bl1);
  randomize(Bl2);
  PrintData(Bl1);
  PrintData(Bl2);
  //cout<<"here"<<endl;
  //auto fl = randomTensor(s1,s2,l);
  auto ts1 = ITensor(s1);
  auto ts2 = ITensor(s2);
      auto ts = ITensor(s1,s2,prime(s1),prime(s2));
  //auto ts1 = ITensor(s1,s2,prime(s1),prime(s2));
  //auto ts2 = ITensor(s1,s2,prime(s1),prime(s2));
  ts.fill(0.0);
  auto dL1 = ITensor(l1);
  auto dL2 = ITensor(l2);
  //auto fl = ITensor(l);
  auto fl1 = ITensor(l1);
  auto fl2 = ITensor(l2);
  auto yn1 = ITensor(s1,s2,l1);
  yn1.fill(0.0);
  auto yn2 = ITensor(s1,s2,l2);
  yn2.fill(0.0);
  cout<<"here0"<<endl;

  //for(int i = 0; i<train_set; i++) {
  for(int i = 0; i<50; i++) {
    //dL.fill(0.0);
    //dL.set(l(pa1[i][2]),1);                     //choose one mode
    dL1.fill(1.0);
    for (int sj = 1; sj<=d; sj++) {
    ts1.set(s1(sj),sam1x[i][sj]);
    ts2.set(s2(sj),sam1y[i][sj]);
    }

    auto psi = ts1*ts2;
    ts = ts + prime(prime(psi,s1),s2)*psi;
    //yn = yn + psi*dL;
    yn1 = yn1 + psi*dL1;

    dL1.fill(0.0);
    for (int sj = 1; sj<=d; sj++) {
    ts1.set(s1(sj),sam2x[i][sj]);
    ts2.set(s2(sj),sam2y[i][sj]);
    }

    psi = ts1*ts2;
    ts = ts + prime(prime(psi,s1),s2)*psi;
    //yn = yn + psi*dL;
    yn1 = yn1 + psi*dL1;

  }

  PrintData(ts);
  PrintData(yn1);

  ITensor S,V;
  auto U = ITensor(s1,s2);
  svd(ts,U,S,V,{"Cutoff",cutoff});        //trick: pesudo inverse
  PrintData(S);
  auto inverse = [](Real r) { if (r < 1E-4) return 1/r; };

  S.apply(inverse);           //how to preserve part of S
  PrintData(S);

  auto ts_inverse = U*S*V;        // how to truncate ts
  auto II = ts*ts_inverse;
  PrintData(II);

  PrintData(ts_inverse);
  PrintData(yn1*ts_inverse);
  Bl1 = yn1*ts_inverse;
  Bl1 = noprime(Bl1);                                     //classifor
  PrintData(Bl1);

  ts.fill(0.0);
  for(int i = 0; i<50; i++) {
    //dL.fill(0.0);
    //dL.set(l(pa2[i][2]),1);
    dL2.fill(1.0);
    for (int sj = 1; sj<=d; sj++) {
    ts1.set(s1(sj),sam2x[i][sj]);
    ts2.set(s2(sj),sam2y[i][sj]);
    }

    auto psi = ts1*ts2;
    ts = ts + prime(prime(psi,s1),s2)*psi;
    //yn =  yn + psi*dL;
    yn2 =  yn2 + psi*dL2;

    dL2.fill(0.0);
    for (int sj = 1; sj<=d; sj++) {
    ts1.set(s1(sj),sam1x[i][sj]);
    ts2.set(s2(sj),sam1y[i][sj]);
    }
    psi = ts1*ts2;
    ts = ts + prime(prime(psi,s1),s2)*psi;
    yn2 = yn2 + psi*dL2;
  }
  PrintData(ts2);
  PrintData(yn2);

  //ITensor S,V;
  U = ITensor(s1,s2);
  svd(ts,U,S,V,{"Cutoff",cutoff});        //trick: pesudo inverse
  PrintData(S);

  S.apply(inverse);
  PrintData(S);

  ts_inverse = U*S*V;
  II = ts*ts_inverse;
  PrintData(II);

  PrintData(yn);
  Bl2 = yn2*ts_inverse;
  Bl2 = noprime(Bl2);                                     //classifor
  PrintData(Bl2);

  auto cn1 = ITensor(s1,s2,l1);
  cn1.fill(0.0);
  for(int i = 0; i<50; i++) {
    //dL.fill(0.0);
    //dL.set(l(pa1[i][2]),1);
    dL1.fill(1.0);
    for (int sj = 1; sj<=d; sj++) {
    ts1.set(s1(sj),sam1x[i][sj]);
    ts2.set(s2(sj),sam1y[i][sj]);
    }
    auto psi = ts1*ts2;
    //cn = cn + (Bl*psi - dL)*psi;
    cn1 = cn1 + (Bl1*psi*psi - psi*dL1);
  }
  PrintData(cn1);

  cout<<"h1"<<endl;

  auto cn2 = ITensor(s1,s2,l2);
  cn2.fill(0.0);
  for(int i = 0; i<50; i++) {
    //dL.fill(0.0);
    //dL.set(l(pa2[i][2]),1);
    dL2.fill(1.0);
    for (int sj = 1; sj<=d; sj++) {
    ts1.set(s1(sj),sam2x[i][sj]);
    ts2.set(s2(sj),sam2y[i][sj]);
    }
    auto psi = ts1*ts2;
    //cn = cn + (Bl*psi - dL)*psi;
    cn2 = cn2 + (Bl2*psi*psi - psi*dL2);
    //PrintData((Bl*psi - dL)*psi);
  }
  PrintData(cn2);


  PrintData(Bl1);
  int bln1 = 0;                                            //saperation line
  for (int bl1 = 1;bl1 <= d;bl1++) {
    for (int bl2 = 1;bl2 <= d;bl2++) {
      //for (int bl3 = 1;bl3 <= 2;bl3++) {
        bln1++;
        outfile_Bl1<<Bl1.real(s1(bl1),s2(bl2),l1(1))<<" ";
      //}
    }
  }
  outfile_Bl1<<endl;

  PrintData(Bl2);
  int bln2 = 0;                                            //saperation line
  for (int bl1 = 1;bl1 <= d;bl1++) {
    for (int bl2 = 1;bl2 <= d;bl2++) {
      //for (int bl3 = 1;bl3 <= 2;bl3++) {
        bln2++;
        outfile_Bl2<<Bl2.real(s1(bl1),s2(bl2),l2(1))<<" ";
      //}
    }
  }
  outfile_Bl2<<endl;

  cout<<"here"<<endl;

  auto cost1 = 0.0;
  //dL1 = ITensor(l1);
  dL1.fill(1.0);                                           //real cost
  for(int j = train_set; j<train_set+50; j++) {
    //dL.fill(0.0);
    //dL.set(l(pa1[j][2]),1);
    ts1.set(s1(1),sam1x[j][1]); ts1.set(s1(2),sam1x[j][2]);
    ts2.set(s2(1),sam1y[j][1]); ts2.set(s2(2),sam1y[j][2]);
    for (int sj = 1; sj<=d; sj++) {
    ts1.set(s1(sj),sam1x[j][sj]);
    ts2.set(s2(sj),sam1y[j][sj]);
    }
    auto fl1 = Bl1*ts1*ts2;
    //cost1 += pow(norm(dL-fl),2);
    cost1 += pow(norm(dL1-fl1),2);
    //PrintData(cost);
  }
  PrintData(cost1);

  cout<<"hh1"<<endl;

  auto cost2 = 0.0;
  //auto dL2 = ITensor(l2);
  dL2.fill(1.0);
  for(int j = train_set; j<train_set+50; j++) {
    //dL.fill(0.0);
    //dL.set(l(pa2[j][2]),1);
    for (int sj = 1; sj<=d; sj++) {
    ts1.set(s1(sj),sam2x[j][sj]);
    ts2.set(s2(sj),sam2y[j][sj]);
  }
    fl2 = Bl2*ts1*ts2;
    //cost += pow(norm(dL-fl),2);
    cost2 += pow(norm(dL2-fl2),2);
  }
  PrintData(cost2);
  cout<<"hh2"<<endl;


  auto outcome = 0;
  auto result = 0;
  //for(int i = train_set; i<set_N; i++) {
  for(int i = train_set; i<train_set+100; i++) {
  //for(int i = 0; i<50; i++) {
  //for(int i = 0; i<100; i++) {

  for (int sj = 1; sj<=d; sj++) {
  ts1.set(s1(sj),sam2x[i][sj]);
  ts2.set(s2(sj),sam2y[i][sj]);
  }

    //fl = Bl*ts1*ts2;
    fl1 = Bl1*ts1*ts2;
    fl2 = Bl2*ts1*ts2;
    outfile_dots<<pa2[i][0]<<" "<<pa2[i][1]<<" "<<pa2[i][2]<<" ";
    //if (fl.real(l(1))>fl.real(l(2))) {
    if (fl1.real(l1(1))>fl2.real(l2(1))) {
      result = 1;
      outfile_dots<<result; }
    else {
      result = 2;
      outfile_dots<<result; }
    if (result == pa2[i][2]) outcome += 1;
    outfile_dots<<endl;

    for (int sj = 1; sj<=d; sj++) {
    ts1.set(s1(sj),sam1x[i][sj]);
    ts2.set(s2(sj),sam1y[i][sj]);
    }

    fl1 = Bl1*ts1*ts2;
    fl2 = Bl2*ts1*ts2;
    outfile_dots<<pa1[i][0]<<" "<<pa1[i][1]<<" "<<pa1[i][2]<<" ";
    result = 0;
    //if (fl.real(l(1))>fl.real(l(2))) {
    if (fl1.real(l1(1))>fl2.real(l2(1))) {
      result = 1;
      outfile_dots<<result; }
    else {
      result = 2;
      outfile_dots<<result; }
    if (result == pa1[i][2]) outcome += 1;
    outfile_dots<<endl;
    }
  outfile_dots<<outcome;
  cout<<"new classification: "<<outcome<<endl;

  outfile_dots.close();
  outfile_Bl1.close();
  outfile_Bl2.close();

    return 0;
  }
