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
  Real sam1x[300][7],sam1y[300][7],sam2x[300][7],sam2y[300][7];
  //Real set1[300][3],set2[300][5];
  int i = 0,train_set = 200, set_N = 300;
  int d = 6;
  Real pi = 3.1415926;
  ifstream file1,file2;
  //char filename[512] = {'s'};
  //string filename ＝ 's.txt';

  //string filename1 = "gaussion1.txt";
  //string filename2 = "gaussion2.txt";
  string filename1 = "class1.txt";
  string filename2 = "class2.txt";
  file1.open(filename1,ios::in);
  file2.open(filename2,ios::in);

  ofstream outfile_Bl,outfile_dots;
  outfile_Bl.open("Bl.txt"),outfile_dots.open("classification.txt");

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
            i++;
            //cout<<i<<" "<<pa1[i]<<endl;
          }
      }
        file2.close(); //close file/

  for(int j=0;j<300;j++)
    {
      for(int sj = 1; sj<=d; sj++) {
        sam1x[j][sj] = sqrt(coeff(d,sj))*pow(cos(pa1[j][0]*(pi)/2),(d-sj))*pow(sin(pa1[j][0]*(pi)/2),(sj-1));
        sam1y[j][sj] = sqrt(coeff(d,sj))*pow(cos(pa1[j][1]*(pi)/2),(d-sj))*pow(sin(pa1[j][1]*(pi)/2),(sj-1));
        cout<<"x: "<<sam1x[j][sj]<<endl;
        cout<<"y: "<<sam1y[j][sj]<<endl;

        sam2x[j][sj] = sqrt(coeff(d,sj))*pow(cos(pa2[j][0]*(pi)/2),(d-sj))*pow(sin(pa2[j][0]*(pi)/2),(sj-1));
        sam2y[j][sj] = sqrt(coeff(d,sj))*pow(cos(pa2[j][1]*(pi)/2),(d-sj))*pow(sin(pa2[j][1]*(pi)/2),(sj-1));
      }
    }

  for(int j=0;j<100;j++)
  {
    cout<<pa1[j][0]<<" "<<pa1[j][1]<<" "<<pa1[j][2]<<" ";
    //cout<<pa1[i][0]<<" "<<pa1[i][1]<<endl;
    //cout<<j<<" "<<pa1[j][0]<<" "<<pa1[j][1]<<" "<<pa1[j][2]<<endl;
    //cout<<j<<" "<<pa2[j][0]<<" "<<pa2[j][1]<<" "<<pa2[j][2]<<endl;
    //cout<<j<<" "<<sam1[j][0]<<" "<<sam1[j][1]<<" "<<sam1[j][2]<<" "<<sam1[j][3]<<" "<<pa1[j][0]<<" "<<pa1[j][1]<<endl;
    //cout<<j<<" "<<sam2[j][0]<<" "<<sam2[j][1]<<" "<<sam2[j][2]<<" "<<sam2[j][3]<<" "<<pa2[j][0]<<" "<<pa2[j][1]<<endl;
  }

  int N = 2;
  Real cutoff = 1E-12;
  auto cost = 0.0;
  auto sites = SpinHalf(N);
  auto state = InitState(sites);          //??
  auto psi = MPS(state);

  auto l = Index("classification",2);
  auto s1 = Index("x1",d);
  auto s2= Index("x1",d);

  auto Bl = ITensor(s1,s2,l);
  randomize(Bl);
  //cout<<"here"<<endl;
  //auto fl = randomTensor(s1,s2,l);
  auto ts1 = ITensor(s1);
  auto ts2 = ITensor(s2);
  auto ts = ITensor(s1,s2,prime(s1),prime(s2));
  ts.fill(0.0);
  auto dL = ITensor(l);
  auto fl = ITensor(l);
  auto yn = ITensor(l,s1,s2);
  yn.fill(0.0);

  cout<<"here0"<<endl;

  //for(int i = 0; i<train_set; i++) {
  for(int i = 0; i<50; i++) {
    dL.fill(0.0);
    dL.set(l(pa1[i][2]),1);

    for (int sj = 1; sj<=d; sj++) {
    ts1.set(s1(sj),sam1x[i][sj]);
    ts2.set(s2(sj),sam1y[i][sj]);
    }

    auto psi = ts1*ts2;
    ts = ts + prime(prime(psi,s1),s2)*psi;
    yn = yn + psi*dL;

    dL.fill(0.0);
    dL.set(l(pa2[i][2]),1);

    for (int sj = 1; sj<=d; sj++) {
    ts1.set(s1(sj),sam2x[i][sj]);
    ts2.set(s2(sj),sam2y[i][sj]);
  }

    psi = ts1*ts2;
    ts = ts + prime(prime(psi,s1),s2)*psi;
    yn =  yn + psi*dL;
  }

  PrintData(ts);
  PrintData(yn);

  ITensor S,V;
  auto U = ITensor(s1,s2);
  svd(ts,U,S,V,{"Cutoff",cutoff});        //trick: pesudo inverse
  PrintData(S);
  auto inverse = [](Real r) { if (r < 1E-4) return 1/r; };

  S.apply(inverse);
  PrintData(S);

  auto ts_inverse = U*S*V;
  auto II = ts*ts_inverse;
  PrintData(II);

  Bl = yn*ts_inverse;
  Bl = noprime(Bl);                                     //classifor
  PrintData(Bl);

  auto cn = ITensor(s1,s2,l);
  cn.fill(0.0);
  for(int i = 0; i<50; i++) {                   //check part

    dL.fill(0.0);
    dL.set(l(pa1[i][2]),1);

    for (int sj = 1; sj<=d; sj++) {
    ts1.set(s1(sj),sam1x[i][sj]);
    ts2.set(s2(sj),sam1y[i][sj]);
    }

    auto psi = ts1*ts2;
    //cout<<"part0"<<endl;
    cn = cn + (Bl*psi - dL)*psi;
    //PrintData((Bl*psi - dL)*psi);
    //cout<<"part1"<<endl;

    dL.fill(0.0);
    dL.set(l(pa2[i][2]),1);

    for (int sj = 1; sj<=d; sj++) {
    ts1.set(s1(sj),sam2x[i][sj]);
    ts2.set(s2(sj),sam2y[i][sj]);
  }

    psi = ts1*ts2;
    cn = cn + (Bl*psi - dL)*psi;
    //PrintData((Bl*psi - dL)*psi);
  }
  PrintData(cn);

  int bln = 0;                                            //saperation line
  for (int bl1 = 1;bl1 <= d;bl1++) {
    for (int bl2 = 1;bl2 <= d;bl2++) {
      for (int bl3 = 1;bl3 <= 2;bl3++) {
        bln++;
        outfile_Bl<<Bl.real(s1(bl1),s2(bl2),l(bl3))<<" ";
      }
    }
  }
  outfile_Bl<<endl;

  cost = 0.0;                                             //real cost
  for(int j = train_set; j<train_set+50; j++) {
    dL.fill(0.0);
    dL.set(l(pa1[j][2]),1);
    ts1.set(s1(1),sam1x[j][1]); ts1.set(s1(2),sam1x[j][2]);
    ts2.set(s2(1),sam1y[j][1]); ts2.set(s2(2),sam1y[j][2]);
    for (int sj = 1; sj<=d; sj++) {
    ts1.set(s1(sj),sam1x[j][sj]);
    ts2.set(s2(sj),sam1y[j][sj]);
    }
    auto fl = Bl*ts1*ts2;
    cost += pow(norm(dL-fl),2);
    //PrintData(cost);

    dL.fill(0.0);
    dL.set(l(pa2[j][2]),1);
    for (int sj = 1; sj<=d; sj++) {
    ts1.set(s1(sj),sam2x[j][sj]);
    ts2.set(s2(sj),sam2y[j][sj]);
  }
    fl = Bl*ts1*ts2;
    cost += pow(norm(dL-fl),2);
  }
  PrintData(cost);


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

    fl = Bl*ts1*ts2;
    outfile_dots<<pa2[i][0]<<" "<<pa2[i][1]<<" "<<pa2[i][2]<<" ";
    if (fl.real(l(1))>fl.real(l(2))) {
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

    fl = Bl*ts1*ts2;
    outfile_dots<<pa1[i][0]<<" "<<pa1[i][1]<<" "<<pa1[i][2]<<" ";
    result = 0;
    if (fl.real(l(1))>fl.real(l(2))) {
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
  outfile_Bl.close();

    return 0;
  }
