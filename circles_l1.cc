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
  int d = 3;
  Real sam1x[300][d+1],sam1y[300][d+1],sam2x[300][d+1],sam2y[300][d+1];
  //Real set1[300][3],set2[300][5];
  int i=0,train_set = 100, set_N = 300;
  Real pi = 3.1415926;
  ifstream file1,file2;
  //char filename[512] = {'s'};
  //string filename ＝ 's.txt';

  //string filename1 = "gaussion1.txt";
  //string filename2 = "gaussion2.txt";

  //string filename1 = "class1.txt";
  //string filename2 = "class2.txt";
  string filename1 = "circle1.txt";
  string filename2 = "circle2.txt";
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

        cout<<"here1"<<endl;
  for(int j=0;j<300;j++)
    {
      for(int sj = 1; sj<=d; sj++) {
        //sam1x[j][sj] = sqrt(coeff(d,sj))*pow(cos(pa1[j][0]*(pi)/2),(d-sj))*pow(sin(pa1[j][0]*(pi)/2),(sj-1));
        //sam1y[j][sj] = sqrt(coeff(d,sj))*pow(cos(pa1[j][1]*(pi)/2),(d-sj))*pow(sin(pa1[j][1]*(pi)/2),(sj-1));

        //sam2x[j][sj] = sqrt(coeff(d,sj))*pow(cos(pa2[j][0]*(pi)/2),(d-sj))*pow(sin(pa2[j][0]*(pi)/2),(sj-1));
        //sam2y[j][sj] = sqrt(coeff(d,sj))*pow(cos(pa2[j][1]*(pi)/2),(d-sj))*pow(sin(pa2[j][1]*(pi)/2),(sj-1));
        sam1x[j][sj] = pow(pa1[j][0],sj-1);
        sam1y[j][sj] = pow(pa1[j][1],sj-1);

        sam2x[j][sj] = pow(pa2[j][0],sj-1);
        sam2y[j][sj] = pow(pa2[j][1],sj-1);
      }
    }

  for(int j=0;j<300;j++)
  {
    //cout<<pa1[j][0]<<" "<<pa1[j][1]<<" "<<pa1[j][2]<<endl;
    //cout<<pa2[j][0]<<" "<<pa2[j][1]<<" "<<pa2[j][2]<<endl;
  }

  int N = 2;
  int nsweeps = 200;
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
  auto dL = ITensor(l);
  auto fl = ITensor(l);

  cout<<"here2"<<endl;

for (int sweep = 0; sweep < nsweeps; sweep++) {

  auto dBl = ITensor(s1,s2,l);
  for(int i = 0; i<train_set; i++) {
  //for(int i = 0; i<50; i++) {
    dL.fill(0.0);
    dL.set(l(pa1[i][2]),1);

    for (int sj = 1; sj<=d; sj++) {
    ts1.set(s1(sj),sam1x[i][sj]);
    ts2.set(s2(sj),sam1y[i][sj]);
    }

    fl = Bl*ts1*ts2;
    dBl += 0.01*ts1*ts2*(dL-fl);              //gradient
    //Bl+=dBl;                            // update, no step?

    dL.fill(0.0);
    dL.set(l(pa2[i][2]),1);

    for (int sj = 1; sj<=d; sj++) {
    ts1.set(s1(sj),sam2x[i][sj]);
    ts2.set(s2(sj),sam2y[i][sj]);
    }

    fl = Bl*ts1*ts2;
    dBl += 0.01*(ts1*ts2*(dL-fl));              //gradient
    //Bl+=dBl;
  }
  Bl += dBl;                //update

  for (int bl1 = 1;bl1 <= d;bl1++) {
    for (int bl2 = 1;bl2 <= d;bl2++) {
      for (int bl3 = 1;bl3 <= 2;bl3++) {
        //outfile_Bl<<Bl.real(s1(bl1),s2(bl2),l(bl3))<<" ";                 // how to fix one index and SVD
      }
    }
  }
  //outfile_Bl<<endl;

  cout<<"here3"<<endl;
  auto cost = 0.0;                                             //real cost
  for(int j = 1; j<train_set; j++) {
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
}

for (int bl1 = 1;bl1 <= d;bl1++) {
  for (int bl2 = 1;bl2 <= d;bl2++) {
    for (int bl3 = 1;bl3 <= 2;bl3++) {
      outfile_Bl<<Bl.real(s1(bl1),s2(bl2),l(bl3))<<" ";                 // how to fix one index and SVD
    }
  }
}
outfile_Bl<<endl;


  cout<<"here4"<<endl;
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
